## =========================================================
## MODULE B -- Omega-structured sim-null for the eMLG climate association
## =========================================================
## Empirical calibration of eMLG-level climate association against population
## structure, used here as a sim-FDR CANDIDATE FILTER (not a per-eMLG p test).
##
## NSIM_TOTAL null "climate" covariates are drawn from the among-population
## covariance Omega (same fixed LD-pruned Omega as the real runs): each is a
## random variable with the population-structure autocorrelation of a real
## climate axis but no causal link to genotypes. For every eMLG we count how many
## of the NSIM_TOTAL structure-null covariates match or beat its REAL climate BF
## (k1 on PC1, k2 on PC2). The candidate set is the eMLGs at the FLOOR: k==0, i.e.
## the real climate BF exceeds ALL NSIM_TOTAL structure nulls for that eMLG.
##
## Why this is a valid (permutation-style) genome-wide FDR, unlike scoping BH to a
## pre-screened set: under the pure-structure null an eMLG's observed BF is
## exchangeable with its NSIM_TOTAL nulls, so P(observed is the strict max) =
## 1/(NSIM_TOTAL+1). Hence the EXPECTED number of floor-survivors under the null
## is N_tested/(NSIM_TOTAL+1), and the FDR of the observed floor set is
## approx [N_tested/(NSIM_TOTAL+1)] / (#observed floor-survivors). Any eMLG ABOVE
## the floor (k>=1) is discarded -- it already loses to >=1 structure null and
## cannot survive FDR. At NSIM=1000 the null expectation (32840/1001~=32.8) equaled
## the observed floor count (32) => FDR~=1 (negative); NSIM=10000 drops the null
## expectation to ~3.3, so a holding survivor count yields a meaningful FDR.
##
## EFFICIENCY / validity: run in BATCHES of BATCH null covariates (each one BayPass
## run on the full 32,840 eMLGs -- validated equivalent to separate runs, and the
## per-marker prior depends on the marker SET not the covariate count). Same MCMC
## -seed each batch => identical core chain reused across batches, so the union of
## batches == one NSIM_TOTAL run. Per-eMLG exceedance counts accumulate across
## batches; each batch's ~4.3 GB raw output is parsed then deleted. A checkpoint
## after every batch makes the ~21 h run resumable. Null covariates are Omega
## draws (axis-agnostic) -> one null set calibrates both PC1 and PC2.
## Approach: Li et al. 2018 Mol Ecol Resour 18:809-824.
##
## Reads : aland_excluded_eMLG/{u_eMLG.geno,omega_mat_omega.out,u_DIEM.size,
##           eMLG_group_order.txt}, aland_excluded_eMLG/PC{1,2}_eMLG_withOmega_summary_betai_reg.out
## Writes: moduleB_BayPass/data/moduleB_eMLG_null.rds  (per-eMLG: real BF1/BF2, exceedance counts
##           k1/k2, empirical p1/p2 = (1+k)/(NSIM+1), null_max; floor flags; + meta)
##         moduleB_BayPass/data/moduleB_eMLG_null_ckpt.rds  (resume checkpoint, removed on success)
##         aland_excluded_eMLG/null/null_b*.env  (the null covariates per batch, kept)
## Run from the formica_hybrid repo root. LONG (~21 h for NSIM_TOTAL=10000).

suppressMessages(library(data.table))
NSIM_TOTAL <- 10000
BATCH      <- 1000
SEED_SIM   <- 2026               # base RNG seed for null-covariate draws (per batch: SEED_SIM+b)
MCMC_SEED  <- 74                 # BayPass MCMC seed (same each batch -> identical core chain)
stopifnot(NSIM_TOTAL %% BATCH == 0); NBATCH <- NSIM_TOTAL / BATCH
D  <- "aland_excluded_eMLG"; ND <- file.path(D, "null"); dir.create(ND, showWarnings = FALSE)
CKPT <- "moduleB_BayPass/data/moduleB_eMLG_null_ckpt.rds"
BAYPASS <- "/Users/petrikem/gitlab/baypass_public-master/sources/g_baypass"
OPT <- sprintf("-nthreads 6 -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed %d", MCMC_SEED)

grp <- readLines(file.path(D, "eMLG_group_order.txt")); NM <- length(grp)
b1  <- fread(file.path(D, "PC1_eMLG_withOmega_summary_betai_reg.out"))$`BF(dB)`
b2  <- fread(file.path(D, "PC2_eMLG_withOmega_summary_betai_reg.out"))$`BF(dB)`
stopifnot(length(b1) == NM, length(b2) == NM)
for (f in c("u_eMLG.geno", "omega_mat_omega.out", "u_DIEM.size"))
  file.copy(file.path(D, f), file.path(ND, f), overwrite = TRUE)

Omega <- as.matrix(fread(file.path(D, "omega_mat_omega.out"))); Omega <- (Omega + t(Omega)) / 2
eg <- eigen(Omega, symmetric = TRUE); vals <- pmax(eg$values, 0); P <- nrow(Omega)
draw_null <- function(seed, n) { set.seed(seed)
  sapply(seq_len(n), function(k) as.numeric(scale(eg$vectors %*% (sqrt(vals) * rnorm(P))))) }

## resume from checkpoint if present
if (file.exists(CKPT)) {
  ck <- readRDS(CKPT); k1 <- ck$k1; k2 <- ck$k2; nmax <- ck$nmax; done <- ck$done
  message("resuming: ", done, " of ", NBATCH, " batches already done")
} else { k1 <- integer(NM); k2 <- integer(NM); nmax <- rep(-Inf, NM); done <- 0L }

for (b in (done + 1L):NBATCH) {
  t0 <- Sys.time()
  Y <- draw_null(SEED_SIM + b, BATCH)
  ef <- file.path(ND, sprintf("null_b%02d.env", b))
  write.table(t(Y), ef, quote = FALSE, row.names = FALSE, col.names = FALSE)
  pref <- file.path(ND, sprintf("b%02d", b))
  st <- system(paste0(BAYPASS, " -countdatafile ", ND, "/u_eMLG.geno -omegafile ", ND,
    "/omega_mat_omega.out -efile ", ef, " -poolsizefile ", ND, "/u_DIEM.size ",
    OPT, " -outprefix ", pref, " > ", pref, "_stdout.log 2>&1"))
  stopifnot("BayPass batch failed" = st == 0)
  of <- paste0(pref, "_summary_betai_reg.out")
  nb <- fread(of, select = c("COVARIABLE", "MRK", "BF(dB)"))
  stopifnot(nrow(nb) == BATCH * NM, all(nb$COVARIABLE[1:NM] == 1),
            all(nb$MRK[1:NM] == seq_len(NM)), nb$COVARIABLE[NM + 1] == 2)
  M <- matrix(nb$`BF(dB)`, nrow = NM, ncol = BATCH); rm(nb)
  k1 <- k1 + rowSums(M >= b1)
  k2 <- k2 + rowSums(M >= b2)
  nmax <- pmax(nmax, apply(M, 1, max)); rm(M); invisible(gc())
  ## drop the batch's large raw outputs, keep the .env
  file.remove(Sys.glob(paste0(pref, "_summary_*.out")))
  done <- b
  saveRDS(list(k1 = k1, k2 = k2, nmax = nmax, done = done, group_id = grp), CKPT)
  message(sprintf("batch %d/%d done in %.1f min (cumulative nulls=%d)",
                  b, NBATCH, as.numeric(difftime(Sys.time(), t0, units = "mins")), b * BATCH))
}

## ---- final per-eMLG summary + floor-survivor candidate set --------------
res <- data.table(group_id = grp, BF1 = b1, BF2 = b2, k1 = k1, k2 = k2,
                  p1 = (1 + k1) / (NSIM_TOTAL + 1), p2 = (1 + k2) / (NSIM_TOTAL + 1),
                  null_max = nmax)
res[, `:=`(floor1 = (BF1 >= 15) & (k1 == 0), floor2 = (BF2 >= 15) & (k2 == 0))]
exp_null <- NM / (NSIM_TOTAL + 1)
attr(res, "meta") <- list(NSIM = NSIM_TOTAL, batch = BATCH, seed_sim = SEED_SIM,
  run = "aland_excluded/withOmega", p_floor = 1 / (NSIM_TOTAL + 1),
  null_exp_floor = exp_null,
  note = "candidate = floor survivor (k==0); FDR ~= null_exp_floor / #survivors")
saveRDS(res, "moduleB_BayPass/data/moduleB_eMLG_null.rds")

n1 <- res[floor1 == TRUE, .N]; n2 <- res[floor2 == TRUE, .N]
cat(sprintf("\n=== NSIM=%d sim-FDR filter ===\n", NSIM_TOTAL))
cat(sprintf("expected floor-survivors under pure-structure null: %.2f (per axis, of BF>=15 sets)\n", exp_null))
cat(sprintf("PC1 floor-survivors: %d  -> set FDR ~= %.3f\n", n1, ifelse(n1 > 0, exp_null / n1, NA)))
cat(sprintf("PC2 floor-survivors: %d  -> set FDR ~= %.3f\n", n2, ifelse(n2 > 0, exp_null / n2, NA)))
cat("\nfloor-survivor candidates:\n")
print(rbind(res[floor1 == TRUE, .(group_id, axis = "PC1", BF = round(BF1, 1), k = k1)],
            res[floor2 == TRUE, .(group_id, axis = "PC2", BF = round(BF2, 1), k = k2)])[order(axis, -BF)])
if (file.exists(CKPT)) file.remove(CKPT)
cat("\nSaved moduleB_BayPass/data/moduleB_eMLG_null.rds\n")
