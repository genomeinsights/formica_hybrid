## =========================================================
## Stage-1-unit BayPass climate scan -- Omega-structured null calibration
## =========================================================
## Adapted from moduleB_climate_GEA/R/moduleB_eMLG_null.R (the original
## Stage-2 eMLG null, 32,854 units, ~21h) for the new Stage-1-direct scan
## (18,361 units, n_snps>=5 -- see moduleB_stage1_prepare_baypass_inputs.R).
## Mechanism unchanged: NSIM_TOTAL null "climate" covariates drawn from the
## Stage-1-derived Omega's eigendecomposition (structure-preserving, no
## causal genotype link), run in BATCHES of 1000 (one BayPass call per
## batch, all 1000 null covariates tested at once), same MCMC seed each
## batch. Per-unit exceedance counts (k1 on PC1, k2 on PC2) accumulate
## across batches; floor-survivor = BF>=15 AND k==0.
##
## Designed to run standalone on a machine with just this directory's files
## (no formica_hybrid repo / LDscnR needed -- pure data.table + BayPass CLI).
##
## Reads (all in the SAME directory as this script, or set D below):
##   u_S1units.geno, omega_mat_omega.out, u_DIEM.size, S1units_group_order.txt,
##   PC1_S1units_withOmega_summary_betai_reg.out,
##   PC2_S1units_withOmega_summary_betai_reg.out
## Writes: moduleB_stage1_S1units_null.rds, moduleB_stage1_S1units_null_ckpt.rds
##   (resume checkpoint, removed on success), null/null_b*.env (kept)
##
## Run:  Rscript moduleB_stage1_S1units_null.R
## =========================================================

suppressMessages(library(data.table))
NSIM_TOTAL <- 10000
BATCH      <- 200    # smaller than the original eMLG null's 1000 -- mini2 has only 16GB RAM,
                     # and peak memory for a batch (BayPass's own MCMC arrays plus the R-side
                     # M <- matrix(...BATCH*NM...) parse) was never measured at the original
                     # batch size on a 16GB machine, so this trades a few more (50 vs 10)
                     # cheaper batches for a much smaller peak footprint. Total null draws
                     # (10,000) and expected total compute are unchanged either way.
SEED_SIM   <- 2026
MCMC_SEED  <- 74
NTHREADS   <- 10
stopifnot(NSIM_TOTAL %% BATCH == 0); NBATCH <- NSIM_TOTAL / BATCH

D  <- "."
ND <- file.path(D, "null"); dir.create(ND, showWarnings = FALSE)
CKPT <- file.path(D, "moduleB_stage1_S1units_null_ckpt.rds")
BAYPASS <- path.expand("~/baypass_public/sources/g_baypass")
stopifnot(file.exists(BAYPASS))
OPT <- sprintf("-nthreads %d -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed %d", NTHREADS, MCMC_SEED)

grp <- readLines(file.path(D, "S1units_group_order.txt")); NM <- length(grp)
b1  <- fread(file.path(D, "PC1_S1units_withOmega_summary_betai_reg.out"))$`BF(dB)`
b2  <- fread(file.path(D, "PC2_S1units_withOmega_summary_betai_reg.out"))$`BF(dB)`
stopifnot(length(b1) == NM, length(b2) == NM)
for (f in c("u_S1units.geno", "omega_mat_omega.out", "u_DIEM.size"))
  file.copy(file.path(D, f), file.path(ND, f), overwrite = TRUE)

Omega <- as.matrix(fread(file.path(D, "omega_mat_omega.out"))); Omega <- (Omega + t(Omega)) / 2
eg <- eigen(Omega, symmetric = TRUE); vals <- pmax(eg$values, 0); P <- nrow(Omega)
draw_null <- function(seed, n) { set.seed(seed)
  sapply(seq_len(n), function(k) as.numeric(scale(eg$vectors %*% (sqrt(vals) * rnorm(P))))) }

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
  st <- system(paste0(BAYPASS, " -countdatafile ", ND, "/u_S1units.geno -omegafile ", ND,
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
  file.remove(Sys.glob(paste0(pref, "_summary_*.out")))
  done <- b
  saveRDS(list(k1 = k1, k2 = k2, nmax = nmax, done = done, group_id = grp), CKPT)
  message(sprintf("batch %d/%d done in %.1f min (cumulative nulls=%d)",
                  b, NBATCH, as.numeric(difftime(Sys.time(), t0, units = "mins")), b * BATCH))
}

res <- data.table(group_id = grp, BF1 = b1, BF2 = b2, k1 = k1, k2 = k2,
                  p1 = (1 + k1) / (NSIM_TOTAL + 1), p2 = (1 + k2) / (NSIM_TOTAL + 1),
                  null_max = nmax)
res[, `:=`(floor1 = (BF1 >= 15) & (k1 == 0), floor2 = (BF2 >= 15) & (k2 == 0))]
exp_null <- NM / (NSIM_TOTAL + 1)
attr(res, "meta") <- list(NSIM = NSIM_TOTAL, batch = BATCH, seed_sim = SEED_SIM,
  run = "Stage-1-direct S1units (n_snps>=5), aland_excluded/withOmega",
  p_floor = 1 / (NSIM_TOTAL + 1), null_exp_floor = exp_null,
  note = "candidate = floor survivor (k==0); FDR ~= null_exp_floor / #survivors")
saveRDS(res, file.path(D, "moduleB_stage1_S1units_null.rds"))

n1 <- res[floor1 == TRUE, .N]; n2 <- res[floor2 == TRUE, .N]
cat(sprintf("\n=== NSIM=%d sim-FDR filter (Stage-1 units, n=%d) ===\n", NSIM_TOTAL, NM))
cat(sprintf("expected floor-survivors under pure-structure null: %.2f (per axis, of BF>=15 sets)\n", exp_null))
cat(sprintf("PC1 floor-survivors: %d  -> set FDR ~= %.3f\n", n1, ifelse(n1 > 0, exp_null / n1, NA)))
cat(sprintf("PC2 floor-survivors: %d  -> set FDR ~= %.3f\n", n2, ifelse(n2 > 0, exp_null / n2, NA)))
cat("\nfloor-survivor candidates:\n")
print(rbind(res[floor1 == TRUE, .(group_id, axis = "PC1", BF = round(BF1, 1), k = k1)],
            res[floor2 == TRUE, .(group_id, axis = "PC2", BF = round(BF2, 1), k = k2)])[order(axis, -BF)])
if (file.exists(CKPT)) file.remove(CKPT)
cat("\nSaved moduleB_stage1_S1units_null.rds\n")
