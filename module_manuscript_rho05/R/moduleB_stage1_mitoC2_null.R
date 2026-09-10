## =========================================================
## Stage-1-unit BayPass mitotype C2 -- Omega-structured null calibration
## =========================================================
## Adapted from moduleB_stage1_S1units_null.R (the PC1/PC2 climate null) to
## calibrate the mitotype C2 contrast (moduleB_stage1_prepare_mito_contrast.R,
## run_baypass_stage1.sh step 4) the same way: floor-survivor / sim-FDR test
## against 10,000 Omega-structured null "contrasts".
##
## KEY DIFFERENCE FROM THE CLIMATE NULL: the real mitotype test is a BINARY
## group contrast (-contrastfile, C2 statistic -- Olazcuaga et al. 2020), not
## a continuous covariate regression (-efile, BF). A structure-preserving null
## contrast therefore can't just be a random continuous draw -- it must be a
## population-level +1/-1 PARTITION with the SAME group sizes as the real
## contrast (7 Faquilonia-like vs 12 Fpolyctena-like populations here; read
## from u.mito_contrast rather than hardcoded).
##
## Reuses the EXISTING 10,000 Omega-eigenvector draws already on disk
## (null/null_b01.env .. null_b50.env, from moduleB_stage1_S1units_null.R) --
## no new random draws. Each row (one structured "climate-like" draw across
## the 19 populations) is rank-thresholded to a matching binary partition:
## the top 7 populations by value -> +1, the remaining 12 -> -1. Because these
## draws come from Omega's own eigendecomposition, populations with high
## covariance (incl. the two documented shared-origin pairs) tend to land on
## the same side of the threshold more often than chance -- the same
## structure-preserving property the continuous PC1/PC2 null already has,
## carried over to a binary partition instead of imposed by relabelling.
##
## BayPass's core model computes C2 for MULTIPLE contrasts in one call when
## -contrastfile has multiple rows (BayPass manual section 2.2.4), exactly
## mirroring how -efile batches multiple covariates as columns -- so this
## batches BATCH=200 null contrasts per call, identical sizing to the climate
## null (same MCMC seed, same mini2 16GB-RAM batch size).
##
## PERSISTS each batch's null log10(1/pval) matrix (bf_matrices/mitoC2_bf_b##.rds)
## -- the climate null (moduleB_stage1_S1units_null.R) did NOT persist its BF
## matrices, which meant moduleC_stage1_null_regen.R had to re-run BayPass from
## scratch (~9-10h) to get anything beyond exceedance counts. Persisting here
## from the start avoids repeating that; ~1.5GB total, trivial next to the
## BayPass compute cost.
##
## Designed to run standalone on a machine with just this directory's files
## (no formica_hybrid repo / LDscnR needed -- pure data.table + BayPass CLI).
## RESOURCE NOTE: shares the machine with moduleC_stage1_null_regen.R (climate
## Module C), which already uses -nthreads 10 (all cores) for ~9-10h. Do NOT
## launch this in parallel with that job -- queue it to start after, or reduce
## NTHREADS for both if running concurrently is unavoidable.
##
## Reads (all in the SAME directory as this script, or set D below):
##   u_S1units.geno, omega_mat_omega.out, u_DIEM.size, u.mito_contrast,
##   S1units_group_order.txt, mito_C2_S1units_summary_contrast.out (observed),
##   null/null_b01.env .. null_b50.env (reused, NOT regenerated)
## Writes: moduleB_stage1_mitoC2_null.rds, moduleB_stage1_mitoC2_null_ckpt.rds
##   (resume checkpoint, removed on success),
##   null/nullc2_b*.contrast (kept), null/bf_matrices/mitoC2_bf_b##.rds (kept)
##
## Run:  Rscript moduleB_stage1_mitoC2_null.R
## =========================================================

suppressMessages(library(data.table))
NSIM_TOTAL <- 10000
BATCH      <- 200     # identical batching to the climate null (mini2 16GB RAM)
MCMC_SEED  <- 74       # SAME seed as the climate null and the real mito_C2 run
NTHREADS   <- 10
stopifnot(NSIM_TOTAL %% BATCH == 0); NBATCH <- NSIM_TOTAL / BATCH

D  <- "."
ND <- file.path(D, "null")
BFDIR <- file.path(ND, "bf_matrices"); dir.create(BFDIR, showWarnings = FALSE, recursive = TRUE)
CKPT <- file.path(D, "moduleB_stage1_mitoC2_null_ckpt.rds")
OUT  <- file.path(D, "moduleB_stage1_mitoC2_null.rds")
BAYPASS <- path.expand("~/baypass_public/sources/g_baypass")
stopifnot(file.exists(BAYPASS))
## NOTE: no -nocovscaling -- the real single mito_C2 run (run_baypass_stage1.sh
## step 4) omitted it too (that flag is covariate-mode scaling, not contrast mode).
OPT <- sprintf("-nthreads %d -nval 500 -burnin 5000 -thin 25 -seed %d", NTHREADS, MCMC_SEED)

grp <- readLines(file.path(D, "S1units_group_order.txt")); NM <- length(grp)
obs <- fread(file.path(D, "mito_C2_S1units_summary_contrast.out"))
stopifnot(nrow(obs) == NM, all(obs$MRK == seq_len(NM)))
obs_C2 <- obs$C2; obs_log10p <- obs$`log10(1/pval)`

mito_contrast <- scan(file.path(D, "u.mito_contrast"), quiet = TRUE)
P <- length(mito_contrast)
N_POS <- sum(mito_contrast == 1L); N_NEG <- sum(mito_contrast == -1L)
stopifnot("u.mito_contrast must be all +-1 (no 0-excluded populations)" = N_POS + N_NEG == P)
message(sprintf("mitotype contrast: %d populations, %d Faquilonia-like (+1) / %d Fpolyctena-like (-1)",
                P, N_POS, N_NEG))

ENV_FILES <- file.path(ND, sprintf("null_b%02d.env", seq_len(NBATCH)))
stopifnot("some climate-null .env files missing -- run moduleB_stage1_S1units_null.R first" =
            all(file.exists(ENV_FILES)))
for (f in c("u_S1units.geno", "omega_mat_omega.out", "u_DIEM.size"))
  if (!file.exists(file.path(ND, f))) file.copy(file.path(D, f), file.path(ND, f), overwrite = TRUE)

## rank-threshold one structured draw (length P) to a +-1 partition matching
## the real contrast's group sizes: top N_POS populations by value -> +1.
mk_null_contrast_row <- function(v) {
  ord <- order(v, decreasing = TRUE)
  out <- rep(-1L, length(v)); out[ord[seq_len(N_POS)]] <- 1L
  out
}

if (file.exists(CKPT)) {
  ck <- readRDS(CKPT); k3 <- ck$k3; nmax3 <- ck$nmax3; done <- ck$done
  message("resuming: ", done, " of ", NBATCH, " batches already done")
} else { k3 <- integer(NM); nmax3 <- rep(-Inf, NM); done <- 0L }

for (b in (done + 1L):NBATCH) {
  t0 <- Sys.time()
  cf_out <- file.path(BFDIR, sprintf("mitoC2_bf_b%02d.rds", b))
  cf     <- file.path(ND, sprintf("nullc2_b%02d.contrast", b))
  pref   <- file.path(ND, sprintf("c2b%02d", b))
  of     <- paste0(pref, "_summary_contrast.out")

  if (file.exists(cf_out)) {
    message(sprintf("  batch %d: reduce from persisted null log10p matrix (no BayPass)", b))
    Mp <- readRDS(cf_out)
    stopifnot(nrow(Mp) == NM, ncol(Mp) == BATCH, all(is.finite(Mp)))
  } else {
    ## structured null contrasts for this batch, from the ALREADY-DRAWN climate
    ## null covariates (same 200 draws null_b##.env used for PC1/PC2) -- rank-
    ## thresholded row-wise to a +-1 partition of matching group size.
    Yb <- as.matrix(fread(ENV_FILES[b], header = FALSE))   # BATCH x P
    stopifnot(nrow(Yb) == BATCH, ncol(Yb) == P)
    Cmat <- t(apply(Yb, 1, mk_null_contrast_row))            # BATCH x P, +-1
    write.table(Cmat, cf, quote = FALSE, row.names = FALSE, col.names = FALSE)

    reuse <- file.exists(of) && file.info(of)$size > 0
    if (reuse) {
      message(sprintf("  reusing existing raw output for batch %d (skipping BayPass)", b))
    } else {
      st <- system(paste0(BAYPASS, " -countdatafile ", ND, "/u_S1units.geno -omegafile ", ND,
        "/omega_mat_omega.out -contrastfile ", cf, " -poolsizefile ", ND, "/u_DIEM.size ",
        OPT, " -outprefix ", pref, " > ", pref, "_stdout.log 2>&1"))
      stopifnot("BayPass batch failed" = st == 0)
    }

    nb <- fread(of, select = c("CONTRAST", "MRK", "C2", "log10(1/pval)"))
    stopifnot(nrow(nb) == BATCH * NM, all(nb$CONTRAST[1:NM] == 1),
              all(nb$MRK[1:NM] == seq_len(NM)),
              if (BATCH > 1) nb$CONTRAST[NM + 1] == 2 else TRUE,
              all(is.finite(nb$C2)), all(is.finite(nb$`log10(1/pval)`)))
    Mp <- matrix(nb$`log10(1/pval)`, nrow = NM, ncol = BATCH); rm(nb)

    saveRDS(Mp, paste0(cf_out, ".tmp")); file.rename(paste0(cf_out, ".tmp"), cf_out)
  }

  k3 <- k3 + rowSums(Mp >= obs_log10p)
  nmax3 <- pmax(nmax3, apply(Mp, 1, max)); rm(Mp); invisible(gc())
  file.remove(Sys.glob(paste0(pref, "_summary_*.out")))
  done <- b
  saveRDS(list(k3 = k3, nmax3 = nmax3, done = done, group_id = grp), CKPT)
  message(sprintf("batch %d/%d done in %.1f min (cumulative nulls=%d)",
                  b, NBATCH, as.numeric(difftime(Sys.time(), t0, units = "mins")), b * BATCH))
}

res <- data.table(group_id = grp, C2 = obs_C2, log10p = obs_log10p, k3 = k3,
                  p3 = (1 + k3) / (NSIM_TOTAL + 1), null_max = nmax3)
res[, floor3 := (log10p >= 3) & (k3 == 0)]
exp_null <- NM / (NSIM_TOTAL + 1)
attr(res, "meta") <- list(NSIM = NSIM_TOTAL, batch = BATCH,
  run = "Stage-1-direct S1units (n_snps>=5) mitotype C2, aland_excluded/withOmega",
  null_source = "rank-thresholded climate-null draws (null_b01..b50.env), group sizes matching u.mito_contrast",
  n_pos = N_POS, n_neg = N_NEG,
  p_floor = 1 / (NSIM_TOTAL + 1), null_exp_floor = exp_null,
  note = "candidate = floor survivor (log10p>=3 AND k3==0); FDR ~= null_exp_floor / #survivors")
saveRDS(res, OUT)

n3 <- res[floor3 == TRUE, .N]
cat(sprintf("\n=== NSIM=%d sim-FDR filter (Stage-1 units mitoC2, n=%d) ===\n", NSIM_TOTAL, NM))
cat(sprintf("expected floor-survivors under pure-structure null: %.2f\n", exp_null))
cat(sprintf("mitoC2 floor-survivors: %d  -> set FDR ~= %.3f\n", n3, ifelse(n3 > 0, exp_null / n3, NA)))
cat("\nfloor-survivor candidates:\n")
print(res[floor3 == TRUE, .(group_id, C2 = round(C2, 1), log10p = round(log10p, 2), k = k3)][order(-log10p)])
if (file.exists(CKPT)) file.remove(CKPT)
cat("\nSaved moduleB_stage1_mitoC2_null.rds\n")
