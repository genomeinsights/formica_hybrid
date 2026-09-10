## =========================================================
## Stage-1-direct Module C -- regenerate the null BF matrix and reduce each
## covariate to genome-wide statistics ON THE FLY (persisting BF this time)
## =========================================================
## Adapted from moduleC_climate_vs_sorting/R/moduleC_null_regen.R for the
## Stage-1-direct universe (18,361 units, n_snps>=5). moduleB_stage1_S1units_null.R
## (already run on mini2) kept only per-unit exceedance counts (k1/k2), same as
## the original eMLG null -- the full BF matrix was deleted after each batch. The
## within-covariate rank statistics this Stage-1 Module C needs are NOT
## recoverable from those counts, but ARE exactly regenerable: BayPass ran with a
## fixed MCMC seed (74) on 50 preserved null covariate files (null/null_b01.env
## .. null_b50.env), so re-running BayPass on those exact files reproduces the
## identical 10,000 nulls (Monte-Carlo equivalent, not bit-identical -- BayPass
## is not bit-reproducible; see the equivalence gate below).
##
## SINGLE UNIVERSE (no min grid). Unlike the canonical eMLG Module C (min in
## {5,10}), the Stage-1-direct scan tests only the natural n_snps>=5 universe --
## there is no second min level to sweep. MINS is fixed to one element (5) so the
## SAME shared grid code (moduleC_stat_functions.R's cell_key/minC_stamp) can be
## reused unmodified; the grid is effectively just the tau series {0.5,0.6,0.8}.
##
## BATCH=200 (not 1000), matching moduleB_stage1_S1units_null.R's mini2 sizing
## (16GB RAM). PERSISTS the per-batch BF matrix this time (BFDIR/cRegen_bf_b##.rds,
## ~30MB x 50 batches =~1.5GB) so any future annotation/threshold change is
## re-reducible without re-running BayPass.
##
## Designed to run standalone on mini2 in the SAME directory as the completed
## moduleB_stage1_S1units_null.R run (baypass_stage1_S1units/), reusing its
## null/null_b*.env files, PC1/PC2 *_summary_betai_reg.out, and
## moduleB_stage1_S1units_null.rds (k1/k2 cross-check). Needs a LOCAL copy of
## moduleC_stat_functions.R and moduleC_stage1_annotations.rds alongside it
## (no formica_hybrid repo / LDscnR needed).
##
## Resumable. PILOT vs FULL via MODC_NBATCH (see moduleC_null_regen.R). LONG (~9-10h).
##
## Reads (all in the SAME directory as this script, or set D below):
##   u_S1units.geno, omega_mat_omega.out, u_DIEM.size, S1units_group_order.txt,
##   PC{1,2}_S1units_withOmega_summary_betai_reg.out,
##   null/null_b01.env .. null_b50.env, moduleB_stage1_S1units_null.rds,
##   moduleC_stat_functions.R, moduleC_stage1_annotations.rds
## Writes: moduleC_stage1_null_stats.rds (final; when done==50)
##         moduleC_stage1_null_ckpt.rds  (resume checkpoint, removed on success)
##         null/bf_matrices/cRegen_bf_b##.rds  (persisted null BF, KEPT)
##
## Run:  Rscript moduleC_stage1_null_regen.R
## =========================================================

suppressMessages({ library(data.table); library(digest) })
source("moduleC_stat_functions.R")

## ---- parameters (must match moduleB_stage1_S1units_null.R) --------------
NSIM_TOTAL <- 10000
BATCH      <- 200
NBATCH     <- NSIM_TOTAL / BATCH                       # 50
MCMC_SEED  <- 74
NTHREADS   <- 10
TOL_BF     <- 1e-6
BAYPASS    <- path.expand("~/baypass_public/sources/g_baypass")
OPT        <- sprintf("-nthreads %d -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed %d", NTHREADS, MCMC_SEED)
D          <- "."
ND         <- file.path(D, "null")
STATFNS    <- file.path(D, "moduleC_stat_functions.R")
OBS_PC1    <- file.path(D, "PC1_S1units_withOmega_summary_betai_reg.out")
OBS_PC2    <- file.path(D, "PC2_S1units_withOmega_summary_betai_reg.out")
ANN_FILE   <- file.path(D, "moduleC_stage1_annotations.rds")
CKPT       <- file.path(D, "moduleC_stage1_null_ckpt.rds")
OUT        <- file.path(D, "moduleC_stage1_null_stats.rds")
BFDIR      <- file.path(ND, "bf_matrices")
NBATCH_RUN <- as.integer(Sys.getenv("MODC_NBATCH", NBATCH))
stopifnot(file.exists(BAYPASS), NBATCH_RUN >= 1, NBATCH_RUN <= NBATCH)
dir.create(BFDIR, showWarnings = FALSE, recursive = TRUE)

## ---- annotations + observed BF -------------------------------------------
ann <- readRDS(ANN_FILE)
grp <- readLines(file.path(D, "S1units_group_order.txt"))
stopifnot("annotation order != BayPass order" = identical(ann$group_id, grp))

## SINGLE universe: MINS fixed to one element (the Stage-1 n_snps>=5 universe
## IS the primary/only min level -- no second min to sweep). TAU series unchanged.
TAUS   <- MODULEC_TAU_SERIES                     # c(0.5, 0.6, 0.8)
MINS   <- 5L
TSTAMP <- tauC_stamp(TAUS)
PRIMARY_CELL <- cell_key(MINS, MODULEC_TAU_PRIMARY)   # "min05_tau06"
NM  <- nrow(ann)
STAT_NAMES <- covariate_stat_names(); NSTAT <- length(STAT_NAMES)

## pdir_diff (proportion directional | differentiated, WITHIN the top BF fraction) is
## DEFINED as NA when a null draw's top fraction contains zero differentiated units --
## moduleC_stat_functions.R's own compute_covariate_stats() documents this as expected
## (0/0), not an error. For the canonical 32,854-eMLG universe (~large differentiated
## top-0.1% pool) this apparently never triggered; for the smaller 18,361-unit Stage-1
## universe it does. Every OTHER statistic remains an exact hard stop -- in particular
## the actual PRIMARY sorting stat (sort_gap_differentiated) is never in this list.
PDIR_COLS   <- grep("_pdir_diff$", STAT_NAMES, value = TRUE)
STRICT_COLS <- setdiff(STAT_NAMES, PDIR_COLS)

stopifnot("annotation missing n_loci" = "n_loci" %in% names(ann),
          "all units must satisfy the single min threshold" = all(ann$n_loci >= MINS))
idx_by_min <- setNames(list(seq_len(NM)), minC_stamp(MINS))
message(sprintf("[regen-S1] single universe min_n_loci >= %d : %d Stage-1 units", MINS, NM))

CELLS <- as.data.table(expand.grid(m = MINS, tau = TAUS))[, key := cell_key(m, tau)]
A_cell   <- setNames(vector("list", nrow(CELLS)), CELLS$key)
cell_idx <- setNames(vector("list", nrow(CELLS)), CELLS$key)
for (i in seq_len(nrow(CELLS))) {
  idx <- idx_by_min[[minC_stamp(CELLS$m[i])]]
  cell_idx[[CELLS$key[i]]] <- idx
  A_cell[[CELLS$key[i]]]   <- prepare_annotation_ranks(ann[idx], dir_col = dir_col_for_tau(CELLS$tau[i]))
}
stopifnot("primary cell missing from grid" = PRIMARY_CELL %in% names(A_cell))

## observed per-unit BF (BayPass row order); reduced per cell over its row-subset
b1 <- fread(OBS_PC1)$`BF(dB)`; b2 <- fread(OBS_PC2)$`BF(dB)`
stopifnot(length(b1) == NM, length(b2) == NM,
          "PC1 BF != annotation eBF1 beyond tolerance" = max(abs(b1 - ann$eBF1)) <= TOL_BF,
          "PC2 BF != annotation eBF2 beyond tolerance" = max(abs(b2 - ann$eBF2)) <= TOL_BF,
          all(is.finite(b1)), all(is.finite(b2)))
obs_list <- setNames(lapply(names(A_cell), function(k) {
  idx <- cell_idx[[k]]
  rbind(PC1 = compute_covariate_stats(b1[idx], A_cell[[k]]),
        PC2 = compute_covariate_stats(b2[idx], A_cell[[k]]))
}), names(A_cell))
obs <- obs_list[[PRIMARY_CELL]]

## saved moduleB_stage1_S1units_null.R exceedance counts (cross-check target)
mbnull <- readRDS(file.path(D, "moduleB_stage1_S1units_null.rds"))
stopifnot("Stage-1 null has duplicate group_id" = !any(duplicated(mbnull$group_id)),
          "Stage-1 null group set != unit universe" = setequal(mbnull$group_id, grp))
mo <- match(grp, mbnull$group_id)
k1_saved <- mbnull$k1[mo]; k2_saved <- mbnull$k2[mo]
stopifnot(all(is.finite(k1_saved)), all(is.finite(k2_saved)))

## ---- input fingerprints (checkpoint integrity) ---------------------------
fp_file <- function(f) unname(tools::md5sum(f))
ENV_FILES <- file.path(ND, sprintf("null_b%02d.env", seq_len(NBATCH)))
stopifnot("some null .env files missing" = all(file.exists(ENV_FILES)))
FP_ANN_COLS <- c("group_id", "DI", "recomb", "prop_fixed", "uni_score",
                 paste0("directional_", TSTAMP), "differentiated", "n_loci")
fingerprint <- list(
  ann      = digest(ann[, ..FP_ANN_COLS], algo = "md5"),
  statfns  = fp_file(STATFNS),
  obs_pc1  = fp_file(OBS_PC1), obs_pc2 = fp_file(OBS_PC2),
  env      = vapply(ENV_FILES, fp_file, character(1)),
  geno     = fp_file(file.path(D, "u_S1units.geno")),
  omega    = fp_file(file.path(D, "omega_mat_omega.out")),
  poolsize = fp_file(file.path(D, "u_DIEM.size")),
  params   = list(NSIM = NSIM_TOTAL, BATCH = BATCH, NM = NM, seed = MCMC_SEED, opt = OPT,
                  tau_series = TAUS, tau_primary = MODULEC_TAU_PRIMARY,
                  min_series = MINS, min_primary = MINS, cells = CELLS$key),
  stat_names = STAT_NAMES)

## ---- resume from checkpoint or initialise --------------------------------
if (file.exists(CKPT)) {
  ck <- readRDS(CKPT)
  stopifnot("checkpoint fingerprint mismatch: inputs/definitions changed since it was written" =
              identical(ck$fingerprint, fingerprint),
            identical(ck$group_id, grp))
  null_list <- ck$null_list; k1r <- ck$k1r; k2r <- ck$k2r; done <- ck$done
  stopifnot("checkpoint cell set != current grid" = identical(names(null_list), names(A_cell)))
  message("[regen-S1] resuming: ", done, "/", NBATCH, " batches already done (fingerprint OK)")
} else {
  null_list <- setNames(replicate(length(A_cell),
                 matrix(NA_real_, nrow = NSIM_TOTAL, ncol = NSTAT,
                        dimnames = list(NULL, STAT_NAMES)), simplify = FALSE), names(A_cell))
  k1r <- integer(NM); k2r <- integer(NM); done <- 0L
}

if (done >= NBATCH_RUN) {
  message("[regen-S1] already at or beyond requested batch ", NBATCH_RUN, "; nothing to run")
} else for (b in (done + 1L):NBATCH_RUN) {
  t0   <- Sys.time()
  pref   <- file.path(ND, sprintf("cRegen_b%02d", b))
  of     <- paste0(pref, "_summary_betai_reg.out")
  bf_out <- file.path(BFDIR, sprintf("cRegen_bf_b%02d.rds", b))

  if (file.exists(bf_out)) {
    message(sprintf("  [regen-S1] batch %d: reduce from persisted BF matrix (no BayPass)", b))
    M <- readRDS(bf_out)
    stopifnot("persisted BF matrix wrong shape" = nrow(M) == NM && ncol(M) == BATCH,
              "non-finite BF in persisted matrix" = all(is.finite(M)))
  } else {
    ef <- ENV_FILES[b]
    reuse <- file.exists(of) && file.info(of)$size > 0
    if (reuse) {
      message(sprintf("  [regen-S1] reusing existing raw output for batch %d (skipping BayPass)", b))
    } else {
      st <- system(paste0(BAYPASS, " -countdatafile ", ND, "/u_S1units.geno -omegafile ", ND,
        "/omega_mat_omega.out -efile ", ef, " -poolsizefile ", ND, "/u_DIEM.size ",
        OPT, " -outprefix ", pref, " > ", pref, "_stdout.log 2>&1"))
      stopifnot("BayPass batch failed" = st == 0)
    }

    nb <- fread(of, select = c("COVARIABLE", "MRK", "BF(dB)"))
    row_cov <- rep(seq_len(BATCH), each = NM)
    covnum  <- suppressWarnings(as.integer(nb$COVARIABLE)); okcov <- !is.na(covnum)
    stopifnot("wrong row count" = nrow(nb) == BATCH * NM,
              "MRK ordering wrong (not covariate-major)" = all(nb$MRK == rep(seq_len(NM), times = BATCH)),
              "COVARIABLE (where parseable) disagrees with covariate-major order" =
                all(covnum[okcov] == row_cov[okcov]),
              "unparseable COVARIABLE outside the expected >=1000 field-overflow" =
                all(row_cov[!okcov] >= 1000),
              "non-finite BF in batch" = all(is.finite(nb$`BF(dB)`)))
    if (any(!okcov)) message(sprintf("  [regen-S1] %d rows had overflowed COVARIABLE ('***', cov idx >= 1000); validated via MRK/row position",
                                     sum(!okcov)))
    M <- matrix(nb$`BF(dB)`, nrow = NM, ncol = BATCH); rm(nb, covnum, okcov, row_cov); invisible(gc())

    saveRDS(M, paste0(bf_out, ".tmp")); file.rename(paste0(bf_out, ".tmp"), bf_out)
  }

  k1r <- k1r + rowSums(M >= b1)
  k2r <- k2r + rowSums(M >= b2)

  rows <- ((b - 1L) * BATCH + 1L):(b * BATCH)
  for (k in names(A_cell)) {
    idx <- cell_idx[[k]]
    Mk <- t(apply(M[idx, , drop = FALSE], 2, compute_covariate_stats, A = A_cell[[k]]))
    null_list[[k]][rows, ] <- Mk
    if (!all(is.finite(Mk[, STRICT_COLS])))
      stop(sprintf("non-finite statistic outside *_pdir_diff in batch %d, cell %s", b, k))
    n_na_pdir <- sum(is.na(Mk[, PDIR_COLS]))
    if (n_na_pdir > 0)
      message(sprintf("  [regen-S1] batch %d, cell %s: %d NA pdir_diff value(s) (empty differentiated top fraction; expected, kept NA)",
                      b, k, n_na_pdir))
  }
  rm(M); invisible(gc())

  file.remove(Sys.glob(paste0(pref, "_summary_*.out")))
  done <- b
  saveRDS(list(null_list = null_list, k1r = k1r, k2r = k2r, done = done,
               fingerprint = fingerprint, group_id = grp), CKPT)
  message(sprintf("[regen-S1] batch %d/%d done in %.1f min (cumulative nulls=%d)",
                  b, NBATCH, as.numeric(difftime(Sys.time(), t0, units = "mins")), b * BATCH))
}

## ---- pilot cross-check ----------------------------------------------------
frac <- done / NBATCH
cat(sprintf("\n=== reproduction cross-check after %d/%d batches ===\n", done, NBATCH))
cat(sprintf("PC1: sum(partial k1)=%d ; %.3f x saved (%.0f)  [batch subset, ~1.0 expected]\n",
            sum(k1r), sum(k1r) / (frac * sum(k1_saved)), frac * sum(k1_saved)))
cat(sprintf("PC2: sum(partial k2)=%d ; %.3f x saved (%.0f)\n",
            sum(k2r), sum(k2r) / (frac * sum(k2_saved)), frac * sum(k2_saved)))

## ---- finalise ONLY when all batches are done ------------------------------
if (done == NBATCH) {
  COR_THR <- 0.99; SUM_TOL <- 0.03
  cor_k1 <- cor(k1r, k1_saved, method = "pearson")
  cor_k2 <- cor(k2r, k2_saved, method = "pearson")
  rel1 <- sum(k1r) / sum(k1_saved) - 1
  rel2 <- sum(k2r) / sum(k2_saved) - 1
  d1 <- max(abs(k1r - k1_saved)); d2 <- max(abs(k2r - k2_saved))
  reproduced <- (cor_k1 > COR_THR && cor_k2 > COR_THR &&
                 abs(rel1) < SUM_TOL && abs(rel2) < SUM_TOL)
  cat(sprintf("\nMonte-Carlo equivalence gate (thresholds: r > %.2f, |sum ratio - 1| < %.2f):\n",
              COR_THR, SUM_TOL))
  cat(sprintf("  PC1: Pearson r = %.5f   sum ratio - 1 = %+.4f   (max|dk1| = %d)\n", cor_k1, rel1, d1))
  cat(sprintf("  PC2: Pearson r = %.5f   sum ratio - 1 = %+.4f   (max|dk2| = %d)\n", cor_k2, rel2, d2))
  if (!reproduced)
    stop(sprintf("MONTE-CARLO EQUIVALENCE GATE FAILED (PC1 r=%.4f rel=%+.4f; PC2 r=%.4f rel=%+.4f; ",
                 cor_k1, rel1, cor_k2, rel2),
         sprintf("thresholds r>%.2f, |rel|<%.2f). ", COR_THR, SUM_TOL),
         "Not writing moduleC_stage1_null_stats.rds. Do NOT relax thresholds -- inputs are ",
         "fingerprint-verified, so investigate covariate wiring / BayPass setup; only ",
         "recalibrate tolerances with additional independent reruns if genuinely warranted.")

  n_na_pdir_by_cell <- setNames(integer(length(null_list)), names(null_list))
  for (k in names(null_list)) {
    if (nrow(null_list[[k]]) != NSIM_TOTAL)                          stop(sprintf("%s null_stats wrong nrow", k))
    if (!all(is.finite(null_list[[k]][, STRICT_COLS])))              stop(sprintf("%s null_stats has non-finite entries outside *_pdir_diff", k))
    if (sum(complete.cases(null_list[[k]][, STRICT_COLS])) != NSIM_TOTAL) stop(sprintf("%s not exactly 10,000 complete rows (excluding *_pdir_diff)", k))
    if (!all(is.finite(obs_list[[k]])))                              stop(sprintf("%s observed has non-finite entries", k))
    if (!setequal(colnames(obs_list[[k]]), STAT_NAMES))              stop(sprintf("%s observed missing a statistic", k))
    n_na_pdir_by_cell[k] <- sum(is.na(null_list[[k]][, PDIR_COLS]))
  }
  cat("\npdir_diff NA counts (expected, empty-differentiated-top-fraction; excluded from FDR family):\n")
  print(n_na_pdir_by_cell)

  by_cell <- setNames(lapply(names(A_cell), function(k)
    list(observed = obs_list[[k]], null_stats = null_list[[k]])), names(A_cell))

  res <- list(
    observed   = obs, null_stats = null_list[[PRIMARY_CELL]],
    by_cell     = by_cell, grid = CELLS,
    tau_series  = TAUS, tau_primary = MODULEC_TAU_PRIMARY,
    min_series  = MINS, min_primary = MINS,
    primary_cell = PRIMARY_CELL,
    n_units_by_min = setNames(lengths(idx_by_min), minC_stamp(MINS)),
    stat_names = STAT_NAMES, pdir_cols = PDIR_COLS,
    n_na_pdir_by_cell = n_na_pdir_by_cell,
    k_check    = list(k1r = k1r, k2r = k2r, k1_saved = k1_saved, k2_saved = k2_saved,
                      cor_k1 = cor_k1, cor_k2 = cor_k2, rel1 = rel1, rel2 = rel2,
                      cor_thr = COR_THR, sum_tol = SUM_TOL,
                      max_abs_dk1 = d1, max_abs_dk2 = d2, reproduced = reproduced),
    params = list(NSIM = NSIM_TOTAL, batch = BATCH, mcmc_seed = MCMC_SEED,
                  config = "Stage-1-direct (n_snps>=5) / aland_excluded / withOmega", top_fracs = TOP_FRACS,
                  primary_stats = PRIMARY_STATS, opt = OPT, tol_bf = TOL_BF,
                  tau_series = TAUS, tau_primary = MODULEC_TAU_PRIMARY,
                  min_series = MINS, min_primary = MINS,
                  bf_matrices = normalizePath(BFDIR, mustWork = FALSE)),
    fingerprint = fingerprint,
    session = sessionInfo()
  )
  saveRDS(res, OUT)
  if (file.exists(CKPT)) file.remove(CKPT)
  cat(sprintf("\n[regen-S1] DONE. wrote %s (observed + %d null covariates x %d statistics x %d cells: %s; primary %s)\n",
              OUT, NSIM_TOTAL, NSTAT, length(A_cell), paste(names(A_cell), collapse = "/"), PRIMARY_CELL))
  cat(sprintf("[regen-S1] persisted %d null BF matrices in %s (kept; re-reducible without BayPass)\n",
              length(Sys.glob(file.path(BFDIR, "cRegen_bf_b*.rds"))), BFDIR))
} else {
  cat(sprintf("\n[regen-S1] PILOT/partial: %d/%d batches done.\n", done, NBATCH))
  cat("        Re-run with MODC_NBATCH=50 to complete and write the final object.\n")
}
