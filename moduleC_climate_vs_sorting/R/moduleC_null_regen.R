## =========================================================
## MODULE C (revised) -- regenerate the null BF matrix and reduce each covariate
##                       to genome-wide statistics ON THE FLY
## =========================================================
## Module B kept only per-eMLG exceedance counts from the 10,000 Omega-structured
## null covariates; the full 32,840 x 10,000 null BF matrix was deleted after each
## batch (moduleB_eMLG_null.R line 92). The within-covariate rank statistics this
## revised Module C needs cannot be recovered from those counts. They ARE exactly
## regenerable, because the 10 null covariate files (aland_excluded_eMLG/null/
## null_b01.env .. null_b10.env) survived and BayPass ran with a fixed MCMC seed:
## re-running BayPass on those exact files reproduces the identical 10,000 nulls.
##
## This script re-runs BayPass batch by batch on those .env files and, for each of
## the 10,000 null covariates, reduces its genome-wide BF vector to the same
## statistics computed for the observed PC1/PC2 (via moduleC_stat_functions.R).
##
## TAU SERIES (one BayPass run, three thresholds). The null BF matrix is climate-only
## and IDENTICAL across the sorting-threshold series; only the directional partition
## changes. So each null BF vector is reduced against ALL of tau = {0.5, 0.6, 0.8} in
## the same pass, producing one null-statistic matrix per tau (`by_tau`). tau06 is the
## primary and is also exposed as the top-level `observed`/`null_stats`.
##
## PERSISTED BF MATRICES. Unlike the first design, the per-batch 32,840 x 1,000 BF
## matrix is SAVED (BFDIR/cRegen_bf_b##.rds) rather than deleted, so any future
## annotation/threshold change can be re-reduced WITHOUT re-running BayPass. These
## files are large (~1.5-2 GB total) and live under the git-ignored aland tree.
## Per-eMLG partial exceedance counts are still accumulated so the run is
## cross-checked against Module B's saved k1/k2 (faithful-reproduction gate).
##
## FAITHFUL REGENERATION IS MANDATORY, via a MONTE-CARLO EQUIVALENCE gate. BayPass
## is NOT bit-reproducible (fresh MCMC realization each run -- moduleC_determinism_
## probe.R), so exact k equality is impossible. Instead, per axis, we require
## Pearson r(k_regen, k_saved) > 0.99 AND |sum(k_regen)/sum(k_saved) - 1| < 0.03; if
## not, the run STOPS and no final object is written. The genome-wide statistics are
## separately shown reproducible despite the per-BF noise (moduleC_stat_stability_
## probe.R). Input identity, ordering, completeness and finiteness stay EXACT hard stops.
##
## CHECKPOINT INTEGRITY. The checkpoint records fingerprints of every input that
## defines the statistics (annotation values, moduleC_stat_functions.R, observed
## BayPass outputs, the ten .env files, geno/omega/poolsize, and the numeric
## parameters). On resume these must match EXACTLY, so batches computed under
## different definitions can never be silently combined.
##
## Resumable. PILOT vs FULL is controlled by MODC_NBATCH:
##   MODC_NBATCH=1  -> run batch 1 only (~2 h), validate pipeline, no final object
##   MODC_NBATCH=10 -> run through batch 10 (resumes from checkpoint), write final
##
## Reads : aland_excluded_eMLG/null/{u_eMLG.geno,omega_mat_omega.out,u_DIEM.size,
##            null_b01.env .. null_b10.env}
##         aland_excluded_eMLG/PC{1,2}_eMLG_withOmega_summary_betai_reg.out (observed)
##         moduleB_climate_GEA/data/moduleB_eMLG_null.rds  (saved k1/k2 for the exact check)
##         moduleC_climate_vs_sorting/data/moduleC_annotations.rds
## Writes: moduleC_climate_vs_sorting/data/moduleC_null_stats.rds  (final; when done==10)
##         moduleC_climate_vs_sorting/data/moduleC_null_ckpt.rds   (resume checkpoint)
##         aland_excluded_eMLG/null/bf_matrices/cRegen_bf_b##.rds  (persisted null BF, KEPT)
##         aland_excluded_eMLG/null/cRegen_b*  (BayPass raw outputs, deleted per batch)
## Run from the formica_hybrid repo root. LONG.

suppressMessages({ library(data.table); library(digest) })
source("moduleC_climate_vs_sorting/R/moduleC_stat_functions.R")

## ---- parameters (must match moduleB_eMLG_null.R) ------------------------
NSIM_TOTAL <- 10000
BATCH      <- 1000
NBATCH     <- NSIM_TOTAL / BATCH                       # 10
MCMC_SEED  <- 74                                       # identical each batch (Module B)
TOL_BF     <- 1e-6                                     # observed-vs-eBF tolerance
BAYPASS    <- "/Users/petrikem/gitlab/baypass_public-master/sources/g_baypass"
OPT        <- sprintf("-nthreads 6 -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed %d", MCMC_SEED)
ND         <- "aland_excluded_eMLG/null"
STATFNS    <- "moduleC_climate_vs_sorting/R/moduleC_stat_functions.R"
OBS_PC1    <- "aland_excluded_eMLG/PC1_eMLG_withOmega_summary_betai_reg.out"
OBS_PC2    <- "aland_excluded_eMLG/PC2_eMLG_withOmega_summary_betai_reg.out"
CKPT       <- "moduleC_climate_vs_sorting/data/moduleC_null_ckpt.rds"
OUT        <- "moduleC_climate_vs_sorting/data/moduleC_null_stats.rds"
BFDIR      <- file.path(ND, "bf_matrices")                    # persisted null BF matrices (KEPT)
NBATCH_RUN <- as.integer(Sys.getenv("MODC_NBATCH", NBATCH))   # how far to run THIS invocation
stopifnot(file.exists(BAYPASS), NBATCH_RUN >= 1, NBATCH_RUN <= NBATCH)
dir.create(BFDIR, showWarnings = FALSE, recursive = TRUE)

## ---- annotations + observed BF ------------------------------------------
ann <- readRDS("moduleC_climate_vs_sorting/data/moduleC_annotations.rds")
grp <- readLines(file.path("aland_excluded_eMLG", "eMLG_group_order.txt"))
stopifnot("annotation order != BayPass order" = identical(ann$group_id, grp))

## (min_n_loci x tau) GRID. The null BF is generated once over the primary (>=min)
## universe; each null covariate is reduced against every grid cell. tau varies the
## directional partition; min varies the analysis universe (a strict ROW-SUBSET of
## eMLGs with n_loci >= min -- valid because Omega is fixed, so per-eMLG BF does not
## depend on the size threshold). DI/recomb/magnitude/orientation ranks depend only
## on min; sort_gap_* depend on both. Cell key = "min05_tau06" etc.
TAUS   <- MODULEC_TAU_SERIES                     # c(0.5, 0.6, 0.8)
MINS   <- MODULEC_MIN_SERIES                     # c(5, 10)
TSTAMP <- tauC_stamp(TAUS)
PRIMARY_CELL <- cell_key(MODULEC_MIN_PRIMARY, MODULEC_TAU_PRIMARY)   # "min05_tau06"
NM  <- nrow(ann)
STAT_NAMES <- covariate_stat_names(); NSTAT <- length(STAT_NAMES)

## row-subset per min threshold (min = smallest MUST be the full universe)
stopifnot("annotation missing n_loci" = "n_loci" %in% names(ann),
          "all eMLGs must satisfy the smallest min threshold" = all(ann$n_loci >= min(MINS)))
idx_by_min <- setNames(lapply(MINS, function(m) which(ann$n_loci >= m)), minC_stamp(MINS))
stopifnot("primary (smallest) min is not the full universe" =
            length(idx_by_min[[minC_stamp(min(MINS))]]) == NM,
          "a min subset is empty" = all(lengths(idx_by_min) > 0))
for (m in MINS) message(sprintf("[regen] min_n_loci >= %d : %d eMLGs", m, length(idx_by_min[[minC_stamp(m)]])))

## per-cell annotation ranks (built on the min-subset, with that tau's directional col)
CELLS <- as.data.table(expand.grid(m = MINS, tau = TAUS))[, key := cell_key(m, tau)]
A_cell   <- setNames(vector("list", nrow(CELLS)), CELLS$key)
cell_idx <- setNames(vector("list", nrow(CELLS)), CELLS$key)
for (i in seq_len(nrow(CELLS))) {
  idx <- idx_by_min[[minC_stamp(CELLS$m[i])]]
  cell_idx[[CELLS$key[i]]] <- idx
  A_cell[[CELLS$key[i]]]   <- prepare_annotation_ranks(ann[idx], dir_col = dir_col_for_tau(CELLS$tau[i]))
}
stopifnot("primary cell missing from grid" = PRIMARY_CELL %in% names(A_cell))

## observed per-eMLG BF (BayPass row order); reduced per cell over its row-subset
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
obs <- obs_list[[PRIMARY_CELL]]                  # primary alias (min05_tau06)

## saved Module B exceedance counts (assert one-to-one, recover canonical order)
mbnull <- readRDS("moduleB_climate_GEA/data/moduleB_eMLG_null.rds")
stopifnot("Module B null has duplicate group_id" = !any(duplicated(mbnull$group_id)),
          "Module B null group set != eMLG universe" = setequal(mbnull$group_id, grp))
mo <- match(grp, mbnull$group_id)
k1_saved <- mbnull$k1[mo]; k2_saved <- mbnull$k2[mo]
stopifnot(all(is.finite(k1_saved)), all(is.finite(k2_saved)))

## ---- input fingerprints (checkpoint integrity) --------------------------
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
  geno     = fp_file(file.path(ND, "u_eMLG.geno")),
  omega    = fp_file(file.path(ND, "omega_mat_omega.out")),
  poolsize = fp_file(file.path(ND, "u_DIEM.size")),
  params   = list(NSIM = NSIM_TOTAL, BATCH = BATCH, NM = NM, seed = MCMC_SEED, opt = OPT,
                  tau_series = TAUS, tau_primary = MODULEC_TAU_PRIMARY,
                  min_series = MINS, min_primary = MODULEC_MIN_PRIMARY, cells = CELLS$key),
  stat_names = STAT_NAMES)

## ---- resume from checkpoint or initialise -------------------------------
if (file.exists(CKPT)) {
  ck <- readRDS(CKPT)
  stopifnot("checkpoint fingerprint mismatch: inputs/definitions changed since it was written" =
              identical(ck$fingerprint, fingerprint),
            identical(ck$group_id, grp))
  null_list <- ck$null_list; k1r <- ck$k1r; k2r <- ck$k2r; done <- ck$done
  stopifnot("checkpoint cell set != current grid" = identical(names(null_list), names(A_cell)))
  message("[regen] resuming: ", done, "/", NBATCH, " batches already done (fingerprint OK)")
} else {
  null_list <- setNames(replicate(length(A_cell),
                 matrix(NA_real_, nrow = NSIM_TOTAL, ncol = NSTAT,
                        dimnames = list(NULL, STAT_NAMES)), simplify = FALSE), names(A_cell))
  k1r <- integer(NM); k2r <- integer(NM); done <- 0L
}

if (done >= NBATCH_RUN) {
  message("[regen] already at or beyond requested batch ", NBATCH_RUN, "; nothing to run")
} else for (b in (done + 1L):NBATCH_RUN) {
  t0   <- Sys.time()
  ef   <- ENV_FILES[b]
  pref <- file.path(ND, sprintf("cRegen_b%02d", b))
  of   <- paste0(pref, "_summary_betai_reg.out")

  ## (re)run BayPass on the preserved null covariates -- unless a complete raw output
  ## for this batch already exists (e.g. a prior invocation crashed AFTER BayPass but
  ## before parsing), in which case reuse it. Completeness is re-asserted at parse.
  reuse <- file.exists(of) && file.info(of)$size > 0
  if (reuse) {
    message(sprintf("  [regen] reusing existing raw output for batch %d (skipping BayPass)", b))
  } else {
    st <- system(paste0(BAYPASS, " -countdatafile ", ND, "/u_eMLG.geno -omegafile ", ND,
      "/omega_mat_omega.out -efile ", ef, " -poolsizefile ", ND, "/u_DIEM.size ",
      OPT, " -outprefix ", pref, " > ", pref, "_stdout.log 2>&1"))
    stopifnot("BayPass batch failed" = st == 0)
  }

  ## parse and validate the 32,840 x 1,000 BF matrix (covariate-major rows).
  ## MRK is the structural guarantee (covariate-major => within each NM-row block the
  ## markers run 1..NM; MRK does not overflow). COVARIABLE is validated only where it
  ## parses: BayPass prints '***' for covariate indices that overflow its fixed-width
  ## integer field (>= 1000), so those rows are checked via MRK + row position instead.
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
  if (any(!okcov)) message(sprintf("  [regen] %d rows had overflowed COVARIABLE ('***', cov idx >= 1000); validated via MRK/row position",
                                   sum(!okcov)))
  M <- matrix(nb$`BF(dB)`, nrow = NM, ncol = BATCH); rm(nb, covnum, okcov, row_cov); invisible(gc())

  ## persist the raw null BF matrix (KEPT) so future re-reductions need no BayPass.
  ## atomic write (tmp then rename) so a crash mid-save cannot leave a partial file.
  bf_out <- file.path(BFDIR, sprintf("cRegen_bf_b%02d.rds", b))
  saveRDS(M, paste0(bf_out, ".tmp")); file.rename(paste0(bf_out, ".tmp"), bf_out)

  ## accumulate per-eMLG exceedance counts (exact cross-check vs Module B; tau-independent)
  k1r <- k1r + rowSums(M >= b1)
  k2r <- k2r + rowSums(M >= b2)

  ## reduce each null covariate to genome-wide statistics, ONCE PER GRID CELL
  ## (identical code path; min = row-subset of M, tau = directional partition)
  rows <- ((b - 1L) * BATCH + 1L):(b * BATCH)
  for (k in names(A_cell)) {
    idx <- cell_idx[[k]]
    null_list[[k]][rows, ] <- t(apply(M[idx, , drop = FALSE], 2, compute_covariate_stats, A = A_cell[[k]]))
    if (!all(is.finite(null_list[[k]][rows, ])))
      stop(sprintf("non-finite statistic in batch %d, cell %s (e.g. empty differentiated top fraction)", b, k))
  }
  rm(M); invisible(gc())

  ## drop the large raw BayPass output; keep only .env + logs + the persisted BF matrix
  file.remove(Sys.glob(paste0(pref, "_summary_*.out")))
  done <- b
  saveRDS(list(null_list = null_list, k1r = k1r, k2r = k2r, done = done,
               fingerprint = fingerprint, group_id = grp), CKPT)
  message(sprintf("[regen] batch %d/%d done in %.1f min (cumulative nulls=%d)",
                  b, NBATCH, as.numeric(difftime(Sys.time(), t0, units = "mins")), b * BATCH))
}

## ---- pilot cross-check (partial counts should track ~ done/NBATCH of saved) ----
frac <- done / NBATCH
cat(sprintf("\n=== reproduction cross-check after %d/%d batches ===\n", done, NBATCH))
cat(sprintf("PC1: sum(partial k1)=%d ; %.3f x saved (%.0f)  [batch subset, ~1.0 expected]\n",
            sum(k1r), sum(k1r) / (frac * sum(k1_saved)), frac * sum(k1_saved)))
cat(sprintf("PC2: sum(partial k2)=%d ; %.3f x saved (%.0f)\n",
            sum(k2r), sum(k2r) / (frac * sum(k2_saved)), frac * sum(k2_saved)))

## ---- finalise ONLY when all batches are done ----------------------------
if (done == NBATCH) {
  ## MANDATORY faithful-reproduction gate -- MONTE-CARLO EQUIVALENCE form.
  ## BayPass is NOT bit-reproducible (a fresh MCMC realization each run; per-BF
  ## differences up to ~25 dB -- moduleC_determinism_probe.R), so exact k equality
  ## is unattainable. Faithfulness is asserted, PER AXIS, as:
  ##   Pearson r(k_regen, k_saved) > COR_THR   AND   |sum(k_regen)/sum(k_saved)-1| < SUM_TOL.
  ## Pearson (not Spearman) so nonlinear discrepancies in the counts cannot hide.
  ## These are pipeline-integrity bounds (the 1000-null pilot agreed to ~0.1-0.2%),
  ## NOT a numerical-reproducibility claim -- if the completed run fails, do NOT
  ## relax them; investigate (inputs are fingerprint-checked) and, if warranted,
  ## recalibrate tolerances with additional independent reruns. Input identity,
  ## ordering, completeness and finiteness remain EXACT hard stops. The primary
  ## genome-wide statistics are separately shown reproducible (~0.0005 MCMC SD,
  ## <=6% of null spread -- moduleC_stat_stability_probe.R).
  COR_THR <- 0.99; SUM_TOL <- 0.03
  cor_k1 <- cor(k1r, k1_saved, method = "pearson")
  cor_k2 <- cor(k2r, k2_saved, method = "pearson")
  rel1 <- sum(k1r) / sum(k1_saved) - 1
  rel2 <- sum(k2r) / sum(k2_saved) - 1
  d1 <- max(abs(k1r - k1_saved)); d2 <- max(abs(k2r - k2_saved))   # reported diagnostic only
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
         "Not writing moduleC_null_stats.rds. Do NOT relax thresholds -- inputs are ",
         "fingerprint-verified, so investigate covariate wiring / BayPass setup; only ",
         "recalibrate tolerances with additional independent reruns if genuinely warranted.")

  ## final structural assertions on EVERY grid cell's null-statistic matrix (EXACT hard stops)
  for (k in names(null_list)) {
    if (nrow(null_list[[k]]) != NSIM_TOTAL)                 stop(sprintf("%s null_stats wrong nrow", k))
    if (!all(is.finite(null_list[[k]])))                   stop(sprintf("%s null_stats has non-finite entries", k))
    if (sum(complete.cases(null_list[[k]])) != NSIM_TOTAL) stop(sprintf("%s not exactly 10,000 complete rows", k))
    if (!all(is.finite(obs_list[[k]])))                    stop(sprintf("%s observed has non-finite entries", k))
    if (!setequal(colnames(obs_list[[k]]), STAT_NAMES))    stop(sprintf("%s observed missing a statistic", k))
  }

  by_cell <- setNames(lapply(names(A_cell), function(k)
    list(observed = obs_list[[k]], null_stats = null_list[[k]])), names(A_cell))

  res <- list(
    ## --- primary (min05_tau06) aliases: keep the flat schema older consumers expect ---
    observed   = obs,                       # 2 x NSTAT (PC1, PC2) at primary cell
    null_stats = null_list[[PRIMARY_CELL]], # 10000 x NSTAT        at primary cell
    ## --- full (min x tau) grid ---
    by_cell     = by_cell,                  # list("min05_tau06", ...) -> {observed, null_stats}
    grid        = CELLS,                    # m, tau, key
    tau_series  = TAUS, tau_primary = MODULEC_TAU_PRIMARY,
    min_series  = MINS, min_primary = MODULEC_MIN_PRIMARY,
    primary_cell = PRIMARY_CELL,
    n_eMLG_by_min = setNames(lengths(idx_by_min), minC_stamp(MINS)),
    stat_names = STAT_NAMES,
    k_check    = list(k1r = k1r, k2r = k2r, k1_saved = k1_saved, k2_saved = k2_saved,
                      cor_k1 = cor_k1, cor_k2 = cor_k2, rel1 = rel1, rel2 = rel2,
                      cor_thr = COR_THR, sum_tol = SUM_TOL,
                      max_abs_dk1 = d1, max_abs_dk2 = d2, reproduced = reproduced),
    params = list(NSIM = NSIM_TOTAL, batch = BATCH, mcmc_seed = MCMC_SEED,
                  config = "aland_excluded / withOmega", top_fracs = TOP_FRACS,
                  primary_stats = PRIMARY_STATS, opt = OPT, tol_bf = TOL_BF,
                  tau_series = TAUS, tau_primary = MODULEC_TAU_PRIMARY,
                  min_series = MINS, min_primary = MODULEC_MIN_PRIMARY,
                  bf_matrices = normalizePath(BFDIR, mustWork = FALSE)),
    fingerprint = fingerprint,
    session = sessionInfo()
  )
  saveRDS(res, OUT)
  if (file.exists(CKPT)) file.remove(CKPT)
  cat(sprintf("\n[regen] DONE. wrote %s (observed + %d null covariates x %d statistics x %d cells: %s; primary %s)\n",
              OUT, NSIM_TOTAL, NSTAT, length(A_cell), paste(names(A_cell), collapse = "/"), PRIMARY_CELL))
  cat(sprintf("[regen] persisted %d null BF matrices in %s (kept; any future (min,tau) re-reducible without BayPass)\n",
              length(Sys.glob(file.path(BFDIR, "cRegen_bf_b*.rds"))), BFDIR))
} else {
  cat(sprintf("\n[regen] PILOT/partial: %d/%d batches done. Checkpoint saved to %s.\n",
              done, NBATCH, CKPT))
  cat("        Re-run with MODC_NBATCH=10 to complete and write the final object.\n")
}
