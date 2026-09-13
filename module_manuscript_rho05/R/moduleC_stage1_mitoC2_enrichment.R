## =========================================================
## module_manuscript_rho05 -- Module C genome-wide enrichment, mitoC2
## =========================================================
## Same DI / recombination / directional-sorting calibration as the
## canonical PC1/PC2/bio_winter Module C, but for the mitotype C2
## population-contrast statistic (mito_C2_S1units_summary_contrast.out),
## NOT a continuous-covariate BF regression -- so it needs its OWN null,
## built from genuine BayPass contrast-mode reruns, not the PC1/PC2 null.
##
## That dedicated contrast null was ALREADY BUILT by
## moduleB_stage1_mitoC2_null.R (2026-09-08/09, on mini2), for the
## Stage-1-direct sim-FDR floor test: 10,000 null contrasts, each a
## population-level +1/-1 partition matching the real 7 Faquilonia-like /
## 12 Fpolyctena-like split, rank-thresholded from the SAME Omega-eigenvector
## draws used for the PC1/PC2 null (see that script's header) -- then
## rerun through BayPass in -contrastfile mode, batched (50 x 200,
## identical MCMC seed 74, identical withOmega Stage-1 setup as the real
## mito_C2 run). Unlike the original PC1/PC2 null run, that script
## PERSISTED every batch's full null log10(1/pval) matrix
## (null/bf_matrices/mitoC2_bf_b##.rds), so recovering the per-unit
## genome-wide reduction needed for Module C required NO new BayPass --
## only pulling those 50 already-computed matrices back from mini2
## (2026-09-13) and reducing them here exactly like
## moduleC_stage1_null_regen.R reduces the climate null's persisted BF.
##
## log10(1/pval) (not C2) is used as the per-unit "association strength"
## fed into compute_covariate_stats() -- same choice moduleB_stage1_mitoC2_
## null.R made for its own exceedance test, and monotonic with C2 within a
## batch, so rank-based statistics (the primary family here) are identical
## either way.
##
## INTEGRITY CHECK: because these are the EXACT persisted matrices already
## used by moduleB_stage1_mitoC2_null.R (not a fresh MCMC realization), the
## per-unit exceedance count recomputed here must equal that script's saved
## k3 EXACTLY (not just Monte-Carlo-equivalent) -- hard stop otherwise.
##
## Reads : module_manuscript_rho05/data/moduleC_stage1_annotations.rds
##         module_manuscript_rho05/data/moduleB_stage1_mitoC2_null.rds  (k3 cross-check)
##         module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/
##           S1units_group_order.txt, mito_C2_S1units_summary_contrast.out,
##           u.mito_contrast, null/bf_matrices/mitoC2_bf_b01..50.rds
## Writes: module_manuscript_rho05/data/moduleC_stage1_mitoC2_null_stats.rds
##         module_manuscript_rho05/data/moduleC_stage1_mitoC2_results.rds
##         module_manuscript_rho05/data/moduleC_stage1_mitoC2_primary_tests.tsv
##         module_manuscript_rho05/data/moduleC_stage1_mitoC2_tau_sensitivity.tsv
##         module_manuscript_rho05/Figures/moduleC_stage1_mitoC2_null_calibration.{png,pdf}
##         module_manuscript_rho05/Figures/moduleC_stage1_mitoC2_threshold_sensitivity.{png,pdf}
##         module_manuscript_rho05/doc/moduleC_stage1_mitoC2_report.md
##
## Run from the repo root:
##   Rscript module_manuscript_rho05/R/moduleC_stage1_mitoC2_enrichment.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2) })
source("moduleC_climate_vs_sorting/R/moduleC_stat_functions.R")

BASE <- "module_manuscript_rho05"
BP_DIR <- file.path(BASE, "baypass_stage1", "aland_excluded_S1units")
ND     <- file.path(BP_DIR, "null"); BFDIR <- file.path(ND, "bf_matrices")
DATA <- file.path(BASE, "data"); FIG <- file.path(BASE, "Figures"); DOC <- file.path(BASE, "doc")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(DOC, showWarnings = FALSE, recursive = TRUE)

AX <- "mitoC2"
NSIM_TOTAL <- 10000L; BATCH <- 200L; NBATCH <- NSIM_TOTAL / BATCH

## ---- inputs ------------------------------------------------------------
ann <- readRDS(file.path(DATA, "moduleC_stage1_annotations.rds"))
grp <- readLines(file.path(BP_DIR, "S1units_group_order.txt"))
stopifnot("annotation order != BayPass order" = identical(ann$group_id, grp))
NM <- nrow(ann)

obs_c2 <- fread(file.path(BP_DIR, "mito_C2_S1units_summary_contrast.out"))
stopifnot("mito_C2 row count != Stage-1 unit count" = nrow(obs_c2) == NM,
          "mito_C2 MRK not 1..N in order" = all(obs_c2$MRK == seq_len(NM)),
          "non-finite mito_C2 log10p" = all(is.finite(obs_c2$`log10(1/pval)`)))
obs_log10p <- obs_c2$`log10(1/pval)`

mbnull <- readRDS(file.path(DATA, "moduleB_stage1_mitoC2_null.rds"))
stopifnot("moduleB_stage1_mitoC2_null group set != unit universe" = setequal(mbnull$group_id, grp))
mo <- match(grp, mbnull$group_id)
k3_saved <- mbnull$k3[mo]
stopifnot(all(is.finite(k3_saved)))

BF_FILES <- file.path(BFDIR, sprintf("mitoC2_bf_b%02d.rds", seq_len(NBATCH)))
stopifnot("some persisted mitoC2 null log10p matrices are missing -- pull them back from mini2 (null/bf_matrices/mitoC2_bf_b##.rds)" =
            all(file.exists(BF_FILES)))

## ---- cell grid (mirrors moduleC_stage1_null_regen.R: single min level) --
TAUS <- MODULEC_TAU_SERIES; MINS <- 5L; TSTAMP <- tauC_stamp(TAUS)
PRIMARY_CELL <- cell_key(MINS, MODULEC_TAU_PRIMARY)
stopifnot("all units must satisfy the single min threshold" = all(ann$n_loci >= MINS))
CELLS <- as.data.table(expand.grid(m = MINS, tau = TAUS))[, key := cell_key(m, tau)]
A_cell <- setNames(lapply(CELLS$tau, function(tau) prepare_annotation_ranks(ann, dir_col = dir_col_for_tau(tau))), CELLS$key)

STAT_NAMES <- covariate_stat_names(); NSTAT <- length(STAT_NAMES)
PDIR_COLS   <- grep("_pdir_diff$", STAT_NAMES, value = TRUE)
STRICT_COLS <- setdiff(STAT_NAMES, PDIR_COLS)

obs_by_cell <- setNames(lapply(names(A_cell), function(k)
  compute_covariate_stats(obs_log10p, A_cell[[k]])), names(A_cell))

## ---- reduce the 50 persisted null log10p matrices (NO BayPass) ---------
null_list <- setNames(replicate(length(A_cell),
               matrix(NA_real_, nrow = NSIM_TOTAL, ncol = NSTAT, dimnames = list(NULL, STAT_NAMES)),
               simplify = FALSE), names(A_cell))
k3_check <- integer(NM)

for (b in seq_len(NBATCH)) {
  t0 <- Sys.time()
  Mp <- readRDS(BF_FILES[b])
  stopifnot("persisted mitoC2 null matrix wrong shape" = nrow(Mp) == NM && ncol(Mp) == BATCH,
            "non-finite log10p in persisted null matrix" = all(is.finite(Mp)))
  k3_check <- k3_check + rowSums(Mp >= obs_log10p)
  rows <- ((b - 1L) * BATCH + 1L):(b * BATCH)
  for (k in names(A_cell)) {
    Mk <- t(apply(Mp, 2, compute_covariate_stats, A = A_cell[[k]]))
    null_list[[k]][rows, ] <- Mk
    if (!all(is.finite(Mk[, STRICT_COLS])))
      stop(sprintf("non-finite statistic outside *_pdir_diff in batch %d, cell %s", b, k))
  }
  rm(Mp); invisible(gc())
  message(sprintf("[mito-enrich] batch %d/%d reduced in %.1f s", b, NBATCH,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

## exact integrity check (deterministic reduction of the SAME persisted
## matrices already used by moduleB_stage1_mitoC2_null.R -- must match
## exactly, not just Monte-Carlo-equivalent)
stopifnot("k3 recomputed from the persisted null matrices != the saved floor-test k3 -- data mismatch, do not proceed" =
            all(k3_check == k3_saved))
message("[mito-enrich] exact integrity check PASSED: recomputed exceedance counts == moduleB_stage1_mitoC2_null.rds$k3 for all units")

for (k in names(null_list)) {
  if (nrow(null_list[[k]]) != NSIM_TOTAL) stop(sprintf("%s null_stats wrong nrow", k))
  if (sum(complete.cases(null_list[[k]][, STRICT_COLS])) != NSIM_TOTAL)
    stop(sprintf("%s not exactly %d complete rows (excluding *_pdir_diff)", k, NSIM_TOTAL))
}

null <- null_list[[PRIMARY_CELL]]
obs  <- obs_by_cell[[PRIMARY_CELL]]
NSIM <- nrow(null)

nres_mito <- list(
  observed = obs, null_stats = null, by_cell = setNames(lapply(names(A_cell), function(k)
    list(observed = obs_by_cell[[k]], null_stats = null_list[[k]])), names(A_cell)),
  tau_series = TAUS, tau_primary = MODULEC_TAU_PRIMARY, min_series = MINS,
  primary_cell = PRIMARY_CELL, n_units_by_min = setNames(NM, minC_stamp(MINS)),
  k3_integrity = list(matches_saved_exactly = TRUE, n_units = NM),
  null_source = "moduleB_stage1_mitoC2_null.R persisted batches (null/bf_matrices/mitoC2_bf_b##.rds), reduced here; no new BayPass",
  built = as.character(Sys.time()))
saveRDS(nres_mito, file.path(DATA, "moduleC_stage1_mitoC2_null_stats.rds"))

## ---- empirical two-sided P -------------------------------------------
emp_p_two_sided <- function(t_obs, t_null) {
  t_null <- t_null[!is.na(t_null)]
  med <- median(t_null)
  (1 + sum(abs(t_null - med) >= abs(t_obs - med))) / (length(t_null) + 1)
}

LAB <- c(rho_DI = "DI (Spearman rho)", rho_rec = "recombination (Spearman rho)",
         sort_gap_differentiated = "sorting, differentiated only (log10p percentile gap)",
         sort_gap_all = "sorting, all units (log10p percentile gap)",
         rho_sort_magnitude = "sorting magnitude / prop_fixed (Spearman rho)",
         rho_sort_orientation = "sorting orientation / uni_score (Spearman rho)",
         pear_DI = "DI (raw log10(1/pval) Pearson)", pear_rec = "recombination (raw log10(1/pval) Pearson)",
         pear_sort_magnitude = "sorting magnitude (raw log10(1/pval) Pearson)",
         bf_gap_differentiated = "sorting, differentiated (raw log10(1/pval) gap)",
         bf_gap_all = "sorting, all units (raw log10(1/pval) gap)")

PRIMARY <- c("rho_DI", "rho_rec", "sort_gap_differentiated")
SUPP    <- c("sort_gap_all", "rho_sort_magnitude", "rho_sort_orientation",
             "pear_DI", "pear_rec", "pear_sort_magnitude",
             "bf_gap_differentiated", "bf_gap_all")

build_row <- function(stat, family) {
  tn <- null[, stat]; to <- obs[[stat]]
  lab <- if (stat %in% names(LAB)) LAB[[stat]] else stat
  data.table(family = family, test = lab, stat = stat, axis = AX,
             observed = to, null_median = median(tn, na.rm = TRUE),
             null_lo = unname(quantile(tn, 0.025, na.rm = TRUE)), null_hi = unname(quantile(tn, 0.975, na.rm = TRUE)),
             p_emp = emp_p_two_sided(to, tn))
}
tab <- rbindlist(c(
  lapply(PRIMARY, build_row, family = "primary"),
  lapply(SUPP,    build_row, family = "supplementary")
))
tab[family == "primary", p_adj := p.adjust(p_emp, method = "BH")]

## ---- tau sensitivity ----------------------------------------------------
GRID_STATS <- c("rho_DI", "rho_rec", "sort_gap_differentiated")
grid_rows <- rbindlist(lapply(names(obs_by_cell), function(k) {
  pc <- CELLS[key == k]
  nn <- null_list[[k]]; oo <- obs_by_cell[[k]]
  rbindlist(lapply(GRID_STATS, function(s) {
    tn <- nn[, s]; to <- oo[[s]]
    data.table(cell = k, tau = pc$tau,
               test = if (s %in% names(LAB)) LAB[[s]] else s, stat = s, axis = AX,
               observed = to, null_median = median(tn, na.rm = TRUE),
               null_lo = unname(quantile(tn, 0.025, na.rm = TRUE)), null_hi = unname(quantile(tn, 0.975, na.rm = TRUE)),
               p_emp = emp_p_two_sided(to, tn))
  }))
}))
chk <- grid_rows[cell == PRIMARY_CELL & stat %in% GRID_STATS][order(stat)]
ref <- tab[stat %in% GRID_STATS][order(stat)]
stopifnot("grid primary-cell rows disagree with flat primary alias" = isTRUE(all.equal(chk$observed, ref$observed)))

## ---- threshold sensitivity table ----------------------------------------
th_stats <- grep("^top", STAT_NAMES, value = TRUE)
th_tab <- rbindlist(lapply(th_stats, build_row, family = "threshold"))
frac_lookup <- c(top0001 = "top 0.1%", top0005 = "top 0.5%", top0010 = "top 1%")
metric_lookup <- c(DI_med = "median DI", rec_med = "median recomb",
                   pdiff = "prop. differentiated", pdir_diff = "prop. directional | diff.")
th_tab[, threshold_tag := sub("_.*", "", stat)]
th_tab[, metric_tag    := sub("^top[0-9]+_", "", stat)]
th_tab[, frac   := unname(frac_lookup[threshold_tag])]
th_tab[, metric := unname(metric_lookup[metric_tag])]
stopifnot(all(!is.na(th_tab$frac)), all(!is.na(th_tab$metric)))

## ---- results object -----------------------------------------------------
results <- list(
  axis = AX,
  primary = tab[family == "primary"], supplementary = tab[family == "supplementary"],
  threshold = th_tab, tau_sensitivity = grid_rows, observed = obs, null_stats = null,
  k3_integrity = nres_mito$k3_integrity, mito_contrast_meta = attr(mbnull, "meta"),
  meta = list(NSIM = NSIM, p_formula = "two-sided vs null median; (1+#{|Tn-med|>=|To-med|})/(NSIM+1)",
              fdr = "BH across the 3 primary tests (single axis)",
              tau_series = TAUS, tau_primary = MODULEC_TAU_PRIMARY,
              primary_cell = PRIMARY_CELL, n_units_by_min = nres_mito$n_units_by_min,
              null_source = nres_mito$null_source,
              built = as.character(Sys.time())))
saveRDS(results, file.path(DATA, "moduleC_stage1_mitoC2_results.rds"))

fwrite(grid_rows[, .(cell, tau, test, axis, observed = round(observed, 4),
        null_median = round(null_median, 4), null_lo = round(null_lo, 4),
        null_hi = round(null_hi, 4), p_emp = signif(p_emp, 3))],
       file.path(DATA, "moduleC_stage1_mitoC2_tau_sensitivity.tsv"), sep = "\t")
fwrite(tab[family == "primary", .(test, axis, observed = round(observed, 4),
        null_median = round(null_median, 4), null_lo = round(null_lo, 4),
        null_hi = round(null_hi, 4), p_emp = signif(p_emp, 3), p_adj = signif(p_adj, 3))],
       file.path(DATA, "moduleC_stage1_mitoC2_primary_tests.tsv"), sep = "\t")

cat("\n=== PRIMARY TESTS (Stage-1-direct structured-null calibration, mitoC2) ===\n")
print(tab[family == "primary", .(test, axis, observed = round(observed, 4),
      null_median = round(null_median, 4),
      null95 = sprintf("[%.3f, %.3f]", null_lo, null_hi),
      p_emp = signif(p_emp, 3), p_adj = signif(p_adj, 3))])

## ========================================================================
## FIGURE 1 -- structured-null calibration (primary tests + magnitude)
## ========================================================================
fig_stats <- c(PRIMARY, "rho_sort_magnitude")
long <- rbindlist(lapply(fig_stats, function(s)
  data.table(label = LAB[[s]], value = null[, s])))
long[, label := factor(label, levels = LAB[fig_stats])]
mk <- rbindlist(lapply(fig_stats, function(s) data.table(
  label = LAB[[s]], median = median(null[, s]),
  lo = quantile(null[, s], 0.025), hi = quantile(null[, s], 0.975),
  obs = obs[[s]])))
mk[, label := factor(label, levels = LAB[fig_stats])]

p1 <- ggplot(long, aes(value)) +
  geom_histogram(bins = 60, fill = "grey80", colour = NA) +
  geom_rect(data = mk, inherit.aes = FALSE, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            fill = "#999999", alpha = 0.18) +
  geom_vline(data = mk, aes(xintercept = median), linetype = 2, colour = "grey40") +
  geom_vline(data = mk, aes(xintercept = obs), colour = "#CC79A7", linewidth = 0.8) +
  facet_wrap(~ label, scales = "free", ncol = 2) +
  labs(x = "genome-wide statistic", y = "count (of 10,000 null contrasts)",
       title = "Module C (Stage-1-direct): observed mitoC2 vs 10,000 structure-matched null population contrasts",
       subtitle = "shaded = null 95% interval; dashed = null median; solid pink = observed mitoC2 (own dedicated contrast null, NOT the PC1/PC2/bio_winter null)") +
  theme_classic(base_size = 9) +
  theme(strip.text = element_text(size = 8), strip.background = element_blank())
ggsave(file.path(FIG, "moduleC_stage1_mitoC2_null_calibration.png"), p1,
       width = 175, height = 140, units = "mm", dpi = 300)
ggsave(file.path(FIG, "moduleC_stage1_mitoC2_null_calibration.pdf"), p1,
       width = 175, height = 140, units = "mm")

## ========================================================================
## FIGURE 2 -- rank-threshold sensitivity (top 0.1% / 0.5% / 1%)
## ========================================================================
lev  <- c("top 0.1%", "top 0.5%", "top 1%")
mlev <- unname(metric_lookup)
setfac <- function(d) { d[, frac := factor(frac, levels = lev)]
                        d[, metric := factor(metric, levels = mlev)]; d }
th_meta <- unique(th_tab[, .(stat, frac, metric)])
thl <- rbindlist(lapply(th_stats, function(s) data.table(stat = s, value = null[, s])))
thl <- setfac(th_meta[thl, on = "stat"])
thm <- setfac(th_tab[, .(stat, frac, metric, median = null_median, lo = null_lo, hi = null_hi, obs = observed)])

p2 <- ggplot(thl, aes(value)) +
  geom_histogram(bins = 45, fill = "grey80", colour = NA) +
  geom_rect(data = thm, inherit.aes = FALSE, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            fill = "#999999", alpha = 0.18) +
  geom_vline(data = thm, aes(xintercept = median), linetype = 2, colour = "grey40") +
  geom_vline(data = thm, aes(xintercept = obs), colour = "#CC79A7", linewidth = 0.6) +
  facet_wrap(vars(metric, frac), scales = "free", ncol = 3,
             labeller = label_wrap_gen(multi_line = FALSE)) +
  labs(x = "value within the top-ranked fraction of Stage-1 units", y = "count (null contrasts)",
       title = "Module C (Stage-1-direct) threshold sensitivity, mitoC2: top rank fractions, observed vs null",
       subtitle = "prop. directional is computed among differentiated units in the fraction") +
  theme_classic(base_size = 8) +
  theme(strip.text = element_text(size = 7), strip.background = element_blank())
ggsave(file.path(FIG, "moduleC_stage1_mitoC2_threshold_sensitivity.png"), p2,
       width = 185, height = 165, units = "mm", dpi = 300)
ggsave(file.path(FIG, "moduleC_stage1_mitoC2_threshold_sensitivity.pdf"), p2,
       width = 185, height = 165, units = "mm")

## ========================================================================
## data-driven Markdown report
## ========================================================================
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
mm <- attr(mbnull, "meta")

rep_lines <- c(
"# Module C, Stage-1-direct: genome-wide, unit-level mitoC2 calibration",
"",
sprintf("*Generated %s. NSIM = %d structure-matched null population contrasts; Stage-1-direct unit universe = %d clusters (n_snps>=5, Aland excluded, 19 pops, Stage-1-derived Omega).*",
        as.character(Sys.Date()), NSIM, NM),
"",
"## Data provenance",
"",
"- **Scope:** the Stage-1 LARGE-CLUSTER scan (18,361 clusters with >=5 SNPs), identical universe to the PC1/PC2/bio_winter Module C runs.",
sprintf("- **Observed mitoC2 association:** per-unit BayPass log10(1/pval) from the C2 population-contrast test (`mito_C2_S1units_summary_contrast.out`), contrasting %d Faquilonia-like vs %d Fpolyctena-like populations (19 populations total, Aland excluded) -- same withOmega Stage-1-direct BayPass setup as PC1/PC2/bio_winter (`-omegafile` pointed at the identical frozen Stage-1 Omega), but in `-contrastfile` mode rather than `-efile` covariate-regression mode.",
        mm$n_pos, mm$n_neg),
"- **Null contrasts: a DEDICATED contrast-mode null, NOT the PC1/PC2/bio_winter null.** A continuous-covariate null (random Omega-structured numeric draws through `-efile`) is not a matched reference for a binary population-contrast statistic, so this reuses the null already built by `moduleB_stage1_mitoC2_null.R` for the Stage-1 sim-FDR floor test: each of 10,000 null contrasts is a population-level +1/-1 partition with the SAME group sizes as the real split (7 vs 12), obtained by rank-thresholding the same Omega-eigenvector draws used for the climate null (so populations with high covariance still tend to land on the same side of the null partition, preserving the structure-matching property), then rerun through BayPass in `-contrastfile` mode with the identical MCMC seed (74) and Stage-1 Omega as the real mitoC2 run.",
"- **No new BayPass run was needed for THIS Module C reduction.** Unlike the PC1/PC2 null (which only kept exceedance counts, forcing a ~9-10h BayPass rerun for `moduleC_stage1_null_regen.R`), `moduleB_stage1_mitoC2_null.R` persisted every batch's full null log10(1/pval) matrix (`null/bf_matrices/mitoC2_bf_b##.rds`, 50 files). Those were pulled back from mini2 and reduced here directly.",
sprintf("- **Exact integrity check (not Monte-Carlo tolerance):** because these are the literal persisted matrices already used for the floor test (not a fresh MCMC realization), the per-unit null-exceedance count recomputed here matches `moduleB_stage1_mitoC2_null.rds`'s saved `k3` EXACTLY for all %d units (hard stop otherwise, enforced in the script).", NM),
"- **Annotations (per-unit, joined by `group_id`):** identical to the PC1/PC2/bio_winter runs -- consensus best-SNP Diagnostic Index, unidirectional-sorting status, sorting magnitude `prop_fixed` and signed orientation `uni_score`, recombination rate (cM/Mb) at the best-SNP marker, cluster size.",
"",
"## Validation checks",
"",
sprintf("- All 50 persisted null matrices present and correctly shaped (%d x 200), all finite.", NM),
"- Exact k3 exceedance-count match against the floor-test run (see above).",
sprintf("- mitoC2 log10(1/pval) vector: N = %d, order verified identical to the Stage-1 BayPass row order, all finite.", NM),
"",
"## Methods",
"",
"Identical statistical framework to the PC1/PC2/bio_winter Module C analyses, with log10(1/pval) substituted for BF(dB) as the per-unit association-strength input (monotonic with C2 within each null draw, so the rank-based primary statistics are unaffected by this substitution): the genome-wide mitoC2 log10(1/pval) vector is reduced to Spearman rho with DI, Spearman rho with recombination, and the directional-sorting percentile-gap among differentiated units (primary), plus supplementary statistics. The observed reduction is compared against the dedicated 10,000-null-contrast distribution by a two-sided empirical P; the three primary tests (mitoC2 x DI/recombination/sorting) are BH-FDR corrected.",
"",
"## Results",
"",
"| test | axis | observed | null median | null 95% | p_emp | p_adj |",
"|---|---|---|---|---|---|---|")
prim_md <- tab[family == "primary"][order(stat)]
rep_lines <- c(rep_lines,
  prim_md[, sprintf("| %s | %s | %s | %s | [%s, %s] | %s | %s |",
    test, axis, fmt(observed), fmt(null_median), fmt(null_lo), fmt(null_hi),
    signif(p_emp, 3), signif(p_adj, 3))])
rep_lines <- c(rep_lines, "",
"### Supplementary statistics (not in the FDR family)",
"",
"| test | axis | observed | null median | null 95% | p_emp |",
"|---|---|---|---|---|---|")
supp_md <- tab[family == "supplementary"][order(stat)]
rep_lines <- c(rep_lines,
  supp_md[, sprintf("| %s | %s | %s | %s | [%s, %s] | %s |",
    test, axis, fmt(observed), fmt(null_median), fmt(null_lo), fmt(null_hi), signif(p_emp, 3))])

ts_series <- paste(sprintf("%.1f", TAUS), collapse = ", ")
rep_lines <- c(rep_lines, "",
  "### Sensitivity to the fixation threshold (tau)",
  "",
  sprintf("Reported over the fixation-threshold tau in {%s} (Stage-1-direct universe fixed at n_snps>=5). DI and recombination do not depend on tau (shown once, at the primary tau); directional sorting is recomputed at each tau. Empirical P only (the FDR family is the three primary tests at the primary tau).",
          ts_series),
  "",
  "**Directional sorting (differentiated-only) across tau:**",
  "",
  "| tau | axis | observed | null 95% | p_emp |",
  "|---|---|---|---|---|")
gsd <- grid_rows[stat == "sort_gap_differentiated"][order(tau)]
rep_lines <- c(rep_lines,
  gsd[, sprintf("| %.1f | %s | %s | [%s, %s] | %s |",
    tau, axis, fmt(observed), fmt(null_lo), fmt(null_hi), signif(p_emp, 3))])

## ---- data-driven interpretation (NOT presupposing the outcome) -----------
gp <- function(st) tab[family == "primary" & stat == st]
sig_primary <- tab[family == "primary" & p_adj < 0.05][order(p_adj)]
s1 <- gp("sort_gap_differentiated"); r1 <- gp("rho_rec"); d1 <- gp("rho_DI")
pd1 <- tab[family == "supplementary" & stat == "pear_DI"]
mag1 <- tab[stat == "rho_sort_magnitude"]
sideword <- function(row) if (row$observed > row$null_hi) "above" else if (row$observed < row$null_lo) "below" else "within"

describe_axis <- function(row, label, extra = "") {
  sig <- row$p_adj < 0.05
  if (sig) {
    sprintf("**%s: a mitoC2 association survives FDR.** Observed %.3f (FDR %.3f, %s the null)%s.",
            label, row$observed, row$p_adj, sideword(row), extra)
  } else {
    sprintf("**%s: no mitoC2 association survives FDR.** Observed %.3f (FDR %.3f), within the structured-null 95%% interval.",
            label, row$observed, row$p_adj)
  }
}

sort_sentence <- describe_axis(s1, "Directional sorting (primary, differentiated-only)")
if (s1$p_adj >= 0.05)
  sort_sentence <- paste0(sort_sentence, sprintf(" Sorting magnitude (`prop_fixed`) is likewise null (supplementary, p_emp %.2f).", mag1$p_emp))
rec_sentence <- describe_axis(r1, "Recombination")
di_sentence  <- describe_axis(d1, "Diagnostic Index",
  extra = sprintf(", corroborated by the raw log10(1/pval) analysis (Pearson %.3f, p %.3g)", pd1$observed, pd1$p_emp))
di_sentence <- paste0(di_sentence, " (DI is a signed index; the sign is reported as-is and should not be read as locus-level adaptation in diagnostic regions.)")
overall_sentence <- if (nrow(sig_primary) > 0) {
  sprintf("**Overall:** of the three primary tests, %d survives FDR (%s).",
          nrow(sig_primary), paste(sig_primary$test, collapse = "; "))
} else {
  "**Overall:** no primary test is exceptional; mitoC2 association evidence is not concentrated in diagnostic, directionally-sorted, or low-recombination Stage-1 units beyond population structure and genomic architecture."
}

interp <- c(sort_sentence, rec_sentence, di_sentence, overall_sentence)

rep_lines <- c(rep_lines, "",
"## Interpretation",
"",
interp,
"",
"## What this analysis can and cannot establish",
"",
"- **Can:** calibrate genome-wide, Stage-1-unit-level mitoC2-association *patterns* against a null built specifically for the contrast-mode statistic (matched group sizes, same Omega-eigenvector structure-preservation as the continuous-covariate null), with the same unit universe and annotation set as the PC1/PC2/bio_winter runs.",
"- **Cannot:** identify individual mitoC2-associated loci beyond what the Stage-1-direct sim-FDR floor test already flags (`moduleB_stage1_mitoC2_null.rds`); be compared directly, panel-for-panel, against the PC1/PC2/bio_winter calibration figure on ONE set of null histograms -- the null distributions are built from genuinely different processes (continuous Omega-eigenvector draws vs rank-thresholded +1/-1 partitions of those same draws) and are reported on separate figures for that reason, even though the primary statistics themselves (Spearman rho, percentile gap) are on a comparable scale.",
"")
writeLines(rep_lines, file.path(DOC, "moduleC_stage1_mitoC2_report.md"))

cat(sprintf("\n[mito-enrich] wrote null_stats, results, TSVs, 2 figures, and doc/moduleC_stage1_mitoC2_report.md\n"))
cat(sprintf("[mito-enrich] dedicated contrast null built from persisted moduleB_stage1_mitoC2_null.R batches; exact k3 integrity check PASSED\n"))
