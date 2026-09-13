## =========================================================
## module_manuscript_rho05 -- Module C genome-wide enrichment, bio_winter
## =========================================================
## Same DI / recombination / directional-sorting calibration as the
## canonical PC1/PC2 Module C (moduleC_stage1_analyse.R), applied to
## bio_winter, the Stage-1-direct (18,361 units, n_snps>=5, Aland excluded)
## BayPass covariate run. bio_winter's per-unit BF(dB) is the SAME
## continuous-covariate-regression statistic as PC1/PC2 (verified: same
## file format, same BayPass row order, same withOmega Stage-1-direct
## setup) -- so it is compared against the EXISTING 10,000-null-covariate
## reference distribution in moduleC_stage1_null_stats.rds. That null
## distribution is covariate-agnostic (it is built purely from random
## Omega-structured covariates' own BF-vs-annotation relationship, never
## from PC1/PC2 values themselves), so no new BayPass/mini2 run is needed --
## this script only regenerates the OBSERVED reduction for bio_winter and
## reuses the null already on disk. An ann-vs-null fingerprint check below
## guards against silently reusing a null built from stale annotations.
##
## mitoC2 is NOT included here (per-2026-09-13 decision): it is a BayPass
## C2 population-contrast statistic, not a continuous-covariate regression,
## so the existing null is not a matched reference for it -- left out
## pending a decision on a dedicated contrast null.
##
## Reads : module_manuscript_rho05/data/moduleC_stage1_annotations.rds
##         module_manuscript_rho05/data/moduleC_stage1_null_stats.rds
##         module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/
##           S1units_group_order.txt, bio_winter_S1units_withOmega_summary_betai_reg.out
## Writes: module_manuscript_rho05/data/moduleC_stage1_bio_winter_results.rds
##         module_manuscript_rho05/data/moduleC_stage1_bio_winter_primary_tests.tsv
##         module_manuscript_rho05/Figures/moduleC_stage1_bio_winter_null_calibration.{png,pdf}
##         module_manuscript_rho05/Figures/moduleC_stage1_bio_winter_threshold_sensitivity.{png,pdf}
##         module_manuscript_rho05/doc/moduleC_stage1_bio_winter_report.md
##
## Run from the repo root:
##   Rscript module_manuscript_rho05/R/moduleC_stage1_bio_winter_enrichment.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(digest) })
source("moduleC_climate_vs_sorting/R/moduleC_stat_functions.R")

BASE <- "module_manuscript_rho05"
BP_DIR <- file.path(BASE, "baypass_stage1", "aland_excluded_S1units")
DATA <- file.path(BASE, "data"); FIG <- file.path(BASE, "Figures"); DOC <- file.path(BASE, "doc")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(DOC, showWarnings = FALSE, recursive = TRUE)

AX <- "bio_winter"

## ---- inputs ----------------------------------------------------------------
ann <- readRDS(file.path(DATA, "moduleC_stage1_annotations.rds"))
grp <- readLines(file.path(BP_DIR, "S1units_group_order.txt"))
stopifnot("annotation order != BayPass order" = identical(ann$group_id, grp))

nres <- readRDS(file.path(DATA, "moduleC_stage1_null_stats.rds"))
TAUS <- nres$tau_series; TSTAMP <- tauC_stamp(TAUS)
PRIMARY_CELL <- nres$primary_cell

bw <- fread(file.path(BP_DIR, "bio_winter_S1units_withOmega_summary_betai_reg.out"))
NM <- nrow(ann)
stopifnot("bio_winter row count != Stage-1 unit count" = nrow(bw) == NM,
          "bio_winter MRK not 1..N in order" = all(bw$MRK == seq_len(NM)),
          "non-finite bio_winter BF" = all(is.finite(bw$`BF(dB)`)))
b_bw <- bw$`BF(dB)`

## ---- fingerprint check: the annotations underlying the existing null must
## exactly match the annotations used here (guards against silently pairing
## bio_winter's observed stats with a null built from a different/stale
## annotation set) ------------------------------------------------------------
FP_ANN_COLS <- c("group_id", "DI", "recomb", "prop_fixed", "uni_score",
                 paste0("directional_", TSTAMP), "differentiated", "n_loci")
ann_fp <- digest(ann[, ..FP_ANN_COLS], algo = "md5")
stopifnot("annotation fingerprint != the fingerprint recorded when the null was built; do not reuse this null" =
            identical(ann_fp, nres$fingerprint$ann))
message("[bw-enrich] annotation fingerprint matches the null's recorded fingerprint -- null reuse is valid")

## ---- observed bio_winter statistics, per (min,tau) cell -------------------
CELLS <- as.data.table(expand.grid(m = nres$min_series, tau = TAUS))[, key := cell_key(m, tau)]
obs_by_cell <- setNames(vector("list", nrow(CELLS)), CELLS$key)
for (i in seq_len(nrow(CELLS))) {
  A <- prepare_annotation_ranks(ann, dir_col = dir_col_for_tau(CELLS$tau[i]))
  obs_by_cell[[CELLS$key[i]]] <- compute_covariate_stats(b_bw, A)
}
stopifnot("primary cell missing from grid" = PRIMARY_CELL %in% names(obs_by_cell),
          "primary cell missing from existing null (by_cell)" = PRIMARY_CELL %in% names(nres$by_cell))
obs <- obs_by_cell[[PRIMARY_CELL]]
null <- nres$by_cell[[PRIMARY_CELL]]$null_stats
NSIM <- nrow(null)

EXPECT_STATS <- covariate_stat_names()
stopifnot("observed stat set != expected" = setequal(names(obs), EXPECT_STATS),
          "null stat set != expected"     = setequal(colnames(null), EXPECT_STATS),
          "observed has non-finite entries" = all(is.finite(obs)))
PDIR_COLS   <- grep("_pdir_diff$", EXPECT_STATS, value = TRUE)
STRICT_COLS <- setdiff(EXPECT_STATS, PDIR_COLS)
stopifnot("null_stats does not have 10,000 rows" = NSIM == 10000,
          "null_stats has non-finite entries outside *_pdir_diff" = all(is.finite(null[, STRICT_COLS])))

## ---- empirical two-sided P (same formula as the PC1/PC2 analysis) --------
emp_p_two_sided <- function(t_obs, t_null) {
  t_null <- t_null[!is.na(t_null)]
  med <- median(t_null)
  (1 + sum(abs(t_null - med) >= abs(t_obs - med))) / (length(t_null) + 1)
}

LAB <- c(rho_DI = "DI (Spearman rho)", rho_rec = "recombination (Spearman rho)",
         sort_gap_differentiated = "sorting, differentiated only (BF percentile gap)",
         sort_gap_all = "sorting, all units (BF percentile gap)",
         rho_sort_magnitude = "sorting magnitude / prop_fixed (Spearman rho)",
         rho_sort_orientation = "sorting orientation / uni_score (Spearman rho)",
         pear_DI = "DI (raw-BF Pearson)", pear_rec = "recombination (raw-BF Pearson)",
         pear_sort_magnitude = "sorting magnitude (raw-BF Pearson)",
         bf_gap_differentiated = "sorting, differentiated (raw-BF gap)",
         bf_gap_all = "sorting, all units (raw-BF gap)")

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

## ---- tau sensitivity (single axis) -----------------------------------------
GRID_STATS <- c("rho_DI", "rho_rec", "sort_gap_differentiated")
grid_rows <- rbindlist(lapply(names(obs_by_cell), function(k) {
  pc <- CELLS[key == k]
  nn <- nres$by_cell[[k]]$null_stats; oo <- obs_by_cell[[k]]
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
stopifnot("grid primary-cell rows disagree with flat primary alias" =
            isTRUE(all.equal(chk$observed, ref$observed)))

## ---- threshold sensitivity table ------------------------------------------
th_stats <- grep("^top", EXPECT_STATS, value = TRUE)
th_tab <- rbindlist(lapply(th_stats, build_row, family = "threshold"))
frac_lookup <- c(top0001 = "top 0.1%", top0005 = "top 0.5%", top0010 = "top 1%")
metric_lookup <- c(DI_med = "median DI", rec_med = "median recomb",
                   pdiff = "prop. differentiated", pdir_diff = "prop. directional | diff.")
th_tab[, threshold_tag := sub("_.*", "", stat)]
th_tab[, metric_tag    := sub("^top[0-9]+_", "", stat)]
th_tab[, frac   := unname(frac_lookup[threshold_tag])]
th_tab[, metric := unname(metric_lookup[metric_tag])]
stopifnot("unresolved threshold fraction label" = all(!is.na(th_tab$frac)),
          "unresolved threshold metric label"   = all(!is.na(th_tab$metric)))

## ---- results object ---------------------------------------------------------
results <- list(
  axis = AX,
  primary = tab[family == "primary"], supplementary = tab[family == "supplementary"],
  threshold = th_tab, tau_sensitivity = grid_rows, observed = obs, null_stats = null,
  ann_fingerprint_check = list(ann_fp = ann_fp, null_ann_fp = nres$fingerprint$ann, match = TRUE),
  params = nres$params, fingerprint = nres$fingerprint, session = sessionInfo(),
  meta = list(NSIM = NSIM, p_formula = "two-sided vs null median; (1+#{|Tn-med|>=|To-med|})/(NSIM+1)",
              fdr = "BH across the 3 primary tests (single axis)",
              tau_series = TAUS, tau_primary = nres$tau_primary,
              primary_cell = PRIMARY_CELL, n_units_by_min = nres$n_units_by_min,
              null_source = "reused from moduleC_stage1_null_stats.rds (PC1/PC2 run); no new BayPass",
              built = as.character(Sys.time())))
saveRDS(results, file.path(DATA, "moduleC_stage1_bio_winter_results.rds"))

fwrite(grid_rows[, .(cell, tau, test, axis, observed = round(observed, 4),
        null_median = round(null_median, 4), null_lo = round(null_lo, 4),
        null_hi = round(null_hi, 4), p_emp = signif(p_emp, 3))],
       file.path(DATA, "moduleC_stage1_bio_winter_tau_sensitivity.tsv"), sep = "\t")
fwrite(tab[family == "primary", .(test, axis, observed = round(observed, 4),
        null_median = round(null_median, 4), null_lo = round(null_lo, 4),
        null_hi = round(null_hi, 4), p_emp = signif(p_emp, 3), p_adj = signif(p_adj, 3))],
       file.path(DATA, "moduleC_stage1_bio_winter_primary_tests.tsv"), sep = "\t")

cat("\n=== PRIMARY TESTS (Stage-1-direct structured-null calibration, bio_winter) ===\n")
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
  geom_vline(data = mk, aes(xintercept = obs), colour = "#009E73", linewidth = 0.8) +
  facet_wrap(~ label, scales = "free", ncol = 2) +
  labs(x = "genome-wide statistic", y = "count (of 10,000 null covariates)",
       title = "Module C (Stage-1-direct): observed bio_winter vs 10,000 Omega-structured null covariates",
       subtitle = "shaded = null 95% interval; dashed = null median; solid green = observed bio_winter") +
  theme_classic(base_size = 9) +
  theme(strip.text = element_text(size = 8), strip.background = element_blank())
ggsave(file.path(FIG, "moduleC_stage1_bio_winter_null_calibration.png"), p1,
       width = 175, height = 140, units = "mm", dpi = 300)
ggsave(file.path(FIG, "moduleC_stage1_bio_winter_null_calibration.pdf"), p1,
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
  geom_vline(data = thm, aes(xintercept = obs), colour = "#009E73", linewidth = 0.6) +
  facet_wrap(vars(metric, frac), scales = "free", ncol = 3,
             labeller = label_wrap_gen(multi_line = FALSE)) +
  labs(x = "value within the top-ranked fraction of Stage-1 units", y = "count (null covariates)",
       title = "Module C (Stage-1-direct) threshold sensitivity, bio_winter: top rank fractions, observed vs null",
       subtitle = "prop. directional is computed among differentiated units in the fraction") +
  theme_classic(base_size = 8) +
  theme(strip.text = element_text(size = 7), strip.background = element_blank())
ggsave(file.path(FIG, "moduleC_stage1_bio_winter_threshold_sensitivity.png"), p2,
       width = 185, height = 165, units = "mm", dpi = 300)
ggsave(file.path(FIG, "moduleC_stage1_bio_winter_threshold_sensitivity.pdf"), p2,
       width = 185, height = 165, units = "mm")

## ========================================================================
## data-driven Markdown report
## ========================================================================
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
N_units <- nres$n_units_by_min[[1]]

rep_lines <- c(
"# Module C, Stage-1-direct: genome-wide, unit-level bio_winter calibration",
"",
sprintf("*Generated %s. NSIM = %d Omega-structured null covariates (REUSED from the PC1/PC2 Module C run -- see Data provenance); Stage-1-direct unit universe = %d clusters (n_snps>=5, Aland excluded, 19 pops, Stage-1-derived Omega).*",
        as.character(Sys.Date()), NSIM, N_units),
"",
"## Data provenance",
"",
"- **Scope:** the Stage-1 LARGE-CLUSTER scan (18,361 clusters with >=5 SNPs), identical universe to the canonical PC1/PC2 Module C run (`moduleC_stage1_report.md`).",
"- **Observed bio_winter association:** per-unit BayPass BF(dB) on the bio_winter covariate (`module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/bio_winter_S1units_withOmega_summary_betai_reg.out`) -- same withOmega Stage-1-direct BayPass setup as PC1/PC2 (identical `.geno`/`omega_mat_omega.out`/poolsize inputs, same continuous-covariate regression model).",
"- **Null covariates: REUSED, not regenerated.** The 10,000 Omega-structured null covariates and their reduced genome-wide statistics come from `moduleC_stage1_null_stats.rds` (the PC1/PC2 Module C run). That null distribution is a property of random-covariate-BF-vs-annotation relationships under this Omega/genotype/unit setup -- it does not depend on which real covariate (PC1, PC2, or bio_winter) is being tested, so no new BayPass/mini2 run was needed here. Validity of this reuse is enforced by an exact md5 fingerprint match between this run's annotations and the annotations recorded when the null was built (checked below, hard stop on mismatch).",
"- **Annotations (per-unit, joined by `group_id`):** identical to the PC1/PC2 run -- consensus best-SNP Diagnostic Index, unidirectional-sorting status, sorting magnitude `prop_fixed` and signed orientation `uni_score`, recombination rate (cM/Mb) at the best-SNP marker, cluster size.",
"",
"## Validation checks",
"",
sprintf("- Annotation fingerprint used here == fingerprint recorded when the reused null was built (md5 `%s`); the stopifnot in this script hard-stops if they ever diverge.", ann_fp),
sprintf("- bio_winter BF vector: N = %d, order verified identical to the Stage-1 BayPass row order (`S1units_group_order.txt`), all finite.", NM),
"",
"## Methods",
"",
"Identical to the PC1/PC2 Module C analysis (see `moduleC_stage1_report.md`): the genome-wide bio_winter BF vector is reduced to Spearman rho with DI, Spearman rho with recombination, and the directional-sorting percentile gap among differentiated units (primary), plus supplementary statistics. The observed reduction is compared against the (reused) 10,000-null distribution by a two-sided empirical P; the three primary tests (bio_winter x DI/recombination/sorting) are BH-FDR corrected.",
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
    sprintf("**%s: a bio_winter association survives FDR.** Observed %.3f (FDR %.3f, %s the null)%s.",
            label, row$observed, row$p_adj, sideword(row), extra)
  } else {
    sprintf("**%s: no bio_winter association survives FDR.** Observed %.3f (FDR %.3f), within the structured-null 95%% interval.",
            label, row$observed, row$p_adj)
  }
}

sort_sentence <- describe_axis(s1, "Directional sorting (primary, differentiated-only)")
if (s1$p_adj >= 0.05)
  sort_sentence <- paste0(sort_sentence, sprintf(" Sorting magnitude (`prop_fixed`) is likewise null (supplementary, p_emp %.2f).", mag1$p_emp))
rec_sentence <- describe_axis(r1, "Recombination")
di_sentence  <- describe_axis(d1, "Diagnostic Index",
  extra = sprintf(", corroborated by the raw-BF analysis (Pearson %.3f, p %.3g)", pd1$observed, pd1$p_emp))
di_sentence <- paste0(di_sentence, " (DI is a signed index; the sign is reported as-is and should not be read as locus-level adaptation in diagnostic regions.)")
overall_sentence <- if (nrow(sig_primary) > 0) {
  sprintf("**Overall:** of the three primary tests, %d survives FDR (%s).",
          nrow(sig_primary), paste(sig_primary$test, collapse = "; "))
} else {
  "**Overall:** no primary test is exceptional; bio_winter association evidence is not concentrated in diagnostic, directionally-sorted, or low-recombination Stage-1 units beyond population structure and genomic architecture."
}

interp <- c(sort_sentence, rec_sentence, di_sentence, overall_sentence)

rep_lines <- c(rep_lines, "",
"## Interpretation",
"",
interp,
"",
"## What this analysis can and cannot establish",
"",
"- **Can:** calibrate genome-wide, Stage-1-unit-level bio_winter-association *patterns* against a structure- and architecture-preserving null, with the same unit universe and statistic as the PC1/PC2 run, at no extra BayPass cost.",
"- **Cannot:** identify individual bio_winter-associated loci beyond what the Stage-1-direct sim-FDR floor test already flags (`moduleB_stage1_bio_winter_null.rds`); XtX/BF share among-population allele-frequency variation, so raw-BF sensitivity variants may absorb genuine bio_winter differentiation as well as confounding.",
"")
writeLines(rep_lines, file.path(DOC, "moduleC_stage1_bio_winter_report.md"))

cat(sprintf("\n[bw-enrich] wrote results, TSVs, 2 figures, and doc/moduleC_stage1_bio_winter_report.md\n"))
cat(sprintf("[bw-enrich] null reused from moduleC_stage1_null_stats.rds (fingerprint-verified); no new BayPass run\n"))
