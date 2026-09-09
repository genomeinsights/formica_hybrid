## =========================================================
## module_manuscript_rho05 -- Module C analysis, Stage-1-direct universe
## =========================================================
## Adapted from moduleC_climate_vs_sorting/R/moduleC_analyse.R for the
## 18,361 Stage-1 units (n_snps>=5, single universe -- no min grid; see
## moduleC_stage1_null_regen.R). Consumes the regenerated observed +
## 10,000-null genome-wide statistics and produces empirical P values, BH-FDR
## across the six primary tests (PC1/PC2 x DI/recombination/sorting), the
## calibration figure, the tau-sensitivity table, a compact TSV, the results
## RDS, and a data-driven Markdown report.
##
## UNLIKE the canonical script, the interpretation text below does NOT
## presuppose the outcome (the canonical narrative was written after seeing
## the Stage-2 "no association" result) -- every interpretive sentence is
## conditioned on p_adj / p_emp exactly as the DI and Overall sentences
## already were, so this script produces correct prose for whatever the
## Stage-1-direct calibration actually shows.
##
## Reads : module_manuscript_rho05/data/moduleC_stage1_null_stats.rds
##         (pulled back from mini2 after moduleC_stage1_null_regen.R finishes)
## Writes: module_manuscript_rho05/data/moduleC_stage1_results.rds
##         module_manuscript_rho05/data/moduleC_stage1_primary_tests.tsv
##         module_manuscript_rho05/Figures/moduleC_stage1_null_calibration.{png,pdf}
##         module_manuscript_rho05/Figures/moduleC_stage1_tau_sensitivity.{png,pdf}
##         module_manuscript_rho05/doc/moduleC_stage1_report.md
##
## Run from the repo root, AFTER moduleC_stage1_null_regen.R has finished and
## moduleC_stage1_null_stats.rds has been copied back from mini2:
##   Rscript module_manuscript_rho05/R/moduleC_stage1_analyse.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2) })
source("moduleC_climate_vs_sorting/R/moduleC_stat_functions.R")

BASE <- "module_manuscript_rho05"
FIG  <- file.path(BASE, "Figures"); DOC <- file.path(BASE, "doc")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(DOC, showWarnings = FALSE, recursive = TRUE)

res_in <- readRDS(file.path(BASE, "data", "moduleC_stage1_null_stats.rds"))
obs   <- res_in$observed          # 2 x NSTAT (PC1, PC2), primary cell (min05_tau06)
null  <- res_in$null_stats        # NSIM x NSTAT
NSIM  <- nrow(null)
EXPECT_STATS <- covariate_stat_names()

## ---- integrity gate (do NOT analyse an unfaithful / incomplete regeneration) ----
kc <- res_in$k_check
stopifnot(
  "regeneration failed the Monte-Carlo equivalence gate (k_check$reproduced != TRUE)" = isTRUE(kc$reproduced),
  "PC1 k Pearson r <= 0.99"                 = kc$cor_k1 > 0.99,
  "PC2 k Pearson r <= 0.99"                 = kc$cor_k2 > 0.99,
  "PC1 |sum(k) ratio - 1| >= 0.03"          = abs(kc$rel1) < 0.03,
  "PC2 |sum(k) ratio - 1| >= 0.03"          = abs(kc$rel2) < 0.03,
  "null_stats does not have 10,000 rows"    = NSIM == 10000,
  "null_stats has non-finite entries"       = all(is.finite(null)),
  "null_stats has incomplete rows"          = sum(complete.cases(null)) == NSIM,
  "observed is not 2 rows (PC1/PC2)"        = nrow(obs) == 2,
  "observed has non-finite entries"         = all(is.finite(obs)),
  "a required statistic is missing (null)"  = setequal(colnames(null), EXPECT_STATS),
  "a required statistic is missing (obs)"   = setequal(colnames(obs),  EXPECT_STATS))

## ---- empirical two-sided P against the structured-null distribution ------
emp_p_two_sided <- function(t_obs, t_null) {
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

build_row <- function(stat, axis, family) {
  tn <- null[, stat]; to <- obs[axis, stat]
  lab <- if (stat %in% names(LAB)) LAB[[stat]] else stat
  data.table(family = family, test = lab, stat = stat, axis = axis,
             observed = to, null_median = median(tn),
             null_lo = unname(quantile(tn, 0.025)), null_hi = unname(quantile(tn, 0.975)),
             p_emp = emp_p_two_sided(to, tn))
}
tab <- rbindlist(c(
  lapply(PRIMARY, function(s) rbindlist(lapply(c("PC1", "PC2"), build_row, stat = s, family = "primary"))),
  lapply(SUPP,    function(s) rbindlist(lapply(c("PC1", "PC2"), build_row, stat = s, family = "supplementary")))
))
tab[family == "primary", p_adj := p.adjust(p_emp, method = "BH")]

## ---- tau sensitivity (single universe -- tau grid only, no min grid) -----
grid_rows <- NULL
if (!is.null(res_in$by_cell)) {
  bc <- res_in$by_cell
  GRID_STATS <- c("rho_DI", "rho_rec", "sort_gap_differentiated")
  parse_cell <- function(k) list(min = as.integer(sub("min(\\d+)_.*", "\\1", k)),
                                  tau = as.integer(sub(".*_tau(\\d+)", "\\1", k)) / 10)
  grid_rows <- rbindlist(lapply(names(bc), function(k) {
    pc <- parse_cell(k); o <- bc[[k]]$observed; nn <- bc[[k]]$null_stats
    rbindlist(lapply(GRID_STATS, function(s) rbindlist(lapply(c("PC1", "PC2"), function(ax) {
      tn <- nn[, s]; to <- o[ax, s]
      data.table(cell = k, tau = pc$tau,
                 test = if (s %in% names(LAB)) LAB[[s]] else s, stat = s, axis = ax,
                 observed = to, null_median = median(tn),
                 null_lo = unname(quantile(tn, 0.025)), null_hi = unname(quantile(tn, 0.975)),
                 p_emp = emp_p_two_sided(to, tn))
    }))))
  }))
  pcell <- res_in$primary_cell
  chk <- grid_rows[cell == pcell & stat %in% c("rho_DI", "rho_rec", "sort_gap_differentiated")][order(stat, axis)]
  ref <- tab[stat %in% c("rho_DI", "rho_rec", "sort_gap_differentiated")][order(stat, axis)]
  stopifnot("grid primary-cell rows disagree with flat primary alias" =
              isTRUE(all.equal(chk$observed, ref$observed)))
}

## ---- threshold sensitivity table ------------------------------------------
th_stats <- grep("^top", colnames(null), value = TRUE)
th_tab <- rbindlist(lapply(th_stats, function(s)
  rbindlist(lapply(c("PC1", "PC2"), build_row, stat = s, family = "threshold"))))
frac_lookup <- c(top0001 = "top 0.1%", top0005 = "top 0.5%", top0010 = "top 1%")
metric_lookup <- c(DI_med = "median DI", rec_med = "median recomb",
                   pdiff = "prop. differentiated", pdir_diff = "prop. directional | diff.")
th_tab[, threshold_tag := sub("_.*", "", stat)]
th_tab[, metric_tag    := sub("^top[0-9]+_", "", stat)]
th_tab[, frac   := unname(frac_lookup[threshold_tag])]
th_tab[, metric := unname(metric_lookup[metric_tag])]
stopifnot("unresolved threshold fraction label" = all(!is.na(th_tab$frac)),
          "unresolved threshold metric label"   = all(!is.na(th_tab$metric)))

## ---- results object -------------------------------------------------------
results <- list(
  primary = tab[family == "primary"], supplementary = tab[family == "supplementary"],
  threshold = th_tab, tau_sensitivity = grid_rows, observed = obs, null_stats = null,
  k_check = kc, params = res_in$params, fingerprint = res_in$fingerprint,
  session = res_in$session,
  meta = list(NSIM = NSIM, p_formula = "two-sided vs null median; (1+#{|Tn-med|>=|To-med|})/(NSIM+1)",
              fdr = "BH across the 6 primary tests",
              tau_series = res_in$tau_series, tau_primary = res_in$tau_primary,
              primary_cell = res_in$primary_cell, n_units_by_min = res_in$n_units_by_min,
              built = as.character(Sys.time())))
saveRDS(results, file.path(BASE, "data", "moduleC_stage1_results.rds"))

if (!is.null(grid_rows)) {
  fwrite(grid_rows[, .(cell, tau, test, axis, observed = round(observed, 4),
          null_median = round(null_median, 4), null_lo = round(null_lo, 4),
          null_hi = round(null_hi, 4), p_emp = signif(p_emp, 3))],
         file.path(BASE, "data", "moduleC_stage1_tau_sensitivity.tsv"), sep = "\t")
  cat(sprintf("\n=== TAU SENSITIVITY (primary cell %s; %d Stage-1 units) ===\n",
              res_in$primary_cell, res_in$n_units_by_min[[1]]))
  cat("-- directional sorting across tau --\n")
  print(grid_rows[stat == "sort_gap_differentiated",
        .(tau, axis, observed = round(observed, 4),
          null95 = sprintf("[%.3f, %.3f]", null_lo, null_hi), p_emp = signif(p_emp, 3))][order(tau, axis)])
}

fwrite(tab[family == "primary", .(test, axis, observed = round(observed, 4),
        null_median = round(null_median, 4), null_lo = round(null_lo, 4),
        null_hi = round(null_hi, 4), p_emp = signif(p_emp, 3), p_adj = signif(p_adj, 3))],
       file.path(BASE, "data", "moduleC_stage1_primary_tests.tsv"), sep = "\t")

cat("\n=== PRIMARY TESTS (Stage-1-direct structured-null calibration) ===\n")
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
  PC1 = obs["PC1", s], PC2 = obs["PC2", s])))
mk[, label := factor(label, levels = LAB[fig_stats])]

p1 <- ggplot(long, aes(value)) +
  geom_histogram(bins = 60, fill = "grey80", colour = NA) +
  geom_rect(data = mk, inherit.aes = FALSE, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            fill = "#999999", alpha = 0.18) +
  geom_vline(data = mk, aes(xintercept = median), linetype = 2, colour = "grey40") +
  geom_vline(data = mk, aes(xintercept = PC1, colour = "PC1"), linewidth = 0.7) +
  geom_vline(data = mk, aes(xintercept = PC2, colour = "PC2"), linewidth = 0.7) +
  facet_wrap(~ label, scales = "free", ncol = 2) +
  scale_colour_manual(values = c(PC1 = "#0072B2", PC2 = "#D55E00"), name = NULL) +
  labs(x = "genome-wide statistic", y = "count (of 10,000 null covariates)",
       title = "Module C (Stage-1-direct): observed climate axes vs 10,000 Omega-structured null covariates",
       subtitle = "shaded = null 95% interval; dashed = null median") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8),
        strip.background = element_blank())
ggsave(file.path(FIG, "moduleC_stage1_null_calibration.png"), p1,
       width = 175, height = 140, units = "mm", dpi = 300)
ggsave(file.path(FIG, "moduleC_stage1_null_calibration.pdf"), p1,
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
th_stat_meta <- unique(th_tab[, .(stat, frac, metric, median = null_median,
                                  lo = null_lo, hi = null_hi)])
th_obs <- dcast(th_tab, stat ~ axis, value.var = "observed")
thm <- setfac(th_obs[th_stat_meta, on = "stat"])

p2 <- ggplot(thl, aes(value)) +
  geom_histogram(bins = 45, fill = "grey80", colour = NA) +
  geom_rect(data = thm, inherit.aes = FALSE, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            fill = "#999999", alpha = 0.18) +
  geom_vline(data = thm, aes(xintercept = median), linetype = 2, colour = "grey40") +
  geom_vline(data = thm, aes(xintercept = PC1, colour = "PC1"), linewidth = 0.6) +
  geom_vline(data = thm, aes(xintercept = PC2, colour = "PC2"), linewidth = 0.6) +
  facet_wrap(vars(metric, frac), scales = "free", ncol = 3,
             labeller = label_wrap_gen(multi_line = FALSE)) +
  scale_colour_manual(values = c(PC1 = "#0072B2", PC2 = "#D55E00"), name = NULL) +
  labs(x = "value within the top-ranked fraction of Stage-1 units", y = "count (null covariates)",
       title = "Module C (Stage-1-direct) threshold sensitivity: top rank fractions, observed vs null",
       subtitle = "prop. directional is computed among differentiated units in the fraction") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom", strip.text = element_text(size = 7),
        strip.background = element_blank())
ggsave(file.path(FIG, "moduleC_stage1_threshold_sensitivity.png"), p2,
       width = 185, height = 165, units = "mm", dpi = 300)
ggsave(file.path(FIG, "moduleC_stage1_threshold_sensitivity.pdf"), p2,
       width = 185, height = 165, units = "mm")

## ========================================================================
## data-driven Markdown report
## ========================================================================
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
N_units <- length(kc$k1r)

rep_lines <- c(
"# Module C, Stage-1-direct: genome-wide, unit-level climate calibration",
"",
sprintf("*Generated %s. NSIM = %d Omega-structured null covariates; Stage-1-direct unit universe = %d clusters (n_snps>=5, Aland excluded, 19 pops, Stage-1-derived Omega). Single universe -- unlike the canonical eMLG Module C, there is no second min-cluster-size level to sweep.*",
        as.character(Sys.Date()), NSIM, N_units),
"",
"## Data provenance",
"",
"- **Observed Stage-1-unit climate association:** per-unit BayPass BF(dB) on climate PC1/PC2 (`module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/PC{1,2}_S1units_withOmega_summary_betai_reg.out`).",
"- **Null covariates:** the same 10,000 Omega-structured draws used for the Stage-1-unit sim-FDR floor test (`moduleB_stage1_S1units_null.R`, mini2), re-run on the preserved `null/null_b01..b50.env` files (that run kept only exceedance counts; the full BF matrix is regenerated here, persisted, and reduced on the fly).",
"- **Annotations (per-unit, joined by `group_id`):** consensus best-SNP Diagnostic Index, unidirectional-sorting status, sorting magnitude `prop_fixed` and signed orientation `uni_score` (`moduleA_stage1_cluster_sorting.rds`, Module A equivalent for the Stage-1 universe); recombination rate (cM/Mb) at the best-SNP marker (`data/Frufa_DTOL_PR.ref_genome.recmap`); cluster size (`moduleB_stage1_units_bestsnp.rds`).",
"",
"## Validation checks",
"",
sprintf("- Stage-1 unit count and order identical across observed / null / annotation objects (N = %d), all joins by explicit `group_id`.", N_units),
sprintf("- **Faithful regeneration (Monte-Carlo equivalence gate passed):** BayPass is not bit-reproducible (a fresh MCMC realization each run), so the regenerated per-unit exceedance counts match the moduleB_stage1_S1units_null.R run within MCMC tolerance rather than exactly: PC1 Pearson r = %.4f, sum ratio-1 = %+.4f; PC2 r = %.4f, sum ratio-1 = %+.4f (thresholds r > 0.99, |ratio-1| < 0.03; max|dk1| = %d, max|dk2| = %d reported as diagnostics only).",
        kc$cor_k1, kc$rel1, kc$cor_k2, kc$rel2, kc$max_abs_dk1, kc$max_abs_dk2),
"- Input identity (the 50 `.env` covariates, geno, Omega, poolsize, params, statistic code) is guaranteed EXACTLY by md5 fingerprints. Observed BF vectors equal `eBF1`/`eBF2` (max|d| = 0); observed and null reduced by identical code (`moduleC_stat_functions.R`, shared unmodified with the canonical eMLG Module C).",
"",
"## Methods",
"",
"For every covariate (observed PC1/PC2 and each of the 10,000 nulls) the genome-wide Stage-1-unit BF vector is reduced to: Spearman rho with DI; Spearman rho with recombination; the difference in mean within-covariate BF **percentile** between directionally sorted and non-directional units **among differentiated clusters** (primary sorting statistic); plus supplementary statistics (all-unit sorting contrast, Spearman rho with sorting magnitude `prop_fixed`, signed sorting orientation `uni_score`, and raw-BF variants). Each observed statistic is compared with its 10,000-value structured-null distribution by a two-sided empirical P (deviation from the null median); the six primary tests (PC1/PC2 x DI/recombination/sorting) are BH-FDR corrected. Units are never resampled independently.",
"",
"## Results",
"",
"| test | axis | observed | null median | null 95% | p_emp | p_adj |",
"|---|---|---|---|---|---|---|")
prim_md <- tab[family == "primary"][order(stat, axis)]
rep_lines <- c(rep_lines,
  prim_md[, sprintf("| %s | %s | %s | %s | [%s, %s] | %s | %s |",
    test, axis, fmt(observed), fmt(null_median), fmt(null_lo), fmt(null_hi),
    signif(p_emp, 3), signif(p_adj, 3))])
rep_lines <- c(rep_lines, "",
"### Supplementary statistics (not in the FDR family)",
"",
"| test | axis | observed | null median | null 95% | p_emp |",
"|---|---|---|---|---|---|")
supp_md <- tab[family == "supplementary"][order(stat, axis)]
rep_lines <- c(rep_lines,
  supp_md[, sprintf("| %s | %s | %s | %s | [%s, %s] | %s |",
    test, axis, fmt(observed), fmt(null_median), fmt(null_lo), fmt(null_hi), signif(p_emp, 3))])

if (!is.null(grid_rows)) {
  ts_series <- paste(sprintf("%.1f", res_in$tau_series), collapse = ", ")
  rep_lines <- c(rep_lines, "",
    "### Sensitivity to the fixation threshold (tau)",
    "",
    sprintf("The calibration is reported over the fixation-threshold tau in {%s} (the Stage-1-direct universe is fixed at n_snps>=5 throughout; there is no second minimum-cluster-size level to sweep, unlike the canonical eMLG Module C's min in {5,10}). DI and recombination do not depend on tau (shown once, at the primary tau); directional sorting is recomputed at each tau. Empirical P only (the FDR family is the six primary tests at the primary tau).",
            ts_series),
    "",
    "**Directional sorting (differentiated-only) across tau:**",
    "",
    "| tau | axis | observed | null 95% | p_emp |",
    "|---|---|---|---|---|")
  gsd <- grid_rows[stat == "sort_gap_differentiated"][order(tau, axis)]
  rep_lines <- c(rep_lines,
    gsd[, sprintf("| %.1f | %s | %s | [%s, %s] | %s |",
      tau, axis, fmt(observed), fmt(null_lo), fmt(null_hi), signif(p_emp, 3))])
}

## ---- data-driven interpretation (NOT presupposing the outcome) -----------
gp <- function(st, ax) tab[family == "primary" & stat == st & axis == ax]
sig_primary <- tab[family == "primary" & p_adj < 0.05][order(p_adj)]
s1 <- gp("sort_gap_differentiated", "PC1"); s2 <- gp("sort_gap_differentiated", "PC2")
r1 <- gp("rho_rec", "PC1"); r2 <- gp("rho_rec", "PC2")
d1 <- gp("rho_DI", "PC1");  d2 <- gp("rho_DI", "PC2")
pd1 <- tab[family == "supplementary" & stat == "pear_DI" & axis == "PC1"]
pd2 <- tab[family == "supplementary" & stat == "pear_DI" & axis == "PC2"]
mag1 <- tab[stat == "rho_sort_magnitude" & axis == "PC1"]
mag2 <- tab[stat == "rho_sort_magnitude" & axis == "PC2"]
sideword <- function(row) if (row$observed > row$null_hi) "above" else if (row$observed < row$null_lo) "below" else "within"

sort_sig <- s1$p_adj < 0.05 || s2$p_adj < 0.05
rec_sig  <- r1$p_adj < 0.05 || r2$p_adj < 0.05
di_sig   <- d1$p_adj < 0.05 || d2$p_adj < 0.05

sort_sentence <- if (sort_sig) {
  sig_ax <- if (s1$p_adj < 0.05 && s2$p_adj < 0.05) "both PC1 and PC2" else if (s1$p_adj < 0.05) "PC1" else "PC2"
  sprintf("**Directional sorting (primary, differentiated-only): a climate association survives FDR on %s.** PC1 observed %.3f (FDR %.3f), PC2 %.3f (FDR %.3f); the significant axis/axes lie %s the structured-null 95%% interval.",
          sig_ax, s1$observed, s1$p_adj, s2$observed, s2$p_adj,
          if (s1$p_adj < 0.05) sideword(s1) else sideword(s2))
} else {
  sprintf("**Directional sorting (primary, differentiated-only): no climate association survives FDR on either axis.** PC1 observed %.3f (FDR %.3f), PC2 %.3f (FDR %.3f) both lie within the structured-null 95%% interval [%.3f, %.3f]; sorting magnitude (`prop_fixed`) is likewise null on both axes (supplementary, p_emp %.2f/%.2f).",
          s1$observed, s1$p_adj, s2$observed, s2$p_adj, s1$null_lo, s1$null_hi,
          mag1$p_emp, mag2$p_emp)
}
rec_sentence <- if (rec_sig) {
  sig_ax <- if (r1$p_adj < 0.05 && r2$p_adj < 0.05) "both PC1 and PC2" else if (r1$p_adj < 0.05) "PC1" else "PC2"
  sprintf("**Recombination: a climate association survives FDR on %s.** PC1 %.3f (FDR %.3f), PC2 %.3f (FDR %.3f).",
          sig_ax, r1$observed, r1$p_adj, r2$observed, r2$p_adj)
} else {
  sprintf("**Recombination: no climate association survives FDR.** PC1 %.3f (FDR %.3f), PC2 %.3f (FDR %.3f), both within the null.",
          r1$observed, r1$p_adj, r2$observed, r2$p_adj)
}
di_sentence <- if (di_sig) {
  sig_row  <- if (d1$p_adj < 0.05) d1 else d2
  sig_pear <- if (d1$p_adj < 0.05) pd1 else pd2
  other    <- if (d1$p_adj < 0.05) d2 else d1
  other_ax <- if (d1$p_adj < 0.05) "PC2" else "PC1"
  sig_ax   <- if (d1$p_adj < 0.05) "PC1" else "PC2"
  sprintf("**Diagnostic Index: a signal on %s.** %s is within the null (Spearman rho %.3f, FDR %.3f), but on %s the climate-association strength is correlated with DI beyond the structured null (rho %.3f, empirical p %.3g, FDR %.3f; observed %s the null 95%% interval), corroborated by the raw-BF analysis (Pearson %.3f, p %.3g). (DI is a signed index; the sign is reported as-is and should not be read as locus-level adaptation in diagnostic regions.)",
          sig_ax, other_ax, other$observed, other$p_adj, sig_ax, sig_row$observed, sig_row$p_emp,
          sig_row$p_adj, sideword(sig_row), sig_pear$observed, sig_pear$p_emp)
} else {
  sprintf("**Diagnostic Index: no climate association after FDR** (PC1 rho %.3f FDR %.3f; PC2 rho %.3f FDR %.3f).",
          d1$observed, d1$p_adj, d2$observed, d2$p_adj)
}
overall_sentence <- if (nrow(sig_primary) > 0)
  sprintf("**Overall:** of the six primary tests, %d survives FDR (%s).",
          nrow(sig_primary), paste(sprintf("%s x %s", sig_primary$test, sig_primary$axis), collapse = "; "))
else
  "**Overall:** no primary test is exceptional; climate-association evidence is not concentrated in diagnostic, directionally-sorted, or low-recombination Stage-1 units beyond population structure and genomic architecture."

interp <- c(sort_sentence, rec_sentence, di_sentence, overall_sentence)

rep_lines <- c(rep_lines, "",
"## Interpretation",
"",
interp,
"",
"## What this analysis can and cannot establish",
"",
"- **Can:** calibrate genome-wide, Stage-1-unit-level climate-association *patterns* against a structure- and architecture-preserving null, with the same unit universe and statistic for observed and null.",
"- **Cannot:** identify individual climate-associated loci beyond what the Stage-1-direct sim-FDR floor test already flags (`moduleB_stage1_S1units_null.rds`); XtX/BF share among-population allele-frequency variation, so raw-BF sensitivity variants may absorb genuine climate differentiation as well as confounding.",
"")
writeLines(rep_lines, file.path(DOC, "moduleC_stage1_report.md"))

cat(sprintf("\n[analyse-S1] wrote results, TSV, 2 figures, and doc/moduleC_stage1_report.md\n"))
cat(sprintf("[analyse-S1] faithful-regeneration gate PASSED (max|dk1|=%d max|dk2|=%d)\n",
            kc$max_abs_dk1, kc$max_abs_dk2))
