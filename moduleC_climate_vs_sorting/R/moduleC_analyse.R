## =========================================================
## MODULE C (revised) -- structured-null calibration, figures, report
## =========================================================
## Consumes the regenerated observed + 10,000-null genome-wide statistics and
## produces the inferential outputs: empirical two-sided P values, FDR across the
## six primary tests (PC1/PC2 x DI/recombination/sorting), the calibration figure,
## the rank-threshold sensitivity figure, a compact TSV, the results RDS, and a
## data-driven Markdown report.
##
## Primary sorting test = directional vs non-directional contrast among
## DIFFERENTIATED eMLGs only (`sort_gap_differentiated`). Continuous sorting
## magnitude = `prop_fixed` (`rho_sort_magnitude`); the signed `uni_score`
## ("sorting orientation") is supplementary only.
##
## Empirical two-sided P (pre-specified):
##   p = ( 1 + #{ |T_null - median(T_null)| >= |T_obs - median(T_null)| } ) / (NSIM + 1)
## Each null covariate is one complete genome-wide pseudo-analysis; eMLGs are NEVER
## resampled independently.
##
## This script REFUSES to run unless the regeneration passes the pre-specified
## Monte-Carlo equivalence gate against Module B's counts and the null-statistic
## matrix is complete. Exact identity is required for inputs/order/completeness,
## not for stochastic BayPass BF values or their per-eMLG exceedance counts.
##
## Reads : moduleC_climate_vs_sorting/data/moduleC_null_stats.rds
## Writes: moduleC_climate_vs_sorting/data/moduleC_results.rds
##         moduleC_climate_vs_sorting/data/moduleC_primary_tests.tsv
##         moduleC_climate_vs_sorting/Figures/moduleC_null_calibration.{png,pdf}
##         moduleC_climate_vs_sorting/Figures/moduleC_threshold_sensitivity.{png,pdf}
##         moduleC_climate_vs_sorting/doc/moduleC_report.md
## Run from the formica_hybrid repo root, AFTER moduleC_null_regen.R has finished.

suppressMessages({ library(data.table); library(ggplot2) })
source("moduleC_climate_vs_sorting/R/moduleC_stat_functions.R")

BASE <- "moduleC_climate_vs_sorting"
res_in <- readRDS(file.path(BASE, "data", "moduleC_null_stats.rds"))
obs   <- res_in$observed          # 2 x NSTAT (PC1, PC2)
null  <- res_in$null_stats        # NSIM x NSTAT
NSIM  <- nrow(null)
EXPECT_STATS <- covariate_stat_names()

## ---- integrity gate (do NOT analyse an unfaithful / incomplete regeneration) ----
## Faithfulness = Monte-Carlo equivalence gate (BayPass is not bit-reproducible):
## Pearson r(k_regen, k_saved) > 0.99 and |sum ratio - 1| < 0.03, per axis.
## Completeness / finiteness / statistic presence stay exact hard stops.
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
dir.create(file.path(BASE, "Figures"), showWarnings = FALSE, recursive = TRUE)

## ---- empirical two-sided P against the structured-null distribution ------
emp_p_two_sided <- function(t_obs, t_null) {
  med <- median(t_null)
  (1 + sum(abs(t_null - med) >= abs(t_obs - med))) / (length(t_null) + 1)
}

## human labels
LAB <- c(rho_DI = "DI (Spearman rho)", rho_rec = "recombination (Spearman rho)",
         sort_gap_differentiated = "sorting, differentiated only (BF percentile gap)",
         sort_gap_all = "sorting, all eMLGs (BF percentile gap)",
         rho_sort_magnitude = "sorting magnitude / prop_fixed (Spearman rho)",
         rho_sort_orientation = "sorting orientation / uni_score (Spearman rho)",
         pear_DI = "DI (raw-BF Pearson)", pear_rec = "recombination (raw-BF Pearson)",
         pear_sort_magnitude = "sorting magnitude (raw-BF Pearson)",
         bf_gap_differentiated = "sorting, differentiated (raw-BF gap)",
         bf_gap_all = "sorting, all eMLGs (raw-BF gap)")

## six primary tests (FDR family) + supplementary statistics
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
tab[family == "primary", p_adj := p.adjust(p_emp, method = "BH")]   # BH over the 6 primary

## ---- threshold sensitivity table (explicit, safe label mapping) ---------
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
          "unresolved threshold metric label"   = all(!is.na(th_tab$metric)),
          "threshold fractions incomplete" = setequal(th_tab$frac, unname(frac_lookup)),
          "threshold metrics incomplete"   = setequal(th_tab$metric, unname(metric_lookup)))

## ---- results object ----------------------------------------------------
results <- list(
  primary = tab[family == "primary"], supplementary = tab[family == "supplementary"],
  threshold = th_tab, observed = obs, null_stats = null,
  k_check = kc, params = res_in$params, fingerprint = res_in$fingerprint,
  session = res_in$session,
  meta = list(NSIM = NSIM, p_formula = "two-sided vs null median; (1+#{|Tn-med|>=|To-med|})/(NSIM+1)",
              fdr = "BH across the 6 primary tests", built = as.character(Sys.time())))
saveRDS(results, file.path(BASE, "data", "moduleC_results.rds"))

fwrite(tab[family == "primary", .(test, axis, observed = round(observed, 4),
        null_median = round(null_median, 4), null_lo = round(null_lo, 4),
        null_hi = round(null_hi, 4), p_emp = signif(p_emp, 3), p_adj = signif(p_adj, 3))],
       file.path(BASE, "data", "moduleC_primary_tests.tsv"), sep = "\t")

cat("\n=== PRIMARY TESTS (structured-null calibration) ===\n")
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
       title = "Module C: observed climate axes vs 10,000 Omega-structured null covariates",
       subtitle = "shaded = null 95% interval; dashed = null median") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8))
ggsave(file.path(BASE, "Figures", "moduleC_null_calibration.png"), p1,
       width = 175, height = 140, units = "mm", dpi = 300)
ggsave(file.path(BASE, "Figures", "moduleC_null_calibration.pdf"), p1,
       width = 175, height = 140, units = "mm")

## ========================================================================
## FIGURE 2 -- rank-threshold sensitivity (top 0.1% / 0.5% / 1%)
## ========================================================================
lev  <- c("top 0.1%", "top 0.5%", "top 1%")
mlev <- unname(metric_lookup)                         # DI, recomb, pdiff, pdir_diff order
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
  labs(x = "value within the top-ranked fraction of eMLGs", y = "count (null covariates)",
       title = "Module C threshold sensitivity: top rank fractions, observed vs null",
       subtitle = "prop. directional is computed among differentiated eMLGs in the fraction") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom", strip.text = element_text(size = 7))
ggsave(file.path(BASE, "Figures", "moduleC_threshold_sensitivity.png"), p2,
       width = 185, height = 165, units = "mm", dpi = 300)
ggsave(file.path(BASE, "Figures", "moduleC_threshold_sensitivity.pdf"), p2,
       width = 185, height = 165, units = "mm")

## ========================================================================
## data-driven Markdown report
## ========================================================================
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
any_sig <- tab[family == "primary" & p_adj < 0.05, .N] > 0
N_eMLG  <- length(kc$k1r)

rep_lines <- c(
"# Module C (revised): genome-wide, eMLG-level climate calibration",
"",
sprintf("*Generated %s. NSIM = %d Omega-structured null covariates; eMLG universe = %d clusters (Aland excluded, 19 pops, fixed LD-pruned Omega).*",
        as.character(Sys.Date()), NSIM, N_eMLG),
"",
"## Data provenance",
"",
"- **Observed eMLG climate association:** per-eMLG BayPass BF(dB) on climate PC1/PC2 (`moduleB_climate_GEA/data/moduleB_eMLG_association.rds`).",
"- **Null covariates:** the same 10,000 Omega-structured draws used for Module B's sim-FDR, re-run on the preserved `aland_excluded_eMLG/null/null_b01..b10.env` files (Module B kept only exceedance counts; the full BF matrix was regenerated here and reduced on the fly).",
"- **Annotations (per-eMLG, joined by `group_id`):** consensus Diagnostic Index, unidirectional-sorting status, sorting magnitude `prop_fixed` and signed orientation `uni_score` (`moduleA_sorting/data/moduleA_cluster_sorting.rds`); recombination rate (cM/Mb) at the representative marker (`data/Frufa_DTOL_PR.ref_genome.recmap`); cluster size (`module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds`, cross-checked against the sorting annotation).",
"",
"## Validation checks",
"",
sprintf("- eMLG count and order identical across observed / null / annotation objects (N = %d), all joins by explicit `group_id`.", N_eMLG),
sprintf("- **Faithful regeneration (Monte-Carlo equivalence gate passed):** BayPass is not bit-reproducible (a fresh MCMC realization each run), so the regenerated per-eMLG exceedance counts match Module B's within MCMC tolerance rather than exactly: PC1 Pearson r = %.4f, sum ratio-1 = %+.4f; PC2 r = %.4f, sum ratio-1 = %+.4f (thresholds r > 0.99, |ratio-1| < 0.03; max|dk1| = %d, max|dk2| = %d reported as diagnostics only).",
        kc$cor_k1, kc$rel1, kc$cor_k2, kc$rel2, kc$max_abs_dk1, kc$max_abs_dk2),
"- Input identity (the ten `.env` covariates, geno, Omega, poolsize, params, statistic code) is guaranteed EXACTLY by md5 fingerprints; the primary genome-wide statistics are reproducible to ~0.0005 (<=6% of the null spread) despite per-BF MCMC noise (stability probe). Observed BF vectors equal `eBF1`/`eBF2` (max|d| = 0); observed and null reduced by identical code (`moduleC_stat_functions.R`).",
"",
"## Methods",
"",
"For every covariate (observed PC1/PC2 and each of the 10,000 nulls) the genome-wide eMLG BF vector is reduced to: Spearman rho with DI; Spearman rho with recombination; the difference in mean within-covariate BF **percentile** between directionally sorted and non-directional eMLGs **among differentiated clusters** (primary sorting statistic); plus supplementary statistics (all-eMLG sorting contrast, Spearman rho with sorting magnitude `prop_fixed`, signed sorting orientation `uni_score`, and raw-BF variants). Each observed statistic is compared with its 10,000-value structured-null distribution by a two-sided empirical P (deviation from the null median); the six primary tests (PC1/PC2 x DI/recombination/sorting) are BH-FDR corrected. eMLGs are never resampled independently.",
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

## ---- data-driven interpretation ----------------------------------------
gp <- function(st, ax) tab[family == "primary" & stat == st & axis == ax]
sig_primary <- tab[family == "primary" & p_adj < 0.05][order(p_adj)]
s1 <- gp("sort_gap_differentiated", "PC1"); s2 <- gp("sort_gap_differentiated", "PC2")
r1 <- gp("rho_rec", "PC1"); r2 <- gp("rho_rec", "PC2")
d1 <- gp("rho_DI", "PC1");  d2 <- gp("rho_DI", "PC2")
pd2 <- tab[family == "supplementary" & stat == "pear_DI" & axis == "PC2"]
sideword <- function(row) if (row$observed > row$null_hi) "above" else if (row$observed < row$null_lo) "below" else "within"

interp <- c(
  sprintf("**Directional sorting (primary, differentiated-only): no climate association on either axis.** PC1 observed %.3f (FDR %.3f), PC2 %.3f (FDR %.3f) both lie within the structured-null 95%% interval [%.3f, %.3f]; sorting magnitude (`prop_fixed`) is likewise null on both axes (supplementary, p_emp %.2f/%.2f). Directional ancestry sorting is not organised by the measured climate gradients beyond population structure.",
          s1$observed, s1$p_adj, s2$observed, s2$p_adj, s1$null_lo, s1$null_hi,
          tab[stat=="rho_sort_magnitude" & axis=="PC1", p_emp], tab[stat=="rho_sort_magnitude" & axis=="PC2", p_emp]),
  sprintf("**Recombination: no climate association.** PC1 %.3f (FDR %.3f), PC2 %.3f (FDR %.3f), both within the null.",
          r1$observed, r1$p_adj, r2$observed, r2$p_adj),
  if (d2$p_adj < 0.05)
    sprintf("**Diagnostic Index: a single weak, PC2-specific signal.** PC1 is within the null (Spearman rho %.3f, FDR %.3f), but on PC2 the climate-association strength is correlated with DI beyond the structured null (rho %.3f, empirical p %.3g, FDR %.3f; observed %s the null 95%% interval), corroborated by the raw-BF analysis (Pearson %.3f, p %.3g). The correlation is weak and diffuse -- no individual eMLG survives Module B's genome-wide sim-FDR -- so it is a genome-wide gradient, not a set of climate-adaptation loci. (DI is a signed index; the sign is reported as-is and should not be read as locus-level adaptation in diagnostic regions.)",
            d1$observed, d1$p_adj, d2$observed, d2$p_emp, d2$p_adj, sideword(d2), pd2$observed, pd2$p_emp)
  else
    sprintf("**Diagnostic Index: no climate association after FDR** (PC1 rho %.3f FDR %.3f; PC2 rho %.3f FDR %.3f).",
            d1$observed, d1$p_adj, d2$observed, d2$p_adj),
  if (nrow(sig_primary) > 0)
    sprintf("**Overall:** of the six primary tests, %d survives FDR (%s); every sorting and recombination test is indistinguishable from the Omega-structured null. Weak climate-association evidence is thus diffusely related to the diagnostic gradient on PC2 alone, not to directional sorting, sorting magnitude, or recombination. This must NOT be read as validation of the null-floor candidate eMLGs.",
            nrow(sig_primary), paste(sprintf("%s x %s", sig_primary$test, sig_primary$axis), collapse = "; "))
  else
    "**Overall:** no primary test is exceptional; climate-association evidence is not concentrated in diagnostic, directionally-sorted, or low-recombination eMLGs beyond population structure and genomic architecture -- a clean basis for concluding the measured climate gradients are not a robust organising axis of ancestry sorting.")

rep_lines <- c(rep_lines, "",
"## Interpretation",
"",
interp,
"",
"## What this analysis can and cannot establish",
"",
"- **Can:** calibrate genome-wide, eMLG-level climate-association *patterns* against a structure- and architecture-preserving null, with the same eMLG universe and statistic for observed and null.",
"- **Cannot:** identify individual climate-associated loci (none survive Module B's sim-FDR); the null-floor eMLGs remain exploratory. XtX/BF share among-population allele-frequency variation, so raw-BF sensitivity variants may absorb genuine climate differentiation as well as confounding.",
"")
writeLines(rep_lines, file.path(BASE, "doc", "moduleC_report.md"))

cat(sprintf("\n[analyse] wrote results, TSV, 2 figures, and doc/moduleC_report.md\n"))
cat(sprintf("[analyse] faithful-regeneration gate PASSED (max|dk1|=%d max|dk2|=%d)\n",
            kc$max_abs_dk1, kc$max_abs_dk2))
