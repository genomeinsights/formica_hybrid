## =========================================================================
## module_population_partitioning -- follow-up 14: synthesis and decision
## table for Analyses 1-3 (Analysis 4 stopped -- no verified parental
## locality metadata, see FOLLOWUP_STATUS.md). Does NOT modify any existing
## pp_*.R script or output, or R/10-12.
##
## Combines: geographic prediction (10), individual influence (11), ancestry-
## run proxy (12) into one decision table + one 4-panel figure + draft text
## (methods/results/main-text/discussion paragraphs -- NOT inserted into the
## manuscript automatically, per the brief).
##
## Run from the formica_hybrid repo root, after R/10, R/11, R/12:
##   Rscript module_population_partitioning/R/14_followup_synthesis.R
## Writes: module_population_partitioning/data/followup/14_followup_synthesis.rds
##         module_population_partitioning/Figures/followup/14_synthesis.png
## =========================================================================
suppressMessages(library(data.table))
t_start <- Sys.time()
OUTDIR <- "module_population_partitioning/data/followup"
FIGDIR <- "module_population_partitioning/Figures/followup"

geo <- readRDS(file.path(OUTDIR, "10_geographic_prediction.rds"))
infl <- readRDS(file.path(OUTDIR, "11_individual_influence.rds"))
run <- readRDS(file.path(OUTDIR, "12_ancestry_run_proxy.rds"))

## ---------------------------------------------------------------------
## 1. combined decision table
## ---------------------------------------------------------------------
decision_table <- data.table(
  item = c("Geographic effect size + permutation",
          "Leave-one-population-out range",
          "Largest individual influence",
          "Effect of excluding influential individuals",
          "Ancestry-run association",
          "Parental differentiation by DI",
          "Limitations and power notes"),
  value = c(
    sprintf("R2=%.3f (adj R2=%.3f); permutation p=%.3f (n=%d reps; observed R2 BELOW null median %.3f)",
            geo$obs_primary["R2"], geo$obs_primary["R2adj"], geo$p_perm, geo$n_perm, median(geo$perm_R2)),
    sprintf("in-sample adj R2: all 20 folds <= %.4f (none positive); out-of-sample CV R2 = %.3f (strongly negative)",
            max(geo$loo_insample$R2adj), geo$R2_cv),
    sprintf("%s (%s, n=%d individuals) influence_rms=%.3f -- driven by small population size, not flagged as biologically unusual",
            infl$top1_overall, infl$indiv_tab[Sample_ID == infl$top1_overall, Population],
            infl$indiv_tab[Sample_ID == infl$top1_overall, n_pop], infl$indiv_tab[Sample_ID == infl$top1_overall, influence_rms]),
    sprintf("max drift in headline FST-vs-|r| statistic across all targeted exclusion scenarios = %.4f (baseline %.4f) -- %s",
            infl$key_stat_drift, infl$scenarios[1, rho_FST_vs_absr], infl$verdict),
    sprintf("pop. long-run fraction vs residual magnitude rho=%.2f [%.2f,%.2f]; vs high-DI-differentiation contribution rho=%.2f [%.2f,%.2f] (flagged as largely tautological, see Analysis 3); vs individual influence rho=%.2f [%.2f,%.2f] (weak); Sielva a clean near-zero outlier",
            run$assoc[target == "residual-profile magnitude", rho], run$assoc[target == "residual-profile magnitude", loo_min], run$assoc[target == "residual-profile magnitude", loo_max],
            run$assoc[target == "contribution to high-DI differentiation", rho], run$assoc[target == "contribution to high-DI differentiation", loo_min], run$assoc[target == "contribution to high-DI differentiation", loo_max],
            run$assoc[target == "individual influence (Analysis 2)", rho], run$assoc[target == "individual influence (Analysis 2)", loo_min], run$assoc[target == "individual influence (Analysis 2)", loo_max]),
    "NOT AVAILABLE -- Analysis 4 stopped: no verified parental Sample_ID -> species/colony/locality/coordinate metadata found in this repository (user confirmed 2026-09-13, chose not to proceed on unverified ID-derived locality)",
    "n=20 populations throughout (low power for weak effects, load-bearing evidence is permutation/cross-validation/LOO-range, not point estimates alone); Analysis 2's scenario recomputation used an approximate (not refit) residualization; the Analysis 3 run proxy cannot distinguish all recent-backcross histories from scattered heterozygosity/genotyping noise; DI was estimated from the same parental individuals any future parental-DI analysis would use (circularity caveat, undischarged)."
  )
)
cat("=== decision table ===\n"); print(decision_table, width = 200)

## ---------------------------------------------------------------------
## 2. verdict against the brief's Outcome A-E framework
## ---------------------------------------------------------------------
verdict <- paste(
  "Evidence from Analyses 1-3 is most consistent with Outcome A (weak support):",
  "no detectable geographic organization of residual ancestry profiles (Analysis 1),",
  "no individual- or small-population-driven artefact (Analysis 2), and no signature of",
  "pervasive, currently-ongoing large-scale backcrossing (Analysis 3's short run-length",
  "distributions). This STRENGTHENS, but does not by itself ESTABLISH, an interpretation of",
  "heterogeneous locus-specific sorting as reflecting locus-specific selection or",
  "incompatibility resolution rather than ongoing geographically-structured gene flow --",
  "selection is not established merely because these alternative explanations were not",
  "supported. Outcome D/E (parental geographic structure) remains COMPLETELY UNTESTED",
  "(Analysis 4 stopped for lack of verified metadata) and must be reported as an open",
  "question, not treated as resolved by Analyses 1-3."
)
cat("\n=== verdict ===\n"); cat(strwrap(verdict, width = 78), sep = "\n")

## ---------------------------------------------------------------------
## 3. synthesis figure (4 panels; panel D is a robustness panel since
##    parental results are unavailable, per the brief's own fallback rule)
## ---------------------------------------------------------------------
suppressMessages({ library(ggplot2); library(patchwork) })
theme_ms <- theme_bw(base_size = 11) + theme(strip.background = element_blank(), panel.grid.minor = element_blank())

## A: geographic organization of residual profiles (population-space PCA)
pcs <- geo$pca_scores
panel_A <- ggplot(pcs, aes(Longitude, Latitude, colour = PC1)) +
  geom_point(size = 3) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b", midpoint = 0) +
  labs(title = "A. residual-profile PC1 by geography", x = NULL, y = NULL) + theme_ms + theme(legend.position = "none")

## B: observed geographic effect vs permutation null
panel_B <- ggplot(data.table(R2 = geo$perm_R2), aes(R2)) +
  geom_histogram(bins = 50, fill = "grey75") +
  geom_vline(xintercept = geo$obs_primary["R2"], colour = "firebrick", linewidth = 1) +
  labs(title = sprintf("B. geographic R2 vs null (p=%.3f)", geo$p_perm), x = "R2", y = NULL) + theme_ms

## C: individual influence vs ancestry-run relationship (population-level)
targets <- run$targets
panel_C <- ggplot(targets, aes(mean_frac_long_2cM, mean_influence)) +
  geom_point(aes(colour = Population == "Sielva"), size = 2.6) +
  scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey30"), guide = "none") +
  geom_smooth(method = "lm", se = TRUE, colour = "#7570b3", linewidth = 0.7) +
  labs(title = sprintf("C. long-run fraction vs influence (rho=%.2f)", run$assoc[target == "individual influence (Analysis 2)", rho]),
      x = "pop. mean long-run fraction (>=2cM)", y = "pop. mean individual influence") + theme_ms

## D (robustness panel, replaces unavailable parental panel): before/after
## the targeted individual/Sielva exclusion scenarios from Analysis 2
sc <- infl$scenarios[scenario %in% c("baseline (none dropped)", sprintf("drop top1 overall (%s)", infl$top1_overall),
                                     "drop top5 overall", "drop all Sielva")]
sc[, scenario_short := c("baseline", "drop top1", "drop top5", "drop Sielva")]
sc[, scenario_short := factor(scenario_short, levels = scenario_short)]
panel_D <- ggplot(sc, aes(scenario_short, rho_FST_vs_absr)) +
  geom_col(fill = "#1b9e77") +
  labs(title = "D. robustness: FST-vs-concordance rho\nunder targeted exclusion (Analysis 2)", x = NULL, y = "rho") +
  theme_ms + theme(axis.text.x = element_text(angle = 30, hjust = 1))

fig_synth <- (panel_A + panel_B) / (panel_C + panel_D)
ggsave(file.path(FIGDIR, "14_synthesis.png"), fig_synth, width = 10, height = 9, dpi = 200)
cat("\n[synth] figure saved: 14_synthesis.png\n")

## ---------------------------------------------------------------------
## 4. draft text (NOT inserted into the manuscript)
## ---------------------------------------------------------------------
draft_methods <- paste(
 "To test alternative explanations for the observed heterogeneous, locus-specific sorting of",
 "ancestry, we conducted three follow-up analyses on the DI25 rho05 population-partitioning",
 "dataset (20,807 LD-reduced units, 20 hybrid populations, 165 individuals). First, we tested",
 "whether geographic distance predicts each population's leave-one-chromosome-out residual",
 "ancestry profile, using multivariate OLS (centred latitude/longitude as predictors, per-unit",
 "standardized profiles as the response), a complete-row population-label permutation test",
 "(10,000 replicates), a chromosome-block bootstrap for uncertainty, and leave-one-population-out",
 "cross-validation. Second, we tested whether the population-level pattern is driven by a small",
 "number of individuals, using a closed-form per-individual influence statistic (the exact",
 "population-mean shift from excluding that individual, normalized by population size) and",
 "recomputing the headline FST-vs-local-concordance statistic under a small set of targeted",
 "exclusion scenarios (most influential individual overall; each population's own most",
 "influential member; the five most influential overall; all of one geographically disjunct",
 "population). Third, we implemented a conservative, unphased ancestry-run proxy: runs of",
 "consecutive, genetic-map-ordered LD-reduced units with identical ancestry-homozygous genotype",
 "state, with missing calls breaking (not bridging) a run in the primary analysis, related to",
 "population-level differentiation and influence metrics via Spearman correlation with",
 "leave-one-population-out uncertainty (n=20 populations, the biological replicate throughout)."
)
draft_results <- paste(
 "Geographic distance did not predict residual ancestry profiles: the observed multivariate R2",
 sprintf("(%.3f) fell below the median of its own permutation null (%.3f, p=%.3f), and the",
        geo$obs_primary["R2"], median(geo$perm_R2), geo$p_perm),
 sprintf("leave-one-population-out cross-validated R2 was strongly negative (%.3f). The", geo$R2_cv),
 "population-level concordance pattern was not driven by individual or small-population",
 "artefacts: excluding the most influential individual, the five most influential",
 "individuals overall, each population's own most influential member, or an entire",
 "geographically disjunct population (a genuine Alpine site among otherwise-Fennoscandian",
 "populations) shifted the headline FST-vs-local-concordance statistic by at most",
 sprintf("%.3f from a baseline of %.3f.", infl$key_stat_drift, infl$scenarios[1, rho_FST_vs_absr]),
 "An unphased ancestry-run proxy found predominantly short runs genome-wide, arguing against",
 "pervasive, currently-ongoing large-scale backcrossing; the one Alpine population showed",
 "essentially no long ancestry-homozygous runs, consistent with its independently-documented",
 "F1-like heterozygosity. Two population-level associations between run length and",
 "differentiation-related statistics were strong but attributed to measurement overlap rather",
 "than an independent recency signal (see Discussion)."
)
draft_maintext <- paste(
 "Three independent follow-up analyses found no evidence that the heterogeneous, locus-specific",
 "sorting of ancestry documented above is attributable to geographic population structure,",
 "individual-level sampling artefacts, or widespread recent backcrossing."
)
draft_discussion <- paste(
 "These results narrow, but do not close, the space of alternative explanations for the observed",
 "pattern. Geography, individual influence, and coarse recent-introgression signatures were each",
 "tested directly and found wanting, which is consistent with -- though does not establish --",
 "locus-specific selection or incompatibility resolution as an explanation; ruling out several",
 "alternatives is not equivalent to confirming the remaining one. A fourth alternative, geographic",
 "structure within the parental reference samples themselves (which could contribute to the",
 "high-DI pattern via ascertainment or genuine parental-source admixture), could not be evaluated:",
 "no verified locality metadata exists for the 15+15 parental individuals underlying this analysis,",
 "and we did not infer locality from sample identifiers alone. This is a genuine, currently",
 "unresolved gap, not a null result, and should be flagged as such wherever these follow-up",
 "analyses are cited. The ancestry-run proxy used here is a conservative, unphased approximation;",
 "it cannot distinguish all recent backcross histories from chance short runs in an F1-like",
 "individual, and a validated local-ancestry tract caller, if adopted later, should supersede it",
 "rather than be treated as confirmatory of the present proxy's conclusions."
)

cat("\n=== draft text (NOT inserted into the manuscript) ===\n")
for (nm in c("draft_methods", "draft_results", "draft_maintext", "draft_discussion")) {
  cat(sprintf("\n--- %s ---\n", nm)); cat(strwrap(get(nm), width = 78), sep = "\n")
}

result <- list(decision_table = decision_table, verdict = verdict,
              draft_methods = draft_methods, draft_results = draft_results,
              draft_maintext = draft_maintext, draft_discussion = draft_discussion,
              session_info = sessionInfo(), run_time = Sys.time(),
              elapsed_secs = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
saveRDS(result, file.path(OUTDIR, "14_followup_synthesis.rds"))
cat(sprintf("\n[synth] saved -> %s (elapsed %.1fs)\n", file.path(OUTDIR, "14_followup_synthesis.rds"), result$elapsed_secs))
