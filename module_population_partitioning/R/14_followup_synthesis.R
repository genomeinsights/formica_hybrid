## =========================================================================
## module_population_partitioning -- follow-up 14: synthesis and decision
## table for Analyses 1-3 plus the FST-vs-local-concordance permutation-null
## diagnostic (Analysis 4 stopped -- no verified parental locality metadata,
## see FOLLOWUP_STATUS.md). Does NOT modify any existing pp_*.R script or
## output, or R/10-12.
##
## AUDIT FIX (item 6): this script previously described Analyses 1-3 as
## "three independent follow-up analyses" that found "no evidence of
## widespread recent backcrossing" and used "ruled out" framing. That
## language is removed throughout. The genotype-state run proxy (Analysis 3)
## cannot be used to exclude recent backcrossing (see 12_ancestry_run_proxy.R
## header) and its synthesis-figure panel is replaced with the new
## FST-vs-local-concordance permutation-null diagnostic (item 8), which
## IS methodologically independent of Analyses 1-3.
##
## Combines: geographic prediction (10), individual influence (11), ancestry-
## run proxy (12, exploratory only), FST-concordance null diagnostic into one
## decision table + one 4-panel figure + draft text (methods/results/
## main-text/discussion paragraphs -- NOT inserted into the manuscript
## automatically, per the brief).
##
## Run from the formica_hybrid repo root, after R/10, R/11, R/12 and
## pp_fst_concordance_null.R:
##   Rscript module_population_partitioning/R/14_followup_synthesis.R
## Writes: module_population_partitioning/data/followup/14_followup_synthesis.rds
##         module_population_partitioning/Figures/followup/14_synthesis.png
## =========================================================================
suppressMessages(library(data.table))
t_start <- Sys.time()
OUTDIR <- "module_population_partitioning/data/followup"
FIGDIR <- "module_population_partitioning/Figures/followup"
DATADIR <- "module_population_partitioning/data"

u_check <- readRDS(file.path(DATADIR, "pp_units_Fmat.rds"))$u
stopifnot("unit table must have exactly the 20,807 DI25 rho05 units (min_r2_rho=0.5) -- check for a legacy/stale input" = nrow(u_check) == 20807L)

geo <- readRDS(file.path(OUTDIR, "10_geographic_prediction.rds"))
infl <- readRDS(file.path(OUTDIR, "11_individual_influence.rds"))
run <- readRDS(file.path(OUTDIR, "12_ancestry_run_proxy.rds"))
fstnull <- readRDS(file.path(DATADIR, "pp_fst_concordance_null.rds"))
frac_single_marker <- mean(run$all_runs$is_single_marker)

## ---------------------------------------------------------------------
## 1. combined decision table
## ---------------------------------------------------------------------
decision_table <- data.table(
  item = c("Geographic effect size + permutation (linear lat/long only)",
          "Leave-one-population-out range",
          "Largest individual influence",
          "Effect of excluding influential individuals/populations",
          "FST-vs-local-concordance null check (within-unit label permutation)",
          "Unphased genotype-state run proxy (exploratory only)",
          "Parental differentiation by DI",
          "Limitations and power notes"),
  value = c(
    sprintf("R2=%.3f (adj R2=%.3f); permutation p=%.3f (n=%d reps; observed R2 BELOW null median %.3f) -- tests only LINEAR prediction from latitude/longitude, not non-linear or historical (e.g. colonization-route) spatial structure",
            geo$obs_primary["R2"], geo$obs_primary["R2adj"], geo$p_perm, geo$n_perm, median(geo$perm_R2)),
    sprintf("in-sample adj R2: all 20 folds <= %.4f (none positive); out-of-sample CV R2 = %.3f (strongly negative)",
            max(geo$loo_insample$R2adj), geo$R2_cv),
    sprintf("%s (%s, n=%d individuals) influence_rms=%.3f -- driven by small population size, not flagged as biologically unusual",
            infl$top1_overall, infl$indiv_tab[Sample_ID == infl$top1_overall, Population],
            infl$indiv_tab[Sample_ID == infl$top1_overall, n_pop], infl$indiv_tab[Sample_ID == infl$top1_overall, influence_rms]),
    sprintf("max drift in headline FST-vs-|r| statistic across all targeted exclusion scenarios = %.4f (baseline %.4f) -- %s; NB the top-5-overall scenario wipes out %s entirely (%s), so it is a stress test, not a pure individual-level exclusion",
            infl$key_stat_drift, infl$scenarios[1, rho_FST_vs_absr], infl$verdict,
            paste(infl$top5_wiped_pops, collapse = ","), infl$top5_label),
    sprintf("observed |r| rho=%.4f (signed r rho=%.4f); null (independent within-unit population-label permutation, n=%d reps) 95%% interval |r| [%.4f, %.4f], signed r [%.4f, %.4f] -- observed value falls far outside the null, so the relationship is not a tautological consequence of computing both statistics from the same population-frequency data",
            fstnull$rho_obs_absr, fstnull$rho_obs_r, fstnull$B,
            fstnull$ci_absr[1], fstnull$ci_absr[2], fstnull$ci_r[1], fstnull$ci_r[2]),
    sprintf("%.1f%% of segments are single-marker (primary max_gap); NOT a validated local-ancestry tract analysis and NOT used to infer tract age or exclude recent backcrossing; population-level associations with residual magnitude/high-DI-differentiation contribution are NOT independent tests (same genotypes reused), reported descriptively only",
            100 * frac_single_marker),
    "NOT AVAILABLE -- Analysis 4 stopped: no verified parental Sample_ID -> species/colony/locality/coordinate metadata found in this repository (user confirmed 2026-09-13, chose not to proceed on unverified ID-derived locality)",
    "n=20 populations throughout (low power for weak effects, load-bearing evidence is permutation/cross-validation/LOO-range, not point estimates alone); Analysis 2's scenario recomputation used an approximate (not refit) residualization; the Analysis 3 run proxy is exploratory and unphased and cannot resolve introgression timing or distinguish recent-backcross histories from scattered heterozygosity/genotyping noise; DI was estimated from the same parental individuals any future parental-DI analysis would use (circularity caveat, undischarged)."
  )
)
cat("=== decision table ===\n"); print(decision_table, width = 200)

## ---------------------------------------------------------------------
## 2. verdict
## ---------------------------------------------------------------------
verdict <- paste(
  "Residual ancestry profiles were not detectably predicted by linear geographic coordinates and",
  "the principal population-partitioning results were stable to targeted individual and population",
  "exclusions. An exploratory unphased genotype-state run summary was dominated by short segments",
  "but was not sufficient to infer tract age or exclude recent introgression.",
  "The geographic result concerns only LINEAR prediction from latitude/longitude (Analysis 1); it",
  "does not test non-linear, discrete, or historical (e.g. postglacial colonization route) forms of",
  "spatial or population structure, which remain unaddressed. A separate, methodologically",
  "independent diagnostic (within-unit population-label permutation) supports the",
  "FST-vs-local-concordance relationship being a genuine cross-unit signal rather than a mechanical",
  "consequence of shared computation. Parental geographic/genetic structure remains COMPLETELY",
  "UNTESTED (Analysis 4 stopped for lack of verified metadata) and must be reported as an open",
  "question, not treated as resolved by Analyses 1-3. Together, failing to find support for these",
  "particular alternative explanations is COMPATIBLE WITH, but does NOT BY ITSELF ESTABLISH, an",
  "interpretation of heterogeneous locus-specific sorting as reflecting locus-specific selection or",
  "incompatibility resolution -- none of Analyses 1-3 were designed to test selection or",
  "incompatibility directly, and ancestry-informative loci are broadly differentiated among hybrid",
  "populations while differentiation is assembled from many partly independent, region-specific",
  "ancestry outcomes (linkage causes neighbouring loci to distinguish populations similarly,",
  "especially in low-recombination regions; the corrected row-centred PCA still shows a real,",
  "non-dominant recurring axis, so this is not a claim that every region is independent)."
)
cat("\n=== verdict ===\n"); cat(strwrap(verdict, width = 78), sep = "\n")

## ---------------------------------------------------------------------
## 3. synthesis figure (4 panels). Panel C is the FST-concordance
## permutation-null diagnostic (item 8), replacing the former genotype-run
## panel: the run proxy cannot test introgression timing and per item 6 is
## not shown in the synthesis figure except as clearly-labelled exploratory
## material in its own script's figures (12_*.png).
## ---------------------------------------------------------------------
suppressMessages({ library(ggplot2); library(patchwork) })
theme_ms <- theme_bw(base_size = 11) + theme(strip.background = element_blank(), panel.grid.minor = element_blank())

## A: geographic organization of residual profiles (population-space PCA)
pcs <- geo$pca_scores
panel_A <- ggplot(pcs, aes(Longitude, Latitude, colour = PC1)) +
  geom_point(size = 3) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b", midpoint = 0) +
  labs(title = "A. residual-profile PC1 by geography", x = NULL, y = NULL) + theme_ms + theme(legend.position = "none")

## B: observed geographic effect vs permutation null (LINEAR lat/long only)
panel_B <- ggplot(data.table(R2 = geo$perm_R2), aes(R2)) +
  geom_histogram(bins = 50, fill = "grey75") +
  geom_vline(xintercept = geo$obs_primary["R2"], colour = "firebrick", linewidth = 1) +
  labs(title = sprintf("B. geographic R2 (linear lat/long) vs null (p=%.3f)", geo$p_perm), x = "R2", y = NULL) + theme_ms

## C: FST-vs-local-concordance permutation-null diagnostic (item 8; replaces
## the former genotype-run association panel -- see header note)
panel_C <- ggplot(data.table(rho = fstnull$null_rho_absr), aes(rho)) +
  geom_histogram(bins = 40, fill = "grey75") +
  geom_vline(xintercept = fstnull$rho_obs_absr, colour = "firebrick", linewidth = 1) +
  labs(title = sprintf("C. FST-vs-local-concordance null (within-unit label perm.)\nobserved abs(r)=%.3f, null 95%% [%.3f, %.3f]",
                       fstnull$rho_obs_absr, fstnull$ci_absr[1], fstnull$ci_absr[2]),
      x = "null Spearman rho", y = NULL) + theme_ms

## D (robustness panel, replaces unavailable parental panel): before/after
## the targeted individual/Sielva exclusion scenarios from Analysis 2
sc <- infl$scenarios[scenario %in% c("baseline (none dropped)", sprintf("drop top1 overall (%s)", infl$top1_overall),
                                     infl$top5_label, "drop all Sielva")]
sc[, scenario_short := c("baseline", "drop top1",
                         sprintf("drop top5\n(wipes %s)", paste(infl$top5_wiped_pops, collapse = ",")),
                         "drop Sielva")]
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
 "To probe alternative explanations for the observed heterogeneous, locus-specific sorting of",
 "ancestry, we conducted several follow-up analyses on the DI25 rho05 population-partitioning",
 "dataset (20,807 LD-reduced units, 20 hybrid populations, 165 individuals). First, we tested",
 "whether geographic distance LINEARLY predicts each population's leave-one-chromosome-out",
 "residual ancestry profile, using multivariate OLS (centred latitude/longitude as predictors,",
 "per-unit standardized profiles as the response), a complete-row population-label permutation",
 "test (10,000 replicates), a chromosome-block bootstrap for uncertainty, and leave-one-",
 "population-out cross-validation; this tests only linear spatial prediction, not non-linear or",
 "historical (e.g. colonization-route) forms of structure. Second, we tested whether the",
 "population-level pattern is driven by a small number of individuals, using a closed-form",
 "per-individual influence statistic (the exact population-mean shift from excluding that",
 "individual, normalized by population size) and recomputing the headline FST-vs-local-",
 "concordance statistic under a small set of targeted exclusion scenarios (most influential",
 "individual overall; each population's own most influential member; the five most influential",
 "overall, which fully removes one small population; all of one geographically disjunct",
 "population). Third, we tested whether the FST-vs-local-concordance relationship is a",
 "mechanical consequence of computing both statistics from the same population-frequency data,",
 "by independently permuting population labels within each unit (destroying genuine cross-unit",
 "population correspondence while leaving each unit's own FST unchanged) and recomputing the",
 "correlation under 500 such permutations. Fourth, as an exploratory, non-confirmatory summary,",
 "we implemented an unphased genotype-state run proxy: runs of consecutive, genetic-map-ordered",
 "LD-reduced units with identical ancestry-homozygous genotype state, with missing calls breaking",
 "(not bridging) a run and a maximum genetic-distance gap chosen from the empirical gap",
 "distribution, related descriptively (not as independent tests) to population-level",
 "differentiation and influence metrics via Spearman correlation with leave-one-population-out",
 "uncertainty (n=20 populations, the biological replicate throughout)."
)
draft_results <- paste(
 "Geographic distance did not linearly predict residual ancestry profiles: the observed",
 sprintf("multivariate R2 (%.3f) fell below the median of its own permutation null (%.3f, p=%.3f), and the",
        geo$obs_primary["R2"], median(geo$perm_R2), geo$p_perm),
 sprintf("leave-one-population-out cross-validated R2 was strongly negative (%.3f). The", geo$R2_cv),
 "population-level concordance pattern was not driven by individual or small-population",
 "artefacts: excluding the most influential individual, the five most influential individuals",
 "overall (which fully removes one small population), each population's own most influential",
 "member, or an entire geographically disjunct population (a genuine Alpine site among otherwise-",
 "Fennoscandian populations) shifted the headline FST-vs-local-concordance statistic by at most",
 sprintf("%.3f from a baseline of %.3f. An independent diagnostic that permutes population labels",
        infl$key_stat_drift, infl$scenarios[1, rho_FST_vs_absr]),
 sprintf("within each unit produced a null Spearman correlation centred near zero (95%% interval [%.3f, %.3f])",
        fstnull$ci_absr[1], fstnull$ci_absr[2]),
 sprintf("versus the observed value of %.3f, indicating the FST-vs-local-concordance relationship is",
        fstnull$rho_obs_absr),
 "not a tautological consequence of shared computation. An exploratory unphased genotype-state",
 "run summary was dominated by short segments genome-wide",
 sprintf("(%.1f%% single-marker) and was not sufficient to infer tract age or exclude recent",
        100 * frac_single_marker),
 "introgression; associations between run-length summaries and residual-profile magnitude or",
 "high-DI-differentiation contribution are not independent tests, since both reuse the same",
 "genotype data, and are reported descriptively only."
)
draft_maintext <- paste(
 "Residual ancestry profiles were not detectably predicted by linear geographic coordinates and",
 "the principal population-partitioning results were stable to targeted individual and population",
 "exclusions. An exploratory unphased genotype-state run summary was dominated by short segments",
 "but was not sufficient to infer tract age or exclude recent introgression."
)
draft_discussion <- paste(
 "These results narrow, but do not close, the space of alternative explanations for the observed",
 "pattern. Linear geographic structure and individual/small-population sampling artefacts were",
 "each tested directly and not supported; a targeted diagnostic further indicates the FST-vs-",
 "local-concordance relationship is not a mechanical artefact of shared computation. None of this",
 "establishes locus-specific selection or incompatibility resolution as the explanation -- it is",
 "compatible with, but does not by itself confirm, that interpretation. Two further caveats limit",
 "how far these results generalize. First, the geographic test only addresses LINEAR prediction",
 "from latitude/longitude; non-linear, discrete, or historical (e.g. postglacial colonization",
 "route) forms of spatial or population structure were not tested and remain open. Second, a",
 "fourth alternative, geographic or genetic structure within the parental reference samples",
 "themselves (which could contribute to the high-DI pattern via ascertainment or genuine parental-",
 "source admixture), could not be evaluated: no verified locality metadata exists for the 15+15",
 "parental individuals underlying this analysis, and we did not infer locality from sample",
 "identifiers alone. This is a genuine, currently unresolved gap, not a null result, and should be",
 "flagged as such wherever these follow-up analyses are cited. The unphased genotype-state run",
 "proxy used here is an exploratory, conservative approximation, not a validated local-ancestry",
 "tract analysis; it cannot resolve introgression timing, cannot exclude recent backcrossing, and",
 "short runs must not be read as evidence against recent introgression. A validated local-ancestry",
 "tract caller, if adopted later, should supersede it rather than be treated as confirmatory of the",
 "present proxy's conclusions. Ancestry-informative loci are broadly differentiated among hybrid",
 "populations, but differentiation is assembled from many partly independent, region-specific",
 "ancestry outcomes: linkage causes neighbouring loci to distinguish populations similarly,",
 "especially in low-recombination regions, whereas distant and unlinked regions generally",
 "distinguish different subsets of populations -- the corrected row-centred PCA still shows a",
 "real, non-dominant recurring axis, so this is not a claim that every region behaves independently."
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
