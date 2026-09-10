# module_population_partitioning

> **Status (2026-09, revised after AUDIT.md).** Exploratory but now audited
> and corrected: an independent review (`AUDIT.md`, root of this module)
> checked the first pass against the saved data, found it directionally
> right but overclaiming in places and one real code bug, and recommended a
> specific next analysis (residualizing against genome-wide ancestry). All
> of that has been addressed below — see "Changes from the first pass". Not
> yet in the manuscript.

Tests whether different high-DI, LD-reduced units **partition the 20 hybrid
populations differently** — i.e. whether two units can both be strongly
differentiated (high FST) while distinguishing *different* subsets of
populations. If so, high marginal differentiation at individual loci does not
by itself imply consistently elevated multilocus covariance, which is one
candidate explanation for why the empirical FST-DI relationship (FST rising
with diagnostic index) exceeds what the LD-preserving neutral simulations
produce (`module_di25/R/di25_fst_vs_di.R`,
`module_di25/doc/di25_neutral_null_supplement.md`).

## Conventions (reused, not re-derived)

DI25 from-scratch clustering, **cM5** merge cap (`module_di25/data/di25_clustering_cM5.rds`,
11,052 units) · unit representation = **best-SNP** for >2-marker clusters
(`eMLG_best_snp(fill=FALSE)`, strictly observed calls) or the cluster's own
representative for 1–2-marker clusters — real single-SNP genotypes throughout
· sorting threshold **τ = 0.6**, φ = 0.85, `sort_rule="binom"`, α = 0.05
(locked Module A convention, re-applied post-hoc, no new threshold) · allele
orientation to *F. aquilonia* via the parental reference samples, matching
`parallelism_stats()`'s internal (unsaved) orientation exactly · FST = Weir &
Cockerham (1984), hybrid populations only, same estimator as
`module_di25/R/di25_fst_vs_di.R::wc_ac()`.

## Pipeline (`R/`)

| order | script | role |
|---|---|---|
| 1 | `pp_prep_units.R` | Data audit; builds the 20-population × 11,052-unit ancestry-oriented allele-frequency matrix (`Fmat`) and per-unit FST. |
| 2 | `pp_local_concordance.R` | All within-chromosome unit-pair similarity (signed r, \|r\|, Euclidean, physical distance); adjacent pairs; distance-binned similarity; FST vs local similarity (adjacent + ≤100kb window), by FST quartile, sort_class, DI decile. |
| 3 | `pp_permutation_null.R` | Two baselines for the large-distance "floor": genuine cross-chromosome unit pairs, and a population-label permutation (pure small-n sampling-noise floor). |
| 4 | `pp_robustness_pca.R` | Big low-recomb blocks (n_loci>50) vs rest; best-SNP vs representative-only units; excluding the 3 named polyctena blocks; row-centered PCA of sorted-unit profiles. |
| 5 | `pp_block_bootstrap.R` | Chromosome-block bootstrap 95% CIs for the distance-bin curve and the FST-decile trend (replaces naive `sd/sqrt(N)`, which treats millions of non-independent pairs as independent). |
| 6 | `pp_pca_refined.R` | Row-permutation null for PC1; PC1/PC2 population loadings; leave-one-population-out PCA; separate aqu- vs pol-sorted PCA; direct test of why the naive (non-row-centered) PCA gives PC1=71%. |
| 7 | `pp_residualize_ancestry.R` | **Key analysis**: residualizes each unit's population profile against leave-one-chromosome-out genome-wide ancestry, then re-runs the full concordance analysis on the residuals. |
| 8 | `pp_extra_robustness.R` | Leave-one-population-out (esp. Sielva) and literal once-per-named-region collapse of the genome-wide FST-vs-concordance statistic. |
| 9 | `pp_figures.R` | The 4 figures (below), revised to lead with signed r and block-bootstrap CIs. |

Run from the repo root, in the order above.

## Changes from the first pass (per AUDIT.md — verified independently before acting)

1. **Bug fixed**: `pp_local_concordance.R`'s vectorised Euclidean distance used
   `colSums()` on a pairs×populations matrix (returns one value per
   population, silently recycled) instead of `rowSums()` (one value per
   pair) — confirmed by hand, and the source of the "26 warnings" glossed
   over in the first pass. Fixed; Euclidean is otherwise redundant with r for
   this near-complete data (`eucl = sqrt(2·19·(1−r))/√20`) and was never used
   in any conclusion.
2. **Signed r is now the primary statistic**, |r| secondary. Orientation to
   *F. aquilonia* makes the sign biologically meaningful (same direction /
   opposite direction / unrelated), and |r| has a positive floor purely from
   sampling 20 populations, inflating apparent concordance for any two units.
3. **Corrected an overclaim**: the first pass said the genome-wide
   FST-vs-concordance association was "almost entirely driven by" the 48
   `n_loci>50` clusters. Checked directly: excluding them moves ρ from 0.081
   to 0.079 — negligible. Those 48 units *do* have a much stronger internal
   association (ρ=0.41) — a distinct high-concordance regime — but they're
   too few (48/11,052) to move the pooled genome-wide statistic. Also
   corrected: only 4 of those 48 units belong to the three previously-named
   polyctena blocks, not "mostly" them; the label "giant low-recombination
   block" for all 48 was dropped since `n_loci>50` identifies large LD
   clusters, not independently-verified low-recombination regions.
4. **Cross-chromosome empirical pairs, not the label-permutation null, are
   now the primary large-distance reference** (both give the same |r| ≈
   0.19, but cross-chromosome pairs preserve real among-population
   covariance and are the more defensible baseline for "do physical
   neighbours show excess concordance"). The permutation null is kept as a
   secondary, explicitly-labelled small-n sampling-noise floor.
5. **Chromosome-block bootstrap 95% CIs** (2000 replicates, resampling the 26
   chromosomes with replacement) replace the naive `sd/sqrt(N)` error bars in
   Figs 2–3.
6. **PCA re-examined**: a row-permutation null shows the observed PC1 (11.0%)
   is above chance (permutation null 95th pct 7.2%) — a real, if modest,
   recurring axis, not nothing. Its loadings are dominated by one population,
   **Sielva** (loading −0.82, next largest −0.32 for Åland, all others
   0.05–0.16); leave-one-population-out PCA confirms this (dropping Sielva:
   PC1 11.0%→9.9%, the largest of any population's removal; every other
   population's removal leaves PC1 within 10.9–11.8%). Separately, the naive
   (non-row-centered) PCA's PC1=71.5% is now directly confirmed, not just
   inferred, to be the aquilonia/polyctena sort-direction split: mean naive
   PC1 score is −0.8 for aquilonia-sorted units vs +2.5 for polyctena-sorted
   (SD ≈0.3 each) — i.e. failing to row-centre just separates the two sorted
   classes by their own baseline level, not a shared population partition.
7. **Residualization against genome-wide ancestry** (the audit's top
   recommendation) — see below.
8. **Leave-one-population-out and once-per-region robustness** for the
   headline FST-vs-concordance statistic (not just PCA): dropping any single
   population, including Sielva, moves ρ only within 0.069–0.085 (|r|) /
   0.103–0.118 (r) of the full-sample 0.081/0.114 — no population
   disproportionately drives this particular statistic (contrast with PC1
   above). Collapsing the 3 named blocks to one representative unit each
   changes ρ not at all (0.081→0.081, 0.114→0.114).

## Key results (revised)

- **Local concordance decays fast and is now formally uncertainty-quantified.**
  Chromosome-block-bootstrap 95% CI, signed r: 0.258 [0.249,0.267] at 0–5kb →
  0.146 [0.139,0.151] at 5–20kb → 0.076 [0.072,0.081] at 20–100kb → 0.044
  [0.041,0.048] at 100–500kb → 0.030 [0.027,0.034] at 0.5–2Mb → 0.027
  [0.024,0.031] at 2–10Mb → 0.024 [0.013,0.034] at >10Mb, approaching (though
  its point estimate stays a little above) the empirical cross-chromosome
  baseline of 0.027. |r| shows the same shape, converging on the
  cross-chromosome baseline (0.193) by ≈0.5Mb.
- **FST vs local concordance genome-wide is weak but, under a proper
  chromosome-block bootstrap, real**: signed r ρ=0.114, |r| ρ=0.081; the
  FST-decile slope's 95% block-bootstrap CI excludes 0 for both (signed r
  slope 0.0042 [0.0036,0.0049]; |r| slope 0.0016 [0.0012,0.0019]) — small,
  but not naive-SE noise.
- **Residualizing against genome-wide ancestry (leave-one-chromosome-out)
  is the key diagnostic.** Per-unit, genome-wide ancestry explains little of
  most units' among-population variation (median R²=0.037, mean R²=0.070;
  weakly related to FST, ρ=0.053) — most of a unit's profile is unit-specific,
  not a shared-ancestry echo. Re-running the full concordance analysis on the
  residuals: **short-range concordance is essentially unchanged** (0–5kb
  signed r 0.244 raw-comparable; 20–100kb 0.053), while **the raw curve's
  slowly-decaying long-range floor collapses to ≈0** (0.5–2Mb: 0.030 raw →
  0.004 residual; 2–10Mb: 0.027 → 0.001; >10Mb: 0.024 → 0.0004). This
  directly distinguishes the two candidate explanations the audit posed:
  the modest long-range floor in the raw curve was mostly the shared
  genome-wide ancestry gradient, not locus-specific signal, whereas the
  short-range (<~100kb) concordance is genuinely locus-specific and survives
  removing that gradient. The FST-vs-concordance association also
  attenuates somewhat on residuals (signed r ρ 0.114→0.092; |r| ρ
  0.081→0.059) but does not vanish — part of it, not all, is ancestry-tracking.
- **No dominant shared partition, but a real modest one, mostly carried by
  Sielva.** Row-centered PCA: PC1=11.0% (vs a row-permutation null 95th
  percentile of 7.2% — real, not chance), decaying gradually (46% cumulative
  by PC6). PC1 is disproportionately driven by Sielva (the F1-like colony
  with elevated heterozygosity) and, to a lesser extent, Åland.
- Best-SNP (`is_emlg`) vs representative-only units, and current_map_DI
  decile, still do not materially change any of the above.

**Interpretation:** outside a handful of physically massive, near-fully-linked
blocks, marginal differentiation (FST) and multilocus partition concordance
are weakly related but not decoupled — the relationship is small and mostly,
though not entirely, attributable to a shared population-level ancestry
gradient rather than a strong, independent, locus-specific mechanism. The
genuinely locus-specific signal that survives ancestry-residualization is
real but short-range (<~100kb) and modest in magnitude. See
`doc/partition_concordance_summary.md` for the full interpretation, revised
draft Methods/Results paragraphs, and remaining open questions.

## Not yet run

- Within-block fine-scale structure of the `n_loci>50` regime (is ρ=0.41
  uniform inside those 48 units, or itself driven by a few).
- Joining the recombination map to test whether the ≤100kb decay length
  varies with local recombination rate (still describe as "consistent with
  linkage", not "LD-driven", until this is done).
- Inspecting which populations contribute most to each individual high-FST
  unit (beyond the aggregate PC1 loadings already shown).
- Comparison against population profiles from the existing simulations,
  needed before attributing the empirical-vs-simulated FST gap to this
  mechanism specifically (out of scope for this pass, per the original
  brief).

## Inputs

`module_di25/data/di25_inputs.rds`, `di25_clustering_cM5.rds`,
`di25_sorting_emlg.rds`, `di25_three_blocks.rds` · `moduleA_sorting/R/parallelism_stats.R`
(`classify_sort()` only) · repo-root `data/hybrids_and_parents_maf005.Rdata`.

## Outputs

`data/`: `pp_units_Fmat.rds`, `pp_concordance_results.rds`, `pp_all_pairs.csv.gz`,
`pp_null_check.rds`, `pp_robustness.rds`, `pp_units_final.rds`,
`pp_block_bootstrap.rds`, `pp_pca_refined.rds`, `pp_residual_ancestry.rds`,
`pp_extra_robustness.rds`.
`Figures/`: `fig1_heatmap_Chr26.png` (+ zoomed named-block panel),
`fig2_similarity_vs_distance.png` (signed primary / \|r\| secondary,
block-bootstrap CI), `fig3_FST_vs_similarity.png` (same treatment),
`fig4_genomewide_summary.png` (distance-decay raw vs ancestry-residualized,
+ residual genome-wide panel).
`doc/`: `partition_concordance_summary.md` — interpretation, draft
Methods/Results paragraphs, open questions. `AUDIT.md` (module root) —
the independent review this revision responds to.
