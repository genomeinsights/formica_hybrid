# module_population_partitioning

> **Status (2026-09-13, migrated to rho05).** Exploratory but audited and
> corrected: an independent review (`AUDIT.md`, module root) checked the
> first pass against the saved data, found it directionally right but
> overclaiming in places and one real code bug — see "Changes from the first
> pass" — and recommended a residualization analysis, now implemented. Then
> migrated its primary dataset from the legacy `module_di25` (fixed
> `min_r2=0.2`) clustering to the corrected `module_di25_rho05`
> (`min_r2_rho=0.5`) one — see "Migration to the rho05 primary dataset". Full
> provenance for both this module's inputs and the parallel
> `module_manuscript_rho05` full-genome universe is in
> `CROSS_MODULE_INPUTS.md`. Not yet in the manuscript.

Tests whether different high-DI, LD-reduced units **partition the 20 hybrid
populations differently** — i.e. whether two units can both be strongly
differentiated (high FST) while distinguishing *different* subsets of
populations. If so, high marginal differentiation at individual loci does not
by itself imply consistently elevated multilocus covariance, which is one
candidate explanation for why the empirical FST-DI relationship (FST rising
with diagnostic index) exceeds what the LD-preserving neutral simulations
produce (`module_di25/R/di25_fst_vs_di.R`,
`module_di25/doc/di25_neutral_null_supplement.md`; the rho05-corrected,
full-range counterpart is `module_manuscript_rho05/data/di25_fst_vs_di_rho05.rds`
— see `CROSS_MODULE_INPUTS.md` Table C).

## Two analysis universes (do not conflate — decided 2026-09-13)

This module and `module_manuscript_rho05` deliberately use **different unit
universes** for different questions; per the user's explicit decision, this
module does not migrate onto the full-genome universe wholesale:

- **Primary (this module, always)**: the DI25-restricted, from-scratch,
  rho05 (`min_r2_rho=0.5`) clustering — units defined purely by
  ancestry-informative (DI>−25) variation, never diluted by linked lower-DI
  markers pulled in from a full-genome clustering. No parental-MAF gate
  (the DI>−25 ascertainment IS the diagnostic gate here — see
  `pp_prep_units.R` header). Strictly observed (`fill=FALSE`) best-SNP/
  representative genotypes. **20,807 units.** Used for population-profile
  concordance, geographic prediction, individual-influence diagnostics,
  sorting profiles, ancestry-run proxies — everything in this module.
- **Secondary (module_manuscript_rho05, reference only)**: the full-genome
  rho05 units (395,996 folded-parental-MAF≥0.15-gated / 661,386 total),
  consensus-filled (`fill=TRUE`) genotypes, fixed DI bins spanning the full
  range. For empirical FST across the full DI range, within-species parental
  differentiation, ancestry-uninformative background comparisons, and
  genome-wide ancestry/structure covariates. **Not read by any script in
  this module.** A full-genome-units-restricted-to-DI>−25 analysis may be
  run later as an explicit *sensitivity check* against the primary universe
  above, never as a replacement for it.

Full per-object provenance (paths, producing scripts, checksums, dimensions,
encoded parameters) for both universes: `CROSS_MODULE_INPUTS.md`.

## Conventions (reused, not re-derived)

DI25 from-scratch clustering, **cM5** merge cap, **rho05** (`min_r2_rho=0.5`,
decay-relative) Stage-2 quality gate (`module_di25_rho05/data/di25_clustering_cM5_rho05.rds`,
20,807 units) · unit representation = **best-SNP** for >2-marker clusters
(`eMLG_best_snp(fill=FALSE)`, strictly observed calls) or the cluster's own
representative for 1–2-marker clusters — real single-SNP genotypes throughout
· sorting threshold **τ = 0.6**, φ = 0.85, `sort_rule="binom"`, α = 0.05
(locked Module A convention, re-applied post-hoc, no new threshold) · allele
orientation to *F. aquilonia* via the parental reference samples, matching
`parallelism_stats()`'s internal (unsaved) orientation exactly · FST = Weir &
Cockerham (1984), hybrid populations only, same estimator as
`module_di25/R/di25_fst_vs_di.R::wc_ac()` · **no parental-MAF gate** (by
design — see "Two analysis universes" above).

## Pipeline (`R/`)

| order | script | role |
|---|---|---|
| 1 | `pp_prep_units.R` | Data audit; builds the 20-population × 20,807-unit ancestry-oriented allele-frequency matrix (`Fmat`), per-unit FST, and resolves the 3 named polyctena blocks into the rho05 unit set by physical position. |
| 2 | `pp_local_concordance.R` | All within-chromosome unit-pair similarity (signed r, \|r\|, Euclidean, physical distance); adjacent pairs; distance-binned similarity; FST vs local similarity (adjacent + ≤100kb window), by FST quartile, sort_class, DI decile. |
| 3 | `pp_permutation_null.R` | Two baselines for the large-distance "floor": genuine cross-chromosome unit pairs, and a population-label permutation (pure small-n sampling-noise floor). |
| 4 | `pp_robustness_pca.R` | Big low-recomb clusters (n_loci>50) vs rest; best-SNP vs representative-only units; excluding the 3 named polyctena blocks; row-centered PCA of sorted-unit profiles. |
| 5 | `pp_block_bootstrap.R` | Chromosome-block bootstrap 95% CIs for the distance-bin curve and the FST-decile trend (replaces naive `sd/sqrt(N)`, which treats millions of non-independent pairs as independent). |
| 6 | `pp_pca_refined.R` | Row-permutation null for PC1; PC1/PC2 population loadings; leave-one-population-out PCA; separate aqu- vs pol-sorted PCA; direct test of why the naive (non-row-centered) PCA gives a large PC1. |
| 7 | `pp_residualize_ancestry.R` | **Key analysis**: residualizes each unit's population profile against leave-one-chromosome-out genome-wide ancestry, then re-runs the full concordance analysis on the residuals. |
| 8 | `pp_extra_robustness.R` | Leave-one-population-out (esp. Sielva) and literal once-per-named-region collapse of the genome-wide FST-vs-concordance statistic. |
| 9 | `pp_figures.R` | The 4 figures (below), leading with signed r and block-bootstrap CIs. |

Run from the repo root, in the order above.

## Migration to the rho05 primary dataset (2026-09-13)

Switched `pp_prep_units.R` from the superseded `module_di25/data/di25_clustering_cM5.rds`
+ `di25_sorting_emlg.rds` (fixed `min_r2=0.2`, 11,052 units) to
`module_di25_rho05/data/di25_clustering_cM5_rho05.rds` + `di25_sorting_emlg_rho05.rds`
(`min_r2_rho=0.5`, 20,807 units) — confirmed by reading both rho05 scripts
directly that every other parameter is unchanged (same DI25 panel, same
Stage-1 partition, same cM=5 cap, same `fill=FALSE` convention; ONLY the
Stage-2 quality gate differs). The 3 named polyctena blocks
(`module_di25/data/di25_three_blocks.rds`) were re-resolved into the new
unit set **by physical position** (chromosome + Mb span), not by reusing
their old `group_id`s, which belong to a different partition and would
silently select the wrong units. The full pipeline (steps 1–9) was re-run
end to end; every qualitative conclusion below is unchanged, several
relationships are now visibly stronger (finer clustering → less
LD-pseudoreplication smoothing).

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
   FST-vs-concordance association was "almost entirely driven by" the
   `n_loci>50` clusters. Checked directly: excluding them barely moves ρ.
   Those clusters *do* have a much stronger internal association (a distinct
   high-concordance regime) — but they're too few to move the pooled
   genome-wide statistic. Also corrected: only a minority of them belong to
   the three previously-named polyctena blocks; the label "giant
   low-recombination block" for all of them was dropped since `n_loci>50`
   identifies large LD clusters, not independently-verified low-recombination
   regions.
4. **Cross-chromosome empirical pairs, not the label-permutation null, are
   now the primary large-distance reference** (both give the same |r| ≈
   0.19, but cross-chromosome pairs preserve real among-population
   covariance and are the more defensible baseline for "do physical
   neighbours show excess concordance"). The permutation null is kept as a
   secondary, explicitly-labelled small-n sampling-noise floor.
5. **Chromosome-block bootstrap 95% CIs** (2000 replicates, resampling the 26
   chromosomes with replacement) replace the naive `sd/sqrt(N)` error bars in
   Figs 2–3.
6. **PCA re-examined**: a row-permutation null shows the observed PC1 is
   above chance — a real, if modest, recurring axis, not nothing. Its
   loadings are dominated by one population, **Sielva**; leave-one-
   population-out PCA confirms this (dropping Sielva causes the single
   largest PC1 drop of any population). Separately, the naive
   (non-row-centered) PCA's large PC1 is directly confirmed, not just
   inferred, to be the aquilonia/polyctena sort-direction split: naive PC1
   score cleanly separates aquilonia-sorted from polyctena-sorted units —
   i.e. failing to row-centre just separates the two sorted classes by their
   own baseline level, not a shared population partition.
7. **Residualization against genome-wide ancestry** (the audit's top
   recommendation) — see below.
8. **Leave-one-population-out and once-per-region robustness** for the
   headline FST-vs-concordance statistic (not just PCA): dropping any single
   population, including Sielva, barely moves ρ — no population
   disproportionately drives this particular statistic (contrast with PC1
   above). Collapsing the 3 named blocks to one representative unit each
   changes ρ not at all.

## Key results (rho05 primary dataset, 20,807 units)

- **Local concordance decays fast and is formally uncertainty-quantified.**
  Chromosome-block-bootstrap 95% CI, signed r: 0.406 [0.400,0.412] at 0–5kb →
  0.202 [0.194,0.209] at 5–20kb → 0.099 [0.093,0.105] at 20–100kb → 0.056
  [0.051,0.061] at 100–500kb → 0.038 [0.034,0.042] at 0.5–2Mb → 0.033
  [0.030,0.037] at 2–10Mb → 0.031 [0.015,0.041] at >10Mb, approaching (though
  its point estimate stays a little above) the empirical cross-chromosome
  baseline of 0.027. |r| shows the same shape, converging on the
  cross-chromosome baseline (0.191) by ≈0.5Mb.
- **FST vs local concordance genome-wide is weak but, under a proper
  chromosome-block bootstrap, clearly real — and stronger under rho05 than
  under the legacy clustering**: signed r ρ=0.144, |r| ρ=0.142 (legacy
  values were 0.114/0.081); the FST-decile slope's 95% block-bootstrap CI
  excludes 0 for both (signed r slope 0.0056 [0.0049,0.0062]; |r| slope
  0.0031 [0.0027,0.0034]) — small in absolute terms, but the finer rho05
  partition (less LD-driven pseudoreplication) sharpens rather than weakens
  this relationship. **Diagnostic (item 8, 2026-09-13,
  `R/pp_fst_concordance_null.R`)**: is this a genuine cross-unit signal or a
  mechanical consequence of computing both statistics from the same
  population-frequency matrix? Independently permuting population labels
  within each unit (destroying real cross-unit correspondence while leaving
  each unit's own FST unchanged) gives a null centred near 0 (500 reps: |r|
  95% interval [-0.013, 0.014], signed r [-0.014, 0.012]) — the observed
  ρ≈0.14 falls far outside it, so the relationship is not tautological. See
  `data/pp_fst_concordance_null.rds`.
- **Residualizing against genome-wide ancestry (leave-one-chromosome-out)
  is the key diagnostic, and confirms the same picture under rho05.**
  Per-unit, genome-wide ancestry explains little of most units' among-
  population variation (median R²=0.043, mean R²=0.077; weakly related to
  FST, ρ=0.075). Re-running the full concordance analysis on the residuals:
  **short-range concordance is largely retained** (0–5kb signed r 0.387 vs
  0.406 raw; 20–100kb 0.068 vs 0.099 raw), while **the raw curve's
  slowly-decaying long-range floor collapses to ≈0** (0.5–2Mb: 0.038 raw →
  0.005 residual; 2–10Mb: 0.033 → 0.001; >10Mb: 0.031 → 0.002). The
  FST-vs-concordance association attenuates only modestly on residuals
  (signed r ρ 0.144→0.119; |r| ρ 0.142→0.121) — most of it is not
  ancestry-tracking.
- **No dominant shared partition, but a real modest one, mostly carried by
  Sielva.** Row-centered PCA: PC1=12.1% (vs a row-permutation null 95th
  percentile of 6.6% — real, not chance), decaying gradually. PC1 loadings
  are dominated by Sielva (0.77, next largest 0.40 for Åland); leave-one-
  population-out confirms it (dropping Sielva: PC1 12.1%→10.1%, still the
  largest single-population effect). **Audit fix (item 3, 2026-09-13)**:
  "sorted units" is now correctly restricted to `sort_class %in%
  c("aquilonia","polyctena")` (previously `!= "unsorted"` also wrongly
  included 46 "unresolved" units, reported separately). This is a small
  correction, not a qualitative change — PC1 moved from 11.82% to 12.05%
  (rounds to 12.1% here); the same populations/loadings dominate.
- **Robustness holds**: excluding the 50 `n_loci>50` clusters (internal
  ρ=0.35 vs 0.14 for the rest) leaves the pooled FST-concordance ρ
  unchanged; excluding the 3 named blocks (now 38 rho05 units, resolved by
  physical position) leaves it unchanged; leave-one-population-out moves ρ
  only within 0.125–0.153 (|r|) / 0.131–0.149 (r) of the full-sample
  values; once-per-region collapse changes nothing; best-SNP vs
  representative-only units and current_map_DI decile still do not
  materially change any of the above.

**Interpretation:** outside a handful of physically massive, near-fully-linked
clusters, marginal differentiation (FST) and multilocus partition concordance
are weakly but genuinely related — largely, though not entirely, independent
of the shared population-level ancestry gradient, i.e. more consistent with a
real, if modest, locus-specific mechanism than under the legacy clustering.
The genuinely locus-specific signal that survives ancestry-residualization is
real but short-range (<~100kb) and modest in magnitude. See
`doc/partition_concordance_summary.md` for the full interpretation, revised
draft Methods/Results paragraphs, and remaining open questions.

## Recombination-map join (`pp_recombination.R`, `pp_fig_recombination.R`)

Joins `data/Frufa_DTOL_PR.ref_genome.recmap` (per-chromosome `approx(rule=2)`
interpolation of cM position and local cM/Mb rate, same convention as
`module_di25/R/di25_ld_clustering.R`). Answers the item flagged in "Not yet
run" below — **upgrades "consistent with linkage" to a directly tested,
confirmed recombination effect**:

- Local similarity decays cleanly and monotonically against GENETIC (cM)
  distance (signed r: 0.355 at <0.001cM → 0.037 beyond 2cM), a tighter,
  more universal relationship than against physical distance alone.
- **At a FIXED physical distance (100–500kb), concordance is significantly
  higher in low- vs high-recombination-rate regions**: mean signed r 0.068
  (low tertile) vs 0.049 (high tertile); single-bin contrast, block-bootstrap
  95% CI [0.011, 0.027], clearly excluding 0. (At 0.5–2Mb the contrast
  shrinks to ~0, CI [-0.005, 0.008] — by then most pairs have already
  decayed to background regardless of local recombination rate.) **Audit fix
  (item 4, 2026-09-13)**: the single 100–500kb bin only broadly controls for
  physical distance. Added a distance-adjusted contrast using 20kb-wide
  strata spanning the same 100–500kb scope, combined across strata and
  chromosome-block-bootstrapped (not pair-level SEs): adjusted contrast
  0.018, CI [0.011, 0.027] — matches the single-bin estimate closely and
  **supports** (does not, on its own, "confirm") a genuine recombination-rate
  effect independent of the coarse-bin distance confound. Also checked (and
  found, as expected): exactly 2 units (F11207 on Chr21, F13909 on Chr27)
  fall outside the genetic map's covered physical range and are excluded
  from cM-dependent analyses.
- A unit's own local recombination rate predicts its local (≤100kb)
  concordance to neighbours (Spearman ρ=-0.171 signed r, -0.152 |r|; clean
  monotonic decile trend, 0.196→0.106 across recombination-rate deciles).

Figure: `Figures/fig5_recombination.png`. Data: `data/pp_recombination.rds`.

## Not yet run

- Within-cluster fine-scale structure of the `n_loci>50` regime (is ρ=0.35
  uniform inside those 50 units, or itself driven by a few).
- Inspecting which populations contribute most to each individual high-FST
  unit (beyond the aggregate PC1 loadings already shown).
- The full-genome-units-restricted-to-DI>−25 sensitivity check against this
  module's primary (DI25-specific) unit set, per "Two analysis universes"
  above.
- Comparison against population profiles from the existing simulations,
  needed before attributing the empirical-vs-simulated FST gap to this
  mechanism specifically (out of scope for this pass, per the original
  brief).

## Inputs

`module_di25/data/di25_inputs.rds` (marker panel + genotypes, unchanged by
the rho05 migration), `module_di25_rho05/data/di25_clustering_cM5_rho05.rds`,
`di25_sorting_emlg_rho05.rds`, `module_di25/data/di25_three_blocks.rds`
(resolved by physical position, not group_id) · `moduleA_sorting/R/parallelism_stats.R`
(`classify_sort()` only) · repo-root `data/hybrids_and_parents_maf005.Rdata`.
Full provenance: `CROSS_MODULE_INPUTS.md`.

## Outputs

`data/`: `pp_units_Fmat.rds` (now includes `blk_rho05`, the physically-resolved
named blocks), `pp_concordance_results.rds`, `pp_all_pairs.csv.gz`,
`pp_null_check.rds`, `pp_robustness.rds`,
`pp_block_bootstrap.rds`, `pp_pca_refined.rds`, `pp_residual_ancestry.rds`,
`pp_extra_robustness.rds`, `pp_fst_concordance_null.rds` (within-unit
population-label permutation null for the FST-vs-local-concordance
relationship).

**`pp_units_final.rds` (audit item 1, resolved 2026-09-13)**: this object
was stale (11,052 rows, the pre-rho05-migration lineage — the `saveRDS()`
call that would have refreshed it to 20,807 rows was dropped during an
earlier revision and never re-added). No script read it. It has been moved,
not deleted, to `data/legacy/pp_units_final_minr2_02_11052.rds`, clearly
named by its actual provenance. **`pp_residual_ancestry.rds$u` (20,807
rows) is the canonical current unit table.** Every script in this module
that loads a unit table now asserts, via `stopifnot()`, that it has exactly
20,807 rows with `Fmat`/`Resid` column identity matching `u$group_id` in
order, so accidental use of the legacy 11,052-unit lineage fails loudly
instead of silently propagating stale numbers. See `FOLLOWUP_STATUS.md`
and `CROSS_MODULE_INPUTS.md` for the same resolution note.
`Figures/`: `fig1_heatmap_Chr26.png` (+ zoomed named-block panel),
`fig2_similarity_vs_distance.png` (signed primary / \|r\| secondary,
block-bootstrap CI), `fig3_FST_vs_similarity.png` (same treatment),
`fig4_genomewide_summary.png` (distance-decay raw vs ancestry-residualized,
+ residual genome-wide panel).
`doc/`: `partition_concordance_summary.md` — interpretation, draft
Methods/Results paragraphs, open questions. `AUDIT.md` (module root) — the
independent review this module's statistics were revised in response to.
`CROSS_MODULE_INPUTS.md` (module root) — full input provenance across both
analysis universes.

## Follow-up analyses

`R/1{0,1,2,4}_*.R` — alternative-explanation follow-ups (geographic
prediction, individual influence, an exploratory genotype-state run proxy,
synthesis) plus `R/pp_fst_concordance_null.R` (a within-unit permutation
null for the core FST-vs-concordance result). `R/15_*.R` through
`R/22_*.R` — a separate follow-up importing `module_manuscript_rho05`'s
frozen Stage-1-direct BayPass candidate-locus scans (PC1, PC2, bio_winter,
mitoC2; 18,361 units, 19 populations, Åland excluded — a different unit
universe and population set from this module's own primary pipeline above)
for population-frequency heatmaps, profile-similarity, structure-adjustment,
and matched-null analysis. See `FOLLOWUP_STATUS.md` for full write-ups
(Analyses 1–6 and the FST-concordance diagnostic) and `CROSS_MODULE_INPUTS.md`
Table G for the candidate-locus follow-up's full input provenance.
