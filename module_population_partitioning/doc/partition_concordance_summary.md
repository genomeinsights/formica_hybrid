# Do different high-DI LD units partition the hybrid populations differently?

Revised 2026-09-13 after independent audit (`../AUDIT.md`) and migration to
the rho05 primary dataset (`../CROSS_MODULE_INPUTS.md`, `../README.md`
"Migration to the rho05 primary dataset"). All numbers below are from the
current primary universe: the DI25-restricted, `min_r2_rho=0.5` clustering,
20,807 units. This version supersedes both earlier summaries.

## Interpretation

- **Local partition concordance is real, short-range, and (mostly)
  locus-specific.** Signed correlation between unit population profiles
  decays from 0.406 at <5kb to background within ~0.5–2Mb (chromosome-block
  bootstrap 95% CI), and — critically — this short-range decay is **largely
  retained after residualizing out each population's genome-wide ancestry**
  (20–100kb: 0.068 residual vs 0.099 raw). That is the strongest single
  piece of evidence that ≤100kb-scale concordance is a genuine,
  locus-specific phenomenon rather than an echo of populations simply
  differing in overall admixture proportion.
- **The (weak) long-range "floor" in the raw distance-decay curve is mostly
  shared ancestry, not signal.** The raw curve does not fully reach the
  empirical cross-chromosome baseline even at >10Mb (0.031 vs 0.027 — close,
  but the point estimate stays on the same side); after residualizing
  against genome-wide ancestry, the long-range floor collapses to ≈0 (0.5–2Mb:
  0.005; 2–10Mb: 0.001; >10Mb: 0.002). This distinguishes the two
  explanations directly: distant units are *not* independently converging on
  a shared, locus-driven partition — what little raw correlation they show
  is population-level ancestry-tracking.
- **FST predicts local concordance weakly but, under proper uncertainty
  quantification, genuinely — and more clearly than under the legacy
  (`min_r2=0.2`) clustering.** A chromosome-block bootstrap confirms the
  FST-decile slope is non-zero for both signed r and |r| (signed r rises
  from ~0.12 to ~0.17 across FST deciles; slope 95% CI [0.0049,0.0062]).
  About 15–17% of this is attributable to genome-wide ancestry-tracking (ρ
  0.144→0.119 signed r on residuals) — most is not.
- **The 50 large (`n_loci>50`) LD clusters are a distinct, more strongly
  concordant regime (internal ρ=0.35), but they do not drive the
  genome-wide result** — excluding them (and separately, excluding the 3
  named blocks, now resolved as 38 rho05 units by physical position rather
  than by reusing legacy group_ids) leaves the pooled ρ unchanged. Avoid
  describing all 50 as "giant low-recombination blocks": `n_loci` identifies
  large LD clusters, not independently-verified low-recombination regions
  (that needs the recombination map, not yet joined in).
- **A modest, real (not chance) shared axis exists among sorted units'
  population profiles (PC1=12.1%, above a row-permutation null's 95th
  percentile of 6.6%), and it is disproportionately carried by one
  population: Sielva**, the F1-like colony with elevated heterozygosity
  (loading 0.77, next-largest 0.40 for Åland; dropping Sielva is the single
  largest mover of PC1 in a leave-one-population-out sweep, PC1 12.1%→10.1%).
  This is a different, more specific claim than "no shared partition" — a
  small number of loci that happen to involve Sielva (and to a lesser extent
  Åland) recur more than chance, while the genome-wide FST-vs-concordance
  statistic itself is not disproportionately driven by any single population
  (leave-one-population-out moves ρ only within 0.125–0.153). **Audit fix
  (item 3, 2026-09-13)**: "sorted units" was corrected from `sort_class !=
  "unsorted"` (which also wrongly included 46 directionally-**unresolved**
  units) to `sort_class %in% c("aquilonia","polyctena")` only; this moved
  PC1 from 11.82% to 12.05% (rounds to 12.1%) — a small correction, not a
  qualitative change in the finding.
- **Diagnostic (item 8, 2026-09-13): is the FST-vs-concordance relationship
  tautological?** Both statistics are computed from the same population-
  allele-frequency matrix, raising the question of whether ρ≈0.14 could be a
  mechanical artefact of shared computation. Independently permuting
  population labels within each unit (destroying genuine cross-unit
  population correspondence while leaving each unit's own FST unchanged)
  gives a null centred near 0 (500 reps: |r| 95% interval [-0.013, 0.014]);
  the observed value falls far outside it, so the relationship is a genuine
  cross-unit signal, not a tautology (`R/pp_fst_concordance_null.R`).
- **Core biological interpretation (unchanged by the audit fixes above,
  explicitly preserved)**: ancestry-informative loci are broadly
  differentiated among hybrid populations, but differentiation is assembled
  from many partly independent, region-specific ancestry outcomes. Linkage
  causes neighbouring loci to distinguish populations similarly, especially
  in low-recombination regions, whereas distant and unlinked regions
  generally distinguish different subsets of populations. This is not a
  claim that every region is independent — the corrected row-centred PCA
  still shows a real, non-dominant recurring axis (see above).
- Net picture: strongly differentiated (high-FST, high-DI) loci in this
  dataset generally do **not** share a consistent multilocus population
  partition beyond the scale of direct physical linkage (~100kb) and beyond
  a modest, Sielva/Åland-linked recurring axis — consistent with, but a more
  qualified version of, the originally proposed mechanism for the
  empirical-vs-simulated FST gap. The rho05 migration strengthened rather
  than weakened this picture: every relationship above is somewhat clearer
  than under the legacy clustering, not an artefact of it.

## Draft Methods (revised)

> For each of the 20,807 DI25 LD-reduced units (from-scratch, hybrids-only
> clustering restricted to DI>−25 markers, 5cM merge cap, decay-relative
> `min_r2_rho=0.5` Stage-2 quality gate; best-SNP representation for
> clusters >2 markers, single-marker representative otherwise, strictly
> observed genotypes), we computed an ancestry-oriented allele-frequency
> profile across the 20 hybrid populations, polarised to *F. aquilonia*
> using the parental reference samples, and per-unit differentiation as the
> Weir & Cockerham (1984) FST among the 20 hybrid populations. Because
> profiles are oriented, the signed Pearson correlation between two units'
> profiles is directly interpretable (same vs. opposite ancestry direction
> across populations) and was treated as the primary concordance statistic,
> with the absolute correlation retained as a secondary, direction-
> independent measure. We related concordance to physical distance
> (within-chromosome pairs, binned) against two references: genuine
> cross-chromosome unit pairs (preserving real among-population ancestry
> covariance) and a population-label permutation null (the small-n,
> 20-population sampling floor). Uncertainty on the distance-decay and
> FST-decile trends was estimated by a chromosome-block bootstrap (2000
> replicates resampling the 26 chromosomes with replacement), avoiding the
> pseudoreplication of treating millions of non-independent unit pairs as
> independent observations. To separate locus-specific concordance from a
> shared population-level ancestry gradient, we additionally computed, for
> each population, a leave-one-chromosome-out estimate of its genome-wide
> ancestry, regressed each unit's population profile on this covariate, and
> repeated the full concordance analysis on the residuals. This unit set is
> deliberately restricted to DI25-ascertained, ancestry-informative markers
> (no additional parental-MAF gate) and kept separate from the full-genome,
> MAF-gated rho05 unit set used for full-DI-range analyses elsewhere in the
> project.

## Draft Results (revised)

> Signed correlation between DI25 LD-reduced units' population profiles
> decayed from 0.406 [95% CI 0.400–0.412] at <5kb to 0.031 [0.015–0.041]
> beyond 10Mb, close to but not fully reaching the empirical cross-chromosome
> baseline (0.027); |r| showed the same shape, converging on its
> cross-chromosome baseline (0.191) by ~0.5Mb. Per-unit FST was weakly but,
> under a chromosome-block bootstrap, genuinely associated with local
> (≤100kb) concordance (signed r ρ=0.144, |r| ρ=0.142; FST-decile slope 95%
> CI excluded zero for both). This association was not materially driven by
> the 50 largest (>50-marker) LD clusters (excluding them left the pooled ρ
> effectively unchanged), though these clusters showed markedly stronger
> internal concordance (ρ=0.35) than the rest of the genome (ρ=0.14).
> Residualizing each unit's profile against a leave-one-chromosome-out
> estimate of genome-wide population ancestry left short-range (<100kb)
> concordance largely intact, but collapsed the raw curve's long-range floor
> to ≈0 (2–10Mb: 0.033→0.001), and only modestly attenuated the
> FST-concordance association (ρ 0.144→0.119). A row-centred PCA of
> population profiles among directionally-sorted (aquilonia/polyctena) units
> found a real, if modest, shared axis (PC1=12.1% of variance, above a
> row-permutation null's 95th percentile of 6.6%), disproportionately loaded
> on one population (Sielva). An independent within-unit population-label
> permutation diagnostic supports the FST-vs-concordance relationship being a
> genuine signal rather than a tautological consequence of shared
> computation (observed ρ=0.14 vs. a null 95% interval of [-0.013, 0.014]).
> These results indicate that short-range (<~100kb) partition concordance
> among strongly-differentiated loci is a genuine, largely
> ancestry-independent, locus-specific signal, but it does not extend to a
> strong shared multilocus partition genome-wide: FST and concordance are
> only weakly linked, and most of what shared structure exists is modest and
> traceable to one population.

## Open questions (not modelled here)

1. Is the ρ=0.35 concordance within the 50 `n_loci>50` clusters itself
   uniform, or concentrated in a subset? Needs within-cluster fine-scale
   analysis.
2. The short-range (≤100kb) residual concordance is real but its magnitude
   has not been related to local recombination rate directly — needs the
   recombination map joined in before using "LD-driven" rather than
   "consistent with linkage."
3. Sielva's disproportionate role in PC1 — is this a real, repeatable
   population-specific signal (co-adaptation, local introgression) or an
   artefact of Sielva's unusual biology (F1-like, elevated heterozygosity,
   n≈1 colony — see `module_di25/doc/di25_pruning_test_summary.md`)? The
   genome-wide FST-vs-concordance statistic is NOT disproportionately driven
   by Sielva, so this is specifically a PC1/shared-axis question, not a
   headline-result robustness concern.
4. Whether the full-genome-units-restricted-to-DI>−25 sensitivity check
   (per README "Two analysis universes") would change any conclusion here —
   not yet run.
5. Ultimately distinguishing "many independent selective events on different
   population subsets" from "one shared event whose signature is masked by
   drift/admixture noise at the population level" needs simulation
   (conditional/neutral-null, in the spirit of Module E) — out of scope for
   this pass per the original brief.
6. Comparison against population profiles from the existing neutral
   simulations, needed before attributing the empirical-vs-simulated FST gap
   to this mechanism specifically.
