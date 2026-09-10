# Do different high-DI LD units partition the hybrid populations differently?

Revised after independent audit (`../AUDIT.md`) — see `../README.md` "Changes
from the first pass" for what changed and why. This version supersedes the
original summary.

## Interpretation

- **Local partition concordance is real, short-range, and (mostly)
  locus-specific.** Signed correlation between unit population profiles
  decays from 0.258 at <5kb to background within ~0.5–2Mb (chromosome-block
  bootstrap 95% CI), and — critically — this short-range decay is
  **essentially unchanged after residualizing out each population's
  genome-wide ancestry**. That is the strongest single piece of evidence
  that ≤100kb-scale concordance is a genuine, locus-specific phenomenon
  rather than an echo of populations simply differing in overall admixture
  proportion.
- **The (weak) long-range "floor" in the raw distance-decay curve is mostly
  shared ancestry, not signal.** The raw curve does not fully reach the
  empirical cross-chromosome baseline even at >10Mb (0.024 vs 0.027 — close,
  but the point estimate stays on the same side); after residualizing
  against genome-wide ancestry, the long-range floor collapses to ≈0 (0.5–2Mb:
  0.004; 2–10Mb: 0.001; >10Mb: 0.0004). This distinguishes the two
  explanations directly: distant units are *not* independently converging on
  a shared, locus-driven partition — what little raw correlation they show
  is population-level ancestry-tracking.
- **FST predicts local concordance weakly but, under proper uncertainty
  quantification, genuinely.** A chromosome-block bootstrap (not the naive
  `sd/sqrt(N)` used originally, which badly understated uncertainty for
  millions of non-independent pairs) confirms the FST-decile slope is
  non-zero for both signed r and |r|. The effect is small (signed r rises
  from ~0.08 to ~0.13 across FST deciles) and about a fifth of it is
  attributable to genome-wide ancestry-tracking (ρ 0.114→0.092 on residuals)
  — most is not.
- **The 48 large (`n_loci>50`) LD clusters are a distinct, more strongly
  concordant regime (internal ρ=0.41), but they do not drive the genome-wide
  result** — excluding them changes the pooled ρ by <0.003. Avoid describing
  all 48 as "giant low-recombination blocks": only 4 belong to the three
  previously-named polyctena regions, and `n_loci` identifies large LD
  clusters, not independently-verified low-recombination regions (that
  needs the recombination map, not yet joined in).
- **A modest, real (not chance) shared axis exists among sorted units'
  population profiles (PC1=11.0%, above a row-permutation null's 95th
  percentile of 7.2%), and it is disproportionately carried by one
  population: Sielva**, the F1-like colony with elevated heterozygosity
  (loading −0.82, next-largest −0.32 for Åland, all others 0.05–0.16;
  dropping Sielva is the single largest mover of PC1 in a
  leave-one-population-out sweep). This is a different, more specific claim
  than "no shared partition" — a small number of loci that happen to involve
  Sielva (and to a lesser extent Åland) recur more than chance, while the
  genome-wide FST-vs-concordance statistic itself is not disproportionately
  driven by any single population (dropping Sielva moves it only as much as
  dropping most other populations).
- Net picture: strongly differentiated (high-FST, high-DI) loci in this
  dataset generally do **not** share a consistent multilocus population
  partition beyond the scale of direct physical linkage (~100kb) and beyond
  a modest, Sielva/Åland-linked recurring axis — consistent with, but a more
  qualified version of, the originally proposed mechanism for the
  empirical-vs-simulated FST gap.

## Draft Methods (revised)

> For each of the 11,052 DI25 LD-reduced units (from-scratch, hybrids-only
> clustering restricted to DI>−25 markers, 5cM merge cap; best-SNP
> representation for clusters >2 markers, single-marker representative
> otherwise, strictly observed genotypes), we computed an ancestry-oriented
> allele-frequency profile across the 20 hybrid populations, polarised to
> *F. aquilonia* using the parental reference samples, and per-unit
> differentiation as the Weir & Cockerham (1984) FST among the 20 hybrid
> populations. Because profiles are oriented, the signed Pearson correlation
> between two units' profiles is directly interpretable (same vs. opposite
> ancestry direction across populations) and was treated as the primary
> concordance statistic, with the absolute correlation retained as a
> secondary, direction-independent measure. We related concordance to
> physical distance (within-chromosome pairs, binned) against two references:
> genuine cross-chromosome unit pairs (preserving real among-population
> ancestry covariance) and a population-label permutation null (the
> small-n, 20-population sampling floor). Uncertainty on the distance-decay
> and FST-decile trends was estimated by a chromosome-block bootstrap (2000
> replicates resampling the 26 chromosomes with replacement), avoiding the
> pseudoreplication of treating millions of non-independent unit pairs as
> independent observations. To separate locus-specific concordance from a
> shared population-level ancestry gradient, we additionally computed, for
> each population, a leave-one-chromosome-out estimate of its genome-wide
> ancestry, regressed each unit's population profile on this covariate, and
> repeated the full concordance analysis on the residuals.

## Draft Results (revised)

> Signed correlation between DI25 LD-reduced units' population profiles
> decayed from 0.258 [95% CI 0.249–0.267] at <5kb to 0.024 [0.013–0.034]
> beyond 10Mb, close to but not fully reaching the empirical cross-chromosome
> baseline (0.027); |r| showed the same shape, converging on its
> cross-chromosome baseline (0.193) by ~0.5Mb. Per-unit FST was weakly but,
> under a chromosome-block bootstrap, genuinely associated with local
> (≤100kb) concordance (signed r ρ=0.114, |r| ρ=0.081; FST-decile slope 95%
> CI excluded zero for both). This association was not materially driven by
> the 48 largest (>50-marker) LD clusters (excluding them changed ρ by
> <0.003), though these clusters showed markedly stronger internal
> concordance (ρ=0.41) than the rest of the genome (ρ=0.08). Residualizing
> each unit's profile against a leave-one-chromosome-out estimate of
> genome-wide population ancestry left short-range (<100kb) concordance
> essentially unchanged, but collapsed the raw curve's long-range floor to
> ≈0 (2–10Mb: 0.027→0.001), and modestly attenuated the FST-concordance
> association (ρ 0.114→0.092). A row-centred PCA of population profiles
> among sorted units found a real, if modest, shared axis (PC1=11.0% of
> variance, above a row-permutation null's 95th percentile of 7.2%),
> disproportionately loaded on one population (Sielva). These results
> indicate that short-range (<~100kb) partition concordance among
> strongly-differentiated loci is a genuine, largely ancestry-independent,
> locus-specific signal, but it does not extend to a strong shared
> multilocus partition genome-wide: FST and concordance are only weakly
> linked, and most of what shared structure exists is modest and traceable
> to one population.

## Open questions (not modelled here)

1. Is the ρ=0.41 concordance within the 48 `n_loci>50` units itself uniform,
   or concentrated in a subset? Needs within-block fine-scale analysis.
2. The short-range (≤100kb) residual concordance is real but its magnitude
   (signed r ~0.05–0.24 depending on distance) has not been related to local
   recombination rate directly — needs the recombination map joined in
   before using "LD-driven" rather than "consistent with linkage."
3. Sielva's disproportionate role in PC1 — is this a real, repeatable
   population-specific signal (co-adaptation, local introgression) or an
   artefact of Sielva's unusual biology (F1-like, elevated heterozygosity,
   n≈1 colony — see `module_di25/doc/di25_pruning_test_summary.md`)? The
   genome-wide FST-vs-concordance statistic is NOT disproportionately driven
   by Sielva, so this is specifically a PC1/shared-axis question, not a
   headline-result robustness concern.
4. Ultimately distinguishing "many independent selective events on different
   population subsets" from "one shared event whose signature is masked by
   drift/admixture noise at the population level" needs simulation
   (conditional/neutral-null, in the spirit of Module E) — out of scope for
   this pass per the original brief.
5. Comparison against population profiles from the existing neutral
   simulations, needed before attributing the empirical-vs-simulated FST gap
   to this mechanism specifically.
