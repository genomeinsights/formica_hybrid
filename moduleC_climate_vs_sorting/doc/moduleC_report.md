# Module C (revised): genome-wide, eMLG-level climate calibration

*Generated 2026-08-06. NSIM = 10000 Omega-structured null covariates; eMLG universe = 32840 clusters (Aland excluded, 19 pops, fixed LD-pruned Omega).*

## Data provenance

- **Observed eMLG climate association:** per-eMLG BayPass BF(dB) on climate PC1/PC2 (`moduleB_climate_GEA/data/moduleB_eMLG_association.rds`).
- **Null covariates:** the same 10,000 Omega-structured draws used for Module B's sim-FDR, re-run on the preserved `aland_excluded_eMLG/null/null_b01..b10.env` files (Module B kept only exceedance counts; the full BF matrix was regenerated here and reduced on the fly).
- **Annotations (per-eMLG, joined by `group_id`):** consensus Diagnostic Index, unidirectional-sorting status, sorting magnitude `prop_fixed` and signed orientation `uni_score` (`moduleA_sorting/data/moduleA_cluster_sorting.rds`); recombination rate (cM/Mb) at the representative marker (`data/Frufa_DTOL_PR.ref_genome.recmap`); cluster size (`module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds`, cross-checked against the sorting annotation).

## Validation checks

- eMLG count and order identical across observed / null / annotation objects (N = 32840), all joins by explicit `group_id`.
- **Faithful regeneration (Monte-Carlo equivalence gate passed):** BayPass is not bit-reproducible (a fresh MCMC realization each run), so the regenerated per-eMLG exceedance counts match Module B's within MCMC tolerance rather than exactly: PC1 Pearson r = 1.0000, sum ratio-1 = +0.0000; PC2 r = 1.0000, sum ratio-1 = +0.0000 (thresholds r > 0.99, |ratio-1| < 0.03; max|dk1| = 212, max|dk2| = 284 reported as diagnostics only).
- Input identity (the ten `.env` covariates, geno, Omega, poolsize, params, statistic code) is guaranteed EXACTLY by md5 fingerprints; the primary genome-wide statistics are reproducible to ~0.0005 (<=6% of the null spread) despite per-BF MCMC noise (stability probe). Observed BF vectors equal `eBF1`/`eBF2` (max|d| = 0); observed and null reduced by identical code (`moduleC_stat_functions.R`).

## Methods

For every covariate (observed PC1/PC2 and each of the 10,000 nulls) the genome-wide eMLG BF vector is reduced to: Spearman rho with DI; Spearman rho with recombination; the difference in mean within-covariate BF **percentile** between directionally sorted and non-directional eMLGs **among differentiated clusters** (primary sorting statistic); plus supplementary statistics (all-eMLG sorting contrast, Spearman rho with sorting magnitude `prop_fixed`, signed sorting orientation `uni_score`, and raw-BF variants). Each observed statistic is compared with its 10,000-value structured-null distribution by a two-sided empirical P (deviation from the null median); the six primary tests (PC1/PC2 x DI/recombination/sorting) are BH-FDR corrected. eMLGs are never resampled independently.

## Results

| test | axis | observed | null median | null 95% | p_emp | p_adj |
|---|---|---|---|---|---|---|
| DI (Spearman rho) | PC1 | 0.017 | -0.022 | [-0.093, 0.047] | 0.302 | 0.494 |
| DI (Spearman rho) | PC2 | -0.125 | -0.022 | [-0.093, 0.047] | 0.001 | 0.006 |
| recombination (Spearman rho) | PC1 | 0.001 | 0.003 | [-0.019, 0.023] | 0.896 | 0.896 |
| recombination (Spearman rho) | PC2 | -0.006 | 0.003 | [-0.019, 0.023] | 0.449 | 0.539 |
| sorting, differentiated only (BF percentile gap) | PC1 | 0.060 | 0.040 | [0.007, 0.070] | 0.241 | 0.494 |
| sorting, differentiated only (BF percentile gap) | PC2 | 0.023 | 0.040 | [0.007, 0.070] | 0.329 | 0.494 |

### Supplementary statistics (not in the FDR family)

| test | axis | observed | null median | null 95% | p_emp |
|---|---|---|---|---|---|
| sorting, all eMLGs (raw-BF gap) | PC1 | 0.959 | 0.522 | [0.125, 1.001] | 0.051 |
| sorting, all eMLGs (raw-BF gap) | PC2 | 0.401 | 0.522 | [0.125, 1.001] | 0.612 |
| sorting, differentiated (raw-BF gap) | PC1 | 1.043 | 0.650 | [0.251, 1.091] | 0.068 |
| sorting, differentiated (raw-BF gap) | PC2 | 0.665 | 0.650 | [0.251, 1.091] | 0.951 |
| DI (raw-BF Pearson) | PC1 | 0.023 | -0.017 | [-0.085, 0.049] | 0.275 |
| DI (raw-BF Pearson) | PC2 | -0.119 | -0.017 | [-0.085, 0.049] | 6e-04 |
| recombination (raw-BF Pearson) | PC1 | -0.010 | 0.002 | [-0.016, 0.021] | 0.234 |
| recombination (raw-BF Pearson) | PC2 | -0.016 | 0.002 | [-0.016, 0.021] | 0.0683 |
| sorting magnitude (raw-BF Pearson) | PC1 | 0.117 | 0.100 | [0.049, 0.149] | 0.525 |
| sorting magnitude (raw-BF Pearson) | PC2 | 0.126 | 0.100 | [0.049, 0.149] | 0.327 |
| sorting magnitude / prop_fixed (Spearman rho) | PC1 | 0.094 | 0.079 | [0.024, 0.123] | 0.579 |
| sorting magnitude / prop_fixed (Spearman rho) | PC2 | 0.104 | 0.079 | [0.024, 0.123] | 0.359 |
| sorting orientation / uni_score (Spearman rho) | PC1 | 0.061 | -0.016 | [-0.109, 0.073] | 0.11 |
| sorting orientation / uni_score (Spearman rho) | PC2 | -0.124 | -0.016 | [-0.109, 0.073] | 0.0148 |
| sorting, all eMLGs (BF percentile gap) | PC1 | 0.056 | 0.032 | [-0.002, 0.063] | 0.168 |
| sorting, all eMLGs (BF percentile gap) | PC2 | 0.001 | 0.032 | [-0.002, 0.063] | 0.0689 |

### Sensitivity to minimum cluster size and fixation threshold

The calibration is reported over a grid of the minimum eMLG cluster size (min_n_loci in {5, 10}: min=5 is the full 32,840-eMLG universe; min=10 restricts to the 14,349 high-local-LD clusters expected to carry recent-selection signal) and the fixation threshold tau in {0.5, 0.6, 0.8}. Because the eMLG BayPass runs use a FIXED (LD-pruned) Omega, each cluster's BF is invariant to the size threshold, so min=10 is a strict row-subset of the SAME 10,000-null BF matrices, re-reduced over the smaller universe -- no extra BayPass. DI and recombination depend only on min (tau-invariant); directional sorting depends on both. Empirical P only (the FDR family is the six primary tests at the primary cell).

**Directional sorting (differentiated-only) across the grid:**

| min_n_loci | tau | axis | observed | null 95% | p_emp |
|---|---|---|---|---|---|
| 5 | 0.5 | PC1 | 0.056 | [0.006, 0.057] | 0.095 |
| 5 | 0.5 | PC2 | 0.023 | [0.006, 0.057] | 0.457 |
| 5 | 0.6 | PC1 | 0.060 | [0.007, 0.070] | 0.241 |
| 5 | 0.6 | PC2 | 0.023 | [0.007, 0.070] | 0.329 |
| 5 | 0.8 | PC1 | 0.052 | [0.036, 0.095] | 0.378 |
| 5 | 0.8 | PC2 | 0.046 | [0.036, 0.095] | 0.198 |
| 10 | 0.5 | PC1 | 0.069 | [0.008, 0.069] | 0.084 |
| 10 | 0.5 | PC2 | 0.036 | [0.008, 0.069] | 0.749 |
| 10 | 0.6 | PC1 | 0.073 | [0.015, 0.079] | 0.172 |
| 10 | 0.6 | PC2 | 0.043 | [0.015, 0.079] | 0.691 |
| 10 | 0.8 | PC1 | 0.041 | [0.030, 0.119] | 0.168 |
| 10 | 0.8 | PC2 | 0.049 | [0.030, 0.119] | 0.293 |

**DI and recombination across min (tau-invariant), at primary tau = 0.6:**

| test | min_n_loci | axis | observed | null 95% | p_emp |
|---|---|---|---|---|---|
| DI (Spearman rho) | 5 | PC1 | 0.017 | [-0.093, 0.047] | 0.302 |
| DI (Spearman rho) | 5 | PC2 | -0.125 | [-0.093, 0.047] | 0.001 |
| DI (Spearman rho) | 10 | PC1 | -0.007 | [-0.101, 0.024] | 0.371 |
| DI (Spearman rho) | 10 | PC2 | -0.138 | [-0.101, 0.024] | 3e-04 |
| recombination (Spearman rho) | 5 | PC1 | 0.001 | [-0.019, 0.023] | 0.896 |
| recombination (Spearman rho) | 5 | PC2 | -0.006 | [-0.019, 0.023] | 0.449 |
| recombination (Spearman rho) | 10 | PC1 | -0.006 | [-0.028, 0.022] | 0.831 |
| recombination (Spearman rho) | 10 | PC2 | -0.004 | [-0.028, 0.022] | 0.933 |

## Interpretation

**Directional sorting (primary, differentiated-only): no climate association on either axis.** PC1 observed 0.060 (FDR 0.494), PC2 0.023 (FDR 0.494) both lie within the structured-null 95% interval [0.007, 0.070]; sorting magnitude (`prop_fixed`) is likewise null on both axes (supplementary, p_emp 0.58/0.36). Directional ancestry sorting is not organised by the measured climate gradients beyond population structure.
**Recombination: no climate association.** PC1 0.001 (FDR 0.896), PC2 -0.006 (FDR 0.539), both within the null.
**Diagnostic Index: a single weak, PC2-specific signal.** PC1 is within the null (Spearman rho 0.017, FDR 0.494), but on PC2 the climate-association strength is correlated with DI beyond the structured null (rho -0.125, empirical p 0.001, FDR 0.006; observed below the null 95% interval), corroborated by the raw-BF analysis (Pearson -0.119, p 0.0006). The correlation is weak and diffuse -- no individual eMLG survives Module B's genome-wide sim-FDR -- so it is a genome-wide gradient, not a set of climate-adaptation loci. (DI is a signed index; the sign is reported as-is and should not be read as locus-level adaptation in diagnostic regions.)
**Overall:** of the six primary tests, 1 survives FDR (DI (Spearman rho) x PC2); every sorting and recombination test is indistinguishable from the Omega-structured null. Weak climate-association evidence is thus diffusely related to the diagnostic gradient on PC2 alone, not to directional sorting, sorting magnitude, or recombination. This must NOT be read as validation of the null-floor candidate eMLGs.

## What this analysis can and cannot establish

- **Can:** calibrate genome-wide, eMLG-level climate-association *patterns* against a structure- and architecture-preserving null, with the same eMLG universe and statistic for observed and null.
- **Cannot:** identify individual climate-associated loci (none survive Module B's sim-FDR); the null-floor eMLGs remain exploratory. XtX/BF share among-population allele-frequency variation, so raw-BF sensitivity variants may absorb genuine climate differentiation as well as confounding.

