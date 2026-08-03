# Module C (revised): genome-wide, eMLG-level climate calibration

*Generated 2026-08-03. NSIM = 10000 Omega-structured null covariates; eMLG universe = 32840 clusters (Aland excluded, 19 pops, fixed LD-pruned Omega).*

## Data provenance

- **Observed eMLG climate association:** per-eMLG BayPass BF(dB) on climate PC1/PC2 (`moduleB_climate_GEA/data/moduleB_eMLG_association.rds`).
- **Null covariates:** the same 10,000 Omega-structured draws used for Module B's sim-FDR, re-run on the preserved `aland_excluded_eMLG/null/null_b01..b10.env` files (Module B kept only exceedance counts; the full BF matrix was regenerated here and reduced on the fly).
- **Annotations (per-eMLG, joined by `group_id`):** consensus Diagnostic Index, unidirectional-sorting status, sorting magnitude `prop_fixed` and signed orientation `uni_score` (`data/moduleC_C3_cl.rds`); recombination rate (cM/Mb) at the representative marker (`data/Frufa_DTOL_PR.ref_genome.recmap`); cluster size (`module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds`, cross-checked against C3_cl).

## Validation checks

- eMLG count and order identical across observed / null / annotation objects (N = 32840), all joins by explicit `group_id`.
- **Faithful regeneration (Monte-Carlo equivalence gate passed):** BayPass is not bit-reproducible (a fresh MCMC realization each run), so the regenerated per-eMLG exceedance counts match Module B's within MCMC tolerance rather than exactly: PC1 Pearson r = 1.0000, sum ratio-1 = +0.0000; PC2 r = 1.0000, sum ratio-1 = +0.0000 (thresholds r > 0.99, |ratio-1| < 0.03; max|dk1| = 259, max|dk2| = 256 reported as diagnostics only).
- Input identity (the ten `.env` covariates, geno, Omega, poolsize, params, statistic code) is guaranteed EXACTLY by md5 fingerprints; the primary genome-wide statistics are reproducible to ~0.0005 (<=6% of the null spread) despite per-BF MCMC noise (stability probe). Observed BF vectors equal `eBF1`/`eBF2` (max|d| = 0); observed and null reduced by identical code (`moduleC_stat_functions.R`).

## Methods

For every covariate (observed PC1/PC2 and each of the 10,000 nulls) the genome-wide eMLG BF vector is reduced to: Spearman rho with DI; Spearman rho with recombination; the difference in mean within-covariate BF **percentile** between directionally sorted and non-directional eMLGs **among differentiated clusters** (primary sorting statistic); plus supplementary statistics (all-eMLG sorting contrast, Spearman rho with sorting magnitude `prop_fixed`, signed sorting orientation `uni_score`, and raw-BF variants). Each observed statistic is compared with its 10,000-value structured-null distribution by a two-sided empirical P (deviation from the null median); the six primary tests (PC1/PC2 x DI/recombination/sorting) are BH-FDR corrected. eMLGs are never resampled independently.

## Results

| test | axis | observed | null median | null 95% | p_emp | p_adj |
|---|---|---|---|---|---|---|
| DI (Spearman rho) | PC1 | 0.017 | -0.022 | [-0.093, 0.047] | 0.303 | 0.571 |
| DI (Spearman rho) | PC2 | -0.125 | -0.022 | [-0.093, 0.047] | 0.001 | 0.006 |
| recombination (Spearman rho) | PC1 | 0.001 | 0.002 | [-0.019, 0.023] | 0.909 | 0.909 |
| recombination (Spearman rho) | PC2 | -0.006 | 0.002 | [-0.019, 0.023] | 0.456 | 0.571 |
| sorting, differentiated only (BF percentile gap) | PC1 | 0.055 | 0.033 | [0.007, 0.057] | 0.107 | 0.32 |
| sorting, differentiated only (BF percentile gap) | PC2 | 0.024 | 0.033 | [0.007, 0.057] | 0.475 | 0.571 |

### Supplementary statistics (not in the FDR family)

| test | axis | observed | null median | null 95% | p_emp |
|---|---|---|---|---|---|
| sorting, all eMLGs (raw-BF gap) | PC1 | 0.797 | 0.418 | [0.103, 0.788] | 0.0315 |
| sorting, all eMLGs (raw-BF gap) | PC2 | 0.326 | 0.418 | [0.103, 0.788] | 0.622 |
| sorting, differentiated (raw-BF gap) | PC1 | 0.895 | 0.562 | [0.246, 0.886] | 0.0409 |
| sorting, differentiated (raw-BF gap) | PC2 | 0.613 | 0.562 | [0.246, 0.886] | 0.763 |
| DI (raw-BF Pearson) | PC1 | 0.023 | -0.017 | [-0.085, 0.049] | 0.275 |
| DI (raw-BF Pearson) | PC2 | -0.119 | -0.017 | [-0.085, 0.049] | 7e-04 |
| recombination (raw-BF Pearson) | PC1 | -0.010 | 0.002 | [-0.017, 0.021] | 0.243 |
| recombination (raw-BF Pearson) | PC2 | -0.016 | 0.002 | [-0.017, 0.021] | 0.0708 |
| sorting magnitude (raw-BF Pearson) | PC1 | 0.117 | 0.100 | [0.048, 0.149] | 0.525 |
| sorting magnitude (raw-BF Pearson) | PC2 | 0.126 | 0.100 | [0.048, 0.149] | 0.33 |
| sorting magnitude / prop_fixed (Spearman rho) | PC1 | 0.094 | 0.079 | [0.024, 0.122] | 0.575 |
| sorting magnitude / prop_fixed (Spearman rho) | PC2 | 0.104 | 0.079 | [0.024, 0.122] | 0.357 |
| sorting orientation / uni_score (Spearman rho) | PC1 | 0.061 | -0.016 | [-0.109, 0.073] | 0.109 |
| sorting orientation / uni_score (Spearman rho) | PC2 | -0.124 | -0.016 | [-0.109, 0.073] | 0.0148 |
| sorting, all eMLGs (BF percentile gap) | PC1 | 0.050 | 0.025 | [-0.002, 0.050] | 0.0594 |
| sorting, all eMLGs (BF percentile gap) | PC2 | 0.000 | 0.025 | [-0.002, 0.050] | 0.0714 |

## Interpretation

**Directional sorting (primary, differentiated-only): no climate association on either axis.** PC1 observed 0.055 (FDR 0.320), PC2 0.024 (FDR 0.571) both lie within the structured-null 95% interval [0.007, 0.057]; sorting magnitude (`prop_fixed`) is likewise null on both axes (supplementary, p_emp 0.58/0.36). Directional ancestry sorting is not organised by the measured climate gradients beyond population structure.
**Recombination: no climate association.** PC1 0.001 (FDR 0.909), PC2 -0.006 (FDR 0.571), both within the null.
**Diagnostic Index: a single weak, PC2-specific signal.** PC1 is within the null (Spearman rho 0.017, FDR 0.571), but on PC2 the climate-association strength is correlated with DI beyond the structured null (rho -0.125, empirical p 0.001, FDR 0.006; observed below the null 95% interval), corroborated by the raw-BF analysis (Pearson -0.119, p 0.0007). The correlation is weak and diffuse -- no individual eMLG survives Module B's genome-wide sim-FDR -- so it is a genome-wide gradient, not a set of climate-adaptation loci. (DI is a signed index; the sign is reported as-is and should not be read as locus-level adaptation in diagnostic regions.)
**Overall:** of the six primary tests, 1 survives FDR (DI (Spearman rho) x PC2); every sorting and recombination test is indistinguishable from the Omega-structured null. Weak climate-association evidence is thus diffusely related to the diagnostic gradient on PC2 alone, not to directional sorting, sorting magnitude, or recombination. This must NOT be read as validation of the null-floor candidate eMLGs.

## What this analysis can and cannot establish

- **Can:** calibrate genome-wide, eMLG-level climate-association *patterns* against a structure- and architecture-preserving null, with the same eMLG universe and statistic for observed and null.
- **Cannot:** identify individual climate-associated loci (none survive Module B's sim-FDR); the null-floor eMLGs remain exploratory. XtX/BF share among-population allele-frequency variation, so raw-BF sensitivity variants may absorb genuine climate differentiation as well as confounding.

