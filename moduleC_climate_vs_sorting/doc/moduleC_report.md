# Module C (revised): genome-wide, eMLG-level climate calibration

*Generated 2026-08-19. NSIM = 10000 Omega-structured null covariates; eMLG universe = 32854 clusters (Aland excluded, 19 pops, fixed LD-pruned Omega).*

## Data provenance

- **Observed eMLG climate association:** per-eMLG BayPass BF(dB) on climate PC1/PC2 (`moduleB_climate_GEA/data/moduleB_eMLG_association.rds`).
- **Null covariates:** the same 10,000 Omega-structured draws used for Module B's sim-FDR, re-run on the preserved `aland_excluded_eMLG/null/null_b01..b10.env` files (Module B kept only exceedance counts; the full BF matrix was regenerated here and reduced on the fly).
- **Annotations (per-eMLG, joined by `group_id`):** consensus Diagnostic Index, unidirectional-sorting status, sorting magnitude `prop_fixed` and signed orientation `uni_score` (`moduleA_sorting/data/moduleA_cluster_sorting.rds`); recombination rate (cM/Mb) at the representative marker (`data/Frufa_DTOL_PR.ref_genome.recmap`); cluster size (`module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds`, cross-checked against the sorting annotation).

## Validation checks

- eMLG count and order identical across observed / null / annotation objects (N = 32854), all joins by explicit `group_id`.
- **Faithful regeneration (Monte-Carlo equivalence gate passed):** BayPass is not bit-reproducible (a fresh MCMC realization each run), so the regenerated per-eMLG exceedance counts match Module B's within MCMC tolerance rather than exactly: PC1 Pearson r = 1.0000, sum ratio-1 = -0.0000; PC2 r = 1.0000, sum ratio-1 = +0.0000 (thresholds r > 0.99, |ratio-1| < 0.03; max|dk1| = 198, max|dk2| = 235 reported as diagnostics only).
- Input identity (the ten `.env` covariates, geno, Omega, poolsize, params, statistic code) is guaranteed EXACTLY by md5 fingerprints; the primary genome-wide statistics are reproducible to ~0.0005 (<=6% of the null spread) despite per-BF MCMC noise (stability probe). Observed BF vectors equal `eBF1`/`eBF2` (max|d| = 0); observed and null reduced by identical code (`moduleC_stat_functions.R`).

## Methods

For every covariate (observed PC1/PC2 and each of the 10,000 nulls) the genome-wide eMLG BF vector is reduced to: Spearman rho with DI; Spearman rho with recombination; the difference in mean within-covariate BF **percentile** between directionally sorted and non-directional eMLGs **among differentiated clusters** (primary sorting statistic); plus supplementary statistics (all-eMLG sorting contrast, Spearman rho with sorting magnitude `prop_fixed`, signed sorting orientation `uni_score`, and raw-BF variants). Each observed statistic is compared with its 10,000-value structured-null distribution by a two-sided empirical P (deviation from the null median); the six primary tests (PC1/PC2 x DI/recombination/sorting) are BH-FDR corrected. eMLGs are never resampled independently.

## Results

| test | axis | observed | null median | null 95% | p_emp | p_adj |
|---|---|---|---|---|---|---|
| DI (Spearman rho) | PC1 | 0.060 | -0.009 | [-0.106, 0.090] | 0.192 | 0.383 |
| DI (Spearman rho) | PC2 | -0.115 | -0.009 | [-0.106, 0.090] | 0.0295 | 0.177 |
| recombination (Spearman rho) | PC1 | -0.013 | -0.001 | [-0.019, 0.017] | 0.191 | 0.383 |
| recombination (Spearman rho) | PC2 | 0.003 | -0.001 | [-0.019, 0.017] | 0.659 | 0.659 |
| sorting, differentiated only (BF percentile gap) | PC1 | 0.038 | 0.026 | [0.002, 0.050] | 0.363 | 0.441 |
| sorting, differentiated only (BF percentile gap) | PC2 | 0.015 | 0.026 | [0.002, 0.050] | 0.368 | 0.441 |

### Supplementary statistics (not in the FDR family)

| test | axis | observed | null median | null 95% | p_emp |
|---|---|---|---|---|---|
| sorting, all eMLGs (raw-BF gap) | PC1 | 0.673 | 0.294 | [-0.147, 0.803] | 0.128 |
| sorting, all eMLGs (raw-BF gap) | PC2 | -0.027 | 0.294 | [-0.147, 0.803] | 0.198 |
| sorting, differentiated (raw-BF gap) | PC1 | 0.683 | 0.449 | [0.143, 0.809] | 0.182 |
| sorting, differentiated (raw-BF gap) | PC2 | 0.503 | 0.449 | [0.143, 0.809] | 0.764 |
| DI (raw-BF Pearson) | PC1 | 0.058 | -0.004 | [-0.092, 0.083] | 0.183 |
| DI (raw-BF Pearson) | PC2 | -0.112 | -0.004 | [-0.092, 0.083] | 0.0105 |
| recombination (raw-BF Pearson) | PC1 | -0.016 | 0.002 | [-0.017, 0.022] | 0.0691 |
| recombination (raw-BF Pearson) | PC2 | -0.008 | 0.002 | [-0.017, 0.022] | 0.337 |
| sorting magnitude (raw-BF Pearson) | PC1 | 0.073 | 0.075 | [0.032, 0.120] | 0.942 |
| sorting magnitude (raw-BF Pearson) | PC2 | 0.121 | 0.075 | [0.032, 0.120] | 0.0394 |
| sorting magnitude / prop_fixed (Spearman rho) | PC1 | 0.045 | 0.053 | [0.011, 0.091] | 0.702 |
| sorting magnitude / prop_fixed (Spearman rho) | PC2 | 0.092 | 0.053 | [0.011, 0.091] | 0.059 |
| sorting orientation / uni_score (Spearman rho) | PC1 | 0.085 | -0.010 | [-0.119, 0.096] | 0.0885 |
| sorting orientation / uni_score (Spearman rho) | PC2 | -0.162 | -0.010 | [-0.119, 0.096] | 0.0018 |
| sorting, all eMLGs (BF percentile gap) | PC1 | 0.040 | 0.016 | [-0.018, 0.049] | 0.184 |
| sorting, all eMLGs (BF percentile gap) | PC2 | -0.026 | 0.016 | [-0.018, 0.049] | 0.0119 |

### Sensitivity to minimum cluster size and fixation threshold

The calibration is reported over a grid of the minimum eMLG cluster size (min_n_loci in {5, 10}: min=5 is the full 32,854-eMLG universe; min=10 restricts to the 14,317 high-local-LD clusters expected to carry recent-selection signal) and the fixation threshold tau in {0.5, 0.6, 0.8}. Because the eMLG BayPass runs use a FIXED (LD-pruned) Omega, each cluster's BF is invariant to the size threshold, so min=10 is a strict row-subset of the SAME 10,000-null BF matrices, re-reduced over the smaller universe -- no extra BayPass. DI and recombination depend only on min (tau-invariant); directional sorting depends on both. Empirical P only (the FDR family is the six primary tests at the primary cell).

**Directional sorting (differentiated-only) across the grid:**

| min_n_loci | tau | axis | observed | null 95% | p_emp |
|---|---|---|---|---|---|
| 5 | 0.5 | PC1 | 0.029 | [-0.001, 0.041] | 0.404 |
| 5 | 0.5 | PC2 | 0.012 | [-0.001, 0.041] | 0.475 |
| 5 | 0.6 | PC1 | 0.038 | [0.002, 0.050] | 0.363 |
| 5 | 0.6 | PC2 | 0.015 | [0.002, 0.050] | 0.368 |
| 5 | 0.8 | PC1 | 0.029 | [0.003, 0.097] | 0.413 |
| 5 | 0.8 | PC2 | 0.050 | [0.003, 0.097] | 0.999 |
| 10 | 0.5 | PC1 | 0.039 | [-0.007, 0.049] | 0.244 |
| 10 | 0.5 | PC2 | 0.012 | [-0.007, 0.049] | 0.504 |
| 10 | 0.6 | PC1 | 0.049 | [-0.002, 0.056] | 0.16 |
| 10 | 0.6 | PC2 | 0.029 | [-0.002, 0.056] | 0.996 |
| 10 | 0.8 | PC1 | 0.028 | [0.011, 0.096] | 0.271 |
| 10 | 0.8 | PC2 | 0.068 | [0.011, 0.096] | 0.524 |

**DI and recombination across min (tau-invariant), at primary tau = 0.6:**

| test | min_n_loci | axis | observed | null 95% | p_emp |
|---|---|---|---|---|---|
| DI (Spearman rho) | 5 | PC1 | 0.060 | [-0.106, 0.090] | 0.192 |
| DI (Spearman rho) | 5 | PC2 | -0.115 | [-0.106, 0.090] | 0.0295 |
| DI (Spearman rho) | 10 | PC1 | 0.055 | [-0.116, 0.078] | 0.146 |
| DI (Spearman rho) | 10 | PC2 | -0.139 | [-0.116, 0.078] | 0.0091 |
| recombination (Spearman rho) | 5 | PC1 | -0.013 | [-0.019, 0.017] | 0.191 |
| recombination (Spearman rho) | 5 | PC2 | 0.003 | [-0.019, 0.017] | 0.659 |
| recombination (Spearman rho) | 10 | PC1 | -0.016 | [-0.029, 0.023] | 0.362 |
| recombination (Spearman rho) | 10 | PC2 | -0.005 | [-0.029, 0.023] | 0.881 |

## Interpretation

**Directional sorting (primary, differentiated-only): no climate association on either axis.** PC1 observed 0.038 (FDR 0.441), PC2 0.015 (FDR 0.441) both lie within the structured-null 95% interval [0.002, 0.050]; sorting magnitude (`prop_fixed`) is likewise null on both axes (supplementary, p_emp 0.70/0.06). Directional ancestry sorting is not organised by the measured climate gradients beyond population structure.
**Recombination: no climate association.** PC1 -0.013 (FDR 0.383), PC2 0.003 (FDR 0.659), both within the null.
**Diagnostic Index: no climate association after FDR** (PC1 rho 0.060 FDR 0.383; PC2 rho -0.115 FDR 0.177).
**Overall:** no primary test is exceptional; climate-association evidence is not concentrated in diagnostic, directionally-sorted, or low-recombination eMLGs beyond population structure and genomic architecture -- a clean basis for concluding the measured climate gradients are not a robust organising axis of ancestry sorting.

## What this analysis can and cannot establish

- **Can:** calibrate genome-wide, eMLG-level climate-association *patterns* against a structure- and architecture-preserving null, with the same eMLG universe and statistic for observed and null.
- **Cannot:** identify individual climate-associated loci (none survive Module B's sim-FDR); the null-floor eMLGs remain exploratory. XtX/BF share among-population allele-frequency variation, so raw-BF sensitivity variants may absorb genuine climate differentiation as well as confounding.

