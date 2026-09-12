# Module C, Stage-1-direct: genome-wide, unit-level climate calibration

*Generated 2026-09-12. NSIM = 10000 Omega-structured null covariates; Stage-1-direct unit universe = 18361 clusters (n_snps>=5, Aland excluded, 19 pops, Stage-1-derived Omega). Single universe -- unlike the canonical eMLG Module C, there is no second min-cluster-size level to sweep.*

## Data provenance

- **Observed Stage-1-unit climate association:** per-unit BayPass BF(dB) on climate PC1/PC2 (`module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/PC{1,2}_S1units_withOmega_summary_betai_reg.out`).
- **Null covariates:** the same 10,000 Omega-structured draws used for the Stage-1-unit sim-FDR floor test (`moduleB_stage1_S1units_null.R`, mini2), re-run on the preserved `null/null_b01..b50.env` files (that run kept only exceedance counts; the full BF matrix is regenerated here, persisted, and reduced on the fly).
- **Annotations (per-unit, joined by `group_id`):** consensus best-SNP Diagnostic Index, unidirectional-sorting status, sorting magnitude `prop_fixed` and signed orientation `uni_score` (`moduleA_stage1_cluster_sorting.rds`, Module A equivalent for the Stage-1 universe); recombination rate (cM/Mb) at the best-SNP marker (`data/Frufa_DTOL_PR.ref_genome.recmap`); cluster size (`moduleB_stage1_units_bestsnp.rds`).

## Validation checks

- Stage-1 unit count and order identical across observed / null / annotation objects (N = 18361), all joins by explicit `group_id`.
- **Faithful regeneration (Monte-Carlo equivalence gate passed):** BayPass is not bit-reproducible (a fresh MCMC realization each run), so the regenerated per-unit exceedance counts match the moduleB_stage1_S1units_null.R run within MCMC tolerance rather than exactly: PC1 Pearson r = 1.0000, sum ratio-1 = +0.0000; PC2 r = 1.0000, sum ratio-1 = +0.0000 (thresholds r > 0.99, |ratio-1| < 0.03; max|dk1| = 0, max|dk2| = 0 reported as diagnostics only).
- Input identity (the 50 `.env` covariates, geno, Omega, poolsize, params, statistic code) is guaranteed EXACTLY by md5 fingerprints. Observed BF vectors equal `eBF1`/`eBF2` (max|d| = 0); observed and null reduced by identical code (`moduleC_stat_functions.R`, shared unmodified with the canonical eMLG Module C).
- The `*_pdir_diff` threshold-sensitivity statistics (proportion directional among differentiated units, within the top BF fraction) are undefined (0/0) whenever a null draw's top fraction contains zero differentiated units -- expected for the smallest fraction (top 0.1%, ~18 units) given this smaller 18,361-unit universe. Those draws are excluded from the corresponding p-value/quantile calculation rather than imputed; NA counts (of 10,000 nulls): top0001_pdir_diff=11, top0005_pdir_diff=0, top0010_pdir_diff=0. No other statistic -- including the primary FDR family -- is ever NA.

## Methods

For every covariate (observed PC1/PC2 and each of the 10,000 nulls) the genome-wide Stage-1-unit BF vector is reduced to: Spearman rho with DI; Spearman rho with recombination; the difference in mean within-covariate BF **percentile** between directionally sorted and non-directional units **among differentiated clusters** (primary sorting statistic); plus supplementary statistics (all-unit sorting contrast, Spearman rho with sorting magnitude `prop_fixed`, signed sorting orientation `uni_score`, and raw-BF variants). Each observed statistic is compared with its 10,000-value structured-null distribution by a two-sided empirical P (deviation from the null median); the six primary tests (PC1/PC2 x DI/recombination/sorting) are BH-FDR corrected. Units are never resampled independently.

## Results

| test | axis | observed | null median | null 95% | p_emp | p_adj |
|---|---|---|---|---|---|---|
| DI (Spearman rho) | PC1 | 0.100 | 0.029 | [-0.018, 0.072] | 0.0036 | 0.00432 |
| DI (Spearman rho) | PC2 | -0.131 | 0.029 | [-0.018, 0.072] | 1e-04 | 3e-04 |
| recombination (Spearman rho) | PC1 | -0.099 | 0.065 | [-0.013, 0.144] | 1e-04 | 3e-04 |
| recombination (Spearman rho) | PC2 | 0.020 | 0.065 | [-0.013, 0.144] | 0.269 | 0.269 |
| sorting, differentiated only (BF percentile gap) | PC1 | 0.036 | 0.123 | [0.078, 0.168] | 5e-04 | 0.00075 |
| sorting, differentiated only (BF percentile gap) | PC2 | 0.019 | 0.123 | [0.078, 0.168] | 2e-04 | 4e-04 |

### Supplementary statistics (not in the FDR family)

| test | axis | observed | null median | null 95% | p_emp |
|---|---|---|---|---|---|
| sorting, all units (raw-BF gap) | PC1 | 0.761 | 0.372 | [-0.097, 0.993] | 0.144 |
| sorting, all units (raw-BF gap) | PC2 | -0.418 | 0.372 | [-0.097, 0.993] | 0.0085 |
| sorting, differentiated (raw-BF gap) | PC1 | 0.640 | 0.826 | [0.364, 1.432] | 0.452 |
| sorting, differentiated (raw-BF gap) | PC2 | 0.588 | 0.826 | [0.364, 1.432] | 0.342 |
| DI (raw-BF Pearson) | PC1 | 0.102 | 0.038 | [0.001, 0.083] | 0.0042 |
| DI (raw-BF Pearson) | PC2 | -0.104 | 0.038 | [0.001, 0.083] | 1e-04 |
| recombination (raw-BF Pearson) | PC1 | -0.079 | 0.032 | [-0.047, 0.112] | 0.0062 |
| recombination (raw-BF Pearson) | PC2 | -0.027 | 0.032 | [-0.047, 0.112] | 0.146 |
| sorting magnitude (raw-BF Pearson) | PC1 | 0.070 | 0.193 | [0.123, 0.272] | 0.0025 |
| sorting magnitude (raw-BF Pearson) | PC2 | 0.115 | 0.193 | [0.123, 0.272] | 0.0433 |
| sorting magnitude / prop_fixed (Spearman rho) | PC1 | 0.051 | 0.307 | [0.231, 0.375] | 1e-04 |
| sorting magnitude / prop_fixed (Spearman rho) | PC2 | 0.119 | 0.307 | [0.231, 0.375] | 1e-04 |
| sorting orientation / uni_score (Spearman rho) | PC1 | 0.108 | -0.018 | [-0.059, 0.024] | 1e-04 |
| sorting orientation / uni_score (Spearman rho) | PC2 | -0.201 | -0.018 | [-0.059, 0.024] | 1e-04 |
| sorting, all units (BF percentile gap) | PC1 | 0.050 | 0.052 | [0.015, 0.092] | 0.904 |
| sorting, all units (BF percentile gap) | PC2 | -0.059 | 0.052 | [0.015, 0.092] | 1e-04 |

### Sensitivity to the fixation threshold (tau)

The calibration is reported over the fixation-threshold tau in {0.5, 0.6, 0.8} (the Stage-1-direct universe is fixed at n_snps>=5 throughout; there is no second minimum-cluster-size level to sweep, unlike the canonical eMLG Module C's min in {5,10}). DI and recombination do not depend on tau (shown once, at the primary tau); directional sorting is recomputed at each tau. Empirical P only (the FDR family is the six primary tests at the primary tau).

**Directional sorting (differentiated-only) across tau:**

| tau | axis | observed | null 95% | p_emp |
|---|---|---|---|---|
| 0.5 | PC1 | 0.050 | [0.071, 0.155] | 0.0033 |
| 0.5 | PC2 | 0.015 | [0.071, 0.155] | 1e-04 |
| 0.6 | PC1 | 0.036 | [0.078, 0.168] | 5e-04 |
| 0.6 | PC2 | 0.019 | [0.078, 0.168] | 2e-04 |
| 0.8 | PC1 | -0.007 | [0.112, 0.242] | 1e-04 |
| 0.8 | PC2 | 0.048 | [0.112, 0.242] | 7e-04 |

## Interpretation

**Directional sorting (primary, differentiated-only): a climate association survives FDR on both PC1 and PC2.** PC1 observed 0.036 (FDR 0.001, below the null); PC2 observed 0.019 (FDR 0.000, below the null).
**Recombination: a climate association survives FDR on PC1 only.** PC2 is within the null (observed 0.020, FDR 0.269); on PC1 the association is beyond the structured null (observed -0.099, FDR 0.000, below the null).
**Diagnostic Index: a climate association survives FDR on both PC1 and PC2.** PC1 observed 0.100 (FDR 0.004, above the null), corroborated by the raw-BF analysis (Pearson 0.102, p 0.0042); PC2 observed -0.131 (FDR 0.000, below the null), corroborated by the raw-BF analysis (Pearson -0.104, p 0.0001). (DI is a signed index; the sign is reported as-is and should not be read as locus-level adaptation in diagnostic regions.)
**Overall:** of the six primary tests, 5 survives FDR (DI (Spearman rho) x PC2; recombination (Spearman rho) x PC1; sorting, differentiated only (BF percentile gap) x PC2; sorting, differentiated only (BF percentile gap) x PC1; DI (Spearman rho) x PC1).

## What this analysis can and cannot establish

- **Can:** calibrate genome-wide, Stage-1-unit-level climate-association *patterns* against a structure- and architecture-preserving null, with the same unit universe and statistic for observed and null.
- **Cannot:** identify individual climate-associated loci beyond what the Stage-1-direct sim-FDR floor test already flags (`moduleB_stage1_S1units_null.rds`); XtX/BF share among-population allele-frequency variation, so raw-BF sensitivity variants may absorb genuine climate differentiation as well as confounding.

