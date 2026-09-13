# Module C, Stage-1-direct: genome-wide, unit-level mitoC2 calibration

*Generated 2026-09-13. NSIM = 10000 structure-matched null population contrasts; Stage-1-direct unit universe = 18361 clusters (n_snps>=5, Aland excluded, 19 pops, Stage-1-derived Omega).*

## Data provenance

- **Scope:** the Stage-1 LARGE-CLUSTER scan (18,361 clusters with >=5 SNPs), identical universe to the PC1/PC2/bio_winter Module C runs.
- **Observed mitoC2 association:** per-unit BayPass log10(1/pval) from the C2 population-contrast test (`mito_C2_S1units_summary_contrast.out`), contrasting 7 Faquilonia-like vs 12 Fpolyctena-like populations (19 populations total, Aland excluded) -- same withOmega Stage-1-direct BayPass setup as PC1/PC2/bio_winter (`-omegafile` pointed at the identical frozen Stage-1 Omega), but in `-contrastfile` mode rather than `-efile` covariate-regression mode.
- **Null contrasts: a DEDICATED contrast-mode null, NOT the PC1/PC2/bio_winter null.** A continuous-covariate null (random Omega-structured numeric draws through `-efile`) is not a matched reference for a binary population-contrast statistic, so this reuses the null already built by `moduleB_stage1_mitoC2_null.R` for the Stage-1 sim-FDR floor test: each of 10,000 null contrasts is a population-level +1/-1 partition with the SAME group sizes as the real split (7 vs 12), obtained by rank-thresholding the same Omega-eigenvector draws used for the climate null (so populations with high covariance still tend to land on the same side of the null partition, preserving the structure-matching property), then rerun through BayPass in `-contrastfile` mode with the identical MCMC seed (74) and Stage-1 Omega as the real mitoC2 run.
- **No new BayPass run was needed for THIS Module C reduction.** Unlike the PC1/PC2 null (which only kept exceedance counts, forcing a ~9-10h BayPass rerun for `moduleC_stage1_null_regen.R`), `moduleB_stage1_mitoC2_null.R` persisted every batch's full null log10(1/pval) matrix (`null/bf_matrices/mitoC2_bf_b##.rds`, 50 files). Those were pulled back from mini2 and reduced here directly.
- **Exact integrity check (not Monte-Carlo tolerance):** because these are the literal persisted matrices already used for the floor test (not a fresh MCMC realization), the per-unit null-exceedance count recomputed here matches `moduleB_stage1_mitoC2_null.rds`'s saved `k3` EXACTLY for all 18361 units (hard stop otherwise, enforced in the script).
- **Annotations (per-unit, joined by `group_id`):** identical to the PC1/PC2/bio_winter runs -- consensus best-SNP Diagnostic Index, unidirectional-sorting status, sorting magnitude `prop_fixed` and signed orientation `uni_score`, recombination rate (cM/Mb) at the best-SNP marker, cluster size.

## Validation checks

- All 50 persisted null matrices present and correctly shaped (18361 x 200), all finite.
- Exact k3 exceedance-count match against the floor-test run (see above).
- mitoC2 log10(1/pval) vector: N = 18361, order verified identical to the Stage-1 BayPass row order, all finite.

## Methods

Identical statistical framework to the PC1/PC2/bio_winter Module C analyses, with log10(1/pval) substituted for BF(dB) as the per-unit association-strength input (monotonic with C2 within each null draw, so the rank-based primary statistics are unaffected by this substitution): the genome-wide mitoC2 log10(1/pval) vector is reduced to Spearman rho with DI, Spearman rho with recombination, and the directional-sorting percentile-gap among differentiated units (primary), plus supplementary statistics. The observed reduction is compared against the dedicated 10,000-null-contrast distribution by a two-sided empirical P; the three primary tests (mitoC2 x DI/recombination/sorting) are BH-FDR corrected.

## Results

| test | axis | observed | null median | null 95% | p_emp | p_adj |
|---|---|---|---|---|---|---|
| DI (Spearman rho) | mitoC2 | 0.015 | 0.021 | [-0.037, 0.069] | 0.821 | 0.963 |
| recombination (Spearman rho) | mitoC2 | 0.006 | 0.008 | [-0.082, 0.114] | 0.963 | 0.963 |
| sorting, differentiated only (log10p percentile gap) | mitoC2 | 0.037 | 0.016 | [-0.030, 0.072] | 0.428 | 0.963 |

### Supplementary statistics (not in the FDR family)

| test | axis | observed | null median | null 95% | p_emp |
|---|---|---|---|---|---|
| sorting, all units (raw log10(1/pval) gap) | mitoC2 | 0.033 | 0.028 | [-0.076, 0.176] | 0.917 |
| sorting, differentiated (raw log10(1/pval) gap) | mitoC2 | 0.060 | 0.019 | [-0.051, 0.141] | 0.376 |
| DI (raw log10(1/pval) Pearson) | mitoC2 | 0.031 | 0.033 | [-0.020, 0.082] | 0.94 |
| recombination (raw log10(1/pval) Pearson) | mitoC2 | -0.015 | 0.004 | [-0.075, 0.103] | 0.644 |
| sorting magnitude (raw log10(1/pval) Pearson) | mitoC2 | 0.050 | 0.031 | [-0.035, 0.169] | 0.626 |
| sorting magnitude / prop_fixed (Spearman rho) | mitoC2 | 0.062 | 0.032 | [-0.044, 0.159] | 0.541 |
| sorting orientation / uni_score (Spearman rho) | mitoC2 | 0.037 | -0.005 | [-0.067, 0.048] | 0.177 |
| sorting, all units (log10p percentile gap) | mitoC2 | 0.012 | 0.017 | [-0.049, 0.095] | 0.858 |

### Sensitivity to the fixation threshold (tau)

Reported over the fixation-threshold tau in {0.5, 0.6, 0.8} (Stage-1-direct universe fixed at n_snps>=5). DI and recombination do not depend on tau (shown once, at the primary tau); directional sorting is recomputed at each tau. Empirical P only (the FDR family is the three primary tests at the primary tau).

**Directional sorting (differentiated-only) across tau:**

| tau | axis | observed | null 95% | p_emp |
|---|---|---|---|---|
| 0.5 | mitoC2 | 0.044 | [-0.019, 0.065] | 0.308 |
| 0.6 | mitoC2 | 0.037 | [-0.030, 0.072] | 0.428 |
| 0.8 | mitoC2 | 0.026 | [-0.092, 0.118] | 0.617 |

## Interpretation

**Directional sorting (primary, differentiated-only): no mitoC2 association survives FDR.** Observed 0.037 (FDR 0.963), within the structured-null 95% interval. Sorting magnitude (`prop_fixed`) is likewise null (supplementary, p_emp 0.54).
**Recombination: no mitoC2 association survives FDR.** Observed 0.006 (FDR 0.963), within the structured-null 95% interval.
**Diagnostic Index: no mitoC2 association survives FDR.** Observed 0.015 (FDR 0.963), within the structured-null 95% interval. (DI is a signed index; the sign is reported as-is and should not be read as locus-level adaptation in diagnostic regions.)
**Overall:** no primary test is exceptional; mitoC2 association evidence is not concentrated in diagnostic, directionally-sorted, or low-recombination Stage-1 units beyond population structure and genomic architecture.

## What this analysis can and cannot establish

- **Can:** calibrate genome-wide, Stage-1-unit-level mitoC2-association *patterns* against a null built specifically for the contrast-mode statistic (matched group sizes, same Omega-eigenvector structure-preservation as the continuous-covariate null), with the same unit universe and annotation set as the PC1/PC2/bio_winter runs.
- **Cannot:** identify individual mitoC2-associated loci beyond what the Stage-1-direct sim-FDR floor test already flags (`moduleB_stage1_mitoC2_null.rds`); be compared directly, panel-for-panel, against the PC1/PC2/bio_winter calibration figure on ONE set of null histograms -- the null distributions are built from genuinely different processes (continuous Omega-eigenvector draws vs rank-thresholded +1/-1 partitions of those same draws) and are reported on separate figures for that reason, even though the primary statistics themselves (Spearman rho, percentile gap) are on a comparable scale.

