# Module C, Stage-1-direct: genome-wide, unit-level bio_winter calibration

*Generated 2026-09-13. NSIM = 10000 Omega-structured null covariates (REUSED from the PC1/PC2 Module C run -- see Data provenance); Stage-1-direct unit universe = 18361 clusters (n_snps>=5, Aland excluded, 19 pops, Stage-1-derived Omega).*

## Data provenance

- **Scope:** the Stage-1 LARGE-CLUSTER scan (18,361 clusters with >=5 SNPs), identical universe to the canonical PC1/PC2 Module C run (`moduleC_stage1_report.md`).
- **Observed bio_winter association:** per-unit BayPass BF(dB) on the bio_winter covariate (`module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/bio_winter_S1units_withOmega_summary_betai_reg.out`) -- same withOmega Stage-1-direct BayPass setup as PC1/PC2 (identical `.geno`/`omega_mat_omega.out`/poolsize inputs, same continuous-covariate regression model).
- **Null covariates: REUSED, not regenerated.** The 10,000 Omega-structured null covariates and their reduced genome-wide statistics come from `moduleC_stage1_null_stats.rds` (the PC1/PC2 Module C run). That null distribution is a property of random-covariate-BF-vs-annotation relationships under this Omega/genotype/unit setup -- it does not depend on which real covariate (PC1, PC2, or bio_winter) is being tested, so no new BayPass/mini2 run was needed here. Validity of this reuse is enforced by an exact md5 fingerprint match between this run's annotations and the annotations recorded when the null was built (checked below, hard stop on mismatch).
- **Annotations (per-unit, joined by `group_id`):** identical to the PC1/PC2 run -- consensus best-SNP Diagnostic Index, unidirectional-sorting status, sorting magnitude `prop_fixed` and signed orientation `uni_score`, recombination rate (cM/Mb) at the best-SNP marker, cluster size.

## Validation checks

- Annotation fingerprint used here == fingerprint recorded when the reused null was built (md5 `0fa3b167ed2da567464d500bd66d6ed4`); the stopifnot in this script hard-stops if they ever diverge.
- bio_winter BF vector: N = 18361, order verified identical to the Stage-1 BayPass row order (`S1units_group_order.txt`), all finite.

## Methods

Identical to the PC1/PC2 Module C analysis (see `moduleC_stage1_report.md`): the genome-wide bio_winter BF vector is reduced to Spearman rho with DI, Spearman rho with recombination, and the directional-sorting percentile gap among differentiated units (primary), plus supplementary statistics. The observed reduction is compared against the (reused) 10,000-null distribution by a two-sided empirical P; the three primary tests (bio_winter x DI/recombination/sorting) are BH-FDR corrected.

## Results

| test | axis | observed | null median | null 95% | p_emp | p_adj |
|---|---|---|---|---|---|---|
| DI (Spearman rho) | bio_winter | -0.055 | 0.029 | [-0.018, 0.072] | 0.001 | 0.0015 |
| recombination (Spearman rho) | bio_winter | 0.142 | 0.065 | [-0.013, 0.144] | 0.0556 | 0.0556 |
| sorting, differentiated only (BF percentile gap) | bio_winter | 0.036 | 0.123 | [0.078, 0.168] | 5e-04 | 0.0015 |

### Supplementary statistics (not in the FDR family)

| test | axis | observed | null median | null 95% | p_emp |
|---|---|---|---|---|---|
| sorting, all units (raw-BF gap) | bio_winter | -0.307 | 0.372 | [-0.097, 0.993] | 0.0202 |
| sorting, differentiated (raw-BF gap) | bio_winter | 0.545 | 0.826 | [0.364, 1.432] | 0.268 |
| DI (raw-BF Pearson) | bio_winter | -0.042 | 0.038 | [0.001, 0.083] | 2e-04 |
| recombination (raw-BF Pearson) | bio_winter | 0.105 | 0.032 | [-0.047, 0.112] | 0.0748 |
| sorting magnitude (raw-BF Pearson) | bio_winter | 0.124 | 0.193 | [0.123, 0.272] | 0.0699 |
| sorting magnitude / prop_fixed (Spearman rho) | bio_winter | 0.112 | 0.307 | [0.231, 0.375] | 1e-04 |
| sorting orientation / uni_score (Spearman rho) | bio_winter | -0.091 | -0.018 | [-0.059, 0.024] | 8e-04 |
| sorting, all units (BF percentile gap) | bio_winter | -0.016 | 0.052 | [0.015, 0.092] | 8e-04 |

### Sensitivity to the fixation threshold (tau)

Reported over the fixation-threshold tau in {0.5, 0.6, 0.8} (Stage-1-direct universe fixed at n_snps>=5). DI and recombination do not depend on tau (shown once, at the primary tau); directional sorting is recomputed at each tau. Empirical P only (the FDR family is the three primary tests at the primary tau).

**Directional sorting (differentiated-only) across tau:**

| tau | axis | observed | null 95% | p_emp |
|---|---|---|---|---|
| 0.5 | bio_winter | 0.010 | [0.071, 0.155] | 1e-04 |
| 0.6 | bio_winter | 0.036 | [0.078, 0.168] | 5e-04 |
| 0.8 | bio_winter | 0.126 | [0.112, 0.242] | 0.167 |

## Interpretation

**Directional sorting (primary, differentiated-only): a bio_winter association survives FDR.** Observed 0.036 (FDR 0.001, below the null).
**Recombination: no bio_winter association survives FDR.** Observed 0.142 (FDR 0.056), within the structured-null 95% interval.
**Diagnostic Index: a bio_winter association survives FDR.** Observed -0.055 (FDR 0.001, below the null), corroborated by the raw-BF analysis (Pearson -0.042, p 0.0002). (DI is a signed index; the sign is reported as-is and should not be read as locus-level adaptation in diagnostic regions.)
**Overall:** of the three primary tests, 2 survives FDR (DI (Spearman rho); sorting, differentiated only (BF percentile gap)).

## What this analysis can and cannot establish

- **Can:** calibrate genome-wide, Stage-1-unit-level bio_winter-association *patterns* against a structure- and architecture-preserving null, with the same unit universe and statistic as the PC1/PC2 run, at no extra BayPass cost.
- **Cannot:** identify individual bio_winter-associated loci beyond what the Stage-1-direct sim-FDR floor test already flags (`moduleB_stage1_bio_winter_null.rds`); XtX/BF share among-population allele-frequency variation, so raw-BF sensitivity variants may absorb genuine bio_winter differentiation as well as confounding.

