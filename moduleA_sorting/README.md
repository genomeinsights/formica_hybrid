# Module A — Ancestry sorting & genomic architecture

Quantifies the **extent and parental direction of ancestry sorting** across the
20 hybrid populations, and its **genomic architecture** (relation to diagnostic
index, recombination, parental summaries, cluster size). Works at two
resolutions — per SNP and per LD-reduced unit — using the clustering produced by
`module0_ld_pruning`.

All scripts run **from the repo root** (`~/gitlab/formica_hybrid`), e.g.
`Rscript moduleA_sorting/R/moduleA_sorting_phenomenon.R`. Expensive steps are
**cached** (read back if present; a per-file parameter stamp rebuilds a stale
cache, and each script has a `RECOMPUTE <- TRUE` switch to force a rebuild).

## Conventions (locked)

Parental-MAF gate ≥ 0.15 · sorting threshold **τ = 0.5** · near-fixation floor
**φ = 0.90** (`fix_th = 0.1`, on the major/fixed-for allele) · **cM1** clustering
· diagnostic index (DI) kept **ungated** as a covariate.

## Pipeline (`R/`)

Shared statistical libraries (sourced, no I/O): `Ohta.R` (`ohta_fast_prepare`),
`parallelism_stats.R` (`parallelism_stats`), `eMLG_parallelism.R`
(`build_sorted_eMLG`, `build_group_consensus`, `cluster_DI`).

| order | script | role |
|---|---|---|
| 1 | `moduleA_sorting_phenomenon.R` | Per-SNP parallelism → sorted markers (A1); companion eMLG for the sorted clusters (A2); one cluster-level test + `score_eMLG` dilution check (A3); the τ×φ threshold-sensitivity sweeps (sorting classification A1b, DI-decile direction A1c). Produces `moduleA_snp.rds`, used by the other two scripts. |
| 2 | `moduleA_architecture.R` | Reproduces the genomic-architecture results (the old manuscript "Module B"): DI–recombination/π/dxy correlations + decile table; sorting depleted in low recombination + magnitude slope; the `P(aquilonia)` direction logistic model + its threshold sweep; assembles the 3-panel architecture figure. |
| 3 | `moduleA_fig_sorting_manhattan.R` | Genome-wide sorting landscape (toward-aquilonia, toward-polyctena, net direction, DI, and marker-level `ld_w`), across thresholds. Pure post-hoc on `moduleA_snp.rds`. |

## Inputs

Shared, repo-root `data/`: `hybrids_and_parents_maf005.Rdata`,
`hybrids_only_maf005.Rdata`, `Frufa_DTOL_PR.ref_genome.recmap`.
From Module 0: `module0_ld_pruning/data/eMLG_5loci_0025_cM1.rds` (canonical clustering).

## Outputs

`data/`:
- `moduleA_snp.rds` — per-SNP `parallelism_stats` output (the hub file)
- `eMLG_sorted_cM1.rds`, `eMLG_sorted_cM1_parents.rds` — A2/A3 companion eMLG (hybrid + matched parent side)
- `moduleA_clusters.rds`, `moduleA_dilution.rds` — cluster-level test & dilution check
- `moduleA_sortth_fixth_sweep.rds`, `moduleA_panelB_sweep.rds` — threshold-sensitivity sweeps
- `moduleA_architecture.rds` — architecture correlations, decile tables, models, direction sweep
- `moduleA_prep_snp.rds`, `moduleA_r_unit.rds`, `moduleA_r_snp.rds` — caches

`Figures/`:
- `moduleA_sorting_sweep.png`, `moduleA_panelB_sweep.png` — sweep supplements S1, S2
- `moduleA_direction_sweep.png` — direction-model coefficients vs τ (supplement S3)
- `moduleA_sorting_manhattan{,_directional}.png` — genome-wide landscape (S4)
- `moduleA_architecture_fig.{pdf,png}` — 3-panel architecture figure (Fig 1)

`doc/`:
- `supplementary_methods_moduleA.{tex,pdf}` — Materials & Methods (sorting + architecture)
- `moduleA_architecture_summary.{tex,pdf}` — architecture results summary
- `moduleA_supp_th_sensitivity.{tex,pdf}` — threshold-sensitivity supplementary figures (S1–S4)

## Not yet in this module

The per-SNP direction "stem-plot" exploration (`R_PK/moduleA_fig_snp_direction*.R`,
`Figures_PK/moduleA_snp_direction*.png`) is held out pending revisit; those
scripts still point at the pre-reorganisation `RDS_data/moduleA_snp.rds` path.
