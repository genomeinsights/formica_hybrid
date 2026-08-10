------------------------------------------------------------------------

editor_options: markdown: wrap: 72 ---

# Module A — Ancestry sorting & genomic architecture

Quantifies the **extent and parental direction of ancestry sorting** across the 20 hybrid populations, and its **genomic architecture** (relation to diagnostic index, recombination, parental summaries, cluster size). Works at two resolutions — per SNP and per LD-reduced unit — using the clustering produced by `module0_ld_pruning`.

All scripts run **from the repo root** (`~/gitlab/formica_hybrid`), e.g. `Rscript moduleA_sorting/R/moduleA_sorting_phenomenon.R`. Expensive steps are **cached** (read back if present; a per-file parameter stamp rebuilds a stale cache, and each script has a `RECOMPUTE <- TRUE` switch to force a rebuild).

## Conventions (locked)

Parental-MAF gate ≥ 0.15 · sorting threshold **τ = 0.6** · near-fixation floor **φ = 0.85** (`fix_th = 0.15`, on the major/fixed-for allele) · **cM05** clustering · diagnostic index (DI) kept **ungated** as a covariate · sorting call `sort_rule = "binom"` (α = 0.05): magnitude gate = total near-fixation (`prop_fixed ≥ τ`), with direction resolved only when the split is significantly biased toward one parent. A non-significant direction is reported as **direction unresolved** (or ambiguous when too few populations are fixed to test); this is not positive evidence of equal sorting toward both parents. The earlier component rule is retained only for provenance.

## Pipeline (`R/`)

`ohta_fast_prepare()` (and the rest of the Ohta D-statistic suite) now comes from the **`LDscnR` package** (`devtools::load_all("~/gitlab/LDscnR/")`), not a local file. Two reviewed stat libraries are still sourced locally (no I/O): `parallelism_stats.R` (`parallelism_stats`), `eMLG_parallelism.R` (`build_sorted_eMLG`, `build_group_consensus`, `cluster_DI`).

| order | script | role |
|----|----|----|
| 1 | `moduleA_sorting_phenomenon.R` | Per-SNP parallelism → sorted markers (A1); companion eMLG for the sorted clusters (A2); one cluster-level test + `score_eMLG` dilution check (A3); the τ×φ threshold-sensitivity sweeps (sorting classification A1b, DI-decile direction A1c). Produces `moduleA_snp.rds`, used by the other two scripts. |
| 2 | `moduleA_architecture.R` | Reproduces the genomic-architecture results (the old manuscript "Module B"): DI–recombination/π/dxy correlations + decile table; sorting depleted in low recombination + magnitude slope; the `P(aquilonia)` direction logistic model + its threshold sweep; assembles the 3-panel architecture figure; and the within-cluster DI-spread vs cluster-size supplementary panel that justifies `di_agg = "max"`. |
| 3 | `moduleA_fig_sorting_manhattan.R` | Genome-wide sorting landscape (toward-aquilonia, toward-polyctena, net direction, DI, and marker-level `ld_w`), across thresholds. Pure post-hoc on `moduleA_snp.rds`. |
| 4 | `moduleA_cluster_sorting.R` | The per-eMLG sorting call over the **full has_eMLG universe** (`colnames($eMLG)`) with the locked conventions — the annotation the climate modules (B, C) join on. Formerly the flat `R/moduleC_sorting_climate.R` checkpoint ("C3_cl"); relocated here because it is a sorting product. Own full-universe parent-side consensus, so it is independent of A3's sorted-focused companion set. |
| 5 | `moduleA_direction_unresolved_peaks.R` | **Diagnostic** (legacy filename; post-hoc on `moduleA_snp.rds`): tests whether the per-SNP direction-unresolved Manhattan peaks at τ=0.5 are LD-cluster artifacts — per-SNP vs LD-reduced recount, `n_fixed` tabulation, and the per-peak population split. They collapse (single high-`ld_w` blocks; 84% have `n_fixed` 10–11). |
| 6 | `moduleA_tau_sensitivity.R` | **Diagnostic only**: compares stored τ-independent counts across {0.5, 0.6, 0.8}. The reported analysis remains τ=0.6; direction-unresolved calls increase at relaxed τ partly because directional power is lower. |

## Inputs

Shared, repo-root `data/`: `hybrids_and_parents_maf005.Rdata`, `hybrids_only_maf005.Rdata`, `Frufa_DTOL_PR.ref_genome.recmap`. From Module 0: `module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds` (canonical clustering).

## Outputs

`data/`: - `moduleA_snp.rds` — per-SNP `parallelism_stats` output (the hub file) - `eMLG_sorted_cM05.rds`, `eMLG_sorted_cM05_parents.rds` — A2/A3 companion eMLG (hybrid + matched parent side) - `moduleA_clusters.rds`, `moduleA_dilution.rds` — cluster-level test & dilution check - `moduleA_cluster_sorting.rds` — **has_eMLG sorting annotation consumed by `moduleB_climate_GEA` / `moduleC_climate_vs_sorting`** (`group_id`, `n_loci`, `differentiated`, `sort_class`, `DI`, `prop_fixed`, `uni_score`, `directional`, `sorted`) - `moduleA_sortth_fixth_sweep.rds`, `moduleA_panelB_sweep.rds` — threshold-sensitivity sweeps - `moduleA_architecture.rds` — architecture correlations, decile tables, models, direction sweep, within-cluster DI-spread table (`di_size_tab`, `di_sd_size_spearman`) - `moduleA_prep_snp.rds`, `moduleA_r_unit.rds`, `moduleA_r_snp.rds` — caches

`Figures/`: - `moduleA_sorting_sweep.png`, `moduleA_panelB_sweep.png` — sweep supplements S1, S2 - `moduleA_direction_sweep.png` — direction-model coefficients vs τ (supplement S3) - `moduleA_sorting_manhattan{,_directional}.png` — genome-wide landscape (S4) - `moduleA_architecture_fig.{pdf,png}` — 3-panel architecture figure (Fig 1) - `moduleA_di_variance.{pdf,png}` — within-cluster DI spread vs LD-cluster size (architecture supplementary panel)

`doc/`: - `supplementary_methods_moduleA.{tex,pdf}` — Materials & Methods (sorting + architecture) - `moduleA_architecture_summary.{tex,pdf}` — architecture results summary (incl. the within-cluster DI-spread supplementary panel) - `moduleA_supp_th_sensitivity.{tex,pdf}` — threshold-sensitivity supplementary figures (S1–S4) - `moduleA_tau_sensitivity.md` — diagnostic notes for the τ sensitivity series. Files or objects retaining “bidirectional” in their names are legacy identifiers; manuscript-facing text uses “direction unresolved.”

## Not yet in this module

The per-SNP direction "stem-plot" exploration (`R_PK/moduleA_fig_snp_direction*.R`, `Figures_PK/moduleA_snp_direction*.png`) is held out pending revisit; those scripts still point at the pre-reorganisation `RDS_data/moduleA_snp.rds` path.
