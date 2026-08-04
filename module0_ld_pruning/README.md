# Module 0 — LD decay & LD-based complexity reduction

Foundation of the `_PK` pipeline. Estimates linkage-disequilibrium (LD) decay,
derives a per-marker local-LD statistic (`ld_w`), and reduces the ~1.1 M markers
to a set of **LD-reduced units** — each represented by a single marker and, for
clusters with ≥5 markers, an expected multi-locus genotype (eMLG) consensus. The
resulting clustering is the input consumed by `moduleA_sorting` (and other
downstream modules).

All scripts are run **from the repo root** (`~/gitlab/formica_hybrid`), e.g.
`Rscript module0_ld_pruning/R/module0_LD_decay_from_DIEM.R`. They use functions
from the `LDscnR` package (`devtools::load_all("~/gitlab/LDscnR/")`).

## Pipeline (`R/`)

| order | script | role |
|---|---|---|
| 1 | `module0_LD_decay_from_DIEM.R` | Parse DIEM output → oriented genotypes; chromosome-specific LD-decay fits; per-marker local-LD support `ld_w`; ROC of `ld_w`/decay-rate vs map recombination. Writes decay + track intermediates. |
| 2 | `module0_ld_pruning_DIEM.R` | Two-stage clustering: Stage 1 (connected components + complete-linkage refinement) and Stage 2 (distance-restricted, quality-gated merging). Produces the LD-reduced partition + eMLG consensus (the canonical clustering). |
| 3 | `module0_fig_ld_tracks.R` | Local-LD / decay-rate / recombination track figures for example chromosomes (Chr26, Chr10). |
| 4 | `module0_fig_fidelity_hist.R` | Consensus round-trip fidelity histogram (`score_eMLG = cor(round(x),x)²`) per ld_w flagging threshold, singleton vs merged. Reads `results_min_loci5.rds`. Writes `Figures/module0_fidelity_hist.*` and the manuscript **Fig S6** (`manuscript/figures/figS06_fidelity_hist.*`). Ported from the deprecated flat `R/fig_supplementary.R`. Lightweight, re-runnable. |

> **⚠ Cold-run generators.** Scripts 1–2 are meant to be run **once, from
> scratch**, to generate the caches and clustering; they are **not idempotent on
> warm re-runs**. `module0_LD_decay_from_DIEM.R` only creates `GTs_hybrids_005`
> inside its DIEM-parse branch, so a warm re-run (with `diem_parsed.rds` present)
> stops at the `save()` on line 192. `module0_ld_pruning_DIEM.R` has no
> `file.exists` guards and recomputes the clustering ~8× at several thresholds;
> lines 181–182 also reference an undefined `map` (should be `map_hyb_005`). The
> outputs below already exist, so downstream modules and the docs need no re-run.
> Script 3 is lightweight and re-runnable.

## Inputs

Shared, in repo-root `data/`:
- `hybrids_only_maf005.Rdata`, `hybrids_and_parents_maf005.Rdata` — genotypes (± parents)
- `Frufa_DTOL_PR.ref_genome.recmap` — genetic/recombination map
- `Formica_hybrids_filtered_diem_output.bed.gz`, `Sample_covariate_info_outlier_analysis_20.txt` — DIEM output & sample table

## Outputs

`data/` (module intermediates & the clustering):
- `eMLG_5loci_0025_cM05.rds` — **canonical clustering (0.5 cM), consumed by all downstream modules (A, B, C)**
- `eMLG_5loci_0025_cM1.rds` — 1 cM variant (sensitivity only)
- `pruned_stage1.rds`, `pruned_markers.rds`, `eMLG_groups.rds` — clustering intermediates
- `diem_parsed.rds`, `ld_decay_DIEM_100w.rds`, `ld_decay_DI.rds` — parsed genotypes & LD-decay fits
- `ld_tracks_a_windows.rds`, `ld_tracks_ldw_persnp.rds`, `ldw_a_recmap_comparison.rds` — track/comparison data
- `results_min_loci5.rds` — threshold sweep of the ≥5-marker eMLG construction

`Figures/`:
- `Chr26_stage1_vs_combined_{high,low}.png` — Stage 1 vs combined clustering (Chr26 block)
- `ld_tracks_chr26_chr10.{pdf,png}`, `ld_tracks.png`, `p_ldw_a_recmap_scatter.png`, `p_roc_low_recombination.png`
- `module0_fidelity_hist.{pdf,png}` — consensus fidelity histogram (also written to `manuscript/figures/figS06_fidelity_hist.*`)

`doc/`:
- `supplementary_methods_module0.{tex,pdf}` — Materials & Methods (LD decay + complexity reduction)

## Key parameters

Hybrid MAF ≥ 0.05 (1,114,340 markers) · LD-decay: 100 nominal windows/chr with
50% overlap, 0.95 `r²`
quantile · Stage 1: relative LD threshold ρ = 0.5 · Stage 2: `ld_w` > 0.025 or
≥5 markers, consensus fidelity ≥ 0.80, `r²` floor 0.2, 0.5 cM distance restriction (canonical; 1 cM built as a sensitivity variant)
· final partition **474,014 clusters; 32,840 with a stored consensus**. The
1-cM sensitivity partition contains 470,035 clusters, of which 32,871 have a
stored consensus.
