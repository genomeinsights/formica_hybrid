# formica_hybrid

Population genomics of 20 hybrid populations between *Formica polyctena* and *F. aquilonia* — climate adaptation and the predictability of ancestry sorting, with allopatric parental references. Haplodiploid ants (haploid males).

This README is the **data-flow map**: for each analysis, what it reads from upstream and what it writes. The *reasoning* behind the methods is in `manuscript_notes/supplementary_methods_pipeline.pdf` (the canonical methods document); this file is for navigating the code.

## Pipeline at a glance

```         
DIEM output ─► [0a] LD decay ─► hybrids_only + hybrids_and_parents (+ ld_decay, ld_w)
                                   │
                                   └► [0b] complexity reduction ─► eMLG_5loci_0025_cM05.rds  (canonical clustering)
                                                                     │
        ┌────────────────────────────────────────────────────────────┼───────────────────────────────┐
   [A] sorting phenomenon                            [B] genomic architecture        [C] climate association
   moduleA_*                                          moduleB_*                        BayPass + moduleC_*
```

All three analysis modules join on the **canonical clustering** and on the two genotype matrices. Module D (intrinsic two-locus arm) also joins on them and is built here (structure-corrected, descriptive; see below); Module E (neutral null) is a separate workstream — see *Status*.

## Key shared objects

Everything downstream reads these. `data/` is git-ignored (large, regenerable).

| file | contents | produced by |
|----|----|----|
| `data/hybrids_only_maf005.Rdata` | `GTs_hybrids_005`, `map_hyb_005` (incl. `ld_w_095`, `DiagnosticIndex`), `ld_decay`, `sample_data` | 0a |
| `data/hybrids_and_parents_maf005.Rdata` | `GTs_with_parents`, `sample_data_with_parents`, `map_hyb_005` | 0a |
| `data/eMLG_5loci_0025_cM05.rds` | `$groups`, `$eMLG`, `$pruned`, `$params` — the keystone clustering | 0b |
| `data/Frufa_DTOL_PR.ref_genome.recmap` | recombination map (cM, cM/Mb) | external |
| `{with_aland,aland_excluded,…}/*_summary_betai_reg.out` | BayPass BF(dB) per marker | BayPass runs |

**Two genotype matrices, by design:** hybrids-only for all LD estimation and clustering (including parents would let parent–hybrid structure dominate LD); hybrids+parents only where parental allele frequencies are needed (ancestry orientation, divergence).

## Stage 0 — data, LD decay, complexity reduction

**`R/LD_decay_from_DIEM.R`** - reads: `data/Formica_hybrids_filtered_diem_output.bed.gz`, `data/Sample_covariate_info_outlier_analysis_20.txt`, `data/Frufa_DTOL_PR.ref_genome.recmap` - writes: `data/diem_parsed.rds`, `data/hybrids_only_maf005.Rdata`, `data/hybrids_and_parents_maf005.Rdata`, `data/ld_decay_DIEM_100w.rds`, `data/ld_tracks_ldw_persnp.rds`, `data/ld_tracks_a_windows.rds`, `Figures/p_roc_low_recombination.png` - parses DIEM, biallelic + MAF≥0.05 filter, ancestry polarisation; fits LD decay and per-SNP `ld_w`; ROC of `ld_w`/`a` vs low recombination.

**`R/ld_pruning_DIEM.R`** - reads: `data/hybrids_only_maf005.Rdata`, `data/Frufa_DTOL_PR.ref_genome.recmap` - writes: `data/eMLG_5loci_0025_cM05.rds` (canonical), `data/pruned_markers.rds`, `data/eMLG_groups.rds` - two-stage LD complexity reduction (LDscnR) → pruning representatives + eMLG consensus genotypes.

## Module A — sorting phenomenon

**`R/moduleA_sorting_phenomenon.R`** - reads: `hybrids_and_parents_maf005.Rdata`, `hybrids_only_maf005.Rdata`, `eMLG_5loci_0025_cM05.rds`; sources `dev/R/{Ohta,parallelism_stats,eMLG_parallelism}.R` - writes: `data/moduleA_snp.rds`, `data/moduleA_clusters.rds`, `data/moduleA_dilution.rds`, `data/eMLG_sorted_cM05.rds` - per-SNP and cluster-level parallelism; gate = pooled-parental MAF ≥ 0.15.

**`R/moduleA_di_asymmetry.R`** - reads: `data/moduleA_snp.rds`, `data/moduleA_clusters.rds` - writes: `data/moduleA_di_asymmetry.rds`, `Figures/moduleA_fig1.{pdf,png}` - DI-governs-direction and cluster-size analyses.

## Module B — genomic architecture

**`R/moduleB_architecture.R`** - reads: `hybrids_and_parents_maf005.Rdata`, `hybrids_only_maf005.Rdata`, `eMLG_5loci_0025_cM05.rds`, recmap; sources `Ohta`, `parallelism_stats` - writes: `data/moduleB_architecture.rds`, `Figures/moduleB_fig2.{pdf,png}` - DI vs recombination / π / d_xy / F_ST; sorting vs recombination; direction × architecture.

**`R/moduleB_eMLG_vs_rep.R`** (validation) - reads: `moduleA_clusters.rds`, `moduleA_snp.rds`, `eMLG_5loci_0025_cM05.rds`, `eMLG_sorted_cM05.rds`, `hybrids_only_maf005.Rdata` - writes: `data/moduleB_eMLG_vs_rep.rds`, `Figures/eMLG_vs_rep_cor.png` - eMLG consensus vs representative SNP: direction robust to unit choice; consensus needed for magnitude/LD.

## Module C — climate association

BayPass inputs and runs are upstream (HPC): `R/moduleC_prepare_{with_aland,aland_excluded,sielva_excluded}.R` → `R/moduleC_write_baypass_inputs.R` → `<set>/run_baypass.sh` → `<set>/*_summary_betai_reg.out`. Covariates PC1/PC2 are per-population climate axes.

**`R/moduleC_analyse_baypass.R`** — reads clustering + BayPass `.out` → `Figures/manhattan_*.png` (outlier definition + Manhattans).

**`R/moduleC_diagnostic_index_enrichment.R`** — reads clustering, `hybrids_only`, BayPass `.out` → `data/diagnostic_index_enrichment<tag>.csv`, `Figures/diagnostic_index_enrichment_{forest,proportions}*.png` (DI-enrichment of outlier clusters).

**`R/moduleC_ancestry_confound.R`** — reads `hybrids_and_parents` → `data/moduleC_ancestry_confound.rds`, `Figures/moduleC_ancestry_confound.{pdf,png}` (PC↔ancestry confound; motivates Ω + Åland-excluded controls).

**`R/moduleC_sorting_climate.R`** — reads `hybrids_and_parents`, `hybrids_only`, clustering, `moduleA_snp.rds`, BayPass `.out`; sources shared stats → `data/moduleC_C3_cl.rds` (**consensus checkpoint, reused across threshold settings**), `data/moduleC_sorting_climate_<tag>.rds`, `Figures/moduleC_fig3_<tag>.{pdf,png}` (sorting × outlier overlap; threshold/binary version).

**`R/moduleC_rate_based.R`** (**primary Module C analysis**) — reads `data/moduleC_C3_cl.rds`, clustering, `hybrids_and_parents`, BayPass `.out`, `data/diagnostic_index_enrichment<tag>.csv` → `data/moduleC_rate_based_<tag>.rds`, `Figures/moduleC_dose_response_<tag>.{pdf,png}` — size-normalised, cluster-level enrichment (replaces the size-gated outlier count).

## Module D — intrinsic (two-locus incompatibility) arm

Screens **unlinked** LD-reduced unit pairs for residual associations beyond shared ancestry and relatedness — candidate two-locus (Dobzhansky–Muller) incompatibilities, complementary to the per-locus sorting of A/B and the climate axis of C. *Unlinked* = different chromosome, or same chromosome \> 10 cM on the genetic map (≈99% admixture-LD decay in \~50 generations). One scan → four filters → an annotation → modules → the neutral null (Module E). Minimal-pipeline spec: `manuscript_notes/moduleD_plan.md`.

**`dev/R/moduleD_emmax.R`** — reads `eMLG_5loci_0025_cM05.rds`, `hybrids_only_maf005.Rdata` (`sample_data`), `moduleC_C3_cl.rds` (differentiated / DI gate), recmap; sources `dev/R/{emmax,moduleD_paralogy}.R` — writes `data/moduleD_emmax.rds`, `Figures/moduleD_emmax_manhattans.pdf` — the structure-corrected two-locus scan: each differentiated eMLG dosage is a trait tested against every other in an EMMAX LMM with a **double-LOCO** VanRaden GRM built from differentiated units, each focal **conditioned on the top 10 genome PCs**; associations signed cis/trans by the structure-corrected GLS coefficient; calibrated (λ ≈ 0.98–1.05). Focal set is targeted (bidirectional + extreme-covariance + random controls).

**`dev/R/moduleD_paralogy.R`** (shared filter) — `flag_paralogy()`: median within-population \|r\| \> 0.9 flags cross-chromosome assembly duplicates that a genetic-distance rule cannot catch; stores genotype concordance + per-unit excess Ho as corroboration.

**`dev/R/moduleD_network_build.R`** — reads `moduleD_emmax.rds`, clustering, `moduleC_C3_cl.rds`, `hybrids_only_maf005.Rdata`, `moduleD_bidirectional.rds`, recmap; sources `dev/R/moduleD_paralogy.R` — writes `data/moduleD_network.rds` (meta-nodes + meta-edges) — the single reproducible read-out: global Benjamini–Hochberg FDR (`q < 0.01`) over all unlinked tests → paralogy filter → third-level within-chromosome **average-linkage** merge (`|r| > 0.5`, ≤ 10 cM; single linkage chains, so avoided) → single-population **leverage** filter (leave-one-population-out `|r| ≥ 0.3`, which subsumes any near-fixed/MAF cutoff) → low-recombination **structure annotation** (recombination percentile \< 0.1 — a label carried into Module E's null, *never* a filter) → cross-chromosome correlation **modules** (`|r| > 0.4`).

**`dev/R/moduleD_module_heatmaps.R`** — reads `moduleD_network.rds`, clustering, `moduleC_C3_cl.rds`, `hybrids_only_maf005.Rdata`, recmap — writes `Figures/moduleD_module_<id>_heatmap.png` — per-module genotype heatmaps (each module = one coherent multi-chromosome co-ancestry block).

**`dev/R/moduleD_network_circos.R`** — reads `moduleD_network.rds`, clustering, `hybrids_only_maf005.Rdata`, `moduleC_C3_cl.rds`, `moduleD_bidirectional.rds`, recmap — writes `Figures/moduleD_{trans,candidate,structure,module}_circos.png` + `moduleD_trans_network.png` — the network and its region-band views (combined, candidate-only, structure-only, within-module).

**`dev/R/moduleD_bidirectional.R`** (annotation only, not a DMI screen) — reads `moduleA_snp.rds`, clustering — writes `data/moduleD_bidirectional.rds` (`reps`) — the set of bidirectionally sorting units for the circos/heatmap sorting-ring annotation.

**`dev/R/moduleD_ohta_dmi.R`** — retained only as **Module E's measurement instrument** (among-population Ohta LD; reusable `moduleD_pop_freq_matrix()` / `moduleD_prefilter()` / `moduleD_scan()`), not a standalone screen. The intrinsic call (excess over neutral, judged at each pair's *local* recombination rate) is made only against Module E's recombination-matched null.

## Shared code and helpers

- `dev/R/Ohta.R` — `ohta_fast_prepare()` (per-population allele-frequency prep).
- `dev/R/parallelism_stats.R` — the core sorting statistic (`parallelism_stats()`).
- `dev/R/eMLG_parallelism.R` — `build_sorted_eMLG()`, `build_group_consensus()`, `cluster_DI()`.
- `R/fig_ld_tracks.R` — reads `data/ld_tracks_{a_windows,ldw_persnp}.rds` + recmap → `Figures/ld_tracks_chr26_chr10.{pdf,png}`.
- `LDscnR` package (`~/gitlab/LDscnR`) — LD decay, complexity reduction, consensus construction.

## Conventions

- **Version-tagged clustering.** The canonical clustering filename encodes its parameters (`eMLG_5loci_0025_cM05.rds`); never overwrite a fixed name — group IDs are not stable across clusterings.
- **Parameter-tagged outputs.** Threshold-dependent outputs carry the settings in the filename (e.g. `_5_15`), so alternative settings coexist.
- **Terminology.** The units are *LD clusters*, not haplotype blocks: a biallelic SNP carries one bipartition, so markers on different genealogical branches form different clusters that can overlap in the same interval.
- **Units.** eMLG consensus for magnitude and LD-based tests; pruning representatives are adequate for direction (see `moduleB_eMLG_vs_rep.R`).
- `data/` and `Figures/` are git-ignored (regenerable); the canonical figures are committed selectively.

## Documentation

- `manuscript_notes/supplementary_methods_pipeline.{tex,pdf}` — **canonical methods** (LD decay → clustering → A/B/C), with the parameter table.
- `manuscript_notes/module{B,C}_results_summary.{tex,pdf}`, `moduleA_results_summary.md` — per-module results.
- `manuscript/supplementary_methods_ld_module_D.{tex,pdf}`, `supplementary_results_module_D.{tex,pdf}` — Module D methods + results (manuscript style); `manuscript_notes/moduleD_plan.md` — the minimal-pipeline spec; `manuscript_notes/moduleD_wip_summary.{tex,pdf}` — working notes.
- `dev/methods_notes.md` — LD-pruning / eMLG design rationale.
- `dev/HANDOFF_SUMMARY_*.md` — historical thread handoffs.

## Status

- **Module D** (structure-corrected two-locus screen for intrinsic incompatibilities; `dev/R/moduleD_{emmax,network_build,…}.R`) — **built, descriptive**. Largely negative: the scan's strongest signal is technical (paralogy) or a low-recombination co-ancestry component (resolving into several distinct multi-chromosome modules; nature unresolved pending E), not pair-specific. The residue is an aquilonia co-sorting module (belongs with A/C) and a single ancestry-independent trans candidate (F11431–F49480). Low recombination is annotated, not filtered; the intrinsic call (excess over neutral, judged at each pair's *local* recombination rate) plugs in when Module E's null is ready.
- **Module E** (recombination-matched haplodiploid neutral null; `dev/R/moduleE_*.R`) — separate workstream. The inference license: until it exists, the descriptive sorting results are "consistent with neutral" only.
