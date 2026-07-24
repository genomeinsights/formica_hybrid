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

All three analysis modules join on the **canonical clustering** and on the two genotype matrices. Modules D (Ohta DMI test) and E (neutral null) are not built here — see *Status* below.

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

**`dev/R/moduleA_sorting_phenomenon.R`** - reads: `hybrids_and_parents_maf005.Rdata`, `hybrids_only_maf005.Rdata`, `eMLG_5loci_0025_cM05.rds`; sources `dev/R/{Ohta,parallelism_stats,eMLG_parallelism}.R` - writes: `data/moduleA_snp.rds`, `data/moduleA_clusters.rds`, `data/moduleA_dilution.rds`, `data/eMLG_sorted_cM05.rds` - per-SNP and cluster-level parallelism; gate = pooled-parental MAF ≥ 0.15.

**`dev/R/moduleA_di_asymmetry.R`** - reads: `data/moduleA_snp.rds`, `data/moduleA_clusters.rds` - writes: `data/moduleA_di_asymmetry.rds`, `Figures/moduleA_fig1.{pdf,png}` - DI-governs-direction and cluster-size analyses.

## Module B — genomic architecture

**`dev/R/moduleB_architecture.R`** - reads: `hybrids_and_parents_maf005.Rdata`, `hybrids_only_maf005.Rdata`, `eMLG_5loci_0025_cM05.rds`, recmap; sources `Ohta`, `parallelism_stats` - writes: `data/moduleB_architecture.rds`, `Figures/moduleB_fig2.{pdf,png}` - DI vs recombination / π / d_xy / F_ST; sorting vs recombination; direction × architecture.

**`dev/R/moduleB_eMLG_vs_rep.R`** (validation) - reads: `moduleA_clusters.rds`, `moduleA_snp.rds`, `eMLG_5loci_0025_cM05.rds`, `eMLG_sorted_cM05.rds`, `hybrids_only_maf005.Rdata` - writes: `data/moduleB_eMLG_vs_rep.rds`, `Figures/eMLG_vs_rep_cor.png` - eMLG consensus vs representative SNP: direction robust to unit choice; consensus needed for magnitude/LD.

## Module C — climate association

BayPass inputs and runs are upstream (HPC): `R/prepare_{with_aland,aland_excluded,sielva_excluded}.R` → `R/write_baypass_inputs.R` → `<set>/run_baypass.sh` → `<set>/*_summary_betai_reg.out`. Covariates PC1/PC2 are per-population climate axes.

**`R/analyse_baypass.R`** — reads clustering + BayPass `.out` → `Figures/manhattan_*.png` (outlier definition + Manhattans).

**`R/diagnostic_index_enrichment.R`** — reads clustering, `hybrids_only`, BayPass `.out` → `data/diagnostic_index_enrichment<tag>.csv`, `Figures/diagnostic_index_enrichment_{forest,proportions}*.png` (DI-enrichment of outlier clusters).

**`dev/R/moduleC_ancestry_confound.R`** — reads `hybrids_and_parents` → `data/moduleC_ancestry_confound.rds`, `Figures/moduleC_ancestry_confound.{pdf,png}` (PC↔ancestry confound; motivates Ω + Åland-excluded controls).

**`dev/R/moduleC_sorting_climate.R`** — reads `hybrids_and_parents`, `hybrids_only`, clustering, `moduleA_snp.rds`, BayPass `.out`; sources shared stats → `data/moduleC_C3_cl.rds` (**consensus checkpoint, reused across threshold settings**), `data/moduleC_sorting_climate_<tag>.rds`, `Figures/moduleC_fig3_<tag>.{pdf,png}` (sorting × outlier overlap; threshold/binary version).

**`dev/R/moduleC_rate_based.R`** (**primary Module C analysis**) — reads `data/moduleC_C3_cl.rds`, clustering, `hybrids_and_parents`, BayPass `.out`, `data/diagnostic_index_enrichment<tag>.csv` → `data/moduleC_rate_based_<tag>.rds`, `Figures/moduleC_dose_response_<tag>.{pdf,png}` — size-normalised, cluster-level enrichment (replaces the size-gated outlier count).

## Shared code and helpers

- `dev/R/Ohta.R` — `ohta_fast_prepare()` (per-population allele-frequency prep).
- `dev/R/parallelism_stats.R` — the core sorting statistic (`parallelism_stats()`).
- `dev/R/eMLG_parallelism.R` — `build_sorted_eMLG()`, `build_group_consensus()`, `cluster_DI()`.
- `dev/R/fig_ld_tracks.R` — reads `data/ld_tracks_{a_windows,ldw_persnp}.rds` + recmap → `Figures/ld_tracks_chr26_chr10.{pdf,png}`.
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
- `dev/methods_notes.md` — LD-pruning / eMLG design rationale.
- `dev/HANDOFF_SUMMARY_*.md` — historical thread handoffs.

## Status

- **Module D** (among-region two-locus Ohta D′2st DMI test on unlinked eMLG pairs) — designed, not built. Needs the null (E) to be interpretable.
- **Module E** (recombination-matched haplodiploid neutral null; `dev/R/moduleE_*.R`) — separate workstream. The inference license: until it exists, the descriptive sorting results are "consistent with neutral" only.
