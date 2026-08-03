# formica_hybrid — climate adaptation and the predictability of ancestry sorting

Population genomics of **20 hybrid populations** between *Formica polyctena* and *F. aquilonia*, with allopatric parental references. Haplodiploid ants (haploid males). The pipeline asks not merely whether parental ancestry sorts in hybrids, but whether it sorts **predictably** — the same parental allele approaching fixation across independent replicate populations — and, if so, whether that predictable directional sorting is driven by **extrinsic climate adaptation**.

This repository is organised as **four sequential modules**, each in its own folder with its own `README.md`, `R/`, `data/`, `Figures/` and `doc/`. Each module builds on the ones before it and reads/writes version-tagged shared objects.

> The intrinsic counterpart (Module D, two-locus incompatibilities), the neutral-null inference licence (Module E, haplodiploid SLiM), and the genome-wide epistasis scan are **downstream** of these four modules and are developed separately (`dev/R/moduleD_*`, `moduleE_slim/`). They are not part of this four-module core.

## Pipeline at a glance

```         
                module0_ld_pruning
   parse DIEM → genotype matrices + LD decay → canonical cM05 clustering
                          │
          ┌───────────────┴───────────────┐
   moduleA_sorting                   moduleB_climate_GEA
   sorting phenomenon +              per-eMLG BayPass climate
   genomic architecture              association + 10k-null sim-FDR
          │                                │
          └───────────────┬───────────────┘
                          ▼
                moduleC_climate_vs_sorting
   does the genome-wide climate-association pattern track ancestry
   sorting / DI / recombination beyond population structure?
```

Modules A and B are **parallel arms** off the shared clustering; module C is where they converge. Everything joins on **one canonical clustering** and the **two genotype matrices**.

## The four modules

| \# | Folder | Question it answers | Reads | Writes (key) |
|---------------|----------------------|---------------|---------------|---------------|
| 1 | [`module0_ld_pruning`](module0_ld_pruning/README.md) | LD decay + LD-based complexity reduction: reduce \~1.1 M markers to LD clusters, each with an eMLG consensus | raw DIEM, recmap, sample table | genotype matrices (shared), LD-decay fits, **`eMLG_5loci_0025_cM05.rds`** (canonical clustering) |
| 2 | [`moduleA_sorting`](moduleA_sorting/README.md) | The genetic architecture of ancestry sorting: extent, parental direction, and its relation to DI / recombination / diversity / cluster size | clustering, genotype matrices, recmap | `moduleA_snp.rds`, `moduleA_clusters.rds`, `moduleA_architecture.rds`, **`moduleA_cluster_sorting.rds`** (per-eMLG sorting annotation for the climate modules) |
| 3 | [`moduleB_climate_GEA`](moduleB_climate_GEA/README.md) | Genotype–environment (climate) association: is any locus associated with climate PC1/PC2 beyond structure? (BayPass BF + a 10 000-draw Ω-structured sim-FDR) | clustering, genotype matrices, BayPass runs | `moduleB_eMLG_association.rds`, `moduleB_eMLG_null.rds`, candidate/outlier tables |
| 4 | [`moduleC_climate_vs_sorting`](moduleC_climate_vs_sorting/README.md) | To what extent does climate adaptation drive ancestry sorting? Calibrate the genome-wide climate-association *pattern* (vs sorting / DI / recombination) against the same structured null | module B's BF + null draws, module A's sorting annotation, clustering | `moduleC_annotations.rds`, `moduleC_null_stats.rds`, `moduleC_results.rds` |

**Headline findings.** Sorting is predictable and directional (A). No individual locus survives a genome-wide climate-association FDR (B). At the pattern level, the climate gradients do **not** organise directional sorting or recombination beyond structure — only a weak, diffuse `DI × PC2` gradient survives (C). Read as an extrinsic-driver test, the answer so far is largely negative; the inferential licence to call any of this non-neutral is Module E's job (downstream).

## How to run

Every script is run **from the repo root** (`~/gitlab/formica_hybrid`) with working-directory-relative paths, e.g.

``` bash
Rscript moduleA_sorting/R/moduleA_sorting_phenomenon.R
```

R analysis uses the `LDscnR` package (`devtools::load_all("~/gitlab/LDscnR/")`) for LD decay / complexity reduction / consensus construction. Module B and Module C's null regeneration also need the **BayPass** v3 binary (path set at the top of those scripts — edit for your machine). Run order within each module is in that module's README; across modules it is 0 → A → B → C (A and B may run in either order).

## Shared objects and layout

`data/` (git-ignored, regenerable) holds the shared inputs every module reads:

| file | contents |
|------------------------------------|------------------------------------|
| `hybrids_only_maf005.Rdata` | `GTs_hybrids_005`, `map_hyb_005` (incl. `ld_w_095`, `DiagnosticIndex`), `sample_data` — hybrids only, for all LD estimation & clustering |
| `hybrids_and_parents_maf005.Rdata` | genotypes incl. parents — only where parental allele frequencies are needed (ancestry orientation, divergence) |
| `Frufa_DTOL_PR.ref_genome.recmap` | recombination map (cM, cM/Mb) |
| `Formica_hybrids_filtered_diem_output.bed.gz`, `Sample_covariate_info_outlier_analysis_20.txt` | raw DIEM output + sample/covariate table (Module 0 inputs) |

**Two genotype matrices, by design:** hybrids-only for all LD estimation and clustering (including parents would let parent–hybrid structure dominate LD); hybrids+parents only where parental allele frequencies are needed.

Each module writes its own products under `<module>/data/` and `<module>/Figures/`, and reads upstream products from the sibling module's `data/` folder (paths are always module-qualified, e.g. `module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds`). `data_legacy/` holds pre-reorganisation artifacts (see its README); nothing in the four-module pipeline reads from it.

## Conventions

- **The LD cluster, not the marker, is the unit of analysis.** Markers in LD carry no independent information; tests treat clusters as observations. A signal present at marker level but absent at cluster level is spatial pseudo-replication. The units are *LD clusters*, never "haplotype blocks" (a biallelic SNP carries one bipartition, so markers on different genealogical branches form different clusters that can overlap the same interval).
- **Canonical clustering = `eMLG_5loci_0025_cM05.rds` (0.5 cM).** Version-tagged by its parameters; never overwrite a fixed name (group IDs are not stable across clusterings). The 1 cM build is kept only as a sensitivity variant.
- **A predictor is never also used to select the data.** In particular the diagnostic index (DI) is a covariate throughout and is kept **ungated**.
- **eMLG consensus for magnitude/LD tests; pruning representatives suffice for direction** (see `moduleB_eMLG_vs_rep`).
- **Parameter-tagged outputs.** Threshold-dependent outputs carry the settings in the filename (e.g. `_5_15`), so alternative settings coexist.
- **Descriptive by design.** Each pattern is, on its own, compatible with neutral drift as well as selection; establishing departure from neutrality (Module E) and discriminating its cause (D vs C) are downstream.

## Documentation status

The per-module `README.md` files and this file are current. Some higher-level documents (`doc/`, `manuscript_notes/`, the flat `R/` script set) still describe the earlier three-module scheme (A = sorting, B = architecture, C = climate) and are being revised; prefer the module READMEs and this file where they disagree.
