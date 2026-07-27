# Pipeline documentation — Stage 0 through Modules A/B/C

This folder documents the *formica_hybrid* analysis pipeline from raw DIEM
genotypes through the three descriptive modules (A, B, C), plus a shared
block-bootstrap sensitivity pass on their main coefficients. It is written to be
**self-contained**: methods, parameters, reasoning, results and figures are all
here, so a reader can understand the workflow and its findings without opening
the code. It doubles as the source material for the manuscript's Materials &
Methods and Results.

Each document maps to one stage of the pipeline and to the corresponding scripts
in `R/`:

| Document | Covers | Scripts (`R/`) |
|---|---|---|
| [01 — LD decay & complexity reduction](01_stage0_ld_and_complexity_reduction.md) | Stage 0: data, LD-decay model, per-marker local-LD support (`ld_w`), two-stage LD clustering → eMLG consensus genotypes | `LD_decay_from_DIEM.R`, `ld_pruning_DIEM.R` |
| [02 — Module A: the sorting phenomenon](02_moduleA_sorting_phenomenon.md) | Whether ancestry sorts *predictably* across replicate populations, and what governs its direction | `moduleA_sorting_phenomenon.R`, `moduleA_di_asymmetry.R` |
| [03 — Module B: genomic architecture](03_moduleB_genomic_architecture.md) | How sorting relates to recombination, diversity, divergence and cluster size | `moduleB_architecture.R`, `moduleB_eMLG_vs_rep.R` |
| [04 — Module C: association with climate](04_moduleC_climate_association.md) | Whether ancestry-diagnostic and directionally sorted regions coincide with genotype–environment (climate) association | `moduleC_analyse_baypass.R`, `moduleC_diagnostic_index_enrichment.R`, `moduleC_ancestry_confound.R`, `moduleC_sorting_climate.R`, `moduleC_rate_based.R` |
| [05 — Sensitivity: block bootstrap](05_sensitivity_block_bootstrap.md) | Dependence-aware confidence intervals for the main A/B/C coefficients (chromosome + 10 cM blocks) | `sensitivity_block_bootstrap.R` |

## The scientific question

Twenty independently formed hybrid populations between *Formica polyctena* and
*F. aquilonia*, with allopatric parental references, let us ask not merely
*whether* parental ancestry sorts in hybrids but whether it sorts **predictably**
— the same parental allele approaching fixation repeatedly across replicates.
The headline question the full pipeline is built to discriminate is whether such
predictable, directional sorting is driven by **intrinsic** genetic
incompatibilities (Dobzhansky–Muller interactions) or by **extrinsic** climate
adaptation.

Modules A–C, documented here, establish and characterise the phenomenon
descriptively. Two further modules complete the discrimination and are **not**
covered here because they are separate workstreams:

- **Module D (intrinsic arm):** among-population two-locus linkage disequilibrium
  between *unlinked* diagnostic clusters (an Ohta decomposition). Experimental /
  under revision; scripts remain in `dev/R/`.
- **Module E (neutral null):** a recombination-matched, haplodiploid demographic
  simulation seeded from real haploid-male haplotypes, which supplies the
  inference licence. Built in a separate session.

Everything in Modules A–C is therefore **descriptive**. Each pattern is, on its
own, compatible with neutral drift as well as with selection; establishing
departure from neutrality (Module E) and discriminating its cause (D vs C) are
deferred to those modules.

## Three organising principles

The same three ideas recur throughout and explain most of the design choices:

1. **The LD cluster, not the marker, is the unit of analysis.** Markers in
   linkage disequilibrium do not carry independent information. Quantities are
   measured per cluster and statistical tests treat clusters — not markers — as
   observations. Reporting a result at both the marker and the cluster level is
   used diagnostically: a signal present at marker level but absent at cluster
   level is spatial pseudo-replication, not a property of the units.
2. **A variable to be tested as a predictor is never also used to select the
   data.** In particular the diagnostic index (DI) is a predictor throughout, so
   the data are gated on *parental allele frequency*, never on DI — otherwise the
   question "does DI predict sorting?" would be circular.
3. **A criterion whose stringency would depend on cluster size is expressed as a
   rate, not a count.** A fixed count of significant members within a cluster
   makes eligibility scale with size; a proportion (significance rate) does not.

## Key shared objects

All modules join on the same two genotype matrices and the same clustering.
`data/` and `Figures/` are git-ignored (large / regenerable).

| File | Contents | Produced by |
|---|---|---|
| `data/hybrids_only_maf005.Rdata` | `GTs_hybrids_005` (164 hybrids × 1,114,340 markers), `map_hyb_005` (incl. `ld_w_095`, `DiagnosticIndex`), `ld_decay`, `sample_data` | Stage 0a |
| `data/hybrids_and_parents_maf005.Rdata` | `GTs_with_parents`, `sample_data_with_parents`, `map_hyb_005` | Stage 0a |
| `data/eMLG_5loci_0025_cM05.rds` | the canonical clustering: `$groups`, `$eMLG`, `$pruned`, `$params` (474,014 clusters; 32,840 with a stored consensus) | Stage 0b |
| `data/Frufa_DTOL_PR.ref_genome.recmap` | recombination map (cM, cM/Mb) | external |

**Why two genotype matrices.** Hybrids-only is used for all LD estimation and
clustering: including parents would let parent–hybrid differentiation dominate
every pairwise LD estimate and obscure the within-hybrid-zone structure the
clustering is meant to capture. Hybrids+parents is used only where parental
allele frequencies are needed (ancestry orientation and between-species
divergence).

## Conventions used in these documents

- **Figures** are referenced by their path under `../Figures/`. Because
  `Figures/` is regenerable and git-ignored, the images may not render on a
  remote host; every figure therefore carries a complete textual caption so the
  document is informative on its own.
- **Parameters** are quoted with their exact values; a consolidated table is at
  the end of document 01.
- **Numbers** are from the runs recorded in `data/module*.rds`
  (Modules A/B ran 2026-07-21; Module C 2026-07-22/23, `_5_15` threshold set).
- The condensed, publication-formatted version of the methods lives separately in
  `manuscript_notes/supplementary_methods_pipeline.{tex,pdf}`; this folder is the
  fuller, descriptive companion and also carries the results.
