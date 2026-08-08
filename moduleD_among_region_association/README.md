# Module D — Among-region association: a two-locus screen for intrinsic incompatibilities

The **intrinsic-side** test of the headline question. Modules A–C characterise the
*per-locus* direction of ancestry sorting and its (non-)link to climate. Module D asks the
complementary, **pair-level** question: between **unlinked** LD-reduced units, is there
association *beyond* what shared ancestry and relatedness predict — the signature of a
two-locus (Dobzhansky–Muller) incompatibility?

**Headline result (descriptive, pre-simulation):** largely negative for widespread
pair-specific epistasis. After a structure-corrected scan and one reproducible read-out, the
strongest signal is **technical** (cross-chromosome assembly paralogy) or a
**low-recombination co-ancestry** component (which resolves into several distinct
multi-chromosome modules), not pair-specific. The candidate residue is small: an **aquilonia
co-sorting module** (F33454, belongs with Modules A/C) and a single **ancestry-independent
trans pair** (F11431–F49480). Low recombination is treated as an **annotation, never a
filter**, so every unit is carried into the null; the intrinsic call (excess over neutral,
judged at each pair's *local* recombination rate) is deferred to **Module E**. Full write-up
in [`doc/supplementary_results_module_D.pdf`](doc/supplementary_results_module_D.pdf); minimal
pipeline spec in [`doc/moduleD_plan.md`](doc/moduleD_plan.md).

*Unlinked* = different chromosome, or same chromosome > 10 cM on the genetic map (≈99%
admixture-LD decay in ~50 generations). Scripts are run **from the repository root**; shared
inputs are read from the repo `data/`, and this module's outputs are written under
`moduleD_among_region_association/{data,Figures}/`.

---

## How it relates to the other modules

```
module0  ──►  data/eMLG_5loci_0025_cM05.rds        (canonical clustering, eMLG consensus)
moduleA  ──►  data/moduleA_snp.rds                 (SNP-level sort_class -> bidirectional annotation)
moduleC  ──►  data/moduleC_C3_cl.rds               (differentiated / DI / sort_class gate)
             + data/Frufa_DTOL_PR...recmap         (genetic map: the unlinked / recomb definitions)
                                   │
                                   ▼
                         Module D (this folder)
                                   │
                                   ▼
moduleE  ◄──  Ohta.R + moduleD_ohta_dmi.R           (reused as E's LD-measurement instrument;
                                                     the recombination-matched neutral null)
```

## Pipeline (`R/`)

Run from the repo root, e.g. `Rscript moduleD_among_region_association/R/moduleD_emmax.R`.

**`R/moduleD_bidirectional.R`** — reads `data/moduleA_snp.rds`, clustering → writes
`data/moduleD_bidirectional.rds` (`reps`) — the set of bidirectionally sorting units, used
only as a descriptive **annotation** (the "bidirectional" category in the sorting rings), not
a DMI screen.

**`R/moduleD_emmax.R`** — reads clustering, `hybrids_only_maf005.Rdata`, `moduleC_C3_cl.rds`,
recmap (and `data/moduleD_ohta.rds` for focal picks); sources `R/{emmax,moduleD_paralogy}.R`
→ writes `data/moduleD_emmax.rds`, `Figures/moduleD_emmax_manhattans.pdf` — the
structure-corrected two-locus scan: each differentiated eMLG dosage as a trait vs all others
in an EMMAX LMM with a **double-LOCO** VanRaden GRM from differentiated units, each focal
**conditioned on the top 10 genome PCs**; signed cis/trans by the GLS coefficient (λ ≈
0.98–1.05).

**`R/moduleD_paralogy.R`** (shared filter) — `flag_paralogy()`: median within-population
|r| > 0.9 flags cross-chromosome assembly duplicates a distance rule cannot catch; stores
concordance + excess Ho as corroboration.

**`R/moduleD_network_build.R`** — reads `data/moduleD_emmax.rds`, clustering,
`moduleC_C3_cl.rds`, `hybrids_only`, `data/moduleD_bidirectional.rds`, recmap; sources
`R/moduleD_paralogy.R` → writes `data/moduleD_network.rds` — the single reproducible
read-out: global BH-FDR (`q<0.01`) → paralogy filter → within-chromosome **average-linkage**
merge (`|r|>0.5`, ≤10 cM) → single-population **leverage** filter (leave-one-population-out
`|r|≥0.3`) → low-recombination **structure annotation** (percentile < 0.1; carried to E, never
a filter) → cross-chromosome correlation **modules** (`|r|>0.4`). Emits meta-nodes + meta-edges.

**`R/moduleD_network_circos.R`** — reads `data/moduleD_network.rds`, clustering, `hybrids_only`,
`moduleC_C3_cl.rds`, `data/moduleD_bidirectional.rds`, recmap → writes
`Figures/moduleD_{trans,candidate,structure,module}_circos.png` + `moduleD_trans_network.png`.

**`R/moduleD_module_heatmaps.R`** — reads `data/moduleD_network.rds`, clustering,
`moduleC_C3_cl.rds`, `hybrids_only`, recmap → writes `Figures/moduleD_module_<id>_heatmap.png`
— per-module genotype heatmaps (each module = one coherent multi-chromosome co-ancestry block).

**`R/moduleD_wip_figures.R`** — reads `data/moduleD_{emmax,ohta}.rds` → writes
`Figures/moduleD_{emmax,paralogy}_writeup.png` (the method-diagnostic figures used in `doc/`).

**`R/moduleD_ohta_dmi.R`** — reads clustering, `hybrids_only`, `moduleC_C3_cl.rds`; sources
`R/Ohta.R` → writes `data/moduleD_ohta.rds`, `Figures/moduleD_fig4.{pdf,png}` — among-population
Ohta LD on unlinked pairs. Retained primarily as **Module E's measurement instrument**
(reusable `moduleD_pop_freq_matrix()` / `moduleD_prefilter()` / `moduleD_scan()`), not a
standalone screen; the descriptive result (`D²st ≫ D²is`) is the structure signature.

**Helpers:** `R/emmax.R` (EMMAX LMM; Li, Kemppainen, Rastas & Merilä), `R/Ohta.R` (Ohta
decomposition). Copied in from `dev/R/` (also used by Module E).

## `doc/`

- [`supplementary_methods_ld_module_D.pdf`](doc/supplementary_methods_ld_module_D.pdf) — methods (manuscript style).
- [`supplementary_results_module_D.pdf`](doc/supplementary_results_module_D.pdf) — results (manuscript style).
- [`moduleD_plan.md`](doc/moduleD_plan.md) — the minimal-pipeline spec (what to keep and why).
- [`moduleD_wip_summary.pdf`](doc/moduleD_wip_summary.pdf) — extended working notes.
