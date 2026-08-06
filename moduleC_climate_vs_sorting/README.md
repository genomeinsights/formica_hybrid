# Module C — Does climate adaptation drive ancestry sorting?

The extrinsic-side test of the headline question. Module B established that **no
individual eMLG survives a genome-wide climate-association FDR**; Module C asks the
complementary, *pattern-level* question: across the whole genome, is the strength
of climate association (BayPass BF on PC1/PC2) organised by **ancestry sorting**,
by the **diagnostic index (DI)**, or by **recombination** — beyond what population
structure alone produces? It reuses Module B's per-eMLG BF and Module B's 10 000
Ω-structured null draws, reduces observed and null to the *same* genome-wide
statistics by identical code, and calibrates each observed statistic against its
structured-null distribution.

**Headline result:** of the six primary tests (PC1/PC2 × DI / recombination /
sorting), only **DI × PC2** survives FDR (a weak, diffuse gradient; ρ ≈ −0.124,
FDR ≈ 0.006). **Directional sorting, sorting magnitude, and recombination are all
indistinguishable from the Ω-structured null on both axes.** So the measured
climate gradients do **not** organise directional ancestry sorting beyond
structure. The calibration is reported over a **(min\_n\_loci × τ) grid**
(min ∈ {5, 10}, τ ∈ {0.5, 0.6, 0.8}); the conclusion is unchanged across it, and
DI × PC2 is if anything **stronger in the min = 10 high-LD / selection-candidate
subset** (ρ ≈ −0.138, p ≈ 0.0003). See [`doc/moduleC_report.md`](doc/moduleC_report.md)
for the full grid table, interpretation, and caveats.

Primary configuration: **Åland-excluded** (19 populations), **fixed LD-pruned Ω**,
primary eMLG universe = **32 840 clusters** (min = 5; min = 10 → 14 349 clusters),
**NSIM = 10 000** structured-null covariates. Because Ω is fixed (passed via
`-omegafile`), each cluster's BF is invariant to the size threshold, so **min is a
cheap analysis knob**: min = 10 is a row-subset of the *same* null BF matrices,
re-reduced without re-running BayPass.

---

## How it relates to the other modules

```
module0  ──►  eMLG_5loci_0025_cM05.rds (canonical clustering, cluster size)
moduleA  ──►  moduleA_cluster_sorting.rds  (per-eMLG DI / sort_class / prop_fixed / uni_score)
moduleB  ──►  moduleB_eMLG_association.rds  (observed eMLG BF)  +  aland_excluded_eMLG/null/  (10k null .env)
                                   │
                                   ▼
                        moduleC_climate_vs_sorting  (this module)
```

Module C is where the sorting arm (A) and the climate arm (B) are joined: every
statistic contrasts Module B's climate BF against Module A's sorting/DI
annotation, on Module 0's clustering.

---

## Dependencies

- **R** packages: `data.table` (+ `ggplot2` for the calibration figure).
- **BayPass** v3 binary — required by `moduleC_null_regen.R`, which re-runs BayPass
  on the preserved null covariates to rebuild the null BF matrix (now **persisted**
  to `aland_excluded_eMLG/null/bf_matrices/`, so future (min, τ) changes re-reduce
  without BayPass). Path is set at the top of that script — **edit it for your machine.**
- The three shared repo-root inputs and the Module B eMLG BayPass run directory
  (`aland_excluded_eMLG/`, incl. `null/` and `eMLG_group_order.txt`).

---

## Inputs (external / upstream)

| Input | Used for | Produced by |
|---|---|---|
| `moduleB_climate_GEA/data/moduleB_eMLG_association.rds` | observed per-eMLG BF (`eBF1/eBF2`) | Module B |
| `aland_excluded_eMLG/null/null_b01..b10.env` (+ geno, Ω, poolsize, params) | the 10 000 Ω-structured null covariates | Module B (kept only exceedance counts; the full BF matrix is regenerated here) |
| `aland_excluded_eMLG/eMLG_group_order.txt` | canonical BayPass eMLG row order | Module B |
| `moduleA_sorting/data/moduleA_cluster_sorting.rds` | per-eMLG DI, `sort_class`, `directional`, `prop_fixed`, `uni_score` | Module A |
| `module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds` | cluster size (`n_loci`), cross-check | Module 0 |
| `data/Frufa_DTOL_PR.ref_genome.recmap` | recombination rate (cM/Mb) at the representative marker | external |

---

## Pipeline & run order

Run scripts with `Rscript moduleC_climate_vs_sorting/R/<script>` from the repo root.

| # | Script | Does | Key output |
|---|---|---|---|
| — | `moduleC_stat_functions.R` | Shared genome-wide statistic reducers (DI/recomb/sorting contrasts). **Sourced**, not run standalone — guarantees observed and null are reduced by identical code. | — |
| 1 | `moduleC_annotations.R` | Build the aligned per-eMLG annotation table (BF + sorting/DI + recomb), ordered to the BayPass rows, with fail-loud one-to-one join assertions. | `data/moduleC_annotations.rds` |
| 2 | `moduleC_null_regen.R` | **[BayPass, 10 threads]** Regenerate the 32 840 × 10 000 null BF matrix from the preserved `.env` draws and reduce each covariate against every (min × τ) grid cell; **persists** the per-batch BF matrices (`aland_excluded_eMLG/null/bf_matrices/`) so future (min, τ) re-reduce without BayPass. Checkpointed; Monte-Carlo equivalence gate vs Module B counts. | `data/moduleC_null_stats.rds` (`by_cell` grid) |
| 3 | `moduleC_analyse.R` | Structured-null calibration: two-sided empirical P per statistic, BH-FDR over the six primary tests (at the primary cell), the min×τ grid-sensitivity table, calibration figure, and the report. | `data/moduleC_results.rds`, `data/moduleC_primary_tests.tsv`, `data/moduleC_grid_sensitivity.tsv`, `doc/moduleC_report.md` |
| — | `moduleC_acceptance.R` | Pre-pilot acceptance checklist (scripts parse; annotation = 32 840 unique eMLGs in BayPass order; joins one-to-one). | console |
| — | `moduleC_determinism_probe.R`, `moduleC_stat_stability_probe.R` | Diagnostics behind the faithful-regeneration gate: BayPass BF are not bit-reproducible (MCMC noise), but the genome-wide statistics are stable to ≪ the null spread. | `data/stat_stability_probe.rds` |

---

## Outputs

**`data/`** — `moduleC_annotations.rds` (aligned per-eMLG table, carries
`directional_tau05/06/08` + `n_loci`), `moduleC_null_stats.rds` (observed + 10k-null
genome-wide statistics per (min × τ) cell in `by_cell`), `moduleC_results.rds`
(empirical P, FDR), `moduleC_primary_tests.tsv` (the six-test table),
`moduleC_grid_sensitivity.tsv` (min × τ grid), `stat_stability_probe.rds`.
The persisted null BF matrices live outside `data/` in `aland_excluded_eMLG/null/bf_matrices/`
(~2.4 GB, git-ignored) — kept so any future (min, τ) is re-reducible without BayPass.

**`doc/`** — `moduleC_report.md` (results write-up + interpretation), `moduleC_audit.md`
(provenance / reproducibility audit), `supplementary_methods_moduleC_climate_genomewide.{tex,pdf}`
(Materials & Methods).

---

## What this module can and cannot establish

- **Can:** calibrate genome-wide, eMLG-level climate-association *patterns* against a
  structure- and architecture-preserving null, with the same eMLG universe and
  statistic for observed and null.
- **Cannot:** identify individual climate-associated loci (that is Module B's job,
  and none survive its sim-FDR); the weak PC2×DI gradient is a diffuse property, not
  a set of adaptation loci, and must not be read as validating any candidate eMLG.

## Reference

LD-cluster / eMLG association approach: Li Z, Kemppainen P, Rastas P, Merilä J.
2018. Linkage disequilibrium clustering-based approach for association mapping with
tightly linked genomewide data. *Mol Ecol Resour* 18:809–824.
