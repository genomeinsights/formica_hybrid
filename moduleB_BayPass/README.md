# Module B — BayPass climate association

Tests whether population allele frequencies are associated with the climate
covariates **PC1** and **PC2**, at two resolutions — individual SNPs (standard
BayPass outlier scan) and **LD-complexity-reduced clusters** (eMLGs; one
consensus genotype per LD cluster, *sensu* Li et al. 2018) — then calibrates the
cluster-level signal against a population-structure null to obtain a valid
genome-wide FDR.

**Headline result:** after Ω-structured population-structure control and a valid
sim-FDR, the climate signal does **not** survive genome-wide. The SNP
member-count candidates are differentiation/structure artefacts; only 5 PC2
"eMLG-only" clusters resist structure (set FDR ≈ 0.66, suggestive only). See
[`doc/moduleB_summary.pdf`](doc/moduleB_summary.pdf) for the full write-up,
figures and tables.

Primary configuration throughout: **Åland-excluded** (19 populations, 154
individuals), **fixed LD-pruned Ω** (`withOmega`); marker significance =
**BF(dB) ≥ 15**. Full (20-pop) and internally-estimated-Ω runs are sensitivity
comparisons.

---

## Directory layout

```
moduleB_BayPass/
├── README.md                 (this file)
├── R/                        pipeline scripts (run from the repo root)
├── data/                     output objects (.rds/.csv) + null_env/ (archived null draws)
├── Figures/                  Manhattans + PC–ancestry confound (.png)
└── doc/                      moduleB_summary.{md,pdf}, render script, header;
                              doc/tables/  standalone candidate tables (.tex/.pdf)
```

**Run scripts from the repo root**, e.g. `Rscript moduleB_BayPass/R/moduleB_climate_candidates.R`.
Scripts read shared/raw inputs at their repo-root locations and write all
products under `moduleB_BayPass/`.

---

## Dependencies

- **R** packages: `data.table`, `ggplot2`, `patchwork`.
- **BayPass** v3 binary. Path is set at the top of the run scripts
  (`/Users/petrikem/gitlab/baypass_public-master/sources/g_baypass`) — **edit
  this for your machine.**
- **pandoc + xelatex** (only to rebuild `doc/moduleB_summary.pdf`).

---

## Inputs (external / shared — not part of this module)

| Input | Used for | Produced by |
|---|---|---|
| `module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds` | LD clustering + per-individual eMLG dosages | Module 0 (LD pruning) |
| `data/hybrids_only_maf005.Rdata` | `map_hyb_005`, `GTs_hybrids_005`, `sample_data` (PC1/PC2) | upstream marker processing |
| `data/hybrids_and_parents_maf005.Rdata` | parents + hybrids (ancestry estimate) | upstream marker processing |
| `data/pruned_markers.rds` | LD-pruned marker set for the fixed Ω | Module 0 |
| `with_aland/`, `aland_excluded/` | SNP-level BayPass run dirs (inputs, `run_baypass.sh`, `*_summary_*.out`) | this module (step 1–2) |
| `aland_excluded_eMLG/` | eMLG-level BayPass run dir (inputs, `*_summary_*.out`, `null/` scratch) | this module (step 5–7, 9) |

The three BayPass run directories are large (genotype count files + `*_summary_*.out`)
and live at the **repo root**, not inside the module.

---

## Pipeline & run order

Steps marked **[BayPass]** are the external BayPass runs (shell scripts that call
the `g_baypass` binary); all other steps are `Rscript … moduleB_BayPass/R/<script>`.

| # | Script | Does | Key output |
|---|---|---|---|
| 1 | `moduleB_prepare_with_aland.R`, `moduleB_prepare_aland_excluded.R` (source `moduleB_write_baypass_inputs.R`) | Build BayPass genotype-count / poolsize / PC-covariate inputs | files in `with_aland/`, `aland_excluded/` |
| 2 | **[BayPass]** `with_aland/run_baypass.sh`, `aland_excluded/run_baypass.sh` | Ω + PC1/PC2 association, with & without fixed Ω | `…/PC{1,2}_DIEM_{noOmega,withOmega}_summary_*.out` |
| 3 | `moduleB_analyse_baypass.R` | SNP-level Manhattans (8 configs; member-count candidate colouring) | `Figures/manhattan_*_min5SigLoci_*.png` |
| 4 | `moduleB_ancestry_confound.R` | PC provenance + PC–ancestry confound check | `data/moduleC_ancestry_confound.rds`, `Figures/moduleC_ancestry_confound.png` |
| 5 | `moduleB_climate_candidates.R` | Member-count candidates → eMLG-BF gate → candidate table | `data/moduleB_climate_candidates.rds`, `doc/tables/climate_candidate_*.{tex,pdf}` |
| 6 | `moduleB_write_eMLG_baypass_inputs.R` | Build eMLG count file (round eMLG dosage → per-pop counts; reuse fixed Ω) | `aland_excluded_eMLG/u_eMLG.geno` (+ staged Ω/covariates) |
| 7 | **[BayPass]** `moduleB_BayPass/R/run_baypass_eMLG.sh` | PC1/PC2 association on the 32,840 eMLGs, fixed Ω | `aland_excluded_eMLG/PC{1,2}_eMLG_withOmega_summary_betai_reg.out` |
| 8 | `moduleB_eMLG_manhattan.R` | eMLG-level Manhattan + eMLG-only outlier set | `data/moduleB_eMLG_association.rds`, `data/moduleB_eMLG_outliers.csv`, `Figures/moduleB_eMLG_manhattan.png` |
| 9 | `moduleB_eMLG_null.R` | 10 000 Ω-structured null covariates (10×1000 batches) → per-eMLG sim-FDR. **~21 h**; checkpointed; parses & deletes the ~43 GB raw output | `data/moduleB_eMLG_null.rds`; null draws archived in `data/null_env/` |
| 10 | re-run `moduleB_eMLG_manhattan.R` | picks up step 9 → adds FDR floor-survivor overlays | `Figures/moduleB_eMLG_manhattan.png` (triangles) + `Figures/moduleB_fdr_snp_manhattan.png` |

Notes:
- Steps 3, 4 and 5 depend on step 2; step 5 also needs step 7 (eMLG `.out`).
- Step 8 reads `data/moduleB_eMLG_null.rds` **if present** (guarded), so step 10 is
  simply re-running step 8 after step 9 to draw the FDR overlays.

---

## Outputs

**`data/`**
| File | Contents |
|---|---|
| `moduleB_climate_candidates.rds` | 35 member-count candidates + eMLG BF, `pass_eMLG` flag |
| `moduleB_eMLG_association.rds` | all 32,840 eMLGs: position, member sig counts, eMLG BF per axis, support category |
| `moduleB_eMLG_outliers.csv` | 51 eMLG outliers (BF ≥ 15) with per-SNP `n_sig`/`pct_sig` retained |
| `moduleB_eMLG_null.rds` | 10k sim-FDR: `BF1/BF2`, `k1/k2` (#nulls ≥ obs), `p1/p2`, `null_max`, `floor1/floor2`; `attr(.,"meta")` |
| `moduleC_ancestry_confound.rds` | per-population PC1/PC2 vs aquilonia ancestry |
| `null_env/null_b01…b10.env` | the 10,000 Ω-structured null covariates (byte-exact reproducibility of step 9) |

**`Figures/`** — 2 primary SNP Manhattans (`manhattan_PC{1,2}_aland_excluded_withOmega_…png`),
`moduleB_eMLG_manhattan.png` (eMLG Manhattan; FDR survivors = triangles),
`moduleB_fdr_snp_manhattan.png` (all SNPs; FDR clusters coloured), `moduleC_ancestry_confound.png`.
(Step 3 also writes the other 6-config sensitivity Manhattans to `Figures/`.)

**`doc/`** — `moduleB_summary.{md,pdf}` (results write-up: main-text conclusions,
supplementary methods, complete Tables S1–S4, figures); `doc/tables/climate_candidate_*.{tex,pdf}`.
Rebuild the PDF with:

```bash
bash moduleB_BayPass/doc/render_moduleB_summary.sh
```

---

## Reference

LD-cluster / eMLG association approach: Li Z, Kemppainen P, Rastas P, Merilä J.
2018. Linkage disequilibrium clustering-based approach for association mapping
with tightly linked genomewide data. *Mol Ecol Resour* 18:809–824.
