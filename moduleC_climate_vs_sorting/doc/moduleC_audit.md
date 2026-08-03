# Module C (revised) — pre-analysis audit of available objects

**Date:** 2026-07-31
**Scope:** revised, genome-wide, eMLG-level Module C calibrated against the same
10,000 Ω-structured null covariates used in Module B. This document records the
audit required before any new analysis code is run.

Primary configuration (fixed throughout, matching Module B): **Åland excluded,
19 populations, fixed LD-pruned Ω, complete eMLG set (32,840 clusters with
`has_eMLG == TRUE`).**

All paths are relative to the `formica_hybrid` repo root.

---

## 1. Objects located

| Role | File | Status |
|---|---|---|
| Observed eMLG association | `moduleB_climate_GEA/data/moduleB_eMLG_association.rds` | present (32,840 × 12; `eBF1`,`eBF2` = per-eMLG BF on PC1/PC2) |
| eMLG null summary | `moduleB_climate_GEA/data/moduleB_eMLG_null.rds` | present (32,840 × 10) — **summaries only, see §2** |
| Null covariate inputs | `aland_excluded_eMLG/null/null_b01.env … null_b10.env` | present (10 × 1000 Ω-structured draws; ~345 KB each) |
| Raw null BayPass BF outputs | `aland_excluded_eMLG/null/b*_summary_betai_reg.out` | **ABSENT — deleted after parsing (see §2)** |
| eMLG BayPass run inputs | `aland_excluded_eMLG/{u_eMLG.geno, omega_mat_omega.out, u_DIEM.size, eMLG_group_order.txt}` | present |
| Observed eMLG BayPass outputs | `aland_excluded_eMLG/PC{1,2}_eMLG_withOmega_summary_betai_reg.out` | present |
| Canonical clustering | `module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds` (`$groups`) | present (474,014 clusters; 32,840 `has_eMLG`); moved here in the 2026-08-03 reorg (was `data/`) |
| Sorting + DI per eMLG | `data/moduleC_C3_cl.rds` | present (32,840 × 9: `DI`, `directional`, `sort_class`, `uni_score`, `n_loci`, …) |
| Recombination map | `data/Frufa_DTOL_PR.ref_genome.recmap` | present (chr, pos, cM, cM/Mb) |
| Null generator (provenance) | `moduleB_climate_GEA/R/moduleB_eMLG_null.R` | present |

---

## 2. Is the full 32,840 × 10,000 null BF matrix retained? — **No.**

`moduleB_climate_GEA/R/moduleB_eMLG_null.R` ran BayPass in 10 batches of 1,000 null
covariates. For each batch it built the 32,840 × 1,000 BF matrix in memory, used
it **only** to accumulate per-eMLG exceedance counts and a running per-eMLG
maximum, and then **deleted the raw BayPass output** (line 92:
`file.remove(Sys.glob(paste0(pref, "_summary_*.out")))`).

The saved object `moduleB_eMLG_null.rds` therefore contains, per eMLG, only:

- `BF1`, `BF2` — the *observed* PC1/PC2 BF (identical to `eBF1`/`eBF2`);
- `k1`, `k2` — count of the 10,000 nulls whose BF ≥ observed (per-eMLG right tail);
- `p1`, `p2` = (1 + k) / 10,001 — per-eMLG empirical p;
- `null_max` — the single largest null BF seen for that eMLG across all 10,000;
- `floor1`, `floor2` — floor-survivor flags.

There is **no** per-null-covariate genome-wide BF vector anywhere on disk. A
filesystem search for any large null-matrix object returned nothing, and the
`null/` directory holds only the input `.env` covariate files plus logs.

**Consequence for the revised Module C.** The requested primary statistic is a
*within-covariate* rank/percentile of eBF across all 32,840 eMLGs, computed
identically for the observed axes and for **each** of the 10,000 null covariates.
That requires, for every null covariate *j*, the full genome-wide BF vector under
*j*. The saved exceedance counts cannot reconstruct this: `k` collapses each
eMLG's 10,000 null BF values to a single count, discarding which covariate
produced which value and the cross-eMLG structure within a covariate. **No
statistic of the requested form can be recovered from the saved summaries.**

## 3. Can the matrix be recovered without rerunning BayPass? — **No, but regeneration is exact and cheap to make faithful.**

The raw batch outputs are deleted, so recovery from existing BayPass outputs is
impossible. However, the **inputs are fully preserved and deterministic**:

- the 10 null covariate files `null_b01.env … null_b10.env` are on disk;
- they are also reproducible from scratch — `moduleB_eMLG_null.R` draws them with
  `set.seed(SEED_SIM + b)`, `SEED_SIM = 2026`, from the fixed Ω eigendecomposition;
- BayPass was run with a fixed MCMC seed (`-seed 74`) **identical across batches**,
  so re-running on the same `.env` files and the same `u_eMLG.geno`/Ω reproduces
  the same core chain and hence the same null BF (up to negligible MC noise).

Regenerating the matrix therefore means **re-running BayPass on the 10 saved
`.env` files** (≈2 h/batch on 6 threads, ≈20 h total — the original run took
~19.5 h). This is unavoidable given the deletion, but it reproduces exactly the
10,000 Ω-structured null covariates already used for Module B's sim-FDR, keeping
the two analyses consistent.

**Design note:** the giant matrix never needs to be *stored*. Exactly as the
original script accumulated `k` on the fly, the revised Module C can compute each
per-covariate genome-wide statistic (Spearman ρ with DI, with recombination; the
sorted-vs-unsorted percentile gap; the top-fraction summaries) directly from each
batch's BF matrix in memory and retain only the 10,000 null *statistics*
(≈10,000 × a handful of numbers). No multi-GB intermediate is written.

---

## 4. eMLG count and ordering across objects

| Object | N | Key |
|---|---|---|
| `moduleB_eMLG_association.rds` | 32,840 | `group_id` |
| `moduleB_eMLG_null.rds` | 32,840 | `group_id` |
| `eMLG_group_order.txt` (BayPass row order) | 32,840 | `group_id` |
| `data/moduleC_C3_cl.rds` (sorting + DI) | 32,840 | `group_id` |
| `module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds$groups[has_eMLG]` | 32,840 | `group_id` |

Verified programmatically:

- association `group_id` **==** null `group_id` **==** `eMLG_group_order.txt`,
  identical order (`F3, F6, F8, …`). The BayPass BF row order is this order, so
  regenerated null BF rows map to eMLGs positionally **within the BayPass side
  only**.
- association `group_id` ⊆ `moduleC_C3_cl.rds$group_id` (all present), **but the
  row order differs** (`all.equal` = FALSE). → sorting/DI **must** be joined by
  `group_id`, never positionally.
- all association `group_id` ∈ `has_eMLG` clusters.

**Rule adopted:** every join is an explicit `group_id` key join with a
one-to-one assertion (`stopifnot` on set equality and no duplication). Positional
alignment is used only where it has been verified identical (the BayPass row order
vs `eMLG_group_order.txt`), and is re-asserted at run time.

---

## 5. Annotation completeness (join by `group_id`, all 32,840 eMLGs)

- **Diagnostic Index:** `moduleC_C3_cl.rds$DI` — per-eMLG (consensus) DI, 0 NA,
  range −124.9 … −4.2. This is the established eMLG diagnostic-content summary and
  is the primary DI measure. (Member-level fraction with DI > −25 is available as
  a sensitivity measure from the marker map.)
- **Directional sorting:** `moduleC_C3_cl.rds$directional` (0/1), 0 NA; 4,094 of
  32,840 directional. `sort_class` (aquilonia / polyctena / unsorted / NA) and
  `differentiated` also available. This is the Module A/C sorting definition.
- **Sorting magnitude (optional continuous):** `moduleC_C3_cl.rds$uni_score`.
- **Cluster size:** `n_loci` (present in both clustering and C3_cl).
- **Recombination rate:** built per eMLG from `Frufa_DTOL_PR.ref_genome.recmap`
  by interpolating cM/Mb at the cluster's **representative** marker position
  (the established per-cluster convention in `moduleB_architecture.R` B2/B3).
  Verified: 32,840/32,840 finite, range 0.06 … 141.16 cM/Mb.

---

## 6. Audit conclusions

1. All annotations (DI, sorting status, sorting magnitude, cluster size,
   recombination) and both observed eMLG association axes are present, complete,
   and joinable one-to-one by `group_id` for the full 32,840-eMLG universe under
   the primary configuration.
2. The 10,000 Ω-structured null covariates are preserved (`.env` files, and
   independently reproducible from the recorded seed), so the null design is
   intact.
3. **The only missing ingredient is the full null BF matrix**, which was
   intentionally discarded after Module B extracted its exceedance counts. It
   cannot be recovered from any saved output; it must be **regenerated by
   re-running BayPass on the 10 preserved `.env` files** (~20 h), with the
   per-covariate statistics computed on the fly so nothing large is stored.
4. Nothing here overwrites Module B or the earlier Module C outputs; the revised
   analysis is written under `moduleC_climate_vs_sorting/`.

**Gate:** the ~20 h BayPass regeneration is the pre-condition for the primary
analysis and is the decision to confirm before proceeding.

---

## 7. Post-review corrections (2026-08-01)

A code review after the first pilot required the following corrections. They were
applied, the stale checkpoint (built under the previous statistic definitions) was
archived, the annotation table was rebuilt, and the one-batch pilot is re-run from
scratch. No Module B result or historical Module C output is modified.

| # | Correction | Files |
|---|---|---|
| 1 | **Primary sorting is now differentiated-only.** The primary sorting statistic (`sort_gap_differentiated`) contrasts directional vs non-directional eMLGs **within `differentiated == TRUE`**, matching the Module A/C definition. The all-eMLG contrast (`sort_gap_all`) is retained as a supplementary sensitivity. The six-test FDR family = PC1/PC2 × {DI, recombination, differentiated-only sorting}. | stat_functions, analyse |
| 2 | **Continuous sorting magnitude uses `prop_fixed`, not `uni_score`.** `uni_score` = (n_aqu − n_pol)/n_obs is signed *orientation*; `prop_fixed` (degree of fixation) is *magnitude*. Verified: `prop_fixed` ≈0.6 for both aquilonia and polyctena (unsigned), `uni_score` = +0.6 / −0.6 (signed). `prop_fixed` imported and validated; primary magnitude stat = `rho_sort_magnitude` = Spearman(BF, prop_fixed) on the fixed complete-case subset. `uni_score` kept as supplementary "sorting orientation". | annotations, stat_functions, analyse |
| 3 | **Faithful regeneration is mandatory — Monte-Carlo equivalence gate** (see §8: BayPass is not bit-reproducible, so exact `max\|dk\|==0` is unattainable). `moduleC_null_regen.R` `stop()`s before writing the final object unless, **per axis**, Pearson `r(k_regen, k_saved) > 0.99` AND `\|Σk_regen/Σk_saved − 1\| < 0.03`. `moduleC_analyse.R` re-asserts these plus `k_check$reproduced == TRUE`, 10,000 complete finite null rows, and presence of all statistics before doing anything; the report line quotes the observed `r` and relative differences. Completeness, finiteness, ordering and input fingerprints remain **exact** hard stops. Thresholds are pipeline-integrity bounds, not to be relaxed on failure. | null_regen, analyse |
| 4 | **Threshold-label parser fixed.** The old `sub("top0*", …)` collapsed `top0001→0.1` etc. (→ NA). Replaced with an explicit `frac_lookup` keyed by tag, asserted to resolve all three fractions and all metrics. | analyse |
| 5 | **Stronger alignment/numeric validation.** Observed BF asserted equal to `eBF1`/`eBF2` (tol 1e-6, actual Δ = 0); Module B null summary asserted one-to-one and re-ordered by `group_id`; full-length `COVARIABLE`/`MRK` batch-ordering assertion; all parsed BF asserted finite; each completed `null_stats` batch asserted finite; exactly 10,000 complete rows asserted at finalisation; `n_loci` cross-checked between clustering and C3_cl (0 mismatches). | null_regen, annotations |
| 6 | **Checkpoint integrity via fingerprints.** The checkpoint stores md5 fingerprints of the annotation values, `moduleC_stat_functions.R`, the observed PC1/PC2 outputs, the ten `.env` files, geno/omega/poolsize, and the numeric parameters (NSIM/BATCH/NM/seed/opts). On resume these must match exactly, so batches under different definitions cannot be combined. | null_regen |
| 7 | **Pilot acceptance harness.** `moduleC_acceptance.R` checks items 1–5, 7–11 (incl. tamper tests: an altered reproduction flag makes `analyse.R` stop; a mismatched fingerprint makes `null_regen.R` stop before BayPass). `moduleC_determinism_probe.R` runs BayPass twice on the same small covariate file to confirm bit-reproducibility before the long run (checks 5-live and 6 are read off the pilot log). | new scripts |

### Sorting-magnitude / orientation verification (by `sort_class`)

| sort_class | median `uni_score` (orientation) | median `prop_fixed` (magnitude) |
|---|---|---|
| unsorted | 0.0 | 0.20 |
| aquilonia | +0.6 | 0.60 |
| polyctena | −0.6 | 0.60 |

Magnitude is symmetric across the two sorted directions (confirming `prop_fixed` is
unsigned magnitude); orientation is antisymmetric (confirming `uni_score` is signed).

### Statistic set (per covariate; identical for observed and null)

Primary: `rho_DI`, `rho_rec`, `sort_gap_differentiated`. Supplementary: `sort_gap_all`,
`rho_sort_magnitude` (prop_fixed), `rho_sort_orientation` (uni_score), raw-BF variants
(`pear_DI`, `pear_rec`, `pear_sort_magnitude`, `bf_gap_differentiated`, `bf_gap_all`).
Threshold (top 0.1/0.5/1%): median DI, median recomb, proportion differentiated
(`pdiff`), and proportion directional among differentiated (`pdir_diff`, NA + warning
if a fraction has no differentiated eMLG).

---

## 8. BayPass non-determinism and the Monte-Carlo equivalence gate (2026-08-01)

Before the full run, two probes were run (both delete their large raw outputs).

**Determinism probe** (`moduleC_determinism_probe.R`): BayPass run **twice on the same
covariate file** with identical `-nthreads 6 -seed 74` gives BF differing by up to
**24.6 dB**, with **0 of 98,520** values identical. BayPass with multiple threads is
**not bit-reproducible** even with a fixed seed (consistent with Module B's noted
"±1.4 dB MC error"). Each regeneration is therefore a fresh MCMC realization, and the
original `max|dk|==0` gate is **unsatisfiable** — not because of any pipeline error.

**Statistic-stability probe** (`moduleC_stat_stability_probe.R`): the same two
realizations were reduced to the Module C statistics for 20 covariates. Although
individual BF are not reproducible, the **genome-wide statistics are**:

| statistic | within-covariate MCMC SD | between-covariate null SD | ratio |
|---|---|---|---|
| `rho_DI` | 0.00041 | 0.046 | 0.9% |
| `rho_rec` | 0.00061 | 0.0099 | 6.2% |
| `sort_gap_differentiated` | 0.00057 | 0.012 | 4.7% |
| `rho_sort_magnitude` | 0.00056 | 0.030 | 1.9% |
| `rho_sort_orientation` | 0.00069 | 0.058 | 1.2% |

MCMC noise is ≤ ~6% of the null spread in SD terms (< 0.4% in variance); the null
distribution of every statistic is governed by real between-covariate variation, not
MCMC noise. The analysis is therefore valid despite per-BF non-determinism.

**Resolution (agreed with PK).** The exact gate is replaced by a **Monte-Carlo
equivalence gate**: per axis, Pearson `r(k_regen, k_saved) > 0.99` and
`|Σk_regen/Σk_saved − 1| < 0.03` (Pearson, so nonlinear count discrepancies cannot
hide). This is a pipeline-integrity bound (the 1,000-null pilot agreed to ~0.1–0.2%),
**not** a numerical-reproducibility claim; input identity is guaranteed exactly by the
fingerprints, and statistic reproducibility by the stability probe. Observed `r` and
relative differences are reported. On failure the thresholds are **not** relaxed —
the cause is investigated and, if needed, tolerances recalibrated with additional
independent reruns.
