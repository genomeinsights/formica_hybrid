# Handoff: authoritative rho05 outputs for population-partitioning follow-up

Audit-only document. No new analyses were implemented here, and no completed
result was modified except where an error is explicitly noted below (the
parental-MAF fold fix, corrected in place with full justification). Nothing
was copied into `module_population_partitioning/` -- read the paths below
directly from `module_manuscript_rho05/`.

Git HEAD at the time of this audit: `d0a3da7e1fbb5716dd0c8d1da60d142a687d810f`
(2026-09-13 14:53:45 +0300).

---

## 1. Authoritative script and saved object: empirical DI vs per-locus Fst

- Script: [R/di25_fst_vs_di_rho05.R](R/di25_fst_vs_di_rho05.R)
- Saved object: `data/di25_fst_vs_di_rho05.rds`
- Figures: `Figures/di25_fst_vs_di_rho05.{png,pdf}`
- `run_fingerprint$script_version = "di25_fst_vs_di_rho05_AUDITFIX_2026-09-12"`

This script computes, per LD-reduced unit (see §5), the empirical
Weir & Cockerham per-locus Fst against the diagnostic index (DI), gated on
folded parental MAF (see §2), and compares each DI-bin's empirical Fst to a
POOLED-per-replicate/per-bin neutral-simulation summary (see §3a/3b).

A companion, richer per-locus violin view lives in a separate script/object
(§2 below) and is the one to use for any analysis needing per-locus Fst
values joined to unit-level covariates (sort_class, prop_fixed, DI, MAF) --
`di25_fst_vs_di_rho05.rds` itself only stores DI-bin summaries, not
per-unit rows.

## 2. Corrected parental MAF definition and retained unit count

Parental allele frequency is `pf <- colMeans(parent_genotypes)/2` (range
`[0,1]`). It must be folded to minor-allele frequency before gating or
binning:

```r
pf_folded <- pmin(pf, 1 - pf)
```

**This fold was missing until 2026-09-13** in both `R/di25_fst_vs_di_rho05.R`
and `R/di25_fst_vs_di_violin_rho05.R`. Effect: units with true MAF as low as
0 could pass a nominal "MAF>=0.15" gate whenever their raw `pf` was >=0.85
(e.g. `pf=0.95` has true MAF 0.05), and a MAF-stratification `cut()` on
unfolded `pf` (breaks assuming `[0,0.5]`) silently produced `NA` for roughly
half of all markers. Fixed in both scripts on 2026-09-13.

- Primary parental-MAF gate: `pmaf >= 0.15`
- Retained units after the fix: **395,996** (previously, with the unfolded
  bug: 636,649 -- do not use the 636,649 figure or anything derived from it)
- Per-DI-bin retained counts (primary, MAF-gated, from the violin script's
  log, identical unit set as `di25_fst_vs_di_rho05.rds`'s primary gate):

  | DI_bin | N |
  |---|---|
  | (-Inf,-90] | 9,122 |
  | (-90,-75] | 23,178 |
  | (-75,-60] | 107,620 |
  | (-60,-50] | 88,587 |
  | (-50,-40] | 86,117 |
  | (-40,-30] | 50,713 |
  | (-30,-25] | 13,290 |
  | (-25,-20] | 8,412 |
  | (-20,-15] | 4,719 |
  | (-15, Inf] | 4,238 |

  (Sums to 395,996.)

- Overall Fst with the corrected gate: 0.245 (MAF-gated primary) vs 0.248
  (ungated). The corrected per-bin pattern and medians are "almost
  unchanged" from the pre-fix (buggy-gate) run -- the fix is a correctness
  fix to the unit *count* and *identity*, not a change to the qualitative
  DI-Fst relationship.
- `pmaf` (folded) is carried alongside `best`/`DI`/`fst_locus` in the units
  table saved inside `data/di25_fst_vs_di_violin_rho05.rds$units` (see §6),
  which is the authoritative per-unit table to join against.

## 3. Authoritative neutral-simulation Fst summaries

**Three distinct neutral-simulation Fst quantities exist across two files.
They are not interchangeable -- pick the one matching your question.**

### (a) Pooled, one value per replicate (genome-wide multilocus pool)
- Location: `data/di25_fst_vs_di_rho05.rds$neutral` (a numeric vector, 1000
  values, one per simulation replicate)
- Definition: for each of the 1000 `diem_boot%d_output.bed` replicates, Fst
  is computed as `sum(a)/sum(abc)` (Weir & Cockerham) pooled across all loci
  in that replicate -- a single genome-wide value per replicate.
- Summary: median 0.0213, 95% interval [0.0187, 0.0248] (2.5/97.5
  percentiles across the 1000 replicates).
- Use case: a single reference value/interval for "neutral genome-wide Fst
  under this demographic model," not appropriate for any per-locus or
  per-DI-bin comparison.

### (b) Pooled, per-DI-bin-per-replicate
- Location: `data/di25_fst_vs_di_rho05.rds$env` (a `data.table`, 10 rows,
  one per DI bin), columns `sim_med`/`sim_lo`/`sim_hi` (plus `bin`,
  `DI_bin`, `emp`, `emp_ungated`, `n_emp_units`, `n_emp_units_ungated`,
  `n_sim_units`)
- Definition: within each DI bin, Fst pooled across only the simulated loci
  that fall in that bin, per replicate, then summarized (median/95%
  interval) across the 1000 replicates.
- **Only populated for DI bins 6-10** (the "ancestry-informative" bins where
  the simulation panel has overlap with real markers -- see §7's
  ascertainment note); bins 1-5 have `sim_med`/`sim_lo`/`sim_hi` = `NA`.
- Use case: the DI-bin-matched neutral reference actually plotted alongside
  the empirical DI-bin Fst in `Figures/di25_fst_vs_di_rho05.{png,pdf}`.

### (c) Per-locus, unpooled (the only one usable for a locus-level comparison)
- Location: `data/di25_fst_vs_di_violin_rho05.rds$simO_locus` (a numeric
  vector, **16,615,717 values** -- every (locus, replicate) pair pooled
  together, non-finite values excluded, 0 excluded)
- Also: `data/di25_fst_vs_di_violin_rho05.rds$n_sim_reps` (1000),
  `$n_sim_overlap_units` (16,616 -- the number of *distinct* simulated
  marker positions, before multiplying by replicates)
- Definition: `a/abc` (Weir & Cockerham numerator/denominator) computed
  **independently per locus**, never pooled/averaged across loci within a
  replicate. This is the genuine locus-level analogue of the empirical
  `fst_locus` column.
- Summary: median 0.0188 (compare to (a)'s pooled-per-replicate median of
  0.0213 -- the difference itself demonstrates why (a)/(b) are NOT
  substitutable for a per-locus comparison: pooling across loci narrows and
  shifts the distribution).
- Built by script `R/di25_fst_vs_di_violin_rho05.R`, cached per-replicate in
  `data/fst_sim_locus_cache_rho05/` (1000 files; fingerprint keyed on
  marker_hash + params_hash + `script_version =
  "di25_fst_vs_di_violin_rho05_locus_v1"`).
- **Use this one** for any figure or test comparing empirical per-locus Fst
  (stratified by sort_class, prop_fixed, DI, etc.) against a neutral
  reference distribution at the same (per-locus) resolution.

Do **not** use `fst_sim_cache_full_rho05_v2/` (the old pooled-per-replicate
cache backing (a)/(b)) as a stand-in for (c) -- it contains no per-locus
values at all, only pre-pooled per-replicate/per-bin numbers.

## 4. Authoritative empirical sorting summaries

**Two independent sorting classifications exist, over two different unit
universes, using two different population sets. They are not
interchangeable.**

### (a) Stage-1-direct tau sweep (BayPass/Omega climate pipeline universe)
- Files: `moduleA_stage1_cluster_sorting_tau05.rds`,
  `moduleA_stage1_cluster_sorting_tau06.rds`,
  `moduleA_stage1_cluster_sorting_tau08.rds`
- Unit universe: **18,361 units**, Stage-1-direct clusters only
- Population set: **Aland EXCLUDED** (19 hybrid pops) -- the BayPass/Omega
  convention (Issue 5 in `AUDIT_FIXES.md`)
- Genuine tau sweep at 0.5 / 0.6 / 0.8 (all mtime 2026-09-12 18:01):

  | tau | directional | unresolved | differentiated (total) |
  |---|---|---|---|
  | 0.5 | 2,961 | 210 | 10,170 |
  | 0.6 | 1,945 | 17 | 10,170 |
  | 0.8 | 384 | 0 | 10,170 |

- Use case: any question tied to the Stage-1/Stage-2 climate-candidate
  pipeline and its Omega/BayPass population universe.

### (b) Full LD-reduced-unit sorting classification (Fst-vs-sort_class figure)
- Script: [R/di25_fst_vs_di_sorting_stratified_rho05.R](R/di25_fst_vs_di_sorting_stratified_rho05.R)
- Saved object: `data/di25_fst_vs_di_sorting_stratified_rho05.rds`
  (fields: `units`, `di_high_cutoff`, `sort_params`)
- Figures: `Figures/di25_fst_vs_di_sorting_stratified_rho05.{png,pdf}`
- Unit universe: **395,996 units** (the same MAF-gated set as §2/§3c) --
  computed fresh via `parallelism_stats()`/`classify_sort()`
  (`moduleA_sorting/R/parallelism_stats.R`), not reused from (a)
- Population set: **all 20 hybrid pops, Aland included** -- a pure
  Fst/sorting question, independent of the Omega/BayPass population
  convention, matching `di25_fst_vs_di_rho05.R`'s own population handling
- **tau=0.6 only** (`SORT_TH=0.6`), no sweep; other locked parameters:
  `MIN_PARENT_MAF=0.15`, `FIX_TH=0.15`, `SORT_RULE="binom"`, `ALPHA=0.05`
- Full-dataset (395,996 units) `sort_label` counts:

  | sort_label | N |
  |---|---|
  | differentiated, not fixed (= `classify_sort()`'s "unsorted") | 292,553 |
  | sorted: polyctena | 57,241 |
  | sorted: aquilonia | 41,502 |
  | differentiated, unclassified (NA `uni_score`/`n_obs==0`) | 4,152 |
  | differentiated, unresolved | 548 |
  | not differentiated | (remainder) |

  Note the label mapping: `classify_sort()`'s raw category "unsorted" (its
  literal name; differentiated and observed but below the fixation
  threshold, and the *largest* category, not an edge case) is relabeled
  here to "differentiated, not fixed" for figure clarity; "ambiguous" is
  relabeled "differentiated, low power". See the script's own bugfix
  comment (lines 87-92) -- an earlier render silently dropped "unsorted" to
  `NA` because it was missing from the label-ordering vector; fixed with an
  explicit `stopifnot` guard against any future unmapped category.
- High-DI subset used in the figure: `DI > -25`, restricted to
  `is.finite(fst_locus)`. Per-sort_label medians in that high-DI subset:
  "differentiated, not fixed" (n=16,080) median Fst ~0.25; "sorted:
  aquilonia"/"sorted: polyctena" (n=914/339) median ~0.28-0.29;
  "differentiated, unresolved" (n=36) median ~0.51 (highest, small n).
- Use case: any question relating per-locus Fst magnitude to sorting
  class/prop_fixed across the full rho05 unit set -- this is the object to
  join against for that purpose, not (a).

## 5. Marker/unit definitions, rho, LD-reduction, and representative-SNP convention

- LD-pruning source object: `module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05_bestsnp.rds`
  (fields: `stats`, `geno`, `rep_snp_all`)
- Build script: [module0_ld_pruning_rho05/R/module0_ld_pruning_rho05_DIEM.R](../module0_ld_pruning_rho05/R/module0_ld_pruning_rho05_DIEM.R)
- Parameters (as run, `ld_prune_and_eMLG(...)`, line ~49-51 of the build
  script): `ld_w_threshold = 0.025`, `score_threshold = 0.80`, `min_r2 =
  NULL`, **`min_r2_rho = 0.5`** (the decay-relative LD-pruning join
  threshold this module family is named for), `min_n_loci_flag = 5`,
  `cM_threshold = 0.5`. All non-rho settings held at
  `module0_ld_pruning`'s canonical values; only the rho-based join
  threshold differs from that canonical (fixed `min_r2=0.2`) run --
  see the repo-wide "rho=0.5 standing decision" convention.
- Total LD-reduced units: `nrow(rep_snp_all) = 661,386`
- Large (`has_eMLG==TRUE`) clusters, each with its own local best-SNP call:
  `nrow(stats) = 17,509`; confirmed `sum(rep_snp_all$has_eMLG) == 17,509`
  (exact match)
- **Representative-SNP convention** (`rep_snp_all$rep_snp`): for the 17,509
  large clusters, `rep_snp = best_marker` (from `eMLG_best_snp()`, chosen
  within-cluster); for the remaining 643,877 small/singleton clusters,
  `rep_snp = representative` (the cluster's centrality representative, not
  an eMLG best-SNP call, since eMLG requires >= `min_n_loci_flag` = 5 loci).
  Verified: `rep_snp == best_marker` exactly when `has_eMLG == TRUE`,
  `rep_snp == representative` exactly when `has_eMLG == FALSE`; no `NA`s,
  no duplicates, across all 661,386 rows.
- All downstream rho05 scripts (§1/§2/§4b) use `rep_snp_all$rep_snp` as
  "the unit's marker" (column `best` in their own tables) -- this is the
  single definition of "unit" and "marker" throughout
  `module_manuscript_rho05/`.

## 6. Paths to source-data objects underlying the relevant figures

| Object | Path | Size | mtime | MD5 |
|---|---|---|---|---|
| DI vs per-locus Fst summary | `data/di25_fst_vs_di_rho05.rds` | 9,075 B | 2026-09-13 14:16:34 | `143e887a75ce3b7bf383a2b318776ada` |
| Per-locus Fst violin (per-unit rows + per-locus sim) | `data/di25_fst_vs_di_violin_rho05.rds` | 127,381,797 B | 2026-09-13 14:29:20 | `00540941488f03eed786563f4dce0b8e` |
| Fst-vs-sorting-class figure data | `data/di25_fst_vs_di_sorting_stratified_rho05.rds` | 10,938,608 B | 2026-09-13 14:38:55 | `a67785903e3655c486a4e7b8d8d83a0b` |
| LD-reduced unit/marker source | `module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05_bestsnp.rds` | 9,006,128 B | 2026-09-07 15:10 | `29bcdfffac6fdc2f11f87940611be22a` |
| Per-locus sim cache (1000 files) | `data/fst_sim_locus_cache_rho05/` | -- | -- | (per-file fingerprinted, not a single checksum) |
| Pooled sim cache (1000 files, for §3a/3b only) | `data/fst_sim_cache_full_rho05_v2/` | -- | -- | (per-file fingerprinted) |
| Stage-1-direct tau sweep (§4a) | `moduleA_stage1_cluster_sorting_{tau05,tau06,tau08}.rds` | -- | 2026-09-12 18:01 | (not individually hashed for this audit) |

`data/di25_fst_vs_di_violin_rho05.rds$units` is the authoritative per-unit
table (395,996 rows) to join against for any new analysis needing
`best`/`DI`/`bin`/`pmaf`/`fst_locus` together; `data/
di25_fst_vs_di_sorting_stratified_rho05.rds$units` is the same 395,996-row
table with `differentiated`/`prop_fixed`/`uni_score`/`sort_class`/
`sort_label` additionally joined on.

## 7. Known limitations

- **Simulation-panel ascertainment is not genome-wide random.** The DI25
  simulation panel (`data/diem_outs_demo/diem_boot%d_output.bed`, 1000
  replicates) overlaps only **16,616 of 661,386** empirical units (~2.5%),
  and that overlap is concentrated in DI bins 6-10 (the
  "ancestry-informative" bins) -- it is not a random genome-wide sample of
  markers. This is why `di25_fst_vs_di_rho05.rds$env`'s per-bin simulation
  summary (§3b) is populated only for bins 6-10, and why any per-locus
  empirical-vs-simulated comparison (§3c) implicitly compares against a
  DI-biased marker subset, not a genome-wide neutral marker set.
- **DI-vintage/gating mismatch risk for cross-module comparison.** DI
  throughout `module_manuscript_rho05/` (including this handoff's §1/§2/§4b
  analyses) comes solely from `map_hyb_005$DiagnosticIndex` (full genome,
  never pre-gated); the `DI > -25` cut in §4b is a fresh cut on that
  never-restricted table. `module_population_partitioning`, by contrast,
  uses `module_di25`'s own DI25 objects, which are **already pre-gated** to
  `DI > -25` upstream. A naive side-by-side comparison of "high-DI" results
  between the two modules could show a spurious difference driven by this
  DI-vintage/gating mismatch rather than a real biological one -- see the
  cross-module note in `R/di25_fst_vs_di_sorting_stratified_rho05.R` (lines
  21-32) for the full reasoning, contributed after cross-session
  consultation on 2026-09-13. This does not affect the internal correctness
  of any object listed in this handoff (each uses one consistent DI source
  throughout), only the validity of comparing across modules.
- **No independent MAF-fold-fix intermediate quarantine exists.** See §9 --
  the pre-fix, post-audit intermediate states were overwritten in place,
  not separately archived.

## 8. Checksums and modification times for the authoritative files

Scripts:

| File | MD5 | mtime | Size |
|---|---|---|---|
| `R/di25_fst_vs_di_rho05.R` | `5cc9e23d68e07dcd9621e95183908ff9` | 2026-09-13 14:08:01 | 21,410 B |
| `R/di25_fst_vs_di_violin_rho05.R` | `147bff7c980566ccc5b0bcc3d52f60e8` | 2026-09-13 14:28:29 | 12,932 B |
| `R/di25_fst_vs_di_sorting_stratified_rho05.R` | `daf4c0eed5d069657645c4cb421633cb` | 2026-09-13 14:53:35 | 9,166 B |

Saved objects: see the table in §6 (identical MD5/mtime/size values).

Source clustering object: see the table in §6 (identical MD5/mtime/size
values for `eMLG_5loci_0025_cM05_rho05_bestsnp.rds`).

Git HEAD for all of the above: `d0a3da7e1fbb5716dd0c8d1da60d142a687d810f`
(2026-09-13 14:53:45 +0300).

## 9. Stale or superseded objects -- do not use

- `stale_pre_fix_20260912/data/di25_fst_vs_di_rho05.rds` (8,517 B, mtime
  2026-09-12 18:41) -- the ORIGINAL pre-audit, pre-MAF-fold version, built
  on the 17,509-unit (large-cluster-only, pre-LD-reduction-to-661,386)
  universe. Do not use for any purpose; superseded by §1/§6.
- `stale_pre_fix_20260912/data/fst_sim_cache_full_rho05_OLD/` (1000 files)
  -- the original simulation cache built against a marker set with only
  1,511 overlapping markers (not the current 16,616-overlap panel). Do not
  use; superseded by `data/fst_sim_cache_full_rho05_v2/` (§3a/b) and
  `data/fst_sim_locus_cache_rho05/` (§3c).
- **Gap in the audit trail (disclosed, not hidden):** there is no
  separately-quarantined intermediate copy between (i) the
  Issue-2-audit-fix state (636,649 MAF-gated units, MAF unfolded) and (ii)
  today's MAF-fold-fix state (395,996 units, MAF folded) for
  `di25_fst_vs_di_rho05.rds` -- state (i) was silently overwritten in place
  when the script was rerun today, without an explicit `cp`-to-quarantine
  step first. Similarly, **no stale copy of `di25_fst_vs_di_violin_rho05.rds`
  exists at all**: the pre-MAF-fold, pooled-per-replicate-`simO` version of
  that object was overwritten in place by the same rerun. Only the
  ORIGINAL pre-audit (17,509-unit, 1,511-marker-sim-overlap) versions
  listed above survive as quarantined stale artifacts. If any historical
  reproduction of the intermediate (636,649-unit, unfolded-MAF) numbers is
  ever needed, it is not recoverable from disk and would have to be
  re-derived by temporarily reverting the MAF-fold fix in git history.
- `moduleA_stage1_cluster_sorting_*.rds` (§4a) are NOT stale -- they remain
  the authoritative sorting summaries for the Stage-1-direct/BayPass-Omega
  universe. They are simply a *different* universe from
  `di25_fst_vs_di_sorting_stratified_rho05.rds` (§4b) and must not be
  substituted for it, or vice versa.
