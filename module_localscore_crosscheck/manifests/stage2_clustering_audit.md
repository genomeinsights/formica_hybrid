# Stage-2 rho05 clustering audit

Status: audit of the canonical Stage-2 rho05 LD-clustering object, recorded
before any downstream use (per instruction: record its params/checksum/size
distribution in a manifest rather than assuming them).

## Source

- Path: `module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05.rds`
- MD5: `3b685f6ad972be13d1e0d7de6a3fb7f3`
- Size: 15,633,340 bytes (~15.6 MB), last built 2026-09-07 14:53

## Params (`$params`, as stored)

| Param | Value |
|---|---|
| `ld_w_col` | `ld_w_095` |
| `ld_w_threshold` | 0.025 |
| `min_n_loci_flag` | 5 |
| `rho` (LD-clustering rho) | 0.95 |
| `distance_threshold` | NULL |
| `min_r2_rho` | 0.5 |
| `min_r2` | NULL (derived, see `min_r2_resolved`) |
| `use_cM` | TRUE |
| `cM_threshold` | 0.5 |
| `min_r2_resolved` | ~0.534 per chromosome (0.531 on Chr26 only) |

"rho05" in the object name refers to `min_r2_rho = 0.5` (the r^2-threshold
derivation rho), not the LD-clustering `rho = 0.95` -- these are two
different parameters, both present, and both recorded here rather than
conflated.

## Structure

`readRDS()` returns a list with 4 elements: `eMLG`, `groups`, `pruned`, `params`.

- `$groups`: a `data.table`, 661,386 rows, one per LD cluster at any size.
  Columns: `group_id, Chr, representative, n_loci, score, has_eMLG, members`.
- `$pruned`: character vector, length 661,386 -- one pruned/representative
  marker ID per cluster (matches `nrow($groups)`).
- `$eMLG`: numeric matrix, 165 rows x 17,509 columns -- the eMLG consensus
  genotype matrix, one column per cluster with `has_eMLG == TRUE`
  (`n_loci >= 5`), one row per individual.

## Partition property (verified directly, not assumed)

The `groups$members` lists are **non-overlapping**: total member entries
across all 661,386 clusters = 1,114,423, and the number of *unique* marker
IDs across all `members` lists is also 1,114,423. Every marker belongs to
exactly one Stage-2 LD cluster -- `$groups` is a genuine partition of the
full marker set at this clustering resolution.

This is a different property from what "Stage-2 clusters are non-disjoint /
must not be forced into a partition" (as specified) is actually about:
clusters do not overlap *each other* here. The relevant non-partition
relationship is between **local-score windows and Stage-2 clusters**: a
window can span markers belonging to several different Stage-2 clusters,
and a single Stage-2 cluster's member markers are not guaranteed to fall
entirely inside one window (LD/genetic-distance-based clustering is not
window-bounded). Downstream per-window annotation (Part B) will therefore
report "unique Stage-2 cluster IDs touching a window" and "distinct Stage-1
best-marker positions inside a window" as separate counts, per instruction
-- not because clusters overlap, but because the window<->cluster mapping
is many-to-many.

## `has_eMLG` / Stage-2-unit subset (`n_loci >= 5`)

- `has_eMLG` agrees with `n_loci >= 5` in every one of the 661,386 rows
  (cross-tabulated: 0 disagreements).
- 17,509 clusters qualify (2.6% of all 661,386 clusters), covering
  `sum(n_loci)` = 1,114,423 - (484,866 singletons) - ... (see size
  distribution below for the bulk of coverage).

### Cluster size (`n_loci`) distribution

| Statistic | All 661,386 clusters | Stage-2 subset (n_loci>=5, N=17,509) |
|---|---|---|
| Min | 1 | 5 |
| Median | 1 | 6 |
| 90th pct | 2 | 16 |
| 95th pct | 3 | 30 |
| 99th pct | 7 | 125 |
| 99.9th pct | 39 | 1302 |
| Max | 5382 | 5382 |
| N singletons (n_loci==1) | 484,866 (73.3%) | -- |

The single largest Stage-2 cluster is `F241709` on Chr25, n_loci = 5382 --
flagged here for later cross-reference against the recombination map in
Part C (candidate centromeric/low-recombination region).

### Chromosomal distribution of Stage-2 (n_loci>=5) clusters

Chr2 (1289), Chr1 (1082), Chr4 (1021), Chr9 (785), Chr3 (877), Chr5 (933),
Chr12 (819), Chr7 (790), Chr6 (832), Chr8 (661), Chr11 (704), Chr14 (645),
Chr15 (669), Chr10 (528), Chr16 (581), Chr21 (581), Chr13 (849), Chr25 (537),
Chr9 (785), Chr20 (448), Chr17 (501), Chr22 (494), Chr18 (443), Chr24 (398),
Chr19 (377), Chr26 (302), Chr27 (363). No chromosome has zero Stage-2
clusters; no gross imbalance beyond what chromosome length would predict
(not further tested here -- a proper length-normalized check belongs in
Part C's recombination-decile analysis, not this basic audit).

## `representative` vs. `best_marker` -- confirmed distinct fields

- `groups$representative` is chosen **at clustering time**, purely on LD
  centrality (highest median r^2 to the rest of the cluster) -- it is not
  informed by genotype-consensus fidelity or any association statistic.
- `best_marker` is a **separate, post-hoc** field, computed only for the
  17,509 `has_eMLG == TRUE` clusters by `eMLG_best_snp()`
  (`~/gitlab/LDscnR/R/eMLG_best_snp.R`), defined as the cluster member SNP
  whose own genotype has the highest `|r|` (Pearson, pairwise-complete) to
  the cluster's eMLG consensus genotype (sign-corrected via a `flipped`
  flag). This criterion uses **no BF, p-value, F_ST, DI, or any
  association/differentiation statistic** -- it is genotype-fidelity-only,
  which is exactly the "representative selection independent of
  association strength" property required for the Part E LD-reduced
  sensitivity analysis.
- For clusters with `has_eMLG == FALSE` (n_loci < 5, 643,877 of 661,386),
  `best_marker` is undefined; the companion object's `rep_snp_all$rep_snp`
  column falls back to `representative` for these, which is the
  centrality-representative marker, not a consensus-fidelity best marker
  (there being no consensus to be faithful to at that size).
- Verified concretely: cluster `U2` (`Chr1:86227`, n_loci=2, has_eMLG=FALSE)
  has `best_marker = NA` in the companion object -- `representative` and
  `best_marker` coincide there only via the fallback, not because the two
  selection criteria agree.

## Pre-computed companion object (source, reused as an authoritative input)

- Path: `module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05_bestsnp.rds`
- Built by: `module0_ld_pruning_rho05/R/module0_ld_pruning_rho05_DIEM.R`,
  calling `eMLG_best_snp()` directly on the object audited above.
- Contents: `$stats` (17,509 rows: `group_id, representative, best_marker,
  n_loci, score, best_r, best_abs_r, rep_abs_r, rep_is_best, flipped, n_obs,
  n_filled, n_resid_na`), `$geno` (165 individuals x 17,509 best-SNP-called
  genotype matrix), and `$rep_snp_all` (661,386 rows, `best_marker` filled
  where available, else falling back to `representative` -- see above).
- This is a deterministic function of already-authoritative upstream inputs
  (the Stage-2 rho05 clustering + genotype matrix), not an artifact of the
  exploratory local-score module -- it is used here as an authoritative
  input, consistent with "use the original BayPass inputs ... as
  authoritative inputs," not copied from
  `module_localscore_crosscheck_exploration/`.

## Distinct from the Stage-1-direct `best_marker` object

A second, unrelated `best_marker` object exists:
`module_manuscript_rho05/data/moduleB_stage1_units_bestsnp.rds`, built by
calling `eMLG_best_snp()` on the **Stage-1** clustering
(`module0_ld_pruning/data/pruned_stage1.rds`, `n_snps >= 5`, 18,361 units),
bypassing Stage-2 merging entirely (documented rationale in that script:
Stage-2 over-merging costs power for outlier scans). This is the object the
existing (non-exploratory) `moduleA_stage1_cluster_sorting.R` and
`moduleC_stage1_annotations.R` scripts consume via
`obj$best$stats$best_marker`. **This is the object the new Stage-1-
resolution scripts in this module must use** for Stage-1-cluster genomic
Position assignment -- not the Stage-2 rho05 `best_marker` audited above,
which is reserved for the separate "Additional objective" Stage-2
redundancy/susceptibility analysis (Parts A-J).
