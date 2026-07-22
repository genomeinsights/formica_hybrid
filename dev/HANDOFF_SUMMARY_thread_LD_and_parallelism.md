# Handoff summary — LD-clustering / eMLG / parallelism-sorting thread

Written for merging with the "predictability of sorting" thread into one
pipeline. Distilled conclusions and validated code, not a conversation
replay. Companion doc: `dev/methods_notes.md` (LD-pruning/eMLG design
rationale specifically — not duplicated here, just referenced).

## Critical, time-sensitive issue

**`data/eMLG_5loci_0025.rds` was silently overwritten in place on
2026-07-21 01:05** by the production run switching to the cM=0.5 backstop
(`params$use_cM=TRUE`, `cM_threshold=0.5`, 474,014 groups — up from
469,361 in the old bp/rho=0.95-backstop file, no `params` field). It is
not versioned; the old file's contents no longer exist anywhere.

**Every cluster-level finding below was computed against the OLD file**
and has NOT been re-verified against the new one — group IDs are not
stable across runs (the same string, e.g. `"F10329"`, refers to a
different marker set with a different `n_loci` in each), so nothing that
joins on `group_id` from before 01:05 to after it will work, and no
error will necessarily warn you — it'll just silently give you wrong
answers (this is exactly what happened here first).

Checked by file mtime: `dev/R/tier1_architecture.R`,
`tier2_sorting_vs_recomb.R`, `supplementary_analyses.R` were all last
touched 2026-07-20 17:39–18:46, **before** the overwrite — worth
double-checking whether their numeric results are still current.
`sorting_baypass_crosscheck.R` / `_memberlevel.R` were touched
2026-07-21 10:10/10:20, **after** — likely safer, but mtime isn't proof
of re-execution, so worth confirming they were actually *run* against
the new file rather than just edited.

**Recommendation for the merged pipeline**: version-tag cluster outputs
(e.g. `eMLG_<params-hash-or-date>.rds`) instead of overwriting a fixed
filename, so this can't recur.

## Validated, reusable code

### LDscnR (`~/gitlab/LDscnR`, package — `devtools::load_all()`)
- `ld_complexity_reduction()` (Stage 1) / `ld_prune_and_eMLG()` (Stage 2,
  `dynamic_cut_eMLG()`-based) — pre-existing, now also takes
  `genetic_map`/`cM_threshold` for cM-based (not flat-bp) distance
  restriction. New helper `interpolate_cM()`. Commit `ed8391d`.
- `plot_genotype_heatmap()` — new, generic ComplexHeatmap genotype
  dosage plot: arbitrary named-vector column/row annotations (not
  hardcoded to CL_id/Population), explicit `row_order`, optional
  `col_annotation_legend`/`row_annotation_legend` toggles (large cluster
  counts make a full legend unreadable), `polarize` toggle. Commits
  `2562bc1`, `4ad44b5`. 109/109 package tests pass.
- Unexported helpers used throughout: `polarize_genotypes()`,
  `consensus_dosage()` / `expected_gt_dosage()`, `score_eMLG()`
  (`cor(round(x), x)^2` — round-trip fidelity of a consensus dosage
  against its hard-called version; **not yet applied** to any of the
  parallelism-test consensus below, see open questions).

### formica_hybrid (`~/gitlab/formica_hybrid`)
- `R/ld_pruning_DIEM.R` — production pruning/eMLG script, cM=0.5 backstop
  (commit `76fdac7`, ran overnight → the overwrite above). Edited again
  2026-07-21 09:55, after I last reviewed it — not something I've seen.
- `R/LD_decay_from_DIEM.R` — added parent-inclusive
  `GTs_with_parents`/`sample_data_with_parents` construction (same
  markers as `map_hyb_005`, same order), saved to
  `data/hybrids_and_parents_maf005.Rdata`. Commit `ffb8c7a`. Needed
  because both the parallelism test and (per the other thread's file
  comments) the BayPass cross-check both require true parental
  genotypes for ancestry orientation, which the hybrids-only pipeline
  deliberately excludes (parent-vs-hybrid structure would badly distort
  LD estimates for everything else).
- `dev/R/parallelism_stats.R` — pre-existing (this repo's own code, not
  written by me). Binomial parallelism test: per locus, classifies
  `sort_class` (`aquilonia`/`polyctena`/`bidirectional`/`unsorted`) via
  `uni_score`/`bi_score` crossing `sort_th` (default 0.5) —
  **`sort_class` is NOT gated by the binomial test's own significance**
  (`p_binom`/`q_binom` are computed but not used to assign the label;
  see open questions).
- `dev/R/eMLG_parallelism.R` — **new this session, uncommitted**.
  - `consensus_dosage_matched(GTs_ref, GTs_new, mk)`: consensus dosage
    for marker set `mk`, using the reference-marker/flip decisions
    derived from `GTs_ref` but applied to `GTs_new`. Needed because
    `polarize_genotypes()` picks its reference/flip fresh from whatever
    matrix it's given — computing it separately on parents (30
    individuals, strongly differentiated) vs. hybrids (164) can flip
    different markers, silently putting the two consensus vectors on
    opposite orientations. Strictly generalizes LDscnR's own
    `consensus_dosage()` (identical when `GTs_ref` and `GTs_new` are the
    same data).
  - `build_group_consensus()`: `mclapply` wrapper, many groups at once.
  - `cluster_level_parallelism()`: complete, gap-free cluster-level
    parallelism test for a marker set — singleton clusters passthrough
    their own SNP-level result (nothing to average), small (2-4 loci)
    clusters get a freshly-built matched-flip consensus (both hybrid and
    parent sides), large (>=5 loci, already in `eMLG_5loci_0025$eMLG`)
    clusters reuse the existing hybrid-side eMLG + a freshly-built
    parent-side consensus. `di_agg = c("representative", "max")` controls
    whether a cluster's DI (for the `min_DI` gate) comes from its
    representative marker only or the max DI across all its members.
    Validated against independently-computed manual numbers (exact
    match). **Built against the OLD (pre-01:05) `eMLG_5loci_0025.rds`** —
    needs re-running against the new file before trusting it further.

## Key findings (⚠️ against the now-superseded OLD clustering)

- **Per-SNP sort_class is heavily pseudo-replicated.** 570 SNPs
  individually `sort_class=="bidirectional"` partition into only 315
  distinct Stage-2 clusters. Actually re-testing each cluster's consensus
  (not just noting it contains a flagged member) found 226/315 (71.7%)
  show **no** sorting evidence at all once aggregated; only 1 cluster
  (`F59859`, Chr24, 257 markers) had strong multi-marker support; the
  other 24 "still bidirectional" clusters are singletons or 2-marker
  clusters with no better evidence than the raw per-SNP noise.
- **The DI gate via "representative marker only" spuriously excludes
  clusters.** 60/315 clusters failed `differentiated==FALSE` (representative
  marker's DI < -50) despite `n_obs=20` (plenty of data) — not an
  insufficient-data problem, purely a DI-lookup-convention artifact.
  Switching to `di_agg="max"` (best DI across all members) resolved all
  60 — 59 to "unsorted", 1 to "polyctena". **Bidirectional count was
  unchanged (still 25)** — broadening the DI gate didn't hide or
  manufacture any bidirectional signal, only reshuffled
  unidirectional/unsorted classifications.
- **Genome-wide** (all 265,046 individually-sorted SNPs — aquilonia +
  polyctena + bidirectional): partitioning + re-testing collapses this to
  151,335 independent units (93,231 singleton / 39,580 small / 18,524
  large). Aquilonia consistently outnumbers polyctena ~1.8-2:1 across
  every version tested (58-61% vs. 31-32%) — likely reflects genome-wide
  admixture skew (this repo's own `parallelism_stats.R` docs flag this;
  `null_prob="pooled"` would correct for it but **hasn't been tried**).
  Bidirectional stayed rare and stable (max 25 units) regardless of
  DI-gating convention.
- **Not yet run**: the `score_eMLG()` dilution check — was about to test
  whether the 226 "washed out" clusters are homogeneous (genuinely no
  signal) or heterogeneous (real signal diluted by averaging in less-
  informative members) when the production data changed underneath the
  analysis. This connects directly to the outlier-region framing in the
  other thread (`sorting_baypass_crosscheck*.R`): a cluster/region can be
  a false or true positive depending on how strongly its members
  correlate with the actual causal signal — exactly what `score_eMLG()`
  is a cheap proxy for.

## Related, pre-existing work in this repo (commit history only — I don't
have full narrative context, the other thread likely does)

- `R/diagnostic_index_enrichment.R` (commit `d0d12f2`): tests whether
  BayPass PC1/PC2 outlier LD clusters (>=10 individually significant,
  BF(dB)>=20 members) are enriched for high-DiagnosticIndex loci vs.
  background, correcting for the mechanical confound that the outlier
  definition requires large clusters and cluster size independently
  predicts DI. PC1 enrichment survives the size-adjusted test in all 4
  population/Omega settings; PC2 splits (2 robust, 2 reverse below 1).
- `R/analyse_baypass.R`, `582ad15` (cluster-colored Manhattan plots),
  `fedbfc3` (Sielva-excluded configs) — the wider BayPass PC1/PC2 outlier
  pipeline this enrichment test sits on top of.

## Open questions for the merged pipeline

1. Re-run `cluster_level_parallelism()` (and anything else keyed on
   `eMLG_5loci_0025$groups`) against the new (post-01:05) production
   file — is `F59859` even the same cluster under cM-based backstop?
2. Run the `score_eMLG()` dilution check once re-run on current data.
3. Try `null_prob="pooled"` to see how much of the aquilonia:polyctena
   skew survives correcting for genome-wide admixture asymmetry.
4. Decide on a version-tagging convention for `eMLG_*.rds` so this
   silent-overwrite problem can't recur across two concurrent sessions.
5. Reconcile with the other thread's sorting-vs-BayPass-outlier
   cross-check and DI enrichment work — same "is a cluster/region real or
   an artifact of how strongly its members tag the causal signal"
   question, approached from two directions (parallelism-sorting here,
   climate-association outliers there).
