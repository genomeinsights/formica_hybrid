# Cross-module input provenance audit

Written 2026-09-13, before any follow-up work on `module_population_partitioning`,
per repository convention: this project is split across `module_manuscript_rho05/`
(authoritative for current rho05 manuscript-run FST/DI/sorting results) and
`module_population_partitioning/` (this module, population-partition
concordance), and neither is self-contained. This document traces every
external input this module currently reads, plus the newer rho05-corrected
counterparts a reader might otherwise reach for by mistake, so the two are
never silently conflated.

> **UPDATE 2026-09-13 (same day, after the user's migration decision):** this
> module was migrated from the legacy `module_di25` (`min_r2=0.2`) DI25
> clustering to the `module_di25_rho05` (`min_r2_rho=0.5`) one — see
> `../README.md` "Migration to the rho05 primary dataset" and "Two analysis
> universes". **Table A below now reflects what this module currently reads
> (rho05); the previously-primary legacy objects are demoted to Table A-legacy.**
> The full-genome `module_manuscript_rho05` universe (Tables C/D) remains a
> deliberately separate, secondary universe this module does not read — per
> the user's explicit decision, NOT migrated onto wholesale.

**Checksums are SHA-256** (`shasum -a 256`), captured 2026-09-13. Sizes in
bytes. All paths repo-root-relative.

## Naming warning: two unrelated parameters are both called "rho"

Before anything else, because this trips up every table below: the pipeline
uses "rho" for two different things.

1. `RHO_LDW = 0.95` — the `compute_ld_w()` neighbourhood threshold. Identical
   in every clustering object below (legacy and rho05 alike); not what the
   "_rho05" module-name suffix refers to.
2. `min_r2_rho = 0.5` — the Stage-2 (`ld_prune_and_eMLG`) merge-quality gate,
   a **decay-relative** r² threshold. This is what "_rho05" means. Since the
   2026-09-13 migration, the objects this module currently reads (Table A)
   DO use this convention; the now-superseded objects (Table A-legacy) used
   a **fixed** `min_r2 = 0.2` instead (not rho-based at all, and not even
   stored as a named field in those objects' own `$params`).

So "rho05" in a module or file name refers exclusively to parameter (2). Do
not infer anything about parameter (1) from it.

## Table A — objects `module_population_partitioning` currently reads (as of the 2026-09-13 rho05 migration)

Same DI25 marker panel, same Stage-1 partition (reused unchanged from the
legacy lineage), same cM=5 cap, same `MIN_N_LOCI_EMLG=3`, same best-SNP
`fill=FALSE` (strictly observed, non-imputed) convention as the superseded
objects (Table A-legacy) — **the only thing that differs is the Stage-2
gate**: `min_r2_rho=0.5` (decay-relative; resolves to ≈0.535 per chromosome,
see `$params$min_r2_resolved`) instead of the fixed `min_r2=0.2`. Confirmed
by reading `module_di25_rho05/R/di25_ld_clustering_rho05.R` and
`di25_sorting_rho05.R` directly (both state this explicitly as a "paired
comparison, min_r2 is the only thing that differs").

| path | producing script | module | dims / key columns | mtime | sha256 (first 12) | role | status |
|---|---|---|---|---|---|---|---|
| `module_di25/data/di25_inputs.rds` | `module_di25/R/di25_ld_clustering.R` | module_di25 | list(map: 51,612×3 [Chr,Pos,marker], GTs_hyb: 165×51,612, GTs_par: 30×51,612) | 2026-08-07 16:29 | `330d270af27e` | the frozen DI>−25 marker panel + hybrid/parent genotypes; UNCHANGED by the rho05 migration — only the Stage-2 clustering differs | **authoritative** — same object used by both the legacy and rho05 lineages |
| `module_di25_rho05/data/di25_clustering_cM5_rho05.rds` | `module_di25_rho05/R/di25_ld_clustering_rho05.R` | module_di25_rho05 | list($groups: **20,807**×7 [group_id,Chr,representative,n_loci,score,has_eMLG,members]) | 2026-09-07 09:44 | `31829fc26214` | the 20,807 DI25-only, cM5, rho05 LD-reduced units — **"the DI25 high-DI LD units" (item 4), now this module's primary unit set**. Unit COUNT differs from the legacy 11,052 — `min_r2_rho=0.5` merges less aggressively than the old fixed 0.2, so this is a different partition, not a reindex | **authoritative, currently read** |
| `module_di25_rho05/data/di25_sorting_emlg_rho05.rds` | `module_di25_rho05/R/di25_sorting_rho05.R` | module_di25_rho05 | 20,807×24, identical column schema to the legacy `di25_sorting_emlg.rds`, incl. `unit_marker` (=rep_snp), `current_map_DI`, `n_aqu/n_pol/n_obs`, `sort_class` | 2026-09-07 10:03 | `6869ae2347d1` | per-unit sorting stats; `unit_marker` = **item 5, "representative-SNP identities,"** for this module's unit set | **authoritative, currently read** |
| `module_di25/data/di25_three_blocks.rds` | `module_di25/R_legacy/di25_pruning_test.R` | module_di25 | 3×9 (chr, anchor, n_units, n_markers, start_Mb, end_Mb, span_Mb, markers_per_Mb, group_ids) | 2026-08-08 12:22 | `1aa9d3575c20` | defines the 3 named polyctena regions (chromosome + Mb span). `pp_prep_units.R` step 5 resolves these spans into the CURRENT rho05 unit set by physical position (38 units) — the object's own `group_ids` column (legacy min_r2=0.2 unit IDs) is read but never used to select units directly | **legacy object, current resolution logic** — producing script lives in `R_legacy/`; the region definitions (chr + Mb span) are still used, the legacy group_ids are not |
| `moduleA_sorting/R/parallelism_stats.R` | (hand-written, not generated) | moduleA_sorting | code only: `parallelism_stats()`, `classify_sort()` | 2026-08-19 15:19 | `632a41349c8b` | the shared sort-classification engine this module calls directly (not a data object) | **authoritative, stable** — shared utility, not versioned by rho |
| `data/hybrids_and_parents_maf005.Rdata` | (repo-root shared build, not in either module) | shared upstream | `sample_data_with_parents` 195×5 [Population,Sample_ID,PC1,PC2,Mitotype]; `GTs_with_parents` 195×1,114,423; `map_hyb_005` 1,114,423×7 [Chr,Pos,marker,Polarity,DiagnosticIndex,maf_hyb,ld_w_095] | 2026-08-14 11:41 | `8a6c271ed3ce` | **item 8, population/individual metadata** (`sample_data_with_parents`), and the source of `current_map_DI` (`map_hyb_005$DiagnosticIndex`, the LATER full-data map's DI, explicitly NOT the DI25 ascertainment DI — see di25_sorting.R header) | **authoritative, shared** — used identically by module_manuscript_rho05's scripts |

## Table A-legacy — the superseded objects this module no longer reads

Kept only as a historical/comparison record (e.g. the AUDIT.md-era numbers in
git history were computed against these). **Do not use for any current or
future analysis in this module.**

| path | producing script | module | dims / key columns | mtime | sha256 (first 12) | status |
|---|---|---|---|---|---|---|
| `module_di25/data/di25_clustering_cM5.rds` | `module_di25/R/di25_ld_clustering.R` | module_di25 | list($groups: 11,052×7, same columns as the rho05 version) | 2026-08-07 16:36 | `51c365a5a3f1` | **superseded 2026-09-13** — fixed `min_r2=0.2`; replaced by `module_di25_rho05/data/di25_clustering_cM5_rho05.rds` above |
| `module_di25/data/di25_sorting_emlg.rds` | `module_di25/R/di25_sorting.R` | module_di25 | 11,052×24, same schema as the rho05 version | 2026-08-19 16:32 | `063a60059199` | **superseded 2026-09-13** — same lineage as above; replaced by `di25_sorting_emlg_rho05.rds` above |

The per-SNP sorting object (`module_di25/data/di25_sorting_snp.rds`) is
identical between the legacy and rho05 lineages by construction (it doesn't
depend on clustering, reused unchanged by `di25_sorting_rho05.R`); this
module does not read it directly either way (it computes its own per-unit
statistics from genotypes).

## This module's own stale output — resolved (audit item 1, 2026-09-13)

`data/pp_units_final.rds` (this module's own output, not an input from
`module_di25`) was found to still have 11,052 rows (the pre-rho05-migration
lineage) — a `saveRDS()` call in `pp_figures.R` that would have refreshed it
to 20,807 rows was dropped during an earlier revision and never re-added.
Grep-verified: no script in this module reads it. Moved, not deleted, to
`data/legacy/pp_units_final_minr2_02_11052.rds`. **`pp_residual_ancestry.rds$u`
(20,807 rows) is the canonical current unit table**; every script in this
module now asserts exact row count (20,807) and `Fmat`/`Resid` column
identity against it via `stopifnot()`. See `README.md` and
`FOLLOWUP_STATUS.md` for the same note.

## Table C — module_manuscript_rho05: corrected full-range DI–FST, neutral summary, empirical sorting (items 1–3)

These are **not currently read** by `module_population_partitioning` (this
module recomputes its own per-unit FST directly on the DI25 panel rather than
joining a full-range object). Traced here because a future synthesis step
will need them, and because they use a **different, incompatible MAF
convention** from Table A (see warning below).

> **Companion document, cross-checked 2026-09-13**:
> `module_manuscript_rho05/POPULATION_PARTITIONING_HANDOFF.md`, written by
> the peer session that owns `module_manuscript_rho05`, specifically to hand
> off these objects for this module's follow-up work. Re-verified its 4
> checksummed objects independently (size/mtime/MD5) — all match what this
> table already records; no staleness. It adds detail this table
> summarizes below rather than duplicates in full; read it directly before
> using any Table C/D object for a new analysis, especially §3 (three
> distinct, non-interchangeable neutral-FST quantities), §7 (known
> limitations) and §9 (stale objects under `stale_pre_fix_20260912/` — do
> not read from that directory).

| path | producing script | module | dims / key columns | mtime | sha256 (first 12) | role | status |
|---|---|---|---|---|---|---|---|
| `module_manuscript_rho05/data/di25_fst_vs_di_rho05.rds` | `module_manuscript_rho05/R/di25_fst_vs_di_rho05.R` | module_manuscript_rho05 | list($env: 10×10 [bin,DI_bin,**emp**,emp_ungated,n_emp_units,n_emp_units_ungated,n_sim_units,sim_med,sim_lo,sim_hi], $strat: 34×4 MAF-stratified, $neutral: length-3 overall summary, $sim_bg: 1000 background-LD replicates, $run_fingerprint) | 2026-09-13 14:16 | `78276ed9bab2` | **items 1 AND 2 combined**: corrected full-range (fixed DI bins, DI −Inf..+Inf) empirical FST-vs-DI (`env$emp`, MAF≥0.15-gated primary + `env$emp_ungated` sensitivity) AND its neutral-simulation comparison. `env$sim_med/lo/hi` (per-DI-bin, pooled-per-replicate) is **populated for DI bins 6–10 only** — the DI25 sim panel overlaps only ~2.5% of empirical units, concentrated in high-DI bins; bins 1–5 are NA. `$neutral` (length-3) is a SEPARATE, coarser quantity: one pooled genome-wide value per replicate (median 0.0213, 95% interval [0.0187,0.0248]), not usable for any per-bin or per-locus comparison — see the handoff's §3 for the full three-way distinction (pooled-per-replicate vs pooled-per-bin vs per-locus, in two different files, not interchangeable) | **authoritative** — post "AUDIT FIX Issue 2" (2026-09-12) + the 2026-09-13 parental-MAF-fold fix (see below), uses the complete 661,386-unit representation (Table D), not the old 17,509-row partial one |
| `module_manuscript_rho05/data/di25_fst_vs_di_violin_rho05.rds` | `module_manuscript_rho05/R/di25_fst_vs_di_violin_rho05.R` (confirmed via its `OUTRDS`/`saveRDS` call) | module_manuscript_rho05 | list($units: **395,996**×8 [group_id,best,rep_type,DI,bin,pmaf,fst_locus,DI_bin], $neutral, **$simO_locus: 16,615,717×1 per-locus unpooled neutral FST — the only neutral quantity usable for a locus-level comparison**, $n_sim_reps, $n_sim_overlap_units, $n_per_bin) | 2026-09-13 14:29 | `ec063ab3109c` | **per-unit** (not per-bin) FST for every MAF≥0.15-gated, full-genome LD-reduced unit — the full-range analogue of this module's per-unit `FST` column, but on a completely different unit set (395,996 full-genome MAF-gated units vs. this module's 20,807 DI25-only ungated units) | **authoritative** for full-range per-unit FST, post the 2026-09-13 parental-MAF-fold fix (previously 636,649 units under an unfolded-MAF bug — do not use that count or anything derived from it) |
| `module_manuscript_rho05/data/di25_fst_vs_di_sorting_stratified_rho05.rds` | `module_manuscript_rho05/R/di25_fst_vs_di_sorting_stratified_rho05.R` | module_manuscript_rho05 | list($units: 395,996×13, adds differentiated/prop_fixed/uni_score/sort_class/**sort_label**/directional to the violin units) | 2026-09-13 14:38 | `de18674aa998` | **item 3**: per-unit FST stratified by `sort_class`/`prop_fixed`, computed FRESH via `parallelism_stats()` on the same 395,996 MAF-gated units, **all 20 hybrid pops, Åland included** — matches this module's own population set (all 20, Åland included). NOT reused from `moduleA_stage1_cluster_sorting.rds` below, which excludes Åland and covers a narrower subset | **authoritative** for the full-range FST~sort_class question |
| `module_manuscript_rho05/data/moduleA_stage1_cluster_sorting.rds` (+ `_tau05/_tau06/_tau08` variants) | `module_manuscript_rho05/R/moduleA_stage1_cluster_sorting.R` | module_manuscript_rho05 | 18,361×9 [group_id,n_loci,differentiated,sort_class,DI,prop_fixed,uni_score,directional,sorted] | 2026-09-12 18:01 | `cb336adac6c6` | sort_class annotation, but **only for the 18,361-unit "Stage-1-direct" subset** (the BayPass/Omega climate-association unit universe, **Åland EXCLUDED, 19 hybrid pops**) — a DIFFERENT population set from both this module (20 pops) and `di25_fst_vs_di_sorting_stratified_rho05.rds` above (20 pops) | **narrower-scope AND different population set, do not treat as "the" full-genome sort_class annotation** — use `di25_fst_vs_di_sorting_stratified_rho05.rds$units$sort_class` for anything full-genome, and never mix its 19-pop calls into a 20-pop comparison |

**A third, older sort_class candidate to not confuse with either of the
above:** `moduleA_sorting/data/moduleA_cluster_sorting.rds` (produced by
`moduleA_sorting/R/moduleA_cluster_sorting.R`, sha256 `4dd2a71aca3d`, mtime
2026-08-19 13:57) is the **pre-rho05** full-genome sort_class annotation
(over `has_eMLG` full universe, Module A's original min_r2=0.2 lineage).
Neither this module nor, per `di25_fst_vs_di_sorting_stratified_rho05.R`'s
own header, `module_manuscript_rho05`'s stratified analysis reads it — it
predates the rho05 correction entirely. Listed only so it isn't mistaken for
either authoritative object above.

**MAF-gating warning (the concrete collision risk):** every object in this
table applies a **folded parental-MAF ≥ 0.15 gate as the PRIMARY filter**
(`MIN_PARENT_MAF <- 0.15`, matching Module A's locked convention pipeline-
wide — see `di25_fst_vs_di_rho05.R` "AUDIT FIX... DECISION" comment). The
DI25 objects this module actually reads (Table A) apply **no parental-MAF
gate at all** (`min_parent_maf = NULL` in `di25_sorting.R`, deliberately: "the
DI≥−25 selection is the ascertainment gate... the pooled-parental MAF≥0.15
gate Module A needs on the FULL set is likewise unnecessary here"). Both
choices are correct for their own scope, but **FST values are not
comparable between the two without accounting for this** — do not join or
overplot Table A's per-unit FST against Table C's without noting the MAF
convention differs.

**DI-vintage warning** (recorded verbatim from a cross-session note left in
`di25_fst_vs_di_sorting_stratified_rho05.R` by the peer session working on
`module_manuscript_rho05`, 2026-09-13): that script's DI comes from
`map_hyb_005$DiagnosticIndex` for the full genome, never pre-gated, and its
own DI>−25 cut is a fresh cut on that table — NOT a re-gating of the
already-DI>−25-restricted DI25 object. This module's DI25 objects (Table A)
must likewise never be re-gated on `map_hyb_005$DiagnosticIndex` (an unfixed
DIEM seed mismatch drops ~7k markers there — see `di25_sorting.R` header).
The two DI vintages can disagree at the margin; a naive side-by-side "high-DI"
comparison between the two modules could show a spurious DI-vintage artefact
rather than a real difference.

## Table D — full-genome rho05 unit / representative-SNP objects (item 5, full-genome scope)

Not the DI25-restricted units this module uses — the full-genome universe
Table C's per-unit FST objects are built from.

**A third parameter difference (beyond MAF-gating and fill-convention,
below): the merge distance cap and the eMLG-size threshold are also
different**, not just the min_r2 relative-threshold value the "rho05"
suffix refers to. Confirmed by reading `module0_ld_pruning_rho05/R/module0_ld_pruning_rho05_DIEM.R`
directly:

| parameter | this module's units (Table A, DI25-specific) | full-genome units (Table D) |
|---|---|---|
| `cM_threshold` (Stage-2 merge cap) | **5** cM | **0.5** cM — note the filename itself encodes this (`...cM05...`), easy to misread as "05 = rho05" rather than "0.5 cM"; it is the latter |
| `min_n_loci_flag` / `min_n_loci_eMLG` | 1 / **3** | 5 / **5** |
| `min_r2_rho` | 0.5 (both) | 0.5 (both) |

A unit's `n_loci` and `has_eMLG` status are therefore not comparable
across Table A and Table D even where marker IDs overlap.

| path | producing script | module | dims / key columns | mtime | sha256 (first 12) | role |
|---|---|---|---|---|---|---|
| `module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05.rds` | `module0_ld_pruning_rho05/R/module0_ld_pruning_rho05_DIEM.R` | module0_ld_pruning_rho05 | raw Stage-2 clustering object (full genome, `min_r2_rho=0.5`) | 2026-09-07 14:53 | `391207d96fde` | the full-genome rho05 clustering (input to the bestsnp object below) |
| `module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05_bestsnp.rds` | `module0_ld_pruning_rho05/R/module0_ld_pruning_rho05_DIEM.R` | module0_ld_pruning_rho05 | list($stats: **17,509**×13 [only eMLG/large-cluster units], $geno: 165×17,509, $rep_snp_all: **661,386**×6 [group_id,representative,n_loci,has_eMLG,best_marker,**rep_snp**]) | 2026-09-07 15:10 | `de0e37f70d97` | full-genome representative-SNP identities. **`$stats` alone is incomplete** (omits 643,877 singleton/small-cluster units) — use `$rep_snp_all` for the complete 661,386-unit representation, per the AUDIT FIX documented in `di25_fst_vs_di_rho05.R` |

**Representative-marker-convention warning (a second, separate collision
risk):** this full-genome object's best-SNP representative is built with
`eMLG_best_snp(..., fill = TRUE, round_fill = TRUE)` — **consensus-filled**
genotypes for missing calls. This module's DI25 units (Table A/A-legacy) use
`fill = FALSE` — **strictly observed** genotypes only, by deliberate choice
(`di25_sorting.R`: "we deliberately keep observed-only genotypes... Set
`fill = TRUE` for literal parity with Module A"). Do not treat a `rep_snp` /
`unit_marker` identity as carrying the same genotype-completeness semantics
across these two families of objects even where the marker ID matches.

## Table E — shared upstream inputs (items 8–9)

| path | dims / key columns | mtime | sha256 (first 12) | role |
|---|---|---|---|---|
| `data/hybrids_and_parents_maf005.Rdata` | see Table A | 2026-08-14 11:41 | `8a6c271ed3ce` | population/individual metadata + full-genome map; identical file used by both modules |
| `data/Frufa_DTOL_PR.ref_genome.recmap` | 27,648×4 [chr,pos,cM,cM/Mb] | 2026-06-04 10:52 | `bc0e6d419b9c` | the physical/genetic map; same file path used by `module_di25/R/di25_ld_clustering.R`, `module_di25_rho05/R/di25_ld_clustering_rho05.R`, and `module0_ld_pruning_rho05/R/module0_ld_pruning_rho05_DIEM.R` — single authoritative map, no versioning conflict found |

## Table F — this module's own generated objects (items 6–7, internal, not external inputs)

Not external inputs — listed because items 6–7 of the required audit name them
directly. Both are derived entirely from Table A objects; neither reads
anything from `module_manuscript_rho05` or the rho05 lineage.

| path | built by | dims / key columns | mtime | sha256 (first 12) | upstream provenance |
|---|---|---|---|---|---|
| `module_population_partitioning/data/pp_units_Fmat.rds` | `R/pp_prep_units.R` | list($u: 20,807-row unit table, $Fmat: **20×20,807** oriented-aquilonia-frequency matrix — **item 6**, $hybrid_pops, $blk_rho05: the 3 named blocks resolved by physical position) | 2026-09-13 15:21 | `bb4022721f8c` | built from Table A's `di25_inputs.rds` + `di25_clustering_cM5_rho05.rds` + `di25_sorting_emlg_rho05.rds` + `di25_three_blocks.rds` exclusively |
| `module_population_partitioning/data/pp_residual_ancestry.rds` | `R/pp_residualize_ancestry.R` | list($H_loco: 20×27 pop×chr LOCO ancestry, $Resid: 20×20,807 — **item 7**, $alpha_u/$beta_u/$R2_u, $u) | 2026-09-13 15:24 | `4fee798f93a5` | built from `pp_units_Fmat.rds` only |

`data/` is not committed to git for this module (matches repo convention —
see README "Outputs"); both are regenerable by re-running `R/pp_prep_units.R`
then `R/pp_residualize_ancestry.R` against Table A's inputs, unchanged.

## Confirmed non-collision (cross-session, 2026-09-13)

Verified directly with the peer session (`ldscnr-12`) working on
`module_manuscript_rho05`: its units come from
`module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05_bestsnp.rds` (Table
D, full-genome), never touching `module_di25`'s (or `module_di25_rho05`'s)
DI25-restricted clustering anywhere. Verified at the time against the
then-current 11,052-unit legacy DI25 clustering; the subsequent migration to
the 20,807-unit rho05 DI25 clustering (Table A) does not change this
conclusion — `module_manuscript_rho05` reads from `module0_ld_pruning_rho05`
(full-genome), not from either `module_di25` or `module_di25_rho05`
(DI25-restricted), regardless of which DI25-restricted lineage this module
itself uses. Same locked conventions (Weir & Cockerham FST, `classify_sort()`
tau=0.6/binom) but fully independent unit sets and different questions (this
module: multilocus partition concordance between units; `module_manuscript_rho05`:
per-locus FST~DI dose-response stratified by sort_class/prop_fixed). No
duplication.

## Consistency-check status for follow-up work

- Marker/unit identifiers: `unit_marker` (Table A) and `rep_snp`/`best_marker`
  (Tables C/D) are both `Chr<n>:<pos>` SNP IDs and so COULD be joined
  directly on marker identity across the DI25 (20,807) and full-genome
  (395,996 / 661,386) unit sets — but a join would need to resolve the
  MAF-gating and fill-convention differences above first; not attempted here
  (no follow-up analysis requested yet).
- rho verification: Table A/A-legacy confirmed to differ ONLY in Stage-2 `min_r2`
  vs `min_r2_rho`, by reading both scripts directly (not inferred from
  filenames).
- Full-range vs frozen high-DI panel: kept fully separate in this document
  (Table A = frozen DI25 panel; Table C/D = full-range). No script in this
  module reads a full-range object.
- Per-SNP vs per-unit vs replicate-level FST: this module computes per-unit
  FST only (Weir & Cockerham, hybrid pops, on the DI25 best-SNP/representative
  genotype). Table C's `di25_fst_vs_di_violin_rho05.rds$units$fst_locus` is
  also per-unit (on its own, different, unit set); `$sim_bg`/replicate-level
  values in `di25_fst_vs_di_rho05.rds` are a separate, bin-level neutral-sim
  quantity, not per-unit.
- Best-observed-SNP vs eMLG: this module (Table A) and its superseded predecessor (Table A-legacy) both use strictly
  observed (`fill=FALSE`) best-SNP/representative genotypes throughout; Table
  D's full-genome object uses consensus-filled (`fill=TRUE`) genotypes — kept
  distinct above, not interchanged anywhere in this module's current code.
- No legacy uncorrected MAF-gating result has entered any object this module
  currently reads (Table A predates the MAF-gate question entirely — it was
  never gated, by design, not "uncorrected").
