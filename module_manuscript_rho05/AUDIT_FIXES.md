# Audit and repair log -- module_manuscript_rho05

Started 2026-09-12. Stale pre-fix outputs are quarantined under
`stale_pre_fix_20260912/` (mirroring the original directory structure)
rather than deleted. Checksums and environment details are in
`AUDIT_MANIFEST_20260912.md`.

Decisions made by the pipeline author before execution began:
- **Parental-MAF gate (Issue 2):** apply MIN_PARENT_MAF=0.15 as the PRIMARY
  gate for the Fst-vs-DI analysis (matching Module A's locked convention
  pipeline-wide); the ungated/all-strata curve is reported as a sensitivity
  check, not the headline result.
- **bio6/bio11 (Issue 7):** reduce to ONE winter-temperature axis (their
  mean, after standardization) before null calibration, rather than
  calibrating both as nominally-independent tests.

---

## Issue 1: BayPass covariate scaling (PC1/PC2/bio6/bio11 vs structured null)

**Status: FIXED, observed Stage-1-unit scans rerun, floor-survivor/Module C
recalibrated. Full-SNP reruns in progress on mini2 (see "Unfinished work").**

**Verification (before any change):** `u.PC1` n=19 mean=0.190 **SD=3.267**;
`u.PC2` SD=2.236; `u.bio6`/`u.bio11` already close to standardized (SD~0.96)
but not exactly, and not computed over precisely this 19-population set.
`-nocovscaling` used in every observed-covariate BayPass call across all
three driver scripts. Null covariates (`moduleB_stage1_S1units_null.R`,
`moduleC_stage1_null_regen.R`) are drawn via `scale(...)` (SD 1) with
`-nocovscaling` also set -- correctly matched to the unit-SD default
beta-prior grid. Observed and null BF were computed under a mismatched
effective covariate scale.

**Files changed:**
- NEW `module_manuscript_rho05/R/moduleB_stage1_standardize_covariates.R` --
  rho05-specific fix (does NOT touch the shared legacy
  `moduleB_climate_GEA/R/moduleB_write_baypass_inputs.R`, which also writes
  unstandardized PC1/PC2 and is used by other, unaudited modules -- flagged
  below, not fixed here). Standardizes PC1/PC2/bio6/bio11 to mean 0 / SD 1
  over the exact 19-population BayPass order in both `aland_excluded/` and
  `aland_excluded_S1units/`; writes provenance to
  `data/moduleB_stage1_covariate_standardization.csv`; asserts mean~0
  (<1e-8), SD~1 (<1e-8), length 19, per covariate.
- NEW `moduleB_stage1_prepare_bio_winter_covariate.R` -- builds the combined
  `u.bio_winter` (mean of standardized bio6/bio11, re-standardized itself).
- NEW `rerun_observed_standardized_S1units.sh` / `rerun_observed_standardized_fullsnp.sh`
  -- rerun the 4 observed covariates at each resolution with corrected
  covariates; same BayPass parameters as the originals (nthreads 10,
  nocovscaling, nval 500, burnin 5000, thin 25, seed 74); pre-fix outputs
  quarantined via `cp` before being overwritten.
- Omega estimation step in `run_baypass_stage1.sh`: added `-seed 74` for
  future reproducibility (Issue 9) -- **Omega itself was NOT re-estimated**;
  frozen, checksum `7abb549341c5c08a7ba238892c7c0740` (unchanged throughout).

**Reruns performed:**
- 4 Stage-1-unit observed BayPass scans (PC1/PC2/bio6/bio11): ~2.8 min each,
  local machine. DONE.
- 1 Stage-1-unit `bio_winter` scan: DONE.
- 4 full-SNP observed BayPass scans: launched on mini2 (~1h each observed so
  far, in progress -- see "Unfinished work").
- Floor-survivor null (`moduleB_stage1_S1units_null.rds`): k1/k2 and
  floor1/floor2 recomputed by re-reducing the ALREADY-PERSISTED
  `null/bf_matrices/cRegen_bf_b##.rds` (50 files, from `moduleC_stage1_null_regen.R`'s
  earlier run) against the corrected observed BF -- **no BayPass rerun**,
  because the null side was already correctly standardized.
- Module C (`moduleC_stage1_null_regen.R` + `moduleC_stage1_analyse.R`):
  re-reduced from the same persisted matrices with the corrected observed
  BF AND the corrected annotations (Issue 5) -- Monte-Carlo equivalence gate
  r=1.00000 exactly (both sides reduce the identical persisted matrices, not
  independent MCMC runs).
- `bio_winter` floor-survivor calibration: new, reduced from the same
  persisted matrices (no BayPass).

**Old vs corrected headline results (PC1/PC2 floor-survivor test, Stage-1-direct, n=18,361):**

| | PC1 floor survivors | PC1 FDR | PC2 floor survivors | PC2 FDR |
|---|---|---|---|---|
| Pre-fix (unstandardized covariates) | 3 | ~0.61 | 3 | ~0.61 |
| **Corrected** | **1** | **1.84** (below null expectation) | **2** | **0.92** |

The pre-fix "3/3 survivors, FDR~0.61" result -- previously read as confirming
a Stage-1-direct signal recovery over the Stage-2 scan -- was substantially a
scaling artifact. **Neither PC1 nor PC2 is a defensible discovery set after
the fix.** bio_winter, by contrast, is notably strong: 10 of 18,361 floor
survivors, FDR~=0.18 (never calibrated before this audit -- see Issue 7).

**Not rerun and why:** the canonical (Stage-2, `moduleB_climate_GEA`) eMLG
climate scan is a SEPARATE, unaudited pipeline that also uses the same
shared `write_baypass_inputs()` function and therefore likely has the SAME
unstandardized-covariate issue -- **out of scope for this rho05 audit**,
flagged here for awareness only. Fixing it would require its own decision
and rerun cycle.

---

## Issue 2: Fst-vs-DI unit universe (di25_fst_vs_di_rho05.R)

**Status: FIXED and fully rerun (1000/1000 replicates).**

**Verification:** `b$stats$best_marker` (previously used) has 17,509 rows
(large clusters with an eMLG). The complete LD-reduced representation,
`b$rep_snp_all`, has 661,386 rows: `rep_snp` = `best_marker` where
`has_eMLG==TRUE` (17,509), else the centrality `representative` for the
643,877 small/singleton clusters. Verified 661,386 unique `group_id`, 661,386
unique `rep_snp`, no NAs, no dups.

**Files changed:** `module_manuscript_rho05/R/di25_fst_vs_di_rho05.R`,
substantially rewritten:
- Uses `b$rep_snp_all` (all 661,386 units), with `rep_type` (best-SNP vs
  representative) tracked and reported separately.
- Parental MAF>=0.15 applied as the PRIMARY gate (decision above); the
  ungated curve is now an explicit sensitivity line on the same figure, plus
  the pre-existing MAF-stratified breakdown.
- DI and parental MAF joined explicitly by marker name (`match()`), not
  positionally.
- Fresh cache directory `data/fst_sim_cache_full_rho05_v2/` (old
  `fst_sim_cache_full_rho05/` quarantined, NOT reused -- different best/rep-SNP
  marker set invalidates it). Every cache entry now carries a fingerprint
  (marker-list hash, DI-bin hash, parameter hash, source `.bed` file path +
  md5, script-version string); aggregation ABORTS if fewer than the
  requested 1000 replicates are present or any fingerprint mismatches.

**Rerun:** piloted at 5 reps (validated: sim-panel overlap exactly 16,616 as
predicted, vs the old 1,511), then the full 1000-rep run completed
successfully (0 failures, all fingerprint-verified).

**Old vs corrected headline result:**

| | Pre-fix (17,509 units) | Corrected (661,386 units, MAF>=0.15 primary) |
|---|---|---|
| Sim-panel overlap | 1,511 | **16,616** (exact match to the audit's prediction) |
| Empirical Fst, neutral bg | ~0.033 | 0.021 |
| Empirical Fst, diagnostic | ~0.279 | 0.255--0.257 |
| Neutral sim Fst | 0.021 [0.018,0.025] | 0.021 [0.019,0.025] |

Pattern preserved (strong elevation of Fst in diagnostic DI bins over the
neutral simulation baseline); the qualitative conclusion is unchanged, now
computed over the complete, correctly-gated unit universe.

---

## Issue 3: 17-lineage-unit Omega ordering (moduleB_ancestry_climate_mitotype_CORRECTED.R)

**Status: FIXED and rerun; figures regenerated; manuscript doc updated.**

**Verification:** `Omega17` was built via `M %*% Omega19 %*% t(M)` with `M`'s
row order = `units` (`unique()`'s first-appearance order), while `dt17` is
built via `setorder(dt17, unit_id)` (alphabetical). Confirmed by direct
comparison these are genuinely different orderings (`identical(units,
sort(units))` is FALSE). Every Omega-null draw (`NullMat17`, from `Omega17`'s
own eigendecomposition) was compared directly against `dt17`'s columns with
no name-based alignment.

**Fix:** `dimnames(Omega17) <- list(units, units)`, then
`Omega17 <- Omega17[match(dt17$unit_id, rownames(Omega17)), ...]`, with a
hard assertion that rownames/colnames now equal `dt17$unit_id` exactly.

**Independent validation:** all 6 corrected Omega-null p-values match an
independently-computed cross-check EXACTLY: ancestry-PC1 0.0744,
ancestry-PC2 0.1898, ancestry-bio6 0.1416, ancestry-bio11 0.2915,
ancestry-mitotype 0.0501, partial PC1-ancestry|mitotype 0.1411.

**Old vs corrected (Omega-null p only; raw correlations and block-permutation
unaffected by construction, verified unchanged):**

| Test | Pre-fix Omega-null p | Corrected |
|---|---|---|
| ancestry-PC1 | 0.044 | 0.074 |
| ancestry-PC2 | 0.110 | 0.190 |
| ancestry-bio6 | 0.091 | 0.142 |
| ancestry-bio11 | 0.210 | 0.292 |
| ancestry-mitotype | 0.037 | 0.050 |
| partial PC1-ancestry\|mitotype | 0.097 | 0.141 |

Every case moved to a WEAKER (larger) p-value -- the pre-fix misalignment
had made the Omega-null look more significant than it actually is. No test
now clears p<0.05 under the Omega-null (ancestry-mitotype at 0.050 is the
closest).

**Rerun:** `moduleB_ancestry_climate_mitotype_CORRECTED.R` and
`_figs.R`, both complete; 3 figures regenerated.

**Superseded scripts marked (task item 7):**
`moduleB_ancestry_vs_winter_climate.R`, `moduleB_mitotype_vs_ancestry_climate.R`
(own printed results/figures superseded, but their `.rds` data-prep outputs
remain valid, actively-read inputs -- NOT deleted), and
`moduleB_PC1_ancestry_partial_mitotype.R` (fully superseded, has the
Omega-null variable-substitution bug from the original Erratum; nothing
downstream reads its output). All three now carry a clear header comment;
their standalone figures/`.rds` output quarantined.

**Manuscript doc:** `ancestry_climate_mitotype.tex` updated throughout
(Sections 1/3/4/5/6, Summary) with all corrected numbers; a "Second Erratum"
note added documenting this fix. Recompiled cleanly (12 pages).

---

## Issue 4: BayPass test-universe framing (documentation only, no code/data change)

**Status: DONE.** This is a scientific SCOPE clarification, not a bug --
recorded as the decision: **the ≥5-SNP universe (18,361 of 698,251 Stage-1
clusters) is retained**; isolated SNPs and clusters with <5 markers were
never tested by the Stage-1-direct climate scan. Added an explicit "Scope
note" to `ancestry_climate_mitotype.tex`'s Motivation section and to
`moduleC_stage1_analyse.R`'s generated report (Data provenance bullet),
both stating the 18,361-of-698,251 framing plainly and noting the complete
661,386-unit set IS used elsewhere (Fst-vs-DI, Issue 2) but not for this
climate calibration. No rerun needed (nothing computational changed).

---

## Issue 5: Aland harmonization (moduleA_stage1_cluster_sorting.R)

**Status: FIXED and rerun; Module C annotations/null-regen re-reduced.**

**Verification:** `hybrid_pops` was derived from `sample_data$Population`
(20 hybrid populations, including Aland), while BayPass/Omega exclude Aland
(19). `best$geno`'s rows (165 hybrid individuals) included 10 Aland
individuals.

**Fix:** `hybrid_pops <- setdiff(hybrid_pops, "Aland")`; Aland individuals
(10 of 165) dropped from the hybrid genotype rows BEFORE
`ohta_fast_prepare()`/`parallelism_stats()`, not just from the population
list argument. Added `attr(..., "meta")` recording the exact population list
to every saved sorting object. The separate 20-population DI25 sorting
analyses (`module_di25/`) are untouched.

**Old vs corrected directional-classification counts (differentiated count,
10,170, is UNCHANGED -- a property of parent MAF/DI, not hybrid composition):**

| tau | Pre-fix (20 pops) directional | Corrected (19 pops, Aland excluded) |
|---|---|---|
| 0.5 | 3,199 | 2,961 |
| 0.6 (primary) | 2,197 | 1,945 |
| 0.8 | 604 | 384 |

**Rerun:** `moduleA_stage1_cluster_sorting.R`,
`moduleC_stage1_annotations.R` (both instant), and
`moduleC_stage1_null_regen.R` re-reduced from the persisted BF matrices with
the corrected annotations (~5 min, no BayPass).

---

## Issue 6: Null-calibrated terminology in outlier figures

**Status: FIXED for the region-level script (rerun/regenerated); SNP-level
script code-fixed, regeneration pending the mini2 full-SNP rerun.**

**Files changed:**
- `moduleB_stage1_region_manhattan.R`: removed the obsolete "structured null
  doesn't exist yet" comment (it does now). Region assembly and primary
  highlighting now use the FLOOR-SURVIVOR set (raw crossing AND beats all
  10,000 null draws), not raw threshold crossings. Raw-crossing-only units
  get a subdued secondary style (pink/magenta -- changed from an initial
  amber choice that was too close to viridis's yellow end, per user
  feedback) rather than being merged into a region. Where zero units survive
  (mitoC2: 0 of 11), this is stated explicitly and no region is built.
- `moduleB_stage1_snp_manhattan_ldmanhattan.R`: same floor-survivor-primary
  logic added for PC1/PC2/mitoC2 (each has its own null object); bio6/bio11
  (no per-variable null -- see Issue 7) remain raw-threshold-only, now
  explicitly labelled as uncalibrated in the subtitle rather than implied
  validated.
- Quarantined the obsolete `*_PREVIEW.*` bio6/bio11 figures and their
  generating script (`moduleB_stage1_bioclim_snp_manhattan_preview.R`).

**Rerun:** region-Manhattan (PC1/PC2/mitoC2) regenerated --
PC1: 45 raw crossings -> 1 floor survivor -> 1 region (Chr24).
PC2: 65 raw crossings -> 2 floor survivors -> 2 regions (Chr3, Chr4).
mitoC2: 11 raw crossings -> 0 floor survivors -> 0 regions (stated explicitly).
SNP-level Manhattan (5 panels): code fixed, NOT yet regenerated -- needs the
full-SNP scan rerun (Issue 1) to complete first (see "Unfinished work").

---

## Issue 7: bio6/bio11 combined winter-temperature calibration

**Status: FIXED (decision applied) and fully calibrated.**

bio6/bio11 correlation: r=0.960 (raw standardized covariates), r=0.94 (BF
outputs) -- essentially one signal. Per the decision, reduced to
`u.bio_winter` = re-standardized mean of standardized bio6/bio11 (built by
NEW `moduleB_stage1_prepare_bio_winter_covariate.R`). Ran its own Stage-1-unit
BayPass scan and calibrated against the persisted null matrices (no new
BayPass on the null side).

**Result: 10 of 18,361 Stage-1 units survive the floor test, set FDR~=0.18**
-- notably lower (more significant) than PC1's 1.84 or PC2's 0.92. This is
the single strongest floor-survivor result of any covariate in this pipeline
and was never calibrated before this audit (the raw ~124/108 bio6/bio11
threshold crossings previously had no null test at all). Saved
`data/moduleB_stage1_bio_winter_null.rds`. Folded into
`ancestry_climate_mitotype.tex` (Section 5 table + Bottom line + Summary).

**Not yet done:** a full-SNP `bio_winter` scan (for a per-SNP Manhattan
panel) was not run -- out of scope for this audit unless requested; the
existing separate bio6/bio11 full-SNP scans (queued on mini2) remain
useful for visualization, captioned as two correlated components of this
one calibrated axis.

---

## Issue 8: Best-SNP/representative consistency (di25 scripts)

**Status: FIXED and rerun.**

- `di25_recomb_tau_sweep_rho05.R`: was assigning recombination via the
  clustering's centrality `representative`; now uses the saved `unit_marker`
  from `di25_sorting_emlg_rho05.rds` (the marker actually analysed for
  sorting) -- verified these differ for 689 of 20,807 units (3.3%). Rerun;
  qualitatively similar coefficients (all still strongly negative across the
  tau grid); Fig 3's panel (b) regenerated via `di25_sorting_multipanel_rho05.R`.
- `di25_population_fixation_emlg_rho05.R`: was using `consensus_dosage()`
  (averaged across cluster members) for clusters with >2 markers, a
  DIFFERENT genotype representation from the sorting analysis itself.
  Replaced with the same `unit_marker`'s raw observed genotype (fill=FALSE
  -- no consensus-filling of missing calls). Relabelled "LD-reduced
  best-SNP/representative-SNP units." Rerun; new result 21.4% aqu-fixed /
  12.5% pol-fixed cells.
- Stage-1 BayPass scan's `fill=TRUE, round_fill=TRUE` (in
  `eMLG_best_snp()`/`moduleB_stage1_prepare_baypass_inputs.R`) is
  RETAINED as-is (not part of this fix's scope) -- already documented in
  that script's own header that genotypes are the selected best-SNP with
  only missing calls filled from the oriented cluster consensus; `n_filled`
  is recorded per unit in `best$stats`.

---

## Issue 9: Reproducibility hardening

**Status: PARTIAL.** High-value, low-risk fixes applied; broader retrofits
across every historical script/RDS were judged disproportionate (see below).

Done:
- `-seed 74` added to the (frozen, not rerun) Omega-estimation BayPass call.
- `AUDIT_MANIFEST_20260912.md`: checksums for the frozen Omega/genotype/
  poolsize files, the corrected covariate files, BayPass executable versions
  (local v3.0, mini2 v3.1), R/data.table versions, and the Formica repo
  commit hash.
- All NEW scripts this audit wrote (`rerun_observed_standardized_*.sh`,
  `moduleB_stage1_standardize_covariates.R`,
  `moduleB_stage1_prepare_bio_winter_covariate.R`) use `set -euo pipefail`
  (shell) / explicit `stopifnot()` assertions (R), explicit exit-status and
  required-output checks, and atomic writes where they persist anything
  reusable.
- `moduleC_stage1_null_regen.R` and the rewritten `di25_fst_vs_di_rho05.R`
  already had (or now have) proper fingerprint-based cache validation
  (reject existence-only caches), safe `done==NBATCH` handling, and atomic
  writes for persisted BF matrices / cache entries.

NOT done (flagged, not retrofitted, to avoid destabilizing already-verified
outputs under time pressure):
- The ORIGINAL `run_baypass_stage1.sh` / `_fullsnp.sh` / `_fullsnp_bioclim.sh`
  driver scripts lack per-call exit-status/output-existence checks (they do
  have `set -euo pipefail` already). Adding checks to already-run, working
  scripts was judged low-value/nonzero-risk; noted here as a gap for any
  FUTURE rerun of those exact scripts.
- `moduleB_stage1_mitoC2_null.R` and `moduleB_stage1_S1units_null.R` do not
  carry input fingerprints in their checkpoints (unlike
  `moduleC_stage1_null_regen.R`) -- both already completed successfully;
  retrofitting would require a full rerun for no correctness benefit.
- Final RDS objects (e.g. `moduleC_stage1_null_stats.rds`,
  `di25_fst_vs_di_rho05.rds`) are written via plain `saveRDS()`, not
  atomic tmp+rename, in most scripts -- a lower-priority robustness gap
  across dozens of scripts, not retrofitted here.
- Formica/LDscnR commit hashes and BayPass version are NOT embedded inside
  every individual result manifest (only in `AUDIT_MANIFEST_20260912.md`
  centrally) -- embedding them per-object would touch every script in the
  pipeline; out of scope for this pass.

**External archive location for the 50 null `.env` files and persisted BF
matrices:** these live on `mini2` (`~/formica_hybrid/baypass_stage1_S1units/null/`),
NOT in this local git-tracked module (BayPass working directories are
`.gitignore`d, ~1.5-5GB). Retrieval: `scp mini2:~/formica_hybrid/baypass_stage1_S1units/null/bf_matrices/*.rds <dest>`.
Same location holds the mitoC2 null's persisted matrices and logs.

---

## Verification checklist (Issue 10)

- [x] All modified/new R scripts parse (`parse()` check run on every one
      before execution).
- [x] BayPass runs finish successfully: all reruns checked for exit status
      and required non-empty output files.
- [x] Genotype/map/BayPass row orders match: asserted in every rerun script
      (`identical(group_order, ...)`, `nrow(s) == nrow(cl5)`, etc.).
- [x] Covariates and Omega match population order BY NAME: the
      standardization script derives `pop_order` from `sample_data` and
      asserts length==19 against `u_DIEM.size`; the Issue-3 fix asserts
      `identical(rownames(Omega17), dt17$unit_id)` explicitly (not just
      dimension).
- [x] Observed and null covariates use identical scaling: verified (Issue 1)
      -- both mean 0, SD 1, `-nocovscaling` on both sides.
- [x] Every null contains exactly 10,000 draws: asserted in
      `moduleC_stage1_null_regen.R` (`nrow(null_list[[k]]) == NSIM_TOTAL`)
      and in the floor-survivor recompute scripts (50 batches x 200).
- [x] Every simulation aggregate contains exactly 1,000 validated replicates:
      `di25_fst_vs_di_rho05.R` aborts otherwise (verified: ran to completion,
      1000/1000, 0 failures, all fingerprints matched).
- [ ] **Figures use corrected data and are newer than dependencies:
      PARTIAL.** All figures EXCEPT the 5-panel SNP-level Manhattan plots
      (`moduleB_stage1_snp_manhattan_ldmanhattan.R`) have been regenerated
      and verified newer than their inputs. The SNP-level Manhattan plots
      are code-corrected but NOT yet regenerated (blocked on the mini2
      full-SNP BayPass rerun, in progress).
- [x] Manuscript/report text contains no stale candidate counts or invalid
      Omega p-values: `ancestry_climate_mitotype.tex` fully swept and
      updated (verified via grep for every old numeric value); recompiled
      cleanly.

## Unfinished work / how to resume

**Full-SNP BayPass reruns (PC1 DONE; PC2 in progress; bio6, bio11 queued;
~1h each)** running on `mini2` in `~/formica_hybrid/baypass_stage1_fullsnp/`,
launched via `nohup ./run_fullsnp_mini2.sh > run_fullsnp.log 2>&1 &`. To
check status: `ssh mini2 "tail -20 ~/formica_hybrid/baypass_stage1_fullsnp/run_fullsnp.log"`.
PC1's corrected output has already been pulled back and used (see below);
copy PC2/bio6/bio11's `*_fullSNP_stage1Omega_summary_*.out` files back to
`module_manuscript_rho05/baypass_stage1/aland_excluded/` as each finishes
(quarantining the stale pre-fix copies first), then rerun:
```
Rscript module_manuscript_rho05/R/moduleB_stage1_snp_manhattan_ldmanhattan.R
```
to regenerate all 5 SNP-level Manhattan panels with corrected data (code
already fixed for Issue 6's floor-survivor terminology).

**Partial reporting already done (2026-09-12):** PC1 (using the freshly
pulled corrected full-SNP scan) and mitoC2 (never affected by the covariate-
scaling bug, so its existing full-SNP data was already valid) were
regenerated via a temporary scratch copy of the script with the PC2/bio6/bio11
calls commented out -- NOT committed to the repo (the committed script still
processes all 5 panels; rerun it as-is once PC2/bio6/bio11 are ready). PC1:
1 of 45 raw crossings survives the floor test (matches the region-Manhattan
result). mitoC2: 0 of 11 survives (also matches).

**Out-of-scope / flagged, not fixed in this pass:**
- The canonical (Stage-2) `moduleB_climate_GEA` eMLG climate scan likely has
  the same Issue-1 covariate-scaling bug (shares `write_baypass_inputs()`) --
  not touched; would need its own decision + rerun cycle.
- Issue 9's broader retrofits (per-historical-script exit-status checks,
  atomic final-RDS writes, per-object embedded commit hashes) -- see Issue 9
  section above for the specific list and rationale for deferring.
