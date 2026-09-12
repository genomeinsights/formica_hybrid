# Audit and repair log -- module_manuscript_rho05

Started 2026-09-12. This log is built incrementally as each issue is verified,
fixed, and (where practical) rerun. Stale pre-fix outputs are quarantined
under `stale_pre_fix_20260912/` rather than deleted.

Decisions made by the pipeline author before execution began:
- **Parental-MAF gate (Issue 2):** apply MIN_PARENT_MAF=0.15 as the PRIMARY
  gate for the Fst-vs-DI analysis (matching Module A's locked convention
  pipeline-wide); the current ungated/all-strata curve becomes a reported
  sensitivity check, not the headline result.
- **bio6/bio11 (Issue 7):** reduce to ONE winter-temperature axis (derived
  from both, not a post-hoc pick of whichever looks stronger) before null
  calibration, rather than calibrating both as nominally-independent tests.

---

## Issue 1: BayPass covariate scaling (PC1/PC2/bio6/bio11 vs structured null)

**Status: VERIFIED, fix in progress.**

**Verification (before any change):**
- `module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/u.PC1`: n=19, mean=0.190, **SD=3.267**
- `.../u.PC2`: n=19, mean=0.283, **SD=2.236**
- Both `run_baypass_stage1.sh`, `run_baypass_stage1_fullsnp.sh`, and
  `run_baypass_stage1_fullsnp_bioclim.sh` pass `-nocovscaling` to every
  observed-covariate BayPass call (6 call sites checked).
- The structured-null covariates (`moduleB_stage1_S1units_null.R`'s
  `draw_null()`, and `moduleC_stage1_null_regen.R`'s null draws, which reuse
  the same `.env` files) are generated via `scale(...)`, i.e. mean 0, SD 1,
  and BayPass is ALSO run with `-nocovscaling` on those files.
- Net effect: BayPass's covariate-mode beta-prior grid (`-minbeta -0.3
  -maxbeta 0.3` by default) is calibrated for a unit-SD covariate; observed
  PC1/PC2 (raw SD 3.27/2.24) violate that scale while every null draw
  satisfies it. Observed and null Bayes factors are NOT computed under a
  comparable effective model -- confirmed real, not a false claim.
- This affects every downstream object built from the observed PC1/PC2/bio6/
  bio11 BF vectors this session: `moduleB_stage1_S1units_null.rds`,
  `moduleC_stage1_*`, both Manhattan-plot scripts, and Section 6 of
  `ancestry_climate_mitotype.tex` (folded in immediately prior to this audit).

**Fix plan:**
1. Standardize u.PC1/u.PC2/u.bio6/u.bio11 (both `aland_excluded` [full-SNP]
   and `aland_excluded_S1units` dirs) to mean 0, SD 1 over the exact 19-pop
   BayPass order; save a provenance table; assert mean~0/SD~1/length 19/
   correct order.
2. Freeze Omega (checksum recorded below) -- not re-estimated.
3. Rerun the 8 observed BayPass scans (PC1/PC2/bio6/bio11 x {Stage-1-unit,
   full-SNP}) with the corrected covariate files.
4. Recalibrate against the EXISTING null BF matrices -- NOT redrawn/rerun,
   because the null side was already correctly standardized:
   - Floor-survivor k1/k2 (`moduleB_stage1_S1units_null.rds`'s role): reduced
     from `moduleC_stage1_null_regen.R`'s PERSISTED per-batch BF matrices
     (`null/bf_matrices/cRegen_bf_b##.rds`, 50 files, still on mini2), which
     were run on the identical `null_b01..50.env` files with the identical
     Omega/geno/poolsize -- valid for PC1/PC2/bio6/bio11 alike, since the null
     BF does not depend on which observed axis it is later compared against.
   - Module C's null-statistics object (`moduleC_stage1_null_stats.rds`) needs
     only its `observed` rows recomputed (`compute_covariate_stats()` on the
     corrected PC1/PC2 BF) -- the null side is untouched.
   - bio6/bio11 floor-survivor calibration (Issue 7, previously never done)
     reuses the same persisted matrices.

(Filled in below as each step completes.)
