#!/usr/bin/env bash
# AUDIT FIX (Issue 1): rerun the 4 observed Stage-1-unit BayPass scans
# (PC1, PC2, bio6, bio11) now that their covariate files have been
# standardized (moduleB_stage1_standardize_covariates.R) to match the
# Omega-structured null convention (mean 0, SD 1, -nocovscaling on both
# sides). Same BayPass parameters as the original runs (nthreads 10,
# nocovscaling, nval 500, burnin 5000, thin 25, seed 74) -- only the
# covariate VALUES changed, not the model. Omega is untouched (frozen).
#
# Pre-fix outputs for these 4 covariates are quarantined into
# stale_pre_fix_20260912/ before being overwritten (copy, not destructive).
#
# Run from the repo root: bash module_manuscript_rho05/baypass_stage1/rerun_observed_standardized_S1units.sh

set -euo pipefail
cd "$(dirname "$0")"

PATH_TO_BAYPASS=/Users/petrikem/gitlab/baypass_public-master/sources/g_baypass
CORES=10
UNIT_DIR=aland_excluded_S1units
STALE=../stale_pre_fix_20260912/baypass_stage1/${UNIT_DIR}
mkdir -p "$STALE"

for v in PC1 PC2 bio6 bio11; do
  echo "=== quarantining pre-fix ${v} Stage-1-unit outputs ==="
  cp -f ${UNIT_DIR}/${v}_S1units_withOmega_summary_*.out ${UNIT_DIR}/${v}_S1units_withOmega_baypass.log "$STALE"/ 2>/dev/null || true

  echo "=== ${v} association, Stage-1 units, WITH Omega (standardized covariate) ==="
  "${PATH_TO_BAYPASS}" \
    -countdatafile "${UNIT_DIR}/u_S1units.geno" \
    -omegafile     "${UNIT_DIR}/omega_mat_omega.out" \
    -efile         "${UNIT_DIR}/u.${v}" \
    -poolsizefile  "${UNIT_DIR}/u_DIEM.size" \
    -nthreads "${CORES}" \
    -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 \
    -outprefix "${UNIT_DIR}/${v}_S1units_withOmega"
  [ -s "${UNIT_DIR}/${v}_S1units_withOmega_summary_betai_reg.out" ] || { echo "FAILED: ${v} S1units output missing/empty"; exit 1; }
done

echo "=== Done. Rerun outputs (standardized covariates) in ${UNIT_DIR}/: {PC1,PC2,bio6,bio11}_S1units_withOmega_* ==="
