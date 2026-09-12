#!/usr/bin/env bash
# AUDIT FIX (Issue 1): rerun the 4 observed full-SNP BayPass scans
# (PC1, PC2, bio6, bio11 -- 1,114,423 SNPs each) with standardized
# covariates. Same parameters as the originals; ~2.5h EACH (~10h total),
# sequential to avoid oversubscribing cores (each already uses -nthreads 10).
# These feed ONLY the per-SNP Manhattan plots -- the statistical
# recalibration (floor-survivors, Module C) uses the much-faster Stage-1-
# unit scans and does not wait on this.
#
# Pre-fix outputs are quarantined into stale_pre_fix_20260912/ first.
#
# Run (e.g. on mini2, detached):
#   nohup bash module_manuscript_rho05/baypass_stage1/rerun_observed_standardized_fullsnp.sh \
#     > module_manuscript_rho05/baypass_stage1/rerun_fullsnp_standardized.log 2>&1 &

set -euo pipefail
cd "$(dirname "$0")"

PATH_TO_BAYPASS=/Users/petrikem/gitlab/baypass_public-master/sources/g_baypass
CORES=10
D=aland_excluded
STALE=../stale_pre_fix_20260912/baypass_stage1/${D}
mkdir -p "$STALE"

for v in PC1 PC2 bio6 bio11; do
  echo "=== quarantining pre-fix ${v} full-SNP outputs ==="
  cp -f ${D}/${v}_fullSNP_stage1Omega_summary_*.out ${D}/${v}_fullSNP_stage1Omega_baypass.log "$STALE"/ 2>/dev/null || true

  echo "=== ${v} association, FULL SNP set, WITH Stage-1 Omega (standardized covariate) ==="
  "${PATH_TO_BAYPASS}" \
    -countdatafile "${D}/u_DIEM.geno" \
    -omegafile     "${D}/omega_mat_omega.out" \
    -efile         "${D}/u.${v}" \
    -poolsizefile  "${D}/u_DIEM.size" \
    -nthreads "${CORES}" \
    -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 \
    -outprefix "${D}/${v}_fullSNP_stage1Omega"
  [ -s "${D}/${v}_fullSNP_stage1Omega_summary_betai_reg.out" ] || { echo "FAILED: ${v} fullSNP output missing/empty"; exit 1; }
done

echo "=== Done. Rerun outputs (standardized covariates) in ${D}/: {PC1,PC2,bio6,bio11}_fullSNP_stage1Omega_* ==="
