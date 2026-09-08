#!/usr/bin/env bash
# Full per-SNP BayPass scan (all 1,114,423 SNPs, aland_excluded/u_DIEM.geno),
# WITH the Stage-1-derived Omega (aland_excluded/omega_mat_omega.out,
# estimated from all 698,251 Stage-1 cluster representatives -- see
# moduleB_stage1_prepare_baypass_inputs.R). This gives genuine per-SNP
# statistics, fixing the earlier plot's problem: every SNP inheriting its
# Stage-1 cluster's single tested value (flat, unit-resolution) rather than
# an independently computed one.
#
# Sequential (not parallel) to avoid oversubscribing this machine's cores --
# each run already uses -nthreads 10 near the machine's 11-core limit.
#
# EDIT PATH_TO_BAYPASS for the machine you're running on.

set -euo pipefail
cd "$(dirname "$0")"

PATH_TO_BAYPASS=/Users/petrikem/gitlab/baypass_public-master/sources/g_baypass
CORES=10
D=aland_excluded

echo "=== 1. PC1 association, FULL SNP set (1,114,423), WITH Stage-1 Omega ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile "${D}/u_DIEM.geno" \
  -omegafile     "${D}/omega_mat_omega.out" \
  -efile         "${D}/u.PC1" \
  -poolsizefile  "${D}/u_DIEM.size" \
  -nthreads "${CORES}" \
  -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 \
  -outprefix "${D}/PC1_fullSNP_stage1Omega"

echo "=== 2. PC2 association, FULL SNP set, WITH Stage-1 Omega ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile "${D}/u_DIEM.geno" \
  -omegafile     "${D}/omega_mat_omega.out" \
  -efile         "${D}/u.PC2" \
  -poolsizefile  "${D}/u_DIEM.size" \
  -nthreads "${CORES}" \
  -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 \
  -outprefix "${D}/PC2_fullSNP_stage1Omega"

echo "=== 3. Mitotype C2 contrast, FULL SNP set, WITH Stage-1 Omega ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile "${D}/u_DIEM.geno" \
  -omegafile     "${D}/omega_mat_omega.out" \
  -contrastfile  "${D}/u.mito_contrast" \
  -poolsizefile  "${D}/u_DIEM.size" \
  -nthreads "${CORES}" \
  -nval 500 -burnin 5000 -thin 25 -seed 74 \
  -outprefix "${D}/mito_C2_fullSNP_stage1Omega"

echo "=== Done. Outputs in ${D}/: PC1_fullSNP_stage1Omega_*, PC2_fullSNP_stage1Omega_*, mito_C2_fullSNP_stage1Omega_summary_contrast.out ==="
