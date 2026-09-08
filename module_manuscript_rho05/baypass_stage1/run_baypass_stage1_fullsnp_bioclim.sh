#!/usr/bin/env bash
# Full per-SNP BayPass scan (all 1,114,423 SNPs, aland_excluded/u_DIEM.geno)
# for the winter-temperature covariates (bio6, bio11), WITH the Stage-1-
# derived Omega -- same rationale/structure as run_baypass_stage1_fullsnp.sh
# (PC1/PC2/mito-C2): genuine per-SNP statistics for the Manhattan plot,
# not values inherited from a SNP's Stage-1 cluster.
#
# EDIT PATH_TO_BAYPASS for the machine you're running on.

set -euo pipefail
cd "$(dirname "$0")"

PATH_TO_BAYPASS=/Users/petrikem/gitlab/baypass_public-master/sources/g_baypass
CORES=10
D=aland_excluded

echo "=== 1. bio6 association, FULL SNP set (1,114,423), WITH Stage-1 Omega ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile "${D}/u_DIEM.geno" \
  -omegafile     "${D}/omega_mat_omega.out" \
  -efile         "${D}/u.bio6" \
  -poolsizefile  "${D}/u_DIEM.size" \
  -nthreads "${CORES}" \
  -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 \
  -outprefix "${D}/bio6_fullSNP_stage1Omega"

echo "=== 2. bio11 association, FULL SNP set, WITH Stage-1 Omega ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile "${D}/u_DIEM.geno" \
  -omegafile     "${D}/omega_mat_omega.out" \
  -efile         "${D}/u.bio11" \
  -poolsizefile  "${D}/u_DIEM.size" \
  -nthreads "${CORES}" \
  -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 \
  -outprefix "${D}/bio11_fullSNP_stage1Omega"

echo "=== Done. Outputs in ${D}/: bio6_fullSNP_stage1Omega_*, bio11_fullSNP_stage1Omega_* ==="
