#!/usr/bin/env bash
# BayPass runs for the Stage-1-direct methodology (see
# module_manuscript_rho05/R/moduleB_stage1_prepare_baypass_inputs.R and
# moduleB_stage1_prepare_mito_contrast.R for how these inputs were built).
#
# Sequence:
#   1. Estimate Omega from ALL Stage-1 cluster representatives (698,251
#      markers, aland_excluded/u_DIEM.geno_pruned)
#   2. PC1 and PC2 association on the Stage-1 units with >=5 loci (18,361,
#      best-SNP represented), WITH the Omega from step 1
#   3. Mitotype C2 contrast on the same Stage-1 units, WITH the same Omega
#      (core model + -contrastfile, no -efile -- see moduleB_stage1_prepare_
#      mito_contrast.R for the contrast coding)
#
# EDIT PATH_TO_BAYPASS for the machine you're running on.

set -euo pipefail
cd "$(dirname "$0")"

PATH_TO_BAYPASS=/Users/petrikem/gitlab/baypass_public-master/sources/g_baypass
CORES=10
OMEGA_DIR=aland_excluded
UNIT_DIR=aland_excluded_S1units

## AUDIT FIX (Issue 9): this Omega estimation is FROZEN as of 2026-09-12 --
## the Omega on disk was NOT re-estimated by the covariate-scaling fix
## (Issue 1) or any other audit change; checksum recorded in AUDIT_FIXES.md.
## -seed added below for reproducibility of any FUTURE re-estimation only --
## it does not apply retroactively to the frozen Omega already on disk.
echo "=== 1. Estimating Omega (Stage-1 pruned, all 698,251 core_snp) [DO NOT RUN -- Omega is frozen; see AUDIT_FIXES.md] ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile "${OMEGA_DIR}/u_DIEM.geno_pruned" \
  -poolsizefile  "${OMEGA_DIR}/u_DIEM.size" \
  -nthreads "${CORES}" \
  -nval 500 -burnin 5000 -thin 10 -seed 74 \
  -outprefix "${OMEGA_DIR}/omega"

cp "${OMEGA_DIR}/omega_mat_omega.out" "${UNIT_DIR}/omega_mat_omega.out"

echo "=== 2. PC1 association, Stage-1 units (n_loci>=5, best-SNP), WITH Omega ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile "${UNIT_DIR}/u_S1units.geno" \
  -omegafile     "${UNIT_DIR}/omega_mat_omega.out" \
  -efile         "${UNIT_DIR}/u.PC1" \
  -poolsizefile  "${UNIT_DIR}/u_DIEM.size" \
  -nthreads "${CORES}" \
  -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 \
  -outprefix "${UNIT_DIR}/PC1_S1units_withOmega"

echo "=== 3. PC2 association, Stage-1 units, WITH Omega ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile "${UNIT_DIR}/u_S1units.geno" \
  -omegafile     "${UNIT_DIR}/omega_mat_omega.out" \
  -efile         "${UNIT_DIR}/u.PC2" \
  -poolsizefile  "${UNIT_DIR}/u_DIEM.size" \
  -nthreads "${CORES}" \
  -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 \
  -outprefix "${UNIT_DIR}/PC2_S1units_withOmega"

echo "=== 4. Mitotype C2 contrast, Stage-1 units, WITH Omega (core model) ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile "${UNIT_DIR}/u_S1units.geno" \
  -omegafile     "${UNIT_DIR}/omega_mat_omega.out" \
  -contrastfile  "${UNIT_DIR}/u.mito_contrast" \
  -poolsizefile  "${UNIT_DIR}/u_DIEM.size" \
  -nthreads "${CORES}" \
  -nval 500 -burnin 5000 -thin 25 -seed 74 \
  -outprefix "${UNIT_DIR}/mito_C2_S1units"

echo "=== Done. Outputs in ${UNIT_DIR}/: PC1_S1units_withOmega_*, PC2_S1units_withOmega_*, mito_C2_S1units_summary_contrast.out ==="
