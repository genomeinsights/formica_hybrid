#!/usr/bin/env bash
## =========================================================
## module_localscore_crosscheck -- one-replicate end-to-end pilot
## =========================================================
## Runs BayPass on the full SNP set (1,114,423 SNPs, u_DIEM.geno) with the
## authoritative Stage-1-derived Omega, for BOTH null types of replicate 1
## (continuous mode only -- the pilot's purpose is to validate the mechanics
## end-to-end, not to pre-empt the full 10x2x2 sweep):
##   - structured draw 1   (seed 20101)
##   - unstructured draw 1 (seed 20201, permutation of the structured draw)
##
## NOT YET RUN. Reviewed by the user before execution, per the explicit
## checkpoint: "stop after completing a one-replicate end-to-end pilot ...
## before launching expensive BayPass work."
##
## EDIT PATH_TO_BAYPASS for the machine you're running on (mini2).
## =========================================================

set -euo pipefail
cd "$(dirname "$0")/../.."   # repo root

PATH_TO_BAYPASS=/Users/petrikem/gitlab/baypass_public-master/sources/g_baypass
CORES=10
D=module_manuscript_rho05/baypass_stage1/aland_excluded   # authoritative geno/Omega/poolsize
OUT=module_localscore_crosscheck/full_snp_null10/raw
PILOT=module_localscore_crosscheck/full_snp_null10/pilot
mkdir -p "${OUT}"

echo "=== Pilot 1/2: structured draw 1, continuous, FULL SNP set ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile "${D}/u_DIEM.geno" \
  -omegafile     "${D}/omega_mat_omega.out" \
  -efile         "${PILOT}/pilot_structured_draw1.env" \
  -poolsizefile  "${D}/u_DIEM.size" \
  -nthreads "${CORES}" \
  -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 20101 \
  -outprefix "${OUT}/pilot_structured_draw1"

echo "=== Pilot 2/2: unstructured draw 1, continuous, FULL SNP set ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile "${D}/u_DIEM.geno" \
  -omegafile     "${D}/omega_mat_omega.out" \
  -efile         "${PILOT}/pilot_unstructured_draw1.env" \
  -poolsizefile  "${D}/u_DIEM.size" \
  -nthreads "${CORES}" \
  -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 20201 \
  -outprefix "${OUT}/pilot_unstructured_draw1"

echo "=== Pilot done. Outputs in ${OUT}/: pilot_structured_draw1_*, pilot_unstructured_draw1_* ==="
