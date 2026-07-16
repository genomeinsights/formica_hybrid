#!/usr/bin/env bash
# BayPass runs for the full dataset (all populations, including Aland) --
# the comparison baseline for ../aland_excluded/run_baypass.sh.
#
# Expects to live in the same folder as its input files (u_DIEM.geno,
# u_DIEM.geno_pruned, u_DIEM.size, u.PC1, u.PC2 -- all written by
# R/prepare_with_aland.R). Run from anywhere; it cds into its own
# directory first.
#
# Sequence, matching baypass.R's parameters exactly:
#   1. Estimate Omega from the PRUNED genotype set (core model, no covariate)
#   2. PC1 and PC2 association WITHOUT a pre-estimated Omega (BayPass
#      estimates its own from the full, unpruned genotype set during the run)
#   3. PC1 and PC2 association WITH the Omega estimated in step 1
#
# "Without omega" and "with omega" runs use different -outprefix values so
# neither overwrites the other's output.

set -euo pipefail
cd "$(dirname "$0")"

## EDIT this for the machine you're running on -- do not assume the path
## used on the original machine is valid here.
PATH_TO_BAYPASS="~/gitlab/baypass_public-master/sources/g_baypass"
CORES=10

echo "=== 1. Estimating Omega (pruned data, Aland excluded) ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile u_DIEM.geno_pruned \
  -poolsizefile u_DIEM.size \
  -nthreads "${CORES}" \
  -nval 500 -burnin 5000 -thin 10 \
  -outprefix omega

echo "=== 2. PC1 association, WITHOUT pre-estimated Omega ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile u_DIEM.geno \
  -efile u.PC1 \
  -poolsizefile u_DIEM.size \
  -nthreads "${CORES}" \
  -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 \
  -outprefix PC1_DIEM_noOmega

echo "=== 3. PC2 association, WITHOUT pre-estimated Omega ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile u_DIEM.geno \
  -efile u.PC2 \
  -poolsizefile u_DIEM.size \
  -nthreads "${CORES}" \
  -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 \
  -outprefix PC2_DIEM_noOmega

echo "=== 4. PC1 association, WITH Omega from step 1 (pruned data) ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile u_DIEM.geno \
  -omegafile omega_mat_omega.out \
  -efile u.PC1 \
  -poolsizefile u_DIEM.size \
  -nthreads "${CORES}" \
  -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 \
  -outprefix PC1_DIEM_withOmega

echo "=== 5. PC2 association, WITH Omega from step 1 (pruned data) ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile u_DIEM.geno \
  -omegafile omega_mat_omega.out \
  -efile u.PC2 \
  -poolsizefile u_DIEM.size \
  -nthreads "${CORES}" \
  -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 \
  -outprefix PC2_DIEM_withOmega

echo "=== Done. Outputs: omega_*, PC1_DIEM_noOmega_*, PC2_DIEM_noOmega_*, PC1_DIEM_withOmega_*, PC2_DIEM_withOmega_* ==="
