#!/usr/bin/env bash
# eMLG-level BayPass climate association (primary config: Aland-excluded,
# fixed LD-pruned Omega). One "marker" per LD-cluster eMLG.
#
# Inputs are built by R_PK/moduleB_write_eMLG_baypass_inputs.R into
# aland_excluded_eMLG/ (u_eMLG.geno + staged omega_mat_omega.out, u.PC1,
# u.PC2, u_DIEM.size). Parameters match the SNP-level withOmega runs exactly
# (aland_excluded/run_baypass.sh steps 4-5): -nocovscaling, nval 500,
# burnin 5000, thin 25, seed 74, Omega supplied (not re-estimated).
#
# Run from the formica_hybrid repo root.
set -euo pipefail

PATH_TO_BAYPASS=/Users/petrikem/gitlab/baypass_public-master/sources/g_baypass
CORES=8
D=aland_excluded_eMLG

echo "=== PC1 eMLG association (withOmega) ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile "${D}/u_eMLG.geno" \
  -omegafile     "${D}/omega_mat_omega.out" \
  -efile         "${D}/u.PC1" \
  -poolsizefile  "${D}/u_DIEM.size" \
  -nthreads "${CORES}" \
  -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 \
  -outprefix "${D}/PC1_eMLG_withOmega"

echo "=== PC2 eMLG association (withOmega) ==="
"${PATH_TO_BAYPASS}" \
  -countdatafile "${D}/u_eMLG.geno" \
  -omegafile     "${D}/omega_mat_omega.out" \
  -efile         "${D}/u.PC2" \
  -poolsizefile  "${D}/u_DIEM.size" \
  -nthreads "${CORES}" \
  -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 \
  -outprefix "${D}/PC2_eMLG_withOmega"

echo "=== Done. Outputs: ${D}/PC{1,2}_eMLG_withOmega_summary_betai_reg.out ==="
