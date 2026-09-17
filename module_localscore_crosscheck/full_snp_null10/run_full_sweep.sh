#!/usr/bin/env bash
## =========================================================
## module_localscore_crosscheck -- full 10-replicate x 2-mode x 2-null-type
## full-SNP BayPass sweep (40 runs total). Run on mini2; raw output goes to
## an external SSD (/Volumes/T9), not mini2's internal disk.
##
## Reuses the exact BayPass settings validated by the one-replicate pilot
## (structured/unstructured draw 1, continuous mode: ~1h24m30s/run, sane BF
## distributions, row counts matching the full SNP set). Continuous runs use
## -nocovscaling (matching the precedent's PC1/PC2 full-SNP calls); mitoC2
## contrast runs do not (matching the precedent's mito_C2 full-SNP call).
##
## Seeds: BayPass -seed reuses each replicate's own covariate-generation
## seed (STRUCTURED_SEEDS[i] / PERMUTE_SEEDS[i] from
## config/null_covariate_design.md) -- documented, not arbitrary, and
## distinct per replicate/null-type.
## =========================================================

set -euo pipefail
mkdir -p ~/module_localscore_crosscheck_full_sweep
cd ~/module_localscore_crosscheck_full_sweep

BAYPASS=~/baypass_public/sources/g_baypass
D=~/formica_hybrid/baypass_stage1_fullsnp        # authoritative geno/Omega/poolsize (MD5-verified)
REP_DIR=~/module_localscore_crosscheck_full_sweep/replicates
RAW=/Volumes/T9/module_localscore_crosscheck_full_snp_null10/raw
CORES=10
mkdir -p "${RAW}"

STRUCTURED_SEEDS=(20101 20102 20103 20104 20105 20106 20107 20108 20109 20110)
PERMUTE_SEEDS=(20201 20202 20203 20204 20205 20206 20207 20208 20209 20210)

run_continuous () {
  local nulltype=$1 efile=$2 seed=$3 outprefix=$4
  echo "=== [$(date '+%F %T')] continuous ${nulltype}: ${outprefix} (seed ${seed}) ==="
  "${BAYPASS}" \
    -countdatafile "${D}/u_DIEM.geno" \
    -omegafile     "${D}/omega_mat_omega.out" \
    -efile         "${efile}" \
    -poolsizefile  "${D}/u_DIEM.size" \
    -nthreads "${CORES}" \
    -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed "${seed}" \
    -outprefix "${RAW}/${outprefix}"
}

run_mitoC2 () {
  local nulltype=$1 cfile=$2 seed=$3 outprefix=$4
  echo "=== [$(date '+%F %T')] mitoC2 ${nulltype}: ${outprefix} (seed ${seed}) ==="
  "${BAYPASS}" \
    -countdatafile "${D}/u_DIEM.geno" \
    -omegafile     "${D}/omega_mat_omega.out" \
    -contrastfile  "${cfile}" \
    -poolsizefile  "${D}/u_DIEM.size" \
    -nthreads "${CORES}" \
    -nval 500 -burnin 5000 -thin 25 -seed "${seed}" \
    -outprefix "${RAW}/${outprefix}"
}

for i in $(seq -w 1 10); do
  idx=$((10#${i}))
  s_seed=${STRUCTURED_SEEDS[$((idx-1))]}
  u_seed=${PERMUTE_SEEDS[$((idx-1))]}

  run_continuous structured   "${REP_DIR}/structured_continuous_draw${i}.txt"   "${s_seed}" "structured_continuous_draw${i}"
  run_continuous unstructured "${REP_DIR}/unstructured_continuous_draw${i}.txt" "${u_seed}" "unstructured_continuous_draw${i}"
  run_mitoC2     structured   "${REP_DIR}/structured_mitoC2_draw${i}.txt"       "${s_seed}" "structured_mitoC2_draw${i}"
  run_mitoC2     unstructured "${REP_DIR}/unstructured_mitoC2_draw${i}.txt"     "${u_seed}" "unstructured_mitoC2_draw${i}"

  echo "=== [$(date '+%F %T')] replicate ${i} done (4/4 runs) ==="
done

echo "=== [$(date '+%F %T')] FULL SWEEP DONE: 40/40 runs ==="
