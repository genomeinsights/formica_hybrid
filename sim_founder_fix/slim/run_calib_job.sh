#!/bin/bash
# =========================================================================
# One calibration run of the neutral mosaic-founder model (phased founders).
#   run_calib_job.sh <SETTING> <RUN> <BASE_DIR>
# SETTING: a row name of slim/calib_settings.txt (columns: name K initialN).
# Each run = one independent hybrid population (own founder pool, no pairing),
# simulated on chromosomes CHROMS only and sampled at several cycles, so a single
# run gives every time point of the grid. F_ST is computed later among the runs
# of one setting. Keeps results/calib/females_ckl<C>_<SETTING>_<RUN>.vcf.gz.
# Idempotent: skipped if the last cycle's file exists.
# =========================================================================
set -euo pipefail
SETTING=$1; RUN=$2; BASE=$3
cd "$BASE"
export PATH=/usr/local/bin:/opt/homebrew/bin:$PATH
CHROMS="c(1:6)"
CYCLES="60 125 250 500 1000"
LAST=1000
PHASED="$BASE/phased/parents_phased.rds"

read -r _ K INITN < <(awk -v s="$SETTING" '$1 == s' sim_founder_fix/slim/calib_settings.txt)
[ -n "${INITN:-}" ] || { echo "unknown setting $SETTING"; exit 1; }
SIDX=$(awk -v s="$SETTING" '$1 == s {print NR}' sim_founder_fix/slim/calib_settings.txt)
SEED=$(( 200000 + SIDX * 1000 + RUN ))
TAG=$(printf "%s_%02d" "$SETTING" "$RUN")
mkdir -p results/calib logs/calib founders
[ -s "results/calib/females_ckl${LAST}_${TAG}.vcf.gz" ] && { echo "[$TAG] done already"; exit 0; }

NF=$(( INITN / 2 ))                      # initialN/2 aquilonia males + initialN/2 polyctena queens
POOL="founders/calib_${TAG}"
if [ ! -s "$POOL/founders_ch27.vcf" ]; then
  Rscript sim_founder_fix/make_mosaic_founders.R "$POOL" 1 "$NF" "$NF" "$SEED" DI25+neutral "$PHASED" > "logs/calib/gen_${TAG}.log" 2>&1
fi

OUT="out/calib_$TAG/"; mkdir -p "$OUT"
SC="c($(echo $CYCLES | tr ' ' ','))"
start=$(date +%s)
/usr/bin/time -l slim -d "rep=$SEED" -d "TAG=\"$TAG\"" -d "FOUNDER_SEED=$SEED" \
     -d "K=$K" -d "initialN=$INITN" -d "CHROMS=$CHROMS" \
     -d "nCycles=$LAST" -d "sampleCycle=$SC" \
     -d "FDIR=\"$BASE/$POOL/\"" -d "RECDIR=\"$BASE/slim_inputs/recombination_maps/\"" \
     -d "CL=\"$BASE/slim_inputs/climate/climate_rep1.txt\"" -d "folder=\"$BASE/$OUT\"" \
     sim_founder_fix/slim/SpecIAnt_rufa_mosaic_founders.slim > "logs/calib/slim_${TAG}.log" 2>&1
for C in $CYCLES; do
  gzip -c "$OUT/females_ckl${C}_${TAG}.vcf" > "results/calib/females_ckl${C}_${TAG}.vcf.gz.tmp"
  mv "results/calib/females_ckl${C}_${TAG}.vcf.gz.tmp" "results/calib/females_ckl${C}_${TAG}.vcf.gz"
done
cp "$OUT"/*.anc "results/calib/${TAG}.anc" 2>/dev/null || true
rm -rf "$OUT" "$POOL"
grep -v "^initializeGenomicElement" "logs/calib/slim_${TAG}.log" | tail -40 > "logs/calib/slim_${TAG}.tail"; rm -f "logs/calib/slim_${TAG}.log"
echo "[$TAG] done in $(( $(date +%s) - start )) s"
