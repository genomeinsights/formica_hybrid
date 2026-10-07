#!/bin/bash
# =========================================================================
# One founding group x run index of the mosaic-founder simulations.
#   run_group_job.sh <GROUP> <RUN> <BASE_DIR>
# GROUP: one of the 18 founding groups (lan = lanR + lanW, bungrund = bun + grund
# share founders and first-cycle mating, as in Beatriz's sim_table.txt).
# Steps: generate the group's mosaic founder pool (seed = group/run specific, panel
# DI25 + near-neutral), run SLiM once per member population with a shared
# FOUNDER_SEED, keep only the hybrid sample used downstream (new queens at cycle 125;
# Sielva at cycle 11), gzip it, delete everything else (founder pool included).
# Idempotent: a population whose kept file exists is skipped.
# =========================================================================
set -euo pipefail
GROUP=$1; RUN=$2; BASE=$3
cd "$BASE"
export PATH=/usr/local/bin:/opt/homebrew/bin:$PATH

case "$GROUP" in
  lan)      MEMBERS="lanR lanW" ;;
  bungrund) MEMBERS="bun grund" ;;
  *)        MEMBERS="$GROUP" ;;
esac
GIDX=$(printf "%s\n" aland katis lan svan1 svan2 tvar bungrund pik nyr1 nyr2 heina pari hiiv vuos kumm karsi jarven sielva | grep -nx "$GROUP" | cut -d: -f1)
SEED=$(( 100000 + GIDX * 1000 + RUN ))          # founder-pool seed = FOUNDER_SEED

todo=""
for POP in $MEMBERS; do
  TAG=$(printf "%s_%02d" "$POP" "$RUN")
  CYCLE=125; [ "$POP" = "sielva" ] && CYCLE=11
  [ -s "results/females_ckl${CYCLE}_${TAG}.vcf.gz" ] || todo="$todo $POP"
done
[ -z "$todo" ] && { echo "[$GROUP $RUN] done already"; exit 0; }

POOL="founders/pool_${GROUP}_$(printf %02d "$RUN")"
if [ ! -s "$POOL/founders_ch27.vcf" ]; then
  Rscript sim_founder_fix/make_mosaic_founders.R "$POOL" 1 50 50 "$SEED" DI25+neutral > "logs/gen_${GROUP}_${RUN}.log" 2>&1
fi

for POP in $todo; do
  TAG=$(printf "%s_%02d" "$POP" "$RUN")
  CYCLE=125; [ "$POP" = "sielva" ] && CYCLE=11
  OUT="out/$TAG/"; mkdir -p "$OUT"
  start=$(date +%s)
  # stop at the sampled cycle (125, or 11 for Sielva): nothing after it is used
  /usr/bin/time -l slim -d "rep=$SEED" -d "TAG=\"$TAG\"" -d "FOUNDER_SEED=$SEED" \
       -d "nCycles=$CYCLE" -d "sampleCycle=c($CYCLE)" \
       -d "FDIR=\"$BASE/$POOL/\"" -d "RECDIR=\"$BASE/slim_inputs/recombination_maps/\"" \
       -d "CL=\"$BASE/slim_inputs/climate/climate_rep1.txt\"" -d "folder=\"$BASE/$OUT\"" \
       sim_founder_fix/slim/SpecIAnt_rufa_mosaic_founders.slim > "logs/slim_${TAG}.log" 2>&1
  gzip -c "$OUT/females_ckl${CYCLE}_${TAG}.vcf" > "results/females_ckl${CYCLE}_${TAG}.vcf.gz.tmp"
  mv "results/females_ckl${CYCLE}_${TAG}.vcf.gz.tmp" "results/females_ckl${CYCLE}_${TAG}.vcf.gz"
  cp "$OUT"/*.anc "results/${TAG}.anc" 2>/dev/null || true
  rm -rf "$OUT"
  # keep a compact log: drop SLiM's echo of the ~40,000 genomic-element lines
  grep -v "^initializeGenomicElement" "logs/slim_${TAG}.log" | tail -40 > "logs/slim_${TAG}.tail"; rm -f "logs/slim_${TAG}.log"
  echo "[$GROUP $RUN] $TAG done in $(( $(date +%s) - start )) s"
done
rm -rf "$POOL"
