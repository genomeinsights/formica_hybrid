#!/bin/bash
# =========================================================================
# One founding group x replicate of the neutral mosaic-founder simulations
# (the per-population design of Beatriz's sim_table.txt).
#   run_group_job.sh <GROUP> <RUN> <BASE_DIR>
# GROUP: one of the 18 founding groups below. lan = lanR + lanW and bungrund = bun + grund
#        share founders and FOUNDER_SEED, as in sim_table.txt.
# RUN  : replicate number (1, 2, ...).
# BASE_DIR must contain sim_founder_fix/ (this repository folder) and is where
#        founders/, out/, results/ and logs/ are written.
#
# Settings (environment variables, defaults in brackets):
#   K      carrying capacity                         [6250]
#   INITN  number of founders (half aquilonia males,
#          half polyctena queens)                    [100]
#   CYCLE  sampled cycle; Sielva always cycle 11     [125]
#   PHASED phased parents                            [$BASE/sim_founder_fix/data/parents_phased.rds]
#   RECDIR SLiM recombination maps (ch_<id>.recmap)  [$BASE/slim_inputs/recombination_maps/]
#   CL     climate vector file                       [$BASE/slim_inputs/climate/climate_rep1.txt]
#   KEEP   females | both (also keep sampled males)  [females]
# Example: K=12500 INITN=100 sim_founder_fix/slim/run_group_job.sh aland 1 $PWD
#
# Steps: generate the group's founder pool (genotype mosaics of the phased empirical
# parents; ancestry-informative + near-neutral SNPs; seed specific to group and run), run
# SLiM once per member population, keep the sampled new queens (gzipped VCF) and the
# ancestry summary, delete everything else. Idempotent: a population whose kept file
# exists is skipped.
# =========================================================================
set -euo pipefail
GROUP=$1; RUN=$2; BASE=$3
cd "$BASE"
export PATH=/usr/local/bin:/opt/homebrew/bin:$PATH
K=${K:-6250}; INITN=${INITN:-100}; CYCLE_DEFAULT=${CYCLE:-125}; KEEP=${KEEP:-females}
PHASED=${PHASED:-$BASE/sim_founder_fix/data/parents_phased.rds}
RECDIR=${RECDIR:-$BASE/slim_inputs/recombination_maps/}
CL=${CL:-$BASE/slim_inputs/climate/climate_rep1.txt}
for f in "$PHASED" "$CL" "${RECDIR%/}/ch_1.recmap"; do [ -s "$f" ] || { echo "missing input: $f"; exit 1; }; done
mkdir -p founders out results logs

case "$GROUP" in
  lan)      MEMBERS="lanR lanW" ;;
  bungrund) MEMBERS="bun grund" ;;
  *)        MEMBERS="$GROUP" ;;
esac
GIDX=$(printf "%s\n" aland katis lan svan1 svan2 tvar bungrund pik nyr1 nyr2 heina pari hiiv vuos kumm karsi jarven sielva | grep -nx "$GROUP" | cut -d: -f1)
[ -n "$GIDX" ] || { echo "unknown group $GROUP"; exit 1; }
SEED=$(( 100000 + GIDX * 1000 + RUN ))          # founder-pool seed = FOUNDER_SEED
cyc() { if [ "$1" = "sielva" ]; then echo 11; else echo "$CYCLE_DEFAULT"; fi; }

todo=""
for POP in $MEMBERS; do
  TAG=$(printf "%s_%02d" "$POP" "$RUN")
  [ -s "results/females_ckl$(cyc $POP)_${TAG}.vcf.gz" ] || todo="$todo $POP"
done
[ -z "$todo" ] && { echo "[$GROUP $RUN] done already"; exit 0; }

NF=$(( INITN / 2 ))
POOL="founders/pool_${GROUP}_$(printf %02d "$RUN")"
if [ ! -s "$POOL/founders_ch27.vcf" ]; then
  Rscript sim_founder_fix/make_mosaic_founders.R "$POOL" 1 "$NF" "$NF" "$SEED" DI25+neutral "$PHASED" > "logs/gen_${GROUP}_${RUN}.log" 2>&1
fi

for POP in $todo; do
  TAG=$(printf "%s_%02d" "$POP" "$RUN"); C=$(cyc $POP)
  OUT="out/$TAG/"; mkdir -p "$OUT"
  start=$(date +%s)
  # stop at the sampled cycle: nothing after it is used
  slim -d "rep=$SEED" -d "TAG=\"$TAG\"" -d "FOUNDER_SEED=$SEED" -d "K=$K" -d "initialN=$INITN" \
       -d "nCycles=$C" -d "sampleCycle=c($C)" \
       -d "FDIR=\"$BASE/$POOL/\"" -d "RECDIR=\"$RECDIR\"" -d "CL=\"$CL\"" -d "folder=\"$BASE/$OUT\"" \
       sim_founder_fix/slim/SpecIAnt_rufa_mosaic_founders.slim > "logs/slim_${TAG}.log" 2>&1
  for SEX in females $( [ "$KEEP" = both ] && echo males ); do
    gzip -c "$OUT/${SEX}_ckl${C}_${TAG}.vcf" > "results/${SEX}_ckl${C}_${TAG}.vcf.gz.tmp"
    mv "results/${SEX}_ckl${C}_${TAG}.vcf.gz.tmp" "results/${SEX}_ckl${C}_${TAG}.vcf.gz"
  done
  cp "$OUT"/*.anc "results/${TAG}.anc" 2>/dev/null || true
  rm -rf "$OUT"
  # compact log: drop SLiM's echo of the genomic-element set-up
  grep -v "^initializeGenomicElement" "logs/slim_${TAG}.log" | tail -40 > "logs/slim_${TAG}.tail"; rm -f "logs/slim_${TAG}.log"
  echo "[$GROUP $RUN] $TAG done in $(( $(date +%s) - start )) s"
done
rm -rf "$POOL"
