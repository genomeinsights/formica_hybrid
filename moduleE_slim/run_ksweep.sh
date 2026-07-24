#!/bin/bash -l
## =============================================================================
## Module E -- K / founder sweep (PORTABLE, bundle-relative)
## =============================================================================
## Runs the neutral real-founder model over a grid of carrying capacity (K) and
## founder number, 20 independent demes each, sampled every 20 generations to 160.
## All paths resolve relative to this bundle. Requires SLiM 5 on PATH.
##
## Filenames encode founding and K (Naq<..>_Npol<..>_K<..>), so different
## parameter sets never collide and can share one output directory.
##
## Config via environment variables (all optional):
##   FLIST   space-separated "N_AQ N_POL" pairs, one per ';'  e.g. "30 13;12 6"
##   KLIST   space-separated K values                          e.g. "6250 12500"
##   FDIR    founder-VCF directory (default founders/maf015_DIstrat4000)
##   OUT     output directory       (default output)
##   SS      sample size per deme    (default 50)
## Usage:  KLIST="6250 12500" bash run_ksweep.sh [CONC]
##
## MEMORY: ~0.19 GB per 1000 markers at K=6250, scaling ~linearly with K. The
## 40k-marker set is ~7.6 GB/deme at K=6250 and ~15 GB at K=12500 -- set CONC so
## CONC * (per-deme GB) stays under the machine's RAM.
## =============================================================================
set -u
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

MODEL="$SD/SpecIAnt_rufa_neutral_realfounders.slim"
RECDIR="$SD/rec_maps/"
FDIR="${FDIR:-$SD/founders/maf015_DIstrat4000/}"
OUT="${OUT:-$SD/output}"
NCYC=160; EVERY=20; SS="${SS:-50}"; NDEMES=20
CONC="${1:-4}"

command -v slim >/dev/null || { echo "ERROR: 'slim' not on PATH (need SLiM 5)."; exit 1; }
mkdir -p "$OUT" "$OUT/logs"

## founder settings: default full pool 30/13; override with FLIST="30 13;12 6;6 3"
if [ -n "${FLIST:-}" ]; then IFS=';' read -ra FOUNDINGS <<< "$FLIST"; else FOUNDINGS=("30 13"); fi
KS=(${KLIST:-6250 12500})

JOBS=$(for f in "${FOUNDINGS[@]}"; do read -r naq npol <<< "$f"
         for K in "${KS[@]}"; do for r in $(seq 1 $NDEMES); do echo "$naq $npol $K $r"; done; done
       done)
echo "runs: $(echo "$JOBS" | wc -l | tr -d ' ')  foundings: ${FOUNDINGS[*]}  K: ${KS[*]}  CONC=$CONC"
echo "founders: $FDIR"; echo "output:   $OUT"

run_one() {
  read -r naq npol K r <<< "$1"
  local tag="Naq${naq}_Npol${npol}_K${K}_rep${r}"; local seed=$(( naq*1000000 + npol*100000 + K*10 + r ))
  slim -s "$seed" -d "folder='$OUT'" -d "FDIR='$FDIR'" -d "RECDIR='$RECDIR'" \
    -d "N_AQ=$naq" -d "N_POL=$npol" -d "K=$K" -d "nCycles=$NCYC" \
    -d "sampleEvery=$EVERY" -d "sampleSize=$SS" -d "rep=$r" \
    "$MODEL" > "$OUT/logs/${tag}.log" 2>&1
  echo "[$(date +%H:%M:%S)] done $tag (exit $?)"
}
export -f run_one; export MODEL RECDIR FDIR OUT NCYC EVERY SS
echo "$JOBS" | xargs -P "$CONC" -I {} bash -c 'run_one "$@"' _ {}
echo "ALL DONE. VCFs: $(ls "$OUT"/*.vcf 2>/dev/null | wc -l | tr -d ' ') ; errors: $(grep -lri error "$OUT/logs/" 2>/dev/null | wc -l | tr -d ' ')"
