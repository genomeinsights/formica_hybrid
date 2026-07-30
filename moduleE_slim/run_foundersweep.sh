#!/bin/bash -l
## =============================================================================
## Module E -- SHALLOW-BOTTLENECK founder-number sweep (fixed K, inflated pool)
## =============================================================================
## Tests the hypothesis that the neutral null's excess F_ST / background LD comes
## from the severe founding bottleneck (~43 real founders -> K in a few gens).
## Each real founder haplotype is copied C times (make_founders_inflated.R); the
## model then draws N_AQ/N_POL WITHOUT replacement from the inflated pool, which is
## a resample WITH replacement from the 30/13 real founders. Larger N_AQ/N_POL =>
## shallower bottleneck => lower founding F_ST and less drift-LD (diversity still
## capped at the real 30/13 haplotypes). 20 independent demes per setting.
##
## Config (env vars, all optional):
##   NPERBIN  markers per DI bin      (default 1500 -> ~15k markers, ~2.9 GB/deme @K6250)
##   COPIES   founder copies C        (default 40   -> pool 1200 aq / 520 pol)
##   FLIST    "N_AQ N_POL" pairs ';'  (default "30 13;150 65;450 195;1200 520")
##   K        carrying capacity       (default 6250)
##   NCYC     generations             (default 80; sampled every 20 -> gens 20..80)
##   OUT      output dir              (default output_foundersweep)
##   SS       sample size per deme    (default 50)
## Usage:  bash run_foundersweep.sh [CONC]      (CONC = parallel demes, default 6)
##
## After: analyse each founding tag, e.g.
##   Rscript scripts/analyze_di_stratified.R output_foundersweep Naq450_Npol195
## and read the NEUTRAL-bin F_ST + pi vs empirical (target F_ST~0, pi~0.448).
## =============================================================================
set -u
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODEL="$SD/SpecIAnt_rufa_neutral_realfounders.slim"
RECDIR="$SD/rec_maps/"
NPERBIN="${NPERBIN:-1500}"; COPIES="${COPIES:-40}"
K="${K:-6250}"; NCYC="${NCYC:-80}"; EVERY=20; SS="${SS:-50}"; NDEMES=20
OUT="${OUT:-$SD/output_foundersweep}"
CONC="${1:-6}"
FDIR="${FDIR:-$SD/founders/maf015_DIstrat${NPERBIN}_x${COPIES}/}"   # override to use e.g. the burn-in founders

command -v slim >/dev/null || { echo "ERROR: 'slim' not on PATH (need SLiM 5)."; exit 1; }
NAQP=$(( 30 * COPIES )); NPOLP=$(( 13 * COPIES ))     # inflated pool sizes

## build the inflated founders once if absent
if [ ! -f "${FDIR}founders_ch1.vcf" ]; then
  echo "building inflated founders (n_per_bin=$NPERBIN, copies=$COPIES) ..."
  Rscript "$SD/scripts/make_founders_inflated.R" "$NPERBIN" "$COPIES" || exit 1
fi

mkdir -p "$OUT" "$OUT/logs"
if [ -n "${FLIST:-}" ]; then IFS=';' read -ra FOUNDINGS <<< "$FLIST"
else FOUNDINGS=("30 13" "150 65" "450 195" "1200 520"); fi
KS=(${KLIST:-$K})                                    # one or more K values

JOBS=$(for f in "${FOUNDINGS[@]}"; do read -r naq npol <<< "$f"
         for Kv in "${KS[@]}"; do for r in $(seq 1 $NDEMES); do echo "$naq $npol $Kv $r"; done; done
       done)
echo "runs: $(echo "$JOBS" | wc -l | tr -d ' ')  foundings: ${FOUNDINGS[*]}  K: ${KS[*]}  pool=${NAQP}/${NPOLP}  CONC=$CONC"
echo "founders: $FDIR"; echo "output:   $OUT"

run_one() {
  read -r naq npol Kv r <<< "$1"
  local tag="Naq${naq}_Npol${npol}_K${Kv}_rep${r}"
  local seed=$(( naq*1000000 + npol*100000 + Kv*10 + r ))
  slim -s "$seed" -d "folder='$OUT'" -d "FDIR='$FDIR'" -d "RECDIR='$RECDIR'" \
    -d "N_AQ_POOL=$NAQP" -d "N_POL_POOL=$NPOLP" -d "N_AQ=$naq" -d "N_POL=$npol" \
    -d "K=$Kv" -d "nCycles=$NCYC" -d "sampleEvery=$EVERY" -d "sampleSize=$SS" -d "rep=$r" \
    "$MODEL" > "$OUT/logs/${tag}.log" 2>&1
  echo "[$(date +%H:%M:%S)] done $tag (exit $?)"
}
export -f run_one; export MODEL RECDIR FDIR OUT NCYC EVERY SS NAQP NPOLP
echo "$JOBS" | xargs -P "$CONC" -I {} bash -c 'run_one "$@"' _ {}
echo "ALL DONE. VCFs: $(ls "$OUT"/*.vcf 2>/dev/null | wc -l | tr -d ' ') ; errors: $(grep -lri error "$OUT/logs/" 2>/dev/null | wc -l | tr -d ' ')"
