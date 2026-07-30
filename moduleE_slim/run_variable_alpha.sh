#!/bin/bash -l
## =============================================================================
## Module E -- VARIABLE-ALPHA neutral null
## =============================================================================
## Tests the one neutral mechanism that can make F_ST depend on DiagnosticIndex:
## VARIATION IN ADMIXTURE PROPORTION among populations. Each of 20 independent demes
## is founded with a DIFFERENT aquilonia:polyctena ratio, set to a real population's
## estimated ancestry (alpha in [0.35,0.64], mean 0.537 -- matched to empirical). The
## founding counts (N_aq, N_pol per deme) come from moduleE_slim/inputs/varalpha_
## founding.csv (dev/R/moduleE_varalpha_setup.R), chosen to hit each alpha while
## holding total founder haplotypes = 840 so the bottleneck depth is constant and only
## alpha varies. Everything else matches the fixed-alpha replicate control
## (output_replicates: burn-in pool, K=6250, gen 60): the ONLY difference is alpha.
##
## Model is UNCHANGED -- it already takes per-run N_AQ/N_POL. Compare the output against
## (i) the fixed-alpha null (output_replicates/rep01) and (ii) empirical, on F_ST-vs-DI
## and adjacent-LD-vs-DI (dev/R/moduleE_varalpha.R).
##
## Usage: bash run_variable_alpha.sh [CONC]     (CONC = demes in parallel, default 2)
## Output: output_varalpha/hyb_neutral_realfounders_Naq<n>_Npol<n>_K6250_rep<deme>_gen60.vcf
## =============================================================================
set -u
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODEL="$SD/SpecIAnt_rufa_neutral_realfounders.slim"; RECDIR="$SD/rec_maps/"
FDIR="$SD/founders/maf015_DIstrat1500_burnin/"
CSV="$SD/inputs/varalpha_founding.csv"
OUT="${OUT:-$SD/output_varalpha}"
K=6250; NCYC=60; EVERY=60; SS=50; NAQP=1200; NPOLP=520
CONC="${1:-2}"

command -v slim >/dev/null || { echo "ERROR: 'slim' not on PATH (need SLiM 5)."; exit 1; }
[ -f "${FDIR}founders_ch1.vcf" ] || { echo "ERROR: burn-in founders missing ($FDIR)."; exit 1; }
[ -f "$CSV" ] || { echo "ERROR: $CSV missing -- run dev/R/moduleE_varalpha_setup.R first."; exit 1; }
mkdir -p "$OUT/logs"

## job list from the CSV (skip header): "deme N_aq N_pol"
JOBS=$(awk -F, 'NR>1{print $1, $4, $5}' "$CSV")
echo "variable-alpha: $(echo "$JOBS" | wc -l | tr -d ' ') demes  CONC=$CONC  K=$K gen=$NCYC"
echo "founders: $FDIR"; echo "output:   $OUT"

run_one() {
  read -r d naq npol <<< "$1"
  local seed=$(( 42*10000000 + d*10000 + 60 ))          # distinct from replicate/shared-source seeds
  slim -s "$seed" -d "folder='$OUT'" -d "FDIR='$FDIR'" -d "RECDIR='$RECDIR'" \
    -d "N_AQ_POOL=$NAQP" -d "N_POL_POOL=$NPOLP" -d "N_AQ=$naq" -d "N_POL=$npol" \
    -d "K=$K" -d "nCycles=$NCYC" -d "sampleEvery=$EVERY" -d "sampleSize=$SS" -d "rep=$d" \
    "$MODEL" > "$OUT/logs/deme${d}.log" 2>&1
  echo "[$(date +%H:%M:%S)] done deme$d (aq=$naq pol=$npol, exit $?)"
}
export -f run_one; export MODEL RECDIR FDIR OUT K NCYC EVERY SS NAQP NPOLP
echo "$JOBS" | xargs -P "$CONC" -I {} bash -c 'run_one "$@"' _ {}
echo "DONE. VCFs: $(ls "$OUT"/hyb_*.vcf 2>/dev/null | wc -l | tr -d ' ') ; errors: $(grep -lri error "$OUT"/logs/ 2>/dev/null | wc -l | tr -d ' ')"
