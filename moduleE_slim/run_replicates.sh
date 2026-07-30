#!/bin/bash -l
## =============================================================================
## Module E -- WHOLE-EXPERIMENT replication for null envelopes
## =============================================================================
## The 20 simulated demes are ONE realization of the study, not 20 independent
## replicates (assessment section 3). To put a null envelope on experiment-level
## statistics (F_ST-vs-DI slope, adjacent-LD-vs-DI, the F_ST/LD ratio, sorted-locus
## counts), we repeat the COMPLETE 20-deme experiment NREP times: each replicate
## draws its own 450/195 founders (resampled from the diversified burn-in pool via
## a fresh seed) and evolves with independent seeds. Analyse per replicate, then
## take the distribution across replicates as the null envelope.
##
## Founders = the calibrated diversified pool (make_founders_burnin.R; 1200 aq /
## 520 pol), so the neutral LD background matches empirical. K=6250, gen 60.
## Config: NREP (default 20), CONC (arg 1, default 2). ~95 s/deme, 4.8 GB/deme.
## Output: output_replicates/rep<NN>/hyb_..._rep<deme>_gen60.vcf
## =============================================================================
set -u
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODEL="$SD/SpecIAnt_rufa_neutral_realfounders.slim"; RECDIR="$SD/rec_maps/"
FDIR="$SD/founders/maf015_DIstrat1500_burnin/"
OUT="${OUT:-$SD/output_replicates}"
NREP="${NREP:-20}"; NDEMES=20; K=6250; NCYC=60; EVERY=60; SS=50
NAQ=450; NPOL=195; NAQP=1200; NPOLP=520
CONC="${1:-2}"

command -v slim >/dev/null || { echo "ERROR: 'slim' not on PATH (need SLiM 5)."; exit 1; }
[ -f "${FDIR}founders_ch1.vcf" ] || { echo "ERROR: burn-in founders missing ($FDIR) -- run make_founders_burnin.R first."; exit 1; }
mkdir -p "$OUT"

JOBS=$(for rep in $(seq 1 $NREP); do for d in $(seq 1 $NDEMES); do echo "$rep $d"; done; done)
echo "replicates: $NREP x $NDEMES demes = $(echo "$JOBS" | wc -l | tr -d ' ') runs  CONC=$CONC  K=$K gen=$NCYC"
echo "founders: $FDIR"; echo "output:   $OUT"

run_one() {
  read -r rep d <<< "$1"
  local odir="$OUT/rep$(printf '%02d' "$rep")"; mkdir -p "$odir/logs"
  local seed=$(( rep*10000000 + d*10000 + 7 ))
  slim -s "$seed" -d "folder='$odir'" -d "FDIR='$FDIR'" -d "RECDIR='$RECDIR'" \
    -d "N_AQ_POOL=$NAQP" -d "N_POL_POOL=$NPOLP" -d "N_AQ=$NAQ" -d "N_POL=$NPOL" \
    -d "K=$K" -d "nCycles=$NCYC" -d "sampleEvery=$EVERY" -d "sampleSize=$SS" -d "rep=$d" \
    "$MODEL" > "$odir/logs/deme${d}.log" 2>&1
  echo "[$(date +%H:%M:%S)] done rep$rep deme$d (exit $?)"
}
export -f run_one; export MODEL RECDIR FDIR OUT K NCYC EVERY SS NAQ NPOL NAQP NPOLP
echo "$JOBS" | xargs -P "$CONC" -I {} bash -c 'run_one "$@"' _ {}
echo "ALL DONE. replicate dirs: $(ls -d "$OUT"/rep* 2>/dev/null | wc -l | tr -d ' ') ; errors: $(grep -lri error "$OUT"/rep*/logs/ 2>/dev/null | wc -l | tr -d ' ')"
