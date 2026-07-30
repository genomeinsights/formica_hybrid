#!/bin/bash -l
## =============================================================================
## Module E -- SHARED-SOURCE (complex-demography) neutral null
## =============================================================================
## Tests whether SHARED pre-split ancestry (a common admixed source that later
## splits into the 20 populations) -- rather than 20 fully independent foundings --
## reproduces the empirical F_ST/LD pattern. Anchored to Nouhaud et al.: ~50
## generations of hybridisation TOTAL. The lever is the split time T1:
##
##   T1 =  0 (T2=50) : fully INDEPENDENT demes            (= existing null, 50 gen)
##   T1 = 25 (T2=25) : half shared, half independent drift
##   T1 = 45 (T2= 5) : near-PANMICTIC (90% shared history) *
##
## (*) T2 must stay > 0: with T2=0 a sub-deme is just the founding draw and never
## produces a sampleable brood. T1=45/T2=5 keeps the total at 50 generations while
## post-split drift is negligible relative to 45 shared -- the panmictic end of the
## axis. The literal single-population point is already covered by the earlier
## single-panmictic-population test.
##
## Two-stage per (T1>0, replicate):
##   STAGE 1  one ANCESTRAL deme: found 450/195 from the diversified burn-in pool,
##            evolve T1 gen, then OUTPUT_FOUNDERS dumps POOL_AQ males + POOL_POL
##            queens as a founder VCF (haploid male cols + diploid queen cols).
##   split    that combined VCF into per-chromosome founders_ch<id>.vcf.
##   STAGE 2  20 SUB-DEMES, MODEL UNCHANGED: found 450/195 from the ancestral pool
##            (FDIR -> the split copy), evolve T2 gen, sample 50 females each.
## T1=0 skips stages 1/2 and just founds 20 demes from the burn-in pool (50 gen).
##
## Config: T1LIST (default "0 25 45"), NREP (default 1; 20 for null envelopes),
## CONC (arg 1, default 2), POOL_AQ/POOL_POL (default 1200/520).
## Output: output_sharedsource/T1_<nn>/rep<NN>/hyb_..._rep<deme>_gen<T2>.vcf
## Resources: ~95 s + 4.8 GB per deme (stage 1 and each sub-deme); ~21 SLiM runs
## per (T1>0, replicate).
## =============================================================================
set -u
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODEL="$SD/SpecIAnt_rufa_neutral_realfounders.slim"; RECDIR="$SD/rec_maps/"
BURNIN="$SD/founders/maf015_DIstrat1500_burnin/"
OUT="${OUT:-$SD/output_sharedsource}"
T1LIST="${T1LIST:-0 25 45}"; TOTAL="${TOTAL:-50}"
NREP="${NREP:-1}"; NDEMES=20; K=6250; SS=50
NAQ=450; NPOL=195; POOL_AQ="${POOL_AQ:-1200}"; POOL_POL="${POOL_POL:-520}"
CONC="${1:-2}"

command -v slim >/dev/null || { echo "ERROR: 'slim' not on PATH (need SLiM 5)."; exit 1; }
[ -f "${BURNIN}founders_ch1.vcf" ] || { echo "ERROR: burn-in founders missing ($BURNIN) -- run make_founders_burnin.R first."; exit 1; }
mkdir -p "$OUT"

## split a combined multi-chromosome VCF into per-chromosome founders_ch<id>.vcf.
## Also blanks the SLiM ID/QUAL/INFO columns to "." : outputIndividualsToVCF writes
## INFO=...;PO=-1;... and readHaplosomesFromVCF rejects PO=-1 (NonnegativeIntegerForString).
## The real founder VCFs use INFO="." (and readHaplosomesFromVCF is told the type = m10),
## so stripping these fields makes the ancestral pool read exactly like the burn-in pool.
split_ancestral() {   # $1 = ancestral vcf ; $2 = target dir
  local vcf="$1" dir="$2"; mkdir -p "$dir"
  awk -v dir="$dir" 'BEGIN{OFS="\t"}
    /^##/     { meta = meta $0 ORS; next }
    /^#CHROM/ { chdr = $0; next }
    { $3="."; $6="."; $8="."; f = dir "/founders_" $1 ".vcf";
      if (!(f in seen)) { printf "%s%s\n", meta, chdr > f; seen[f]=1 }
      print >> f }' "$vcf"
}

## run one sub-deme (stage 2) or one independent deme (T1=0), model unchanged
run_deme() {          # $1 = "rep d T2 fdir naqp npolp odir"
  read -r rep d T2 fdir naqp npolp odir <<< "$1"
  mkdir -p "$odir/logs"
  local seed=$(( rep*10000000 + d*10000 + T2 ))
  slim -s "$seed" -d "folder='$odir'" -d "FDIR='$fdir'" -d "RECDIR='$RECDIR'" \
    -d "N_AQ_POOL=$naqp" -d "N_POL_POOL=$npolp" -d "N_AQ=$NAQ" -d "N_POL=$NPOL" \
    -d "K=$K" -d "nCycles=$T2" -d "sampleEvery=$T2" -d "sampleSize=$SS" -d "rep=$d" \
    "$MODEL" > "$odir/logs/deme${d}.log" 2>&1
  echo "[$(date +%H:%M:%S)] done T2=$T2 rep$rep deme$d (exit $?)"
}
export -f run_deme; export MODEL RECDIR NAQ NPOL K SS

for T1 in $T1LIST; do
  T2=$(( TOTAL - T1 )); tag=$(printf '%02d' "$T1")
  echo "=================  T1=$T1  T2=$T2  (total $TOTAL gen)  ================="
  for rep in $(seq 1 "$NREP"); do
    odir="$OUT/T1_${tag}/rep$(printf '%02d' "$rep")"; mkdir -p "$odir/logs"

    if [ "$T1" -eq 0 ]; then
      fdir="$BURNIN"; naqp=1200; npolp=520                 # independent: burn-in pool
    else
      ## ---- STAGE 1: ancestral deme -> founder VCF -----------------------------
      aseed=$(( rep*10000000 + 9990000 + T1 ))
      slim -s "$aseed" -d "folder='$odir'" -d "FDIR='$BURNIN'" -d "RECDIR='$RECDIR'" \
        -d "N_AQ_POOL=1200" -d "N_POL_POOL=520" -d "N_AQ=$NAQ" -d "N_POL=$NPOL" \
        -d "K=$K" -d "nCycles=$T1" -d "sampleEvery=999" -d "sampleSize=$SS" -d "rep=0" \
        -d "OUTPUT_FOUNDERS=1" -d "POOL_AQ=$POOL_AQ" -d "POOL_POL=$POOL_POL" \
        "$MODEL" > "$odir/logs/ancestral.log" 2>&1
      avcf=$(ls "$odir"/ancestral_*_gen${T1}.vcf 2>/dev/null | head -1)
      [ -f "$avcf" ] || { echo "ERROR: ancestral VCF missing for T1=$T1 rep$rep -- see $odir/logs/ancestral.log"; continue; }
      fdir="$odir/ancestral_founders/"; split_ancestral "$avcf" "$fdir"
      naqp="$POOL_AQ"; npolp="$POOL_POL"
      echo "[$(date +%H:%M:%S)] T1=$T1 rep$rep ancestral pool ready ($POOL_AQ/$POOL_POL) -> $fdir"
    fi

    ## ---- STAGE 2: 20 sub-demes from the (ancestral or burn-in) pool ----------
    for d in $(seq 1 $NDEMES); do echo "$rep $d $T2 $fdir $naqp $npolp $odir"; done \
      | xargs -P "$CONC" -I {} bash -c 'run_deme "$@"' _ {}
  done
done
echo "ALL DONE. tree:"; find "$OUT" -name 'hyb_*.vcf' | sed 's#/rep[0-9]*/.*##' | sort | uniq -c
echo "errors: $(grep -lri error "$OUT"/T1_*/rep*/logs/ 2>/dev/null | wc -l | tr -d ' ')"
