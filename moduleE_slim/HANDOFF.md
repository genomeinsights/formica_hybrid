# Module E — handoff for a new session / machine

Self-contained bundle for continuing the Module E neutral-null analysis on another
computer. **No external data needed** — the 1 GB empirical dataset is replaced by
the 19.5 MB `inputs/empirical_bundle.rds`. Everything runs from a clone of this
directory. `results/moduleE_results_summary.pdf` is the written-up science; this
file is the operational context.

## Prerequisites
- **SLiM 5** on `PATH` (`slim -v`) — messerlab.org/slim
- **R** with `data.table`, `ggplot2`, `patchwork`
- **Path note:** the *portable* scripts below (`run_ksweep.sh`,
  `scripts/analyze_di_stratified.R`, `scripts/make_founders.R`) resolve all paths
  relative to this bundle and need no editing. The original `dev/R/moduleE_*.R`
  scripts in the parent repo have hardcoded `/Users/petrikem/...` paths and the
  1 GB dataset dependency — prefer the portable versions here.

## The question and the result so far

Twenty replicate *F. polyctena* × *F. aquilonia* hybrid populations show directional
ancestry sorting. Module E is the recombination-matched **haplodiploid neutral
null** that tests whether that exceeds drift. Founders are seeded from **real phased
parental haplotypes**; recombination follows the empirical cM map; all selection is
removed.

**Key result (DI-stratified, the design that worked):** markers are binned by
DiagnosticIndex (DI) over the parent_maf≥0.15 universe. Low-DI markers cannot be
driven by sorting, so they are an **internal demographic control**; the gap vs DI
is a dose-response.

| DI bin | F_ST obs | F_ST null | sorted obs | sorted null |
|---|---|---|---|---|
| ≤ −90 (neutral) | 0.051 | ~0.14 | 0.026 | 0.010 |
| −40…−30 | 0.305 | 0.123 | 0.260 | 0.059 |
| > −15 (diagnostic) | 0.315 | 0.112 | 0.116 | ~0 |

Observed differentiation and sorting **scale with DI**; the neutral null is
**flat**. Drift is DI-blind, so a DI-scaled excess can't be a demographic artifact.
This is the strongest statement so far.

## The immediate next step (why we migrated)

**The null does not yet reproduce the neutral background:** at DI ≤ −90 it gives
F_ST ~0.14 vs observed 0.051 (over-drifts) and π ~0.40 vs 0.448. Both say the
**simulated effective size is too small**. The plan:

1. **Sweep K upward** past 6250 (e.g. 12500, 25000) at full founding (30/13). More
   K → less drift → lower neutral F_ST toward 0.051, higher π toward 0.448.
2. **Founder-number test** at fixed K (e.g. (12 6), (6 3)) — expected to *worsen*
   the neutral fit, constraining founding to be large (rules out a tight-bottleneck
   neutral explanation on the trustworthy markers).
3. **Anchor on the neutral bin** (F_ST + π), NOT on LD-decay (LD favours the
   opposite, over-drifted direction and is unreliable here) and NOT on high-DI
   markers (those are the test).
4. Read the high-DI excess at the neutral-matched cell. Because the null is flat
   across DI, matching neutral F_ST at ~0.05 drops the null's high-DI F_ST to ~0.05
   too, so the excess (obs 0.315) **grows** — the current numbers are conservative.

The migration reason: K=12500+ at 40k markers needs ~15+ GB/deme; the previous
machine was RAM-limited. Use `make_founders.R` to cut markers per bin if needed.

## How to run

```bash
cd moduleE_slim

# (optional) lighter founder set for high K -- 2500/bin = 25k markers, ~half RAM
Rscript scripts/make_founders.R 2500        # -> founders/maf015_DIstrat2500/

# sweep K upward (full founding). Set CONC so CONC*(GB/deme) < RAM.
KLIST="12500 25000" bash run_ksweep.sh 4     # -> output/  (default founders = 4000/bin)
#   to use the lighter set:  FDIR=$PWD/founders/maf015_DIstrat2500/ KLIST="12500 25000" bash run_ksweep.sh 8

# founder-number test at K=6250
FLIST="12 6;6 3" KLIST="6250" bash run_ksweep.sh 6

# analyse (auto-detects K/gens present; anchors on the neutral bin)
Rscript scripts/analyze_di_stratified.R output Naq30_Npol13
#   -> prints neutral-bin fit ranking + the dose-response table; saves output/di_stratified_results.rds
```
`analyze_di_stratified.R` was validated on this machine against the K∈{1500,6250}
runs and reproduces the table above.

## Memory guide
~0.19 GB per 1000 markers per deme at K=6250, scaling ~linearly with K.
40k markers: 7.6 GB (K=6250), ~15 GB (K=12500), ~30 GB (K=25000).
25k markers (2500/bin): 4.7 / 9.5 / 19 GB. Pick CONC accordingly.

## Contents
```
moduleE_slim/
  SpecIAnt_rufa_neutral_realfounders.slim   the neutral model (real founders, no selection)
  run_ksweep.sh                             portable K/founder sweep harness
  rec_maps/                                 27 empirical recombination maps (cM/Mb)
  founders/maf015_DIstrat4000/              40k-marker founder VCFs (4000 per DI bin)
  inputs/empirical_bundle.rds               COMPACT empirical (per-pop freqs + low-DI genotypes + DI + parents)
  inputs/moduleE_founder_haplotypes.rds     phased founder pool (for make_founders.R)
  scripts/analyze_di_stratified.R           portable analysis (reads the bundle)
  scripts/make_founders.R                   portable founder-VCF generator
  scripts/parallelism_stats.R              the sorting statistic (locked conventions)
  scripts/moduleE_*.R                        original (non-portable) analysis scripts, for reference
  results/                                   di_stratified_results.rds, figures, LaTeX + PDF
  HANDOFF.md                                 this file
```

## Model facts and gotchas (hard-won)
- **Founders read as ONE combined VCF per chromosome** (30 haploid aq samples +
  13 diploid pol) via `readHaplosomesFromVCF`; reading species separately makes
  two mutations per site. VCFs must be **uncompressed**.
- **Recombination `scale=1e-8`**: the maps are cM/Mb; `scale=1.0` (in the original
  upstream model) errors because rates exceed 0.5.
- **Filenames encode `_K<K>`** so parameter sets never collide in one output dir.
- **Sample-size matching is mandatory**: sim demes are subsampled to the empirical
  per-population sizes (3–20 diploids) before any statistic; composite r² and
  sort_class are both sample-size sensitive.
- **F_ST estimator** is Nei-style ratio-of-averages; **orientation** and the
  sorting statistic use the locked Module A conventions (parent_maf≥0.15,
  fix_th=0.15, sort_th=0.5, null_prob=0.5), with parents = the phased pool.
- **`empirical_bundle.rds`**: `$emp_mean` is per-pop dosage ×1000 (integer, divide
  by 1000); `$emp_geno_anchor` is individuals × low-DI markers for an optional LD
  anchor; `$universe` has marker/Chr/Pos/cM/DI/parent_maf; `$f_aq_par`/`$f_pol_par`
  are the parental reference. Rebuilt from the 1 GB dataset only if the marker
  universe changes.
- A **queenless-cycle bug** in the founding logic (empty `nmates_all`) is fixed
  with a `size>0` guard.

## Neutral alternatives already excluded (do not re-litigate)
- **Founder heterogeneity** — even the absolute ceiling (each deme from one aq
  haplotype + one pol female) gives at-founding F_ST 0.090 vs a ~0.11 gap; realistic
  sub-pools give 0.01–0.02. At diagnostic markers a 50/50 admixture washes out
  within-species structure.
- **Colony/nest sampling** — measured within-nest excess genotype correlation is
  only 0.047 (polygyny); minor. An earlier full-sib emulation was wrong and retracted.
- **Mixed population ages** — age mixtures interpolate along the LD-vs-drift
  trade-off, no cell matches all summaries.
- **Among-population ancestry variance** — contributes ~0.017 of the high-DI gap.
- **Migration** — not pursued (no evidence for it); the low neutral F_ST is treated
  as a large-Ne / low-drift signal to be matched by the K sweep instead.
```
