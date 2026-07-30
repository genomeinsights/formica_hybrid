# Shared-source (complex-demography) neutral null — design

**Question.** Our current neutral null treats the 20 hybrid populations as *fully
independent* foundings. Reviewers (simulation-design assessment) note this is only one
neutral demography. The natural richer comparator is a **common admixed source that
later splits** into the 20 populations, so the demes share ancestry instead of being
independent. Does that shared ancestry reproduce the empirical F_ST/LD dissociation
(high diagnostic F_ST, low diagnostic adjacent-LD) that independent foundings overshoot
on F_ST?

**Anchor (Nouhaud et al.).** ~50 generations of hybridisation **total**; per-population
Ne ≈ K = 6250. Total generations are held at 50; the free lever is the **split time T1**.

## Model
Two-stage, conditional-null spirit preserved throughout (real phased founders, μ=0,
empirical cM map, haplodiploid life cycle):

- **Stage 1 — ancestral pool.** Found one admixed deme (450 aq + 195 pol from the
  diversified burn-in pool), evolve **T1** generations at K=6250. Then dump
  `POOL_AQ` haploid males + `POOL_POL` diploid queens (default 1200/520) as a founder
  VCF — the `OUTPUT_FOUNDERS` block, males as haploid columns, queens as diploid.
- **Split.** Split that combined VCF into per-chromosome `founders_ch<id>.vcf`
  (blanking the SLiM `ID/QUAL/INFO` columns so it reads exactly like the burn-in pool).
- **Stage 2 — 20 sub-demes, model UNCHANGED.** Each sub-deme founds 450/195 from the
  ancestral pool (`FDIR` → the split copy) and drifts **T2 = 50 − T1** generations,
  sampled to 50 females. The 20 demes now share all pre-split ancestry; their
  among-deme differentiation is only post-split drift + founding draw.

`T1 = 0` skips stages 1–2 and founds 20 demes directly from the burn-in pool
(= the existing independent null, run at 50 gen).

## Sweep (all total 50 generations)
| T1 | T2 | interpretation |
|----|----|----------------|
| 0  | 50 | fully **independent** demes (existing null) |
| 25 | 25 | half shared, half independent drift |
| 45 | 5  | near-**panmictic** (90% of history shared)* |

\* T2 must stay > 0: with T2 = 0 a sub-deme is only the founding draw and never
produces a sampleable brood. T1 = 45 / T2 = 5 keeps the total at 50 while post-split
drift is negligible vs 45 shared — the panmictic end of the axis. (The literal
single-population point is already covered by the earlier single-panmictic test.)

## Expectation / what it decides
As T1 rises, shared ancestry drives neutral (low-DI) F_ST → 0 — which *matches*
empirical and cures the independent null's F_ST overshoot. The decisive question is
whether the **diagnostic** F_ST/LD dissociation emerges too, or stays flat. Prior: our
ancestry decomposition shows a single genome-wide ancestry axis explains only ~6% of
the diagnostic differentiation, so shared ancestry should fix the neutral background but
**not** reproduce the locus-specific diagnostic excess. If so, the "needle doesn't move"
and the non-neutral reading strengthens; if it does move, we've found a neutral
explanation to take seriously. The simulation decides it.

## Files
- `SpecIAnt_rufa_neutral_realfounders.slim` — one additive, guarded block
  (`OUTPUT_FOUNDERS`) + a robustness guard so `sampleEvery > nCycles` yields no
  sampling. Existing runs are byte-for-byte unaffected (block is inert without
  `OUTPUT_FOUNDERS`).
- `run_sharedsource.sh` — the two-stage harness. Config: `T1LIST` (default "0 25 45"),
  `NREP` (default 1; 20 for null envelopes), `CONC` (arg 1), `POOL_AQ/POOL_POL`.
  Output: `output_sharedsource/T1_<nn>/rep<NN>/hyb_..._rep<deme>_gen<T2>.vcf`.
- Analysis reuses `dev/R/moduleE_replicate_envelopes.R` logic, run per T1, overlaying
  the three split times against the empirical F_ST-DI / adjacent-LD-DI / ratio curves.

## Validation (smoke test, K=1500, T1=T2=8, pool 300/150)
Full stage1 → split → stage2 round trip runs clean and produces a 15,000-marker ×
50-sample deme (all diagnostic markers preserved). Mixed haploid-male / diploid-queen
VCF round-trips exactly through `outputIndividualsToVCF` → `readHaplosomesFromVCF`.
Pool 1200/520 availability at K=6250 is confirmed on the first real stage-1 run (the
model errors cleanly and the harness reports it if a pool is short).
