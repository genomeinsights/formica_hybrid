# Sensitivity — block bootstrap of the main coefficients

**Question.** The headline A/B/C coefficients are fitted across very large numbers
of observations, so their model-based *z*-statistics are enormous (|z| up to
~180). Even at the LD-cluster / representative level — where marker
pseudo-replication is already removed — neighbouring units still share ancestry
tracts, local recombination and local diversity, so residual spatial dependence
remains and the model-based standard errors overstate precision. Do the main
effects survive a confidence interval that respects that dependence?

**Headline result.** The **sorting / architecture** effects (DI → direction,
recombination → direction, recombination → magnitude) are **robust**: their
intervals widen 2–4× but still exclude 0 under both block definitions. The
**Module C climate → diagnostic enrichment is not robust**: its interval widens
~5× and includes 0 on both climate axes. Its apparent precision was largely an
artifact of treating spatially-correlated clusters as independent.

This is **one shared sensitivity procedure over the existing coefficients**, not a
new analysis: point estimates are untouched; only the uncertainty is
re-quantified.

### Code & data

| Script (`R/`) | Reads | Writes |
|---|---|---|
| `sensitivity_block_bootstrap.R` | `moduleA_snp.rds`, `hybrids_only_maf005.Rdata` (map only), recmap, `eMLG_5loci_0025_cM05.rds`, the two primary BayPass `.out` files | `data/sensitivity_block_bootstrap.rds`, `Figures/sensitivity_bootstrap_forest.{pdf,png}` |

The sorting frame is `moduleA_snp.rds` subset to the cluster representatives
(`parallelism_stats` is per-marker deterministic, so this equals Module B's
representative-level statistics — no recomputation). The climate frame is the
has_eMLG clusters, primary config (Åland-excluded × withOmega).

---

## Method

**Mechanism.** A non-parametric **block bootstrap**: resample whole genomic blocks
with replacement to their original number, refit the model on the pooled resample,
collect the target coefficient; **B = 10,000** replicates → central-95% percentile
interval. A block drawn *k* times contributes its units *k* times.

**Two block definitions**, both reported:
- **chromosome** — 27 blocks; most conservative, assumption-light about the
  dependence length-scale.
- **~10 cM contiguous blocks** — several hundred blocks; finer, tighter intervals.
  10 cM exceeds the ~10–15 cM admixture-LD decay in these populations.

**Level.** Unit level only, matching the pipeline's unit-of-analysis principle:
LD-pruning representatives for the sorting coefficients (Module B's units;
`du_m` = 252,402 differentiated representatives for magnitude, `cu` = 103,409
unidirectional units for direction), has_eMLG clusters for the climate coefficient
(878 clusters with ≥ 50 loci).

**Self-validation.** Before bootstrapping, each model is refit on the full data and
the script **stops** unless the estimate reproduces its published value — so the
reconstructed frames are known faithful. All five reproduced (1.462, −0.091, 0.052,
4.457, 3.717 vs published 1.46, −0.09, 0.05, 4.5, 3.7).

**Reporting.** Each coefficient is shown on its published scale, with the
model-based, chromosome-block, and 10 cM-block 95% intervals side by side.

---

## Results

| Coefficient (scale) | Estimate | Model-based 95% | **Chromosome 95%** | **10 cM 95%** | Excludes 0? |
|---|---|---|---|---|---|
| DI → direction (log-odds / +1 SD DI) | +1.462 | [1.44, 1.49] | [1.41, 1.52] | [1.41, 1.51] | **yes** |
| recomb → direction (log-odds / +1 SD log-rec) | −0.091 | [−0.115, −0.068] | [−0.145, −0.037] | [−0.137, −0.045] | **yes** |
| recomb → magnitude (Δ prop_fixed / +1 SD log-rec) | +0.052 | [0.051, 0.053] | [0.046, 0.057] | [0.043, 0.058] | **yes** |
| climate rate → diagnostic, **PC1** (pp / +10 pp rate) | +4.46 | [1.87, 7.04] | **[−3.98, 13.06]** | **[−4.16, 12.54]** | **no** |
| climate rate → diagnostic, **PC2** (pp / +10 pp rate) | +3.72 | [0.81, 6.63] | **[−0.71, 17.79]** | **[−0.42, 19.34]** | **no** |

(0 bootstrap failures in all cells; results stable between B = 1,000 and 10,000.)

- **Sorting / architecture effects are robust.** All three widen 2–4× under block
  resampling but still cleanly exclude 0 under *both* block definitions. Notably
  the *weak* recombination → direction lead (−0.091) survives — its
  chromosome-block interval [−0.145, −0.037] still excludes 0 — so the "modest
  low-recombination → aquilonia lead" of Module B is not a large-*N* precision
  artifact (it remains descriptive pending the Module E null).
- **The climate → diagnostic enrichment does not survive.** The model-based
  intervals exclude 0 (PC1 [1.9, 7.0], PC2 [0.8, 6.6]), but the block-bootstrap
  intervals widen ~5× and **include 0 on both axes and both block definitions**.
  They are strongly right-skewed (upper tails to +13 / +19), meaning the apparent
  enrichment is carried by a few chromosomes/blocks. With 878 clusters spread over
  27 chromosomes, the effective information is far below the naive count, and the
  DI-enrichment is not distinguishable from 0 once that dependence is respected.

### Figure

**`../Figures/sensitivity_bootstrap_forest.png`** (also `.pdf`). Per coefficient
(free *x*-axis), the model-based, chromosome-block and 10 cM-block 95% intervals,
with a dashed line at 0. The three sorting/architecture panels sit clearly off 0
under all three methods; the two climate panels cross 0 once block resampling is
applied.

---

## Interpretation

The sensitivity pass sharpens the descriptive picture in two directions. It
**confirms** that the sorting-direction structure — DI governing direction, and the
weak independent recombination lead — is real in the data and not merely a
large-sample precision effect. It **tempers** Module C: the *dynamic*
sorting–climate link was already null, and now the *static* DI-enrichment, its one
positive climate signal, is shown not to survive a dependence-aware interval. The
climate association in these data is, at most, weak and spatially concentrated.

As with the modules themselves, this is about precision, not causal licence: wider
intervals qualify how confidently an effect is distinguished from 0, while the
neutral-null comparison (Module E) remains what any *selection* reading requires.

---

*This is the last of the descriptive documents (Stage 0 → A → B → C → sensitivity).
Modules D (intrinsic) and E (neutral null) are separate workstreams — see the
[index](README.md#the-scientific-question).*
