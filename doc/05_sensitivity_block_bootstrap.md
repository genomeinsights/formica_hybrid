# Sensitivity analyses

Two sensitivity analyses qualify the descriptive A/B/C results without changing any
point estimate: **(1)** a block bootstrap giving dependence-aware confidence
intervals for the main coefficients (*precision*), and **(2)** a Module C check of
whether its climate → diagnostic enrichment is a BayPass-power / allele-frequency
artifact (*confounding*).

## 1 — Block bootstrap of the main coefficients (precision)

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

## 2 — Module C: is the DI-enrichment a power / allele-frequency artifact? (confounding)

**Question.** BayPass evidence scales with among-population variation in hybrid allele
frequency: a cluster near-fixed genome-wide in the hybrids has little among-population
variance and little power (low `p_sig`); one segregating at intermediate frequency has
more. If those higher-power clusters are also more ancestry-diagnostic, the enrichment
`frac_hi ~ p_sig` could reflect a shared dependence on allele-frequency variation rather
than a climate–diagnostic link. This is a *confounding* check, complementary to the
*precision* check above.

**Code:** `R/moduleC_maf_power_sensitivity.R` (primary config, n_loci ≥ 50, PC1 & PC2;
leaves the finalized `moduleC_rate_based.R` untouched).

**Method.** Add simple per-cluster power / variation covariates to the primary
cluster-level model and watch the `p_sig` (enrichment) coefficient:
`M0: frac_hi ~ p_sig + log(n_loci)` (published) → `+ maf_hyb` (hybrid MAF) → `+ recomb`
→ `+ xtx` (BayPass **XtX**, the direct among-population-differentiation / power measure).
The coefficient is reported per +10 pp significance rate; M0 self-validates against the
published +4.5 / +3.7. The fully-adjusted (M3) coefficient additionally carries a
chromosome block-bootstrap 95% CI.

**Results.** (`p_sig` coefficient = pp high-DI per +10 pp climate rate; 878 clusters)

| | M0 baseline | + maf_hyb | + recomb | + xtx (M3) |
|---|---|---|---|---|
| **PC1** | +4.46 (p=0.0008) | +3.48 (p=0.005) | +3.44 (p=0.005) | **+1.82 (p=0.16)** |
| **PC2** | +3.72 (p=0.013) | +2.46 (p=0.078) | +2.77 (p=0.046) | **+1.69 (p=0.23)** |

M3 chromosome block-bootstrap 95% CI: PC1 **[−3.9, 10.0]**, PC2 **[−2.5, 15.5]** (both
include 0).

Hybrid MAF alone removes ~22% (PC1) / ~34% (PC2) of the effect; recombination adds
little; **XtX roughly halves it** (to +1.8 / +1.7) and removes significance on both
axes. So a large share of the apparent enrichment is generic among-population
differentiation: differentiated clusters carry both more BayPass signal and more
diagnostic content, and the *climate-covariate-specific* association explains little
diagnostic content beyond that differentiation.

**Caveat.** XtX and the covariate Bayes factor come from the same BayPass run and are
not independent (the climate association is a structured projection of differentiation),
so M3 is a strong, arguably conservative control and its residual is a lower bound. But
the simpler MAF-only control (M1) already attenuates the effect and drops PC2 below
significance, so the direction of the conclusion is robust to how hard one controls.

### Figure

**`../Figures/moduleC_maf_power_sensitivity.png`** (also `.pdf`). The `p_sig` coefficient
across the covariate ladder, per PC; blue = naive 95% CI, orange (M3) = chromosome
block-bootstrap 95% CI. The estimate falls from ~+4.5 / +3.7 to ~+1.7 and its interval
crosses 0 once differentiation (XtX) is controlled.

**Together with the block bootstrap**, Module C's static DI-enrichment is best read as
weak and largely a differentiation / power byproduct, not a clean climate–diagnostic
link.

---

*This is the last of the descriptive documents (Stage 0 → A → B → C → sensitivity).
Modules D (intrinsic) and E (neutral null) are separate workstreams — see the
[index](README.md#the-scientific-question).*
