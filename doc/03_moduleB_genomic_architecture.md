# Module B — genomic architecture of sorting

**Question.** How does ancestry sorting relate to the genomic landscape —
recombination rate, within-species diversity (π), absolute divergence (d_xy),
relative differentiation (F_ST), and cluster size? In particular, is the
diagnostic-index polarity from Module A just a proxy for recombination or block
structure?

**Headline result.** The diagnostic index (DI) is **orthogonal to recombination**
and reflects **reduced within-species diversity**, not elevated divergence. Sorting
is **depleted, not enriched, in low-recombination regions** once pseudo-replication
is removed. Direction is governed by DI; recombination adds only a weak independent
pull toward *aquilonia*.

### Code & data

| Script (`R/`) | Reads | Writes |
|---|---|---|
| `moduleB_architecture.R` | `hybrids_and_parents_maf005.Rdata`, `hybrids_only_maf005.Rdata`, `eMLG_5loci_0025_cM05.rds`, recmap; sources `Ohta`, `parallelism_stats` | `moduleB_architecture.rds`, `moduleB_fig2.{pdf,png}` |
| `moduleB_eMLG_vs_rep.R` (validation) | `moduleA_clusters.rds`, `moduleA_snp.rds`, `eMLG_5loci_0025_cM05.rds`, `eMLG_sorted_cM05.rds`, `hybrids_only_maf005.Rdata` | `moduleB_eMLG_vs_rep.rds`, `eMLG_vs_rep_cor.png` |

Run 2026-07-21. Structure: **B1** architecture (DI vs covariates) → **B2**
sorting-vs-recombination at unit vs SNP level → **B3** direction vs architecture.

---

## Methods

### Genomic covariates

- **Recombination rate** is assigned per marker by linear interpolation of the
  linkage map (cM/Mb).
- From parental allele frequencies `p_aqu`, `p_pol`:
  - **within-species diversity** = mean of the two species' expected
    heterozygosities;
  - **absolute divergence** `d_xy = p_aqu(1 − p_pol) + p_pol(1 − p_aqu)`;
  - **relative differentiation** `F_ST = (H_T − H_S)/H_T`.

Reporting relative *and* absolute measures alongside diversity is what lets a
differentiation peak caused by *reduced within-species diversity* be distinguished
from one caused by *genuinely elevated divergence* — a distinction a relative measure
alone cannot make.

### Marker-level vs unit-level summaries

Sorting is summarised across recombination deciles **twice**: once treating each
independent unit as one observation (pruning representatives), and once on a random
sample of individual markers. This comparison is **diagnostic**: adjacent markers in a
low-recombination region report the same underlying event, so a relationship present
at marker level but absent at unit level reflects the spatial redundancy of markers,
not a property of the units.

### Cluster size — collider in one question, confounder in another

The same variable is deliberately treated in opposite ways, so this is stated
explicitly:

- **Relating recombination to sorting**, cluster size is a **collider** (a consequence
  of both recombination and local diversity). It is **excluded** from those models:
  conditioning on it would induce an apparent recombination–sorting association that is
  not otherwise there.
- **Relating climate-outlier status to DI** (Module C), cluster size is a
  **confounder** (it independently predicts DI genome-wide and is required by any
  member-count criterion). It is **adjusted for** there.

### Direction of sorting vs architecture

To ask whether architecture predicts the *direction* of sorting beyond DI, the
direction of unidirectionally sorted units is modelled by **logistic regression on
standardised DI, log recombination rate, pooled parental MAF, and log cluster size**,
fitted across units. Any recombination term here is an **association**: under
neutrality a non-recombining autosomal region does *not* experience stronger drift
than its surroundings, so attributing such a term to linked selection requires the
explicit demographic null (Module E).

---

## Results

All findings are **descriptive**, pending the Module E null.

**DI is separable from recombination and reflects reduced diversity.** Across the
genome DI is essentially orthogonal to map-interpolated recombination rate (Spearman
**ρ = −0.03**) and to pooled parental MAF (**ρ = −0.02**), making it a genuinely
independent predictor of sorting. Its dominant correlate is **low within-species
diversity** (**ρ = −0.46 with π**), not elevated absolute divergence (**ρ = +0.13 with
d_xy**) — so DI is largely a *reduced-diversity* differentiation signal, not a
*divergence-island* one. A classical differentiation-island signature (jointly
elevated F_ST and d_xy, reduced π, very large linkage blocks) is confined to the
**lowest-recombination decile**, yet DI is **flat across recombination** and does not
track it. DI and recombination are therefore separable axes.

**Sorting is depleted, not enriched, in low-recombination regions.** At the
independent-unit level, the fraction of units that sort is **lowest in the
least-recombining decile**:

> frac_sorted ≈ **0.06** in the lowest decile vs **~0.42** elsewhere.

The higher SNP-level value in that decile (**0.17**) relative to the unit level (0.06)
is **spatial pseudo-replication** — adjacent, redundant SNPs in a single
low-recombination block. Consistently, sorting *magnitude* **increases** with
recombination (`prop_fixed ~ recombination`, standardised slope **+0.05, z = 117**).
Sorting is thus **not** concentrated where recombination is low — contrary to a naive
linked-selection expectation, and consistent with neutral book-keeping (fewer
independent units per low-recombination block) pending the null.

**Direction is governed by DI; recombination adds only a weak independent pull.** In
the logistic direction model among unidirectional units:

> `P(aquilonia) ~ DI + recombination + parental MAF + log(cluster size)`:
> **DI +1.46, z = 125** (high-DI → *aquilonia*, low-DI → *polyctena*, reproducing the
> Module A polarity); **recombination −0.09, z = −7.6**.

The recombination effect is small, significant, and in the *opposite-to-noise*
direction: **low recombination weakly leans *aquilonia*** after DI, parental MAF and
size are controlled — a modest linked-selection lead, **~16× weaker than DI**. At the
raw unit level, direction is nearly flat across recombination; the stronger
low-recombination *aquilonia* excess at the SNP level is largely pseudo-replication
that flattens once independent units are counted.

**Validation — unit choice does not change direction.** The `moduleB_eMLG_vs_rep.R`
check confirms that sorting **direction** is robust to whether a cluster is
represented by its consensus (eMLG) or its representative marker (**~99% agreement**).
The eMLG's added value is the **magnitude** call for ≥ 5-marker clusters and the
LD-based tests (Module D); representatives are adequate for Module B's direction
analysis.

### Figures

- **Fig. 2 — `../Figures/moduleB_fig2.png`** (also `.pdf`). Three panels:
  - **(a)** Standardised F_ST, d_xy and π across recombination deciles — the
    differentiation-island signature (high F_ST + d_xy, low π) is confined to the
    lowest-recombination decile.
  - **(b)** Fraction of loci sorted vs recombination at the SNP level *and* the
    independent-unit level — the unit-level collapse in low recombination (0.17 → 0.06)
    is the removal of spatial pseudo-replication.
  - **(c)** Direction of unidirectional sorting (fraction toward *aquilonia*) vs
    recombination at the unit level — the DI polarity is **not** mirrored by a
    recombination gradient.
- **`../Figures/eMLG_vs_rep_cor.png`** — correlation of the sorting call between eMLG
  consensus and representative marker (~99% direction agreement).

---

## Interpretation and caveats

Module B establishes that the predictable polarity of sorting is organised by locus
**diagnosticity** — a reduced-diversity axis orthogonal to recombination — and that
sorting is **depleted rather than enriched** in low-recombination regions once
pseudo-replication is removed. A **weak residual low-recombination pull toward
*aquilonia*** survives control for DI: a lead for the intrinsic
(linked-selection / incompatibility) hypothesis. Adjudicating it requires the
recombination-matched, haplodiploid neutral null (Module E); the intrinsic and
extrinsic drivers are then discriminated, respectively, by among-population LD between
unlinked diagnostic units (Module D) and by overlap with climate genotype–environment
outliers (Module C, next).

Both recombination effects are **robust to a dependence-aware confidence interval**: under
the shared block bootstrap ([doc 05](05_sensitivity_block_bootstrap.md)) the magnitude
slope (+0.052) and — notably — the weak direction lead (−0.091) both still exclude 0
(chromosome-block 95% [0.046, 0.057] and [−0.145, −0.037]). So the low-recombination
→ *aquilonia* lead is real in the data, not a large-*N* precision artifact (still
descriptive pending the Module E null).

---

*Next: [Module C — association with climate](04_moduleC_climate_association.md).*
