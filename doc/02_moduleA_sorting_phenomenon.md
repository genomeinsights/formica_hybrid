# Module A — the sorting phenomenon

**Question.** Across twenty independently formed hybrid populations, does parental
ancestry sort **predictably** — the same parental allele approaching fixation
repeatedly — or in directions that vary among populations? And what governs the
*direction* of any predictable sorting?

**Headline result.** Sorting is largely **unidirectional** rather than
bidirectional, and its **direction is governed by the diagnostic index (DI)**:
low-DI loci fix toward *polyctena*, high-DI loci toward *aquilonia*. The signal
resides in many small, largely unlinked units, not in large linkage blocks.

### Code & data

| Script (`R/`) | Reads | Writes |
|---|---|---|
| `moduleA_sorting_phenomenon.R` | `hybrids_and_parents_maf005.Rdata`, `hybrids_only_maf005.Rdata`, `eMLG_5loci_0025_cM05.rds`; sources `dev/R/{Ohta,parallelism_stats,eMLG_parallelism}.R` | `moduleA_snp.rds`, `moduleA_clusters.rds`, `moduleA_dilution.rds`, `eMLG_sorted_cM05.rds` |
| `moduleA_di_asymmetry.R` | `moduleA_snp.rds`, `moduleA_clusters.rds` | `moduleA_di_asymmetry.rds`, `moduleA_fig1.{pdf,png}` |

Run 2026-07-21, ~34 min, 8 cores. Structure: **A1** per-SNP sorting statistic →
**A2** build the sorted-eMLG companion (`eMLG_sorted_cM05.rds`) → **A3** one
cluster-level sorting statistic.

---

## The sorting statistic

The core statistic is implemented in `dev/R/parallelism_stats.R` and applied
identically to individual markers (A1) and to cluster consensus genotypes (A3).

### Orientation and near-fixation

For each locus, dosage is oriented into an ***aquilonia*-allele frequency** using
the parental reference samples: the sign of the between-species allele-frequency
difference defines which allele is counted, and loci where the parents do not
differ have no defined orientation. Each hybrid population is then classified as

- **near-fixed toward *aquilonia*** if its oriented frequency ≥ 1 − θ,
- **near-fixed toward *polyctena*** if ≤ θ,
- **unsorted** otherwise,

with tolerance **θ = 0.15**. θ is held constant across all analyses (not tuned per
analysis) because it interacts with population sample size — a population of few
diploid individuals can only realise a coarse grid of frequencies.

### Decomposition and classification

Let `n_obs` be the number of populations with data and a defined orientation, and
`n_aqu`, `n_pol` the numbers near-fixed in each direction. Sorting decomposes as

$$ \underbrace{\frac{n_{aqu}+n_{pol}}{n_{obs}}}_{\text{prop\_fixed}}
 = \Bigl|\underbrace{\tfrac{n_{aqu}-n_{pol}}{n_{obs}}}_{\text{uni}}\Bigr|
 + \underbrace{\frac{2\min(n_{aqu},n_{pol})}{n_{obs}}}_{\text{bi}}, $$

separating populations fixing in a **consistent** direction (`uni`) from those
fixing in **opposing** directions (`bi`). A single threshold **τ = 0.5** classifies
each locus:

- **unidirectional** when `|uni| ≥ bi` and `|uni| ≥ τ`,
- **bidirectional** when `bi > |uni|` and `bi ≥ τ`,
- **unsorted** otherwise.

Because all three quantities are fractions of the population panel, the same
threshold is directly comparable between individual markers and aggregated clusters,
and genome-wide summaries are simple tallies over units.

### Test of directional predictability

Among the `n_fixed` populations that near-fixed in either direction, predictability
is tested as

$$ n_{aqu} \sim \text{Binomial}(n_{fixed},\, 0.5), $$

two-sided and exact, computed only where **`n_fixed ≥ 5`**, with **Benjamini–Hochberg**
adjustment across tested loci. The null probability **0.5** encodes the symmetric
expectation (populations fix in either direction with equal chance) and is also the
empirically correct one: genome-wide mean *aquilonia* ancestry at diagnostic loci is
≈ 0.5.

A per-locus alternative — testing each locus against its own mean ancestry across
populations — is **not** used: that expectation is a direct function of the outcome
being tested and would absorb the very signal of interest (the "pooled" null is
~0.96-correlated with the outcome and yields zero significant loci).

**The classification (`sort_class`) and the test (`sig` from the BH-adjusted
binomial q-value) are kept as separate quantities** and never combined into one
label. The classification says how much and in which direction a unit sorted at a
fixed threshold; the test says whether the direction departs from the
random-direction expectation.

### The polymorphism gate

Near-fixation is evidence of *sorting* only if the locus was segregating in the
founding admixture. Analyses are therefore restricted to loci **polymorphic in the
parents**, using the **pooled parental minor allele frequency** (from pooled allele
counts across all parental individuals, folded to the minor allele) at a threshold
of **0.15**.

The gate is on parental allele frequency, **not on DI** — deliberately. DI is a
predictor to be tested later (against recombination, cluster size, climate
association, and the direction of sorting). Selecting loci on DI would truncate its
range, attenuate every one of those relationships, and make "does DI predict
sorting?" circular. Parental allele frequency instead removes loci already
near-monomorphic in the founding pool (where hybrid near-fixation reflects the
starting condition, not a sorting event) while leaving DI free to vary fully.

*(Label symmetry: the sorting magnitude and the unsorted/bi/uni classification are
invariant to which allele is called "aquilonia"; only the parental **direction**
depends on orientation, and is interpretable only where the parents differ.)*

### Testing at the cluster level

The same statistic is applied to the eMLG **consensus** genotypes so each LD cluster
contributes one test, not one per marker. Stored consensus genotypes cover clusters
of ≥ 5 markers; for the remaining clusters that contain a sorted marker, consensus
genotypes are built by the identical procedure, so every cluster is tested at
whatever resolution its size permits.

Two subtleties, both handled explicitly:

- **Parental consensus reuses the hybrid choices.** Testing a cluster needs parental
  consensus genotypes (orientation is defined from the parents). These are built
  using the reference-marker choice and per-marker flip decisions **derived from the
  hybrid genotypes**, not recomputed on the parents alone — recomputing on a small,
  strongly differentiated parental sample could pick a different reference marker and
  flip a different subset, silently putting hybrid and parental consensus on opposite
  orientations.
- **Cluster DI is the maximum across members** (`di_agg = "max"`), so a cluster is
  excluded by a DI criterion only if *none* of its members is ancestry-informative.
  (Using the representative marker alone spuriously excluded clusters — 60/315 in a
  test — whose informative markers lie elsewhere in the cluster.)

### Absence of signal vs dilution

A cluster may contain a marker that sorts individually yet not sort once aggregated.
Two situations produce this and are separated using the same round-trip fidelity
measure that governs merging: a **coherent** cluster whose consensus is faithful and
genuinely does not sort, versus a **heterogeneous** cluster in which averaging across
dissimilar members dilutes a real signal. Clusters below the fidelity floor are
flagged for separate consideration (`moduleA_dilution.rds`).

---

## Results

All findings are **descriptive** — each is compatible with neutral drift as well as
with selection (see the [index](README.md#the-scientific-question)).

**Replicated sorting is testable across 20 replicates.** Under drift, or an
equal-fitness incompatibility, populations near-fix ancestry in *random* directions;
under directional selection, or an unequal-fitness Dobzhansky–Muller incompatibility,
the *same* parental allele fixes repeatedly. This is the per-locus, twenty-replicate
generalisation of the three-population parallelism of Nouhaud et al. (2022).

**Most differentiated loci do not sort; those that do are directional, not
bidirectional.** Restricting to parent-polymorphic loci (pooled-parental MAF ≥ 0.15;
**~660,000 loci** at the SNP level):

| Class | Fraction of gated loci |
|---|---|
| unsorted | **67.4%** |
| unidirectional → *polyctena* | **17.3%** |
| unidirectional → *aquilonia* | **14.3%** |
| bidirectional | **0.1%** |

Bidirectional near-fixation — the signature expected under drift or an equal-fitness
incompatibility — is **vanishingly rare**. The consistency of direction across
replicates, not its mere magnitude, is the informative signal. (Overall the SNP-level
tally leans *polyctena* — a flip from earlier DI-gated analyses that leaned
*aquilonia*; the next result reconciles this.)

**The polarity of sorting reverses along the DI axis.** Among unidirectionally sorted
loci, the fraction fixing toward *aquilonia* rises **monotonically from ~0.15 in the
least-diagnostic decile to ~0.74 in the most-diagnostic decile**, crossing parity near
the middle of the DI range (Spearman **ρ = 0.31**; crossover near DI ≈ −50). In a
logistic model,

> `P(aquilonia | unidirectional) ~ DI + parental MAF`:
> DI coefficient **+0.082, z = 181**; parental MAF **−1.92, z = −28**.

DI dominates; higher parental MAF leans *polyctena*. Thus **aquilonia ancestry fixes
preferentially at the most diagnostic loci and polyctena ancestry at the least
diagnostic** — restricting to strongly diagnostic loci (DI above the median) recovers
a ~70% *aquilonia* majority, reconciling the apparent excess of either parent as two
ends of a single diagnostic axis.

**The signal resides in many small units, not in large linkage blocks.** Collapsing
LD-correlated markers into independent units and re-testing each once, sorting is
concentrated in numerous small units:

> **74.7%** of large clusters (≥ 5 markers) show **no** aggregate sorting once
> summarised, versus only **5.0%** of the smallest units.

This collapse is **genuine, not an averaging artifact**: the washed-out clusters
retain high internal fidelity (median `score_eMLG` = **0.856**; only **7 of 16,597**
multi-marker unsorted clusters fall below the 0.80 fidelity floor). The per-marker
impression of pervasive sorting in large blocks is therefore **spatial
pseudo-replication**; the independent-unit signal is carried by many small, largely
unlinked loci. (127,241 clusters passed the gate and were tested at the cluster level.)

**Cluster size predicts direction only through DI.** Larger clusters lean *aquilonia*
and are more diagnostic, so direction *appears* to track genomic architecture — but
this is entirely mediated by DI. Controlling for DI, the size term **reverses sign**:

> `P(aquilonia) ~ log₂(cluster size) + DI + parental MAF`:
> log₂(size) **−0.20, z = −22**; DI **+0.083, z = 120**.

Diagnostic index, not linkage-block size or recombination architecture, is the
governing axis of directional sorting.

### Figure

- **Fig. 1 — `../Figures/moduleA_fig1.png`** (also `.pdf`). Three panels:
  - **(a)** Genome-wide sort-class proportions among parent-polymorphic loci
    (unsorted 67.4%; polyctena 17.3%; aquilonia 14.3%; bidirectional 0.1%).
  - **(b)** Fraction of unidirectionally sorted loci fixing toward *aquilonia* across
    DI deciles (rises ~0.15 → ~0.74; crosses parity mid-range).
  - **(c)** Proportion of independent units that are unsorted as a function of
    linkage-block size (5.0% for the smallest → 74.7% for ≥ 5-marker clusters).

---

## Interpretation and caveats

Module A establishes a **predictable, directionally structured** pattern of ancestry
sorting whose polarity is organised by locus diagnosticity rather than by genomic
architecture. These patterns are **descriptive**: near-fixation in a consistent
direction is expected under directional selection *or* an unequal-fitness
incompatibility, but a formal neutral expectation (Module E) is required before either
can be inferred. The two mechanistic hypotheses the null cannot itself distinguish are
pursued separately — **intrinsic** incompatibilities via among-population LD between
unlinked diagnostic units (Module D), and **extrinsic**, climate-driven selection via
overlap with genotype–environment outliers (Module C).

The governing DI → direction effect is **robust to a dependence-aware confidence
interval**: under the shared block bootstrap ([doc 05](05_sensitivity_block_bootstrap.md)),
the unit-level DI coefficient (+1.46 per SD) has a chromosome-block 95% interval of
[1.41, 1.52] — far from 0 despite a 2–4× widening — so the diagnosticity axis is not an
artifact of the large marker count.

---

*Next: [Module B — genomic architecture](03_moduleB_genomic_architecture.md).*
