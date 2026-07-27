# Stage 0 — LD decay estimation and LD-based complexity reduction

**Purpose.** Turn ~1.1 M ancestry-polarised markers into a set of statistically
independent *LD clusters*, each summarised by one representative marker and (where
large enough) one consensus genotype. Every downstream module operates on these
units. This stage also produces a per-marker measure of local LD support (`ld_w`)
that is reused as a low-recombination indicator.

### Code & data

| Script (`R/`) | Reads | Writes |
|---|---|---|
| `LD_decay_from_DIEM.R` | `Formica_hybrids_filtered_diem_output.bed.gz`, `Sample_covariate_info_outlier_analysis_20.txt`, `Frufa_DTOL_PR.ref_genome.recmap` | `diem_parsed.rds`, `hybrids_only_maf005.Rdata`, `hybrids_and_parents_maf005.Rdata`, `ld_decay_DIEM_100w.rds`, `ld_tracks_ldw_persnp.rds`, `ld_tracks_a_windows.rds`, `p_roc_low_recombination.png` |
| `ld_pruning_DIEM.R` | `hybrids_only_maf005.Rdata`, `Frufa_DTOL_PR.ref_genome.recmap` | `eMLG_5loci_0025_cM05.rds` (canonical), `pruned_markers.rds`, `eMLG_groups.rds` |

The heavy lifting (LD-decay fitting, clustering, consensus construction) is
implemented in the **`LDscnR`** package (`~/gitlab/LDscnR`); these two scripts are
the project-level drivers.

---

## 1. Data and ancestry polarisation

Genotypes and per-marker ancestry information come from the **DIEM** output for
this dataset. Sites are restricted to those called biallelic. For each marker DIEM
supplies a *polarity* and a *diagnostic index* (DI); the DI measures how
informative the marker is about parental ancestry. Genotype dosages at markers
whose polarity indicated the opposite orientation were replaced by
`2 − dosage`, so that after this step dosage is on a single, consistent
ancestry-oriented scale genome-wide. All later use of "DI" refers to this
per-marker diagnostic index.

Minor allele frequencies were computed **among the hybrid individuals only** and
markers retained at **MAF ≥ 0.05**, giving **1,114,340 markers**. Two genotype
matrices were built over this same marker set — hybrids-only (164 individuals; all
LD/clustering) and hybrids+parents (parental allele frequencies only) — for the
reason given in the [index](README.md#key-shared-objects).

The study comprises **20 hybrid populations, 164 hybrid individuals total**, plus
allopatric reference individuals of both parental species. *Formica* wood ants are
haplodiploid (males haploid); this matters for the demographic modelling (Module E)
and for interpreting within-population variation.

---

## 2. LD decay

### Model and fitting

LD decay is modelled per chromosome as

$$ r^2(d) = b + \frac{c - b}{1 + a\,d}, $$

where `d` is physical distance (bp), `a` is the decay rate, `b` is the long-range
background `r²`, and `c` is short-range `r²`. Background LD is estimated from
**inter-chromosomal** marker pairs, which are unlinked by construction and so give
the appropriate baseline. Within each chromosome the curve is fitted across **20
sliding windows with 50% overlap**, sampling **≤ 5,000 marker pairs per window**
stratified across **20 distance classes** so short and long distances are
represented comparably; `r²` is summarised at its **0.95 quantile**, and the
per-window fits are aggregated over the **central 95%** of windows to limit the
influence of outlying windows.

**Allele-frequency filter during fitting.** Only pairs in which *both* markers
exceed **MAF 0.1** contribute to the background and decay fits, because pairs
involving rare alleles mechanically depress `r²` regardless of true linkage and
would bias `a`, `b`, `c`. This applies to curve fitting only — the stored pairwise
LD edge lists retain all markers, so no marker is lost from clustering.

### From a relative threshold to a physical distance

Decay rates vary systematically with chromosome length, so a robust regression of
`a` on chromosome size is fitted genome-wide. For a target relative LD threshold
`ρ`, the corresponding physical distance follows from the fitted curve as

$$ d = \frac{\rho}{a\,(1-\rho)}. $$

This is how every distance-based window below is obtained. Anchoring distances in
each chromosome's *own* fitted curve — rather than a fixed bp window or a fixed
`r²` cutoff — is what makes a single threshold behave comparably across chromosomes
with very different background LD.

### Local LD support (`ld_w`)

For each marker, `ld_w` is the **median `r²` between that marker and its neighbours
within the window corresponding to ρ = 0.95**. It summarises how strongly a marker
is embedded in local LD. Its practical value: it is a *per-marker* quantity read
directly from the LD data — unlike a recombination rate (which must be interpolated
from a genetic map) and unlike the decay rate `a` (fitted per window), `ld_w` needs
no interpolation and is defined at single-marker resolution.

Because reduced recombination elevates local LD, `ld_w` reflects the local
recombination environment. Taking the least-recombining regions (lowest quantiles
of map-based recombination rate) as the target, both `ld_w` and the window-level
decay rate `a` discriminate them well:

> Across the lowest **5%, 10% and 25%** of windows, `a` achieves **AUC 0.91–0.98**
> and `ld_w` **AUC 0.73–0.86** (Fig. S0.1).

The decay rate is the stronger predictor, as expected of a window-level summary,
but `ld_w` keeps the resolution and directness needed for a per-marker decision.
Such a directly-measured statistic is also useful where a pedigree-based
recombination estimate is not: in this cross roughly **two-thirds of
adjacent-marker intervals contain no observed recombinant**, so the map-based
estimate is uninformative over most of its range. `ld_w` is used in Stage 2 below
to decide which clusters to reconsider.

**Directions matter.** `a` is a decay *rate*, so `a` and recombination rate are
positive correlates of each other, whereas `ld_w` is their *inverse*: `ld_w` spikes
exactly where both `a` and recombination are low (Fig. S0.2).

### Figures

- **Fig. S0.1 — `../Figures/p_roc_low_recombination.png`.** ROC curves for
  discriminating low-recombination regions using `ld_w` and the decay rate `a`
  (validated against map-based recombination rate), at three definitions of "low
  recombination" (lowest 5%, 10%, 25% of windows). Both are effective; `a` is
  stronger (AUC 0.91–0.98), `ld_w` weaker (0.73–0.86) but per-marker and
  interpolation-free.
- **Fig. S0.2 — `../Figures/ld_tracks_chr26_chr10.png`.** Per-marker `ld_w` (grey)
  along two example chromosomes (Chr26, Chr10), each with a pronounced
  low-recombination region, plotted against the window-level decay rate `a`
  (orange, left column) and the map-based recombination rate (blue, right column);
  all series min–max scaled. `ld_w` is elevated across each low-recombination
  region, where both `a` and recombination are reduced.

---

## 3. Complexity reduction: one partition, two projections

LD pruning and consensus-genotype construction are usually treated as separate
procedures, but they answer the same question: **what is the finest partition of
markers that can no longer be distinguished statistically, and therefore
constitutes an independent unit in these data?** Once that partition is fixed, two
projections are useful, and both are produced from a **single clustering pass** so
they always refer to the same underlying units:

- **Pruning representative** — one representative marker per cluster (an LD-pruned
  marker set).
- **eMLG consensus** — one consensus dosage per cluster (an *expected multi-locus
  genotype*), for clusters of ≥ 5 markers.

### What a cluster is — and why it is not a haplotype block

A biallelic marker carries exactly **one bipartition** of the sampled chromosomes.
Where the local genealogy has more than two segregating lineages, mutations on
*different branches* define *different* bipartitions. Markers on the same branch are
perfectly correlated (absent recombination, recurrent mutation, missing data);
markers on different branches are only partially correlated, according to how the
branches relate.

Clustering on LD therefore groups markers **by the partition of individuals they
report, not by position**, with two consequences:

1. **Clusters need not partition the chromosome.** Where several branches segregate
   in one region, the markers tagging each are physically interleaved, so the
   resulting clusters *overlap along the sequence*. The distance restriction keeps
   clusters local, but locality is not spatial exclusivity.
2. Because `r²` is frequency-sensitive, two markers tagging **nested** branches can
   be completely associated on a frequency-independent measure yet show low `r²`;
   clustering on `r²` accordingly *separates* nested branches — exactly the
   situation where one interval contains several overlapping clusters.

We therefore avoid "haplotype block" (a contiguous stretch of limited haplotype
diversity with defined boundaries that partitions the sequence). The units here are
**sets of markers reporting the same partition of individuals**, several of which
may occupy the same interval. This is also what makes a consensus genotype
meaningful: averaging markers that report the same partition reinforces that
partition's dosage, whereas averaging across partitions would blur distinct splits.
The Stage-1 complete-linkage requirement and the Stage-2 correlation floor both
exist to prevent that.

### Stage 1 — local LD clusters

The clustering threshold is derived, **per chromosome, from that chromosome's own
fitted decay curve at ρ = 0.5** (its half-decay point), not from a fixed `r²`
cutoff. It is deliberately *not* allowed to vary further *within* a chromosome: a
region with unusually slow local decay (centromere, inversion) is precisely where
markers are most redundant and pruning should be most aggressive, so adapting the
threshold to tolerate high local LD there would defeat the purpose.

Clustering proceeds in two steps on the precomputed pairwise LD edge list:

1. Markers are grouped into **connected components** of the thresholded sparse graph
   (inexpensive).
2. Each component is refined into **complete-linkage sub-clusters** by hierarchical
   clustering on the dense `r²` sub-matrix within it.

Step 2 is necessary because single linkage alone permits **chaining**: two markers
uncorrelated with each other can be joined through an intermediate correlated with
both, and would be wrongly treated as redundant. Complete linkage requires *every*
pairwise `r²` within a final cluster to clear the threshold.

The **representative** of each cluster is the marker with the **highest median `r²`
to the remaining members**. This is preferred over "most correlated with the
cluster's mean genotype" because the latter is sensitive to allele-coding direction
(a genuinely central marker coded on the opposite allele correlates negatively with
the mean and would be passed over); `r²` ignores sign, so the median-`r²` criterion
has no such failure mode.

### Why Stage 1 is conservative, and how Stage 2 resolves it

The pairwise edge list is built within a sliding window, so pairs farther apart than
that window were never evaluated and behave as `r² = 0` in the clustering. Single
linkage tolerates this (a chain of overlapping windows still connects distant
markers); **complete linkage does not** — a genuine single cluster wider than the
window is fragmented purely because unevaluated pairs default to zero. The window is
sized from *average* decay behaviour, which covers typical regions comfortably, but
slowly-recombining regions decay far more slowly and are where fragmentation is
expected.

Stage 2 addresses this *after the fact* rather than by changing Stage 1. Candidate
clusters are compared directly from genotypes with **no window restriction**, and
are selected for reconsideration using Stage 1's own cluster boundaries. Markers are
deliberately **not** pre-split by their local LD before clustering: doing so would
sever real clusters exactly at the selection boundary, whereas selecting whole
clusters leaves genuine clusters intact.

### Stage 2 — quality-gated merging

A Stage-1 cluster is **flagged** for reconsideration if any member exceeds
**`ld_w > 0.025`**, *or* if the cluster contains **≥ 5 markers**. The `ld_w`
criterion follows from what the statistic measures — a marker with low local support
is already near-independent and would not merge with anything — so restricting the
(comparatively expensive) reconsideration to markers with appreciable support
concentrates it exactly where a merge is possible, without omitting merges elsewhere.
The ≥ 5-marker route catches clusters that demonstrably *do* cluster but whose local
support is only moderate.

Flagged clusters are merged by a **distance-restricted, average-linkage dynamic tree
cut** that accepts each candidate merge only while **two quality conditions** hold:

- **Consensus fidelity floor:** the consensus must survive hard-calling, quantified
  as `score_eMLG(x) = cor(round(x), x)² ≥ 0.80` — the property downstream LD
  statistics require of a consensus.
- **Pairwise correlation floor:** the two sides being merged must themselves be
  correlated, at **`r² ≥ 0.2`**.

The two conditions are **not redundant**. Fidelity is a property of the *resulting*
consensus, and once a cluster is large that consensus is dominated by the markers
already in it, so appending a few more changes the score very little whether or not
they report the same partition — the fidelity criterion is least discriminating
exactly where clusters are largest, and the pairwise floor is what stops an
essentially uncorrelated marker being absorbed there. Because merging proceeds only
while *both* conditions hold, **this stage cannot over-merge regardless of how
permissively clusters were flagged.**

*(Note: `r² ≥ 0.2` is a genuine merge gate, not a file-size convenience — it bounds
the residual correlation that can survive between distinct clusters.)*

**Residual independence.** The two floors bound the correlation that can survive
between distinct clusters. Measured after the fact on the merged units — and after
correcting for relatedness, so shared structure does not inflate it — the largest
residual is of order **`r² ≈ 0.3`**, comparable to the fixed `r² = 0.2` that
conventional LD pruning applies uniformly, but here on a relatedness-corrected scale
between consensus genotypes. Any residual is moreover *local* (neighbours within a
short genetic distance); where independent observations are needed across long
distances (e.g. two-locus epistasis screens) a genetic-map distance threshold
removes it directly.

**Distance restriction.** Applied in **genetic units, 0.5 cM**, using the linkage
map. Physical distance is a poor proxy for recombination distance in this genome:
cluster pairs > 1 Mb apart can be at ≈ 0 genetic distance, while others only a few
hundred kb apart span several cM. A genetic threshold also behaves correctly on
exceptional chromosomes, where a slower local recombination rate gives 0.5 cM a
proportionally larger physical span.

### Outputs

The procedure returns both projections: **one representative marker per cluster**,
and a **consensus dosage for clusters of ≥ 5 markers** (`min_n_loci_eMLG = 5`). The
size threshold is a computational limit on consensus construction, *not* a
redefinition of the unit: smaller clusters remain valid units and keep their pruning
representative, and the pruned marker set is unchanged by this setting.

For the canonical clustering `eMLG_5loci_0025_cM05.rds`:

> **474,014 LD clusters** in total, of which **32,840 contain ≥ 5 markers** and
> therefore carry a stored consensus (eMLG) genotype.

Cluster group IDs (e.g. `F59859`) are **not stable across clusterings** — never join
on group ID across different clustering files. The filename encodes the parameters
(`5loci` = min 5 markers for a consensus; `0025` = `ld_w` flag 0.025; `cM05` = 0.5 cM
backstop) so that an analysis cannot silently run against a different clustering.

### Figures

- **Fig. S0.3 — `../Figures/eMLG_score_histograms_by_threshold.png`.** Distribution
  of consensus round-trip fidelity (`score_eMLG`) across flagged clusters at five
  `ld_w` flagging thresholds, separating single-marker clusters (no merge possible)
  from genuine multi-marker merges. The separation matters for reading the mass near
  the upper limit, which at permissive thresholds is largely isolated markers with no
  correlated neighbour, not merges stopping early.
- **Fig. S0.4 — `../Figures/Chr26_stage1_vs_combined_high.png`.** Stage-1 clustering
  vs the combined two-stage result for Chr26 (markers with high local LD support).
  Fragmented Stage-1 clusters within each low-recombination region are consolidated
  into single coherent clusters per physical LD peak — the intended effect of
  Stage 2. (Companion panels for Chr1/2/3/5/10 exist as `Chr*_stage1_vs_combined_*.png`.)

---

## Parameter table (Stage 0)

| Group | Parameter | Value |
|---|---|---|
| Marker set | site type | biallelic |
| | MAF among hybrids | ≥ 0.05 (1,114,340 markers) |
| LD decay | background pairs | inter-chromosomal, 5,000 markers |
| | windows / chromosome | 20, 50% overlap |
| | pairs / window | ≤ 5,000, 20 distance strata |
| | `r²` quantile | 0.95 |
| | MAF for curve fitting | > 0.1 (both members) |
| | robust aggregation | central 95% of windows |
| | local support `ld_w` | median `r²` within ρ = 0.95 window |
| Stage 1 | clustering threshold | per chromosome at ρ = 0.5 |
| | linkage | components, then complete linkage |
| | representative | highest median `r²` to cluster |
| Stage 2 | flagging | `ld_w > 0.025` OR ≥ 5 markers |
| | consensus fidelity floor | `score_eMLG ≥ 0.80` |
| | pairwise correlation floor | `r² ≥ 0.2` |
| | distance restriction | 0.5 cM (genetic map) |
| | consensus size threshold | ≥ 5 markers |

**Output:** `data/eMLG_5loci_0025_cM05.rds` → 474,014 clusters; 32,840 with a stored
consensus.

---

*Next: [Module A — the sorting phenomenon](02_moduleA_sorting_phenomenon.md).*
