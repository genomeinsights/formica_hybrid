# The direction of predictable ancestry sorting is governed by locus diagnosticity

*Draft results summary (main text). All methods are deferred to Supplementary
Materials. Figure/Table callouts are placeholders. Findings are descriptive:
each pattern below is, on its own, compatible with neutral drift as well as
with selection — establishing departure from neutrality and discriminating its
cause are the subject of the companion analyses noted at the end.*

---

## Results

**Replicated ancestry sorting across twenty hybrid populations.** Twenty
independently formed hybrid populations between *Formica polyctena* and *F.
aquilonia*, together with allopatric parental references, allow the direction
of ancestry sorting to be tested locus by locus across replicates. Under drift,
or under an incompatibility in which the two parental combinations are equally
fit, populations near-fix ancestry in *random* directions; under directional
selection, or an unequal-fitness Dobzhansky–Muller incompatibility, the *same*
parental allele fixes repeatedly. This is the per-locus, twenty-replicate
generalisation of the three-population parallelism of Nouhaud et al. (2022).

**Most differentiated loci do not sort; those that do are directional, not
bidirectional.** Restricting to loci that are polymorphic in the parents
(pooled-parental minor-allele frequency ≥ 0.15; ~660,000 loci — a gate that
prevents alleles already near-monomorphic in the founding pool from producing
spurious fixation), the majority of loci remain unsorted across the twenty
replicates (67.4%). Where sorting occurs it is overwhelmingly *unidirectional*:
17.3% of loci near-fix toward *polyctena* and 14.3% toward *aquilonia*, whereas
bidirectional near-fixation — populations fixing in opposing directions, the
signature expected under drift or equal-fitness incompatibility — is
vanishingly rare (0.1%) (Fig. 1a). The consistency of direction across
replicates, rather than its mere magnitude, is the informative signal.

**The polarity of sorting reverses along the diagnostic-index axis.** The net
excess of *polyctena*-directed sorting is not uniform: the direction of
unidirectional sorting reverses systematically with a locus's diagnostic index
(DI), the degree to which it distinguishes the parental species. Among
unidirectionally sorted loci, the fraction fixing toward *aquilonia* rises
monotonically from ~0.15 in the least-diagnostic decile to ~0.74 in the
most-diagnostic decile, crossing parity near the middle of the DI range
(Fig. 1b; Spearman ρ = 0.31). In a logistic model the DI effect is large and
independent of parental allele frequency (P(aquilonia | sorted) ~ DI +
parental MAF; DI z = 181; parental MAF z = −28, the latter indicating that the
most polymorphic loci lean *polyctena*). Thus **aquilonia ancestry is fixed
preferentially at the most diagnostic loci, and polyctena ancestry at the
least diagnostic** — and restricting attention to strongly diagnostic loci
(DI above the median) recovers a ~70% *aquilonia* majority, reconciling the
apparent excess of either parent as two ends of a single diagnostic axis.

**The signal resides in many small units, not in large linkage blocks.**
Because linked markers report the same underlying event, we collapsed
LD-correlated markers into independent units (consensus expected multi-locus
genotypes) and re-tested each unit once. Sorting is concentrated in numerous
small units: 74.7% of large linkage blocks (≥ 5 constituent markers) show *no*
aggregate sorting once their members are summarised, versus only 5.0% of the
smallest units (Fig. 1c). This collapse is genuine and not an artifact of
averaging heterogeneous markers — the consensus genotypes of the washed-out
blocks retain high internal fidelity (median round-trip score 0.86; fewer than
0.05% fall below the fidelity threshold used to construct them). The
per-marker impression of pervasive sorting in large blocks is therefore
spatial pseudo-replication; the independent-unit signal is carried by many
small, largely unlinked loci.

**Cluster size predicts direction only through diagnostic index.** Larger
linkage blocks lean *aquilonia* and are more diagnostic, so sorting direction
appears to track genomic architecture. However, this is entirely mediated by
DI: in a model including block size, DI and parental MAF, the size term
*reverses* sign once DI is controlled (P(aquilonia) ~ log₂(block size) + DI +
parental MAF; size z = −22, DI z = 120). Diagnostic index, not linkage-block
size or recombination architecture, is the governing axis of directional
sorting.

**Interpretation and outlook.** Together these results establish a
*predictable, directionally structured* pattern of ancestry sorting whose
polarity is organised by locus diagnosticity rather than by genomic
architecture. We emphasise that these patterns are descriptive: near-fixation
in a consistent direction is expected under directional selection or an
unequal-fitness incompatibility, but a formal neutral expectation is required
before either can be inferred. We therefore evaluate the sorting against a
recombination-matched, haplodiploid coalescent-with-selection null seeded from
observed parental haplotypes (Supplementary), and separately test the two
mechanistic hypotheses it cannot distinguish on its own: *intrinsic*
incompatibilities, through elevated among-population linkage disequilibrium
between unlinked diagnostic units (a two-locus Ohta decomposition), and
*extrinsic*, climate-driven selection, through the overlap of directionally
sorted units with genotype–environment association outliers. These analyses,
reported next, are designed to convert the descriptive polarity documented here
into a test of whether predictable sorting is driven by genetic incompatibility
or by ecological adaptation.

---

*Fig. 1 (placeholder). (a) Genome-wide sort-class proportions among
parent-polymorphic loci. (b) Fraction of unidirectionally sorted loci fixing
toward aquilonia across diagnostic-index deciles. (c) Proportion of independent
units unsorted as a function of linkage-block size.*
