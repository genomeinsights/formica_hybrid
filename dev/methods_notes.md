# Methods notes: two-stage LD-pruning / eMLG complexity reduction

Running notes on the rationale behind `ld_complexity_reduction()` +
`ld_prune_and_eMLG()` (both part of the LDscnR package -- see
`~/gitlab/LDscnR/R/ld_prune_and_eMLG.R`, `~/gitlab/LDscnR/R/dynamic_cut_eMLG.R`
-- and this repo's `R/ld_pruning_DIEM.R`). Written for eventual reuse in
supporting information.

## The unifying idea

LD-pruning and eMLG (expected multi-locus genotype) complexity reduction
are not two unrelated procedures — they are two different reductions
applied to the *same* partition of markers. Both start from the same
question: what is the finest partition of markers that we can no longer
statistically distinguish, i.e. what constitutes an independent test unit
in this data set? Once that partition is fixed:

- **LD-pruning** picks one representative marker per unit.
- **eMLG complexity reduction** collapses each unit into a single
  consensus dosage genotype for downstream pairwise LD / Ohta-statistic
  analysis.

Both are answering the same underlying question, just projecting the
answer differently. This is why a single combined function
(`ld_prune_and_eMLG()`) can produce both outputs from one clustering pass.

## Stage 1: `ld_complexity_reduction()` (complete linkage, rho = 0.5)

Stage 1 clusters markers genome-wide using complete linkage on pairwise LD
decay, anchored at `rho = 0.5` — the half-decay point of the fitted LD
decay curve (`r² ~ b + (c-b)/(1+a·d)`), rather than an arbitrary fixed r²
cutoff. This choice is deliberate: complete linkage is the strictest
possible aggregation rule (every pairwise comparison within a cluster must
clear the bar), so it is structurally biased toward fragmenting any block
that has real internal heterogeneity in LD decay. Anchoring the per-edge
criterion to the decay curve's own half-point, instead of a stricter fixed
threshold, avoids compounding that structural bias unnecessarily.

In a young, low-recombination hybrid population, true haplotype blocks can
be large (few recombination events have had the chance to break them up).
Complete linkage's strictness means Stage 1 will tend to over-split these
genuinely large blocks into more, smaller clusters than the population
history would justify. This over-splitting is expected and is corrected
in Stage 2, not avoided in Stage 1.

## Stage 2: quality-gated dynamic cut (average linkage)

`ld_prune_and_eMLG()` re-merges Stage 1 clusters that were flagged as
having high local LD support (`ld_w_col > ld_w_threshold`, evaluated per
Stage-1 cluster, not per marker — see below) using a distance-restricted,
average-linkage dynamic tree cut (`dynamic_cut_eMLG()`). The cut walks the
average-linkage dendrogram bottom-up and accepts each candidate merge only
while two quality gates hold:

- `score_eMLG(x) = cor(round(x), x)^2 ≥ score_threshold` — does the
  consensus genotype survive hard-calling (round-trip fidelity), which is
  what downstream LD/Ohta statistics need.
- `pair_r2 ≥ min_r2` — are the two sides actually being merged correlated
  with each other (needed because `score_eMLG` alone does not verify
  this; a merge with pair r² ≈ 0 can still pass a `score_eMLG` gate on its
  own).

Because merging only ever proceeds while quality holds, Stage 2 cannot
over-merge: `score_threshold` is a hard floor on every accepted merge,
independent of how permissively clusters were flagged for reconsideration
in the first place.

### Cluster-level, not marker-level, flagging

Clusters are flagged for Stage 2 by whether *any* member exceeds
`ld_w_threshold`, not by discarding individual high-ld_w markers. Flagging
at the marker level, checked directly, artificially severs real LD blocks
right at the threshold boundary — verified on Chr26: 226 of 15,524 Stage-1
clusters straddled the threshold, together 4,595 markers of which 975 were
low-ld_w members that a marker-level pre-split would have cut loose from
their true cluster. Flagging at the cluster level means the flagged
("high") bucket legitimately contains some sub-threshold "boundary"
markers (cluster-mates of a flagged marker) by design, while the
unflagged ("low") bucket is guaranteed by construction to contain zero
markers above threshold (verified empirically: 0/14,316 in one real-data
test).

### Reading the score histogram

The distribution of final `score_eMLG` values across flagged clusters is a
useful diagnostic precisely because the dynamic cut is quality-gated and
only stops merging when a further merge would break the `score_threshold`
floor. Empirical genome-wide runs at five thresholds (`ld_w_threshold` =
0.025, 0.05, 0.1, 0.15, 0.2, all with `min_n_loci_eMLG = 5`) produced a
figure for this (`R/plot_eMLG_score_histograms.R`, saved to
`Figures/eMLG_score_histograms_by_threshold.png`, intended for
supplementary results) that clarified two things:

- **A real qualitative shift with threshold, not just a shift in mean.**
  At low thresholds (0.025, 0.05) the score distribution has a genuine
  interior mode around 0.85–0.90, with only a modest tail near 1.0. At
  higher thresholds (0.1–0.2) the distribution instead increases
  monotonically toward a large spike at the 1.0 ceiling. Mean score
  actually moves in the *opposite* direction to what the ceiling-spike
  framing would suggest (mean 0.894 at th=0.025 vs 0.937 at th=0.2) —
  the mean is not the useful summary statistic here, the shape is.
- **Mass near the ceiling is NOT uniformly "under-merging" — its
  composition must be checked.** Splitting the ceiling bin (score > 0.99)
  by whether the cluster is a true singleton (`n_loci == 1`, no merge was
  ever possible) or a genuine multi-locus merge changes the story
  substantially:
    - At low thresholds, the ceiling bin is 85–91% trivial singletons —
      clusters that were flagged but had no viable correlated neighbour to
      merge with at all. This is not "wasted headroom," it simply reflects
      that low-`ld_w_threshold` flagging pulls in many isolated markers
      that have nothing nearby to merge into.
    - At `ld_w_threshold = 0.2`, the ceiling bin is mostly genuine merges,
      and the multi-locus ones landing there have a median of 35 loci —
      i.e. real, large haplotype blocks getting fully consolidated at
      near-perfect fidelity. That is a success case, not evidence of
      stopping early.
  The original framing ("mass near 1.0 always means under-merging") was
  too simple; whether it does depends on what fraction of that mass is
  singleton vs. genuine merge, which the figure now shows directly via a
  singleton/merged fill split.
- **The cleanest unambiguous signal for the diminishing-returns argument
  is the low tail, not the ceiling**: the fraction of flagged clusters
  sitting close to the 0.80 floor (score < 0.85) drops from 10.6% at
  th=0.1 to 6.9% at th=0.2 — i.e. raising the threshold leaves fewer
  clusters stuck near-floor, unmerged. This tracks threshold monotonically
  and isn't confounded by the singleton-composition issue above.

### Effect of `ld_w_threshold`

Lowering `ld_w_threshold` flags more Stage-1 clusters into the
distance-restricted runs Stage 2 operates on, giving the dynamic cut more
raw material to merge across — and, per the histogram composition finding
above, also pulls in many more clusters that turn out to have no viable
merge partner at all (pure singletons). The number of flagged clusters
scales steeply as threshold drops (637 at th=0.2 → 64,039 at th=0.025 in
the genome-wide run, roughly 5x per halving), so this is also the main
computational cost lever: because `score_threshold` is a hard floor
regardless of how many clusters are flagged for reconsideration, a lower
`ld_w_threshold` can never cause over-merging, only more (and slower)
opportunity to correct Stage 1's over-splitting. Practically, `th=0.025`
(retaining a real, discriminating cutoff) was the value settled on for
production use, on top of visual inspection of Stage 1 vs. combined
per-chromosome diagnostic plots (`plot_pruning_comparison()`,
`direction="high"`): on both Chr1 (largest) and Chr26 (one of the
smallest) at th=0.025, the previously-fragmented rainbow-noise clusters in
each low-recombination region collapse into one coherent, single-colored
cluster per distinct physical LD peak, which is the target shape for
downstream eMLG-based LD/Ohta analysis.

### A threshold of exactly 0 is not "no flagging cutoff"

It's tempting to think `ld_w_threshold = 0` disables the ld_w-based
flagging criterion, leaving only some other, cheaper criterion (like a
minimum cluster size) to do the restricting. It doesn't: checked directly
on the full marker set, only 23 of 1,114,340 markers (~0.002%) have
`ld_w_095` exactly 0 — effectively 100% of markers have `ld_w_095 > 0`. So
`ld_w_threshold = 0` floods the flagged set to nearly the entire genome on
its own, defeating the whole point of restricting Stage 2's expensive
O(n²) all-pairs-correlation dynamic cut to a manageable subset (see
LDscnR's `R/ld_prune_and_eMLG.R` header comments on why this scaling is
load-bearing, not a workaround to relax).

### `min_n_loci_flag`: pulling substantial low-ld_w clusters into merging, cheaply

Checking the *unflagged* side at th=0.025 turned up 2,804 Stage-1 clusters
with ≥5 loci that were passing straight through unmerged (their internal
members are strongly linked enough for Stage 1's complete linkage to have
grouped them, but their local `ld_w` support wasn't high enough to get
flagged). Two ways to get eMLGs for these:

1. `compute_unflagged_eMLG = TRUE` + `min_n_loci_eMLG = 5`: computes a
   consensus dosage directly from each unflagged cluster's *existing*
   Stage-1 members, with no re-clustering. Cheap (linear pass, no O(n²)
   step) but these clusters never get a chance to merge with a
   physically-nearby, correlated neighbour.
2. `min_n_loci_flag = 5` (new parameter, default `Inf` i.e. off, so
   existing calls/saved results are unaffected): additionally flags any
   Stage-1 cluster with ≥5 loci for the full distance-restricted dynamic
   cut treatment, on top of the `ld_w_threshold` criterion, giving these
   clusters a chance to merge. Verified on a real Chr26+Chr10 subset: 173
   previously-unflagged clusters with ≥5 loci were correctly pulled into
   the flagged pathway, total markers covered stayed identical before and
   after (66,080), and the unflagged bucket's max `n_loci` dropped to 4 as
   expected (nothing ≥5 loci remained unflagged). Because this ORs with a
   real `ld_w_threshold` rather than replacing it, the added flagged-set
   size stays small and controlled (2,804 additional clusters genome-wide
   at th=0.025, vs. 64,039 already flagged from ld_w alone) — unlike
   `ld_w_threshold = 0`, which floods the flagged set regardless of any
   size condition ORed onto it (see above).

### `min_n_loci_eMLG`: a compute-cost gate, not a redefinition of the unit

`min_n_loci_eMLG` sets a minimum cluster size below which an eMLG
consensus genotype is not computed, purely because pairwise LD/Ohta
statistic calculations downstream still scale expensively with the number
of eMLG units and a line has to be drawn somewhere. Critically, this does
not change the underlying partition of markers, and therefore does not
change what counts as an independent test unit for pruning purposes: a
cluster below the size floor is still a fully valid unit and still gets a
pruning representative — it simply is not worth the downstream compute
cost to also produce a consensus eMLG for it (a singleton's dosage is
already directly usable with no consensus step needed anyway). This was
verified directly: the pruned marker set produced by
`ld_prune_and_eMLG()` is byte-identical regardless of the
`min_n_loci_eMLG` setting — only the `eMLG` matrix and `groups$has_eMLG`
column are affected.
