# Structured vs. unstructured null: local-score window comparison

Status: independently regenerated local-score computation on all 40
full-SNP null sweep outputs, using `compute.local.scores()` (Fariello et
al. 2017 / Bonhomme et al. 2019) exactly as called elsewhere in this
pipeline (`xi=1, pval.local.score.thres=0.01, min.maf=0.2, min.nsnp=1e4`),
via `R/03_compute_local_scores_fullsnp.R` (fresh implementation, not
sourced from the exploratory module). Per-SNP detail (`res.local.scores`,
~19MB/run, 751MB total for 40 runs) is kept on disk under
`full_snp_null10/localscores/` (gitignored, regenerable from
`full_snp_null10/raw/` + this script); only the compact
`significant.windows` extract (`all_significant_windows.rds`/`.tsv`, 2,219
rows total across all 40 runs) and this summary are committed.

## Window counts per run

| mode | nulltype | mean | sd | min | max |
|---|---|---|---|---|---|
| continuous | structured | 4.6 | 3.47 | 0 | 10 |
| continuous | unstructured | 4.5 | 2.72 | 1 | 8 |
| mitoC2 | structured | 112.9 | 18.23 | 86 | 147 |
| mitoC2 | unstructured | 99.9 | 13.75 | 74 | 122 |

Per-draw values:

- continuous structured:   3, 9, 3, 9, 2, 4, 0, 2, 10, 4
- continuous unstructured: 7, 3, 8, 3, 2, 2, 4, 8, 7, 1
- mitoC2 structured:   98, 114, 118, 134, 104, 123, 86, 100, 147, 105
- mitoC2 unstructured: 74, 97, 106, 91, 93, 104, 122, 95, 118, 99

## Paired comparison (structured vs. unstructured, matched by replicate)

Each replicate's unstructured draw is a within-pair permutation of its
structured counterpart (same value multiset, see
`config/null_covariate_design.md`), so a paired test on the 10 replicates
is the correct comparison -- it isolates the effect of Omega-linked
population structure while holding the covariate's value distribution
fixed.

| mode | mean difference (structured - unstructured) | paired t-test p | paired Wilcoxon p |
|---|---|---|---|
| continuous | +0.1 windows | 0.946 | 1.000 |
| mitoC2 | +13.0 windows | 0.078 | 0.064 |

## Interpretation

- **Continuous mode**: no detectable difference between structured and
  unstructured nulls (mean difference 0.1 windows, both tests p>0.9). At
  this resolution and threshold, Omega-linked structure in a continuous
  covariate does not measurably inflate local-score window detection
  beyond what an unstructured covariate with the same value distribution
  already produces.
- **mitoC2 mode**: structured consistently produces more windows than its
  paired unstructured counterpart in most replicates (structured > 
  unstructured in 8/10 pairs), with a mean excess of 13 windows (112.9 vs.
  99.9, ~13% relative increase). This does not reach conventional
  significance at n=10 replicates (paired t p=0.078, Wilcoxon p=0.064) --
  an honest, appropriately-powered null result at this replicate count,
  not a confirmed effect. It is directionally consistent with the
  hypothesis that Omega-covariance structure inflates local-score false
  positives specifically for the binary contrast statistic, but this
  sweep alone does not establish that at a conventional threshold.
- The mode asymmetry (mitoC2 far more windows than continuous, ~20-25x, in
  both null types) reflects the underlying statistics' different scales
  and detection thresholds (log10(1/pval) contrast vs. BF(dB) regression),
  not the structured/unstructured contrast itself -- this magnitude
  difference is expected and not itself evidence of anything.

## What remains for the confirmatory document

- This is a null-vs-null comparison (does population structure in the
  covariate inflate detection); it says nothing yet about the OBSERVED
  scans (PC1/bio_winter/mitoC2) relative to either null. That comparison
  (observed window count/positions vs. this null distribution) is the next
  step and belongs in `doc/`.
- The Stage-1 1,000-replicate calibration (using `best_marker` positions,
  see `data/stage1_unit_positions.tsv`) is separate, not yet run.
- Whether to retain a larger replicate count for the mitoC2 comparison
  (to resolve the p=0.06-0.08 result one way or the other) is worth
  deciding explicitly rather than silently extending the sweep.
