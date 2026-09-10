# Audit of the first-pass population-partitioning analyses

## Overall assessment

This module is promising and already provides descriptive support for the core
hypothesis: high differentiation at one LD-reduced unit does not imply that
nearby units partition the hybrid populations in the same way. However, several
statements in the current `README.md` and draft Results are stronger than the
analyses presently justify and should be revised before manuscript use.

The clearest result is the rapid decay of similarity between population
profiles with physical distance:

- Mean signed Pearson correlation falls from 0.258 at <5 kb to 0.076 at
  20–100 kb and 0.024 beyond 10 Mb.
- Mean absolute correlation falls from 0.296 at <5 kb to approximately 0.191
  beyond 10 Mb.
- Random cross-chromosome pairs have mean absolute correlation 0.193.

Thus, there is clear excess concordance at short distances, but little evidence
that units separated by longer distances repeatedly distinguish the same
population subsets. This pattern is consistent with local physical linkage.
Physical distance alone does not demonstrate that the pattern is LD-driven, so
the current wording should say “consistent with LD” unless pairwise LD or local
recombination is analysed directly.

The relationship between per-unit FST and local profile concordance is positive
but weak:

- FST versus local mean absolute correlation: Spearman rho = 0.081.
- FST versus local mean signed correlation: Spearman rho = 0.114.
- Mean local absolute correlation increases from 0.212 in the lowest FST
  quartile to 0.224 in the highest quartile.
- Mean local signed correlation increases from 0.087 to 0.118 across those
  quartiles.

A suitable current interpretation is:

> More highly differentiated units showed slightly greater local concordance,
> but the effect was weak and most neighbouring high-FST units did not share
> strongly concordant population profiles.

The results are compatible with the proposed explanation for the empirical
FST–DI pattern, but do not yet demonstrate that independent, population-specific
sorting explains the discrepancy between empirical and simulated FST.

## Data preparation and conventions

The data preparation appears internally consistent with the cleaned DI25
pipeline:

- The analysis uses the 11,052 LD-reduced units from the DI25-specific cM5
  clustering.
- Units with more than two markers are represented by the best SNP selected by
  `eMLG_best_snp(fill = FALSE)`; strictly observed SNP genotypes are used, not
  eMLG consensus-filled genotypes.
- One- and two-marker clusters use their stored representative SNP.
- Alleles are oriented toward *F. aquilonia* using the parental samples.
- The population-frequency matrix contains 20 hybrid populations × 11,052
  units and has negligible missingness.
- FST is calculated among the hybrid populations using the same Weir and
  Cockerham estimator as the DI25 pipeline.
- Sorting classifications reuse the locked Module A settings rather than
  introducing a new threshold.

These choices should be retained and stated clearly.

## Findings requiring correction

### 1. Large clusters do not drive the genome-wide FST–concordance association

The current summaries say that the weak genome-wide association is “almost
entirely driven” by 48 clusters containing more than 50 SNPs. The reported
statistics do not support this statement:

- overall rho = 0.081;
- rho after excluding clusters with `n_loci_g > 50` = 0.079.

Removing these clusters therefore changes the genome-wide association
negligibly. The 48 large clusters have a stronger FST–concordance relationship
within that subset (rho = 0.41), but they are a distinct subset rather than the
cause of the weak genome-wide trend.

In addition, only 4 of the 48 `n_loci_g > 50` clusters belong to the three named
polyctena blocks. The statement that these clusters are “mostly” the three named
blocks is incorrect. The three named regions contain 10 LD units in total, only
four of which exceed 50 markers.

Finally, `n_loci_g > 50` identifies large LD clusters, not genomic regions and
not low-recombination regions directly. Avoid calling all 48 units “giant
low-recombination blocks” without joining the recombination information or an
independent region definition.

### 2. Signed correlation should be primary

Because all alleles are oriented toward *F. aquilonia*, signed correlation has
direct biological meaning:

- positive r: the same populations have relatively more aquilonia ancestry at
  both units;
- negative r: the same division of populations occurs in opposing ancestry
  directions;
- r near zero: the units distinguish different population subsets, or their
  profiles are too noisy to show a common division.

Absolute correlation treats r = +0.8 and r = -0.8 as equally concordant. This
is useful only for the explicitly direction-independent question “do the units
divide the same populations, regardless of ancestry direction?”

Report signed r as the primary biological statistic and absolute r as a
secondary measure of direction-independent partition similarity. Signed r also
avoids the positive sampling floor inherent in absolute correlations estimated
from only 20 observations.

### 3. Use empirical cross-chromosome pairs as the main baseline

The population-label permutation gives mean absolute r = 0.187, whereas genuine
cross-chromosome pairs give 0.193. Cross-chromosome comparisons preserve the
real genome-wide ancestry differences and covariance among populations and are
therefore the more appropriate baseline for asking whether physical neighbours
show excess concordance.

Use cross-chromosome or otherwise distance/location-matched empirical pairs as
the primary reference. The label-permuted distribution can remain as a
secondary demonstration of the finite-sample absolute-correlation floor.

The statement that the 100–500 kb bin is “statistically indistinguishable” from
the permutation floor has not yet been formally tested. Until a block-aware
comparison is made, describe it as close to the empirical background rather
than statistically indistinguishable.

### 4. Current confidence intervals treat non-independent observations as independent

The error bars in the distance-bin and FST-decile plots are based on
`sd / sqrt(N)`, with genomic pairs or units treated as independent. Units occur
in multiple pairs, and neighbouring units are correlated. With millions of
pairs, these standard errors are consequently much too small.

Estimate uncertainty by resampling chromosomes or genomic blocks. Apply the
same block-bootstrap principle to both the distance-decay curve and the
FST-decile trend.

### 5. Euclidean-distance calculation appears incorrect

In `R/pp_local_concordance.R`, Euclidean distances for all chromosome pairs are
calculated with `colSums()` on an object in which a distance is required for
each row/pair. This appears to require `rowSums()` instead. The resulting means
around 100 are implausible for standardized 20-population profiles and confirm
that the current Euclidean values are not valid.

The Euclidean results are not used in the present biological conclusions, so
this does not affect the reported correlation results. Correct the calculation
and verify its dimensions, or remove Euclidean distance from the module.

### 6. PCA interpretation requires refinement

Row-centring each unit before PCA is appropriate for studying the shape of its
population partition. PC1 explaining 11% of the variance is inconsistent with
one overwhelmingly dominant shared partition, but it does not establish that
the partitions are independent or that no recurring partition exists.

Add:

- a row-wise population-label permutation reference for PC1 variance;
- inspection and plotting of the population loadings for PC1 and subsequent
  axes;
- leave-one-population-out PCA, especially excluding Sielva;
- separate summaries for aquilonia-sorted and polyctena-sorted units.

The existing explanation of the 71% PC1 obtained from the uncorrected PCA also
appears backwards. With units as rows and populations as columns,
`prcomp(center = TRUE)` centres each population across units. Failure to
row-centre the units leaves variation in their overall allele-frequency levels,
which can create a dominant intercept-like axis. Verify this explicitly before
describing the 71% axis as the between-population hybrid-index baseline.

## Most informative next analysis: remove genome-wide ancestry

The current correlations are calculated from raw population allele-frequency
profiles. Those profiles contain both genome-wide differences in ancestry
among populations and locus-specific departures from that background. The
most direct next analysis is therefore to residualize each unit against
genome-wide ancestry.

For each population, calculate genome-wide ancestry using markers on other
chromosomes, preferably as a leave-one-chromosome-out estimate. For each unit
u, fit:

```text
f_pu = alpha_u + beta_u H_p + residual_pu
```

where `H_p` is the population’s leave-one-chromosome-out genome-wide ancestry.
Then calculate:

1. the proportion of each unit’s among-population variation explained by
   genome-wide ancestry;
2. its relationship with FST;
3. signed and absolute local concordance between the residual profiles;
4. whether high-FST units still identify many different residual population
   partitions.

This distinguishes two biologically different explanations:

1. high-FST units repeatedly reflect the same genome-wide ancestry gradient;
2. high-FST units are driven by different populations departing from their
   genome-wide ancestry at different loci.

The second result would provide substantially stronger empirical support for
the proposed population-partitioning mechanism.

## Other priority robustness checks

After correcting the points above, proceed in this order:

1. Leave-one-population-out analysis, with particular attention to Sielva and
   any established shared-origin population pairs.
2. Chromosome- or region-block bootstrap confidence intervals.
3. Literal once-per-region treatment of the three named large sorting regions,
   rather than using `n_loci_g > 50` as a proxy.
4. Focused within-region analysis of the named Chr5, Chr25 and Chr26 regions.
5. Join the recombination map and test whether concordance decay varies with
   recombination rate; until then use “consistent with linkage,” not
   “LD-driven.”
6. Inspect which populations contribute most strongly to each high-FST unit
   and whether a small subset of populations recurs more often than expected.

No new evolutionary simulations are required for this stage. A comparison
with population profiles from the existing simulations will eventually be
needed before claiming that the observed partitioning explains the empirical
versus simulated FST difference, but that comparison can be deferred.

## Figure assessment

### Figure 1: Chr26 heatmap

The heatmap visually demonstrates heterogeneous population profiles, but all
177 Chr26 units are shown as equal-width columns. This obscures physical
distances and does not clearly identify the four-unit named block. Add a focused
panel for the named region with physical coordinates and unit/cluster-size
annotation. Retain physical order rather than clustering genomic columns.

### Figure 2: similarity versus distance

This is currently the clearest result. Revise it to:

- show signed correlation as the primary panel and absolute correlation as the
  secondary panel;
- use the empirical cross-chromosome background as the main reference;
- show chromosome/block-bootstrap uncertainty rather than pairwise standard
  errors.

### Figure 3: FST versus local similarity

This is useful, but describe the increase as small. Lead with signed
correlation, retain absolute correlation as a sensitivity measure, and obtain
uncertainty by block bootstrap. Do not state that the pattern is driven by the
48 large clusters.

### Figure 4: genome-wide summary

The present version adds relatively little because most observations occupy a
narrow range and the colour scale does not isolate the important regions. A
region-focused plot or a plot of ancestry-residual concordance may be more
informative.

## Recommended manuscript-level conclusion at this stage

> Similarity between the ancestry-oriented population profiles of high-DI
> LD-reduced units declined rapidly with physical distance. More highly
> differentiated units showed slightly greater local concordance, but the
> relationship was weak, and most high-FST units did not share strongly
> concordant population profiles with nearby units. The results are therefore
> consistent with high-DI units frequently capturing different population
> partitions, although analyses controlling for genome-wide ancestry are
> required before attributing the empirical–simulation FST difference to this
> mechanism.

