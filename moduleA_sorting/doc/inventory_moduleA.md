# Module A documentation inventory

Status checked against the current scripts and RDS objects on 4 August 2026.
The scripts and parameter-stamped RDS objects are treated as authoritative when
they disagree with older prose.

## Scope

Module A asks how often parental ancestry approaches fixation across the 20
hybrid populations, whether the same parental direction is repeated, and how
sorting magnitude and direction relate to diagnostic index (DI), recombination
and LD redundancy.

## Locked primary definitions

- pooled parental MAF at least 0.15;
- near-fixation floor phi = 0.85;
- sorting-magnitude threshold tau = 0.60;
- DI retained as an ungated predictor;
- canonical Module-0 0.5-cM partition;
- direction classified by an exact binomial random-direction test
  (null probability 0.5, alpha = 0.05), not by the earlier component rule.

## Authoritative headline outputs

- Of 661,312 parental-polymorphic SNPs, 62,504 (9.45%) sort toward
  *F. aquilonia*, 75,416 (11.40%) toward *F. polyctena*, and 785 (0.119%) are
  bidirectional at the locked operating point; 138,705 (20.97%) enter these
  sorted classes.
- Direction reverses sharply along DI. The fraction of unidirectionally sorted
  SNPs assigned to *F. aquilonia* rises from about 0.14--0.19 in the lowest DI
  deciles to about 0.73--0.77 in the highest deciles.
- DI is nearly orthogonal to recombination (Spearman rho = -0.029) and parental
  MAF (-0.018). It is associated mainly with lower within-parent diversity
  (rho = -0.464), and more weakly with higher between-parent allelic difference
  (rho = 0.134).
- The apparent concentration of sorted SNPs in the lowest-recombination decile
  weakens after LD reduction: the sorted fraction is 0.138 in the SNP sample
  but 0.039 among LD-reduced units in that decile. Across the other deciles the
  unit-level fraction is approximately 0.30--0.34.
- In the unit-level magnitude model, sorting increases with recombination
  (standardised coefficient 0.052) after accounting for DI.
- In the binomial-rule direction model, DI is dominant (coefficient 1.567),
  whereas recombination contributes a much smaller negative association
  (-0.091); parental MAF and unit size have still smaller effects.
- Among 13,191 candidate-containing clusters that become unsorted after
  aggregation, only nine have consensus fidelity below 0.80. Loss of the SNP
  signal is therefore rarely explained by poor consensus fidelity.

## Current scripts and outputs

- `moduleA_sorting_phenomenon.R`: SNP sorting, threshold sweeps and DI-direction
  summaries.
- `moduleA_cluster_sorting.R`: sorting at stored eMLG consensuses.
- `moduleA_architecture.R`: marker-versus-unit summaries and magnitude/direction
  models.
- `moduleA_fig_sorting_manhattan.R`: genomic sorting landscape; currently
  requires correction before use (see below).
- `moduleA_clusters.rds` is candidate-conditional, not a genome-wide cluster
  table. `moduleA_cluster_sorting.rds` contains the stored-eMLG analysis.

## Documentation status

- `supplementary_methods_moduleA.tex/pdf`: Methods source; primary thresholds
  and the binomial direction rule have been updated.
- `moduleA_architecture_summary.tex/pdf`: earlier main-text-style draft; values
  and provenance labels are stale and should be replaced, not incrementally
  polished.
- `moduleA_supp_th_sensitivity.tex/pdf`: earlier supplementary draft; its first
  two threshold summaries can be regenerated from current outputs, but the
  direction-sweep and Manhattan panels must wait for the code correction below.

## Items to resolve before freezing results documents

- `moduleA_architecture.R` recomputes its threshold direction sweep with the
  superseded component rule although the primary model uses the binomial rule.
- `moduleA_fig_sorting_manhattan.R` also reconstructs sorting with the component
  rule and labels the old near-fixation floor. Its current figure should not be
  cited as primary evidence.
- `moduleA_snp.rds` is named as a downstream input, but its save call is
  disabled; a clean run therefore does not currently supply the object expected
  by the Manhattan script.
- Several script headers still state the old tau = 0.5, phi = 0.90 and/or cM1
  settings even though the executable constants use the locked values.
- Confirm how ambiguous magnitude-passing units will be tabulated in the final
  reporting denominator.

## Proposed document split

- **Methods:** orientation; near-fixation; magnitude gate; binomial direction
  test; parental-MAF gate; SNP/eMLG aggregation; marker-versus-unit comparison;
  architecture models.
- **Headline Results/Discussion:** three claims only: sorting direction reverses
  with DI; DI is largely independent of recombination; LD-aware summaries show
  that low-recombination SNP peaks are strongly affected by redundancy.
- **Supplementary Results:** prevalence and threshold sensitivity; DI deciles;
  candidate-cluster aggregation and fidelity; full genomic architecture;
  coefficient tables; corrected direction sensitivity and genomic landscape.
