# module_allele_specific_sorting

> **Status (2026-10-06): scripts written for review, NOT run.** Empirical-only
> part. The empirical-vs-simulation comparisons wait for the corrected
> simulations (founder-frequency bug, see `../sim_founder_fix/NOTE_for_Beatriz.md`).

Streamlined successor to the exploratory `module_population_partitioning`
(kept frozen as the record). Goal: show that ancestry sorting among the 20
hybrid populations is **allele-specific** — high F_ST at individual loci
without correspondingly high LD — and test whether any long-range
(unlinked) associations of the kind epistatic selection (BDMIs) predicts
exist.

## The argument in one formula

For two loci with population allele-frequency profiles F_i, F_j (oriented to
*F. aquilonia*), the among-population component of LD is

    rST_ij = cov_pops(F_i, F_j) / sqrt(Fbar_i(1-Fbar_i) Fbar_j(1-Fbar_j)) = conc_ij · sqrt(G_i G_j)

where G is the among-population F_ST (Nei-type) of each locus and conc is
the correlation of the two profiles. High F_ST raises the **ceiling**
sqrt(G_i G_j), but LD among populations materialises only to the extent that
the loci partition the populations the same way (conc). "High F_ST but not
high LD" = high ceiling, low realised fraction except over linkage distances.

## Pipeline (run from the formica_hybrid repo root, in order)

| script | question | key output |
|---|---|---|
| `R/00_utils.R` | shared loading (provenance-guarded), oriented genotypes, standardised matrices, bootstrap, BDMI intervals | — |
| `R/01_pair_stats.R` | all within- and cross-chromosome unit pairs: conc, conc_resid, rST², ceiling G_iG_j, within-pop r_w, r_w_adj | `data/01_*` |
| `R/02_distance_decay.R` | **A1 headline**: ceiling vs realised among-population LD vs within-population LD along bp / cM / unlinked, chromosome-block bootstrap; by pair-F_ST class | `Figures/02_fst_vs_ld_decay.png` |
| `R/03_sorted_anchor_decay.R` | **A2**: do sorted units drag their neighbourhood (sweep-like) or not (allele-specific)? Sorted anchors vs F_ST/recombination/parental-differentiation/density-matched unsorted controls | `Figures/03_sorted_anchor_decay.png` |
| `R/04_longrange_tail.R [B_NULL]` | **A3**: tails of cross-chromosome associations (among-pop residual concordance; within-pop hybrid-index-adjusted LD) vs a chromosome-wise permutation null; BDMI/sorted enrichment of tail endpoints; LOO-population robustness of top pairs | `Figures/04_longrange_tail.png` |

Recombination (A4) is not re-done: `module_population_partitioning/R/pp_recombination.R`
already shows concordance decays with cM and is higher in low-recombination
regions; script 02 adds the cM axis to the new statistics.

Cost notes: 01 and 04 compute all ~2.2e8 cross-chromosome pairs by matrix
cross-products in per-chromosome blocks (peak ~0.5 GB per block). 04 repeats
the scan for each of B_NULL permutations × 3 statistics (default 50; expect
tens of minutes; use mini1/mini2 for larger B).

## Universe and conventions (unchanged from module_population_partitioning)

20,807 DI25 rho05 units (`module_di25_rho05`, cM5 cap, `min_r2_rho=0.5`),
20 hybrid populations, best-SNP/representative genotypes (strictly observed),
orientation to *F. aquilonia* from the parental reference samples (asserted
to reproduce `Fmat` exactly), sorting at τ=0.6/φ=0.85/binom/α=0.05,
W&C F_ST, leave-one-chromosome-out ancestry residuals from
`pp_residualize_ancestry.R`. No parental-MAF gate (DI>−25 is the gate).

## New statistics vs the exploratory module

- **rST² and its ceiling G_iG_j** (Ohta-style among-population LD, normalised) —
  turns "concordance" into the LD statement the paper needs.
- **within-population LD r_w / r_w_adj** on the same units and pairs, so the
  within/among split is computed on one unit set. r_w_adj removes admixture LD
  from variation in individual hybrid index (LOCO), leaving locus-specific
  conspecific-ancestry association.
- **Tail tests with a chromosome-wise permutation null** instead of
  genome-wide means: epistasis predicts sparse specific unlinked pairs.

## Dropped from the exploratory module (deliberately)

Genotype-run proxy, parental-structure-by-DI (no metadata), the BayPass
candidate-locus follow-up (different unit universe), |r| and Euclidean
similarity, the population-label permutation floor. Geographic prediction,
individual influence and the Sielva-driven PC1 stay as supplement-level
statements in `module_population_partitioning`.

## Known limitations

- Mean imputation of the 2.7% missing genotypes after within-population
  centring shrinks r_w slightly toward 0 (same for observed and permuted data).
- Raw cross-chromosome concordance (conc_raw) cannot separate a demographic
  ancestry gradient from a polygenic, aquilonia-biased BDMI network; only the
  residualised and within-population statistics address pair-specific
  epistasis, and those lose any signal aligned with genome-wide ancestry.
- BDMI regions are nodes without partner information: only endpoint
  enrichment is testable, not whether a tail pair is a reported BDMI edge.
- The empirical null answers "more cross-chromosome association than
  chromosomes shuffled against each other", not "more than neutral
  admixture would produce"; the latter needs the corrected simulations run
  through `01`–`04` (the functions take genotypes + populations, so a
  simulation replicate can be passed through unchanged).
