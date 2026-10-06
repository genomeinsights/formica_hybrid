# module_allele_specific_sorting

> **Status (2026-10-06): run and audited.** Empirical analyses complete
> (scripts 01–05, logs in `logs/`); revised after an external audit (unpruned-SNP
> sensitivity, matched-set inference with calipers, permutation baseline,
> long-range analysis dropped, wording). Comparisons with the neutral
> simulations wait for the corrected founder frequencies
> (`../sim_founder_fix/NOTE_for_Beatriz.md`).

Streamlined successor to the exploratory `module_population_partitioning`
(kept as the record). Question: is the strong differentiation among the 20
hybrid populations at ancestry-informative loci accompanied by correspondingly
strong LD, and do sorted loci share their sorting with nearby markers more than
equally differentiated unsorted loci?

## The argument in one formula

For two loci with population allele-frequency profiles F_i, F_j (oriented to
*F. aquilonia*), the among-population component of LD (Ohta 1982) is

    rST_ij = cov_pops(F_i, F_j) / sqrt(Fbar_i(1-Fbar_i) Fbar_j(1-Fbar_j)) = c_ij · sqrt(G_i G_j)

G (Nei-type among-population differentiation) sets a **ceiling** G_iG_j; the
profile correlation c sets the **realised fraction** c². High F_ST produces high
among-population LD only where loci partition populations the same way.

## Pipeline (run from the formica_hybrid repo root, in order)

| script | role | output |
|---|---|---|
| `R/00_utils.R` | provenance-guarded loading, oriented genotypes, standardised matrices, chromosome bootstrap draws | — |
| `R/01_pair_stats.R` | all 9.6 M within- and 206.8 M cross-chromosome unit pairs: concordance (raw/residualised), rST², ceiling, within-pop LD | `data/01_*` (`01_pairs_le2Mb.rds` git-ignored, ~150 MB) |
| `R/02_label_perm_baseline.R [N_PERM] [N_CORES]` | population-label permutation baseline (500) for the realised fraction and concordance per distance bin; exact algebraic shortcut for cross-chromosome pairs, validated against 01 | `data/02_baseline.rds` |
| `R/03_distance_decay.R` | **headline**: ceiling vs among- vs within-population LD along genetic distance (+ physical as sensitivity), chromosome-block bootstrap, permutation band | `Figures/03_decay_cM.*`, `03_decay_bp.*` |
| `R/04_sorted_neighbourhood.R` | sorted loci vs matched unsorted loci (K ≤ 5, 1-SD caliper on F_ST, recombination, parental differentiation, density): ancestry-profile similarity and both-segregating within-pop LD to neighbours, for LD-reduced units AND all unpruned SNPs; own-unit size/span; matched-set differences, anchor-chromosome bootstrap with sets intact | `Figures/04_neighbourhood.*`, `04_balance.*` |
| `R/05_doc_tables.R` | regenerates `doc_manuscript/tables/` and `figures/` from the outputs | — |

## Results (2026-10-06)

- **Ceiling vs realised**: ceiling G_iG_j ≈ 0.09 at every distance; realised
  fraction 0.217 at 0.001–0.01 cM → 0.065 at 1–5 cM → 0.058 beyond 5 cM and 0.0571
  for unlinked pairs, against a permutation baseline of 0.0526. The unlinked
  excess (0.004) is the shared ancestry gradient (residualised unlinked
  concordance ≈ 0). Pairs < 0.001 cM (cold-spot pairs) are less concordant than
  0.001–0.01 cM; monotonic on the physical scale.
- **Sorted neighbourhoods**: 1,552 / 1,577 sorted loci matched (all |SMD| ≤ 0.04).
  Similarity to neighbours within 5 kb is LOWER around sorted loci (≈ −0.06 to −0.07,
  both marker sets), equal beyond 20 kb; both-segregating within-pop LD equal or
  slightly lower, never higher; own LD units not larger (span shorter). Pooled
  r_w is diluted at sorted loci — never use it for this comparison.
- **Unmatched**: 25 sorted loci with no unsorted counterpart are giant LD blocks
  (mean 121 SNPs, 907 kb; multi-Mb blocks on Chr16, 1, 17, 26) — haplotype-scale
  sorting the matched test cannot assess.
- **Interpretation**: consistent with locus-specific sorting outside a few large
  blocks; not yet shown to differ from neutral admixture (simulations pending).

## Dropped (deliberately)

The long-range (cross-chromosome tail / BDMI) analysis — removed after audit
(tested only same-parent ancestry coupling; null anti-conservative for tail
counts; separate question). Last version: commit a73ddc5. Also dropped:
F_ST-tertile figure, unadjusted and pooled within-pop LD, the
aquilonia/polyctena split, the large distance table.

## Universe and conventions

20,807 DI25 rho05 units (`module_di25_rho05`, cM5 cap, `min_r2_rho=0.5`), 51,612
DI25 SNPs for the unpruned sensitivity, 20 hybrid populations, strictly
observed genotypes, orientation from the parental reference samples (asserted
to reproduce `Fmat`), sorting at τ=0.6/φ=0.85/binom/α=0.05, W&C F_ST,
leave-one-chromosome-out ancestry residuals from `pp_residualize_ancestry.R`. No
parental-MAF gate.
