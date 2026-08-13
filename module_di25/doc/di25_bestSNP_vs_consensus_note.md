# Representative SNP vs eMLG consensus: fidelity and sorting

*module_di25 — high-DI (DI > −25) diagnostic markers, 5 cM LD clustering*

## Motivation

LD-complexity reduction summarises each LD block by an **eMLG consensus** (a
polarised average of the block's markers). Averaging is ideal for one value per
block, but some genotype signal is genuinely SNP-specific and can differ even
between strongly-linked markers; there the consensus dilutes it and a single
representative SNP is preferable. We asked (i) how faithfully the consensus is
reproduced by a single SNP, (ii) whether the block **fidelity score** predicts
this, and (iii) whether the choice changes a downstream result (ancestry
sorting).

## Data and units

The from-scratch DI25 clustering (`di25_ld_clustering.R`, 5 cM cap) reduces
51,612 diagnostic markers to LD units; **4,010** of these are eMLG blocks
(> 2 markers, `has_eMLG = TRUE`) where the consensus can differ from any single
member. All correlations are across the **165 hybrids** used for LD estimation;
the sorting analysis uses the 194 hybrids + parents of `di25_sorting.R`.

## Methods

1. **Consensus vs representative.** Per block, Pearson *r* across hybrids between
   the eMLG consensus and its clustering **representative** SNP (chosen for
   cluster centrality).
2. **Best single-SNP proxy.** Per block, the maximum `|r|` between the consensus
   and *any* member SNP, and the gain over the representative. Orientation is
   free (consensus polarity is arbitrary), so `|r|` is used.
3. **Filled best-SNP genotype.** The best-`|r|` member's genotype with its
   *missing* calls filled from the (orientation-matched, rounded) consensus;
   observed calls are never overwritten. Packaged as
   `LDscnR::eMLG_best_snp()`.
4. **Sorting comparison.** Module-A sorting (`sort_rule = "binom"`, φ = 0.85,
   τ ∈ {0.5,…,0.8}, DI/pMAF ungated) on the same 4,010 blocks, once with the
   consensus unit and once with the filled best-SNP unit.
5. **Discriminating check.** For units the best-SNP sorts but the consensus does
   not, whether they are broadly supported (`n_obs`, `n_fixed` across the 20
   hybrid populations) or thin near-fixation calls.

## Results

**1. The consensus is well but imperfectly reproduced by one SNP.** Median
*r*(consensus, representative) = **0.87**; ~6% of blocks fall below 0.7, rising
to ~10–12% among 6–20-marker blocks. Taking the *best* member instead: median
`|r|` = **0.92**, and **no block** falls below **0.71** — a single SNP can always
stand in for the consensus reasonably well. The representative is the best member
in only **45%** of blocks; the median gain from switching is ~0, but reaches
0.17 in the tail (concentrated in mid-sized blocks).

**2. Fidelity score predicts single-SNP faithfulness.** Best-member `|r|`
correlates with the eMLG fidelity score at **Spearman 0.90** (Pearson 0.87),
and this is not a size artifact (partial *r* controlling for log₁₀ n_loci =
0.88). High-score blocks (≥ 0.95) are near-perfectly captured by one SNP
(best `|r|` ≥ 0.98); low-score blocks (< 0.85) are where the consensus does real
work (best `|r|` ~0.85). The score is therefore a principled selector for when to
fall back to a representative SNP.

**3. Filling.** The best-SNP matrix reduces missingness **2.8% → 0.3%** (residual
only where the consensus was itself missing); 80% of blocks become complete; 2%
required orientation flipping.

**4. The best-SNP recovers markedly more sorting than the consensus.** At τ = 0.6
the best-SNP sorts **5.5%** of blocks vs **3.6%** for the consensus (≈1.5×; ≈1.9×
at τ = 0.5). The disagreement is one-directional — 82 units sorted only by the
best-SNP vs 5 only by the consensus at τ = 0.6 — and sorting *direction* is
essentially never in conflict (≤ 7 units), with the best-SNP a few points more
*aquilonia*-leaning.

**5. The extra units are real, not thin.** Of the 82 best-SNP-only sorted units
(τ = 0.6): every one has data in all 20 populations (median `n_obs` = 20), a
median of **12/20** populations near-fixed (vs 15 for both-sorted), **0%** thin
(`n_fixed` ≤ 3 or `n_obs` < 12), **100%** robust (`n_obs` ≥ 15 & `n_fixed` ≥ 9).
The recovery is **not** from the missing-data fill (`n_obs` identical between
representations for 100% of them) but from consensus **dilution**: for every one
of these units the best SNP is near-fixed in more populations than the consensus
(median **+4**), because averaging softens genotypes below the φ = 0.85
near-fixation floor. These are exactly the mid-fidelity blocks (best `|r|` 0.91,
score 0.88) flagged in results 1–2.

## Conclusion

The eMLG consensus is a faithful block summary on average, but for a non-trivial
minority of (mid-sized, mid-fidelity) blocks it erases population-broad,
SNP-specific signal — demonstrably so on ancestry sorting, where the filled
best-SNP surfaces ~50% more genuinely-supported sorted units at τ = 0.6. When the
signal of interest is SNP-specific, the consensus-optimal member SNP (with gaps
filled from the consensus) is the better unit, and the eMLG **fidelity score**
is a usable rule for deciding per block. This selection is now available as
`LDscnR::eMLG_best_snp()`.

## Artifacts

| script (`module_di25/R/`) | output |
|---|---|
| `di25_emlg_vs_representative.R` | `data/di25_emlg_vs_representative.rds`, `Figures/di25_emlg_vs_representative.png` |
| `di25_emlg_best_member.R` | `data/di25_emlg_best_member.rds`, `Figures/di25_emlg_best_member.png` |
| `di25_bestSNP_vs_score.R` | `Figures/di25_bestSNP_vs_score.png` |
| `di25_bestSNP_consensus_filled.R` | `data/di25_bestSNP_consensus_filled.rds` |
| `di25_sorting_bestSNP_vs_consensus.R` | `data/di25_sorting_bestSNP_vs_consensus.rds`, `Figures/di25_sorting_bestSNP_vs_consensus.png` |
| `di25_sorting_bestSNP_discriminate.R` | `data/di25_sorting_bestSNP_discriminate.rds`, `Figures/di25_sorting_bestSNP_discriminate.png` |

Packaged function: `LDscnR::eMLG_best_snp(result, GTs, fill = TRUE)` — returns
per-block `stats` (best_marker, best_abs_r, rep_abs_r, score, fill counts) and,
optionally, the consensus-filled best-SNP genotype matrix.
