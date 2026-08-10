# SNPRelate LD-pruning vs LD-complexity reduction

*Methods-comparison note (Module 0). Worked example: the Chr26 centromeric block.*

## Question

Standard LD-pruning (e.g. `SNPRelate::snpgdsLDpruning`) thins a marker set to an
approximately independent subset by discarding SNPs in LD. Our LD-complexity
reduction (LDCR) instead **clusters** correlated markers into LD-reduced units,
keeps every marker in a cluster, and produces an expected multi-locus genotype
(eMLG) consensus per cluster. How do the two compare where it matters most — a
strong linkage block?

## Setup

We used the **Chr26 centromeric block** (3.40–5.72 Mb, 2.32 Mb, **3,844 markers**,
median `ld_w` = 0.37), the most pronounced low-recombination block in the genome.
SNPRelate `snpgdsLDpruning` was run on the block's genotypes (164 hybrids) at
several r² thresholds (`method="corr"`, whole-block sliding window unless noted).
LDCR is the canonical cM05 clustering (`eMLG_5loci_0025_cM05.rds`); its "units"
are the clusters the block's markers belong to (represented by the cluster
representative or the eMLG consensus). Both are compared on the same block.

## Results

### 1. Degree of reduction — comparable

Both collapse the block by ~90%. LDCR's 473 clusters sit between SNPRelate at
r²≤0.2 and r²≤0.5, i.e. it behaves like LD-pruning at **r² ≈ 0.3–0.4**
(Fig. `snprelate_vs_ldcr_units.png`).

| method (block markers → units) | units | % of 3,844 |
|---|---|---|
| SNPRelate r²≤0.1 | 258 | 6.7% |
| SNPRelate r²≤0.2 | 333 | 8.7% |
| **LD-complexity reduction (cM05)** | **473** | **12.3%** |
| SNPRelate r²≤0.5 | 561 | 14.6% |
| SNPRelate r²≤0.8 | 925 | 24.1% |

LDCR additionally resolves the block's structure: one dominant cluster of **900
markers** (the core linkage block, summarised by a single eMLG) plus a tail of
311 singletons and ~160 small clusters — not a uniform thinning.

### 2. Residual LD — comparable (slightly favours SNPRelate, by design)

Pairwise r² among the retained units (r²≤0.2 pruning). All three strip ~99% of
the LD; differences are confined to the tail.

| method | units | median r² | mean r² | max r² | % pairs r²>0.2 |
|---|---|---|---|---|---|
| SNPRelate, 500 kb window (default) | 374 | 0.0041 | 0.0104 | 0.961 | 0.1% |
| SNPRelate, whole-block window | 333 | 0.0039 | 0.0090 | 0.199 | 0.0% |
| LD-complexity reduction (representatives) | 473 | 0.0041 | 0.0115 | 1.000 | 0.5% |

Two points. First, SNPRelate at its **default 500 kb window** leaves a long-range
residual (max r² = 0.96): the window is narrower than the 2.32 Mb block, so it
never compares the far-apart-but-linked SNPs — a common, silent failure of naive
LD-pruning in low-recombination regions. Second, LDCR retains a slightly larger
tail (a few perfectly-linked pairs) **by design**: Stage 2 merges only within
0.5 cM, so two markers in perfect LD but >0.5 cM apart are kept as *separate
units*, because long-range LD in a young hybrid population is treated as a
between-unit relationship (structure/admixture/selection), not a physical block
to collapse. On raw residual LD the methods are therefore comparable.

### 3. Reproducibility — LDCR wins

SNPRelate starts from a random SNP and greedily prunes, so its output is
seed-dependent. Across 5 random starts (whole-block, r²≤0.2):

- retained counts: 333, 333, 346, 335, 332 — **mean pairwise Jaccard 0.886**
  (~11% of the retained set changes run to run).
- LDCR is **deterministic** (473 clusters; Jaccard = 1 by construction).

Different runs/analysts get the *same* LDCR units, but *different* SNPRelate sets.

### 4. Coverage / retention of informative markers — LDCR wins

SNPRelate keeps **333/3,844 (9%)** and discards the other 91%; LDCR assigns every
marker to a cluster (**100% represented**). Critically, of the **371 top-decile
diagnostic-index markers** in the block, SNPRelate retains **0** — it prunes every
diagnostic marker away, keeping arbitrary LD-tag neighbours instead. LDCR keeps
all of them, retrievable in their clusters.

*(The diagnostic signal is still tagged by the retained neighbours via LD; the
distinction is that LDCR keeps the diagnostic markers themselves, which matters
whenever the loci — not just their LD footprint — are the object of analysis.)*

### 5. Unit completeness — LDCR wins

A pruned unit is one raw SNP and carries its own missingness/genotyping noise; an
eMLG consensus pools the whole cluster.

- raw block SNPs: **5.0% missing** per unit
- eMLG consensus: **1.0% missing** (~5× lower; 38 block clusters have a consensus)

### 6. Interpretability — LDCR wins (qualitative)

LDCR returns interpretable objects — cluster membership, size, and physical/
genetic extent — so the block can be *seen* and audited (the
`Chr26_stage1_vs_combined` figures) and eMLGs built from it. SNPRelate returns
only a list of retained SNP IDs.

## Summary

| axis | winner |
|---|---|
| degree of reduction | comparable (~90% both) |
| residual LD | comparable (SNPRelate marginally lower, by our 0.5 cM design choice) |
| reproducibility | **LDCR** (deterministic vs Jaccard 0.886) |
| coverage of informative markers | **LDCR** (100% vs 9%; 0/371 diagnostic markers kept by SNPRelate) |
| unit completeness | **LDCR** (1.0% vs 5.0% missing) |
| interpretability / eMLG | **LDCR** |

On the metric SNPRelate optimises for — low pairwise residual LD — the two are
essentially equivalent. On every other axis relevant to a downstream ancestry
analysis (a stable, complete, interpretable set of units that keeps the diagnostic
markers and yields eMLG consensuses), LD-complexity reduction is the stronger tool.

---

*Reproduced by the scratch script `chr26_snprelate_vs_ldcr.R` /
`chr26_residual_ld.R` / `chr26_other_comparisons.R` (SNPRelate + the cM05
clustering); numbers above are for the Chr26 centromeric block.*
