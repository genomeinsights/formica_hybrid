# Structured vs. paired-unstructured null covariate design

Status: design, not yet executed against BayPass. This document is the
specification that `R/01_generate_null_covariates.R` implements; read it
before reading the script.

## Goal

The full-SNP null comparison needs two contrasting null covariate designs,
run through the identical BayPass covariate-mode pipeline, so that any
difference in local-score behaviour between them isolates one thing: whether
the covariate's population-to-population pattern of variation is tied to the
population relatedness/covariance structure (Ω) or not.

- **Structured**: a null covariate whose among-population variance pattern
  is drawn to respect Ω (the same "Ω-eigenvector MVN" construction used
  everywhere else in this pipeline -- see Provenance below).
- **Unstructured**: a *paired* counterpart to each structured draw, built by
  randomly permuting that draw's 19 realized population values across the 19
  populations.

Pairing (rather than drawing the unstructured set completely independently,
e.g. from iid N(0,1)) is deliberate: a permutation of a structured draw has,
by construction, the exact same mean, SD, and empirical value distribution
as its structured counterpart (both are already z-scored: mean 0, SD 1,
identical value multiset). The only thing a within-pair permutation changes
is *which population gets which value* -- i.e. whether the covariate's
pattern correlates with Ω. This removes value-distribution differences as a
possible confound when comparing structured vs. unstructured local-score
behaviour; an iid-redraw design would not offer that guarantee as cleanly.

## Provenance of the "structured" mechanism (reused formula, fresh code)

The Ω-eigenvector MVN draw is the same mechanism used throughout
`module_manuscript_rho05` and its Module C (e.g.
`module_manuscript_rho05/R/moduleB_stage1_S1units_null.R`,
`module_manuscript_rho05/R/moduleB_ancestry_climate_mitotype_CORRECTED.R`):

```
eg    <- eigen(Omega, symmetric = TRUE)
vals  <- pmax(eg$values, 0)
draw  <- as.numeric(scale(eg$vectors %*% (sqrt(vals) * rnorm(P))))
```

This module reimplements that formula from scratch (see the script), with
its own validation, rather than sourcing or copying the existing scripts --
per the instruction to reimplement cleanly and not inherit assumptions.
**No unstructured/permuted counterpart exists anywhere in the codebase**
(confirmed by search); the permutation design below is new.

## Inputs (authoritative, from `module_manuscript_rho05`)

| Quantity | Source | Value |
|---|---|---|
| Ω (population covariance matrix) | `module_manuscript_rho05/baypass_stage1/aland_excluded/omega_mat_omega.out` | 19x19, BayPass `omega_mat_omega.out` format |
| Population order | `unique(sample_data[Population != "Aland"]$Population)` from `data/hybrids_only_maf005.Rdata` | 19 populations, Aland already excluded |
| Population sizes (cross-check only) | `module_manuscript_rho05/baypass_stage1/aland_excluded/u_DIEM.size` | 19 integers |
| Real mitoC2 contrast (for the mitoC2 split size) | `module_manuscript_rho05/baypass_stage1/aland_excluded/u.mito_contrast` | 7 populations = +1 (Faquilonia-like), 12 = -1 (Fpolyctena-like) |

Ω is verified byte-identical (MD5) between `aland_excluded/` and
`aland_excluded_S1units/` -- it is a single population-level matrix,
independent of SNP/unit resolution, so it is valid for a full-SNP scan.

**Known fragility, checked explicitly (not assumed):** `omega_mat_omega.out`
carries no population names -- its row/column order is BayPass's own
positional convention (same order as the geno file's population columns,
which is `pop_order` by construction of the original input-writing script).
This exact kind of silent mis-order was previously caught and fixed in
`moduleB_ancestry_climate_mitotype_CORRECTED.R` ("AUDIT FIX Issue 3"). Here
that risk is mitigated by (a) using dimension/order cross-checks against
`u_DIEM.size` and `u.mito_contrast`, both of which are independently
confirmed to already be in `pop_order`, and (b) recording all of this in the
manifest so it is auditable rather than assumed.

## Construction

Let P = 19.

### Continuous nulls (10 structured + 10 paired unstructured)

For replicate `i` in 1..10, with prespecified, disjoint seed streams
(fixed before any pilot result is seen):

- `STRUCTURED_SEEDS  <- 20100 + 1:10`  (i.e. 20101 .. 20110)
- `PERMUTE_SEEDS     <- 20200 + 1:10`  (i.e. 20201 .. 20210)

```
set.seed(STRUCTURED_SEEDS[i])
raw_i         <- Omega_eigvecs %*% (sqrt(Omega_vals) * rnorm(P))
structured_i  <- as.numeric(scale(raw_i))          # mean 0, sd 1

set.seed(PERMUTE_SEEDS[i])
perm_i          <- sample.int(P)
unstructured_i  <- structured_i[perm_i]            # same value multiset, permuted assignment
```

### mitoC2 (contrast) nulls (10 structured + 10 paired unstructured)

Derived from the continuous nulls above by rank-thresholding to the real
group sizes (7 vs 12), exactly mirroring how the real mitoC2 contrast is a
population bipartition and how the precedent calibration null derived
contrast nulls from continuous draws -- reimplemented fresh here, at
10-replicate scale instead of a 10,000-draw calibration pool:

```
top7 <- function(x) { p <- rep(-1L, P); p[order(x, decreasing = TRUE)[1:7]] <- 1L; p }

structured_mitoC2_i    <- top7(structured_i)
unstructured_mitoC2_i  <- top7(unstructured_i)
```

Both are checked to have exactly 7 populations = +1 and 12 = -1.

## Validation checks (run by the generator script, recorded in the manifest)

1. Ω symmetry: report `max(abs(Omega - t(Omega)))` before forcing
   `(Omega + t(Omega)) / 2`.
2. Ω eigenvalues: report how many original eigenvalues were negative /
   clamped to 0 by `pmax(eg$values, 0)`.
3. Dimension/order cross-check: `nrow(Omega) == length(pop_order) ==
   length(u_DIEM.size) == 19`.
4. Each `structured_i`: `mean` == 0 and `sd` == 1 to floating tolerance.
5. Each `unstructured_i`: identical *sorted* value vector to `structured_i`
   (proves it's a true permutation, not a fresh draw) and confirmed *not*
   identical in original order (permutation actually moved something).
6. Each `structured_mitoC2_i` / `unstructured_mitoC2_i`: exactly 7 populations
   coded +1 and 12 coded -1.
7. Seed streams (`STRUCTURED_SEEDS`, `PERMUTE_SEEDS`) do not overlap with each
   other or across replicates.
8. Full reproducibility: script is deterministic given the fixed seeds --
   re-running produces byte-identical output (checked via checksum in the
   manifest).

## What this design deliberately excludes

- No use of `ld_w x statistic` weighting or any LD-weighted null.
- No circular permutation as a calibration mechanism.
- Population order is never re-derived from Ω itself (Ω carries no names) --
  it always comes from the independent, authoritative `pop_order` source
  above, consistent with how the real covariate/contrast files were built.
