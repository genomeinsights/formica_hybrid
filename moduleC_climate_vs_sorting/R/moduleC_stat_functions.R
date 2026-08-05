## =========================================================
## MODULE C (revised) -- shared per-covariate statistic functions
## =========================================================
## The WHOLE POINT of the structured-null calibration is that the observed climate
## axes and every one of the 10,000 Omega-structured null covariates are reduced to
## a genome-wide statistic by IDENTICAL code on the IDENTICAL eMLG universe. This
## file is the single definition of that code; it is sourced by both the
## null-regeneration script (per null covariate) and the analysis script (for the
## observed PC1/PC2). Do not duplicate the statistic logic anywhere else.
##
## Unit of analysis: the eMLG (LD-cluster consensus genotype). N = 32,840.
## Association strength for one covariate = the per-eMLG BayPass BF(dB) vector, in
## the fixed BayPass row order (== eMLG_group_order.txt == annotation table order).
##
## Because different covariates give different genome-wide BF distributions, the
## PRIMARY score is the WITHIN-COVARIATE rank / percentile of BF, not raw BF.
##
## SORTING (primary): directional ancestry sorting is biologically defined
## CONDITIONAL ON PARENTAL DIFFERENTIATION (Module A/C convention). The primary
## sorting statistic therefore contrasts directionally-sorted vs non-directional
## eMLGs WITHIN differentiated == TRUE only (`sort_gap_differentiated`). The
## all-eMLG contrast (`sort_gap_all`) is retained only as a supplementary
## sensitivity.
## SORTING MAGNITUDE (continuous): `prop_fixed` (the degree of fixation), NOT the
## signed `uni_score`. `uni_score` is signed sorting ORIENTATION and is kept only
## as a supplementary "sorting orientation" statistic.

## ---------------------------------------------------------------------------
## SORTING-THRESHOLD (tau) SERIES. Only the DIRECTIONAL classification depends on
## the fixation-fraction threshold tau (Module A's binom rule); the continuous
## annotations DI, recomb, prop_fixed (magnitude), uni_score (orientation) and the
## `differentiated` gate are tau-INDEPENDENT (verified: identical across the stamped
## moduleA_cluster_sorting_tau05/06/08.rds). The annotation table therefore carries
## ONE directional column per tau (directional_tau05/06/08) plus a `directional`
## alias == the primary (tau06); `prepare_annotation_ranks(ann, dir_col=...)` selects
## which one drives the directional-contrast statistics. All tau share the same
## DI/recomb/magnitude/orientation ranks, so only sort_gap_* / bf_gap_* / *_pdir_diff
## differ across the series.
MODULEC_TAU_SERIES  <- c(0.5, 0.6, 0.8)
MODULEC_TAU_PRIMARY <- 0.6
tauC_stamp      <- function(tau) sprintf("tau%02d", round(tau * 10))   # 0.6 -> "tau06"
dir_col_for_tau <- function(tau) paste0("directional_", tauC_stamp(tau))

## ---------------------------------------------------------------------------
## MINIMUM-CLUSTER-SIZE (min_n_loci) SERIES. Requiring >= min markers per eMLG
## targets the high-local-LD clusters expected under recent selection. Because the
## eMLG BayPass runs use a FIXED (LD-pruned) Omega passed via -omegafile, each
## cluster's BF is computed conditionally on Omega and is INVARIANT to the size
## threshold; a larger min is therefore a pure ROW-SUBSET of the same per-eMLG BF
## vectors, re-reduced over the smaller universe (ranks/percentiles recomputed on
## the subset). The null BF matrix is generated once over the primary (>=5) universe
## and reduced against every (min, tau) cell -- min=5 is the full universe, min=10 a
## strict subset (identical consensus genotypes). The genome-wide calibration is thus
## reported over a min x tau grid without any extra BayPass.
MODULEC_MIN_SERIES  <- c(5L, 10L)
MODULEC_MIN_PRIMARY <- 5L
minC_stamp <- function(m)      sprintf("min%02d", as.integer(m))       # 10 -> "min10"
cell_key   <- function(m, tau) paste0(minC_stamp(m), "_", tauC_stamp(tau))  # -> "min05_tau06"

## Precompute the fixed annotation ranks/masks once (constant across covariates).
## `ann` must be ordered to match the BF vector row order and carry:
##   DI, recomb, prop_fixed, uni_score (numeric); <dir_col> (0/1); differentiated (logical).
## `dir_col` names the directional column to use (default "directional" == primary tau);
## pass dir_col_for_tau(tau) to score a specific threshold in the series.
prepare_annotation_ranks <- function(ann, dir_col = "directional") {
  need <- c("DI", "recomb", "prop_fixed", "uni_score", dir_col, "differentiated")
  stopifnot(all(need %in% names(ann)))
  mag_ok <- which(!is.na(ann$prop_fixed))       # sorting magnitude complete-case subset
  ori_ok <- which(!is.na(ann$uni_score))        # sorting orientation complete-case subset
  list(
    DI          = ann$DI,
    recomb      = ann$recomb,
    directional = ann[[dir_col]] == 1L,
    dir_col     = dir_col,
    differ      = ann$differentiated == TRUE,
    rDI         = rank(ann$DI,     ties.method = "average"),
    rRec        = rank(ann$recomb, ties.method = "average"),
    ## sorting magnitude (prop_fixed): fixed complete-case subset, same for all covariates
    mag_ok      = mag_ok,
    mag_vals    = ann$prop_fixed[mag_ok],
    rMag        = rank(ann$prop_fixed[mag_ok], ties.method = "average"),
    ## sorting orientation (uni_score, signed): supplementary
    ori_ok      = ori_ok,
    ori_vals    = ann$uni_score[ori_ok],
    rOri        = rank(ann$uni_score[ori_ok], ties.method = "average"),
    N           = nrow(ann)
  )
}

## The top-fraction thresholds for the (secondary) rank-threshold sensitivity.
TOP_FRACS <- c(0.001, 0.005, 0.01)

## short tag for a fraction, e.g. 0.001 -> "top0001" (used for stat names + labels)
frac_tag <- function(f) sub("\\.", "", sprintf("top%0.3f", f))

## ---------------------------------------------------------------------------
## Primary + threshold statistics for a SINGLE covariate's genome-wide BF vector.
##   bf : numeric length N, aligned to `A` (output of prepare_annotation_ranks).
## Returns a NAMED numeric vector (same names every call -> stackable to a matrix).
## ---------------------------------------------------------------------------
compute_covariate_stats <- function(bf, A) {
  N <- A$N
  stopifnot(length(bf) == N)
  r   <- rank(bf, ties.method = "average")
  pct <- (r - 0.5) / N                       # within-covariate BF percentile in (0,1)
  dir <- A$directional; dd <- A$differ

  ## --- primary rank-based association statistics ---
  rho_DI  <- cor(r, A$rDI)                    # Spearman(BF, DI)  (ranks -> Pearson)
  rho_rec <- cor(r, A$rRec)                   # Spearman(BF, recombination)
  ## PRIMARY sorting: percentile gap, directional vs non-directional, DIFFERENTIATED only
  sort_gap_differentiated <- mean(pct[dd & dir]) - mean(pct[dd & !dir])
  ## supplementary sorting: same contrast over ALL eMLGs
  sort_gap_all <- mean(pct[dir]) - mean(pct[!dir])
  ## continuous sorting magnitude (prop_fixed), complete-case subset
  rho_sort_magnitude   <- cor(rank(bf[A$mag_ok]), A$rMag)
  ## signed sorting orientation (uni_score), supplementary, complete-case subset
  rho_sort_orientation <- cor(rank(bf[A$ori_ok]), A$rOri)

  ## raw-BF sensitivity (captured here: the null BF matrix is not retained and
  ## cannot be revisited).
  pear_DI              <- cor(bf, A$DI)
  pear_rec             <- cor(bf, A$recomb)
  pear_sort_magnitude  <- cor(bf[A$mag_ok], A$mag_vals)
  bf_gap_differentiated <- mean(bf[dd & dir]) - mean(bf[dd & !dir])
  bf_gap_all            <- mean(bf[dir]) - mean(bf[!dir])

  out <- c(rho_DI = rho_DI, rho_rec = rho_rec,
           sort_gap_differentiated = sort_gap_differentiated,
           sort_gap_all = sort_gap_all,
           rho_sort_magnitude = rho_sort_magnitude,
           rho_sort_orientation = rho_sort_orientation,
           pear_DI = pear_DI, pear_rec = pear_rec,
           pear_sort_magnitude = pear_sort_magnitude,
           bf_gap_differentiated = bf_gap_differentiated, bf_gap_all = bf_gap_all)

  ## --- threshold sensitivity: top fraction of eMLGs BY RANK within covariate ---
  ## per fraction: median DI, median recomb, proportion DIFFERENTIATED, and the
  ## proportion DIRECTIONAL AMONG DIFFERENTIATED (so the "is differentiated" and
  ## "is directional | differentiated" effects are separable). If a fraction
  ## contains no differentiated eMLG, pdir_diff is NA with a validation warning.
  ord <- order(bf, decreasing = TRUE)        # highest BF first
  for (f in TOP_FRACS) {
    k   <- max(1L, ceiling(f * N))
    idx <- ord[seq_len(k)]
    tag <- frac_tag(f)
    nd  <- sum(dd[idx])
    if (nd == 0L)
      warning(sprintf("threshold %s: no differentiated eMLG in top fraction -> pdir_diff = NA", tag))
    out[paste0(tag, "_DI_med")]    <- median(A$DI[idx])
    out[paste0(tag, "_rec_med")]   <- median(A$recomb[idx])
    out[paste0(tag, "_pdiff")]     <- mean(dd[idx])                       # prop. differentiated
    out[paste0(tag, "_pdir_diff")] <- if (nd == 0L) NA_real_ else mean(dir[idx][dd[idx]])
  }
  out
}

## Names produced by compute_covariate_stats(), in order (for pre-allocating).
covariate_stat_names <- function() {
  nm <- c("rho_DI", "rho_rec", "sort_gap_differentiated", "sort_gap_all",
          "rho_sort_magnitude", "rho_sort_orientation",
          "pear_DI", "pear_rec", "pear_sort_magnitude",
          "bf_gap_differentiated", "bf_gap_all")
  for (f in TOP_FRACS)
    nm <- c(nm, paste0(frac_tag(f), c("_DI_med", "_rec_med", "_pdiff", "_pdir_diff")))
  nm
}

## The six PRIMARY tests that form the FDR family (PC1/PC2 x DI/recomb/sorting).
## Sorting = differentiated-only directional contrast.
PRIMARY_STATS <- c(DI = "rho_DI", recombination = "rho_rec",
                   sorting = "sort_gap_differentiated")
