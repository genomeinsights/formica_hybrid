## =========================================================
## Cross-chromosome paralogy / duplication filter for the Module D two-locus screens.
## =========================================================
## A genetic-distance rule (LINK_CM) cannot catch cross-chromosome duplicates: two
## clusters on different chromosomes are at infinite recombinational distance yet can
## be genotype-near-identical when the SAME reads map to two assembly locations. Those
## pairs are the loudest inter-chromosomal "LD" and are technical, not epistasis
## (confirmed in the EMMAX arm: the |within-pop r|>0.9 hits are all cross-chromosome,
## up to concordance = 1.000). Distance filtering is orthogonal to this, so both arms
## also apply the signal-based filter below.
##
## DISCRIMINANT = individual-level WITHIN-population correlation of the two clusters'
## consensus dosages. Within-pop removes the among-population / ancestry-axis component
## (which is not what paralogy is about). Real long-range coupling is moderate
## (|r| <~ 0.6 in these data); a duplication sits at |r| ~ 1 (a flipped/complementary
## duplicate gives r ~ -1, which |r| still catches). Genotype concordance and
## per-cluster excess heterozygosity (a duplicated region shows inflated Ho) are
## reported alongside as corroboration -- the flag itself is deliberately simple
## (|within-pop r| > thr) so it can be re-thresholded from the stored values.
## NB |within-pop r| ~ 1 is a paralogy *candidate*: an extremely strong genuine
## within-deme incompatibility could also reach it, so treat it as a screen to be
## confirmed (concordance ~ 1 and/or excess Ho), not a proof.

suppressMessages(library(parallel))

## median |within-population r| of two clusters' dosages
moduleD_within_pop_r <- function(a, b, GT, pop_idx, min_n = 4) {
  x <- GT[, a]; y <- GT[, b]
  v <- vapply(pop_idx, function(ix) {
    xi <- x[ix]; yi <- y[ix]; ok <- is.finite(xi) & is.finite(yi)
    if (sum(ok) < min_n) return(NA_real_)
    if (sd(xi[ok]) == 0 || sd(yi[ok]) == 0) return(NA_real_)
    cor(xi[ok], yi[ok])
  }, numeric(1))
  median(abs(v), na.rm = TRUE)
}

## fraction of individuals with identical hard-called dosage (positive-duplicate signal)
moduleD_concordance <- function(a, b, GT) {
  x <- round(GT[, a]); y <- round(GT[, b]); ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  mean(x[ok] == y[ok])
}

## per-cluster mean observed heterozygosity from raw member genotypes (duplication ->
## excess Ho). marker_Ho = colMeans(GT_raw == 1). Returns a named vector over `ids`.
moduleD_cluster_het <- function(groups, ids, marker_Ho) {
  ml <- groups[.(ids), on = "group_id", .(marker = unlist(members)), by = group_id]
  ml[, ho := marker_Ho[marker]]
  ml[, .(het = mean(ho, na.rm = TRUE)), by = group_id][, setNames(het, group_id)]
}

## annotate a pairs table with the paralogy diagnostics + flag. `GT` is the continuous
## eMLG consensus (individuals x clusters); a_col/b_col name the two cluster-id columns.
flag_paralogy <- function(dt, a_col, b_col, GT, pops, het_of = NULL, thr = 0.9,
                          cores = 1, min_n = 4) {
  dt <- copy(dt); a <- dt[[a_col]]; b <- dt[[b_col]]
  pop_idx <- split(seq_along(pops), as.character(pops))
  idx <- seq_len(nrow(dt))
  dt[, within_pop_r := unlist(mclapply(idx, function(i) moduleD_within_pop_r(a[i], b[i], GT, pop_idx, min_n), mc.cores = cores))]
  dt[, concordance := unlist(mclapply(idx, function(i) moduleD_concordance(a[i], b[i], GT), mc.cores = cores))]
  if (!is.null(het_of)) dt[, `:=`(het1 = het_of[a], het2 = het_of[b])]
  dt[, paralog := is.finite(within_pop_r) & within_pop_r > thr]
  dt[]
}
