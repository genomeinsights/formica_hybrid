## =========================================================================
## module_allele_specific_sorting -- 01: pair statistics for all unit pairs
##
## For every pair of the 20,807 DI25 rho05 units:
##   conc, conc_resid       among-population concordance (raw / ancestry-residualised)
##   rST2, ceil2            among-population LD^2 (= conc^2 G_i G_j) and its
##                          shared-partition ceiling G_i G_j
##   r_w, r_w_adj (and ^2)  within-population LD (pop-centred / + hybrid-index-residualised)
## (definitions: R/00_utils.R header)
##
## Outputs (module_allele_specific_sorting/data/):
##   01_agg_within.rds   per Chr x distance bin x pair-F_ST class sums and counts
##                       (bp and cM bins) -> chromosome-block bootstrap in 02
##   01_agg_cross.rds    per (ChrA, ChrB) x pair-F_ST class sums and counts for
##                       cross-chromosome pairs (the unlinked reference)
##   01_pairs_le2Mb.rds  pair-level table for within-chromosome pairs <= 2 Mb (03)
##   01_units.rds        unit table (+ cM, recombination rate, G_nei, F_ST class)
##
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/01_pair_stats.R
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")

PAIR_TABLE_MAXD <- 2e6

U  <- load_units(); u <- U$u
GG <- load_oriented_genotypes(U)
cat(sprintf("[01] %d units, %d hybrid individuals in %d populations\n",
            nrow(u), nrow(GG$G), length(unique(GG$pop))))

## pair F_ST class: tertile of min(F_ST_i, F_ST_j), tertile breaks from unit F_ST
fst_br <- quantile(u$FST, c(1/3, 2/3), na.rm = TRUE)
u[, fst_tert := cut(FST, c(-Inf, fst_br, Inf), labels = c("low", "mid", "high"))]
cat(sprintf("[01] unit F_ST tertile breaks: %.3f, %.3f\n", fst_br[1], fst_br[2]))

## standardised matrices: crossprod of any two gives the pairwise correlations
Za  <- among_pop_Z(U$Fmat)
Zr  <- among_pop_Z(U$Resid)
Zw  <- within_pop_Z(GG, u, adj = FALSE)
Zwa <- within_pop_Z(GG, u, adj = TRUE)
rm(GG); invisible(gc())

stat_names <- c("conc", "conc_resid", "rST2", "ceil2", "r_w", "r_w_adj", "r_w2", "r_w_adj2")

## pair-level stats for index vectors i (rows) x j (cols); returns a long table
pair_block <- function(i, j, upper_only) {
  A  <- crossprod(Za[, i, drop = FALSE],  Za[, j, drop = FALSE])
  R  <- crossprod(Zr[, i, drop = FALSE],  Zr[, j, drop = FALSE])
  W  <- crossprod(Zw[, i, drop = FALSE],  Zw[, j, drop = FALSE])
  WA <- crossprod(Zwa[, i, drop = FALSE], Zwa[, j, drop = FALSE])
  keep <- if (upper_only) upper.tri(A) else matrix(TRUE, nrow(A), ncol(A))
  ii <- i[row(A)[keep]]; jj <- j[col(A)[keep]]
  gg <- u$G_nei[ii] * u$G_nei[jj]
  data.table(i = ii, j = jj,
             conc = A[keep], conc_resid = R[keep],
             rST2 = A[keep]^2 * gg, ceil2 = gg,
             r_w = W[keep], r_w_adj = WA[keep], r_w2 = W[keep]^2, r_w_adj2 = WA[keep]^2,
             fst_class = pmin(as.integer(u$fst_tert[ii]), as.integer(u$fst_tert[jj])))
}

aggregate_stats <- function(dt, by) {
  dt[, c(lapply(.SD, function(x) sum(x, na.rm = TRUE)),
         lapply(.SD, function(x) sum(!is.na(x)))), by = by, .SDcols = stat_names] |>
    setnames(c(by, paste0("s_", stat_names), paste0("n_", stat_names)))
}

## ---- within-chromosome pairs -------------------------------------------------
chrs <- unique(u$Chr)
agg_bp <- agg_cm <- pairs_near <- vector("list", length(chrs))
n_zero_cM <- 0L
for (k in seq_along(chrs)) {
  idx <- u[Chr == chrs[k], idx]
  if (length(idx) < 2) next
  dt <- pair_block(idx, idx, upper_only = TRUE)
  dt[, Chr := chrs[k]]
  dt[, dist_bp := abs(u$Pos[j] - u$Pos[i])]
  dt[, dist_cM := abs(u$cM_pos[j] - u$cM_pos[i])]
  dt[, bp_bin := cut(dist_bp, BP_BREAKS, labels = BP_LABELS, right = FALSE)]
  dt[, cm_bin := cut(dist_cM, CM_BREAKS, labels = CM_LABELS, right = FALSE)]
  dt[dist_cM == 0, cm_bin := NA]                     # undefined genetic distance (00_utils.R)
  n_zero_cM <<- n_zero_cM + dt[dist_cM == 0, .N]
  agg_bp[[k]] <- aggregate_stats(dt, c("Chr", "bp_bin", "fst_class"))
  agg_cm[[k]] <- aggregate_stats(dt[!is.na(cm_bin)], c("Chr", "cm_bin", "fst_class"))
  pairs_near[[k]] <- dt[dist_bp <= PAIR_TABLE_MAXD,
                        .(Chr, i, j, dist_bp, dist_cM, conc, conc_resid, r_w, r_w_adj)]
  cat(sprintf("[01] %s: %d units, %d pairs (%d <= 2 Mb)\n",
              chrs[k], length(idx), nrow(dt), nrow(pairs_near[[k]])))
  rm(dt); invisible(gc())
}
agg_within <- list(bp = rbindlist(agg_bp), cm = rbindlist(agg_cm), fst_breaks = fst_br, n_zero_cM = n_zero_cM)
cat(sprintf("[01] within-chromosome pairs with zero map distance (excluded from cM bins only): %s\n",
            format(n_zero_cM, big.mark = ",")))

## ---- cross-chromosome pairs (unlinked reference) -----------------------------
agg_cross <- list()
for (a in seq_along(chrs)) for (b in seq_along(chrs)) if (a < b) {
  dt <- pair_block(u[Chr == chrs[a], idx], u[Chr == chrs[b], idx], upper_only = FALSE)
  dt[, `:=`(ChrA = chrs[a], ChrB = chrs[b])]
  agg_cross[[length(agg_cross) + 1L]] <- aggregate_stats(dt, c("ChrA", "ChrB", "fst_class"))
}
agg_cross <- rbindlist(agg_cross)
cat(sprintf("[01] cross-chromosome pairs: %s\n", format(sum(agg_cross$n_conc), big.mark = ",")))

saveRDS(agg_within, file.path(OUT_DATA, "01_agg_within.rds"))
saveRDS(agg_cross,  file.path(OUT_DATA, "01_agg_cross.rds"))
saveRDS(rbindlist(pairs_near), file.path(OUT_DATA, "01_pairs_le2Mb.rds"))
saveRDS(u, file.path(OUT_DATA, "01_units.rds"))
cat("[01] done\n")
