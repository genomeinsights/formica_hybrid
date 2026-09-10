## =========================================================================
## module_population_partitioning -- 09: remaining AUDIT.md robustness items
## not covered by pp_pca_refined.R:
##   (a) leave-Sielva-out (and leave-any-one-population-out) recompute of the
##       headline genome-wide FST-vs-local-similarity statistic and the
##       distance-decay curve -- does the F1-like colony (elevated
##       heterozygosity, already flagged as the dominant PC1 driver in
##       pp_pca_refined.R) drive the main concordance results too, or only PC1?
##   (b) literal once-per-region treatment of the 3 named polyctena blocks
##       (module_di25/data/di25_three_blocks.rds): collapse each region's
##       units to ONE representative before computing genome-wide summaries,
##       instead of letting a multi-unit region contribute several
##       (non-independent) points.
##
## Run from the formica_hybrid repo root, after pp_prep_units.R and
## pp_local_concordance.R:
##   Rscript module_population_partitioning/R/pp_extra_robustness.R
## Writes: module_population_partitioning/data/pp_extra_robustness.rds
## =========================================================================
suppressMessages(library(data.table))
OUTDIR <- "module_population_partitioning/data"
obj <- readRDS(file.path(OUTDIR, "pp_units_Fmat.rds"))
res <- readRDS(file.path(OUTDIR, "pp_concordance_results.rds"))
u <- res$u; Fmat <- obj$Fmat; setDT(u); setorder(u, ChrNum, Pos)

## helper: recompute the near(<=100kb) local |r|/r per unit and the genome-wide
## FST-vs-similarity Spearman rho, from an arbitrary population x unit matrix
recompute_near <- function(Fm, u) {
  chrs <- unique(u$Chr)
  near_r <- near_absr <- rep(NA_real_, nrow(u))
  for (ch in chrs) {
    idx <- which(u$Chr == ch); if (length(idx) < 2) next
    sub <- Fm[, idx, drop = FALSE]; pos <- u$Pos[idx]
    Rm <- suppressWarnings(cor(sub, use = "pairwise.complete.obs"))
    for (k in seq_along(idx)) {
      d <- abs(pos - pos[k]); near <- which(d > 0 & d <= 1e5)
      if (!length(near)) next
      rs <- Rm[k, near]
      near_r[idx[k]] <- mean(rs, na.rm = TRUE); near_absr[idx[k]] <- mean(abs(rs), na.rm = TRUE)
    }
  }
  list(near_r = near_r, near_absr = near_absr,
      rho_absr = cor(u$FST, near_absr, use = "pairwise.complete.obs", method = "spearman"),
      rho_r    = cor(u$FST, near_r,    use = "pairwise.complete.obs", method = "spearman"))
}

## ---------------------------------------------------------------------
## (a) leave-one-population-out
## ---------------------------------------------------------------------
cat("=== (a) leave-one-population-out: genome-wide FST vs local-similarity rho ===\n")
pops20 <- rownames(Fmat)
loo <- rbindlist(lapply(pops20, function(p) {
  Fm <- Fmat[rownames(Fmat) != p, , drop = FALSE]
  out <- recompute_near(Fm, u)
  data.table(dropped = p, rho_absr = round(out$rho_absr, 3), rho_r = round(out$rho_r, 3))
}))
full <- recompute_near(Fmat, u)
cat(sprintf("full 20-population rho: |r| = %.3f, signed r = %.3f\n", full$rho_absr, full$rho_r))
print(loo[order(-abs(rho_r - full$rho_r))])

## ---------------------------------------------------------------------
## (b) once-per-region collapse of the 3 named polyctena blocks
## ---------------------------------------------------------------------
cat("\n=== (b) once-per-region: collapse the 3 named blocks to 1 representative unit each ===\n")
blk <- readRDS("module_di25/data/di25_three_blocks.rds")
blk_list <- strsplit(blk$group_ids, ",")
names(blk_list) <- paste0(blk$anchor, "_Chr", blk$chr)
## representative = the region's own single best-FST unit (avoids inventing a new average)
rep_ids <- vapply(blk_list, function(ids) u[group_id %in% ids][which.max(FST), group_id], character(1))
cat("region representatives (highest-FST unit per named block):\n"); print(data.table(region = names(rep_ids), rep_unit = rep_ids))

all_block_ids <- unlist(blk_list)
drop_ids <- setdiff(all_block_ids, rep_ids)                 # the non-representative units to drop
u_collapsed <- u[!group_id %in% drop_ids]
cat(sprintf("units: %d full -> %d once-per-region-collapsed (dropped %d redundant block units)\n",
            nrow(u), nrow(u_collapsed), length(drop_ids)))
rho_full <- cor(u$FST, u$near_absr, use = "pairwise.complete.obs", method = "spearman")
rho_full_r <- cor(u$FST, u$near_r, use = "pairwise.complete.obs", method = "spearman")
rho_collapsed <- cor(u_collapsed$FST, u_collapsed$near_absr, use = "pairwise.complete.obs", method = "spearman")
rho_collapsed_r <- cor(u_collapsed$FST, u_collapsed$near_r, use = "pairwise.complete.obs", method = "spearman")
cat(sprintf("Spearman FST vs near_absr: full = %.3f, once-per-region-collapsed = %.3f\n", rho_full, rho_collapsed))
cat(sprintf("Spearman FST vs near_r:    full = %.3f, once-per-region-collapsed = %.3f\n", rho_full_r, rho_collapsed_r))

saveRDS(list(loo = loo, full_rho = full, rep_ids = rep_ids, drop_ids = drop_ids,
            rho_full = rho_full, rho_full_r = rho_full_r,
            rho_collapsed = rho_collapsed, rho_collapsed_r = rho_collapsed_r),
        file.path(OUTDIR, "pp_extra_robustness.rds"))
cat("\n[extra-robustness] saved -> pp_extra_robustness.rds\n")
