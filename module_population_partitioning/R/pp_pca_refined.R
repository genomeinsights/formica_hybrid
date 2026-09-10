## =========================================================================
## module_population_partitioning -- 07: PCA refinement (AUDIT.md finding 6).
##
## The row-centered PCA in pp_robustness_pca.R found PC1 ~= 11% of variance
## among sorted units -- not one dominant shared partition, but on its own
## this does not rule out a real (just not overwhelming) recurring partition,
## nor explain WHY the naive column-centered PCA gave PC1 = 71%. This script:
##   (a) a row-wise population-label permutation reference for PC1 variance
##       (is 11% more than expected by chance alone?);
##   (b) the PC1/PC2 population LOADINGS, plus scores split by sort_class, to
##       check the audit's proposed mechanism for the naive 71% axis: that
##       failing to row-center leaves each unit's own overall level (higher
##       for aquilonia-classified units, lower for polyctena-classified ones)
##       as an uncontrolled, dominant source of variance -- i.e. the naive
##       PC1 mostly just separates aqu-sorted from pol-sorted units, not a
##       shared population partition;
##   (c) leave-one-population-out PCA (esp. Sielva, the F1-like colony);
##   (d) separate PCA summaries for aquilonia-sorted vs polyctena-sorted units.
##
## Run from the formica_hybrid repo root, after pp_prep_units.R and
## pp_local_concordance.R:
##   Rscript module_population_partitioning/R/pp_pca_refined.R
## Writes: module_population_partitioning/data/pp_pca_refined.rds
## =========================================================================
suppressMessages(library(data.table))
OUTDIR <- "module_population_partitioning/data"
obj <- readRDS(file.path(OUTDIR, "pp_units_Fmat.rds"))
res <- readRDS(file.path(OUTDIR, "pp_concordance_results.rds"))
u <- res$u; Fmat <- obj$Fmat; setDT(u)

sorted_ids <- u[sort_class != "unsorted", group_id]
M <- t(Fmat[, sorted_ids, drop = FALSE])            # units (rows) x populations (cols)
M_rc <- M - rowMeans(M, na.rm = TRUE)
keep_rows <- stats::complete.cases(M_rc)
M_rc <- M_rc[keep_rows, ]
sorted_ids_kept <- sorted_ids[keep_rows]
pc <- prcomp(M_rc, center = FALSE, scale. = FALSE)
ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
cat(sprintf("[pca] row-centered PCA, %d sorted units x %d pops: PC1-PC6 var.exp. (%%): %s\n",
            nrow(M_rc), ncol(M_rc), paste(round(ve[1:6], 1), collapse = ", ")))

## ---------------------------------------------------------------------
## (a) row-wise population-label permutation reference for PC1
##     -- for each unit (row), independently permute its 20 population
##     values, so each row keeps its own value distribution/variance but any
##     real between-unit covariance in WHICH populations are high/low is
##     destroyed; re-run PCA; repeat -> null distribution of PC1 %.
## ---------------------------------------------------------------------
set.seed(11)
B <- 500
pc1_null <- vapply(seq_len(B), function(b) {
  Mp <- t(apply(M_rc, 1, sample))                 # permute each row's 20 values independently
  pcb <- prcomp(Mp, center = FALSE, scale. = FALSE)
  100 * pcb$sdev[1]^2 / sum(pcb$sdev^2)
}, numeric(1))
cat(sprintf("[pca] row-permutation null for PC1: mean = %.1f%%, 95th pct = %.1f%% (observed PC1 = %.1f%%)\n",
            mean(pc1_null), quantile(pc1_null, 0.95), ve[1]))

## ---------------------------------------------------------------------
## (b) PC1/PC2 loadings (per population) + scores by sort_class
##     -- tests the audit's proposed mechanism for the naive (unrowcentered)
##     71% axis directly: does it separate aqu- vs pol-sorted units?
## ---------------------------------------------------------------------
loadings <- data.table(population = colnames(M_rc), PC1 = pc$rotation[, 1], PC2 = pc$rotation[, 2])
cat("\n[pca] PC1/PC2 population loadings (row-centered PCA):\n"); print(loadings[order(-PC1)])

scores <- data.table(group_id = sorted_ids_kept, PC1 = pc$x[, 1], PC2 = pc$x[, 2])
scores <- merge(scores, u[, .(group_id, sort_class)], by = "group_id")
cat("\n[pca] row-centered PC1 score by sort_class (mean, sd):\n")
print(scores[, .(n = .N, mean_PC1 = round(mean(PC1), 3), sd_PC1 = round(sd(PC1), 3)), by = sort_class])

## naive (column-centered only) PCA, for direct comparison -- does IT separate aqu/pol?
pc_naive <- prcomp(M, center = TRUE, scale. = FALSE)
ve_naive <- 100 * pc_naive$sdev^2 / sum(pc_naive$sdev^2)
scores_naive <- data.table(group_id = sorted_ids, PC1n = pc_naive$x[, 1])
scores_naive <- merge(scores_naive, u[, .(group_id, sort_class)], by = "group_id")
cat(sprintf("\n[pca] naive column-centered PCA: PC1 = %.1f%% of variance\n", ve_naive[1]))
cat("[pca] naive PC1 score by sort_class (mean, sd) -- tests whether naive PC1 mostly\n")
cat("      separates aqu-sorted from pol-sorted units (their differing row-level baseline):\n")
print(scores_naive[, .(n = .N, mean_PC1n = round(mean(PC1n), 1), sd_PC1n = round(sd(PC1n), 1)), by = sort_class])

## ---------------------------------------------------------------------
## (c) leave-one-population-out PCA (esp. Sielva)
## ---------------------------------------------------------------------
cat("\n[pca] leave-one-population-out: PC1 %% (row-centered PCA, 19 remaining populations)\n")
pops20 <- rownames(Fmat)
loo <- rbindlist(lapply(pops20, function(p) {
  Mp <- M[, colnames(M) != p, drop = FALSE]
  Mp_rc <- Mp - rowMeans(Mp, na.rm = TRUE)
  keep <- stats::complete.cases(Mp_rc)
  pcb <- prcomp(Mp_rc[keep, ], center = FALSE, scale. = FALSE)
  veb <- 100 * pcb$sdev^2 / sum(pcb$sdev^2)
  data.table(dropped = p, n_units = sum(keep), PC1 = round(veb[1], 1), PC2 = round(veb[2], 1))
}))
print(loo[order(-PC1)])
cat(sprintf("[pca] full-data (20 pop) PC1 for comparison: %.1f%%\n", ve[1]))

## ---------------------------------------------------------------------
## (d) separate PCA for aquilonia-sorted vs polyctena-sorted units
## ---------------------------------------------------------------------
cat("\n[pca] separate PCA by sort direction (row-centered):\n")
by_dir <- rbindlist(lapply(c("aquilonia", "polyctena"), function(cls) {
  ids <- u[sort_class == cls, group_id]
  Md <- t(Fmat[, ids, drop = FALSE]); Md_rc <- Md - rowMeans(Md, na.rm = TRUE)
  keep <- stats::complete.cases(Md_rc)
  pcd <- prcomp(Md_rc[keep, ], center = FALSE, scale. = FALSE)
  ved <- 100 * pcd$sdev^2 / sum(pcd$sdev^2)
  data.table(sort_class = cls, n_units = sum(keep), PC1 = round(ved[1], 1), PC2 = round(ved[2], 1), PC3 = round(ved[3], 1))
}))
print(by_dir)

saveRDS(list(ve = ve, pc1_null = pc1_null, loadings = loadings, scores = scores,
            ve_naive = ve_naive, scores_naive = scores_naive, loo = loo, by_dir = by_dir),
        file.path(OUTDIR, "pp_pca_refined.rds"))
cat("\n[pca] saved -> pp_pca_refined.rds\n")
