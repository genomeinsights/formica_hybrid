## =========================================================================
## module_population_partitioning -- 04: robustness checks + PCA of sorted-
## unit population profiles (task 6/7 of the original brief -- a first pass,
## not exhaustive; see the module README "not yet run" list).
##
## (a) is the (weak) genome-wide FST-vs-local-similarity association driven by
##     the handful of giant, near-fully-linked low-recombination clusters?
## (b) best-SNP (eMLG, >2-marker clusters) vs representative-only (1-2 marker)
##     units -- does the unit-representation choice change the conclusion?
## (c) same, excluding the 3 previously identified dominant polyctena blocks
##     (module_di25/data/di25_three_blocks.rds: Chr5/25/26).
## (d) PCA of population profiles among SORTED units only, to ask whether a
##     small number of recurring partitions explains most of them. IMPORTANT:
##     must row-center (each unit's own profile, centred on ITS OWN mean
##     across populations -- matching the cor()-based concordance analysis in
##     pp_local_concordance.R) before the SVD. prcomp()'s own `center=TRUE`
##     centers COLUMNS (populations), which instead captures the between-
##     population hybrid-index baseline (the same thing that produces the
##     ~0.19 permutation-null floor in pp_permutation_null.R) and inflates
##     PC1 to ~71% -- a different, less interesting quantity. Row-centering
##     first gives PC1 ~= 11%, consistent with the weak pairwise |r| found
##     genome-wide.
##
## Run from the formica_hybrid repo root, after pp_prep_units.R and
## pp_local_concordance.R:
##   Rscript module_population_partitioning/R/pp_robustness_pca.R
## Reads : module_population_partitioning/data/pp_units_Fmat.rds
##         module_population_partitioning/data/pp_concordance_results.rds
##         module_di25/data/di25_three_blocks.rds
## Writes: module_population_partitioning/data/pp_robustness.rds
## =========================================================================
suppressMessages(library(data.table))
OUTDIR <- "module_population_partitioning/data"
obj <- readRDS(file.path(OUTDIR, "pp_units_Fmat.rds"))
res <- readRDS(file.path(OUTDIR, "pp_concordance_results.rds"))
u <- res$u; Fmat <- obj$Fmat
setDT(u)

## ---------------------------------------------------------------------
## (a) big low-recombination blocks (n_loci > 50) vs the rest
## ---------------------------------------------------------------------
u[, big := n_loci_g > 50]
cat("=== (a) big low-recomb blocks (n_loci>50) vs rest ===\n")
print(u[, .N, by = big])
tab_big <- u[, .(rho = cor(FST, near_absr, use = "pairwise.complete.obs", method = "spearman"),
                 median_FST = median(FST, na.rm = TRUE), mean_near_absr = mean(near_absr, na.rm = TRUE)), by = big]
cat("FST vs near_absr, Spearman rho, by block-size class:\n"); print(tab_big)

## ---------------------------------------------------------------------
## (b) best-SNP eMLG units (is_emlg) vs representative-only (1-2 marker)
## ---------------------------------------------------------------------
cat("\n=== (b) best-SNP eMLG (is_emlg) vs representative-only units ===\n")
tab_emlg <- u[, .(rho = cor(FST, near_absr, use = "pairwise.complete.obs", method = "spearman"),
                  n = .N, mean_near_absr = mean(near_absr, na.rm = TRUE)), by = is_emlg]
print(tab_emlg)

## ---------------------------------------------------------------------
## (c) excluding the 3 previously identified dominant polyctena blocks
## ---------------------------------------------------------------------
blk <- readRDS("module_di25/data/di25_three_blocks.rds")
blk_ids <- unlist(strsplit(blk$group_ids, ","))
cat(sprintf("\n=== (c) excluding the 3 named polyctena blocks (%d units: %s) ===\n",
            length(blk_ids), paste(blk_ids, collapse = ", ")))
u[, in_block3 := group_id %in% blk_ids]
print(u[, .N, by = in_block3])
tab_blk3 <- u[, .(rho = cor(FST, near_absr, use = "pairwise.complete.obs", method = "spearman")), by = in_block3]
print(tab_blk3)
rho_excl <- cor(u[in_block3 == FALSE]$FST, u[in_block3 == FALSE]$near_absr, use = "pairwise.complete.obs", method = "spearman")
cat(sprintf("overall Spearman FST vs near_absr excluding the 3 named blocks: %.3f (n=%d)\n",
            rho_excl, u[in_block3 == FALSE & !is.na(FST) & !is.na(near_absr), .N]))

## ---------------------------------------------------------------------
## (d) PCA of population profiles among SORTED units (row-centered; see
##     header note on why this differs from naive prcomp(center=TRUE))
## ---------------------------------------------------------------------
cat("\n=== (d) PCA of sorted-unit population profiles (row-centered) ===\n")
sorted_ids <- u[sort_class != "unsorted", group_id]
M <- t(Fmat[, sorted_ids, drop = FALSE])                 # units (rows) x populations (cols)
M_rc <- M - rowMeans(M, na.rm = TRUE)                     # each unit centred on ITS OWN mean
keep_rows <- stats::complete.cases(M_rc)                  # drop any unit with residual NA (rare)
pc <- prcomp(M_rc[keep_rows, ], center = FALSE, scale. = FALSE)
ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
cat(sprintf("row-centered PCA on %d/%d sorted units x %d populations:\n", sum(keep_rows), length(sorted_ids), ncol(M)))
cat("  PC1-PC6 variance explained (%):", round(ve[1:6], 1), "\n")
cat("  cumulative PC1-PC6 (%):", round(sum(ve[1:6]), 1), "\n")

## for comparison, the naive (column-centered) version -- captures the
## between-population baseline instead, kept only to document the contrast
pc_naive <- prcomp(M, center = TRUE, scale. = FALSE)
ve_naive <- 100 * pc_naive$sdev^2 / sum(pc_naive$sdev^2)
cat(sprintf("[for reference] naive column-centered PCA PC1 = %.1f%% -- this is the between-population\n", ve_naive[1]))
cat("  hybrid-index baseline, NOT evidence of a shared unit-level partition; see header note.\n")

saveRDS(list(big_block = tab_big, emlg_vs_rep = tab_emlg, exclude_block3 = tab_blk3,
            rho_excl_block3 = rho_excl, pca_var_explained = ve[1:6],
            pca_var_explained_naive = ve_naive[1:6]),
        file.path(OUTDIR, "pp_robustness.rds"))
cat("\n[robustness] saved -> pp_robustness.rds\n")
