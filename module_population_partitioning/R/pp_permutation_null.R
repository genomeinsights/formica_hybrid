## =========================================================================
## module_population_partitioning -- 03: is the ~0.19-0.20 "background" |r|
## floor seen at large physical distance (pp_local_concordance.R, distance-bin
## table) real shared population-level structure, or just small-n (20
## populations) sampling noise for two literally independent profiles?
##
## Two checks: (1) genuine cross-chromosome unit pairs (unlinked by
## construction); (2) a population-label permutation null (keeps each unit's
## own value distribution, destroys any real between-unit covariance).
##
## Run from the formica_hybrid repo root, after pp_prep_units.R:
##   Rscript module_population_partitioning/R/pp_permutation_null.R
## Reads : module_population_partitioning/data/pp_units_Fmat.rds
## Writes: module_population_partitioning/data/pp_null_check.rds
## =========================================================================
suppressMessages(library(data.table))
OUTDIR <- "module_population_partitioning/data"
obj <- readRDS(file.path(OUTDIR, "pp_units_Fmat.rds"))
u <- obj$u; Fmat <- obj$Fmat
setDT(u)

set.seed(1)
## (1) cross-chromosome pairs: sample 20000 random unit pairs on DIFFERENT chromosomes
n <- ncol(Fmat)
i <- sample.int(n, 20000, replace = TRUE); j <- sample.int(n, 20000, replace = TRUE)
keep <- u$Chr[i] != u$Chr[j]
i <- i[keep][1:10000]; j <- j[keep][1:10000]
r_cross <- suppressWarnings(mapply(function(a, b) cor(Fmat[, a], Fmat[, b], use = "pairwise.complete.obs"), i, j))
cat(sprintf("[null] cross-chromosome random pairs (n=%d): mean r = %.3f, mean |r| = %.3f\n",
            length(r_cross), mean(r_cross, na.rm = TRUE), mean(abs(r_cross), na.rm = TRUE)))

## (2) permutation null: independently permute one profile's population labels
##     (so the 20 values are the same multiset, but any real between-unit or
##     between-population covariance structure is destroyed), then correlate.
##     20 hybrid populations -> this is the pure small-n sampling-noise floor.
set.seed(2)
B <- 10000
samp_units <- sample.int(n, 2 * B, replace = TRUE)
a_idx <- samp_units[1:B]; b_idx <- samp_units[(B + 1):(2 * B)]
r_perm <- vapply(seq_len(B), function(k) {
  x <- Fmat[, a_idx[k]]; y <- Fmat[, b_idx[k]]
  y <- sample(y)   # permute population labels on y only -> breaks real covariance
  suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
}, numeric(1))
cat(sprintf("[null] population-label permutation (n=%d, label-shuffled): mean r = %.3f, mean |r| = %.3f, 95th pct |r| = %.3f\n",
            B, mean(r_perm, na.rm = TRUE), mean(abs(r_perm), na.rm = TRUE), quantile(abs(r_perm), 0.95, na.rm = TRUE)))

## For reference: same random UNPERMUTED unrelated-unit pairs (any chromosome, any distance)
r_rand <- suppressWarnings(mapply(function(a, b) cor(Fmat[, a], Fmat[, b], use = "pairwise.complete.obs"), a_idx, b_idx))
cat(sprintf("[null] same random unit pairs, NO permutation (real profiles): mean r = %.3f, mean |r| = %.3f\n",
            mean(r_rand, na.rm = TRUE), mean(abs(r_rand), na.rm = TRUE)))

saveRDS(list(cross_chr = r_cross, perm = r_perm, rand_unpermuted = r_rand,
            null_floor_absr = mean(abs(r_perm), na.rm = TRUE)),
        file.path(OUTDIR, "pp_null_check.rds"))
cat("\n[null] saved -> pp_null_check.rds\n")
