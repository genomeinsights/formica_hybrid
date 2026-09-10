## =========================================================================
## module_population_partitioning -- 06: chromosome-block bootstrap CIs
## (AUDIT.md finding 4: the sd/sqrt(N) error bars in the original figures
## treat millions of non-independent unit PAIRS -- and units that recur in
## many pairs -- as independent observations, so they are far too narrow).
##
## Resamples CHROMOSOMES with replacement (matching the convention already
## used elsewhere in this repo, e.g. module_di25/R/di25_recomb_tau_sweep.R),
## so every pair/unit on a resampled chromosome moves together. Applied to
## the two curves that carried naive CIs before:
##   (i)  local |r| / signed r vs physical-distance bin (all within-chromosome
##        pairs, pp_local_concordance.R section C)
##  (ii) local |r| / signed r vs FST decile (pp_local_concordance.R section D,
##        near_r/near_absr per unit; FST deciles are fixed from the full
##        dataset, only chromosome membership is resampled)
##
## Implementation note: rather than literally resampling and re-summing
## millions of raw pairs/units per replicate, this precomputes small
## per-chromosome summary tables (chromosome x distance-bin; chromosome x
## FST-decile) and resamples CHROMOSOME DRAWS over those -- identical to a
## literal block bootstrap (every pair/unit on a chromosome is weighted by
## how many times that chromosome is drawn) but orders of magnitude faster.
##
## Run from the formica_hybrid repo root, after pp_local_concordance.R:
##   Rscript module_population_partitioning/R/pp_block_bootstrap.R
## Reads : module_population_partitioning/data/pp_concordance_results.rds
##         module_population_partitioning/data/pp_all_pairs.csv.gz
## Writes: module_population_partitioning/data/pp_block_bootstrap.rds
##   list(dist_boot = <distance-bin x {mean_r,mean_absr,lo,hi} block-bootstrap CIs>,
##        fst_boot  = <FST-decile x {mean_r,mean_absr,lo,hi} block-bootstrap CIs>)
## =========================================================================
suppressMessages(library(data.table))
OUTDIR <- "module_population_partitioning/data"
res <- readRDS(file.path(OUTDIR, "pp_concordance_results.rds"))
u <- res$u; setDT(u)
pairs <- fread(file.path(OUTDIR, "pp_all_pairs.csv.gz"))
B <- 2000
set.seed(42)

## ---------------------------------------------------------------------
## (i) distance-bin curve: per-chromosome x per-bin sufficient statistics
## ---------------------------------------------------------------------
chrs <- unique(pairs$Chr)
BRK <- c(0, 5e3, 2e4, 1e5, 5e5, 2e6, 1e7, Inf)
LAB <- c("0-5kb","5-20kb","20-100kb","100-500kb","0.5-2Mb","2-10Mb",">10Mb")
pairs[, dbin := cut(dist_bp, BRK, labels = LAB)]
cell_dist <- pairs[!is.na(dbin), .(n = .N, sum_r = sum(r, na.rm = TRUE), sum_absr = sum(absr, na.rm = TRUE)),
                   by = .(Chr, dbin)]
setkey(cell_dist, Chr)

boot_stat <- function(cell, group_col, chrs, B) {
  levs <- levels(cell[[group_col]]); if (is.null(levs)) levs <- sort(unique(cell[[group_col]]))
  out <- matrix(NA_real_, B, length(levs) * 2, dimnames = list(NULL, c(paste0(levs, "_r"), paste0(levs, "_absr"))))
  for (b in seq_len(B)) {
    draw <- sample(chrs, length(chrs), replace = TRUE)
    dt <- rbindlist(lapply(draw, function(ch) cell[Chr == ch]))
    agg <- dt[, .(n = sum(n), sum_r = sum(sum_r), sum_absr = sum(sum_absr)), by = c(group_col)]
    agg[, mean_r := sum_r / n]; agg[, mean_absr := sum_absr / n]
    mr <- setNames(agg$mean_r, agg[[group_col]]); ma <- setNames(agg$mean_absr, agg[[group_col]])
    out[b, paste0(levs, "_r")]    <- mr[levs]
    out[b, paste0(levs, "_absr")] <- ma[levs]
  }
  out
}
cat(sprintf("[bootstrap] distance-bin curve: %d replicates over %d chromosomes ...\n", B, length(chrs)))
boot_dist <- boot_stat(cell_dist, "dbin", chrs, B)

obs_dist <- pairs[!is.na(dbin), .(mean_r = mean(r, na.rm = TRUE), mean_absr = mean(absr, na.rm = TRUE)), by = dbin][order(dbin)]
dist_ci <- rbindlist(lapply(LAB, function(lv) {
  data.table(dbin = lv,
             mean_r = obs_dist[dbin == lv, mean_r], lo_r = quantile(boot_dist[, paste0(lv, "_r")], 0.025, na.rm = TRUE),
             hi_r = quantile(boot_dist[, paste0(lv, "_r")], 0.975, na.rm = TRUE),
             mean_absr = obs_dist[dbin == lv, mean_absr], lo_absr = quantile(boot_dist[, paste0(lv, "_absr")], 0.025, na.rm = TRUE),
             hi_absr = quantile(boot_dist[, paste0(lv, "_absr")], 0.975, na.rm = TRUE))
}))
dist_ci[, dbin := factor(dbin, levels = LAB)]; setorder(dist_ci, dbin)
cat("\n[bootstrap] distance-bin, chromosome-block bootstrap 95% CI:\n"); print(dist_ci)

## ---------------------------------------------------------------------
## (ii) FST-decile trend: per-chromosome x per-decile sufficient statistics
##      (deciles are FIXED, computed once on the full dataset; only
##      chromosome membership of the CONTRIBUTING units is resampled)
## ---------------------------------------------------------------------
u[, FST_dec := cut(FST, quantile(FST, seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)]
cell_fst <- u[!is.na(FST_dec) & !is.na(near_r), .(n = .N, sum_r = sum(near_r, na.rm = TRUE), sum_absr = sum(near_absr, na.rm = TRUE)),
             by = .(Chr, FST_dec)]
cell_fst[, FST_dec := factor(FST_dec, levels = 1:10)]
setkey(cell_fst, Chr)
cat(sprintf("\n[bootstrap] FST-decile trend: %d replicates over %d chromosomes ...\n", B, length(chrs)))
boot_fst <- boot_stat(cell_fst, "FST_dec", chrs, B)

obs_fst <- u[!is.na(FST_dec) & !is.na(near_r), .(mean_FST = mean(FST), mean_r = mean(near_r, na.rm = TRUE),
                                                 mean_absr = mean(near_absr, na.rm = TRUE)), by = FST_dec][order(FST_dec)]
fst_ci <- rbindlist(lapply(1:10, function(lv) {
  lvc <- as.character(lv)
  data.table(FST_dec = lv, mean_FST = obs_fst[FST_dec == lv, mean_FST],
             mean_r = obs_fst[FST_dec == lv, mean_r], lo_r = quantile(boot_fst[, paste0(lvc, "_r")], 0.025, na.rm = TRUE),
             hi_r = quantile(boot_fst[, paste0(lvc, "_r")], 0.975, na.rm = TRUE),
             mean_absr = obs_fst[FST_dec == lv, mean_absr], lo_absr = quantile(boot_fst[, paste0(lvc, "_absr")], 0.025, na.rm = TRUE),
             hi_absr = quantile(boot_fst[, paste0(lvc, "_absr")], 0.975, na.rm = TRUE))
}))
cat("\n[bootstrap] FST-decile, chromosome-block bootstrap 95% CI:\n"); print(fst_ci)

## overall trend test: is the FST-decile slope (mean_r ~ decile) significant under
## the block bootstrap? (simple bootstrap p-value on the OLS slope of mean_r on decile)
obs_slope_r    <- coef(lm(mean_r ~ FST_dec, obs_fst))[2]
obs_slope_absr <- coef(lm(mean_absr ~ FST_dec, obs_fst))[2]
boot_slopes_r    <- vapply(1:10, function(lv) boot_fst[, paste0(lv, "_r")], numeric(B))
boot_slopes_absr <- vapply(1:10, function(lv) boot_fst[, paste0(lv, "_absr")], numeric(B))
slope_r    <- apply(boot_slopes_r, 1, function(y) coef(lm(y ~ seq_len(10)))[2])
slope_absr <- apply(boot_slopes_absr, 1, function(y) coef(lm(y ~ seq_len(10)))[2])
cat(sprintf("\n[bootstrap] FST-decile slope (signed r):   observed = %.4f, block-bootstrap 95%% CI [%.4f, %.4f]\n",
            obs_slope_r, quantile(slope_r, 0.025, na.rm=TRUE), quantile(slope_r, 0.975, na.rm=TRUE)))
cat(sprintf("[bootstrap] FST-decile slope (absolute |r|): observed = %.4f, block-bootstrap 95%% CI [%.4f, %.4f]\n",
            obs_slope_absr, quantile(slope_absr, 0.025, na.rm=TRUE), quantile(slope_absr, 0.975, na.rm=TRUE)))

saveRDS(list(dist_ci = dist_ci, fst_ci = fst_ci,
            slope_r = slope_r, slope_absr = slope_absr,
            obs_slope_r = obs_slope_r, obs_slope_absr = obs_slope_absr),
        file.path(OUTDIR, "pp_block_bootstrap.rds"))
cat("\n[bootstrap] saved -> pp_block_bootstrap.rds\n")
