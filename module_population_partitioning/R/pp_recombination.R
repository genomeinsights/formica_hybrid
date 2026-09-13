## =========================================================================
## module_population_partitioning -- 11: does local partition concordance
## track genetic (recombination) distance rather than merely physical
## distance? (README "Not yet run" / AUDIT.md priority item 5.)
##
## Everywhere else in this module, "the decay is short-range" has been
## described as "consistent with linkage" rather than "LD-driven", because
## physical distance is only a proxy for recombination distance -- the same
## bp span can span very different genetic distances depending on local
## recombination rate. This script joins the genetic map
## (data/Frufa_DTOL_PR.ref_genome.recmap, the same map module_di25's own
## clustering scripts use) and asks directly:
##   (a) does plotting concordance against GENETIC (cM) distance instead of
##       physical (bp) distance collapse the decay curve, i.e. is cM
##       distance the better predictor?
##   (b) at a FIXED physical distance, is concordance higher in low-
##       recombination regions than high-recombination regions (the direct,
##       physical-distance-controlled test of the recombination hypothesis)?
##   (c) does a unit's own local recombination rate predict its local
##       (<=100kb) concordance to neighbours?
##
## Run from the formica_hybrid repo root, after pp_prep_units.R and
## pp_local_concordance.R:
##   Rscript module_population_partitioning/R/pp_recombination.R
## Reads : module_population_partitioning/data/pp_units_Fmat.rds
##         module_population_partitioning/data/pp_concordance_results.rds
##         module_population_partitioning/data/pp_all_pairs.csv.gz
##         data/Frufa_DTOL_PR.ref_genome.recmap
## Writes: module_population_partitioning/data/pp_recombination.rds
## =========================================================================
suppressMessages(library(data.table))
OUTDIR <- "module_population_partitioning/data"
obj <- readRDS(file.path(OUTDIR, "pp_units_Fmat.rds"))
res <- readRDS(file.path(OUTDIR, "pp_concordance_results.rds"))
u <- res$u; setDT(u); setorder(u, ChrNum, Pos)
stopifnot(identical(u$group_id, obj$u$group_id))   # row order must match pairs$i/j (see below)

## ---------------------------------------------------------------------
## 1. per-unit genetic position (cM) and local recombination rate (cM/Mb),
##    interpolated per chromosome (approx(..., rule=2), same convention as
##    module_di25/R/di25_ld_clustering.R's own genetic_map construction)
## ---------------------------------------------------------------------
gm <- fread("data/Frufa_DTOL_PR.ref_genome.recmap")
gm[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
setnames(gm, "cM/Mb", "cMMb")
stopifnot(all(unique(u$Chr) %in% unique(gm$Chr)))

u[, cM_pos := NA_real_][, recomb_rate := NA_real_]
for (ch in unique(u$Chr)) {
  gmc <- gm[Chr == ch]; setorder(gmc, pos)
  idx <- which(u$Chr == ch)
  f_cM <- approxfun(gmc$pos, gmc$cM, rule = 2)
  f_rate <- approxfun(gmc$pos, gmc$cMMb, rule = 2)
  u$cM_pos[idx] <- f_cM(u$Pos[idx])
  u$recomb_rate[idx] <- f_rate(u$Pos[idx])
}
cat(sprintf("[recomb] local recombination rate (cM/Mb) across %d units: median %.2f, range [%.2f, %.2f]\n",
            nrow(u), median(u$recomb_rate), min(u$recomb_rate), max(u$recomb_rate)))
cat(sprintf("[recomb] Spearman FST vs local recomb rate: rho = %.3f (sanity check, not the main question)\n",
            cor(u$FST, u$recomb_rate, use = "pairwise.complete.obs", method = "spearman")))

## ---------------------------------------------------------------------
## 2. join genetic distance + mean local recombination rate into the
##    within-chromosome pairs table (pairs$i/j are row indices into THIS u,
##    confirmed identical row order to the table pp_local_concordance.R built)
## ---------------------------------------------------------------------
pairs <- fread(file.path(OUTDIR, "pp_all_pairs.csv.gz"))
pairs[, cM_dist := abs(u$cM_pos[i] - u$cM_pos[j])]
pairs[, mean_recomb := (u$recomb_rate[i] + u$recomb_rate[j]) / 2]
cat(sprintf("[recomb] %d pairs joined; median cM distance %.4f, median physical/genetic ratio check: %.0f bp/cM (should be very large for close pairs, ~Mb/cM genome-wide)\n",
            nrow(pairs), median(pairs$cM_dist), median(pairs$dist_bp / pmax(pairs$cM_dist, 1e-6))))

## ---------------------------------------------------------------------
## (a) genetic-distance decay curve (signed r / |r|), vs the physical-
##     distance one already on record -- does cM distance collapse the curve?
## ---------------------------------------------------------------------
CM_BRK <- c(0, 0.001, 0.005, 0.02, 0.1, 0.5, 2, Inf)
CM_LAB <- c("0-0.001cM","0.001-0.005cM","0.005-0.02cM","0.02-0.1cM","0.1-0.5cM","0.5-2cM",">2cM")
pairs[, cmbin := cut(cM_dist, CM_BRK, labels = CM_LAB)]
cm_decay <- pairs[!is.na(cmbin), .(n = .N, mean_r = mean(r, na.rm = TRUE), mean_absr = mean(absr, na.rm = TRUE),
                                   mean_dist_bp = mean(dist_bp)), by = cmbin][order(cmbin)]
cat("\n[recomb] (a) local similarity vs GENETIC distance bin:\n"); print(cm_decay)

## ---------------------------------------------------------------------
## (b) at FIXED physical distance, does concordance differ by local
##     recombination rate? Stratify the existing bp distance bins
##     (100kb-500kb and 0.5-2Mb chosen: enough pairs, far enough that
##     recombination-rate variation has room to matter, short enough that a
##     real signal should still be visible above the null floor)
## ---------------------------------------------------------------------
rec_tertiles <- quantile(u$recomb_rate, c(1/3, 2/3), na.rm = TRUE)
pairs[, recomb_tertile := cut(mean_recomb, c(-Inf, rec_tertiles, Inf), labels = c("low", "mid", "high"))]
cat(sprintf("\n[recomb] (b) recombination-rate tertile cutoffs (cM/Mb): low<%.2f, mid, high>%.2f\n",
            rec_tertiles[1], rec_tertiles[2]))

FOCUS_BINS <- c("100-500kb", "0.5-2Mb")
strat <- pairs[dbin %in% FOCUS_BINS & !is.na(recomb_tertile), .(n = .N, mean_r = mean(r, na.rm = TRUE),
                                                                 mean_absr = mean(absr, na.rm = TRUE)), by = .(dbin, recomb_tertile)]
setorder(strat, dbin, recomb_tertile)
cat("[recomb] (b) local similarity by physical-distance bin x recombination-rate tertile:\n"); print(strat)

## chromosome-block bootstrap CI on the low-vs-high recombination contrast
## (mirrors pp_block_bootstrap.R's convention: resample the 26 chromosomes)
chrs <- unique(pairs$Chr)
cell <- pairs[dbin %in% FOCUS_BINS & !is.na(recomb_tertile),
             .(n = .N, sum_r = sum(r, na.rm = TRUE), sum_absr = sum(absr, na.rm = TRUE)),
             by = .(Chr, dbin, recomb_tertile)]
set.seed(7); B <- 2000
boot_contrast <- matrix(NA_real_, B, length(FOCUS_BINS), dimnames = list(NULL, FOCUS_BINS))
for (b in seq_len(B)) {
  draw <- sample(chrs, length(chrs), replace = TRUE)
  dt <- rbindlist(lapply(draw, function(ch) cell[Chr == ch]))
  agg <- dt[, .(sum_r = sum(sum_r), n = sum(n)), by = .(dbin, recomb_tertile)]
  agg[, mean_r := sum_r / n]
  for (fb in FOCUS_BINS) {
    lo <- agg[dbin == fb & recomb_tertile == "low", mean_r]
    hi <- agg[dbin == fb & recomb_tertile == "high", mean_r]
    if (length(lo) && length(hi)) boot_contrast[b, fb] <- lo - hi
  }
}
obs_contrast <- sapply(FOCUS_BINS, function(fb) {
  strat[dbin == fb & recomb_tertile == "low", mean_r] - strat[dbin == fb & recomb_tertile == "high", mean_r]
})
cat("\n[recomb] (b) low-minus-high recombination contrast in signed r, block-bootstrap 95% CI:\n")
for (fb in FOCUS_BINS) {
  ci <- quantile(boot_contrast[, fb], c(0.025, 0.975), na.rm = TRUE)
  cat(sprintf("  %s: observed = %.4f, 95%% CI [%.4f, %.4f]\n", fb, obs_contrast[fb], ci[1], ci[2]))
}

## ---------------------------------------------------------------------
## (c) per-unit local (<=100kb) concordance vs the unit's OWN local
##     recombination rate
## ---------------------------------------------------------------------
cat(sprintf("\n[recomb] (c) Spearman local(<=100kb) |r| vs local recomb rate: rho = %.3f (n=%d)\n",
            cor(u$near_absr, u$recomb_rate, use = "pairwise.complete.obs", method = "spearman"),
            sum(!is.na(u$near_absr) & !is.na(u$recomb_rate))))
cat(sprintf("[recomb] (c) Spearman local(<=100kb) signed r vs local recomb rate: rho = %.3f (n=%d)\n",
            cor(u$near_r, u$recomb_rate, use = "pairwise.complete.obs", method = "spearman"),
            sum(!is.na(u$near_r) & !is.na(u$recomb_rate))))
u[, recomb_dec := cut(recomb_rate, quantile(recomb_rate, seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)]
recomb_dec_tab <- u[!is.na(recomb_dec), .(n = .N, mean_recomb = mean(recomb_rate), mean_near_r = mean(near_r, na.rm = TRUE),
                                          mean_near_absr = mean(near_absr, na.rm = TRUE)), by = recomb_dec][order(recomb_dec)]
cat("[recomb] (c) local concordance by local-recombination-rate decile:\n"); print(recomb_dec_tab)

saveRDS(list(u = u, cm_decay = cm_decay, strat = strat, boot_contrast = boot_contrast, obs_contrast = obs_contrast,
            rec_tertiles = rec_tertiles, recomb_dec_tab = recomb_dec_tab),
        file.path(OUTDIR, "pp_recombination.rds"))
cat("\n[recomb] saved -> pp_recombination.rds\n")
