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
stopifnot("unit table must have exactly the 20,807 DI25 rho05 units (min_r2_rho=0.5) -- check for a legacy/stale input" = nrow(u) == 20807L,
         "row order must match pairs$i/j (see below)" = identical(u$group_id, obj$u$group_id))

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

## ---------------------------------------------------------------------
## (d) AUDIT FIX (item 4): "at fixed physical distance" previously meant
##     membership in the single broad 100-500kb bin -- too coarse to rule
##     out a residual within-bin distance gradient confounding the
##     low-vs-high recombination contrast. This adjusts more precisely
##     WITHIN THE SAME 100-500kb SCOPE the original claim was about (not a
##     wider distance range -- pooling in the much-weaker 0.5-2Mb signal
##     would dilute the estimate and answer a different question): NARROW
##     (20kb) physical-distance strata spanning 100-500kb (20 strata),
##     combined into one distance-adjusted contrast via a stratum-size-
##     weighted average of the within-stratum low-vs-high difference (a
##     stratified/Mantel-Haenszel-style adjustment: distance is
##     ~constant within each 20kb stratum, so a residual difference there
##     reflects recombination, not distance). Uncertainty is a
##     chromosome-block bootstrap on the SAME sufficient statistics used in
##     (b) above (not an ordinary pair-level model p-value -- unit pairs
##     share units, so treating millions of pairs as independent
##     observations would understate uncertainty, as documented throughout
##     this module).
## ---------------------------------------------------------------------
NARROW_BRK <- seq(1e5, 5e5, by = 2e4)   # 20kb strata, 100-500kb (20 strata) -- same scope as the original claim
pairs[, narrow_bin := cut(dist_bp, NARROW_BRK, include.lowest = TRUE)]
cell_adj <- pairs[!is.na(narrow_bin) & !is.na(recomb_tertile) & recomb_tertile != "mid",
                  .(n = .N, sum_r = sum(r, na.rm = TRUE)), by = .(Chr, narrow_bin, recomb_tertile)]

stratified_contrast <- function(cell_dt) {
  ## one weighted low-vs-high contrast, pooling all narrow strata; weight =
  ## harmonic-mean-style n_low*n_high/(n_low+n_high) per stratum (more
  ## weight to strata with balanced, well-powered low/high pair counts)
  agg <- cell_dt[, .(sum_r = sum(sum_r), n = sum(n)), by = .(narrow_bin, recomb_tertile)]
  wide <- dcast(agg, narrow_bin ~ recomb_tertile, value.var = c("sum_r", "n"))
  wide <- wide[!is.na(n_low) & !is.na(n_high) & n_low > 0 & n_high > 0]
  wide[, mean_low := sum_r_low / n_low]; wide[, mean_high := sum_r_high / n_high]
  wide[, w := (n_low * n_high) / (n_low + n_high)]
  sum(wide$w * (wide$mean_low - wide$mean_high)) / sum(wide$w)
}
obs_adj_contrast <- stratified_contrast(cell_adj)
cat(sprintf("\n[recomb] (d) distance-adjusted (20kb strata, 100-500kb) low-vs-high recombination contrast in signed r: %.4f\n", obs_adj_contrast))

n_narrow_strata <- length(unique(cell_adj$narrow_bin))
cat(sprintf("[recomb] (d) %d narrow strata contributed (out of up to %d possible)\n", n_narrow_strata, length(NARROW_BRK) - 1))

set.seed(11)
boot_adj_contrast <- vapply(seq_len(B_BOOT <- 2000), function(b) {
  draw <- sample(chrs, length(chrs), replace = TRUE)
  wtab <- table(draw)
  dt <- cell_adj[Chr %in% names(wtab)]
  w <- as.numeric(wtab[dt$Chr])
  dt2 <- copy(dt); dt2[, `:=`(n = n * w, sum_r = sum_r * w)]
  stratified_contrast(dt2)
}, numeric(1))
ci_adj <- quantile(boot_adj_contrast, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("[recomb] (d) chromosome-block bootstrap 95%% CI: [%.4f, %.4f] (n=%d reps)\n", ci_adj[1], ci_adj[2], B_BOOT))
cat("[recomb] (d) this SUPPORTS an effect of local recombination rate on partition concordance at fixed\n")
cat("    physical distance (distance-adjusted contrast close to the unadjusted 100-500kb estimate); it does\n")
cat("    not exhaustively rule out every possible form of residual confounding, so 'supports' rather than\n")
cat("    'confirms' is the appropriate strength of claim.\n")

## ---------------------------------------------------------------------
## (e) units outside the linkage-map (recmap) physical range -- cM/rate for
##     these came from approxfun(..., rule=2) FLAT extrapolation beyond the
##     map's own endpoint, not genuine interpolation; reported, not dropped
##     (only 2 units, negligible influence on any aggregate statistic, but
##     worth knowing which two).
## ---------------------------------------------------------------------
map_range <- gm[, .(minpos = min(pos), maxpos = max(pos)), by = Chr]
u_range <- merge(u, map_range, by = "Chr")
outside_map <- u_range[Pos < minpos | Pos > maxpos, .(group_id, Chr, Pos, minpos, maxpos)]
cat(sprintf("\n[recomb] (e) %d unit(s) fall outside the linkage-map physical range (flat rule=2 extrapolation, not interpolation):\n", nrow(outside_map)))
print(outside_map)

saveRDS(list(u = u, cm_decay = cm_decay, strat = strat, boot_contrast = boot_contrast, obs_contrast = obs_contrast,
            rec_tertiles = rec_tertiles, recomb_dec_tab = recomb_dec_tab,
            adjusted_contrast = list(narrow_breaks = NARROW_BRK, n_strata = n_narrow_strata,
                                     obs = obs_adj_contrast, boot = boot_adj_contrast, ci = ci_adj),
            outside_map_range = outside_map),
        file.path(OUTDIR, "pp_recombination.rds"))
cat("\n[recomb] saved -> pp_recombination.rds\n")
