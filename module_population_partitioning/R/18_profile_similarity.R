## =========================================================================
## module_population_partitioning -- candidate-locus follow-up 18: items 2
## and 3, locus-by-locus profile similarity and physical-distance decay.
##
## IMPORTANT LABELLING: the locus x locus statistic here is the Pearson
## correlation between two loci's POPULATION-LEVEL oriented allele-frequency
## PROFILES (one value per of the 19 hybrid populations) -- it is NOT
## within-population gametic linkage disequilibrium (which would require
## individual-level phased or dosage-covariance data within a single
## population). Every output/figure states this explicitly.
##
## For each target/status candidate set: locus x locus correlation heatmap;
## same-chromosome-near (<=100kb) vs same-chromosome-far vs cross-chromosome
## comparison; chromosome-block-bootstrap 95% CI (mirrors pp_block_
## bootstrap.R's per-chromosome-cell resampling); a decay-slope test only
## when enough same-chromosome pairs exist (MIN_PAIRS_FOR_DECAY guard --
## small candidate sets, esp. the 1-2 PC1/PC2 and 10 bio_winter floor
## survivors, are expected to have too few same-chromosome pairs for a
## meaningful slope, and this is reported as "not estimable", not silently
## computed on 0-1 pairs).
##
## Run from the formica_hybrid repo root, after 15_candidate_import.R:
##   Rscript module_population_partitioning/R/18_profile_similarity.R
## Reads : module_population_partitioning/data/followup/15_candidate_data.rds
## Writes: module_population_partitioning/data/followup/18_similarity_raw.rds
##         module_population_partitioning/Figures/followup/18_simheatmap_<target>_<status>.png
##         module_population_partitioning/Figures/followup/18_distance_decay_<target>.png
## =========================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
t_start <- Sys.time()
OUTDIR <- "module_population_partitioning/data/followup"
FIGDIR <- "module_population_partitioning/Figures/followup"
r15 <- readRDS(file.path(OUTDIR, "15_candidate_data.rds"))
u <- r15$u; Fmat_all <- r15$Fmat_all
stopifnot("candidate unit table must have exactly 18,361 rows -- rerun 15_candidate_import.R" = nrow(u) == 18361L,
         "Fmat_all columns must match u$group_id in order" = identical(colnames(Fmat_all), u$group_id))

MIN_PAIRS_FOR_DECAY <- 10L
NEAR_BP <- 1e5
B_BOOT  <- 2000
set.seed(20260914)
theme_ms <- theme_bw(base_size = 11) + theme(strip.background = element_blank(), panel.grid.minor = element_blank())

## ---------------------------------------------------------------------
## all-pairs correlation + distance classification for one locus set.
## Loci with UNDEFINED orientation (fixed for the same allele in both
## parental species -- see 15_candidate_import.R) have an all-NA profile
## and cannot be correlated with anything; they are dropped here (reported,
## not silently discarded) rather than plotted as an uninformative grey
## row/column.
## ---------------------------------------------------------------------
build_pairs <- function(ids) {
  all_na <- vapply(ids, function(id) all(is.na(Fmat_all[, id])), logical(1))
  n_dropped <- sum(all_na)
  ids <- ids[!all_na]
  sub <- u[match(ids, group_id)]
  Fsub <- Fmat_all[, ids, drop = FALSE]
  Rm <- suppressWarnings(cor(Fsub, use = "pairwise.complete.obs"))   # locus x locus, correlation of POPULATION PROFILES
  n <- length(ids)
  ii <- rep(seq_len(n), times = n); jj <- rep(seq_len(n), each = n)
  keep <- ii < jj; ii <- ii[keep]; jj <- jj[keep]
  same_chr <- sub$Chr[ii] == sub$Chr[jj]
  dist_bp <- ifelse(same_chr, abs(sub$Pos[ii] - sub$Pos[jj]), NA_real_)
  dcat <- ifelse(!same_chr, "cross-chr", ifelse(dist_bp <= NEAR_BP, "same-chr near (<=100kb)", "same-chr far"))
  list(Rm = Rm, ids = ids, n_dropped_undef = n_dropped,
      pairs = data.table(i = ii, j = jj, group_id_i = ids[ii], group_id_j = ids[jj],
                         Chr_i = sub$Chr[ii], Chr_j = sub$Chr[jj],
                         dist_bp = dist_bp, dcat = dcat, r = Rm[cbind(ii, jj)], absr = abs(Rm[cbind(ii, jj)])))
}

## ---------------------------------------------------------------------
## chromosome-block bootstrap CI on mean r/|r| by distance category
## (mirrors pp_block_bootstrap.R's boot_stat(): per-chromosome-PAIR
## sufficient statistics, resample the 26 chromosomes with replacement --
## a pair contributes to whichever chromosome its FIRST locus is on, for
## cross-chr pairs; same-chr pairs are keyed by their shared chromosome)
## ---------------------------------------------------------------------
block_boot_dcat <- function(pairs, chrs) {
  pairs <- copy(pairs)
  pairs[, boot_chr := ifelse(dcat == "cross-chr", Chr_i, Chr_i)]   # key every pair by locus i's chromosome
  cell <- pairs[, .(n = .N, sum_r = sum(r, na.rm = TRUE), sum_absr = sum(absr, na.rm = TRUE)), by = .(boot_chr, dcat)]
  levs <- sort(unique(pairs$dcat))
  out <- matrix(NA_real_, B_BOOT, length(levs) * 2, dimnames = list(NULL, c(paste0(levs, "_r"), paste0(levs, "_absr"))))
  for (b in seq_len(B_BOOT)) {
    draw <- sample(chrs, length(chrs), replace = TRUE)
    dt <- rbindlist(lapply(draw, function(ch) cell[boot_chr == ch]))
    if (!nrow(dt)) next
    agg <- dt[, .(n = sum(n), sum_r = sum(sum_r), sum_absr = sum(sum_absr)), by = dcat]
    agg[, mean_r := sum_r / n]; agg[, mean_absr := sum_absr / n]
    mr <- setNames(agg$mean_r, agg$dcat); ma <- setNames(agg$mean_absr, agg$dcat)
    out[b, paste0(levs, "_r")]    <- mr[levs]
    out[b, paste0(levs, "_absr")] <- ma[levs]
  }
  out
}

## ---------------------------------------------------------------------
## build one target/status result: heatmap figure + distance-category
## comparison + block-bootstrap CI + decay-slope test (if estimable)
## ---------------------------------------------------------------------
run_one <- function(tgt, status, ids) {
  bp <- build_pairs(ids)
  ids <- bp$ids   ## undefined-orientation loci dropped inside build_pairs()
  pairs <- bp$pairs
  n_near <- sum(pairs$dcat == "same-chr near (<=100kb)")
  n_far  <- sum(pairs$dcat == "same-chr far")
  n_cross <- sum(pairs$dcat == "cross-chr")
  cat(sprintf("[similarity] %-10s %-8s: n=%3d loci (%d undefined-orientation dropped), %4d pairs (near=%d, far=%d, cross-chr=%d)\n",
              tgt, status, length(ids), bp$n_dropped_undef, nrow(pairs), n_near, n_far, n_cross))
  if (length(ids) < 2L) {
    cat(sprintf("    -> <2 loci with defined orientation remain, skipping similarity analysis for this set\n"))
    return(list(target = tgt, status = status, ids = ids, n_dropped_undef = bp$n_dropped_undef, skipped = TRUE))
  }

  obs <- pairs[, .(n = .N, mean_r = mean(r, na.rm = TRUE), mean_absr = mean(absr, na.rm = TRUE)), by = dcat]
  n_same_chr_pairs <- n_near + n_far
  estimable <- n_same_chr_pairs >= MIN_PAIRS_FOR_DECAY
  ci <- NULL; slope <- NULL
  if (nrow(pairs) >= MIN_PAIRS_FOR_DECAY) {
    boot <- block_boot_dcat(pairs, r15$chrs)
    levs <- sort(unique(pairs$dcat))
    ci <- rbindlist(lapply(levs, function(lv) {
      data.table(dcat = lv, mean_r = obs[dcat == lv, mean_r],
                lo_r = quantile(boot[, paste0(lv, "_r")], 0.025, na.rm = TRUE),
                hi_r = quantile(boot[, paste0(lv, "_r")], 0.975, na.rm = TRUE),
                mean_absr = obs[dcat == lv, mean_absr],
                lo_absr = quantile(boot[, paste0(lv, "_absr")], 0.025, na.rm = TRUE),
                hi_absr = quantile(boot[, paste0(lv, "_absr")], 0.975, na.rm = TRUE))
    }))
  }
  if (estimable) {
    same_chr_pairs <- pairs[dcat != "cross-chr"]
    fit <- lm(absr ~ dist_bp, same_chr_pairs)
    slope <- list(obs_slope = coef(fit)[2], n_pairs = nrow(same_chr_pairs))
  } else {
    cat(sprintf("    -> decay slope NOT ESTIMABLE (%d same-chromosome pairs, need >= %d)\n", n_same_chr_pairs, MIN_PAIRS_FOR_DECAY))
  }

  ## figure: locus x locus correlation heatmap, chromosome-ordered
  sub <- u[match(ids, group_id)]; setorder(sub, ChrNum, Pos)
  ord_ids <- sub$group_id
  Rm_ord <- bp$Rm[ord_ids, ord_ids]
  dm <- as.data.table(as.table(Rm_ord)); setnames(dm, c("locus_i", "locus_j", "r"))
  dm[, locus_i := factor(locus_i, levels = ord_ids)][, locus_j := factor(locus_j, levels = rev(ord_ids))]
  fig <- ggplot(dm, aes(locus_i, locus_j, fill = r)) + geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "grey95", high = "#b2182b", midpoint = 0, limits = c(-1, 1),
                        name = "Pearson r\n(population\nprofiles)") +
    labs(title = sprintf("%s (%s, n=%d loci): locus-by-locus correlation\namong POPULATION PROFILES (NOT gametic LD)", tgt, status, length(ord_ids)),
        subtitle = sprintf("correlation of oriented allele-frequency profiles across the 19 hybrid populations%s",
                           if (bp$n_dropped_undef > 0) sprintf("; %d undefined-orientation loci excluded", bp$n_dropped_undef) else ""),
        x = NULL, y = NULL) +
    theme_ms + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                     plot.title = element_text(size = 11), plot.subtitle = element_text(size = 7.5))
  fname <- sprintf("18_simheatmap_%s_%s.png", tgt, status)
  ggsave(file.path(FIGDIR, fname), fig, width = 7, height = 6.6, dpi = 200)

  list(target = tgt, status = status, ids = ord_ids, n_dropped_undef = bp$n_dropped_undef, pairs = pairs, obs = obs, ci = ci, slope = slope,
      n_same_chr_pairs = n_same_chr_pairs, estimable = estimable, figure = fname, skipped = FALSE)
}

## ---------------------------------------------------------------------
## loop over targets x status (skip empty sets, matching script 17)
## ---------------------------------------------------------------------
sim_results <- list()
for (tgt in names(r15$TARGETS)) {
  sets <- list(raw = r15$cand_ids[[tgt]], floor = r15$floor_ids[[tgt]], topN = r15$topN_ids[[tgt]])
  for (status in names(sets)) {
    ids <- sets[[status]]
    if (length(ids) < 2L) { cat(sprintf("[similarity] %-10s %-8s: <2 loci, skipped (no pairs possible)\n", tgt, status)); next }
    sim_results[[paste(tgt, status, sep = "_")]] <- run_one(tgt, status, ids)
  }
}

## ---------------------------------------------------------------------
## one distance-decay summary figure per target (raw-status pairs, the
## largest/most-informative set per target)
## ---------------------------------------------------------------------
for (tgt in names(r15$TARGETS)) {
  key <- paste(tgt, "raw", sep = "_")
  if (is.null(sim_results[[key]]) || is.null(sim_results[[key]]$ci)) next
  ci <- sim_results[[key]]$ci
  fig <- ggplot(ci, aes(dcat, mean_absr, ymin = lo_absr, ymax = hi_absr)) +
    geom_pointrange() +
    labs(title = sprintf("%s raw candidates: profile similarity by distance category", tgt),
        subtitle = sprintf("chromosome-block bootstrap 95%% CI (n=%d reps); %s",
                           B_BOOT, if (sim_results[[key]]$estimable) sprintf("decay slope estimable (n=%d same-chr pairs)", sim_results[[key]]$n_same_chr_pairs)
                                    else sprintf("decay slope NOT estimable (only %d same-chr pairs)", sim_results[[key]]$n_same_chr_pairs)),
        x = NULL, y = "mean |r| (population-profile correlation)") +
    theme_ms + theme(plot.subtitle = element_text(size = 8))
  ggsave(file.path(FIGDIR, sprintf("18_distance_decay_%s.png", tgt)), fig, width = 6.5, height = 5, dpi = 200)
}

result <- list(sim_results = sim_results, MIN_PAIRS_FOR_DECAY = MIN_PAIRS_FOR_DECAY, NEAR_BP = NEAR_BP, B_BOOT = B_BOOT,
              session_info = sessionInfo(), run_time = Sys.time(),
              elapsed_secs = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
saveRDS(result, file.path(OUTDIR, "18_similarity_raw.rds"))
cat(sprintf("\n[similarity] saved -> %s (elapsed %.1fs)\n", file.path(OUTDIR, "18_similarity_raw.rds"), result$elapsed_secs))
