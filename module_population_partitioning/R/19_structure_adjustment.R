## =========================================================================
## module_population_partitioning -- candidate-locus follow-up 19: item 4,
## structure-adjusted versions of the heatmap and profile-similarity
## analyses.
##
## Each candidate's raw oriented-frequency profile (across the 19 hybrid
## populations) is regressed on the LEAVE-ONE-CHROMOSOME-OUT reference
## structure PC1 score (16_structure_reference.rds$H_loco_pc1, estimated
## from the CANONICAL genome-wide reference panel, never from candidate
## loci) and the RESIDUALS are carried forward. PC1+PC2 is saved as a named
## sensitivity variant. This mirrors pp_residualize_ancestry.R's purpose
## (separate genome-wide background structure from locus-specific signal)
## but adapts its mechanism: pp_residualize_ancestry.R resums an additive
## per-population ancestry MEAN, whereas the structure covariate here is a
## PCA score, which is not additive -- the LOCO refit already happened in
## script 16, this script only fits the per-locus OLS regression + keeps
## residuals.
##
## Both RAW (scripts 17/18) and ADJUSTED (this script) results are always
## reported side by side, never adjusted-only -- the difference between
## them is itself the biologically informative quantity requested.
##
## Run from the formica_hybrid repo root, after 16_structure_reference.R:
##   Rscript module_population_partitioning/R/19_structure_adjustment.R
## Reads : module_population_partitioning/data/followup/15_candidate_data.rds
##         module_population_partitioning/data/followup/16_structure_reference.rds
## Writes: module_population_partitioning/data/followup/19_structure_adjustment.rds
##         module_population_partitioning/Figures/followup/19_heatmap_adj_<target>_<status>.png
##         module_population_partitioning/Figures/followup/19_simheatmap_adj_<target>_<status>.png
##         module_population_partitioning/Figures/followup/19_distance_decay_adj_<target>.png
## =========================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
t_start <- Sys.time()
OUTDIR <- "module_population_partitioning/data/followup"
FIGDIR <- "module_population_partitioning/Figures/followup"
r15 <- readRDS(file.path(OUTDIR, "15_candidate_data.rds"))
r16 <- readRDS(file.path(OUTDIR, "16_structure_reference.rds"))
u <- r15$u; Fmat_all <- r15$Fmat_all
stopifnot("candidate unit table must have exactly 18,361 rows -- rerun 15_candidate_import.R" = nrow(u) == 18361L,
         "Fmat_all columns must match u$group_id in order" = identical(colnames(Fmat_all), u$group_id),
         "structure-reference population order must match Fmat_all's" = identical(rownames(r16$H_loco_pc1), rownames(Fmat_all)))

AQU <- "#21918C"; POL <- "#D3C93B"
theme_ms <- theme_bw(base_size = 11) + theme(strip.background = element_blank(), panel.grid.minor = element_blank())
NEAR_BP <- 1e5; MIN_PAIRS_FOR_DECAY <- 10L; B_BOOT <- 2000
set.seed(20260914)

## ---------------------------------------------------------------------
## 1. per-candidate-unit residualization: f_u ~ H_loco_pc1[, chr(u)]
##    (primary), and f_u ~ H_loco_pc1 + H_loco_pc2 (sensitivity)
## ---------------------------------------------------------------------
chr_of_unit <- u$Chr[match(colnames(Fmat_all), u$group_id)]
n_units <- ncol(Fmat_all); n_pops <- nrow(Fmat_all)
Resid1 <- matrix(NA_real_, n_pops, n_units, dimnames = dimnames(Fmat_all))
Resid2 <- matrix(NA_real_, n_pops, n_units, dimnames = dimnames(Fmat_all))
R2_pc1 <- R2_pc12 <- rep(NA_real_, n_units); names(R2_pc1) <- names(R2_pc12) <- colnames(Fmat_all)

for (ch in r15$chrs) {
  idx <- which(chr_of_unit == ch)
  if (!ch %in% colnames(r16$H_loco_pc1)) next
  H1 <- r16$H_loco_pc1[, ch]; H2 <- r16$H_loco_pc2[, ch]
  for (k in idx) {
    y <- Fmat_all[, k]; ok <- !is.na(y) & !is.na(H1)
    if (sum(ok) < 4) next
    fit1 <- lm(y[ok] ~ H1[ok])
    Resid1[ok, k] <- resid(fit1); R2_pc1[k] <- summary(fit1)$r.squared
    if (!anyNA(H2[ok])) {
      fit2 <- lm(y[ok] ~ H1[ok] + H2[ok])
      Resid2[ok, k] <- resid(fit2); R2_pc12[k] <- summary(fit2)$r.squared
    }
  }
}
u[, R2_structure_pc1 := R2_pc1[group_id]]
u[, R2_structure_pc12 := R2_pc12[group_id]]
cat(sprintf("[adjust] per-unit R2 (variance explained by LOCO reference PC1): median %.3f, mean %.3f (n=%d units fit)\n",
            median(R2_pc1, na.rm = TRUE), mean(R2_pc1, na.rm = TRUE), sum(!is.na(R2_pc1))))
cat(sprintf("[adjust] per-unit R2 (PC1+PC2 sensitivity): median %.3f, mean %.3f\n",
            median(R2_pc12, na.rm = TRUE), mean(R2_pc12, na.rm = TRUE)))

## ---------------------------------------------------------------------
## 2. adjusted heatmap (mirrors script 17's build_heatmap, generalized to
##    take the matrix as a parameter -- same helper-per-script convention)
## ---------------------------------------------------------------------
build_heatmap_adj <- function(tgt, status, ids, Fmat) {
  cfg <- r15$TARGETS[[tgt]]
  sub <- u[match(ids, group_id)]; setorder(sub, ChrNum, Pos)
  ord_ids <- sub$group_id
  Fsub <- Fmat[, ord_ids, drop = FALSE]
  keep <- colSums(!is.na(Fsub)) > 0
  if (sum(keep) < 1) return(NULL)
  Fsub <- Fsub[, keep, drop = FALSE]; ord_ids <- ord_ids[keep]
  Zsub <- suppressWarnings(scale(Fsub, center = TRUE, scale = TRUE))

  dm <- as.data.table(as.table(Zsub)); setnames(dm, c("pop", "group_id", "f_adj_z"))
  dm[, group_id := factor(group_id, levels = ord_ids)]
  p_hm <- ggplot(dm, aes(group_id, pop, fill = f_adj_z)) + geom_tile() +
    scale_fill_gradient2(low = POL, mid = "grey95", high = AQU, midpoint = 0, name = "standardized\nresidual freq.\n(z, within-locus)") +
    labs(title = sprintf("%s: %s (%s, n=%d loci) -- STRUCTURE-ADJUSTED (LOCO ref. PC1 removed)", tgt, r17_status_label(status), "adjusted", length(ord_ids)),
        x = NULL, y = NULL) +
    theme_ms + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(),
                     plot.title = element_text(size = 9.5))
  fname <- sprintf("19_heatmap_adj_%s_%s.png", tgt, status)
  ggsave(file.path(FIGDIR, fname), p_hm, width = max(7, min(20, length(ord_ids) * 0.15 + 4)), height = 6, dpi = 200, limitsize = FALSE)
  fname
}
r17_status_label <- function(status) c(raw = "raw candidates", floor = "floor survivors", topN = sprintf("top-%d ranked", r15$N_TOP))[status]

## ---------------------------------------------------------------------
## 3. adjusted similarity + distance-decay (mirrors script 18's
##    build_pairs()/block_boot_dcat(), generalized to take the matrix)
## ---------------------------------------------------------------------
build_pairs_adj <- function(ids, Fmat) {
  all_na <- vapply(ids, function(id) all(is.na(Fmat[, id])), logical(1))
  n_dropped <- sum(all_na)
  ids <- ids[!all_na]
  if (length(ids) < 2L) return(list(ids = ids, n_dropped_undef = n_dropped, pairs = NULL))
  sub <- u[match(ids, group_id)]
  Fsub <- Fmat[, ids, drop = FALSE]
  Rm <- suppressWarnings(cor(Fsub, use = "pairwise.complete.obs"))
  n <- length(ids); ii <- rep(seq_len(n), times = n); jj <- rep(seq_len(n), each = n)
  keep <- ii < jj; ii <- ii[keep]; jj <- jj[keep]
  same_chr <- sub$Chr[ii] == sub$Chr[jj]
  dist_bp <- ifelse(same_chr, abs(sub$Pos[ii] - sub$Pos[jj]), NA_real_)
  dcat <- ifelse(!same_chr, "cross-chr", ifelse(dist_bp <= NEAR_BP, "same-chr near (<=100kb)", "same-chr far"))
  list(Rm = Rm, ids = ids, n_dropped_undef = n_dropped,
      pairs = data.table(i = ii, j = jj, Chr_i = sub$Chr[ii], Chr_j = sub$Chr[jj],
                         dist_bp = dist_bp, dcat = dcat, r = Rm[cbind(ii, jj)], absr = abs(Rm[cbind(ii, jj)])))
}
block_boot_dcat_adj <- function(pairs, chrs) {
  pairs <- copy(pairs); pairs[, boot_chr := Chr_i]
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
    out[b, paste0(levs, "_r")] <- mr[levs]; out[b, paste0(levs, "_absr")] <- ma[levs]
  }
  out
}
run_one_adj <- function(tgt, status, ids, Fmat) {
  bp <- build_pairs_adj(ids, Fmat)
  if (is.null(bp$pairs)) return(list(target = tgt, status = status, n_dropped_undef = bp$n_dropped_undef, skipped = TRUE))
  pairs <- bp$pairs
  obs <- pairs[, .(n = .N, mean_r = mean(r, na.rm = TRUE), mean_absr = mean(absr, na.rm = TRUE)), by = dcat]
  n_same_chr <- sum(pairs$dcat != "cross-chr")
  estimable <- n_same_chr >= MIN_PAIRS_FOR_DECAY
  ci <- NULL
  if (nrow(pairs) >= MIN_PAIRS_FOR_DECAY) {
    boot <- block_boot_dcat_adj(pairs, r15$chrs)
    levs <- sort(unique(pairs$dcat))
    ci <- rbindlist(lapply(levs, function(lv) data.table(
      dcat = lv, mean_r = obs[dcat == lv, mean_r], lo_r = quantile(boot[, paste0(lv, "_r")], 0.025, na.rm = TRUE),
      hi_r = quantile(boot[, paste0(lv, "_r")], 0.975, na.rm = TRUE), mean_absr = obs[dcat == lv, mean_absr],
      lo_absr = quantile(boot[, paste0(lv, "_absr")], 0.025, na.rm = TRUE), hi_absr = quantile(boot[, paste0(lv, "_absr")], 0.975, na.rm = TRUE))))
  }
  hm_fig <- build_heatmap_adj(tgt, status, bp$ids, Fmat)
  ord_ids <- bp$ids; sub2 <- u[match(ord_ids, group_id)]; setorder(sub2, ChrNum, Pos); ord_ids <- sub2$group_id
  Rm_ord <- bp$Rm[ord_ids, ord_ids]
  dm <- as.data.table(as.table(Rm_ord)); setnames(dm, c("locus_i", "locus_j", "r"))
  dm[, locus_i := factor(locus_i, levels = ord_ids)][, locus_j := factor(locus_j, levels = rev(ord_ids))]
  fig <- ggplot(dm, aes(locus_i, locus_j, fill = r)) + geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "grey95", high = "#b2182b", midpoint = 0, limits = c(-1, 1), name = "Pearson r\n(residual\nprofiles)") +
    labs(title = sprintf("%s (%s, n=%d loci): STRUCTURE-ADJUSTED locus-by-locus\ncorrelation (NOT gametic LD)", tgt, status, length(ord_ids)),
        subtitle = sprintf("residuals after removing LOCO reference-panel PC1%s", if (bp$n_dropped_undef > 0) sprintf("; %d undefined-orientation loci excluded", bp$n_dropped_undef) else ""),
        x = NULL, y = NULL) +
    theme_ms + theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 10.5), plot.subtitle = element_text(size = 7.5))
  simfig <- sprintf("19_simheatmap_adj_%s_%s.png", tgt, status)
  ggsave(file.path(FIGDIR, simfig), fig, width = 7, height = 6.6, dpi = 200)
  cat(sprintf("[adjust] %-10s %-8s: n=%3d loci (%d undef dropped), %4d pairs, %s\n",
              tgt, status, length(ord_ids), bp$n_dropped_undef, nrow(pairs), if (estimable) sprintf("decay estimable (n=%d)", n_same_chr) else "decay NOT estimable"))
  list(target = tgt, status = status, ids = ord_ids, n_dropped_undef = bp$n_dropped_undef, pairs = pairs, obs = obs, ci = ci,
      n_same_chr_pairs = n_same_chr, estimable = estimable, heatmap_figure = hm_fig, simheatmap_figure = simfig, skipped = FALSE)
}

## ---------------------------------------------------------------------
## 4. loop targets x status x {pc1, pc1pc2} -- pc1 is primary, saved fully;
##    pc1pc2 saved as a lighter sensitivity check (stats only, no figures)
## ---------------------------------------------------------------------
adj_results_pc1 <- list()
for (tgt in names(r15$TARGETS)) {
  sets <- list(raw = r15$cand_ids[[tgt]], floor = r15$floor_ids[[tgt]], topN = r15$topN_ids[[tgt]])
  for (status in names(sets)) {
    ids <- sets[[status]]
    if (length(ids) < 2L) next
    adj_results_pc1[[paste(tgt, status, sep = "_")]] <- run_one_adj(tgt, status, ids, Resid1)
  }
}

## PC1+PC2 sensitivity: stats only (raw|r| by dcat), no figures, to keep scope bounded
sens_pc12 <- list()
for (tgt in names(r15$TARGETS)) {
  sets <- list(raw = r15$cand_ids[[tgt]], floor = r15$floor_ids[[tgt]], topN = r15$topN_ids[[tgt]])
  for (status in names(sets)) {
    ids <- sets[[status]]
    if (length(ids) < 2L) next
    bp <- build_pairs_adj(ids, Resid2)
    if (is.null(bp$pairs)) next
    sens_pc12[[paste(tgt, status, sep = "_")]] <- bp$pairs[, .(mean_absr = mean(absr, na.rm = TRUE), n = .N), by = dcat]
  }
}

## ---------------------------------------------------------------------
## 5. save
## ---------------------------------------------------------------------
result <- list(Resid1 = Resid1, Resid2 = Resid2, R2_pc1 = R2_pc1, R2_pc12 = R2_pc12,
              adj_results_pc1 = adj_results_pc1, sens_pc12 = sens_pc12,
              session_info = sessionInfo(), run_time = Sys.time(),
              elapsed_secs = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
saveRDS(result, file.path(OUTDIR, "19_structure_adjustment.rds"))
cat(sprintf("\n[adjust] saved -> %s (elapsed %.1fs)\n", file.path(OUTDIR, "19_structure_adjustment.rds"), result$elapsed_secs))
