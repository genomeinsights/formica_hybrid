## =========================================================================
## module_population_partitioning -- 10: the four figures, REVISED per
## AUDIT.md. Changes from the first-pass version (see git history / README
## changelog for the original):
##   - signed r is now the PRIMARY concordance statistic (biologically
##     interpretable given aquilonia orientation: + = same direction,
##     - = opposing, ~0 = different populations or noise), |r| secondary;
##   - chromosome-block bootstrap 95% CIs (pp_block_bootstrap.R) replace the
##     naive sd/sqrt(N) error bars, which badly understated uncertainty for
##     millions of non-independent pairs;
##   - the empirical cross-chromosome baseline (pp_permutation_null.R) is the
##     primary reference band; the label-permutation floor is kept as a
##     secondary, explicitly-labelled sampling-noise reference;
##   - Fig 1 adds a zoomed panel on the named F7174 block so it isn't lost
##     among 177 equal-width columns;
##   - Fig 3 no longer attributes the (still weak) FST-concordance trend to
##     the 48 large clusters (verified: excluding them changes rho by <0.003);
##   - Fig 4 is replaced: genome-wide near_r before vs after residualizing
##     out genome-wide ancestry (pp_residualize_ancestry.R) is far more
##     informative than the original undifferentiated scatter.
##
## Run from the formica_hybrid repo root, after pp_prep_units.R,
## pp_local_concordance.R, pp_permutation_null.R, pp_block_bootstrap.R,
## pp_residualize_ancestry.R:
##   Rscript module_population_partitioning/R/pp_figures.R
## Writes: module_population_partitioning/Figures/fig{1-4}_*.png
## =========================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
DATADIR <- "module_population_partitioning/data"
FIGDIR  <- "module_population_partitioning/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

obj  <- readRDS(file.path(DATADIR, "pp_units_Fmat.rds"))
res  <- readRDS(file.path(DATADIR, "pp_concordance_results.rds"))
bb   <- readRDS(file.path(DATADIR, "pp_block_bootstrap.rds"))
nullc<- readRDS(file.path(DATADIR, "pp_null_check.rds"))
rez  <- readRDS(file.path(DATADIR, "pp_residual_ancestry.rds"))
u <- res$u; Fmat <- obj$Fmat; setDT(u)

AQU <- "#21918C"; POL <- "#D3C93B"
COL_SIGNED <- "#1b9e77"; COL_ABS <- "#d95f02"; COL_RESID <- "#7570b3"
theme_ms <- theme_bw(base_size = 12) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank())

## empirical cross-chromosome baseline, with an (approximate, pair-level)
## bootstrap band -- the AUDIT-recommended PRIMARY reference (preserves real
## among-population covariance, unlike the label-permutation null)
set.seed(99)
cross_boot_r    <- vapply(1:1000, function(b) mean(sample(nullc$cross_chr, replace = TRUE), na.rm = TRUE), numeric(1))
cross_boot_absr <- vapply(1:1000, function(b) mean(abs(sample(nullc$cross_chr, replace = TRUE)), na.rm = TRUE), numeric(1))
CROSS_R    <- mean(nullc$cross_chr, na.rm = TRUE); CROSS_R_LO <- quantile(cross_boot_r, 0.025); CROSS_R_HI <- quantile(cross_boot_r, 0.975)
CROSS_ABSR <- mean(abs(nullc$cross_chr), na.rm = TRUE); CROSS_ABSR_LO <- quantile(cross_boot_absr, 0.025); CROSS_ABSR_HI <- quantile(cross_boot_absr, 0.975)
PERM_FLOOR <- nullc$null_floor_absr   # secondary reference: pure small-n sampling-noise floor for |r|

## =========================================================================
## Fig 1: population x unit heatmap, Chr26, + zoomed panel on the named
## F7174 block
## =========================================================================
CHR <- "Chr26"
idx <- which(u$Chr == CHR)
sub <- Fmat[, idx, drop = FALSE]
pop_mean <- rowMeans(Fmat, na.rm = TRUE)
pop_ord  <- names(sort(pop_mean, decreasing = TRUE))
dm <- as.data.table(sub, keep.rownames = "pop")
dm <- melt(dm, id.vars = "pop", variable.name = "group_id", value.name = "f_aqu")
dm[, group_id := as.character(group_id)]
dm <- merge(dm, u[idx, .(group_id, Pos, FST, sort_class, n_loci_g)], by = "group_id")
dm[, pop := factor(pop, levels = rev(pop_ord))]
unit_ord <- u[idx][order(Pos), group_id]
dm[, group_id := factor(group_id, levels = unit_ord)]

ann <- u[idx][order(Pos)]; ann[, group_id := factor(group_id, levels = unit_ord)]
p_ann <- ggplot(ann, aes(group_id, 1, fill = FST)) + geom_tile() +
  scale_fill_gradient(low = "grey90", high = "firebrick", name = "FST") +
  theme_void() + theme(legend.position = "top")
p_hm <- ggplot(dm, aes(group_id, pop, fill = f_aqu)) + geom_tile() +
  scale_fill_gradient2(low = POL, mid = "grey95", high = AQU, midpoint = 0.5,
                       name = "oriented\naquilonia\nfreq.", limits = c(0, 1)) +
  labs(x = sprintf("%s units, ordered by position (n=%d, span %.2f-%.2f Mb)",
                   CHR, length(idx), min(ann$Pos)/1e6, max(ann$Pos)/1e6), y = NULL) +
  theme_ms + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                   panel.grid = element_blank())

## zoomed panel: the named F7174 block (module_di25/data/di25_three_blocks.rds) +/- flanking units
blk <- readRDS("module_di25/data/di25_three_blocks.rds")
blk_ids <- unlist(strsplit(blk[chr == 26]$group_ids, ","))
blk_pos <- u[group_id %in% blk_ids, range(Pos)]
zoom_idx <- which(u$Chr == CHR & u$Pos >= blk_pos[1] - 3e5 & u$Pos <= blk_pos[2] + 3e5)
zsub <- Fmat[, zoom_idx, drop = FALSE]
zdm <- as.data.table(zsub, keep.rownames = "pop")
zdm <- melt(zdm, id.vars = "pop", variable.name = "group_id", value.name = "f_aqu")
zdm[, group_id := as.character(group_id)]
zann <- u[zoom_idx][order(Pos)]
zann[, is_block := group_id %in% blk_ids]
zorder <- zann$group_id
zdm <- merge(zdm, zann[, .(group_id, Pos, FST, n_loci_g, is_block)], by = "group_id")
zdm[, pop := factor(pop, levels = rev(pop_ord))]
zdm[, xlab := sprintf("%s\n(%.2fMb, n_loci=%d)%s", group_id, Pos/1e6, n_loci_g, ifelse(is_block, "*", ""))]
xlab_ord <- zann[, sprintf("%s\n(%.2fMb, n_loci=%d)%s", group_id, Pos/1e6, n_loci_g, ifelse(is_block, "*", ""))]
zdm[, xlab := factor(xlab, levels = xlab_ord)]
p_zoom <- ggplot(zdm, aes(xlab, pop, fill = f_aqu)) + geom_tile() +
  scale_fill_gradient2(low = POL, mid = "grey95", high = AQU, midpoint = 0.5, limits = c(0, 1), guide = "none") +
  labs(x = sprintf("F7174 block +/- 300kb flanking (* = in the named block, n=%d)", length(blk_ids)), y = NULL,
      title = "Zoom: named F7174 polyctena block (Chr26)") +
  theme_ms + theme(axis.text.x = element_text(size = 7), panel.grid = element_blank())

fig1 <- p_ann + p_hm + p_zoom + plot_layout(ncol = 1, heights = c(1, 12, 9))
ggsave(file.path(FIGDIR, "fig1_heatmap_Chr26.png"), fig1, width = 11, height = 9, dpi = 200)
cat("[fig1] saved fig1_heatmap_Chr26.png (whole chromosome + zoomed named block)\n")

## =========================================================================
## Fig 2: local similarity vs physical distance -- signed r primary,
## |r| secondary; chromosome-block bootstrap 95% CI; empirical
## cross-chromosome baseline as the main reference band.
## =========================================================================
dci <- copy(bb$dist_ci); dci[, dbin := factor(dbin, levels = dbin)]
fig2a <- ggplot(dci, aes(dbin, mean_r, group = 1)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = CROSS_R_LO, ymax = CROSS_R_HI, fill = COL_SIGNED, alpha = 0.15) +
  geom_hline(yintercept = CROSS_R, linetype = 2, colour = COL_SIGNED) +
  annotate("text", x = dci$dbin[1], y = CROSS_R, label = "empirical cross-chromosome baseline",
           hjust = 0, vjust = -0.6, size = 3.1, colour = COL_SIGNED) +
  geom_line(colour = COL_SIGNED) + geom_point(size = 2.2, colour = COL_SIGNED) +
  geom_errorbar(aes(ymin = lo_r, ymax = hi_r), width = 0.1, colour = COL_SIGNED) +
  labs(x = NULL, y = "mean SIGNED r (primary)") + theme_ms
fig2b <- ggplot(dci, aes(dbin, mean_absr, group = 1)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = CROSS_ABSR_LO, ymax = CROSS_ABSR_HI, fill = COL_ABS, alpha = 0.15) +
  geom_hline(yintercept = CROSS_ABSR, linetype = 2, colour = COL_ABS) +
  geom_hline(yintercept = PERM_FLOOR, linetype = 3, colour = "grey50") +
  annotate("text", x = dci$dbin[1], y = PERM_FLOOR, label = "label-permutation sampling floor",
           hjust = 0, vjust = 1.4, size = 3, colour = "grey40") +
  geom_line(colour = COL_ABS) + geom_point(size = 2.2, colour = COL_ABS) +
  geom_errorbar(aes(ymin = lo_absr, ymax = hi_absr), width = 0.1, colour = COL_ABS) +
  labs(x = "physical distance between unit pair", y = "mean |r| (secondary)") + theme_ms
fig2 <- fig2a / fig2b
ggsave(file.path(FIGDIR, "fig2_similarity_vs_distance.png"), fig2, width = 7.5, height = 7.5, dpi = 200)
cat("[fig2] saved fig2_similarity_vs_distance.png (signed primary + |r| secondary, block-bootstrap CI)\n")

## =========================================================================
## Fig 3: FST vs local partition similarity -- signed r primary, |r|
## secondary; block-bootstrap CI on decile means; explicitly NOT attributed
## to the 48 large clusters (verified negligible effect, see README).
## =========================================================================
fci <- copy(bb$fst_ci)
fig3a <- ggplot(u[!is.na(FST) & !is.na(near_r)], aes(FST, near_r)) +
  geom_point(alpha = 0.05, size = 0.6, colour = "grey60") +
  geom_hline(yintercept = CROSS_R, linetype = 2, colour = COL_SIGNED) +
  geom_point(data = fci, aes(mean_FST, mean_r), colour = COL_SIGNED, size = 2.3, inherit.aes = FALSE) +
  geom_line(data = fci, aes(mean_FST, mean_r), colour = COL_SIGNED, inherit.aes = FALSE) +
  geom_errorbar(data = fci, aes(mean_FST, ymin = lo_r, ymax = hi_r), colour = COL_SIGNED, width = 0, inherit.aes = FALSE) +
  labs(x = expression(F[ST]~"(unit)"), y = "mean local SIGNED r, <=100kb (primary)") + theme_ms
fig3b <- ggplot(fci, aes(mean_FST, mean_absr)) +
  geom_hline(yintercept = CROSS_ABSR, linetype = 2, colour = COL_ABS) +
  geom_point(size = 2.3, colour = COL_ABS) + geom_line(colour = COL_ABS) +
  geom_errorbar(aes(ymin = lo_absr, ymax = hi_absr), width = 0, colour = COL_ABS) +
  labs(x = expression(F[ST]~"decile mean"), y = "mean local |r|, <=100kb (secondary)\ndecile mean +/- block-bootstrap 95% CI") + theme_ms
fig3 <- fig3a + fig3b
ggsave(file.path(FIGDIR, "fig3_FST_vs_similarity.png"), fig3, width = 11, height = 5, dpi = 200)
cat(sprintf("[fig3] saved fig3_FST_vs_similarity.png (rho_r=%.3f, rho_absr=%.3f; block-bootstrap slope CI excludes 0 but is small)\n",
            cor(u$FST, u$near_r, use="pairwise.complete.obs", method="spearman"),
            cor(u$FST, u$near_absr, use="pairwise.complete.obs", method="spearman")))

## =========================================================================
## Fig 4 (REPLACED): distance-decay of SIGNED r, raw profile vs residualized
## against leave-one-chromosome-out genome-wide ancestry, overlaid. This is
## the clearest view of the residualization result: short-range concordance
## (<~100kb) is essentially unchanged (genuine locus-specific signal), while
## the raw curve's slowly-decaying long-range "floor" collapses to ~0 once
## genome-wide ancestry is removed -- i.e. that floor WAS mostly shared
## ancestry, not a locus-specific signal. A genome-position panel is added
## below for the residualized profile alone, coloured by FST, since with the
## long-range floor removed it is now legible (unlike the raw near-100kb
## scatter, which looks the same before/after residualizing and was dropped).
## =========================================================================
dbr <- copy(rez$pairs_resid_summary)[!is.na(dbin)]
dbr[, dbin := factor(dbin, levels = levels(dci$dbin))]
dcomp <- rbindlist(list(
  dci[, .(dbin, mean_r, lo = lo_r, hi = hi_r, which = "raw profile")],
  dbr[, .(dbin, mean_r, lo = NA_real_, hi = NA_real_, which = "residualized (genome-wide ancestry removed)")]
))
fig4a <- ggplot(dcomp, aes(dbin, mean_r, colour = which, group = which)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.1, na.rm = TRUE) +
  scale_colour_manual(values = c("raw profile" = COL_SIGNED, "residualized (genome-wide ancestry removed)" = COL_RESID), name = NULL) +
  labs(x = "physical distance between unit pair", y = "mean SIGNED r",
      title = "Short-range concordance survives residualizing out genome-wide ancestry;\nthe raw curve's long-range floor does not") +
  theme_ms + theme(legend.position = "top")

## residualized-profile local signed r along the genome, coloured by FST (now
## legible since the long-range shared-ancestry floor has been removed)
ru <- rez$u; setDT(ru)
ru[, ChrNum := as.integer(sub("Chr", "", Chr))]; setorder(ru, ChrNum, Pos)
chr_len <- ru[, .(len = max(Pos)), by = .(Chr, ChrNum)][order(ChrNum)]
chr_len[, offset := cumsum(shift(len, fill = 0)) + cumsum(shift(rep(2e6, .N), fill = 0))]
ru2 <- merge(ru, chr_len[, .(Chr, offset)], by = "Chr"); ru2[, gpos := Pos + offset]; setorder(ru2, gpos)
chr_mid <- merge(chr_len, ru2[, .(mn = min(gpos), mx = max(gpos)), by = Chr], by = "Chr")
chr_mid[, mid := (mn + mx) / 2]; setorder(chr_mid, ChrNum)
fig4b <- ggplot(ru2[!is.na(near_r_resid)], aes(gpos, near_r_resid, colour = FST)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
  geom_point(size = 0.6, alpha = 0.5) +
  scale_colour_gradient(low = "grey80", high = "firebrick", name = expression(F[ST])) +
  scale_x_continuous(breaks = chr_mid$mid, labels = chr_mid$ChrNum, expand = c(0.01, 0.01)) +
  labs(x = "chromosome", y = "residual local SIGNED r, <=100kb window") +
  theme_ms + theme(panel.grid.major.x = element_blank())

fig4 <- fig4a / fig4b + plot_layout(heights = c(5, 4))
ggsave(file.path(FIGDIR, "fig4_genomewide_summary.png"), fig4, width = 11, height = 8, dpi = 200)
cat("[fig4] saved fig4_genomewide_summary.png (distance-decay raw vs residualized + residual genome-wide panel, replaces the original Fig 4)\n")

cat("\nAll figures saved to", FIGDIR, "\n")
