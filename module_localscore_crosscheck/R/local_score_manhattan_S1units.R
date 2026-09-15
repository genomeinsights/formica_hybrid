## =========================================================
## module_localscore_crosscheck -- Manhattan plot of BayPass's local-score
## outlier regions computed at STAGE-1-CLUSTER resolution (18,361 units)
## =========================================================
## Companion to local_score_manhattan.R (full-SNP
## resolution) and local_score_regions_S1units.R (the
## analysis this plots). As in moduleB_stage1_region_manhattan.R, every
## genome-wide member SNP of a TESTED Stage-1 cluster (n_snps>=5) is
## plotted inheriting its own cluster's statistic; a SNP is coloured if its
## cluster's representative position (core_snp) falls inside a significant
## local-score window at Stage-1-cluster resolution (one colour per window,
## recycled palette).
##
## Reads : module_localscore_crosscheck/data/localscore_S1units_<tag>.rds
##         module0_ld_pruning/data/pruned_stage1.rds (cl5, members, core_snp)
##         module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/
##           <tag>_S1units_withOmega_summary_{betai_reg,contrast}.out
## Writes: module_localscore_crosscheck/Figures/local_score_manhattan_S1units.{png,pdf}
##
## Run from the repo root:
##   Rscript module_localscore_crosscheck/R/local_score_manhattan_S1units.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
devtools::load_all("~/gitlab/LDscnR/")

UNIT_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded_S1units"
DATA     <- "module_localscore_crosscheck/data"
FIGDIR   <- "module_localscore_crosscheck/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

load("data/hybrids_only_maf005.Rdata")   # map_hyb_005
stage1 <- readRDS("module0_ld_pruning/data/pruned_stage1.rds")
cl <- as.data.table(stage1$clusters)
cl5 <- cl[n_snps >= 5]
group_order <- readLines(file.path(UNIT_DIR, "S1units_group_order.txt"))
stopifnot(identical(group_order, paste0("S1_", cl5$CL_id)))
cl5[, group_id := group_order]
cl5[, core_pos := as.integer(sub(".*:", "", core_snp))]

chr_lens <- map_hyb_005[, .(len = max(Pos)), by = Chr]
chr_lens[, chr_num := as.integer(sub("Chr", "", Chr))]
setorder(chr_lens, chr_num)
chr_lens[, offset := cumsum(shift(len, fill = 0)) + (seq_len(.N) - 1) * 3e5]
chr_lens[, mid := offset + len / 2]
chr_lens[, band := seq_len(.N) %% 2 == 0]
add_gpos <- function(dt) dt[chr_lens, on = "Chr", `:=`(gpos = Pos + i.offset)]

PAL_ALL <- default_cluster_colours()
lum <- { rgb <- t(col2rgb(PAL_ALL)) / 255; 0.299 * rgb[, 1] + 0.587 * rgb[, 2] + 0.114 * rgb[, 3] }
PAL <- PAL_ALL[lum <= 0.75]

make_panel <- function(tag, stat_file, stat_col, thresh, y_lab) {
  ls_res <- readRDS(file.path(DATA, sprintf("localscore_S1units_%s.rds", tag)))
  win <- as.data.table(ls_res$significant.windows)
  n_win <- if (is.null(win) || nrow(win) == 0) 0L else nrow(win)
  if (n_win > 0) win[, win_id := paste0(tag, "_W", seq_len(.N))]

  s <- fread(stat_file)
  stopifnot(nrow(s) == nrow(cl5))
  cl5s <- copy(cl5)
  cl5s[, stat := s[[stat_col]][match(seq_len(.N), s$MRK)]]
  cl5s[, win_id := NA_character_]
  if (n_win > 0) for (i in seq_len(nrow(win))) {
    idx <- cl5s$Chr == win$chr[i] & cl5s$core_pos >= win$beg[i] & cl5s$core_pos <= win$end[i]
    cl5s[idx, win_id := win$win_id[i]]
  }

  snp_dt <- cl5s[, .(marker = unlist(members)), by = .(CL_id, stat, win_id)]
  snp_dt <- map_hyb_005[, .(marker, Chr, Pos)][snp_dt, on = "marker"]
  add_gpos(snp_dt)
  snp_dt[chr_lens, on = "Chr", band := i.band]
  snp_dt[, is_win := !is.na(win_id)]

  win_cols <- if (n_win > 0) setNames(rep(PAL, length.out = n_win), win$win_id) else character(0)
  message(sprintf("[%s] %d significant window(s) at Stage-1-cluster resolution, %d member SNPs coloured",
                  tag, n_win, sum(snp_dt$is_win)))

  ggplot() +
    geom_point(data = snp_dt[is_win == FALSE & band == FALSE], aes(gpos, stat), colour = "grey75", size = 0.3, alpha = 0.6) +
    geom_point(data = snp_dt[is_win == FALSE & band == TRUE], aes(gpos, stat), colour = "grey50", size = 0.3, alpha = 0.6) +
    geom_point(data = snp_dt[is_win == TRUE], aes(gpos, stat, colour = win_id), size = 0.9) +
    geom_hline(yintercept = thresh, linetype = 2, colour = "black", linewidth = 0.3) +
    scale_colour_manual(values = win_cols, guide = "none") +
    scale_x_continuous(breaks = chr_lens$mid, labels = chr_lens$chr_num, expand = c(0.005, 0)) +
    labs(x = "Chromosome", y = y_lab, title = tag,
         subtitle = sprintf("%d local-score window(s), Stage-1-cluster resolution (18,361 units; min.nsnp=100)", n_win)) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 11, face = "bold", hjust = 0),
          plot.subtitle = element_text(size = 8, colour = "grey30"))
}

p_pc1 <- make_panel("PC1", file.path(UNIT_DIR, "PC1_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15, "BF (dB)")
p_pc2 <- make_panel("PC2", file.path(UNIT_DIR, "PC2_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15, "BF (dB)")
p_bw  <- make_panel("bio_winter", file.path(UNIT_DIR, "bio_winter_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15, "BF (dB)")
p_c2  <- make_panel("mitoC2", file.path(UNIT_DIR, "mito_C2_S1units_summary_contrast.out"), "log10(1/pval)", 3, expression(C2~-log[10](p)))

combined <- (p_pc1 / p_pc2 / p_bw / p_c2)

outpng <- file.path(FIGDIR, "local_score_manhattan_S1units.png")
outpdf <- file.path(FIGDIR, "local_score_manhattan_S1units.pdf")
ggsave(outpng, combined, width = 12, height = 16, dpi = 300, limitsize = FALSE)
ggsave(outpdf, combined, width = 12, height = 16, limitsize = FALSE)
cat("wrote", outpng, "and", outpdf, "\n")
