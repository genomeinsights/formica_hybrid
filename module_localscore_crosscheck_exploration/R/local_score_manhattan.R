## =========================================================
## module_localscore_crosscheck -- Manhattan plot of BayPass's own local-score
## outlier regions (Fariello et al. 2017 / Bonhomme et al. 2019), all four
## covariates
## =========================================================
## Visualises the significant.windows produced by
## local_score_regions.R (compute.local.scores(), BayPass's
## own defaults: xi=1, pval.local.score.thres=0.01, min.maf=0.2,
## min.nsnp=1e4) against the full-SNP genome-wide scan. Every SNP whose
## position falls INSIDE a significant local-score window is coloured (one
## colour per window, recycled palette, alternating grey bands elsewhere)
## -- this is the analytic-threshold, spatial-profile-based region
## definition, NOT the Stage-2 LD-cluster + null-calibrated floor-survivor
## regions plotted in moduleB_stage1_snp_manhattan_combined.R. Exploratory
## comparison figure, not yet manuscript-styled (keeps per-panel counts in
## the subtitle rather than stripping them).
##
## Reads : module_localscore_crosscheck/data/localscore_<tag>.rds
##         (from local_score_regions.R)
##         module_manuscript_rho05/baypass_stage1/aland_excluded/
##           <tag>_fullSNP_stage1Omega_summary_{betai_reg,contrast}.out
## Writes: module_localscore_crosscheck/Figures/local_score_manhattan.{png,pdf}
##
## Run from the repo root:
##   Rscript module_localscore_crosscheck/R/local_score_manhattan.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
devtools::load_all("~/gitlab/LDscnR/")

BP_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded"
DATA   <- "module_localscore_crosscheck/data"
FIGDIR <- "module_localscore_crosscheck/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

load("data/hybrids_only_maf005.Rdata")   # map_hyb_005 (Chr, Pos, marker), full-SNP MRK order

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
  ls_res <- readRDS(file.path(DATA, sprintf("localscore_%s.rds", tag)))
  win <- as.data.table(ls_res$significant.windows)
  n_win <- if (is.null(win) || nrow(win) == 0) 0L else nrow(win)
  if (n_win > 0) win[, win_id := paste0(tag, "_W", seq_len(.N))]

  s <- fread(stat_file)
  stopifnot(nrow(s) == nrow(map_hyb_005), all(s$MRK == seq_len(nrow(s))))
  snp_dt <- map_hyb_005[, .(marker, Chr, Pos)]
  snp_dt[, stat := s[[stat_col]]]
  add_gpos(snp_dt)
  snp_dt[chr_lens, on = "Chr", band := i.band]
  snp_dt[, win_id := NA_character_]
  if (n_win > 0) for (i in seq_len(nrow(win))) {
    idx <- snp_dt$Chr == win$chr[i] & snp_dt$Pos >= win$beg[i] & snp_dt$Pos <= win$end[i]
    snp_dt[idx, win_id := win$win_id[i]]
  }
  snp_dt[, is_win := !is.na(win_id)]

  win_cols <- if (n_win > 0) setNames(rep(PAL, length.out = n_win), win$win_id) else character(0)
  message(sprintf("[%s] %d significant window(s), %d SNPs coloured", tag, n_win, sum(snp_dt$is_win)))

  ggplot() +
    geom_point(data = snp_dt[is_win == FALSE & band == FALSE], aes(gpos, stat), colour = "grey75", size = 0.3, alpha = 0.6) +
    geom_point(data = snp_dt[is_win == FALSE & band == TRUE], aes(gpos, stat), colour = "grey50", size = 0.3, alpha = 0.6) +
    geom_point(data = snp_dt[is_win == TRUE], aes(gpos, stat, colour = win_id), size = 0.9) +
    geom_hline(yintercept = thresh, linetype = 2, colour = "black", linewidth = 0.3) +
    scale_colour_manual(values = win_cols, guide = "none") +
    scale_x_continuous(breaks = chr_lens$mid, labels = chr_lens$chr_num, expand = c(0.005, 0)) +
    labs(x = "Chromosome", y = y_lab, title = tag,
         subtitle = sprintf("%d local-score significant window(s), BayPass defaults (xi=1, alpha=0.01)", n_win)) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 11, face = "bold", hjust = 0),
          plot.subtitle = element_text(size = 8, colour = "grey30"))
}

p_pc1 <- make_panel("PC1", file.path(BP_DIR, "PC1_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)", 15, "BF (dB)")
p_pc2 <- make_panel("PC2", file.path(BP_DIR, "PC2_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)", 15, "BF (dB)")
p_bw  <- make_panel("bio_winter", file.path(BP_DIR, "bio_winter_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)", 15, "BF (dB)")
p_c2  <- make_panel("mitoC2", file.path(BP_DIR, "mito_C2_fullSNP_stage1Omega_summary_contrast.out"), "log10(1/pval)", 3, expression(C2~-log[10](p)))

combined <- (p_pc1 / p_pc2 / p_bw / p_c2)

outpng <- file.path(FIGDIR, "local_score_manhattan.png")
outpdf <- file.path(FIGDIR, "local_score_manhattan.pdf")
ggsave(outpng, combined, width = 12, height = 16, dpi = 300, limitsize = FALSE)
ggsave(outpdf, combined, width = 12, height = 16, limitsize = FALSE)
cat("wrote", outpng, "and", outpdf, "\n")
