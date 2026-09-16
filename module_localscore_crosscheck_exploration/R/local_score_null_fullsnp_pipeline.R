## =========================================================
## module_localscore_crosscheck -- FULL-SNP-RESOLUTION false-positive
## demonstration: worst-of-20-random-draws null covariates, run through
## the ACTUAL full-SNP BayPass pipeline (not Stage-1-cluster-inherited
## values), then the same local-score + Stage-2-styled Manhattan
## treatment as the real observed data (Figure fig:obs-fullsnp).
## =========================================================
## Picks: 20 draws were randomly sampled (seed 42) from each null pool
## (continuous, shared by PC1/PC2/bio_winter; mitoC2's own dedicated
## pool), using the window counts already scanned by
## local_score_null_check_S1units.R (Stage-1-cluster resolution, the only
## resolution with 1,000 pre-scanned draws to sample from). The worst of
## each 20 -- continuous-null draw #74 (1 window). mitoC2-null draw #212
## (9 windows) -- was then run through the SAME full-SNP BayPass
## invocation used for the real covariates (run_baypass_fullsnp_null_
## examples.sh: identical geno/Omega/poolsize/MCMC settings, only the
## -efile/-contrastfile swapped for the null draw's raw covariate values,
## pulled from mini2's null_b01.env row 74 / nullc2_b02.contrast row 12).
##
## This closes the full-SNP-resolution gap flagged in Section~8/limitations:
## previously no null had ever been run through full-SNP BayPass (10,000
## draws at that scale was never feasible); this is the first and only one.
##
## Reads : module_manuscript_rho05/baypass_stage1/aland_excluded/
##           nullcont_draw74_fullSNP_stage1Omega_summary_betai_reg.out
##           nullmito_draw212_fullSNP_stage1Omega_summary_contrast.out
## Writes: module_localscore_crosscheck/data/localscore_null_fullsnp_<tag>.rds
##         module_localscore_crosscheck/Figures/local_score_null_fullsnp_examples.{png,pdf}
##
## Run from the repo root:
##   Rscript module_localscore_crosscheck/R/local_score_null_fullsnp_pipeline.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })
devtools::load_all("~/gitlab/LDscnR/")
source("~/gitlab/baypass_public-master/utils/baypass_utils.R")

BP_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded"
DATA   <- "module_localscore_crosscheck/data"
FIGDIR <- "module_localscore_crosscheck/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

load("data/hybrids_only_maf005.Rdata")
pos <- data.frame(Chr = map_hyb_005$Chr, Pos = map_hyb_005$Pos)

s2 <- readRDS("module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05.rds")
g2 <- as.data.table(s2$groups)
marker2s2 <- g2[, .(marker = unlist(members)), by = .(s2_group = group_id)]
setkey(marker2s2, marker)

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

## ---- 1. run compute.local.scores() on the genuine full-SNP null draws ----
pi_pc1 <- fread(file.path(BP_DIR, "PC1_fullSNP_stage1Omega_summary_beta_params.out"))   # not Pi; use pi_xtx below
pi_pc1 <- fread("module_manuscript_rho05/baypass_stage1/aland_excluded/PC1_fullSNP_stage1Omega_summary_pi_xtx.out", select = c("MRK", "M_P"))$M_P
pi_c2  <- fread("module_manuscript_rho05/baypass_stage1/aland_excluded/mito_C2_fullSNP_stage1Omega_summary_pi_xtx.out", select = c("MRK", "M_P"))$M_P

set.seed(1)   # reproducibility, see local_score_regions.R header
cont_bf <- fread(file.path(BP_DIR, "nullcont_draw74_fullSNP_stage1Omega_summary_betai_reg.out"), select = c("MRK", "BF(dB)"))
stopifnot(all(cont_bf$MRK == seq_len(nrow(cont_bf))))
res_cont <- compute.local.scores(snp.position = pos, snp.pi = pi_pc1, snp.bf = cont_bf$`BF(dB)`,
                                 xi = 1, pval.local.score.thres = 0.01, min.maf = 0.2, min.nsnp = 1e4)

mito_p <- fread(file.path(BP_DIR, "nullmito_draw212_fullSNP_stage1Omega_summary_contrast.out"), select = c("MRK", "log10(1/pval)"))
stopifnot(all(mito_p$MRK == seq_len(nrow(mito_p))))
res_mito <- compute.local.scores(snp.position = pos, snp.pi = pi_c2, snp.pvalue = mito_p$`log10(1/pval)`,
                                 xi = 1, pval.local.score.thres = 0.01, min.maf = 0.2, min.nsnp = 1e4)

saveRDS(res_cont, file.path(DATA, "localscore_null_fullsnp_continuous_draw74.rds"))
saveRDS(res_mito, file.path(DATA, "localscore_null_fullsnp_mitoC2_draw212.rds"))

n_win_cont <- if (is.null(res_cont$significant.windows)) 0L else nrow(res_cont$significant.windows)
n_win_mito <- if (is.null(res_mito$significant.windows)) 0L else nrow(res_mito$significant.windows)
cat(sprintf("[full-SNP null] continuous draw #74: %d significant window(s)\n", n_win_cont))
cat(sprintf("[full-SNP null] mitoC2 draw #212: %d significant window(s)\n", n_win_mito))

## ---- 2. manuscript-styled panel (identical logic/arrow convention to
## local_score_snp_manhattan_manuscript_style.R) ------------------------
make_panel <- function(tag, draw_label, stat_file, stat_col, thresh, y_lab, win, is_bf) {
  n_win <- if (is.null(win) || nrow(win) == 0) 0L else nrow(win)
  win <- as.data.table(win)

  s <- fread(stat_file)
  stopifnot(nrow(s) == nrow(map_hyb_005), all(s$MRK == seq_len(nrow(s))))
  snp_dt <- map_hyb_005[, .(marker, Chr, Pos)]
  snp_dt[, stat := s[[stat_col]]]
  add_gpos(snp_dt)
  snp_dt[chr_lens, on = "Chr", band := i.band]
  setkey(snp_dt, Chr, Pos)

  raw <- snp_dt[stat >= thresh]
  raw[marker2s2, on = .(marker), s2_group := i.s2_group]
  raw_regions <- sort(unique(raw$s2_group[!is.na(raw$s2_group)]))

  arrows_dt <- data.table(win_id = character(0), s2_group = character(0), gpos_mid = numeric(0))
  if (n_win > 0) {
    peak_col <- if (is_bf) "BF (dB) peak pos" else "-log10(p-val) peak pos"
    win[, peak_pos := as.integer(get(peak_col))]
    win[, mid_pos := as.integer(round((beg + end) / 2))]
    win[, win_id := paste0(tag, "_W", seq_len(.N))]
    peak_marker <- snp_dt[.(Chr = win$chr, Pos = win$peak_pos), on = .(Chr, Pos), marker]
    win[, peak_s2 := marker2s2[.(peak_marker), on = "marker", s2_group]]
    win[, gpos_mid := mid_pos + chr_lens$offset[match(chr, chr_lens$Chr)]]
    arrows_dt <- win[!is.na(peak_s2), .(win_id, s2_group = peak_s2, gpos_mid)]
  }

  all_regions <- sort(union(raw_regions, unique(arrows_dt$s2_group)))
  region_cols <- setNames(rep(PAL, length.out = length(all_regions)), all_regions)
  region_snps <- if (length(all_regions) > 0) g2[group_id %in% all_regions,
                                      .(marker = unlist(members)), by = .(s2_group = group_id)] else
                                      data.table(marker = character(0), s2_group = character(0))
  snp_dt[region_snps, on = "marker", s2_group := i.s2_group]
  snp_dt[, is_region := !is.na(s2_group)]

  y_top <- max(snp_dt$stat, na.rm = TRUE)
  arrow_y_head <- y_top * 1.35
  arrow_y_tail <- y_top * 1.55

  message(sprintf("[%s null draw %s, FULL-SNP] %d raw crossings -> %d Stage-2 regions coloured; %d local-score window(s) -> %d window(s) arrowed",
                  tag, draw_label, nrow(raw), length(all_regions), n_win, nrow(arrows_dt)))

  ggplot() +
    geom_point(data = snp_dt[is_region == FALSE & band == FALSE], aes(gpos, stat), colour = "grey75", size = 0.3, alpha = 0.6) +
    geom_point(data = snp_dt[is_region == FALSE & band == TRUE], aes(gpos, stat), colour = "grey50", size = 0.3, alpha = 0.6) +
    geom_point(data = snp_dt[is_region == TRUE], aes(gpos, stat, colour = s2_group), size = 0.9) +
    { if (nrow(arrows_dt) > 0) geom_label_repel(
        data = arrows_dt, aes(x = gpos_mid, y = arrow_y_head, label = s2_group, colour = s2_group),
        fill = "white", fontface = "bold", size = 2.6, label.padding = unit(0.12, "lines"),
        nudge_y = arrow_y_tail - arrow_y_head, direction = "y",
        arrow = arrow(length = unit(0.14, "cm"), type = "closed"), segment.colour = "black", segment.size = 0.5,
        box.padding = 0.3, point.padding = 0.1, min.segment.length = 0, max.overlaps = Inf, seed = 1,
        show.legend = FALSE) } +
    geom_hline(yintercept = thresh, linetype = 2, colour = "black", linewidth = 0.3) +
    scale_colour_manual(values = region_cols, guide = "none") +
    scale_x_continuous(breaks = chr_lens$mid, labels = chr_lens$chr_num, expand = c(0.005, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.45))) +
    labs(x = "Chromosome", y = y_lab,
         title = sprintf("%s -- NULL DRAW %s, FULL-SNP BayPass (NOT REAL DATA)", tag, draw_label)) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 11, face = "bold", hjust = 0, colour = "#B03A2E"))
}

p_cont <- make_panel("PC1/PC2/bio_winter pool", "#74", file.path(BP_DIR, "nullcont_draw74_fullSNP_stage1Omega_summary_betai_reg.out"),
                     "BF(dB)", 15, "BF (dB)", res_cont$significant.windows, is_bf = TRUE)
p_mito <- make_panel("mitoC2 pool", "#212", file.path(BP_DIR, "nullmito_draw212_fullSNP_stage1Omega_summary_contrast.out"),
                     "log10(1/pval)", 3, expression(C2~-log[10](p)), res_mito$significant.windows, is_bf = FALSE)

combined <- p_cont / p_mito
outpng <- file.path(FIGDIR, "local_score_null_fullsnp_examples.png")
outpdf <- file.path(FIGDIR, "local_score_null_fullsnp_examples.pdf")
ggsave(outpng, combined, width = 12, height = 8, dpi = 300, limitsize = FALSE)
ggsave(outpdf, combined, width = 12, height = 8, limitsize = FALSE)
cat("wrote", outpng, "and", outpdf, "\n")
