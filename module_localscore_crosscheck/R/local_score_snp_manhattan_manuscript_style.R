## =========================================================
## module_localscore_crosscheck -- local-score outlier regions, OBSERVED data,
## full-SNP resolution, styled EXACTLY like the manuscript's Figure 5
## (moduleB_stage1_snp_manhattan_combined.R)
## =========================================================
## Figure 5 (the real floor-survivor combined figure) plots the genuine
## full-SNP BF/C2 value for every SNP, but decides "raw crossing" /
## "significant" at STAGE-1-CLUSTER resolution (18,361 units) -- the
## full-SNP scan is used only for the Y-axis display, not the test itself.
## This script instead runs the significance test (local-score) AT the
## same full-SNP resolution the points are drawn at: every individual
## raw-crossing SNP (not just Stage-1 cluster representatives) is mapped to
## its canonical Stage-2 (rho05) region for COLOURING (all touched regions).
##
## ARROWS ARE PER LOCAL-SCORE WINDOW, NOT PER TOUCHED STAGE-2 REGION. A
## first version of this script arrowed every Stage-2 region touching a
## significant window and got 138-910 arrows per covariate (8,055 for
## mitoC2) -- because local-score windows are often hundreds of kb wide and
## routinely span dozens of the much finer-grained (rho=0.5, 0.5cM cap)
## Stage-2 clusters. That inflates the apparent candidate count and isn't
## what the method actually claims: compute.local.scores() calls a WINDOW
## significant, not every pre-existing LD cluster it happens to overlap.
## Each window gets exactly one arrow, positioned at its midpoint and
## labelled/coloured by the Stage-2 group of its own peak SNP (the
## strongest raw signal inside that window) -- the closest faithful
## analogue of Figure 5's one-arrow-per-floor-survivor-region convention.
##
## Reads : module_localscore_crosscheck/data/localscore_<tag>.rds
##         (from local_score_regions.R; full-SNP resolution)
##         module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05.rds
##         module_manuscript_rho05/baypass_stage1/aland_excluded/
##           <tag>_fullSNP_stage1Omega_summary_{betai_reg,contrast}.out
## Writes: module_localscore_crosscheck/Figures/local_score_snp_manhattan_combined.{png,pdf}
##
## Run from the repo root:
##   Rscript module_localscore_crosscheck/R/local_score_snp_manhattan_manuscript_style.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })
devtools::load_all("~/gitlab/LDscnR/")

BP_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded"
DATA   <- "module_localscore_crosscheck/data"
FIGDIR <- "module_localscore_crosscheck/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

load("data/hybrids_only_maf005.Rdata")

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

make_panel <- function(tag, stat_file, stat_col, thresh, y_lab, is_bf) {
  ls_res <- readRDS(file.path(DATA, sprintf("localscore_%s.rds", tag)))
  win <- as.data.table(ls_res$significant.windows)
  n_win <- if (is.null(win) || nrow(win) == 0) 0L else nrow(win)

  s <- fread(stat_file)
  stopifnot(nrow(s) == nrow(map_hyb_005), all(s$MRK == seq_len(nrow(s))))
  snp_dt <- map_hyb_005[, .(marker, Chr, Pos)]
  snp_dt[, stat := s[[stat_col]]]
  add_gpos(snp_dt)
  snp_dt[chr_lens, on = "Chr", band := i.band]
  setkey(snp_dt, Chr, Pos)

  ## every INDIVIDUAL raw-crossing SNP -> its Stage-2 region (full-SNP
  ## resolution "raw", unlike Figure 5's Stage-1-cluster-level raw)
  raw <- snp_dt[stat >= thresh]
  raw[marker2s2, on = .(marker), s2_group := i.s2_group]
  raw_regions <- sort(unique(raw$s2_group[!is.na(raw$s2_group)]))

  ## one arrow per WINDOW (not per touched Stage-2 region -- see header),
  ## positioned at the window midpoint, labelled/coloured by the Stage-2
  ## group of the window's OWN peak SNP (the peak-pos columns compute.local.
  ## scores() reports are genomic positions, not marker names -- resolved
  ## back to a marker via an exact Chr+Pos join, then to its Stage-2 group).
  ## NOTE: a window's peak SNP is the strongest RAW value inside the window,
  ## but the whole point of the local-score method is detecting windows
  ## where evidence accumulates WITHOUT any single SNP crossing the raw
  ## threshold -- so the peak SNP is often not itself a raw crossing, and
  ## its Stage-2 group may not be one of `raw_regions`. Looked up directly
  ## from the canonical clustering (marker2s2), independent of raw_regions,
  ## and unioned in below so every arrowed region still gets its own colour
  ## even when none of its members individually crossed the raw threshold.
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

  message(sprintf("[%s] %d raw-crossing SNPs -> %d Stage-2 regions coloured (%d from raw crossings, %d added for arrow-only windows); %d local-score window(s) -> %d window(s) arrowed",
                  tag, nrow(raw), length(all_regions), length(raw_regions),
                  length(all_regions) - length(raw_regions), n_win, nrow(arrows_dt)))

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
    labs(x = "Chromosome", y = y_lab, title = tag) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 11, face = "bold", hjust = 0))
}

p_pc1 <- make_panel("PC1", file.path(BP_DIR, "PC1_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)", 15, "BF (dB)", is_bf = TRUE)
p_pc2 <- make_panel("PC2", file.path(BP_DIR, "PC2_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)", 15, "BF (dB)", is_bf = TRUE)
p_bw  <- make_panel("bio_winter", file.path(BP_DIR, "bio_winter_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)", 15, "BF (dB)", is_bf = TRUE)
p_c2  <- make_panel("mitoC2", file.path(BP_DIR, "mito_C2_fullSNP_stage1Omega_summary_contrast.out"), "log10(1/pval)", 3, expression(C2~-log[10](p)), is_bf = FALSE)

combined <- (p_pc1 / p_pc2 / p_bw / p_c2)

outpng <- file.path(FIGDIR, "local_score_snp_manhattan_combined.png")
outpdf <- file.path(FIGDIR, "local_score_snp_manhattan_combined.pdf")
ggsave(outpng, combined, width = 12, height = 16, dpi = 300, limitsize = FALSE)
ggsave(outpdf, combined, width = 12, height = 16, limitsize = FALSE)
cat("wrote", outpng, "and", outpdf, "\n")
