## =========================================================
## module_manuscript_rho05 -- manuscript-ready combined Manhattan figure
## (PC1 / PC2 / mitoC2 / bio_winter, Stage-2-region-coloured, arrowed)
## =========================================================
## Same data/plotting logic as moduleB_stage1_snp_manhattan_by_region.R
## (every genome-wide SNP; alternating grey chromosome bands; every
## Stage-2 (rho05) region touched by a raw BF/C2 threshold crossing
## coloured with the recycled, near-white-filtered palette; floor-survivor
## regions additionally arrowed and labelled with their Stage-2 group_id),
## but with all descriptive title/subtitle/caption text stripped from the
## plot itself -- only a minimal per-panel covariate-name label remains.
## Every number normally printed in-plot (raw crossings, regions coloured,
## floor survivors, regions arrowed, per-region Stage-1-cluster counts) is
## reported in the LaTeX caption instead.
##
## Reads : identical inputs to moduleB_stage1_snp_manhattan_by_region.R.
## Writes: module_manuscript_rho05/Figures/moduleB_stage1_snp_manhattan_combined.{png,pdf}
##
## Run from the repo root:
##   Rscript module_manuscript_rho05/R/moduleB_stage1_snp_manhattan_combined.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })
devtools::load_all("~/gitlab/LDscnR/")

BP_DIR   <- "module_manuscript_rho05/baypass_stage1"
UNIT_DIR <- file.path(BP_DIR, "aland_excluded_S1units")
SNP_DIR  <- file.path(BP_DIR, "aland_excluded")
DATA_DIR <- "module_manuscript_rho05/data"
FIGDIR   <- "module_manuscript_rho05/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

## ---- shared inputs (identical to moduleB_stage1_snp_manhattan_by_region.R) ----
load("data/hybrids_only_maf005.Rdata")
stage1 <- readRDS("module0_ld_pruning/data/pruned_stage1.rds")
cl <- as.data.table(stage1$clusters)
cl5 <- cl[n_snps >= 5]
group_order <- readLines(file.path(UNIT_DIR, "S1units_group_order.txt"))
stopifnot(identical(group_order, paste0("S1_", cl5$CL_id)))
cl5[, group_id := group_order]

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

## ---- one covariate -> a STRIPPED panel + its caption-worthy stats ---------
process_one <- function(unit_stat_file, unit_stat_col, thresh, snp_stat_file, snp_stat_col,
                        tag, title_stat, null_file = NULL, null_flag_col = NULL) {
  s <- fread(unit_stat_file)
  stopifnot(nrow(s) == nrow(cl5))
  cl5s <- copy(cl5)
  cl5s[, stat := s[[unit_stat_col]][match(seq_len(.N), s$MRK)]]

  raw <- cl5s[stat >= thresh]
  raw[marker2s2, on = .(core_snp = marker), s2_group := i.s2_group]
  all_regions <- sort(unique(raw$s2_group))

  nullobj <- readRDS(null_file)
  floor_ids <- nullobj$group_id[as.logical(nullobj[[null_flag_col]])]
  floor_units <- cl5s[group_id %in% floor_ids]
  floor_units[marker2s2, on = .(core_snp = marker), s2_group := i.s2_group]
  region_counts <- floor_units[, .N, by = s2_group]
  setorder(region_counts, -N)
  n_region <- nrow(region_counts)

  region_cols <- setNames(rep(PAL, length.out = length(all_regions)), all_regions)
  region_snps <- if (length(all_regions) > 0) g2[group_id %in% all_regions,
                                      .(marker = unlist(members)), by = .(s2_group = group_id)] else
                                      data.table(marker = character(0), s2_group = character(0))

  snp_s <- fread(snp_stat_file)
  stopifnot(nrow(snp_s) == nrow(map_hyb_005))
  snp_dt <- map_hyb_005[, .(marker, Chr, Pos)]
  snp_dt[, stat := snp_s[[snp_stat_col]][match(seq_len(.N), snp_s$MRK)]]
  add_gpos(snp_dt)
  snp_dt[chr_lens, on = "Chr", band := i.band]
  snp_dt[region_snps, on = "marker", s2_group := i.s2_group]
  snp_dt[, is_region := !is.na(s2_group)]

  y_top <- max(snp_dt$stat, na.rm = TRUE)
  arrow_y_head <- y_top * 1.35
  arrow_y_tail <- y_top * 1.55
  arrows_dt <- if (n_region > 0) snp_dt[is_region == TRUE & s2_group %in% region_counts$s2_group,
                                        .(gpos_mid = mean(range(gpos))), by = s2_group] else
    data.table(s2_group = character(0), gpos_mid = numeric(0))

  message(sprintf("[%s] %d raw crossings -> %d regions coloured; %d floor survivors -> %d regions arrowed",
                  tag, nrow(raw), length(all_regions), nrow(floor_units), n_region))
  if (n_region > 0)
    message("       ", paste(sprintf("%s=%d", region_counts$s2_group, region_counts$N), collapse = ", "))

  p <- ggplot() +
    geom_point(data = snp_dt[is_region == FALSE & band == FALSE], aes(gpos, stat), colour = "grey75", size = 0.3, alpha = 0.6) +
    geom_point(data = snp_dt[is_region == FALSE & band == TRUE], aes(gpos, stat), colour = "grey50", size = 0.3, alpha = 0.6) +
    geom_point(data = snp_dt[is_region == TRUE], aes(gpos, stat, colour = s2_group), size = 0.9) +
    { if (n_region > 0) geom_label_repel(
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
    labs(x = "Chromosome", y = title_stat, title = tag) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 11, face = "bold", hjust = 0))
  p
}

p_pc1 <- process_one(file.path(UNIT_DIR, "PC1_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15,
                     file.path(SNP_DIR, "PC1_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)",
                     "PC1", "BF (dB)",
                     file.path(DATA_DIR, "moduleB_stage1_S1units_null.rds"), "floor1")
p_pc2 <- process_one(file.path(UNIT_DIR, "PC2_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15,
                     file.path(SNP_DIR, "PC2_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)",
                     "PC2", "BF (dB)",
                     file.path(DATA_DIR, "moduleB_stage1_S1units_null.rds"), "floor2")
p_c2  <- process_one(file.path(UNIT_DIR, "mito_C2_S1units_summary_contrast.out"), "log10(1/pval)", 3,
                     file.path(SNP_DIR, "mito_C2_fullSNP_stage1Omega_summary_contrast.out"), "log10(1/pval)",
                     "mitoC2", expression(C2~-log[10](p)),
                     file.path(DATA_DIR, "moduleB_stage1_mitoC2_null.rds"), "floor3")
p_bw  <- process_one(file.path(UNIT_DIR, "bio_winter_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15,
                     file.path(SNP_DIR, "bio_winter_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)",
                     "bio_winter", "BF (dB)",
                     file.path(DATA_DIR, "moduleB_stage1_bio_winter_null.rds"), "floor")

combined <- (p_pc1 / p_pc2 / p_c2 / p_bw)

outpng <- file.path(FIGDIR, "moduleB_stage1_snp_manhattan_combined.png")
outpdf <- file.path(FIGDIR, "moduleB_stage1_snp_manhattan_combined.pdf")
ggsave(outpng, combined, width = 12, height = 15, dpi = 300, limitsize = FALSE)
ggsave(outpdf, combined, width = 12, height = 15, limitsize = FALSE)
cat("wrote", outpng, "and", outpdf, "\n")
