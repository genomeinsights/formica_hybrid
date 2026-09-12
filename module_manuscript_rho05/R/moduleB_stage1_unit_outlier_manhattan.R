## =========================================================
## module_manuscript_rho05 -- Stage-1-unit outlier Manhattan (BF>=15 raw
## crossings only), coloured by Stage-2 (rho05) cluster, floor survivors
## called out with an arrow
## =========================================================
## One point per Stage-1 unit crossing the raw threshold (BF(dB)>=15 for
## PC1/PC2/bio_winter; -log10(p)>=3 for mitoC2), at its own tested statistic
## and its core_snp's genomic position. Colour = which canonical Stage-2
## (rho05) cluster contains that unit's core_snp (LDscnR::default_cluster_
## colours(), recycled -- not unique per cluster). Floor survivors (beat all
## 10,000 Omega-structured null draws) are overplotted in BLACK at a LARGER
## point size, with an arrow (ggrepel) pointing to each one, labelled by
## group_id.
##
## Unlike moduleB_stage1_region_manhattan.R (every SNP in every TESTED
## cluster) or moduleB_stage1_snp_manhattan_ldmanhattan.R (every SNP in the
## genome), this plots only the raw-crossing Stage-1 UNITS themselves --
## one point each, not their member SNPs.
##
## Run from the repo root:
##   Rscript module_manuscript_rho05/R/moduleB_stage1_unit_outlier_manhattan.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(ggrepel) })
devtools::load_all("~/gitlab/LDscnR/")

BP_DIR   <- "module_manuscript_rho05/baypass_stage1"
UNIT_DIR <- file.path(BP_DIR, "aland_excluded_S1units")
DATA_DIR <- "module_manuscript_rho05/data"
FIGDIR   <- "module_manuscript_rho05/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

## ---- shared inputs ---------------------------------------------------------
load("data/hybrids_only_maf005.Rdata")   # map_hyb_005 (Chr, Pos, marker)
stage1 <- readRDS("module0_ld_pruning/data/pruned_stage1.rds")
cl <- as.data.table(stage1$clusters)
cl5 <- cl[n_snps >= 5]                                  # the 18,361 TESTED clusters
group_order <- readLines(file.path(UNIT_DIR, "S1units_group_order.txt"))
stopifnot(identical(group_order, paste0("S1_", cl5$CL_id)))
cl5[, group_id := group_order]
cl5 <- map_hyb_005[, .(marker, Chr, Pos)][cl5, on = c(marker = "core_snp")]
setnames(cl5, "marker", "core_snp")
stopifnot("core_snp missing Chr/Pos after join" = all(!is.na(cl5$Chr) & !is.na(cl5$Pos)))

## ---- canonical Stage-2 (rho05) clustering: core_snp -> s2_group -----------
message("Loading canonical Stage-2 (rho05) clustering...")
s2 <- readRDS("module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05.rds")
g2 <- as.data.table(s2$groups)
marker2s2 <- g2[, .(marker = unlist(members)), by = .(s2_group = group_id)]
setkey(marker2s2, marker)
cl5[marker2s2, on = .(core_snp = marker), s2_group := i.s2_group]
stopifnot("some core_snp not found in Stage-2 (rho05) clustering" = all(!is.na(cl5$s2_group)))

## ---- genome-wide x-axis (chromosome offsets) -------------------------------
chr_lens <- map_hyb_005[, .(len = max(Pos)), by = Chr]
chr_lens[, chr_num := as.integer(sub("Chr", "", Chr))]
setorder(chr_lens, chr_num)
chr_lens[, offset := cumsum(shift(len, fill = 0)) + (seq_len(.N) - 1) * 3e6]
chr_lens[, mid := offset + len / 2]
add_gpos <- function(dt) dt[chr_lens, on = "Chr", `:=`(gpos = Pos + i.offset)]
add_gpos(cl5)

PAL <- default_cluster_colours()

## ---- one covariate: threshold, colour by Stage-2 cluster, flag floor survivors --
process_one <- function(stat_file, stat_col, thresh, thresh_label, tag, title_stat,
                        null_file, null_flag_col) {
  message("[", tag, "] loading ", stat_file)
  s <- fread(stat_file)
  stopifnot(nrow(s) == nrow(cl5))
  cl5s <- copy(cl5)
  cl5s[, stat := s[[stat_col]][match(seq_len(.N), s$MRK)]]

  raw <- cl5s[stat >= thresh]
  message("[", tag, "] ", nrow(raw), " RAW threshold crossings (", thresh_label, ")")

  nullobj <- readRDS(null_file)
  stopifnot(null_flag_col %in% names(nullobj))
  floor_ids <- nullobj$group_id[as.logical(nullobj[[null_flag_col]])]
  raw[, is_floor := group_id %in% floor_ids]
  n_floor <- sum(raw$is_floor)
  message("[", tag, "] ", n_floor, " of those survive the null-calibrated floor test")

  ## recycle the default cluster-colour palette across however many distinct
  ## Stage-2 groups appear among the raw crossings (not unique per group)
  s2_levels <- sort(unique(raw$s2_group))
  cols <- setNames(rep(PAL, length.out = length(s2_levels)), s2_levels)

  p <- ggplot(raw, aes(gpos, stat)) +
    geom_point(data = raw[is_floor == FALSE], aes(colour = s2_group), size = 2.2) +
    geom_point(data = raw[is_floor == TRUE], colour = "black", size = 4.2) +
    { if (n_floor > 0) geom_label_repel(
        data = raw[is_floor == TRUE], aes(label = group_id),
        colour = "black", fill = "white", size = 3.2, fontface = "bold",
        arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
        segment.size = 0.6, segment.colour = "black",
        box.padding = 0.6, point.padding = 0.3, min.segment.length = 0,
        max.overlaps = Inf, seed = 1) } +
    geom_hline(yintercept = thresh, linetype = 2, colour = "grey40", linewidth = 0.3) +
    scale_colour_manual(values = cols, guide = if (length(s2_levels) > 20) "none" else "legend",
                        name = "Stage-2 (rho05) cluster") +
    scale_x_continuous(breaks = chr_lens$mid, labels = chr_lens$chr_num, expand = c(0.01, 0)) +
    labs(x = "Chromosome", y = title_stat,
         title = sprintf("Stage-1-direct %s: raw threshold crossings (%s), coloured by Stage-2 (rho05) cluster", tag, thresh_label),
         subtitle = sprintf("%d Stage-1 units tested; %d cross the raw threshold; %d survive the null-calibrated floor test (black, arrowed)",
                            nrow(cl5s), nrow(raw), n_floor)) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
  outpng <- file.path(FIGDIR, sprintf("moduleB_stage1_%s_unit_outlier_manhattan.png", tag))
  ggsave(outpng, p, width = 14, height = 5.5, dpi = 200)
  message("[", tag, "] wrote ", outpng)
  invisible(p)
}

process_one(file.path(UNIT_DIR, "PC1_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15, "BF(dB)>=15",
           "PC1", "BF(dB)",
           file.path(DATA_DIR, "moduleB_stage1_S1units_null.rds"), "floor1")
process_one(file.path(UNIT_DIR, "PC2_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15, "BF(dB)>=15",
           "PC2", "BF(dB)",
           file.path(DATA_DIR, "moduleB_stage1_S1units_null.rds"), "floor2")
process_one(file.path(UNIT_DIR, "mito_C2_S1units_summary_contrast.out"), "log10(1/pval)", 3, "-log10(p)>=3",
           "mitoC2", "C2 -log10(p)",
           file.path(DATA_DIR, "moduleB_stage1_mitoC2_null.rds"), "floor3")
process_one(file.path(UNIT_DIR, "bio_winter_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15, "BF(dB)>=15",
           "bio_winter", "BF(dB)",
           file.path(DATA_DIR, "moduleB_stage1_bio_winter_null.rds"), "floor")

message("\n[moduleB-stage1-unit-outlier-manhattan] done")
