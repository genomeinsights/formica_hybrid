## =========================================================
## module_manuscript_rho05 -- per-SNP Manhattan, coloured by SIGNIFICANT
## Stage-2 (rho05) region (floor-survivor-derived), with a region legend
## (Stage-1 floor-survivor count per region) and a downward arrow + label
## over each region
## =========================================================
## Every SNP in the genome is plotted at its own full-SNP BayPass value
## (like moduleB_stage1_snp_manhattan_ldmanhattan.R). Colouring differs: a
## Stage-2 (rho05) cluster is "significant" here iff it contains the
## core_snp of >=1 FLOOR-SURVIVOR Stage-1 unit (beats all 10,000
## Omega-structured null draws, not just a raw BF/C2 threshold crossing).
## Every member SNP of a significant Stage-2 cluster is coloured by that
## cluster's OWN colour (not a single "black" highlight) -- one distinct
## colour per significant region, small enough in number to legend
## explicitly. All other SNPs (untested, tested-but-not-crossing, or
## crossing-but-not-floor-surviving) are grey.
##
## Legend: one entry per significant Stage-2 region, labelled
## "<s2_group> (<n>)" where n = the number of Stage-1 floor-survivor units
## whose core_snp falls in that region (>1 only when two nearby Stage-1
## units both survive and were merged into the same Stage-2 cluster).
##
## Arrow: one per significant region, a vertical arrow from above pointing
## straight down onto the region's genomic midpoint, labelled with the
## Stage-2 group_id.
##
## Run from the repo root:
##   Rscript module_manuscript_rho05/R/moduleB_stage1_snp_manhattan_by_region.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2) })
devtools::load_all("~/gitlab/LDscnR/")

BP_DIR   <- "module_manuscript_rho05/baypass_stage1"
UNIT_DIR <- file.path(BP_DIR, "aland_excluded_S1units")
SNP_DIR  <- file.path(BP_DIR, "aland_excluded")
DATA_DIR <- "module_manuscript_rho05/data"
FIGDIR   <- "module_manuscript_rho05/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

## ---- shared inputs ---------------------------------------------------------
load("data/hybrids_only_maf005.Rdata")   # map_hyb_005 (Chr, Pos, marker)
stage1 <- readRDS("module0_ld_pruning/data/pruned_stage1.rds")
cl <- as.data.table(stage1$clusters)
cl5 <- cl[n_snps >= 5]
group_order <- readLines(file.path(UNIT_DIR, "S1units_group_order.txt"))
stopifnot(identical(group_order, paste0("S1_", cl5$CL_id)))
cl5[, group_id := group_order]

message("Loading canonical Stage-2 (rho05) clustering...")
s2 <- readRDS("module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05.rds")
g2 <- as.data.table(s2$groups)
marker2s2 <- g2[, .(marker = unlist(members)), by = .(s2_group = group_id)]
setkey(marker2s2, marker)

## genome-wide x-axis (chromosome offsets) -- AESTHETIC UPDATE (2026-09-13):
## much smaller inter-chromosome gap (3e6 -> 3e5) and alternating light/dark
## background bands (classic Manhattan-plot convention) instead of a uniform
## white background.
chr_lens <- map_hyb_005[, .(len = max(Pos)), by = Chr]
chr_lens[, chr_num := as.integer(sub("Chr", "", Chr))]
setorder(chr_lens, chr_num)
chr_lens[, offset := cumsum(shift(len, fill = 0)) + (seq_len(.N) - 1) * 3e5]
chr_lens[, mid := offset + len / 2]
## USER FIX (2026-09-13): alternate by POSITIONAL order (chr_lens' row order,
## after setorder above), not chr_num's own parity -- chr_num skips 23
## (22 -> 24), so `chr_num %% 2 == 0` put two "even" chromosomes in a row and
## broke the alternating pattern at that point.
chr_lens[, band := seq_len(.N) %% 2 == 0]
add_gpos <- function(dt) dt[chr_lens, on = "Chr", `:=`(gpos = Pos + i.offset)]

## AUDIT/USER CLARIFICATION (2026-09-13): "significant" for COLOURING purposes
## means every Stage-2 region touched by ANY raw BF>=15 (or C2 threshold)
## Stage-1 crossing -- not just floor survivors. There can be many such
## regions, so they are coloured with the recycled LDscnR default palette
## (not unique per region, no legend entry per region -- matching the
## convention used elsewhere in this pipeline). The floor-survivor subset is
## then called out ON TOP via arrows (region name) and a legend restricted to
## just those regions (name + count of contributing Stage-1 floor survivors),
## using scale_colour_manual's `breaks` to show only a subset of an
## otherwise-recycled discrete colour scale in the legend.
## USER FIX (2026-09-13): drop near-white/pale palette entries -- as labels
## (arrow text/border, coloured to match a floor-survivor region), these are
## hard to read against the plot's white background. Perceived luminance
## (0.299R+0.587G+0.114B) > 0.75 excludes the palest ~10 of 41 default
## colours (e.g. #B3E2CD, a pale mint that was flagged as unreadable);
## still 31 colours to recycle across regions.
PAL_ALL <- default_cluster_colours()
lum <- { rgb <- t(col2rgb(PAL_ALL)) / 255; 0.299 * rgb[, 1] + 0.587 * rgb[, 2] + 0.114 * rgb[, 3] }
PAL <- PAL_ALL[lum <= 0.75]

process_one <- function(unit_stat_file, unit_stat_col, thresh, thresh_label,
                        snp_stat_file, snp_stat_col, tag, title_stat,
                        null_file = NULL, null_flag_col = NULL) {
  message("\n[", tag, "] loading Stage-1-unit scan: ", unit_stat_file)
  s <- fread(unit_stat_file)
  stopifnot(nrow(s) == nrow(cl5))
  cl5s <- copy(cl5)
  cl5s[, stat := s[[unit_stat_col]][match(seq_len(.N), s$MRK)]]

  ## ALL raw threshold crossings (BF>=15 / C2 threshold) -> ALL touched Stage-2 regions
  raw <- cl5s[stat >= thresh]
  raw[marker2s2, on = .(core_snp = marker), s2_group := i.s2_group]
  stopifnot("a raw-crossing core_snp is missing from the Stage-2 (rho05) clustering" =
              nrow(raw) == 0 || all(!is.na(raw$s2_group)))
  all_regions <- sort(unique(raw$s2_group))
  message("[", tag, "] ", nrow(raw), " RAW threshold crossings (", thresh_label,
          ") -> ", length(all_regions), " Stage-2 (rho05) region(s) coloured")

  ## no per-variable floor test available for this covariate (e.g. bio6/
  ## bio11 individually -- only their combined bio_winter axis is calibrated,
  ## Issue 7 decision) -- colour all raw-crossing regions, but no arrows/
  ## floor-survivor subset, and say so explicitly in the caption.
  no_floor_test <- is.null(null_file)
  if (no_floor_test) {
    floor_units <- cl5s[0]
    region_counts <- data.table(s2_group = character(0), N = integer(0))
    n_region <- 0L
    message("[", tag, "] no per-variable floor test -- raw crossings coloured, no arrows")
  } else {
    nullobj <- readRDS(null_file)
    stopifnot(null_flag_col %in% names(nullobj))
    floor_ids <- nullobj$group_id[as.logical(nullobj[[null_flag_col]])]
    floor_units <- cl5s[group_id %in% floor_ids]
    message("[", tag, "] ", nrow(floor_units), " floor-survivor Stage-1 unit(s)")

    ## floor-survivor Stage-1 units -> their Stage-2 (rho05) region(s), with a
    ## per-region count of how many floor-survivor Stage-1 units map into it --
    ## this subset drives the arrows + restricted legend, NOT the colouring
    floor_units[marker2s2, on = .(core_snp = marker), s2_group := i.s2_group]
    stopifnot("a floor-survivor core_snp is missing from the Stage-2 (rho05) clustering" =
                nrow(floor_units) == 0 || all(!is.na(floor_units$s2_group)))
    region_counts <- floor_units[, .N, by = s2_group]
    setorder(region_counts, -N)
    n_region <- nrow(region_counts)
    message("[", tag, "] -> ", n_region, " of those are floor-survivor region(s): ",
            paste(sprintf("%s(n=%d)", region_counts$s2_group, region_counts$N), collapse = ", "))
  }

  ## recycled colour for EVERY raw-crossing region; legend restricted to the
  ## floor-survivor subset via `breaks` (all regions still coloured in-plot)
  region_cols <- setNames(rep(PAL, length.out = length(all_regions)), all_regions)

  ## every member SNP of ANY raw-crossing Stage-2 region, for colouring
  region_snps <- if (length(all_regions) > 0) g2[group_id %in% all_regions,
                                      .(marker = unlist(members)), by = .(s2_group = group_id)] else
                                      data.table(marker = character(0), s2_group = character(0))

  ## ---- REAL per-SNP y-value, from the full-genome BayPass scan -----------
  message("[", tag, "] loading full-SNP scan: ", snp_stat_file)
  snp_s <- fread(snp_stat_file)
  stopifnot(nrow(snp_s) == nrow(map_hyb_005))
  snp_dt <- map_hyb_005[, .(marker, Chr, Pos)]
  snp_dt[, stat := snp_s[[snp_stat_col]][match(seq_len(.N), snp_s$MRK)]]
  add_gpos(snp_dt)
  snp_dt[chr_lens, on = "Chr", band := i.band]
  snp_dt[region_snps, on = "marker", s2_group := i.s2_group]
  snp_dt[, is_region := !is.na(s2_group)]

  ## arrow/label position: genomic midpoint of each region's member SNPs,
  ## pointing straight down from a shared height above the tallest point
  y_top <- max(snp_dt$stat, na.rm = TRUE)
  arrow_y_head <- y_top * 1.35
  arrow_y_tail <- y_top * 1.55
  arrows_dt <- if (n_region > 0) snp_dt[is_region == TRUE & s2_group %in% region_counts$s2_group,
                                        .(gpos_mid = mean(range(gpos))), by = s2_group] else
    data.table(s2_group = character(0), gpos_mid = numeric(0))

  message("[", tag, "] ", nrow(snp_dt), " SNPs plotted (", sum(snp_dt$is_region),
          " in ", n_region, " significant region(s))")

  ## AESTHETIC UPDATE (2026-09-13): alternating chromosomes distinguished via
  ## TWO GREY SHADES on the background (non-significant) points themselves,
  ## not a rectangle band behind them. Arrow LABEL text/border coloured to
  ## match its region (same manual scale as the points, via the shared
  ## `s2_group` aesthetic) instead of black -- the colour match is what
  ## identifies a floor survivor now, so the legend is dropped entirely
  ## (guide = "none"); the per-region Stage-1-cluster count moves into the
  ## plot caption instead.
  caption_txt <- if (no_floor_test)
    "No per-variable null-calibrated floor test exists for this covariate (only its combined bio_winter axis is calibrated) -- coloured regions are RAW threshold crossings only, not validated candidates."
  else if (n_region > 0)
    sprintf("Floor-survivor regions (arrowed, label coloured to match): %s.",
            paste(sprintf("%s = %d Stage-1 cluster%s", region_counts$s2_group, region_counts$N,
                          ifelse(region_counts$N == 1, "", "s")), collapse = "; "))
  else "No Stage-2 region contains a floor-survivor Stage-1 unit for this covariate."

  p <- ggplot() +
    geom_point(data = snp_dt[is_region == FALSE & band == FALSE], aes(gpos, stat), colour = "grey75", size = 0.4, alpha = 0.6) +
    geom_point(data = snp_dt[is_region == FALSE & band == TRUE], aes(gpos, stat), colour = "grey50", size = 0.4, alpha = 0.6) +
    geom_point(data = snp_dt[is_region == TRUE], aes(gpos, stat, colour = s2_group), size = 1.1) +
    { if (n_region > 0) geom_segment(
        data = arrows_dt, aes(x = gpos_mid, xend = gpos_mid, y = arrow_y_tail, yend = arrow_y_head),
        arrow = arrow(length = unit(0.18, "cm"), type = "closed"), linewidth = 0.6, colour = "black") } +
    { if (n_region > 0) geom_label(
        data = arrows_dt, aes(x = gpos_mid, y = arrow_y_tail, label = s2_group, colour = s2_group),
        vjust = 0, size = 3.2, fontface = "bold", fill = "white", label.padding = unit(0.15, "lines"),
        show.legend = FALSE) } +
    geom_hline(yintercept = thresh, linetype = 2, colour = "black", linewidth = 0.3) +
    scale_colour_manual(values = region_cols, guide = "none") +
    scale_x_continuous(breaks = chr_lens$mid, labels = chr_lens$chr_num, expand = c(0.005, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.45))) +
    labs(x = "Chromosome", y = title_stat,
         title = sprintf("Stage-1-direct %s: every SNP, coloured by Stage-2 (rho05) region containing a raw BF/C2 threshold crossing", tag),
         subtitle = sprintf("%d Stage-1 units tested; %d raw crossing(s) -> %d region(s) coloured; %d floor-survivor unit(s) -> %d region(s) arrowed; %d SNPs coloured",
                            nrow(cl5s), nrow(raw), length(all_regions), nrow(floor_units), n_region, sum(snp_dt$is_region)),
         caption = caption_txt) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                                      plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"))
  outpng <- file.path(FIGDIR, sprintf("moduleB_stage1_%s_snp_manhattan_by_region.png", tag))
  ggsave(outpng, p, width = 16, height = 5.5, dpi = 200, limitsize = FALSE)
  message("[", tag, "] wrote ", outpng)
  invisible(p)
}

process_one(file.path(UNIT_DIR, "PC1_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15, "BF(dB)>=15",
           file.path(SNP_DIR, "PC1_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)",
           "PC1", "BF(dB)",
           file.path(DATA_DIR, "moduleB_stage1_S1units_null.rds"), "floor1")
process_one(file.path(UNIT_DIR, "PC2_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15, "BF(dB)>=15",
           file.path(SNP_DIR, "PC2_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)",
           "PC2", "BF(dB)",
           file.path(DATA_DIR, "moduleB_stage1_S1units_null.rds"), "floor2")
process_one(file.path(UNIT_DIR, "mito_C2_S1units_summary_contrast.out"), "log10(1/pval)", 3, "-log10(p)>=3",
           file.path(SNP_DIR, "mito_C2_fullSNP_stage1Omega_summary_contrast.out"), "log10(1/pval)",
           "mitoC2", "C2 -log10(p)",
           file.path(DATA_DIR, "moduleB_stage1_mitoC2_null.rds"), "floor3")

process_one(file.path(UNIT_DIR, "bio6_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15, "BF(dB)>=15",
           file.path(SNP_DIR, "bio6_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)",
           "bio6", "BF(dB)")   # no per-variable floor test -- see bio_winter

## bio11 full-SNP scan still running on mini2 (Issue 1 audit fix rerun) --
## add once available (no per-variable floor test either -- see bio_winter).

message("\n[moduleB-stage1-snp-manhattan-by-region] done")
