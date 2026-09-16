## =========================================================
## module_localscore_crosscheck -- FALSE-POSITIVE DEMONSTRATION: a high-scoring
## NULL draw for each of the four covariates, styled EXACTLY like the real
## Stage-2-region-coloured, arrow-labelled Manhattan figure
## (moduleB_stage1_snp_manhattan_combined.R)
## =========================================================
## THESE ARE NOT REAL RESULTS. Each panel plots one Omega-structured NULL
## draw (no causal link to climate or mitotype) chosen because it happens
## to produce an unusually high number of "significant" local-score
## windows (local_score_null_check_S1units.R already showed
## ~19% of continuous-null draws and the great majority of mitoC2-null
## draws produce >=1 such window purely by chance). The purpose is
## explicitly to illustrate how easily the local-score method's analytic
## threshold, run at Stage-1-cluster resolution and coloured/arrowed with
## the same Stage-2 (rho05) region convention as the real covariate
## figures, manufactures what LOOKS like a compelling multi-region
## candidate signal out of pure population structure.
##
## PC1/PC2/bio_winter panels each use a DIFFERENT high-scoring draw from
## the shared continuous-covariate null pool (draws #536, #6, #204 -- 4, 3,
## 3 significant windows respectively, out of the top 5 scanned) so the
## three panels are not literal repeats of each other. mitoC2 uses its own
## dedicated contrast-mode null pool's worst draw (#480, 16 windows).
##
## "Raw crossing" / "significant" here mean exactly what they mean for the
## real data (raw BF/C2 threshold crossing; local-score window membership
## in place of the floor-survivor null test) -- but applied to a null
## draw, so a "significant" region below is, by construction, a false
## positive, not a candidate locus.
##
## Reads : module_localscore_crosscheck/data/localscore_nullcheck_S1units.rds
##         module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/null/bf_matrices/
##           cRegen_bf_b01..05.rds, mitoC2_bf_b01..05.rds
##         module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05.rds (canonical Stage-2)
## Writes: module_localscore_crosscheck/Figures/local_score_null_examples.{png,pdf}
##
## Run from the repo root:
##   Rscript module_localscore_crosscheck/R/local_score_manhattan_null_examples.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })
devtools::load_all("~/gitlab/LDscnR/")
source("~/gitlab/baypass_public-master/utils/baypass_utils.R")

UNIT_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded_S1units"
BFDIR    <- file.path(UNIT_DIR, "null", "bf_matrices")
DATA     <- "module_localscore_crosscheck/data"
FIGDIR   <- "module_localscore_crosscheck/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

## ---- shared inputs (identical to moduleB_stage1_snp_manhattan_combined.R) --
load("data/hybrids_only_maf005.Rdata")
stage1 <- readRDS("module0_ld_pruning/data/pruned_stage1.rds")
cl <- as.data.table(stage1$clusters)
cl5 <- cl[n_snps >= 5]
group_order <- readLines(file.path(UNIT_DIR, "S1units_group_order.txt"))
stopifnot(identical(group_order, paste0("S1_", cl5$CL_id)))
cl5[, group_id := group_order]
cl5[, core_pos := as.integer(sub(".*:", "", core_snp))]

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

base_pos <- data.table(row = seq_len(nrow(cl5)), Chr = cl5$Chr, Pos = cl5$core_pos)
setorder(base_pos, Chr, Pos)
ord <- base_pos$row
pos_df <- as.data.frame(base_pos[, .(Chr, Pos)])
MIN_NSNP <- 100L

pi_pc1 <- fread(file.path(UNIT_DIR, "PC1_S1units_withOmega_summary_pi_xtx.out"), select = c("MRK", "M_P"))$M_P
pi_c2  <- fread(file.path(UNIT_DIR, "mito_C2_S1units_summary_pi_xtx.out"), select = c("MRK", "M_P"))$M_P

## ---- pull one null draw's Stage-1-cluster stat vector (cl5/MRK order) -----
get_null_draw <- function(batch_pattern, batch, col) {
  M <- readRDS(sprintf(batch_pattern, batch))
  M[, col]   # length nrow(cl5), in ORIGINAL cl5/MRK order
}

## ---- build one panel: Stage-2-region-coloured, arrow-labelled, styled
##      EXACTLY like moduleB_stage1_snp_manhattan_combined.R -----------------
make_null_panel <- function(tag, draw_label, draw_idx, stat_vec_mrk, thresh, y_lab, pi_vec, is_bf) {
  cl5s <- copy(cl5)
  cl5s[, stat := stat_vec_mrk]   # already in cl5/MRK order

  ## compute.local.scores() draws a RANDOM p-value (-log10(runif())) for
  ## every NEGATIVE BF value (BayPass's own code, not this script's) --
  ## re-running the SAME BF-based null draw can therefore give a DIFFERENT
  ## window count each call. Seeded here (per draw index) purely so this
  ## demonstration figure is reproducible on rerun; the underlying
  ## instability is itself a real, separately-noted caveat of the method.
  set.seed(10000L + draw_idx)

  ## local-score windows for THIS null draw (Stage-1-cluster resolution)
  out <- if (is_bf) {
    invisible(capture.output(res <- compute.local.scores(
      snp.position = pos_df, snp.pi = pi_vec[ord], snp.bf = stat_vec_mrk[ord],
      xi = 1, pval.local.score.thres = 0.01, min.maf = 0.2, min.nsnp = MIN_NSNP)))
  } else {
    invisible(capture.output(res <- compute.local.scores(
      snp.position = pos_df, snp.pi = pi_vec[ord], snp.pvalue = stat_vec_mrk[ord],
      xi = 1, pval.local.score.thres = 0.01, min.maf = 0.2, min.nsnp = MIN_NSNP)))
  }
  win <- res$significant.windows
  n_win <- if (is.null(win) || nrow(win) == 0) 0L else nrow(win)

  ## ALL raw threshold crossings -> ALL touched Stage-2 regions (identical
  ## logic to the real-data script)
  raw <- cl5s[stat >= thresh]
  raw[marker2s2, on = .(core_snp = marker), s2_group := i.s2_group]
  raw_regions <- sort(unique(raw$s2_group))

  ## ONE ARROW PER WINDOW (not per touched Stage-2 region -- see
  ## local_score_snp_manhattan_manuscript_style.R's header for why: a
  ## window's peak unit is the strongest RAW value inside it, and is often
  ## not itself a raw crossing, so its Stage-2 group is looked up directly
  ## from the canonical clustering, independent of raw_regions, and unioned
  ## in below so the arrowed region still gets its own colour either way.
  arrows_dt <- data.table(win_id = character(0), s2_group = character(0), gpos_mid = numeric(0))
  if (n_win > 0) {
    win <- as.data.table(win)
    peak_col <- if (is_bf) "BF (dB) peak pos" else "-log10(p-val) peak pos"
    win[, peak_pos := as.integer(get(peak_col))]
    win[, mid_pos := as.integer(round((beg + end) / 2))]
    win[, win_id := paste0(tag, "_W", seq_len(.N))]
    peak_snp <- cl5s[.(Chr = win$chr, core_pos = win$peak_pos), on = .(Chr, core_pos), core_snp]
    win[, peak_s2 := marker2s2[.(peak_snp), on = "marker", s2_group]]
    win[, gpos_mid := mid_pos + chr_lens$offset[match(chr, chr_lens$Chr)]]
    arrows_dt <- win[!is.na(peak_s2), .(win_id, s2_group = peak_s2, gpos_mid)]
  }

  all_regions <- sort(union(raw_regions, unique(arrows_dt$s2_group)))
  region_cols <- setNames(rep(PAL, length.out = length(all_regions)), all_regions)
  region_snps <- if (length(all_regions) > 0) g2[group_id %in% all_regions,
                                      .(marker = unlist(members)), by = .(s2_group = group_id)] else
                                      data.table(marker = character(0), s2_group = character(0))

  ## expand to every genome-wide member SNP of a tested Stage-1 cluster,
  ## inheriting that cluster's null stat (identical to moduleB_stage1_
  ## region_manhattan.R's approach)
  snp_dt <- cl5s[, .(marker = unlist(members)), by = .(CL_id, stat)]
  snp_dt <- map_hyb_005[, .(marker, Chr, Pos)][snp_dt, on = "marker"]
  add_gpos(snp_dt)
  snp_dt[chr_lens, on = "Chr", band := i.band]
  snp_dt[region_snps, on = "marker", s2_group := i.s2_group]
  snp_dt[, is_region := !is.na(s2_group)]

  y_top <- max(snp_dt$stat, na.rm = TRUE)
  arrow_y_head <- y_top * 1.35
  arrow_y_tail <- y_top * 1.55

  message(sprintf("[%s null draw %s] %d raw crossings -> %d regions coloured; %d local-score window(s) -> %d window(s) arrowed",
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
         title = sprintf("%s -- NULL DRAW %s (NOT REAL DATA)", tag, draw_label),
         subtitle = sprintf("%d raw crossings -> %d Stage-2 regions coloured; %d local-score window(s) -> %d window(s) arrowed as \"significant\"",
                            nrow(raw), length(all_regions), n_win, nrow(arrows_dt))) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 11, face = "bold", hjust = 0, colour = "#B03A2E"),
          plot.subtitle = element_text(size = 7.5, colour = "grey30"))
}

## ---- picks: 3 different high-scoring continuous-null draws (PC1/PC2/
## bio_winter labels; from local_score_null_check_S1units.R's
## top-5 by window count), 1 high-scoring mitoC2-null draw (its top-1) -------
batch_col <- function(k, batch_size = 200L) c(batch = ((k - 1L) %/% batch_size) + 1L,
                                              col   = ((k - 1L) %% batch_size) + 1L)

## NOTE: compute.local.scores() draws a random p-value for every negative
## BF value, so the raw (unseeded) window counts recorded by
## local_score_null_check_S1units.R for the continuous-null
## pool (draws #536/#6/#204: 4/3/3 windows) are not exactly reproducible
## call-to-call. These three picks were reselected from the same top
## candidates AFTER fixing the per-draw seed used below (see make_null_panel),
## keeping only draws that still show >=2 windows under that seed.
PICKS <- list(
  PC1        = list(pool = "cRegen_bf_b%02d.rds", draw = 624, pi = pi_pc1, is_bf = TRUE,  thresh = 15, y_lab = "BF (dB)"),
  PC2        = list(pool = "cRegen_bf_b%02d.rds", draw = 536, pi = pi_pc1, is_bf = TRUE,  thresh = 15, y_lab = "BF (dB)"),
  bio_winter = list(pool = "cRegen_bf_b%02d.rds", draw = 107, pi = pi_pc1, is_bf = TRUE,  thresh = 15, y_lab = "BF (dB)"),
  mitoC2     = list(pool = "mitoC2_bf_b%02d.rds", draw = 480, pi = pi_c2,  is_bf = FALSE, thresh = 3,  y_lab = expression(C2~-log[10](p)))
)

panels <- list()
for (tag in names(PICKS)) {
  p <- PICKS[[tag]]
  bc <- batch_col(p$draw)
  stat_vec <- get_null_draw(file.path(BFDIR, p$pool), bc["batch"], bc["col"])
  panels[[tag]] <- make_null_panel(tag, sprintf("#%d", p$draw), p$draw, stat_vec, p$thresh, p$y_lab, p$pi, p$is_bf)
}

combined <- (panels$PC1 / panels$PC2 / panels$bio_winter / panels$mitoC2)

outpng <- file.path(FIGDIR, "local_score_null_examples.png")
outpdf <- file.path(FIGDIR, "local_score_null_examples.pdf")
ggsave(outpng, combined, width = 12, height = 16, dpi = 300, limitsize = FALSE)
ggsave(outpdf, combined, width = 12, height = 16, limitsize = FALSE)
cat("wrote", outpng, "and", outpdf, "\n")
