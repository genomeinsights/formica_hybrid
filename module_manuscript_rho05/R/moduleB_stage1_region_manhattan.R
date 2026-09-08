## =========================================================
## module_manuscript_rho05 -- SNP-level Manhattan, coloured by Stage-2-
## assembled outlier region (Stage-1-direct BayPass scan, all 3 covariates)
## =========================================================
## The Stage-1 scan (PC1, PC2, mito C2) tests one best-SNP per Stage-1
## cluster (n_snps>=5, 18,361 units) -- deliberately NOT merged into Stage-2,
## per the module_3sp Goldilocks finding (over-merging dilutes outlier
## power). This script does the POST-HOC region assembly module_3sp
## describes for exactly this situation (ld_outlier_test(assembly=
## "stage2_discovered")'s pattern, reimplemented directly here since our
## statistic is BayPass BF/C2, not EMMAX p-values): take only the
## SIGNIFICANT Stage-1 clusters and Stage-2-merge THOSE into regions, purely
## to describe their physical extent -- this cannot create, remove, or
## re-rank a significant unit, it only groups already-significant ones that
## are close and correlated enough to plausibly be one signal.
##
## Every member SNP of a TESTED Stage-1 cluster (n_snps>=5) is then plotted,
## inheriting its own cluster's tested statistic (BF(dB) or C2) -- SNPs in
## untested clusters (n_snps<5) have no test result and are excluded, not
## shown as flat/zero. Colour = which assembled Stage-2 region a SNP's
## cluster belongs to (grey = tested, not in any significant region).
##
## Descriptive thresholds (NOT calibrated against a structured null -- that
## step doesn't exist yet for this Stage-1 scan): BF(dB) >= 15 for PC1/PC2
## (matches the manuscript's own established "descriptive outlier
## threshold"); p < 0.001 (-log10(pval) >= 3) for mito C2, the closest
## analogue.
##
## Region assembly: min_r2_rho = 0.5 (today's Stage-2 default for
## LD-reduced-units use), score_threshold = 0.80, genetic_map + cM_threshold
## = 0.5cM (matches the LD-reduced-units convention used elsewhere).
##
## Run from the repo root:  Rscript module_manuscript_rho05/R/moduleB_stage1_region_manhattan.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2) })
devtools::load_all("~/gitlab/LDscnR/")

BP_DIR <- "module_manuscript_rho05/baypass_stage1"
UNIT_DIR <- file.path(BP_DIR, "aland_excluded_S1units")
FIGDIR <- "module_manuscript_rho05/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

## ---- shared inputs ---------------------------------------------------------
load("data/hybrids_only_maf005.Rdata")   # map_hyb_005 (Chr, Pos, marker), ld_decay
stage1 <- readRDS("module0_ld_pruning/data/pruned_stage1.rds")
cl <- as.data.table(stage1$clusters)
cl5 <- cl[n_snps >= 5]                                  # the 18,361 TESTED clusters
group_order <- readLines(file.path(UNIT_DIR, "S1units_group_order.txt"))
stopifnot(identical(group_order, paste0("S1_", cl5$CL_id)))   # confirms MRK row order == cl5 row order
cl5[, group_id := group_order]

rec_map <- fread("./data/Frufa_DTOL_PR.ref_genome.recmap")
rec_map[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
genetic_map <- rec_map[, .(Chr, Pos = pos, cM)]

## chromosome plotting order/offsets (by number)
chr_lens <- map_hyb_005[, .(len = max(Pos)), by = Chr]
chr_lens[, chr_num := as.integer(sub("Chr", "", Chr))]
setorder(chr_lens, chr_num)
chr_lens[, offset := cumsum(shift(len, fill = 0)) + (seq_len(.N) - 1) * 3e6]
chr_lens[, mid := offset + len / 2]

add_gpos <- function(dt) {
  dt[chr_lens, on = "Chr", `:=`(gpos = Pos + i.offset)]
}

## ---- one covariate: threshold, assemble Stage-2 regions, plot -------------
process_one <- function(stat_file, stat_col, thresh, thresh_label, tag, title_stat) {
  message("[", tag, "] loading ", stat_file)
  s <- fread(stat_file)
  stopifnot(nrow(s) == nrow(cl5))
  cl5s <- copy(cl5)
  cl5s[, stat := s[[stat_col]][match(seq_len(.N), s$MRK)]]   # MRK is a 1..n row index matching cl5's order

  sig_ids <- cl5s[stat >= thresh, CL_id]
  message("[", tag, "] ", length(sig_ids), " significant Stage-1 units (", thresh_label, ")")

  ## ---- assemble significant clusters into Stage-2 regions ------------------
  regions <- NULL
  if (length(sig_ids) >= 1) {
    sub_stage1 <- list(map_snp = stage1$map_snp[CL_id %in% sig_ids],
                       clusters = stage1$clusters[CL_id %in% sig_ids])
    if (length(sig_ids) == 1) {
      ## nothing to merge -- one cluster is its own region
      grp <- data.table(group_id = "R1", Chr = sub_stage1$clusters$Chr,
                        members = sub_stage1$clusters$members)
    } else {
      res <- ld_prune_and_eMLG(
        GTs = GTs_hybrids_005, stage1 = sub_stage1, ld_w_col = "ld_w_095",
        ld_w_threshold = 0, min_n_loci_flag = 1, score_threshold = 0.80,
        min_r2 = NULL, min_r2_rho = 0.5, LD_decay = ld_decay,
        genetic_map = genetic_map, cM_threshold = 0.5,
        compute_unflagged_eMLG = FALSE
      )
      grp <- as.data.table(res$groups)[, .(group_id, Chr, members)]
    }
    grp[, region_id := paste0(tag, "_R", seq_len(.N))]
    grp[, n_snp_region := lengths(members)]
    setorder(grp, -n_snp_region)
    regions <- grp[, .(region_id, Chr, n_snp_region,
                       from = vapply(members, function(mk) min(map_hyb_005[marker %in% mk, Pos]), numeric(1)),
                       to   = vapply(members, function(mk) max(map_hyb_005[marker %in% mk, Pos]), numeric(1)))]
    setorder(regions, Chr, from)
    print(regions)
    saveRDS(list(regions = regions, groups = grp),
            file.path("module_manuscript_rho05/data", sprintf("moduleB_stage1_%s_regions.rds", tag)))

    ## SNP -> region_id map (every member SNP of every assembled region)
    snp2region <- grp[, .(marker = unlist(members)), by = region_id]
  } else {
    snp2region <- data.table(marker = character(0), region_id = character(0))
  }

  ## ---- every SNP in a TESTED cluster, inheriting its cluster's stat --------
  snp_dt <- cl5s[, .(marker = unlist(members)), by = .(CL_id, stat)]
  snp_dt <- map_hyb_005[, .(marker, Chr, Pos)][snp_dt, on = "marker"]
  snp_dt <- add_gpos(snp_dt)
  snp_dt[snp2region, on = "marker", region_id := i.region_id]
  snp_dt[, is_sig := !is.na(region_id)]

  n_region <- if (!is.null(regions)) nrow(regions) else 0L
  message("[", tag, "] ", nrow(snp_dt), " member SNPs plotted (", sum(snp_dt$is_sig),
          " in ", n_region, " assembled region(s))")

  p <- ggplot() +
    geom_point(data = snp_dt[is_sig == FALSE], aes(gpos, stat), colour = "grey75", size = 0.4, alpha = 0.6) +
    geom_point(data = snp_dt[is_sig == TRUE], aes(gpos, stat, colour = region_id), size = 1.1) +
    geom_hline(yintercept = thresh, linetype = 2, colour = "black", linewidth = 0.3) +
    scale_x_continuous(breaks = chr_lens$mid, labels = chr_lens$chr_num, expand = c(0.01, 0)) +
    scale_colour_viridis_d(guide = if (n_region > 15) "none" else "legend") +
    labs(x = "Chromosome", y = title_stat,
         title = sprintf("Stage-1-direct %s: every SNP in a tested cluster (n_snps>=5), coloured by Stage-2-assembled region", tag),
         subtitle = sprintf("%d Stage-1 units tested; %d significant (%s) -> %d assembled region(s); descriptive threshold, not null-calibrated",
                            nrow(cl5s), length(sig_ids), thresh_label, n_region)) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
  outpng <- file.path(FIGDIR, sprintf("moduleB_stage1_%s_region_manhattan.png", tag))
  ggsave(outpng, p, width = 14, height = 5, dpi = 200)
  message("[", tag, "] wrote ", outpng)
  invisible(list(snp_dt = snp_dt, regions = regions))
}

r_pc1 <- process_one(file.path(UNIT_DIR, "PC1_S1units_withOmega_summary_betai_reg.out"),
                     "BF(dB)", 15, "BF(dB)>=15", "PC1", "BF(dB)")
r_pc2 <- process_one(file.path(UNIT_DIR, "PC2_S1units_withOmega_summary_betai_reg.out"),
                     "BF(dB)", 15, "BF(dB)>=15", "PC2", "BF(dB)")
r_c2  <- process_one(file.path(UNIT_DIR, "mito_C2_S1units_summary_contrast.out"),
                     "log10(1/pval)", 3, "-log10(p)>=3", "mitoC2", "C2 -log10(p)")

message("[moduleB-stage1-region-manhattan] done")
