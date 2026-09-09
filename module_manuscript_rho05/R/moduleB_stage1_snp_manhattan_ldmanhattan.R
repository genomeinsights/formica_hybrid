## =========================================================
## module_manuscript_rho05 -- SNP-level Manhattan via LDscnR's ld_manhattan(),
## coloured by the FULL canonical Stage-2 cluster (rho05), default palette
## =========================================================
## y-axis is now GENUINE per-SNP data: the full-genome BayPass scan
## (run_baypass_stage1_fullsnp.sh, all 1,114,423 SNPs, Stage-1-derived Omega),
## not a value inherited from a SNP's Stage-1 cluster. The earlier version
## gave every SNP in a Stage-1 cluster the SAME y (that cluster's one tested
## value), so the plot was still Stage-1-unit-resolution under the hood --
## flat horizontal steps of identical-height points, no real per-SNP
## variation, even though colouring already correctly used the full Stage-2
## grouping.
##
## Significance decision is UNCHANGED (Stage-1-direct, per the module_3sp
## Goldilocks finding: testing individual SNPs or full-Stage-2-merged units
## dilutes/inflates the test count vs the right-sized Stage-1 units) -- only
## the plotted y-value moved to real per-SNP data. Concretely:
##   1. Threshold the STAGE-1-UNIT scan (PC1/PC2 BF(dB)>=15, mito C2
##      -log10(p)>=3, descriptive, not null-calibrated) to get significant
##      Stage-1 clusters.
##   2. Find which FULL canonical Stage-2 (rho05) group each significant
##      Stage-1 cluster's core_snp falls in -- flags that group.
##   3. Every member SNP of a flagged Stage-2 group is coloured by that
##      group (LDscnR::default_cluster_colours(), recycled -- not unique per
##      region); its y-value is its OWN row in the full-SNP scan.
##   4. Grey (non-flagged) points also use their own real per-SNP value now.
##
## Run from the repo root:  Rscript module_manuscript_rho05/R/moduleB_stage1_snp_manhattan_ldmanhattan.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2) })
devtools::load_all("~/gitlab/LDscnR/")

BP_DIR <- "module_manuscript_rho05/baypass_stage1"
UNIT_DIR <- file.path(BP_DIR, "aland_excluded_S1units")
SNP_DIR  <- file.path(BP_DIR, "aland_excluded")   # full-SNP scan output
FIGDIR <- "module_manuscript_rho05/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

## ---- shared inputs ---------------------------------------------------------
load("data/hybrids_only_maf005.Rdata")   # map_hyb_005 (Chr, Pos, marker) -- same row order as u_DIEM.geno
stage1 <- readRDS("module0_ld_pruning/data/pruned_stage1.rds")
cl <- as.data.table(stage1$clusters)
cl5 <- cl[n_snps >= 5]
group_order <- readLines(file.path(UNIT_DIR, "S1units_group_order.txt"))
stopifnot(identical(group_order, paste0("S1_", cl5$CL_id)))
cl5[, group_id := group_order]

## ---- FULL canonical Stage-2 clustering (rho05): marker -> Stage-2 group ---
message("Loading canonical Stage-2 (rho05) clustering...")
s2 <- readRDS("module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05.rds")
g2 <- as.data.table(s2$groups)
marker2s2 <- g2[, .(marker = unlist(members)), by = .(s2_group = group_id)]
setkey(marker2s2, marker)
cat("Stage-2 (rho05) groups: ", nrow(g2), " covering ", nrow(marker2s2), " markers\n", sep = "")

## ---- one covariate: threshold Stage-1 units -> flag Stage-2 parent(s) ->
## plot every SNP at its OWN full-SNP-scan value, coloured by flagged group -
process_one <- function(unit_stat_file, unit_stat_col, thresh, thresh_label,
                        snp_stat_file, snp_stat_col, tag, title_stat) {
  message("\n[", tag, "] loading Stage-1-unit scan: ", unit_stat_file)
  s <- fread(unit_stat_file)
  stopifnot(nrow(s) == nrow(cl5))
  cl5s <- copy(cl5)
  cl5s[, stat := s[[unit_stat_col]][match(seq_len(.N), s$MRK)]]

  sig <- cl5s[stat >= thresh]
  message("[", tag, "] ", nrow(sig), " significant Stage-1 units (", thresh_label, ")")

  sig_s2 <- unique(marker2s2[.(sig$core_snp), on = "marker", nomatch = NULL]$s2_group)
  message("[", tag, "] -> ", length(sig_s2), " Stage-2 (rho05) group(s) flagged")

  regions <- if (length(sig_s2)) g2[group_id %in% sig_s2, members] else list()
  n_region_snps <- sum(lengths(regions))
  message("[", tag, "] ", n_region_snps, " member SNPs across those Stage-2 groups will be coloured")

  ## ---- REAL per-SNP y-value, from the full-genome BayPass scan -----------
  message("[", tag, "] loading full-SNP scan: ", snp_stat_file)
  snp_s <- fread(snp_stat_file)
  stopifnot(nrow(snp_s) == nrow(map_hyb_005))
  y <- setNames(snp_s[[snp_stat_col]][match(seq_len(nrow(map_hyb_005)), snp_s$MRK)], map_hyb_005$marker)

  p <- ld_manhattan(
    map = map_hyb_005, value = y, value_label = title_stat,
    regions = regions, hline = thresh,
    title = sprintf("%s: every SNP (own full-SNP BF/C2), coloured by Stage-2 (rho05) cluster containing a significant Stage-1 unit", tag),
    point_size = 0.6
  ) + labs(subtitle = sprintf(
    "Significance decided at Stage-1-unit resolution (%d tested, %d significant, %s); Stage-2 used only to describe physical extent -> %d group(s), %d SNPs coloured; y-axis is the full per-SNP scan",
    nrow(cl5s), nrow(sig), thresh_label, length(sig_s2), n_region_snps))

  outpng <- file.path(FIGDIR, sprintf("moduleB_stage1_%s_snp_manhattan.png", tag))
  ggsave(outpng, p, width = 22, height = 4.5, dpi = 200, limitsize = FALSE)
  message("[", tag, "] wrote ", outpng)
  invisible(p)
}

process_one(file.path(UNIT_DIR, "PC1_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15, "BF(dB)>=15",
           file.path(SNP_DIR, "PC1_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)",
           "PC1", "BF(dB)")
process_one(file.path(UNIT_DIR, "PC2_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15, "BF(dB)>=15",
           file.path(SNP_DIR, "PC2_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)",
           "PC2", "BF(dB)")
process_one(file.path(UNIT_DIR, "mito_C2_S1units_summary_contrast.out"), "log10(1/pval)", 3, "-log10(p)>=3",
           file.path(SNP_DIR, "mito_C2_fullSNP_stage1Omega_summary_contrast.out"), "log10(1/pval)",
           "mitoC2", "C2 -log10(p)")
process_one(file.path(UNIT_DIR, "bio6_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15, "BF(dB)>=15",
           file.path(SNP_DIR, "bio6_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)",
           "bio6", "BF(dB)")
process_one(file.path(UNIT_DIR, "bio11_S1units_withOmega_summary_betai_reg.out"), "BF(dB)", 15, "BF(dB)>=15",
           file.path(SNP_DIR, "bio11_fullSNP_stage1Omega_summary_betai_reg.out"), "BF(dB)",
           "bio11", "BF(dB)")

message("\n[moduleB-stage1-snp-manhattan] done")
