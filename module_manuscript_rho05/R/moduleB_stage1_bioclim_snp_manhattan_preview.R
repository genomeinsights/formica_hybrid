## =========================================================
## module_manuscript_rho05 -- bio6/bio11 Manhattan, Stage-1-unit PREVIEW
## =========================================================
## Quick look using only the already-completed Stage-1-unit scan (no full-
## SNP scan needed yet -- that's running separately, ~2.5h/covariate). Same
## caveat as the original PC1/PC2/mito-C2 attempt before the fix: every SNP
## in a Stage-1 cluster is plotted at that cluster's ONE tested value (not an
## independent per-SNP statistic), so the y-axis is Stage-1-unit-resolution,
## not truly per-SNP. Colouring (which Stage-2 (rho05) cluster) IS already
## correct/final -- that part doesn't depend on which scan the y-axis uses.
## Superseded by moduleB_stage1_snp_manhattan_ldmanhattan.R's approach once
## the full-SNP bio6/bio11 scan finishes.
##
## Run from the repo root: Rscript module_manuscript_rho05/R/moduleB_stage1_bioclim_snp_manhattan_preview.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2) })
devtools::load_all("~/gitlab/LDscnR/")

UNIT_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded_S1units"
FIGDIR <- "module_manuscript_rho05/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

load("data/hybrids_only_maf005.Rdata")   # map_hyb_005
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

process_one <- function(unit_stat_file, tag, title_stat) {
  message("\n[", tag, "] loading ", unit_stat_file)
  s <- fread(unit_stat_file)
  stopifnot(nrow(s) == nrow(cl5))
  cl5s <- copy(cl5)
  cl5s[, stat := s[["BF(dB)"]][match(seq_len(.N), s$MRK)]]

  sig <- cl5s[stat >= 15]
  message("[", tag, "] ", nrow(sig), " significant Stage-1 units (BF(dB)>=15)")

  sig_s2 <- unique(marker2s2[.(sig$core_snp), on = "marker", nomatch = NULL]$s2_group)
  message("[", tag, "] -> ", length(sig_s2), " Stage-2 (rho05) group(s) flagged")
  regions <- if (length(sig_s2)) g2[group_id %in% sig_s2, members] else list()
  n_region_snps <- sum(lengths(regions))
  message("[", tag, "] ", n_region_snps, " member SNPs across those Stage-2 groups will be coloured")

  ## PREVIEW y-value: inherited from each SNP's own Stage-1 cluster (not yet
  ## the genuine per-SNP full-scan value -- see header caveat)
  y <- setNames(rep(cl5s$stat, cl5s$n_snps), unlist(cl5s$members))

  p <- ld_manhattan(
    map = map_hyb_005, value = y, value_label = title_stat,
    regions = regions, hline = 15,
    title = sprintf("%s (PREVIEW, Stage-1-unit resolution): coloured by Stage-2 (rho05) cluster containing a significant Stage-1 unit", tag),
    point_size = 0.8
  ) + labs(subtitle = sprintf(
    "%d Stage-1 units tested; %d significant (BF(dB)>=15) -> %d Stage-2 group(s), %d SNPs coloured; y-axis is each SNP's Stage-1-cluster value (unit-resolution) -- full-SNP scan running separately",
    nrow(cl5s), nrow(sig), length(sig_s2), n_region_snps))

  outpng <- file.path(FIGDIR, sprintf("moduleB_stage1_%s_snp_manhattan_PREVIEW.png", tag))
  ggsave(outpng, p, width = 22, height = 4.5, dpi = 200, limitsize = FALSE)
  message("[", tag, "] wrote ", outpng)
  invisible(p)
}

process_one(file.path(UNIT_DIR, "bio6_S1units_withOmega_summary_betai_reg.out"), "bio6", "BF(dB)")
process_one(file.path(UNIT_DIR, "bio11_S1units_withOmega_summary_betai_reg.out"), "bio11", "BF(dB)")

message("\n[moduleB-stage1-bioclim-manhattan-preview] done")
