## =========================================================
## module_manuscript_rho05 -- Fig [complexity_reduction]: Stage 1 vs combined
## Stage 1+2, chromosome 26, full data set, min_r2_rho = 0.5
## =========================================================
## Mirrors module0_ld_pruning_DIEM.R's canonical diagnostic call exactly
## (ld_w_threshold = 0.025, score_threshold = 0.80, min_n_loci_flag = 5,
## cM_threshold = 1 -- this cM=1 variant, not the downstream cM=0.5 one, is
## what the canonical figure was built from), varying only min_r2.
##
## Restricted to a Chr26-only subset of Stage 1 for speed -- ld_prune_and_eMLG()
## merges strictly per-chromosome (confirmed in LDscnR/R/ld_prune_and_eMLG.R:
## `chr_levels <- unique(flagged$Chr)`), so this reproduces exactly what a
## genome-wide run would give for Chr26. Stage 1 + decay are reused UNCHANGED
## from module0_ld_pruning/data/ (min_r2 doesn't affect them) -- a paired
## comparison against the canonical (min_r2 = 0.2) figure.
##
## Run from the repo root: Rscript module_manuscript_rho05/R/module0_fig_complexity_reduction_rho05.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
devtools::load_all("~/gitlab/LDscnR/")

FIGDIR <- "module_manuscript_rho05/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

load("data/hybrids_only_maf005.Rdata")              # GTs_hybrids_005, map_hyb_005, ld_decay
pruned_stage1 <- readRDS("module0_ld_pruning/data/pruned_stage1.rds")
s1_chr26 <- list(map_snp = pruned_stage1$map_snp[Chr == "Chr26"],
                 clusters = pruned_stage1$clusters[Chr == "Chr26"])

rec_map <- fread("./data/Frufa_DTOL_PR.ref_genome.recmap")
rec_map[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
genetic_map <- rec_map[, .(Chr, Pos = pos, cM)]

res_rho05 <- ld_prune_and_eMLG(
  GTs = GTs_hybrids_005, stage1 = s1_chr26, ld_w_col = "ld_w_095",
  ld_w_threshold = 0.025, score_threshold = 0.80,
  min_r2 = NULL, min_r2_rho = 0.5, LD_decay = ld_decay,
  min_n_loci_flag = 5, genetic_map = genetic_map, cM_threshold = 1,
  compute_unflagged_eMLG = TRUE, min_n_loci_eMLG = 5
)

cat("min_r2 resolved (Chr26, cM=1):", signif(res_rho05$params$min_r2_resolved[["Chr26"]], 4), "\n")

map <- copy(map_hyb_005)
plot_pruning_comparison("Chr26", list(map_snp = pruned_stage1$map_snp, clusters = pruned_stage1$clusters),
                         res_rho05, map, out_folder = FIGDIR)

g <- res_rho05$groups
message(sprintf("Chr26 groups=%d eMLG=%d singletons=%d max_n_loci=%d markers_in=%d reduction=%.2f%%",
                 nrow(g), sum(g$has_eMLG), sum(g$n_loci==1), max(g$n_loci), sum(g$n_loci),
                 100*(1-nrow(g)/sum(g$n_loci))))
message("[complexity-reduction-rho05] done -> ", file.path(FIGDIR, "Chr26_stage1_vs_combined_high.png"))
