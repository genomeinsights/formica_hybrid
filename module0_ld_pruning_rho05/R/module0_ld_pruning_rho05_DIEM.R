## =========================================================
## module0_ld_pruning_rho05 -- genome-wide Stage 2 (ld_prune_and_eMLG) under
## the decay-relative min_r2_rho = 0.5 default, as a paired alternative to
## module0_ld_pruning's fixed min_r2 = 0.2 (canonical eMLG_5loci_0025_cM05.rds).
##
## Stage 1 (ld_complexity_reduction, rho = 0.5 join) and the LD-decay fit do
## NOT depend on min_r2 at all -- both are reused UNCHANGED from
## module0_ld_pruning/data/, so the two arms are a genuinely paired
## comparison (one Stage 1 partition, one decay fit -- min_r2 is the only
## thing that differs), not two independent draws.
##
## All other settings (ld_w_threshold, score_threshold, min_n_loci_flag,
## cM_threshold = 0.5) are held at module0_ld_pruning's canonical values --
## this varies min_r2 alone.
##
## Run from the repo root: Rscript module0_ld_pruning_rho05/R/module0_ld_pruning_rho05_DIEM.R
## =========================================================

library(ggplot2)
library(igraph)
library(data.table)
library(SNPRelate)
library(parallel)
devtools::load_all("~/gitlab/LDscnR/")

FORCE_COLD <- FALSE
OUT <- "module0_ld_pruning_rho05/data"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create("module0_ld_pruning_rho05/Figures", showWarnings = FALSE, recursive = TRUE)

OUT_RDS <- file.path(OUT, "eMLG_5loci_0025_cM05_rho05.rds")
if (!FORCE_COLD && file.exists(OUT_RDS)) {
  stop("module0_ld_pruning_rho05_DIEM.R's output already exists (", OUT_RDS, ");\n",
       "  set FORCE_COLD <- TRUE to regenerate.", call. = FALSE)
}

message("=== Loading canonical Stage 1 + decay (SHARED, unaffected by min_r2) ===")
load("data/hybrids_only_maf005.Rdata")     # GTs_hybrids_005, map_hyb_005, ld_decay, sample_data
pruned_stage1 <- readRDS("module0_ld_pruning/data/pruned_stage1.rds")

rec_map <- fread("./data/Frufa_DTOL_PR.ref_genome.recmap")
rec_map[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
genetic_map <- rec_map[, .(Chr, Pos = pos, cM)]

message("=== Stage 2: genome-wide ld_prune_and_eMLG, min_r2_rho = 0.5 ===")
t0 <- Sys.time()
eMLG_5loci_0025_rho05 <- ld_prune_and_eMLG(
  GTs = GTs_hybrids_005, stage1 = pruned_stage1, ld_w_col = "ld_w_095",
  ld_w_threshold = 0.025, score_threshold = 0.80,
  min_r2 = NULL, min_r2_rho = 0.5, LD_decay = ld_decay,
  min_n_loci_flag = 5, genetic_map = genetic_map, cM_threshold = 0.5,
  compute_unflagged_eMLG = TRUE, min_n_loci_eMLG = 5
)
message("Stage 2 elapsed: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
saveRDS(eMLG_5loci_0025_rho05, OUT_RDS)

message("=== best-SNP representative companion (canonical eMLG representation) ===")
eMLG_best_rho05 <- eMLG_best_snp(eMLG_5loci_0025_rho05, GTs_hybrids_005,
                                  fill = TRUE, round_fill = TRUE)
local({
  g  <- data.table::as.data.table(eMLG_5loci_0025_rho05$groups)
  g[, has_eMLG := as.logical(has_eMLG)]
  allrep <- g[, .(group_id, representative, n_loci, has_eMLG)]
  allrep[data.table::as.data.table(eMLG_best_rho05$stats)[, .(group_id, best_marker)],
         on = "group_id", best_marker := i.best_marker]
  allrep[, rep_snp := data.table::fifelse(has_eMLG & !is.na(best_marker),
                                          best_marker, representative)]
  eMLG_best_rho05$rep_snp_all <<- allrep[]
})
saveRDS(eMLG_best_rho05, file.path(OUT, "eMLG_5loci_0025_cM05_rho05_bestsnp.rds"))

message("=== Chr26 Stage1+2 diagnostic figure ===")
map <- copy(map_hyb_005)
plot_pruning_comparison("Chr26", pruned_stage1, eMLG_5loci_0025_rho05, map,
                         out_folder = "module0_ld_pruning_rho05/Figures/")

g <- eMLG_5loci_0025_rho05$groups
message(sprintf(
  "groups=%d eMLG=%d singletons=%d max_n_loci=%d markers_in=%d reduction=%.2f%%",
  nrow(g), sum(g$has_eMLG), sum(g$n_loci == 1), max(g$n_loci), sum(g$n_loci),
  100 * (1 - nrow(g) / sum(g$n_loci))
))
message("Done.")
