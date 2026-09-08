## =========================================================
## module_di25_rho05 -- Stage 2 (ld_prune_and_eMLG) for the high-DI
## (DI > -25) markers under min_r2_rho = 0.5, as a paired alternative to
## module_di25's fixed min_r2 = 0.2 (canonical di25_clustering_cM5.rds).
##
## Stage 1 (rho = 0.5 join) and the LD-decay fit do NOT depend on min_r2 --
## both are reused UNCHANGED from module_di25/data/, so this is a paired
## comparison against the canonical run (min_r2 is the only thing that
## differs). All other settings match di25_ld_clustering.R's canonical
## cM = 5 call exactly.
##
## Run from the repo root: Rscript module_di25_rho05/R/di25_ld_clustering_rho05.R
## =========================================================

suppressMessages({
  library(data.table)
  library(igraph)
  library(SNPRelate)
})
devtools::load_all("~/gitlab/LDscnR/")

OUT <- "module_di25_rho05/data"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

CM_THRESHOLD    <- 5      # canonical cap (di25_clustering_cM5.rds)
LD_W_THRESHOLD  <- 0      # with min_n_loci_flag = 1 this flags every cluster
MIN_N_LOCI_FLAG <- 1
MIN_N_LOCI_EMLG <- 3
SCORE_THRESHOLD <- 0.80
CORES           <- 4

OUT_RDS <- file.path(OUT, "di25_clustering_cM5_rho05.rds")
if (file.exists(OUT_RDS)) {
  stop("di25_ld_clustering_rho05.R's output already exists (", OUT_RDS, ").", call. = FALSE)
}

message("[di25-rho05] loading cached Stage 1 + decay (shared, unaffected by min_r2)")
di25_stage1_obj <- readRDS("module_di25/data/di25_stage1.rds")   # $stage1, $map
di25_decay      <- readRDS("module_di25/data/di25_ld_decay.rds")
di25_inputs     <- readRDS("module_di25/data/di25_inputs.rds")   # $GTs_hyb

rec_map <- fread("./data/Frufa_DTOL_PR.ref_genome.recmap")
rec_map[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
genetic_map <- rec_map[, .(Chr, Pos = pos, cM)]

message(sprintf("[di25-rho05] Stage 2 @ %.1f cM, min_r2_rho = 0.5", CM_THRESHOLD))
t0 <- Sys.time()
res <- ld_prune_and_eMLG(
  GTs = di25_inputs$GTs_hyb, stage1 = di25_stage1_obj$stage1, ld_w_col = "ld_w_095",
  ld_w_threshold = LD_W_THRESHOLD, min_n_loci_flag = MIN_N_LOCI_FLAG,
  score_threshold = SCORE_THRESHOLD, min_r2 = NULL, min_r2_rho = 0.5, LD_decay = di25_decay,
  genetic_map = genetic_map, cM_threshold = CM_THRESHOLD,
  compute_unflagged_eMLG = TRUE, min_n_loci_eMLG = MIN_N_LOCI_EMLG,
  cores = CORES
)
message("[di25-rho05] Stage 2 elapsed: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
saveRDS(res, OUT_RDS)

g <- res$groups
message(sprintf(
  "groups=%d eMLG=%d singletons=%d max_n_loci=%d markers_in=%d reduction=%.2f%%",
  nrow(g), sum(g$has_eMLG), sum(g$n_loci == 1), max(g$n_loci), sum(g$n_loci),
  100 * (1 - nrow(g) / sum(g$n_loci))
))
message("[di25-rho05] done")
