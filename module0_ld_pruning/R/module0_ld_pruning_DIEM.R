library(ggplot2)
library(igraph)
library(data.table)
library(SNPRelate)
library(parallel)
devtools::load_all("~/gitlab/LDscnR/")

## --- warm-run safeguard (COLD-RUN GENERATOR) --------------------------------
## Heavy from-scratch clustering: no internal caching, recomputes
## ld_prune_and_eMLG() ~8x over ~1.1M markers. Refuse to re-run when the
## canonical clustering already exists.
FORCE_COLD <- FALSE   # set TRUE only for a deliberate from-scratch regeneration
if (!FORCE_COLD && file.exists("module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds")) {
  stop("module0_ld_pruning_DIEM.R is a COLD-RUN generator and its clustering output\n",
       "  (module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds) already exists; a re-run would\n",
       "  recompute the two-stage clustering from scratch (~hours). Remove that file and\n",
       "  set FORCE_COLD <- TRUE to regenerate.", call. = FALSE)
}

# ------------------------------------------------------------
# Read in data (generated in LD_decay_from_DIEM.R)
# ------------------------------------------------------------
message("=== Loading data and Creating gds ======")
# loads GTs_hybrids_005,map_hyb_005, ld_decay and sample_data from LD_decay_from_DIEM.R
# includes maf and ld_w_095
load("data/hybrids_only_maf005.Rdata")

## alias used by the diagnostic plot_pruning_comparison() calls at the end of
## this script; the loaded .Rdata provides map_hyb_005, not `map`.
map <- copy(map_hyb_005)

gds_hyb <- create_gds_from_geno(geno=GTs_hybrids_005, map_hyb_005, "gds_hybrids.gds")


message("=== Pruning SNPs ======")

## Stage 1 on the WHOLE marker set together -- cheap regardless of scale
## (full Chr26: 8.3s for 26,846 markers), and crucially uses LD_decay$el's
## real window-covered edges to group markers correctly. Splitting markers
## by ld_w BEFORE clustering (an earlier version of this) artificially
## severs real blocks right at the threshold boundary: checked on Chr26
## alone, 226 of 15,524 Stage-1 clusters (formed from a combined run) mixed
## low- and high-ld_w members, together 4,595 markers of which 975 were
## low-ld_w members that a pre-split would have cut loose from their real
## cluster.
ld_w_threshold <- 0.2
pruned_stage1 <- ld_complexity_reduction(
  map = map_hyb_005, LD_decay = ld_decay, rho = 0.5, cores = 1
)

## Combined LD-pruning + eMLG generation in one pass, replacing the earlier
## cluster-flagging -> merge_ld_clusters() -> separate eMLG-extraction ->
## separate dynamic-cut sequence (two independent all-pairs correlation
## passes over overlapping data) with ONE distance-restricted, quality-
## gated dynamic cut directly on Stage 1's flagged clusters. Now part of
## LDscnR (ld_prune_and_eMLG()/dynamic_cut_eMLG()) -- see their roxygen
## docs there for the full rationale: average vs single vs complete
## linkage comparison, the score_eMLG/pair_r2 quality gate (with both bugs
## it caught along the way), why the distance restriction uses
## consecutive-gap (not total-span) semantics, and why this is defensible
## for LD-pruning specifically in this young, low-recombination hybrid
## population -- not just for eMLG summarization.
result <- ld_prune_and_eMLG(
  GTs = GTs_hybrids_005, stage1 = pruned_stage1, ld_w_col = "ld_w_095",
  ld_w_threshold = ld_w_threshold, score_threshold = 0.80, min_r2 = 0.2,
  distance_threshold = 5e5,compute_unflagged_eMLG = FALSE
)

pruned_markers <- result$pruned
eMLG           <- result$eMLG
eMLG_groups    <- result$groups

message(
  "Keeping ", length(pruned_markers), "/", map_hyb_005[TRUE,.N], " (",
  round(100 * length(pruned_markers) / map_hyb_005[TRUE,.N], 2), "%) SNPs"
)
saveRDS(pruned_markers, "module0_ld_pruning/data/pruned_markers.rds")
saveRDS(pruned_stage1, "module0_ld_pruning/data/pruned_stage1.rds")
saveRDS(list(eMLG = eMLG, groups = eMLG_groups), "module0_ld_pruning/data/eMLG_groups.rds")

# ------------------------------------------------------------
# Per-chromosome diagnostic: Stage 1 (fragmented) vs the combined
# ld_prune_and_eMLG() result, for the ld_w>threshold markers -- where
# fragmented centromeric/inversion blocks are most likely. Reuses the
# whole-genome objects already computed above, no independent
# recomputation. plot_pruning_comparison() below works for any chromosome;
# Chr26 was the original worked example throughout development (it has
# the most pronounced single low-recombination block of all 26
# chromosomes, ~23% of markers exceed ld_w_095>0.2).
# ------------------------------------------------------------
library(patchwork)

## cycle a modest, visually-distinct palette across CL_id -- with hundreds
## of clusters no palette gives every one a unique color, but neighbouring
## clusters (what we care about here) will very likely differ

# test with min_loci>=5

results_min_loci5 <- lapply(c(0.1,0.15,0.2),function(th){
  out <- ld_prune_and_eMLG(
    GTs = GTs_hybrids_005, stage1 = pruned_stage1, ld_w_col = "ld_w_095",
    ld_w_threshold = th, score_threshold = 0.80, min_r2 = 0.2,
    distance_threshold = 5e5,compute_unflagged_eMLG = FALSE,min_n_loci_eMLG = 5
  )
  out$groups[,th:=th]
  return(out)
})

results_min_loci5_2 <- lapply(c(0.05,0.025),function(th){
  out <- ld_prune_and_eMLG(
    GTs = GTs_hybrids_005, stage1 = pruned_stage1, ld_w_col = "ld_w_095",
    ld_w_threshold = th, score_threshold = 0.80, min_r2 = 0.2,
    distance_threshold = 5e5,compute_unflagged_eMLG = FALSE,min_n_loci_eMLG = 5
  )
  out$groups[,th:=th]
  return(out)
})

#c(results_min_loci5,results_min_loci5_2)
saveRDS(c(results_min_loci5,results_min_loci5_2),"module0_ld_pruning/data/results_min_loci5.rds")

results_min_loci5_2[[2]]$groups[TRUE,hist(score)]

result_01_min_loci5 <- ld_prune_and_eMLG(
  GTs = GTs_hybrids_005, stage1 = pruned_stage1, ld_w_col = "ld_w_095",
  ld_w_threshold = 0.1, score_threshold = 0.80, min_r2 = 0.2,
  distance_threshold = 5e5,compute_unflagged_eMLG = FALSE,min_n_loci_eMLG = 5
)
#result_01_min_loci5$groups

pruned_markers <- results_min_loci5_2[[2]]$pruned
eMLG           <- results_min_loci5_2[[2]]$eMLG
eMLG_groups    <- results_min_loci5_2[[2]]$groups


#result_01_min_loci5$groups[,hist(score)]
## plot_pruning_comparison() is now part of LDscnR (see its roxygen docs
## there for the full parameter rationale, including why min_n_loci_flag
## must be passed through to keep the Stage 1 and Combined panels in sync).
plot_pruning_comparison("Chr1", pruned_stage1, results_min_loci5_2[[2]], map_hyb_005, ld_w_threshold = 0.025)
plot_pruning_comparison("Chr1", pruned_stage1, results_min_loci5_2[[2]], map_hyb_005, ld_w_threshold = 0.025, direction = "low")

results_min_loci5_2[[2]]$group[grepl("U",group_id) & n_loci>=5]

# ------------------------------------------------------------
# Final filtering and eMLG construction
# ------------------------------------------------------------
# distance_threshold is genetic (cM), not physical (bp): physical distance
# is a poor, sometimes wildly misleading proxy for recombination distance
# -- checked directly on this data, some flagged-cluster pairs over 1 Mb
# apart physically sat at cM distance 0 (fully linked), while others only
# a few hundred kb apart already spanned several cM. cM_threshold = 0.5 is
# the settled choice: it reproduces the old flat distance_threshold = 5e5's
# behaviour almost exactly on typical chromosomes (real, adaptive
# protection where a flat bp value would just be a no-op), while
# additionally -- correctly -- staying out of the way of genuinely
# exceptional, very-low-recombination chromosomes like Chr26 (its own
# slower recombination rate gives 0.5 cM a proportionally larger physical
# span there, rather than clipping its one real block down to match a
# threshold calibrated for a typical chromosome). Checked directly against
# three alternatives before settling here: a flat bp constant (adequate
# down to ~100kb, starts degrading the largest blocks by 50kb), an
# a_pred/rho-derived bp threshold (rho=0.95 was far too strict, ~11kb,
# collapsed every large block tested), and predicting cM from
# locally-windowed decay rate a instead of using the real map (worse than
# just using a flat bp constant -- summing a noisy per-window prediction
# over a long block amplifies error badly). See LDscnR's
# ld_prune_and_eMLG()/interpolate_cM() roxygen docs for the parameter
# mechanics.
rec_map <- fread("./data/Frufa_DTOL_PR.ref_genome.recmap")
rec_map[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
genetic_map <- rec_map[, .(Chr, Pos = pos, cM)]

## two variants one with 0.5cM
eMLG_5loci_0025 <- ld_prune_and_eMLG(
  GTs = GTs_hybrids_005, stage1 = pruned_stage1, ld_w_col = "ld_w_095",
  ld_w_threshold = 0.025, score_threshold = 0.80, min_r2 = 0.2, min_n_loci_flag = 5,
  genetic_map = genetic_map, cM_threshold = 0.5,
  compute_unflagged_eMLG = TRUE, min_n_loci_eMLG = 5
)

saveRDS(eMLG_5loci_0025,"module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds")

## --- best-SNP representative companion (canonical eMLG representation) --------
## Instead of the averaged consensus, represent every cluster by a single real
## SNP: eMLG_best_snp() picks, per has_eMLG block (>=5 loci), the member whose
## genotype best matches the consensus (|r|), returning that SNP's real calls
## (missing filled from consensus) in $geno and its id in $stats$best_marker --
## so cluster DI/position come from a real SNP, never averaged. Clusters below
## the eMLG threshold (n_loci < 5) have no consensus and are already single
## SNPs, so they keep their centrality representative. $rep_snp_all gives the
## unified per-cluster real-SNP mapping (best_marker where has_eMLG, else
## representative) that downstream modules (A cluster-level, B/C/D) join on.
eMLG_best_cM05 <- eMLG_best_snp(eMLG_5loci_0025, GTs_hybrids_005,
                               fill = TRUE, round_fill = TRUE)
local({
  g  <- data.table::as.data.table(eMLG_5loci_0025$groups)
  g[, has_eMLG := as.logical(has_eMLG)]
  allrep <- g[, .(group_id, representative, n_loci, has_eMLG)]
  allrep[data.table::as.data.table(eMLG_best_cM05$stats)[, .(group_id, best_marker)],
         on = "group_id", best_marker := i.best_marker]
  allrep[, rep_snp := data.table::fifelse(has_eMLG & !is.na(best_marker),
                                          best_marker, representative)]
  eMLG_best_cM05$rep_snp_all <<- allrep[]
})
saveRDS(eMLG_best_cM05, "module0_ld_pruning/data/eMLG_5loci_0025_cM05_bestsnp.rds")

## one with 1cM
eMLG_5loci_0025_cM1 <- ld_prune_and_eMLG(
  GTs = GTs_hybrids_005, stage1 = pruned_stage1, ld_w_col = "ld_w_095",
  ld_w_threshold = 0.025, score_threshold = 0.80, min_r2 = 0.2, min_n_loci_flag = 5,
  genetic_map = genetic_map, cM_threshold = 1,
  compute_unflagged_eMLG = TRUE, min_n_loci_eMLG = 5
)

saveRDS(eMLG_5loci_0025_cM1,"module0_ld_pruning/data/eMLG_5loci_0025_cM1.rds")

## plotting
eMLG_5loci_0025_cM1 <- readRDS("module0_ld_pruning/data/eMLG_5loci_0025_cM1.rds")
pruned_stage1 <- readRDS("module0_ld_pruning/data/pruned_stage1.rds")
plot_pruning_comparison("Chr26", pruned_stage1, eMLG_5loci_0025_cM1, map, out_folder = "module0_ld_pruning/Figures/")
plot_pruning_comparison("Chr26", pruned_stage1, eMLG_5loci_0025_cM1, map, out_folder = "module0_ld_pruning/Figures/",direction = "low")

eMLG_5loci_0025_cM1$groups[,.N]
eMLG_5loci_0025_cM1$groups[n_loci>=5,.N]
