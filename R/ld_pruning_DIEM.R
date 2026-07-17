library(ggplot2)
library(igraph)
library(data.table)
library(SNPRelate)
library(parallel)
devtools::load_all("~/gitlab/LDscnR/")

# ------------------------------------------------------------
# Read in data (generated in LD_decay_from_DIEM.R)
# ------------------------------------------------------------

message("=== Loading data and Creating gds ======")
# loads GTs_hybrids_005,map_hyb_005, ld_decay and sample_data from LD_decay_from_DIEM.R
# includes maf and ld_w_095
load("./data/hybrids_only_maf005.Rdata")

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
saveRDS(pruned_markers, "./data/pruned_markers.rds")
saveRDS(list(eMLG = eMLG, groups = eMLG_groups), "./data/eMLG_groups.rds")

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
saveRDS(c(results_min_loci5,results_min_loci5_2),"./data/results_min_loci5.rds")

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

eMLG_5loci_0025 <- ld_prune_and_eMLG(
  GTs = GTs_hybrids_005, stage1 = pruned_stage1, ld_w_col = "ld_w_095",
  ld_w_threshold = 0.025, score_threshold = 0.80, min_r2 = 0.2,min_n_loci_flag=5,
  distance_threshold = 5e5,compute_unflagged_eMLG = TRUE,min_n_loci_eMLG = 5
)
eMLG_5loci_0025$groups
saveRDS(eMLG_5loci_0025,"./data/eMLG_5loci_0025.rds")

plot_pruning_comparison("Chr1", pruned_stage1, eMLG_5loci_0025, map_hyb_005, ld_w_threshold = 0.025, min_n_loci_flag = 5)
plot_pruning_comparison("Chr1", pruned_stage1, eMLG_5loci_0025, map_hyb_005, ld_w_threshold = 0.025, min_n_loci_flag = 5, direction = "low")

eMLG_5loci_0025$groups[n_loci>=10,table(has_eMLG,ifelse(grepl("F",group_id),"Flagged","Unflagged"))]

eMLG_5loci_0025$group

