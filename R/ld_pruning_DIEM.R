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
## gated dynamic cut directly on Stage 1's flagged clusters. See
## dev/R/ld_prune_and_eMLG.R for the full rationale: average vs single vs
## complete linkage comparison, the score_eMLG/pair_r2 quality gate (with
## both bugs it caught along the way), why the distance restriction uses
## consecutive-gap (not total-span) semantics, and why this is defensible
## for LD-pruning specifically in this young, low-recombination hybrid
## population -- not just for eMLG summarization.
source("./dev/R/ld_prune_and_eMLG.R")

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
# Chr26 diagnostic: Stage 1 (fragmented) vs the combined ld_prune_and_eMLG()
# result, for the ld_w>0.2 markers -- where fragmented centromeric/
# inversion blocks are most likely. Reuses the whole-genome objects already
# computed above, no independent recomputation.
# ------------------------------------------------------------
library(patchwork)

## cycle a modest, visually-distinct palette across CL_id -- with hundreds
## of clusters no palette gives every one a unique color, but neighbouring
## clusters (what we care about here) will very likely differ

# test with min_loci>=5
result_min_loci5 <- ld_prune_and_eMLG(
  GTs = GTs_hybrids_005, stage1 = pruned_stage1, ld_w_col = "ld_w_095",
  ld_w_threshold = ld_w_threshold, score_threshold = 0.80, min_r2 = 0.2,
  distance_threshold = 5e5,compute_unflagged_eMLG = FALSE,min_n_loci_eMLG = 5
)

pruned_markers <- result_min_loci5$pruned
eMLG           <- result_min_loci5$eMLG
eMLG_groups    <- result_min_loci5$groups

pal_cluster <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33",
                  "#A65628","#F781BF","#1B9E77","#D95F02","#7570B3","#66A61E")

plot_clusters_chr26 <- function(map_snp, title) {
  dt <- copy(map_snp)
  dt[, cl_rank := match(CL_id, unique(CL_id))]
  dt[, col := pal_cluster[(cl_rank - 1) %% length(pal_cluster) + 1]]
  ggplot(dt, aes(Pos / 1e6, ld_w_095, color = col)) +
    geom_point(size = 1.2, alpha = 0.85) +
    scale_color_identity() +
    theme_bw(base_size = 13) +
    labs(x = "Chr26 position (Mbp)", y = expression(ld["w,"*rho*"=0.95"]), title = title)
}

needs_merge_ids <- pruned_stage1$map_snp[ld_w_095 > ld_w_threshold, unique(CL_id)]
chr26_stage1_snp <- pruned_stage1$map_snp[Chr == "Chr26" & CL_id %in% needs_merge_ids]

message("=== Chr26 diagnostic, Stage 1: ld_complexity_reduction (flagged clusters only) ===")
message("Stage 1: ", uniqueN(chr26_stage1_snp$CL_id), " clusters")

p_chr26_stage1 <- plot_clusters_chr26(
  chr26_stage1_snp,
  sprintf("Stage 1: ld_complexity_reduction() -- %d clusters", uniqueN(chr26_stage1_snp$CL_id))
)
ggsave("./Figures/chr26_stage1_ld_complexity_reduction.png", p_chr26_stage1, width = 10, height = 5)

## combined result restricted to Chr26's FLAGGED groups only (group_id
## prefix "F", vs "U" for unflagged pass-through) -- matches the Stage 1
## diagnostic's scope above (also flagged-clusters-only); eMLG_groups
## otherwise includes thousands of tiny unflagged low-ld_w singletons
## genome-wide that would make this comparison apples-to-oranges
chr26_groups <- eMLG_groups[Chr == "Chr26" & startsWith(group_id, "F")]
chr26_final_snp <- chr26_groups[TRUE, .(marker = unlist(members)), by = group_id]
setnames(chr26_final_snp, "group_id", "CL_id")
chr26_final_snp <- map_hyb_005[Chr == "Chr26", .(marker, Pos, ld_w_095)][chr26_final_snp, on = "marker"]

message("=== Chr26 diagnostic, ld_prune_and_eMLG (combined) ===")
message("Combined: ", uniqueN(chr26_final_snp$CL_id), " groups")

p_chr26_combined <- plot_clusters_chr26(
  chr26_final_snp,
  sprintf("ld_prune_and_eMLG() -- %d groups", uniqueN(chr26_final_snp$CL_id))
)
ggsave("./Figures/chr26_combined_ld_prune_and_eMLG.png", p_chr26_combined, width = 10, height = 5)

p_chr26_compare <- p_chr26_stage1 / p_chr26_combined
ggsave("./Figures/chr26_stage1_vs_combined.png", p_chr26_compare, width = 10, height = 9)
p_chr26_compare

## eMLG correlation heatmap for Chr26's flagged/merged groups -- subset of
## the whole-genome eMLG matrix already computed above, no recomputation
chr26_group_ids <- chr26_groups$group_id
chr26_eMLG <- eMLG[, chr26_group_ids, drop = FALSE]
message("Chr26 eMLG matrix: ", nrow(chr26_eMLG), " individuals x ", ncol(chr26_eMLG), " groups")

R2_chr26 <- suppressWarnings(cor(chr26_eMLG, use = "pairwise.complete.obs")^2)
R2_chr26[!is.finite(R2_chr26)] <- 0

plot_eMLG_heatmap <- function(R2_mat, title) {
  dt <- as.data.table(as.table(R2_mat))
  setnames(dt, c("Var1", "Var2", "r2"))
  dt[, Var1 := factor(Var1, levels = colnames(R2_mat))]
  dt[, Var2 := factor(Var2, levels = colnames(R2_mat))]
  ggplot(dt, aes(Var1, Var2, fill = r2)) +
    geom_tile() +
    scale_fill_viridis_c(name = expression(r^2), limits = c(0, 1)) +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
    labs(x = "Group", y = "Group", title = title)
}

pos_lookup <- setNames(map_hyb_005$Pos, map_hyb_005$marker)
group_pos <- vapply(chr26_groups$members, function(mk) mean(pos_lookup[mk]), numeric(1))
ord_pos <- order(group_pos)

p_eMLG_pos <- plot_eMLG_heatmap(
  R2_chr26[ord_pos, ord_pos],
  sprintf("Chr26 eMLG correlation -- %d groups, position order", ncol(R2_chr26))
)
ggsave("./Figures/chr26_eMLG_heatmap_posorder.png", p_eMLG_pos, width = 8, height = 7, dpi = 150)
