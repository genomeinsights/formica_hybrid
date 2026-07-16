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



result_01_min_loci5 <- ld_prune_and_eMLG(
  GTs = GTs_hybrids_005, stage1 = pruned_stage1, ld_w_col = "ld_w_095",
  ld_w_threshold = 0.1, score_threshold = 0.80, min_r2 = 0.2,
  distance_threshold = 5e5,compute_unflagged_eMLG = FALSE,min_n_loci_eMLG = 5
)


pruned_markers <- result_01_min_loci5$pruned
eMLG           <- result_01_min_loci5$eMLG
eMLG_groups    <- result_01_min_loci5$groups

eMLG_groups[,hist(score)]

#' Plot Stage 1 vs combined (ld_prune_and_eMLG) cluster diagnostic for one
#' chromosome, saving only the stacked comparison figure (p_chrXX_compare).
#'
#' @param chr Chromosome name (e.g. "Chr26"), used in plot labels/filename.
#' @param pruned_stage1 Output of ld_complexity_reduction() (whole-genome).
#' @param result Output of ld_prune_and_eMLG() (whole-genome).
#' @param map Map data.table with Chr, Pos, marker, and `ld_w_col`.
#' @param ld_w_col Name of the ld_w column in `map`/`pruned_stage1$map_snp`.
#' @param ld_w_threshold Threshold used to flag clusters for the Stage-1
#'   comparison view -- should match what was used to build `result`.
#' @param direction "high" (default) shows the FLAGGED clusters/groups
#'   (ld_w_col > threshold, group_id prefix "F") -- the main diagnostic view.
#'   "low" shows the UNFLAGGED side (ld_w_col < threshold, group_id prefix
#'   "U") as a sanity check: this view should never contain markers with
#'   ld_w_col > threshold, since "unflagged" is defined as no cluster member
#'   exceeding threshold (verified empirically -- see conversation with
#'   Claude on 2026-07-16). The "high" view, by contrast, is expected to
#'   contain some sub-threshold "boundary" markers (cluster-mates of a
#'   flagged marker) since flagging happens at the cluster level, not the
#'   marker level.
#' @param out_folder Folder to write the figure to (created if missing).
#'
#' @return The stacked comparison plot, invisibly (also saved to disk as
#'   `<out_folder>/<chr>_stage1_vs_combined_<direction>.png`).
plot_pruning_comparison <- function(chr, pruned_stage1, result, map,
                                     ld_w_col = "ld_w_095", ld_w_threshold = 0.2,
                                     direction = c("high", "low"),
                                     out_folder = "./Figures/", width = 10, height = 9) {
  direction <- match.arg(direction)
  dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)

  pal_cluster <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33",
                    "#A65628","#F781BF","#1B9E77","#D95F02","#7570B3","#66A61E")

  ## cycle a modest, visually-distinct palette across CL_id -- with hundreds
  ## of clusters no palette gives every one a unique color, but neighbouring
  ## clusters (what we care about here) will very likely differ
  plot_clusters <- function(map_snp, title) {
    dt <- data.table::copy(map_snp)
    dt[, cl_rank := match(CL_id, unique(CL_id))]
    dt[, col := pal_cluster[(cl_rank - 1) %% length(pal_cluster) + 1]]
    ggplot(dt, aes(Pos / 1e6, .data[[ld_w_col]], color = col)) +
      geom_point(size = 1.2, alpha = 0.85) +
      scale_color_identity() +
      theme_bw(base_size = 13) +
      labs(x = paste(chr, "position (Mbp)"), y = expression(ld["w,"*rho*"=0.95"]), title = title)
  }

  ## "high"/"low" pick the cluster set (Stage 1) and the matching group_id
  ## prefix (Combined) together, so the two panels always stay in sync --
  ## see @param direction above. NOTE: "low" must be the COMPLEMENT of the
  ## "high" (any member > threshold) cluster set, not a separate "any member
  ## < threshold" condition -- those are not the same thing, since a flagged
  ## cluster can (and often does) contain sub-threshold "boundary" members
  ## too. Using "any member < threshold" directly would pull flagged
  ## clusters right back into the "low" view.
  group_prefix <- if (direction == "high") "F" else "U"
  flagged_ids <- pruned_stage1$map_snp[get(ld_w_col) > ld_w_threshold, unique(CL_id)]
  chr_ids <- pruned_stage1$map_snp[Chr == chr, unique(CL_id)]
  needs_merge_ids <- if (direction == "high") flagged_ids else setdiff(chr_ids, flagged_ids)
  stage1_snp <- pruned_stage1$map_snp[Chr == chr & CL_id %in% needs_merge_ids]
  message(chr, " Stage 1 (", direction, "): ", uniqueN(stage1_snp$CL_id), " clusters")
  p_stage1 <- plot_clusters(
    stage1_snp,
    sprintf("Stage 1: ld_complexity_reduction() -- %d clusters (%s)", uniqueN(stage1_snp$CL_id), direction)
  )

  ## combined result restricted to this chromosome's matching-direction
  ## groups only (group_id prefix "F" for flagged/high, "U" for
  ## unflagged/low) -- matches the Stage 1 scope above; result$groups
  ## otherwise includes both buckets genome-wide, making the comparison
  ## apples-to-oranges
  groups_chr <- result$groups[Chr == chr & startsWith(group_id, group_prefix)]
  final_snp <- groups_chr[, .(marker = unlist(members)), by = group_id]
  setnames(final_snp, "group_id", "CL_id")
  final_snp <- map[Chr == chr, c("marker", "Pos", ld_w_col), with = FALSE][final_snp, on = "marker"]
  message(chr, " combined (", direction, "): ", uniqueN(final_snp$CL_id), " groups")
  p_combined <- plot_clusters(
    final_snp,
    sprintf("ld_prune_and_eMLG() -- %d groups (%s)", uniqueN(final_snp$CL_id), direction)
  )

  p_compare <- p_stage1 / p_combined
  fname <- paste0(out_folder, chr, "_stage1_vs_combined_", direction, ".png")
  ggsave(fname, p_compare, width = width, height = height)
  message("Saved: ", fname)

  invisible(p_compare)
}

plot_pruning_comparison("Chr26", pruned_stage1, result, map_hyb_005, ld_w_threshold = ld_w_threshold)
plot_pruning_comparison("Chr26", pruned_stage1, result, map_hyb_005, ld_w_threshold = ld_w_threshold, direction = "low")
