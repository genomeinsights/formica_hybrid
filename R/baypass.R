library(ggplot2)
library(igraph)
library(data.table)
library(SNPRelate)
library(parallel)
devtools::load_all("~/gitlab/LDscnR/")
# ------------------------------------------------------------
# Permutation functions not currently used
# ------------------------------------------------------------

simulate_pop_null_continuous <- function(Omega, mu = NULL, sigma2 = 1, scale_y = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  P <- nrow(Omega)
  stopifnot(ncol(Omega) == P)

  if (is.null(mu)) mu <- rep(0, P)
  stopifnot(length(mu) == P)

  Omega <- (Omega + t(Omega)) / 2
  eig <- eigen(Omega, symmetric = TRUE)
  vals <- pmax(eig$values, 0)

  z <- rnorm(P)
  y <- as.numeric(mu + eig$vectors %*% (sqrt(sigma2 * vals) * z))

  if (scale_y) y <- as.numeric(scale(y))
  y
}

simulate_pop_null_binary <- function(Omega, prevalence, mu = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  P <- nrow(Omega)
  stopifnot(ncol(Omega) == P)
  stopifnot(prevalence > 0 && prevalence < 1)

  if (is.null(mu)) mu <- rep(0, P)
  stopifnot(length(mu) == P)

  Omega <- (Omega + t(Omega)) / 2
  eig <- eigen(Omega, symmetric = TRUE)
  vals <- pmax(eig$values, 0)

  liability <- as.numeric(mu + eig$vectors %*% (sqrt(vals) * rnorm(P)))
  thr <- as.numeric(quantile(liability, probs = 1 - prevalence, type = 8))
  y <- as.integer(liability > thr)

  list(
    y = y,
    liability = liability,
    threshold = thr
  )
}

simulate_null_continuous <- function(K, h2 = 0.5, scale_y = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(K)
  stopifnot(ncol(K) == n)

  # Numerical stabilization
  K <- (K + t(K)) / 2
  eig <- eigen(K, symmetric = TRUE)
  vals <- pmax(eig$values, 0)

  # g ~ N(0, h2 K)
  z_g <- rnorm(n)
  g <- eig$vectors %*% (sqrt(h2 * vals) * z_g)

  # e ~ N(0, (1-h2) I)
  e <- rnorm(n, mean = 0, sd = sqrt(1 - h2))

  y <- as.numeric(g + e)

  if (scale_y) {
    y <- as.numeric(scale(y))
  }

  y
}


# ------------------------------------------------------------
# Read in data (generated in LD_decay_from_DIEM.R)
# ------------------------------------------------------------

message("=== Loading data and Creating gds ======")
# loads GTs_hybrids_005,map_hyb_005, ld_decay and sample_data from LD_decay_from_DIEM.R
# includes maf and ld_w_095
load("./data/hybrids_only_maf005.Rdata")

## TEMPORARY safeguard, not a real fix: GTs_hybrids_005 contains 1,811 of
## 3,557,980 non-NA values (0.05%) outside {0,1,2} -- specifically -5 (517),
## -2 (105), -1 (206), 3 (308), 4 (85), 7 (590), presumably DIEM-specific
## sentinel/state codes from the parsing step in LD_decay_from_DIEM.R, not
## valid genotype dosages. Left untreated, these silently corrupt anything
## downstream that assumes hard 0/1/2 calls (expected_gt_dosage()'s
## consensus signals ended up as high as 7 instead of the valid [0,2]
## range). What these codes actually mean (missing calls? valid states
## needing remapping instead of NA?) is still an open question -- this
## line is a stopgap (treat them as missing) so downstream eMLG work isn't
## silently wrong in the meantime; the real fix belongs in
## LD_decay_from_DIEM.R's DIEM-parsing step once that's resolved.
GTs_hybrids_005[!(GTs_hybrids_005 %in% c(0, 1, 2))] <- NA

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

## Classify CLUSTERS, not markers: a cluster needs the expensive
## all-pairwise merge step if ANY of its members individually reads above
## the ld_w threshold -- i.e. sits in a denser-than-average local LD
## neighbourhood, which is where the sliding window is most likely to have
## further fragmented the true block into more than one Stage-1 cluster
## (see ld_complexity_reduction()'s "Known limitation", and
## merge_ld_clusters()'s documented example for this exact pattern). This
## reuses Stage 1's already-correct cluster boundaries instead of
## re-deriving them from a threshold, so no real block gets split by the
## classification itself -- only ~6% of markers genome-wide exceed 0.2, so
## the flagged subset stays small and merge_ld_clusters() stays cheap
## (even the largest single-chromosome high-ld_w subset, Chr10 at 8,875
## markers, finishes both stages in ~7s).
## Filtering map_snp directly (one vectorized column comparison) instead
## of vapply-ing over every cluster's members is ~1,600x faster in
## practice (3.15s -> 0.002s on a 15,524-cluster chromosome, checked
## directly) -- map_snp already has ld_w_095 and CL_id on the same
## per-marker row, so no per-cluster R-level loop or named-vector lookup
## is needed at all.
##
## Matching on CL_id with %in% (not positional/negative indexing like
## clusters[needs_merge] / clusters[-needs_merge]) is deliberate: CL_id is
## a value (ld_complexity_reduction()'s .GRP, assigned once across the
## whole call), and while it currently happens to equal row position in
## this specific (fresh, unfiltered) clusters table, that's an
## implementation detail, not a documented guarantee -- relying on it is
## exactly the kind of silent, no-error-thrown bug already hit once this
## session (merge_ld_clusters()'s R2[sub$CL_id,...] indexing). %in% costs
## nothing extra here and doesn't depend on that coincidence holding.
needs_merge_ids <- pruned_stage1$map_snp[ld_w_095 > ld_w_threshold, unique(CL_id)]

flagged   <- pruned_stage1$clusters[CL_id %in% needs_merge_ids]
unflagged <- pruned_stage1$clusters[!CL_id %in% needs_merge_ids]

ld_result_flagged <- list(
  map_snp  = pruned_stage1$map_snp[marker %in% unlist(flagged$members)],
  clusters = flagged,
  pruned   = flagged$core_snp
)

pruned_merged <- merge_ld_clusters(
  GTs = GTs_hybrids_005[TRUE, unlist(flagged$members)],
  ld_result = ld_result_flagged, LD_decay = ld_decay, rho = 0.5, cores = 1
)

pruned_markers <- c(unflagged$core_snp, pruned_merged$pruned)

message(
  "Keeping ", length(pruned_markers), "/", map_hyb_005[TRUE,.N], " (",
  round(100 * length(pruned_markers) / map_hyb_005[TRUE,.N], 2), "%) SNPs"
)
saveRDS(pruned_markers,"./data/pruned_markers.rds")
# ------------------------------------------------------------
# Chr26 diagnostic: visualise clusters after each pruning stage
# (Stage 1: ld_complexity_reduction, Stage 2: merge_ld_clusters), for the
# ld_w>0.2 markers -- where fragmented centromeric/inversion blocks are
# most likely. This is the same split used above, worked through for one
# chromosome so the effect of merging is visible.
# ------------------------------------------------------------
library(patchwork)

## cycle a modest, visually-distinct palette across CL_id -- with hundreds
## of clusters no palette gives every one a unique color, but neighbouring
## clusters (what we care about here) will very likely differ
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

## Pulled straight from the whole-genome objects the main pipeline already
## computed above (pruned_stage1, flagged, pruned_merged) instead of an
## independent idx-filtered recomputation. The previous version called
## ld_complexity_reduction(idx = idx_chr26) on markers pre-filtered to
## ld_w_095 > ld_w_threshold -- exactly the "filter markers before
## clustering" approach shown above to sever real blocks at the threshold
## boundary (226 of 15,524 Chr26 clusters straddled it in a combined run,
## cutting 975 low-ld_w markers loose). Reusing the flagged-cluster
## results here means the plots show what the pipeline actually does, and
## costs nothing extra since Stage 1/2 aren't recomputed.
chr26_stage1_snp <- pruned_stage1$map_snp[Chr == "Chr26" & CL_id %in% flagged$CL_id]

message("=== Chr26 diagnostic, Stage 1: ld_complexity_reduction (flagged clusters only) ===")
message("Stage 1: ", uniqueN(chr26_stage1_snp$CL_id), " clusters")

p_chr26_stage1 <- plot_clusters_chr26(
  chr26_stage1_snp,
  sprintf("Stage 1: ld_complexity_reduction() -- %d clusters", uniqueN(chr26_stage1_snp$CL_id))
)
ggsave("./Figures/chr26_stage1_ld_complexity_reduction.png", p_chr26_stage1, width = 10, height = 5)

## every flagged cluster's full marker set (including any lower-ld_w
## members merged in) is already in pruned_merged$map_snp -- no unflagged
## cluster can contain a Chr26 marker above ld_w_threshold by construction
chr26_stage2_snp <- pruned_merged$map_snp[Chr == "Chr26"]

message("=== Chr26 diagnostic, Stage 2: merge_ld_clusters (flagged clusters only) ===")
message("Stage 2: ", uniqueN(chr26_stage2_snp$CL_id), " clusters")

p_chr26_stage2 <- plot_clusters_chr26(
  chr26_stage2_snp,
  sprintf("Stage 2: merge_ld_clusters() -- %d clusters", uniqueN(chr26_stage2_snp$CL_id))
)
ggsave("./Figures/chr26_stage2_merge_ld_clusters.png", p_chr26_stage2, width = 10, height = 5)

p_chr26_compare <- p_chr26_stage1 / p_chr26_stage2
ggsave("./Figures/chr26_stage1_vs_stage2.png", p_chr26_compare, width = 10, height = 9)
p_chr26_compare

# ------------------------------------------------------------
# Chr26 diagnostic: eMLG extraction and correlation heatmap, plus the
# experimental dynamic tree cut (dev/R/dynamic_cut_eMLG.R) that further
# consolidates Stage-2 clusters into coherent eMLG groups -- average
# linkage + a score_eMLG/pair_r2 quality gate walked bottom-up per merge,
# rather than a single global r2 threshold. See that file for the
# validated rationale (complete linkage fragments diffuse-but-real blocks;
# single linkage chains catastrophically, even just as the tree structure
# feeding the quality gate).
# ------------------------------------------------------------
message("=== Chr26 diagnostic, eMLG extraction and correlation heatmap ===")

chr26_clusters <- pruned_merged$clusters[Chr == "Chr26"]

chr26_eMLG <- do.call(cbind, lapply(chr26_clusters$members, function(mk) {
  expected_gt_dosage(GTs_hybrids_005[, mk, drop = FALSE])
}))
colnames(chr26_eMLG) <- chr26_clusters$CL_id
message("Chr26 eMLG matrix: ", nrow(chr26_eMLG), " individuals x ", ncol(chr26_eMLG), " clusters")

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
    labs(x = "Cluster", y = "Cluster", title = title)
}

## position-ordered heatmap -- shows raw structure before any grouping
pos_lookup <- setNames(map_hyb_005$Pos, map_hyb_005$marker)
cluster_pos <- vapply(chr26_clusters$members, function(mk) mean(pos_lookup[mk]), numeric(1))
ord_pos <- order(cluster_pos)

p_eMLG_pos <- plot_eMLG_heatmap(
  R2_chr26[ord_pos, ord_pos],
  sprintf("Chr26 eMLG correlation -- %d clusters, position order", ncol(R2_chr26))
)
ggsave("./Figures/chr26_eMLG_heatmap_posorder.png", p_eMLG_pos, width = 8, height = 7, dpi = 150)

message("=== Chr26 diagnostic, dynamic tree cut (experimental) ===")
source("./dev/R/dynamic_cut_eMLG.R")
n_loci_chr26 <- setNames(chr26_clusters$n_snps, chr26_clusters$CL_id)
dyn_score_threshold <- 0.80

dyn_groups <- dynamic_cut_eMLG(chr26_eMLG, n_loci_chr26, threshold = dyn_score_threshold, min_r2 = 0.2)
message("Dynamic cut (score_eMLG>=", dyn_score_threshold, "): ",
        ncol(chr26_eMLG), " clusters -> ", length(dyn_groups), " groups")

## heatmap reordered by dynamic-cut group assignment -- should show the
## soft/diffuse blocks the position-ordered view only hints at as one
## clean group each, unlike the Stage-2 (complete-linkage) result above
grp_of <- setNames(
  rep(seq_along(dyn_groups), vapply(dyn_groups, function(g) length(g$members), integer(1))),
  unlist(lapply(dyn_groups, `[[`, "members"))
)
ord_dyn <- order(grp_of[colnames(chr26_eMLG)])

p_eMLG_dyn <- plot_eMLG_heatmap(
  R2_chr26[ord_dyn, ord_dyn],
  sprintf("Chr26 eMLG correlation -- reordered by dynamic cut (%d groups)", length(dyn_groups))
)
ggsave("./Figures/chr26_eMLG_heatmap_dynamiccut.png", p_eMLG_dyn, width = 8, height = 7, dpi = 150)

## same ld_w-vs-position view as the Stage 1/Stage 2 plots above, but
## colored by dynamic-cut group instead of Stage-2 CL_id -- expand each
## group's Stage-2 members back down to raw markers via chr26_clusters
chr26_clusters[, dyn_group := grp_of[as.character(CL_id)]]
chr26_dyn_snp <- chr26_clusters[, .(marker = unlist(members)), by = dyn_group]
setnames(chr26_dyn_snp, "dyn_group", "CL_id")
chr26_dyn_snp <- map_hyb_005[Chr == "Chr26", .(marker, Pos, ld_w_095)][chr26_dyn_snp, on = "marker"]

p_chr26_dyn <- plot_clusters_chr26(
  chr26_dyn_snp,
  sprintf("Dynamic cut (score_eMLG>=%.2f) -- %d groups", dyn_score_threshold, length(dyn_groups))
)
ggsave("./Figures/chr26_dynamiccut_ldw.png", p_chr26_dyn, width = 10, height = 5)

message("=== Estimating Omega ===")
path_to_baypass = "~/gitlab/baypass_public-master/sources/g_baypass" ## wherever your baypass is located
baypass_folder <- "./out_baypass/"
pop <- sample_data$Population
baypass_pruned <- do.call(cbind, lapply(unique(pop), function(y){
  t(apply(GTs_hybrids_005[pop==y,map_hyb_005$marker %in% pruned_markers], 2, function(x){
    c(length(which(x==0))*2+length(which(x==1)),
      length(which(x==2))*2+length(which(x==1))
    )
  }))
}))


write.table(baypass_pruned, file=paste0(baypass_folder,"u_DIEM.geno_pruned"), quote = FALSE, row.names = FALSE, col.names = FALSE)
poolsize <- t(as.vector(table(pop)))
write.table(poolsize, file=paste0(baypass_folder,"u_DIEM.size"), quote = FALSE, row.names = FALSE, col.names = FALSE)
cores=10
#paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno_pruned -nthreads ",cores," -nval 500 -burnin 5000 -thin 10 -poolsizefile ",baypass_folder,"u_DIEM.size -outprefix ", baypass_folder, "omega")

Omega <- as.matrix(read.table("./out_baypass/omega_mat_omega.out"))

baypass <- do.call(cbind, lapply(unique(pop), function(y){
  t(apply(GTs_hybrids_005[pop==y,], 2, function(x){
    c(length(which(x==0))*2+length(which(x==1)),
      length(which(x==2))*2+length(which(x==1))
    )
  }))
}))

write.table(baypass, file=paste0(baypass_folder,"u_DIEM.geno"), quote = FALSE, row.names = FALSE, col.names = FALSE)

pc1_env <- sample_data[!duplicated(Population),PC1]
write.table(t(as.matrix(pc1_env)), file=paste0(baypass_folder,"u.PC1"), quote = FALSE, row.names = FALSE, col.names = FALSE)

pc2_env <- sample_data[!duplicated(Population),PC2]
write.table(t(as.matrix(pc2_env)), file=paste0(baypass_folder,"u.PC2"), quote = FALSE, row.names = FALSE, col.names = FALSE)

## run from command line preferrably
## no pruning
paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno -efile ", baypass_folder, "u.PC1 -poolsizefile ",baypass_folder,"u_DIEM.size -nthreads ",cores," -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ", baypass_folder, "PC1_DIEM")
paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno -efile ", baypass_folder, "u.PC2 -poolsizefile ",baypass_folder,"u_DIEM.size -nthreads ",cores," -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ", baypass_folder, "PC2_DIEM")


## wit omega on pruned data
paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno -omegafile ./out_baypass/omega_mat_omega.out -efile ", baypass_folder, "u.PC2 -poolsizefile ",baypass_folder,"u_DIEM.size -nthreads ",cores," -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ", baypass_folder, "PC2_DIEM")
paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno -omegafile ./out_baypass/omega_mat_omega.out -efile ", baypass_folder, "u.PC1 -poolsizefile ",baypass_folder,"u_DIEM.size -nthreads ",cores," -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ", baypass_folder, "PC1_DIEM")

#message("=== Writing baypass ===")


#rep <- 1
# cores=12
# message("=== Running baypass ===")
# for(rep in 1:10){
#   ## get simulated phenotype
#   phe_sim <- simulate_null_continuous(GRM,h2 = 0.5)
#   write.table(t(as.matrix(phe_sim)), file=paste0(baypass_folder,rep,".env"), quote = FALSE, row.names = FALSE, col.names = FALSE)
# 
#   ## running Baypass
#   system(paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno -efile ", paste0(baypass_folder,rep,".env"), " -poolsizefile ",baypass_folder,"u_DIEM.size -nthreads ",cores," -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ", paste0(baypass_folder,rep)))
# }
# q("no")
#
# cores = 12
# pc1_env <- sample_info[!duplicated(Population),PC1]
# write.table(t(as.matrix(pc1_env)), file=paste0(baypass_folder,"u.PC1"), quote = FALSE, row.names = FALSE, col.names = FALSE)
# paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno -efile ", baypass_folder, "u.PC1 -poolsizefile ",baypass_folder,"u_DIEM.size -nthreads ",cores," -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ", baypass_folder, "PC1_DIEM")
#
# pc2_env <- sample_info[!duplicated(Population),PC2]
# write.table(t(as.matrix(pc2_env)), file=paste0(baypass_folder,"u.PC2"), quote = FALSE, row.names = FALSE, col.names = FALSE)
# paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno -efile ", baypass_folder, "u.PC2 -poolsizefile ",baypass_folder,"u_DIEM.size -nthreads ",cores," -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ", baypass_folder, "PC2_DIEM")
#

#
#
# simulate_null_binary <- function(K, prevalence, h2 = 0.5, seed = NULL) {
#   if (!is.null(seed)) set.seed(seed)
#
#   n <- nrow(K)
#   stopifnot(ncol(K) == n)
#   stopifnot(prevalence > 0 && prevalence < 1)
#
#   # Numerical stabilization
#   K <- (K + t(K)) / 2
#   eig <- eigen(K, symmetric = TRUE)
#   vals <- pmax(eig$values, 0)
#
#   # latent liability
#   z_g <- rnorm(n)
#   g <- eig$vectors %*% (sqrt(h2 * vals) * z_g)
#   e <- rnorm(n, mean = 0, sd = sqrt(1 - h2))
#   liability <- as.numeric(g + e)
#
#   # Threshold chosen to match target prevalence
#   # prevalence = proportion with phenotype 1
#   thr <- as.numeric(quantile(liability, probs = 1 - prevalence, type = 8))
#
#   y_bin <- as.integer(liability > thr)
#
#   list(
#     y = y_bin,
#     liability = liability,
#     threshold = thr
#   )
# }
