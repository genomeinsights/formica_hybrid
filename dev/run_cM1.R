library(ggplot2)
library(igraph)
library(data.table)
library(SNPRelate)
library(parallel)
devtools::load_all("~/gitlab/LDscnR/")

load("./data/hybrids_only_maf005.Rdata")
rec_map <- fread("./data/Frufa_DTOL_PR.ref_genome.recmap")
rec_map[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
genetic_map <- rec_map[, .(Chr, Pos = pos, cM)]

# pruned_stage1 <- ld_complexity_reduction(
#   map = map_hyb_005, LD_decay = ld_decay, rho = 0.5, cores = 1
# )

#saveRDS(pruned_stage1,"data/pruned_stage1.rds")

pruned_stage1 <- readRDS("data/pruned_stage1.rds")

eMLG_5loci_0025_cM1 <- ld_prune_and_eMLG(
  GTs = GTs_hybrids_005, stage1 = pruned_stage1, ld_w_col = "ld_w_095",
  ld_w_threshold = 0.025, score_threshold = 0.80, min_r2 = 0.2, min_n_loci_flag = 5,
  genetic_map = genetic_map, cM_threshold = 1,
  compute_unflagged_eMLG = TRUE, min_n_loci_eMLG = 5
)

saveRDS(eMLG_5loci_0025_cM1,"./data/eMLG_5loci_0025_cM1.rds")
eMLG_5loci_0025_cM1 <- readRDS("./data/eMLG_5loci_0025_cM1.rds")
# 
plot_pruning_comparison("Chr10", pruned_stage1, eMLG_5loci_0025_cM1, map_hyb_005)
plot_pruning_comparison("Chr10", pruned_stage1, eMLG_5loci_0025_cM1, map_hyb_005, direction = "low")
# 
# 
# eMLG_5loci_0025_cM1 <- readRDS("./data/eMLG_5loci_0025_cM1.rds")
# eMLG_5loci_0025_cM05 <- readRDS("./data/eMLG_5loci_0025_cM05.rds")
# plot_pruning_comparison("Chr10", pruned_stage1, eMLG_5loci_0025_cM1, map_hyb_005)
# plot_pruning_comparison("Chr10", pruned_stage1, eMLG_5loci_0025_cM1, map_hyb_005, direction = "low")
# 
# 

## data for Bea
# groups <- eMLG_5loci_0025_cM1$groups
# 
# g <- groups[, .(group_id, members)]
# 
# long <- data.table(
#   chromosome = gsub("Chr","chromosome_",rep(groups$Chr, lengths(groups$members))),
#   marker = unlist(groups$members),
#   group_id = rep(groups$group_id, lengths(groups$members))
# )
# long <- map[,.(marker,DiagnosticIndex,position=Pos)][long,on="marker"]
# long[,LG:=as.numeric(gsub("chromosome_","",chromosome))]
# setorder(long,LG,position)
# long <- long[,.(chromosome,position,DiagnosticIndex,group_id)]
# 
# AIMs <-  fread("data/species_diagnostic_marker_fixation_levels.tsv")
# 
# table(AIMs[,paste(chromosome,position)] %in% long[,paste(chromosome,position)])
# 
# saveRDS(long[DiagnosticIndex> -25],"data/cluster_info.rds")
