##################################################
## (expected) Multi-locus genotype complexity reduction
##################################################
## Finds clusters of highly correlated loci from population genomic data and extracts
## multi-locus-genotypes from them
##
## Needs a genotype matrix (GTs), a map file with locus info and sample_info
## file with data for sampled individuals
##
## LD-decay for empiricla data is estimated by function LD_decay in file ./R/LD_decay.R
## We will read in previously prepared data for all these files
##
## Some functions come from the LDscaR package must be installed from github
## Functions for LD-complexity reduction must be sourced from "./R/Functions_MLG.R"

## load libraries and source functions ------------------------------
library(ggplot2)
library(data.table)
library(parallel)
library(patchwork)
library(LDscnR)
source("./R_04062026/Functions_MLG.R")
source("./R_04062026/Functions_misc.R")

## Read in example data ------------------------------------

data <- readRDS("./example_data/data_for_sorting.rds")
GTs_1mb <- data$GTs
sample_info <- data$sample_info
map_1mb <- data$map
GTs_1mb_hyb <- GTs_1mb[sample_info[,which(Species=="hybrid" & Population != "Sielva")],]
rm(data)
gc()

ld_decay <- readRDS("./example_data/ld_decay_1mb_100w_hybrids_wo_sielva.rds")
gds_1mb <- create_gds_from_geno(geno=GTs_1mb, map_1mb, "gds_formicia")

# estimate mean LD in window sizes corresponding to the distance at which
# LD has decayed to 95% towards background LD rate.


map_1mb[,ld_w_095:=compute_ld_w(rho=0.95,ld_decay)]


## LD complexity reduction ------------------------------------

## first step takes LD_decay object and edge lists of pairwise r^2 values from el_folder
## to use fast single linkage clustering to find correlated sets of loci along
## chromosomes. For this LD can be estimated in sliding windows

map_CL_1mb_09_5e5  <-  LD_clustering(ld_decay, # decay object from ./R/LD_decay.R
                            map=map_1mb,
                            ld_th = 0.9, # minimum LD for single linkage clustering
                            d_th=5e5, # cannot be joined beyond this distance threshold
                            l_min=5, # minmum number of loci required for clusters
                            cores = 4)

#alternatives
map_CL_1mb <- map_CL_1mb_06_5e6
map_CL_1mb <- map_CL_1mb_06_5e6
map_CL_1mb <- map_CL_1mb_08_5e5
map_CL_1mb <- map_CL_1mb_08_1e5
map_CL_1mb <- map_CL_1mb_08_1e4
map_CL_1mb <- map_CL_1mb_08_5e4
map_CL_1mb <- map_CL_1mb_09_5e5

## check results
ggplot(map_CL_1mb, aes(Pos/1e6,ld_w_095,col=CL_col)) +
  geom_point(data=map_CL_1mb[CL_col!="grey"],size=0.25,alpha=0.5)+
  scale_color_identity() +
  facet_wrap(~ Chr, scales = "free", ncol = 5) +
  labs(
    x = "Chromosome position (Mbp)",
    y = expression("Decay rate (a) | "*ld["w,"*rho*"=0.95"]/100)
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

## extract the expected multi-locus genotypes
eMLGs_1mb_09_5e5 <- make_eMLGs(GTs=GTs_1mb_hyb,
                    map_cl = map_CL_1mb,
                    l_min = 5,
                    cores = 10)

## We are interested in sets of highly correlated SNPs that can easily be
## collpased into multi locus genotypes.

CeMLGs_1mb_09_5e5  <- collapse_eMLGs_by_chr_window(MLGs=eMLGs_1mb_09_5e5,
                             distance_threshold = 5e5,
                             r2_threshold = 0.8,
                             method = "complete",
                             prefix = "C_eMLG")

CeMLGs_1mb_09_5e5$map_eMLG[,hist(r2_eMLG)]
## join data
CeMLGs <- CeMLGs_1mb_09_5e5

#saveRDS(CeMLGs_1mb_09_5e5,"./example_data/CeMLGs_1mb_09_5e5.rds")


map_CeMLG <- CeMLGs$map_SNP[
  marker %in% map_1mb$marker,
  .(marker, C_eMLG_id)
][map_1mb, on = "marker"]

map_CeMLG <- CeMLGs$map_eMLG[,
  .(C_eMLG_id, n_loci, mean_r2_to_members, min_r2_to_members, r2_eMLG)
][map_CeMLG, on = "C_eMLG_id"]

ids <- sort(unique(map_CeMLG$C_eMLG_id))

col_map <- setNames(
  rep(col_vector, length.out = length(ids)),
  ids
)

# color for each LD-cluster
map_CeMLG[, color := col_map[as.character(C_eMLG_id)]]

map_CeMLG[is.na(color), color := "grey"]

map_CeMLG[,Chr := factor(Chr, levels = paste0("Chr",1:27))]

## Supplementary Figure of LD-complexity reduction --------------------------
p_all_chr <- ggplot(map_CeMLG, aes(Pos/1e6,ld_w_095,col=color)) +
  geom_point(data=map_CeMLG[color!="grey" ],size=0.5,alpha=0.5)+
  scale_color_identity() +
  facet_wrap(~ Chr, scales = "free", ncol = 3) +
  labs(
    x = "Chromosome position (Mbp)",
    y = expression("Decay rate (a) | "*ld["w,"*rho*"=0.95"]/100)
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(plot = p_all_chr,filename  = "./Figures/LD_complexity_reduction.png",height = 16,width = 12, units = "in")

## focusing on Chr 10
p1 <- ggplot(map_CeMLG[color!="grey" & Chr=="Chr10"], aes(Pos/1e6,ld_w_095,col=color)) +
  #geom_point(data=map_CeMLG[color=="grey" & DiagnosticIndex>-15],size=0.1,alpha=0.5)+
  geom_point(data=map_CeMLG[color!="grey" & Chr=="Chr10"])+
  scale_color_identity() +
  labs(
    x = "Chromosome position (Mbp)",
    y = expression(ld["w,"*rho*"=0.95"])
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("(A) LD network complexity reduction (all loci)")

p2 <- ggplot(map_CeMLG[color!="grey" & Chr=="Chr10"], aes(Pos/1e6,ld_w_095,col=color)) +
  geom_point(data=map_CeMLG[color!="grey" & Chr=="Chr10" & DiagnosticIndex>-15])+
  scale_color_identity() +
  labs(
    x = "Chromosome position (Mbp)",
    y = expression(ld["w,"*rho*"=0.95"])
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    pa