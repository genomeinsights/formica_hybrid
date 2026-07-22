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

pop <- sample_data$Population

gds_hyb <- create_gds_from_geno(geno=GTs_hybrids_005, map_hyb_005, "gds_hybrids.gds")

# ------------------------------------------------------------
## LD complexity reduction 
# ------------------------------------------------------------

## first step takes LD_decay object and edge lists of pairwise r^2 values from el_folder
## to use fast single linkage clustering to find correlated sets of loci along
## chromosomes. For this LD can be estimated in sliding windows

map_CL_1mb_09_5e5  <-  LD_clustering(ld_decay, # decay object from ./R/LD_decay.R
                                     map=map_1mb,
                                     ld_th = 0.9, # minimum LD for single linkage clustering
                                     d_th=5e5, # cannot be joined beyond this distance threshold
                                     l_min=5, # minmum number of loci required for clusters
                                     cores = 4)


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