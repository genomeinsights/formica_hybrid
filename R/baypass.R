library(ggplot2)
library(igraph)
library(data.table)
library(LDscnR)
library(SNPRelate)
library(parallel)

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

message("=== Reading in data ===")

path_to_baypass = "~/gitlab/baypass_public-master/sources/g_baypass" ## wherever your baypass is located
baypass_folder <- "./out_baypass/"

tmp <- readRDS("./data/diem_parsed.rds")
GTs <- tmp$GTs
map <- tmp$map
sample_data <- fread("data/Sample_covariate_info_outlier_analysis_20.txt")
rm(tmp)
gc()


pop <- sample_data$Population

ld_decay <- readRDS("./data/ld_decay_DIEM_100w.rds")

## get GRM

message("=== Keeping hybrids and Creating gds ======")
GTs_hybrids <- GTs[sample_data$Sample_ID,]

## filter by maf easies through SNP relate
#showfile.gds(closeall = TRUE)
gds_hyb <- create_gds_from_geno(geno=GTs_hybrids, map, "gds_hybrids.gds")
map[,maf_hyb:= snpgdsSNPRateFreq(gds_hyb)$MinorFreq]
map_hyb_005 <- map[maf_hyb>=0.05]
GTs_hybrids_005 <- GTs_hybrids[,map_hyb_005$marker]

## Redo with maf>0.05
gds_hyb <- create_gds_from_geno(geno=GTs_hybrids_005, map_hyb_005, "gds_hybrids.gds")
#showfile.gds(closeall = TRUE)
#
ld_decay_light <- compute_LD_decay(
  gds_hyb,n_win_decay = 100,
  el_data_folder = NULL, # too large to keep in memory
  keep_el = TRUE,
  slide=100,
  ld_method = "corr"
)


message("=== Estimating ld_w ===")
map_hyb_005[,ld_w_095:=compute_ld_w(ld_decay,rho = 0.95)]
map_hyb_005[is.na(ld_w_095),ld_w_095:=0]
map_hyb_005[,ld_w_05:=compute_ld_w(ld_decay,rho = 0.5)]
map_hyb_005[is.na(ld_w_05),ld_w_05:=0]

map_hyb_005[,plot(ld_w_05,ld_w_095)]

message("=== Pruning SNPs ======")
pruned_res <- ld_complexity_reduction(GTs=GTs_hybrids_005, map=map_hyb_005, LD_decay=ld_decay, rho = 0.5, cores = 1)
pruned_markers <- pruned_res$pruned
message("Keeping ",length(pruned_markers),"/",map_hyb_005[,.N], "SNPs")

message("=== Estimating Omega ===")
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
system(paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno_pruned -nthreads ",cores," -nval 500 -burnin 5000 -thin 10 -poolsizefile ",baypass_folder,"u_DIEM.size -outprefix ", baypass_folder, "omega"))

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
paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno -efile ", baypass_folder, "u.PC1 -poolsizefile ",baypass_folder,"u_DIEM.size -nthreads ",cores," -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ", baypass_folder, "PC1_DIEM")

pc2_env <- sample_data[!duplicated(Population),PC2]
write.table(t(as.matrix(pc2_env)), file=paste0(baypass_folder,"u.PC2"), quote = FALSE, row.names = FALSE, col.names = FALSE)
paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno -efile ", baypass_folder, "u.PC2 -poolsizefile ",baypass_folder,"u_DIEM.size -nthreads ",cores," -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ", baypass_folder, "PC2_DIEM")


# ~/gitlab/baypass_public-master/sources/g_baypass -countdatafile ./out_baypass/u_DIEM.geno -efile ./out_baypass/u.PC1 -poolsizefile ./out_baypass/u_DIEM.size -nthreads 1 -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ./out_baypass/PC1_DIEM
# ~/gitlab/baypass_public-master/sources/g_baypass -countdatafile ./out_baypass/u_DIEM.geno -efile ./out_baypass/u.PC2 -poolsizefile ./out_baypass/u_DIEM.size -nthreads 1 -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ./out_baypass/PC2_DIEM


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
