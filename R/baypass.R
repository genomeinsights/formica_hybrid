library(ggplot2)
library(data.table)
library(LDscnR)
library(SNPRelate)
source("./R/prune_snps_for_grm.R")
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

message("=== Reading in data ===")
path_to_baypass = "/Users/petrikemppainen/baypass_public-master/sources/g_baypass" ## wherever your baypass is located
baypass_folder <- "./out_baypass/"


data <- readRDS("./data/filtered_DIEM_1mb.rds")
map <- data$map
GTs <- data$GTs
sample_info <- data$sample_info
rm(data)

pop <- sample_info$Population

ld_decay <- readRDS("./data/ld_decay_corr.rds")
## get GRM

message("=== Creating gds ======")
gds <- create_gds_from_geno(geno=GTs, map, "gds_formicia")


message("=== Estimating ld_w ===")
map[,ld_w_95:=compute_ld_w(ld_decay,rho = 0.95)]
map[is.na(ld_w_95),ld_w_95:=0]

message("=== Pruning SNPs ======")
pruned_SNPs <- prune_snps_for_grm(map,
                   ld_decay,
                   rho_ld = 0.5,
                   ld_w_col = "ld_w_95",
                   edge_ld_col = "r2",
                   pos1_col = "pos1",
                   pos2_col = "pos2",
                   marker_sep = ":",
                   block_prefix = NULL,
                   show_progress = FALSE)

message("=== Estimating Omega ===")
baypass_pruned <- do.call(cbind, lapply(unique(pop), function(y){
  t(apply(GTs[pop==y,map$marker %in% pruned_SNPs$core_snp], 2, function(x){
    c(length(which(x==0))*2+length(which(x==1)),
      length(which(x==2))*2+length(which(x==1))
    )
  }))
}))


write.table(baypass_pruned, file=paste0(baypass_folder,"u_DIEM.geno_pruned"), quote = FALSE, row.names = FALSE, col.names = FALSE)
poolsize <- t(as.vector(table(pop)))
write.table(poolsize, file=paste0(baypass_folder,"u_DIEM.size"), quote = FALSE, row.names = FALSE, col.names = FALSE)

system(paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno_pruned -nthreads 12 -nval 500 -burnin 5000 -thin 10 -poolsizefile ",baypass_folder,"u_DIEM.size -outprefix ", baypass_folder, "omega"))

Omega <- as.matrix(read.table("omega_omega.out"))

baypass <- do.call(cbind, lapply(unique(pop), function(y){
  t(apply(GTs[pop==y,], 2, function(x){
    c(length(which(x==0))*2+length(which(x==1)),
      length(which(x==2))*2+length(which(x==1))
    )
  }))
}))

write.table(baypass, file=paste0(baypass_folder,"u_DIEM.geno"), quote = FALSE, row.names = FALSE, col.names = FALSE)

pc1_env <- sample_info[!duplicated(Population),PC1]
write.table(t(as.matrix(pc1_env)), file=paste0(baypass_folder,"u.PC1"), quote = FALSE, row.names = FALSE, col.names = FALSE)
paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno -efile ", baypass_folder, "u.PC1 -poolsizefile ",baypass_folder,"u_DIEM.size -nthreads ",cores," -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ", baypass_folder, "PC1_DIEM")

pc2_env <- sample_info[!duplicated(Population),PC2]
write.table(t(as.matrix(pc2_env)), file=paste0(baypass_folder,"u.PC2"), quote = FALSE, row.names = FALSE, col.names = FALSE)
paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno -efile ", baypass_folder, "u.PC2 -poolsizefile ",baypass_folder,"u_DIEM.size -nthreads ",cores," -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ", baypass_folder, "PC2_DIEM")


message("=== Writing baypass ===")


#rep <- 1
cores=12
message("=== Running baypass ===")
for(rep in 1:10){
  ## get simulated phenotype
  phe_sim <- simulate_null_continuous(GRM,h2 = 0.5)
  write.table(t(as.matrix(phe_sim)), file=paste0(baypass_folder,rep,".env"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  ## running Baypass
  system(paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno -efile ", paste0(baypass_folder,rep,".env"), " -poolsizefile ",baypass_folder,"u_DIEM.size -nthreads ",cores," -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ", paste0(baypass_folder,rep)))
}
q("no")
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
