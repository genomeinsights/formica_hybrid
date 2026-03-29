library(data.table)
library(wesanderson)
library(circlize)
library(ggplot2)
library(SNPRelate)
library(LDscnR)

data <- readRDS("./data/filtered_DIEM_1mb.rds")
map <- data$map
GTs <- data$GTs
sample_info <- data$sample_info
rm(data)
gc()

#snpgdsClose(gds); unlink(gds)
gds <- create_gds_from_geno(geno=GTs, map, "gds_formicia")

rm(GTs) #not needed anymores
gc()map
## preliminary estimate of LD-decay based on a subset of the data and large sliding window for LD
ld_decay_pre <- compute_LD_decay(
  gds,
  ## for LD-decay and bg
  q = 0.95,
  ## for bg
  n_sub_bg = 5000,
  ## for decay
  n_win_decay = 20,
  max_pairs = 5000,
  max_SNPs_decay = 5000,
  n_strata = 20,
  overlap = 0.5,
  prob_robust = 0.95,
  keep_el = FALSE,
  slide=1000,
  cores = 10
)

if(FALSE) ld_decay_pre ## slide~400 is needed
rm(ld_decay_pre) ## don't need it anymore
gc()
## Now do the same with all sampled but a sliding window that covers 99% of LD-decay
## r2 are saved to folder "./EL/"
ld_decay <- compute_LD_decay(
  gds,
  el_data_folder = "./EL/",
  ## for LD-decay and bg
  q = 0.95,
  ## for bg
  n_sub_bg = 5000,
  ## for decay
  n_win_decay = 20,
  max_pairs = 5000,
  max_SNPs_decay = Inf,
  n_strata = 20,
  overlap = 0.5,
  prob_robust = 0.95,
  keep_el = FALSE,
  slide=400, ## slide~400 is needed
  cores = 10
)
saveRDS(ld_decay,"./data/ld_decay.rds")
gc()

## check that ld_w works fine
if(FALSE){
  ld_w_99 <- compute_ld_w(ld_decay,rho = 0.99,cores = 10)
  plot(ld_w_99)
}

q("no")
