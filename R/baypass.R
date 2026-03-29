library(ggplot2)
library(data.table)
path_to_baypass = "/Users/petrikemppainen/baypass_public-master/sources/g_baypass" ## wherever your baypass is located
baypass_folder <- "./out_baypass/"


data <- readRDS("./data/filtered_DIEM_1mb.rds")
map_DIEM <- data$map
GTs_DIEM <- data$GTs
sample_info <- data$sample_info
rm(data)
gc()

pop <- sample_info$Population

baypass <- do.call(cbind, lapply(unique(pop), function(y){
  t(apply(GTs_DIEM[pop==y,], 2, function(x){
    c(length(which(x==0))*2+length(which(x==1)),
      length(which(x==2))*2+length(which(x==1))
    )
  }))
}))

write.table(baypass, file=paste0(baypass_folder,"u_DIEM.geno"), quote = FALSE, row.names = FALSE, col.names = FALSE)
poolsize <- t(as.vector(table(pop)))
write.table(poolsize, file=paste0(baypass_folder,"u_DIEM.size"), quote = FALSE, row.names = FALSE, col.names = FALSE)
cores = 12
pc1_env <- sample_info[!duplicated(Population),PC1]
write.table(t(as.matrix(pc1_env)), file=paste0(baypass_folder,"u.PC1"), quote = FALSE, row.names = FALSE, col.names = FALSE)
paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno -efile ", baypass_folder, "u.PC1 -poolsizefile ",baypass_folder,"u_DIEM.size -nthreads ",cores," -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ", baypass_folder, "PC1_DIEM")

pc2_env <- sample_info[!duplicated(Population),PC2]
write.table(t(as.matrix(pc2_env)), file=paste0(baypass_folder,"u.PC2"), quote = FALSE, row.names = FALSE, col.names = FALSE)
paste0(path_to_baypass, " -countdatafile ", baypass_folder, "u_DIEM.geno -efile ", baypass_folder, "u.PC2 -poolsizefile ",baypass_folder,"u_DIEM.size -nthreads ",cores," -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74 -outprefix ", baypass_folder, "PC2_DIEM")
