library(ggplot2)
library(data.table)
library(parallel)
library(patchwork)
library(LDscnR)
source("./R_04062026/Ohta.R")

CeMLGs_1mb <- readRDS("./example_data/CeMLGs_1mb_09_5e5.rds")

keep <- sorting_res_C_eMLG[,maf_parents>0.05 & parent_missing_pol<=0.5 & parent_missing_aqu<=0.5]
dt <- sorting_res_C_eMLG[keep & !is.na(C_eMLG_id),.(mean_ldw=mean(ld_w_095),mean_Bi=mean(bi_index),mean_cl_size=mean(n_loci),mean_aqu=mean(delta_pol_hybrid),mean_pol=mean(delta_aqu_hybrid),mean_DI=mean(DiagnosticIndex)),by=.(C_eMLG_id,color)]

CMLGs <- round(CeMLGs_1mb$eMLG)

keep <- sample_info[,which(Species=="hybrid" & Population != "Sielva")]
prep <- ohta_fast_prepare(data_set = CMLGs, pops=sample_info[keep,Population])

pairs <- t(combn(ncol(CMLGs),2))

table(apply(CMLGs,2,function(x) length(unique(na.omit(x)))))
# get ld first and prefilter by r2>0.01
m <- cor(CMLGs,use="pair")^2
m[diag(m)] <- NA
hist(m[m>0.1],breaks=100)
idx <- which(upper.tri(m) & !is.na(m), arr.ind = TRUE)

el <- data.table(
  idx1 = idx[,1],
  idx2 = idx[,2],
  MLG1 = colnames(m)[idx[,1]],
  MLG2 = colnames(m)[idx[,2]],
  r2 = m[idx]
)
CeMLGs_pos <- map_CeMLG[!is.na(C_eMLG_id),.(mid=trunc(median(Pos)),start=min(Pos),end=max(Pos)),by=.(Chr,C_eMLG_id)]
el <- cbind(CeMLGs_pos[match(el$MLG1,C_eMLG_id)][,.(Chr1=Chr, mid1=mid,   start1=start,     end1=end)],CeMLGs_pos[match(el$MLG2,C_eMLG_id)][,.(Chr2=Chr, mid2=mid,   start2=start,     end2=end)],el)

qqplot(el[Chr1!=Chr2,r2],el[Chr1==Chr2,r2])
abline(0,1)

el <- el[r2>0.2]

pairs <- as.matrix(el[,.(idx1, idx2)])

ohta <- dstat_unphased_scan(
  pairs = pairs,
  prep = prep,
  tot_maf = 0,
  pop_maf = 0,
  cores = 1
)
ohta[is.na(D2st)]
ohta[,Ohta_D := D2is-D2st]
pairs(ohta[,.(r2,Fst_mean,D2is,D2st,Ohta_D)])

ohta <- cbind(el,ohta)

ohta[, min_dist := fifelse(
  Chr1 != Chr2,
  NA_real_,
  pmax(0, pmax(start1, start2) - pmin(end1, end2))
)]

ohta[Chr1 == Chr2]
ohta[Chr1 != Chr2]

pairs(ohta[,.(r2,Fst_mean,D2is,D2st,Ohta_D)])

ohta[,.(r2,Fst_mean,D2is,D2st,Ohta_D)]
#ohta_D_q
ohta[,ohta_D_q := rank(Ohta_D,ties.method = "random")/.N]
ohta[,D2is_q := rank(D2is,ties.method = "random")/.N]
ohta[,D2st_q := rank(D2st,ties.method = "random")/.N]
ohta[,r2_q := rank(r2,ties.method = "random")/.N]
qt <- 0.9
out <- rbindlist(lapply(seq(0,0.95,by=0.05),function(qt){
  tmp <- unlist(unique(ohta[r2_q>qt,.(MLG1,MLG2)]))
  data.table(qt,mean_Bi = dt[C_eMLG_id %in% tmp,mean_Bi],
             mean_aqu = dt[C_eMLG_id %in% tmp,mean_aqu],
             mean_pol = dt[C_eMLG_id %in% tmp,mean_pol],
             mean_DI = dt[C_eMLG_id %in% tmp,mean_DI])

}))

ggplot(out, aes(factor(qt),mean_Bi))+
  geom_boxplot()

ggplot(out, aes(factor(qt),mean_pol))+
  geom_boxplot()

ggplot(out, aes(factor(qt),mean_aqu))+
  geom_boxplot()

ggplot(out, aes(factor(qt),mean_DI))+
  geom_boxplot()



# ohta[,hist(Ohta_D)]
# ohta[,plot(D2is_q,D2st_q)]
keep <- dt[mean_DI>-15,C_eMLG_id]
g <- igraph::graph_from_data_frame(ohta[ohta_D_q>0.9 & Chr1!=Chr2 & MLG1 %in% keep & MLG2 %in% keep ,.(MLG1,MLG2)], directed = FALSE)
comps <- igraph::components(g)
cls <- split(names(comps$membership), comps$membership)
#cls <- cls[vapply(cls, length, integer(1)) >= l_min]

# data.table(t(dstat_unphased_fast(
#   index = which(colnames(CMLGs) %in% cls$`1`)[1:2],
#   prep = prep,
#   tot_maf = 0,
#   pop_maf = 0
# )))

candidate_CMLGs <- cls$`2`
image(t(polarize_genotypes(CeMLGs_1mb$eMLG[,candidate_CMLGs])))
heatmap(polarize_genotypes(CeMLGs_1mb$eMLG[,candidate_CMLGs]),scale="none")

markers <- map_CeMLG[C_eMLG_id %in% candidate_CMLGs,marker]

image(t(polarize_genotypes(GTs_1mb[,markers])))
heatmap(GTs_1mb[,markers],scale="none")
