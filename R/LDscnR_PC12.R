#pak::pak("genomeinsights/LDscnR")

library(ggplot2)
library(data.table)
library(SNPRelate)
library(LDscnR)

#### --- read in parsed data --- ####
data <- readRDS("./data/filtered_DIEM_1mb.rds")
map <- data$map
GTs <- data$GTs
sample_info <- data$sample_info
rm(data)
gc()



## read and transform BF to a pseudo-F scale
## from ./R/baypass.R
BP_PC1 <- fread("./out_baypass/PC1_summary_betai_reg.out")

r = BP_PC1[,rank(`BF(dB)`) / (length(`BF(dB)`) + 1)]
map[,F_BF_PC1 :=  qf(r, df1 = 1, df2 = 18)]

BP_PC2 <- fread("./out_baypass/PC2_summary_betai_reg.out")
r = BP_PC2[,rank(`BF(dB)`) / (length(`BF(dB)`) + 1)]
map[,F_BF_PC2 :=  qf(r, df1 = 1, df2 = 18)]
#map[,cor.test(F_BF_PC2,F_BF_PC1)]

## read and transform BF to a pseudo-F scale
xtx_PC1 <- fread("./out_baypass/PC1_DIEM_summary_pi_xtx.out")
xtx_PC2 <- fread("./out_baypass/PC2_DIEM_summary_pi_xtx.out")
#plot(xtx_PC1$XtXst,xtx_PC2$XtXst)
r = BP_PC1[,rank(`BF(dB)`) / (length(`BF(dB)`) + 1)]
map[,F_BF_PC1 :=  qf(r, df1 = 1, df2 = 18)]

BP_PC2 <- fread("./out_baypass/PC2_summary_betai_reg.out")
r = BP_PC2[,rank(`BF(dB)`) / (length(`BF(dB)`) + 1)]
map[,F_BF_PC2 :=  qf(r, df1 = 1, df2 = 18)]
map[,cor.test(F_BF_PC2,F_BF_PC1)]
## read in LD-decay data: from ./R/LD_decay.R
ld_decay <- readRDS("./data/ld_decay_corr.rds")
ld_decay$decay_sum

#snpgdsClose(gds); unlink(gds)
gds <- create_gds_from_geno(geno=GTs, map, "gds_formicia")

if(FALSE){
  ld_ws <- precalculate_ld_w(rho=c(seq(0.75,0.95,by=0.05),0.99),ld_decay)
  saveRDS(ld_ws,"./data/ld_ws.rds")
}else{
  ld_ws <- readRDS("./data/ld_ws.rds")
}

draws_075_03_2_30 <- ld_rho_draws(gds,
                                        ld_decay  = ld_decay,
                                        F_vals     = map[,.(F_BF_PC1,F_BF_PC2)],
                                        n_draws    = 25,
                                        stat_type  = "q",
                                        rho        = NULL,
                                        ld_ws      = ld_ws,
                                        rho_d_lim  = list(min=0.9,max=0.99),
                                        rho_ld_lim = list(min=0.9,max=0.99),
                                        alpha_lim  = list(min=0.3,max=2),
                                        lmin_lim   = list(min=1,max=30),
                                        cores      = 8,
                                        mode       = c("per_method")
)


map[,C_PC1:= add_consistency_to_map(map, consistency_obj = consistency_score(draws_075_03_2_30$draws[method=="F_BF_PC1_prime"]))$F_BF_PC1_prime_C]
map[,C_PC2:= add_consistency_to_map(map, consistency_obj = consistency_score(draws_075_03_2_30$draws[method=="F_BF_PC2_prime"]))$F_BF_PC2_prime_C]

map[,C_joint := pmax(C_PC1,C_PC2,na.rm = TRUE)]

plot(ld_ws[,colnames(ld_ws)=="1"],cex=0.25)

map[,ld_w_95 := ld_ws[,colnames(ld_ws)=="0.95"]]
map[,ld_w_99 := ld_ws[,colnames(ld_ws)=="0.99"]]
#saveRDS(list(map=map,ld_ws=ld_ws,draws=draws_075_03_2_30),"LDscnR_data.rds")


map_manh <- add_ORs(gds, ld_decay, map, stat = c("C_joint"),
                    sign_th = 0.2, sign_if = "greater", mode = "joint", rho_d = 0.995,
                    rho_ld = 0.995)


layout <- prep_manhattan(map_manh[,])

p_PC1 <- plot_manhattan_gg(layout, y_vars = c("C_joint","C_PC1","C_PC2"), y_labels = c("C_joint","C_PC1","C_PC2"),
                         thresholds = c(0.2,0.2,0.2), col_var = "OR_id",
                         point_size = 1, ncol = 1)
ggsave("./Figures/manh_075_03_2_30.png",p_PC1,width = 8,height = 6,units = "in",dpi=300)



#
# map_manh <- add_ORs(gds, ld_decay, map, stat = "C_PC2",
#                     sign_th = 0.2, sign_if = "greater", mode = "joint", rho_d = 0.995,
#                     rho_ld = 0.995)
#
# setnames(map_manh,"Pos","bp")
# layout <- prep_manhattan(map_manh[,])
#
# p_PC2 <- plot_manhattan_gg(layout, y_vars = c("C_PC2"), y_labels = c("C_PC2"),
#                            thresholds = 0.2, col_var = "OR_id",
#                            point_size = 1, ncol = 1)
# p_PC1 / p_PC2
#
