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
BP_PC1 <- fread("./out_baypass/PC1_DIEM_summary_betai_reg.out")

r = BP_PC1[,rank(`BF(dB)`) / (length(`BF(dB)`) + 1)]
map[,F_BF_PC1 :=  qf(r, df1 = 1, df2 = 15)]

BP_PC2 <- fread("./out_baypass/PC2_DIEM_summary_betai_reg.out")
r = BP_PC2[,rank(`BF(dB)`) / (length(`BF(dB)`) + 1)]
map[,F_BF_PC2 :=  qf(r, df1 = 1, df2 = 15)]

## read in LD-decay data
ld_decay <- readRDS("./data/ld_decay.rds")

gds <- create_gds_from_geno(geno=GTs, map, "gds_formicia")

ld_ws <- precalculate_ld_w(c(seq(0.75,1,by=0.05),0.99),ld_decay)

draws_075_03_4_30 <- ld_rho_draws(gds,
                                        ld_decay  = ld_decay,
                                        F_vals     = map[,.(F_BF_PC1,F_BF_PC2)],
                                        q_vals     = NULL,
                                        C_scores   = NULL,
                                        n_draws    = 25,
                                        stat_type  = "q",
                                        rho        = NULL,
                                        ld_ws      = ld_ws,
                                        rho_d_lim  = list(min=0.9,max=0.99),
                                        rho_ld_lim = list(min=0.9,max=0.99),
                                        alpha_lim  = list(min=0.3,max=4),
                                        lmin_lim   = list(min=1,max=30),
                                        cores      = 10,
                                        mode       = c("per_method")
)

map[,C_PC1:= add_consistency_to_map(map, consistency_obj = consistency_score(draws_075_03_4_30$draws[method=="F_BF_PC1_prime"]))$F_BF_PC1_prime_C]
map[,C_PC2:= add_consistency_to_map(map, consistency_obj = consistency_score(draws_075_03_4_30$draws[method=="F_BF_PC2_prime"]))$F_BF_PC2_prime_C]

map[,C_joint := pmax(C_PC1,C_PC2,na.rm = TRUE)]

map_manh <- add_ORs(gds, ld_decay, map, stat = c("C_joint"),
                    sign_th = 0.2, sign_if = "greater", mode = "joint", rho_d = 0.995,
                    rho_ld = 0.995)

setnames(map_manh,"Pos","bp")
layout <- prep_manhattan(map_manh[,])

p_PC1 <- plot_manhattan_gg(layout, y_vars = c("C_joint","C_PC1","C_PC2"), y_labels = c("C_joint","C_PC1","C_PC2"),
                         thresholds = c(0.2,0.2,0.2), col_var = "OR_id",
                         point_size = 1, ncol = 1)
p_PC1
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
