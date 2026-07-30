## =============================================================================
## Module E -- VARIABLE-ALPHA null: set up per-deme founding counts
## =============================================================================
## The fixed null founds every deme 450 aq / 195 pol -> the SAME admixture proportion
## alpha = 450/(450+2*195) = 0.536 for all 20 demes (zero spread). But the empirical
## populations vary in aquilonia ancestry (alpha in [0.35,0.64], mean 0.537, sd 0.071).
## Variation in admixture proportion is the one NEUTRAL mechanism that can make F_ST
## depend on DiagnosticIndex (high-DI loci respond most to alpha), so it must be tested.
##
## Here each simulated deme is assigned a real population's estimated alpha, and the
## founding counts are chosen to hit that alpha while holding the total founder
## haplotype count fixed at T_HAP = 840 (= the fixed null's 450 + 2*195), so the
## founding BOTTLENECK DEPTH is constant and only alpha varies:
##     N_aq  = round(alpha * T_HAP)                 (haploid aq males)
##     N_pol = round((T_HAP - N_aq)/2)              (diploid pol queens)
## Writes moduleE_slim/inputs/varalpha_founding.csv (deme, pop, alpha, N_aq, N_pol).
## =============================================================================
suppressMessages(library(data.table))
F<-"/Users/petrikem/gitlab/formica_hybrid"; T_HAP<-840L
POOL_AQ<-1200L; POOL_POL<-520L
hi<-readRDS(file.path(F,"data/moduleE_sim/empirical_hybrid_index.rds"))
al<-hi$HI                                            # 20 empirical alphas, named by pop
dt<-data.table(deme=seq_along(al), pop=names(al), alpha=as.numeric(al))
dt[,N_aq := as.integer(round(alpha*T_HAP))]
dt[,N_pol:= as.integer(round((T_HAP-N_aq)/2))]
dt[,alpha_realized := N_aq/(N_aq+2*N_pol)]
stopifnot(all(dt$N_aq<=POOL_AQ), all(dt$N_pol<=POOL_POL), all(dt$N_aq>0), all(dt$N_pol>0))
fwrite(dt, file.path(F,"moduleE_slim/inputs/varalpha_founding.csv"))
cat("wrote varalpha_founding.csv\n"); print(dt[order(alpha)])
cat(sprintf("\nalpha spread: mean=%.3f sd=%.3f range=[%.3f,%.3f]  (N_aq %d-%d, N_pol %d-%d)\n",
    mean(dt$alpha_realized),sd(dt$alpha_realized),min(dt$alpha_realized),max(dt$alpha_realized),
    min(dt$N_aq),max(dt$N_aq),min(dt$N_pol),max(dt$N_pol)))
