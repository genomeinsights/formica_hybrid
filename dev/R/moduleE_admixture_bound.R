## =============================================================================
## Module E -- ANALYTIC BOUND: how much of the diagnostic F_ST can admixture-
## proportion variation (a single genome-wide ancestry axis) actually explain?
## =============================================================================
## Empirical-data-only argument (no simulation). If population i has aquilonia ancestry
## alpha_i, its expected frequency at a locus with parental difference delta=faq-fpol is
## p_i = fpol + alpha_i*delta, so among-population variance from admixture ALONE is
## delta^2 * var(alpha). We compare the observed among-population F_ST vs DI to the F_ST
## predicted from this pure-admixture model using the EMPIRICAL alpha spread (sd~0.07).
## Result: admixture explains <5% of the diagnostic F_ST; ~95% is locus-specific.
## Output: Figures/moduleE_admixture_bound.{pdf,png} + data/moduleE_sim/moduleE_admixture_bound.rds
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2) })
F<-"/Users/petrikem/gitlab/formica_hybrid"
b<-readRDS(file.path(F,"moduleE_slim/inputs/empirical_bundle.rds")); uni<-as.data.table(b$universe); setkey(uni,marker)
FVCF<-file.path(F,"moduleE_slim/founders/maf015_DIstrat1500_burnin")
fm<-rbindlist(lapply(list.files(FVCF,"founders_ch.*[.]vcf$",full.names=TRUE),
   function(f) fread(f,skip="#CHROM",select=3,header=TRUE,col.names="marker")))
markers<-fm$marker; DI<-uni[markers,DI]; faq<-b$f_aq_par[markers]; fpol<-b$f_pol_par[markers]
mm<-match(markers,colnames(b$emp_mean)); P<-b$emp_mean[,mm]/2000        # 20 x M per-pop freq
HI<-readRDS(file.path(F,"data/moduleE_sim/empirical_hybrid_index.rds"))$HI

BR<-c(-Inf,-90,-75,-60,-50,-40,-30,-25,-20,-15,Inf); dib<-cut(DI,BR)
mids<-sapply(levels(dib),function(l) mean(DI[dib==l],na.rm=TRUE))
finite<-apply(P,2,function(c)all(is.finite(c))) & is.finite(faq)&is.finite(fpol)
Fst_bin<-function(M,keep){pb<-colMeans(M);Ht<-2*pb*(1-pb);Hs<-colMeans(2*M*(1-M));ok<-Ht>0&keep
  data.table(dib,Ht,Hs,ok)[ok==TRUE,.(v=sum(Ht-Hs)/sum(Ht)),keyby=dib]$v}
Ppred<-outer(HI,faq)+outer(1-HI,fpol)                                    # pure-admixture prediction, per locus
obs<-Fst_bin(P,finite); pred<-Fst_bin(Ppred,finite)

## headline numbers on diagnostic loci
dg<-which(DI> -25 & finite & abs(faq-fpol)>0.3)
Fst_g<-function(M){pb<-colMeans(M[,dg]);Ht<-2*pb*(1-pb);Hs<-colMeans(2*M[,dg]*(1-M[,dg]));ok<-Ht>0;sum((Ht-Hs)[ok])/sum(Ht[ok])}
o<-Fst_g(P); p<-Fst_g(Ppred)
cat(sprintf("diagnostic F_ST: observed=%.3f  admixture-predicted=%.3f  -> admixture explains %.1f%%\n",o,p,100*p/o))
saveRDS(list(mids=mids,obs=obs,pred=pred,diag_obs=o,diag_pred=p,alpha_sd=sd(HI)),
        file.path(F,"data/moduleE_sim/moduleE_admixture_bound.rds"))

dt<-rbindlist(list(data.table(DI=mids,v=obs,src="observed empirical"),
                   data.table(DI=mids,v=pred,src="admixture-proportion variation\n(single ancestry axis, empirical sd_alpha=0.07)")))
cols<-c("observed empirical"="#d95f02","admixture-proportion variation\n(single ancestry axis, empirical sd_alpha=0.07)"="#1b6ca8")
g<-ggplot(dt,aes(DI,v,colour=src))+geom_line(linewidth=1)+geom_point(size=2)+scale_colour_manual(values=cols,name=NULL)+
  annotate("text",x=-20,y=o,label=sprintf("observed %.2f",o),hjust=1.1,vjust=-0.6,colour="#d95f02",size=3.3)+
  annotate("text",x=-20,y=p,label=sprintf("admixture ceiling %.3f (%.0f%%)",p,100*p/o),hjust=1.1,vjust=1.6,colour="#1b6ca8",size=3.3)+
  labs(x="DiagnosticIndex",y=expression(among-population~F[ST]),
       title="Admixture-proportion variation cannot produce the diagnostic F_ST",
       subtitle="observed rise vs the ceiling a single genome-wide ancestry axis (empirical spread) allows")+
  theme_bw(base_size=11)+theme(legend.position="bottom",legend.text=element_text(size=8))
ggsave(file.path(F,"Figures/moduleE_admixture_bound.pdf"),g,width=8,height=5)
ggsave(file.path(F,"Figures/moduleE_admixture_bound.png"),g,width=8,height=5,dpi=140)
cat("saved figure + rds\n")
