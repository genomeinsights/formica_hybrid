## =============================================================================
## Module E -- VALIDATION: does raising T1 actually increase realized relatedness
## among the simulated populations? (manipulation check for the shared-source null)
## =============================================================================
## Every locus is neutral, so the shared-source sweep is only meaningful if larger T1
## produces progressively MORE shared ancestry among the 20 populations. We test this
## WITHOUT the empirical data, on the FULL 50-female sample per deme (no downsampling,
## so the small-n sampling floor that flattened the downsampled F_ST is removed), using
## neutral loci (DI <= -90):
##   (1) among-population F_ST (full samples)         -- should DECREASE with T1
##   (2) mean pairwise F_ST between demes             -- should DECREASE with T1
##   (3) allele-frequency covariance among pops:      -- should INCREASE with T1
##         mean pairwise correlation of per-deme allele-freq deviations (coancestry),
##         and the 20x20 covariance-matrix mean off-diagonal.
##   (4) divergence FROM the shared ancestral pool (T1>0): mean F_ST(deme vs ancestral
##         pool) -- the post-split divergence, should DECREASE with T1 (= T2 shrinks).
## Output: data/moduleE_sim/moduleE_sharedsource_relatedness.rds + console table +
##         Figures/moduleE_sharedsource_relatedness.{pdf,png}
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
F      <- "/Users/petrikem/gitlab/formica_hybrid"
SSDIR  <- file.path(F, "moduleE_slim/output_sharedsource")
FVCF   <- file.path(F, "moduleE_slim/founders/maf015_DIstrat1500_burnin")
TAG<-"Naq450_Npol195"; K<-6250
T1S <- c(0L,25L,45L); T2S <- 50L - T1S; names(T2S)<-T1S
LAB <- c("0"="T1=0 (independent)","25"="T1=25 (half-shared)","45"="T1=45 (near-panmictic)")

b <- readRDS(file.path(F,"moduleE_slim/inputs/empirical_bundle.rds")); uni<-as.data.table(b$universe); setkey(uni,marker)
fm <- rbindlist(lapply(list.files(FVCF,"founders_ch.*[.]vcf$",full.names=TRUE),
       function(f) fread(f,skip="#CHROM",select=c(1,2,3),header=TRUE,col.names=c("chs","Pos","marker"))))
fm[,`:=`(Chr=sub("ch","Chr",chs),chr_id=as.integer(sub("ch","",chs)))]; setorder(fm,chr_id,Pos)
markers<-fm$marker; DI<-uni[markers,DI]
NEU <- which(DI <= -90)                                  # neutral loci: demography only
cat(sprintf("neutral loci (DI<=-90): %d of %d\n", length(NEU), length(markers)))

read_dose<-function(f){dt<-fread(f,skip="#CHROM",header=TRUE,sep="\t",showProgress=FALSE);mk<-paste(sub("ch","Chr",dt[[1]]),dt[[2]],sep=":");G<-as.matrix(dt[,10:ncol(dt)]);gt<-sub(":.*$","",G);d<-matrix(as.integer(sub("[|/].*$","",gt))+as.integer(sub("^.*[|/]","",gt)),nrow=nrow(G));m<-match(markers,mk);keep<-!is.na(m);D<-matrix(0L,length(markers),ncol(d));D[keep,]<-d[m[keep],,drop=FALSE];t(D)}
# ancestral pool freqs (haploid males + diploid queens): read as haplosome dosage/allele count
anc_freq<-function(vcf){dt<-fread(vcf,skip="#CHROM",header=TRUE,sep="\t",showProgress=FALSE);mk<-paste(sub("ch","Chr",dt[[1]]),dt[[2]],sep=":");G<-as.matrix(dt[,10:ncol(dt)]);
  # males haploid ("0"/"1"), queens diploid ("a|b"): allele count per column, summed / total alleles
  cnt<-matrix(0L,nrow(G),ncol(G)); nall<-integer(ncol(G))
  hap<- !grepl("[|/]",G[1,]);
  cnt[,hap]<-matrix(as.integer(G[,hap]),nrow(G)); nall[hap]<-1L
  di<-!hap; cnt[,di]<-matrix(as.integer(sub("[|/].*$","",G[,di]))+as.integer(sub("^.*[|/]","",G[,di])),nrow(G)); nall[di]<-2L
  fr<-rowSums(cnt)/sum(nall); m<-match(markers,mk); out<-rep(NA_real_,length(markers)); out[!is.na(m)]<-fr[m[!is.na(m)]]; out}

fst_global<-function(Fr,idx){pb<-colMeans(Fr[,idx,drop=FALSE]);Ht<-2*pb*(1-pb);Hs<-colMeans(2*Fr[,idx,drop=FALSE]*(1-Fr[,idx,drop=FALSE]));sum(Ht-Hs)/sum(Ht)}
fst_pair<-function(p,q){pb<-(p+q)/2;Ht<-2*pb*(1-pb);Hs<-(2*p*(1-p)+2*q*(1-q))/2;sum(Ht-Hs,na.rm=T)/sum(Ht,na.rm=T)}

summ <- data.table()
COVmats <- list()
for (t1 in T1S) {
  gen<-T2S[as.character(t1)]; dir<-file.path(SSDIR,sprintf("T1_%02d",t1),"rep01")
  Fr <- t(vapply(1:20,function(i){D<-read_dose(file.path(dir,sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf",TAG,K,i,gen)));colMeans(D)/2},numeric(length(markers)))) # 20 x M full-sample freqs
  Frn <- Fr[,NEU,drop=FALSE]
  fst_g <- fst_global(Fr,NEU)
  pw <- combn(20,2,function(ij) fst_pair(Frn[ij[1],],Frn[ij[2],])); mean_pw<-mean(pw)
  # allele-freq covariance among pops on neutral loci (coancestry): center per locus
  Dev <- sweep(Frn,2,colMeans(Frn)); COV <- Dev %*% t(Dev) / ncol(Frn)   # 20x20
  offdiag <- COV[upper.tri(COV)]; diagv<-diag(COV)
  Corr <- COV / sqrt(outer(diagv,diagv)); mean_corr <- mean(Corr[upper.tri(Corr)])
  COVmats[[as.character(t1)]] <- COV
  # divergence from the shared ancestral pool (T1>0)
  fst_vs_anc <- NA_real_
  if (t1>0) { av<-list.files(dir,"^ancestral_.*[.]vcf$",full.names=TRUE)[1]; pa<-anc_freq(av)[NEU]
    fst_vs_anc <- mean(vapply(1:20,function(i) fst_pair(Frn[i,],pa),numeric(1))) }
  summ <- rbind(summ, data.table(T1=t1,T2=gen,label=LAB[as.character(t1)],
        neutral_Fst_full=fst_g, mean_pairwise_Fst=mean_pw,
        mean_coancestry_corr=mean_corr, mean_offdiag_cov=mean(offdiag),
        Fst_vs_ancestral=fst_vs_anc))
  cat(sprintf("T1=%2d (T2=%2d): neutral F_ST(full)=%.4f  mean-pairwise F_ST=%.4f  coancestry r=%.3f  F_ST(vs anc pool)=%s\n",
      t1,gen,fst_g,mean_pw,mean_corr, ifelse(is.na(fst_vs_anc),"--",sprintf("%.4f",fst_vs_anc))))
}
saveRDS(list(summ=summ,COVmats=COVmats), file.path(F,"data/moduleE_sim/moduleE_sharedsource_relatedness.rds"))
cat("\n=== manipulation check (full 50-sample data, neutral loci) ===\n"); print(summ[,.(T1,T2,neutral_Fst_full,mean_pairwise_Fst,mean_coancestry_corr,Fst_vs_ancestral)])

## ---- figure ------------------------------------------------------------------
mlt <- melt(summ[,.(T1,`neutral F_ST (full)`=neutral_Fst_full,`mean coancestry r`=mean_coancestry_corr)],id.vars="T1")
p1 <- ggplot(summ,aes(factor(T1),neutral_Fst_full,group=1))+geom_line(colour="#1b9e77")+geom_point(colour="#1b9e77",size=3)+
  geom_text(aes(label=sprintf("%.3f",neutral_Fst_full)),vjust=-1,size=3)+
  labs(x="split time T1 (shared generations)",y=expression(neutral~F[ST]~"(full 50-sample)"),
       title="A  among-population F_ST is FLAT across T1\n(expected: should fall as T1 rises -- it does not)")+
  ylim(0,max(summ$neutral_Fst_full)*1.25)+theme_bw(base_size=10)
p2 <- ggplot(summ,aes(factor(T1),mean_coancestry_corr,group=1))+geom_line(colour="#d95f02")+geom_point(colour="#d95f02",size=3)+
  geom_text(aes(label=sprintf("%.3f",mean_coancestry_corr)),vjust=-1.4,size=3)+
  geom_hline(yintercept=-1/19,linetype=2,colour="grey50")+
  annotate("text",x=1,y=-1/19,label="-1/(n-1): independent demes",vjust=-0.6,hjust=0,size=2.8,colour="grey40")+
  labs(x="split time T1 (shared generations)",y="mean pairwise coancestry (allele-freq r)",
       title="B  coancestry pinned at -1/(n-1) for all T1\n(no shared covariance = no realized relatedness gained)")+
  theme_bw(base_size=10)
comp <- p1|p2
comp <- comp + plot_annotation(title="Manipulation check FAILS: raising T1 does NOT increase realized relatedness among populations",
   subtitle="neutral loci (DI<=-90), full 50-female samples per deme, no downsampling; among-deme divergence is founding-draw-dominated, T1-invariant",
   theme=theme(plot.title=element_text(face="bold",size=12)))
ggsave(file.path(F,"Figures/moduleE_sharedsource_relatedness.pdf"),comp,width=11,height=5)
ggsave(file.path(F,"Figures/moduleE_sharedsource_relatedness.png"),comp,width=11,height=5,dpi=140)
cat("saved figure + rds\n")
