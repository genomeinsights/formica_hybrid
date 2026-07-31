## =============================================================================
## Module E -- Chr14 2-4 Mb heatmaps: individual genotypes + population allele freqs
## =============================================================================
## The raw data behind moduleE_fst_vs_adjld_emp_vs_sim (Chr14 2-4 Mb, 141 diagnostic-
## rich SNPs from the 15k set), shown as heatmaps for empirical vs the neutral null.
## UNSORTED: columns in genomic-position order, rows in natural sample order (no
## clustering, no DI sorting). Oriented to aquilonia dosage so colour = ancestry
## (blue = polyctena allele, red = aquilonia allele), as in Module A.
##   (1) individual genotypes: rows = individuals (grouped by population), cols = SNPs
##   (2) population allele frequencies: rows = 20 populations, cols = SNPs
## The point: empirical adjacent columns change WHICH populations are red/blue
## (locus-specific mosaic); the null shows smoother ancestry blocks (admixture LD).
## Output: Figures/moduleE_chr14_genotypes.{pdf,png}, Figures/moduleE_chr14_popfreq.{pdf,png}
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
F<-"/Users/petrikem/gitlab/formica_hybrid"
FVCF<-file.path(F,"moduleE_slim/founders/maf015_DIstrat1500_burnin"); FXDIR<-file.path(F,"moduleE_slim/output_replicates/rep01")
TAG<-"Naq450_Npol195"; K<-6250; GEN<-60L
b<-readRDS(file.path(F,"moduleE_slim/inputs/empirical_bundle.rds"))
fm<-rbindlist(lapply(list.files(FVCF,"founders_ch.*[.]vcf$",full.names=TRUE),
   function(f) fread(f,skip="#CHROM",select=c(1,2,3),header=TRUE,col.names=c("chs","Pos","marker"))))
fm[,chr_id:=as.integer(sub("ch","",chs))]; setorder(fm,chr_id,Pos)
markers<-fm$marker; faq<-b$f_aq_par[markers]; fpol<-b$f_pol_par[markers]
sel<-which(fm$chr_id==14 & fm$Pos>=2e6 & fm$Pos<=4e6); sel<-sel[order(fm$Pos[sel])]
smk<-markers[sel]; spos<-fm$Pos[sel]; aqsign<-ifelse(faq[sel]>=fpol[sel],1L,-1L)   # +1: ALT=aq allele
cat(sprintf("Chr14 2-4Mb SNPs: %d\n",length(sel)))
orient<-function(D){ sweep(D,2,ifelse(aqsign>0,0,2),`-`)*rep(aqsign,each=nrow(D)) } # -> aquilonia dosage 0..2

read_dose<-function(f){dt<-fread(f,skip="#CHROM",header=TRUE,sep="\t",showProgress=FALSE);mk<-paste(sub("ch","Chr",dt[[1]]),dt[[2]],sep=":");G<-as.matrix(dt[,10:ncol(dt)]);gt<-sub(":.*$","",G);d<-matrix(as.integer(sub("[|/].*$","",gt))+as.integer(sub("^.*[|/]","",gt)),nrow=nrow(G));m<-match(markers,mk);keep<-!is.na(m);D<-matrix(0L,length(markers),ncol(d));D[keep,]<-d[m[keep],,drop=FALSE];t(D)}

## ---- empirical genotypes (rows grouped by population) -----------------------
e<-new.env();load(file.path(F,"data/hybrids_only_maf005.Rdata"),envir=e);sdh<-as.data.table(e$sample_data);maph<-as.data.table(e$map_hyb_005)
hyb<-sdh[!Population%in%c("aquilonia_parent","polyctena_parent"),Sample_ID]
Ge<-e$GTs_hybrids_005[match(hyb,rownames(e$GTs_hybrids_005)),match(smk,maph$marker),drop=FALSE];Ge[is.na(Ge)]<-0L
popv<-sdh[match(hyb,Sample_ID),Population];rm(e);gc()
poplev<-unique(popv); ord<-order(match(popv,poplev)); Ge<-Ge[ord,,drop=FALSE]; popv<-popv[ord]
Ea<-orient(Ge)                                                   # individuals x SNP, aq dosage
emp_ns<-b$emp_ns

## ---- null genotypes (downsampled to matched sizes; rows grouped by deme) -----
sim_list<-lapply(1:20,function(i){D<-read_dose(file.path(FXDIR,sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf",TAG,K,i,GEN)))[,sel,drop=FALSE];set.seed(400+i);D[sample(nrow(D),min(emp_ns[i],nrow(D))),,drop=FALSE]})
Sa<-orient(do.call(rbind,sim_list)); simpop<-rep(paste0("deme",1:20),sapply(sim_list,nrow))

## ---- population allele frequencies ------------------------------------------
Pe<-t(vapply(split(seq_along(popv),factor(popv,poplev)),function(ix) colMeans(Ea[ix,,drop=FALSE])/2,numeric(length(sel))))
Ps<-t(vapply(sim_list,function(D) colMeans(orient(D))/2,numeric(length(sel))))
rownames(Ps)<-paste0("deme",1:20)

## ---- long tables for plotting -----------------------------------------------
xb<-pretty(spos/1e6,5); xbreaks<-sapply(xb,function(v) which.min(abs(spos/1e6-v)))
melt_mat<-function(A){dt<-CJ(ci=seq_len(ncol(A)),ri=seq_len(nrow(A))); dt[,val:=as.vector(A)]; dt[]}  # column-major: ri fastest
hm_indiv<-function(A,poprows,title){
  dt<-melt_mat(A); setnames(dt,"val","dose")
  bnd<-cumsum(rle(as.character(poprows))$lengths); mids<-bnd-diff(c(0,bnd))/2; labs<-rle(as.character(poprows))$values
  ggplot(dt,aes(ci,ri,fill=dose))+geom_raster()+
    scale_fill_gradient2(low="#2166ac",mid="#f7f7f7",high="#b2182b",midpoint=1,limits=c(0,2),name="aq dosage")+
    scale_x_continuous(breaks=xbreaks,labels=sprintf("%.1f",spos[xbreaks]/1e6),expand=c(0,0))+
    scale_y_continuous(breaks=mids,labels=labs,expand=c(0,0))+
    geom_hline(yintercept=bnd[-length(bnd)]+.5,colour="grey60",linewidth=.2)+
    labs(x="position on Chr14 (Mb)",y=NULL,title=title)+theme_minimal(base_size=8)+
    theme(panel.grid=element_blank(),axis.text.y=element_text(size=5),legend.position="bottom",legend.key.height=unit(.3,"cm"))
}
hm_freq<-function(P,title){
  dt<-melt_mat(P); setnames(dt,"val","freq")                          # ri=1 (first pop) at bottom
  ggplot(dt,aes(ci,ri,fill=freq))+geom_raster()+
    scale_fill_gradient2(low="#2166ac",mid="#f7f7f7",high="#b2182b",midpoint=.5,limits=c(0,1),name="aq allele freq")+
    scale_x_continuous(breaks=xbreaks,labels=sprintf("%.1f",spos[xbreaks]/1e6),expand=c(0,0))+
    scale_y_continuous(breaks=1:nrow(P),labels=rownames(P),expand=c(0,0))+
    labs(x="position on Chr14 (Mb)",y=NULL,title=title)+theme_minimal(base_size=8)+
    theme(panel.grid=element_blank(),axis.text.y=element_text(size=6),legend.position="bottom",legend.key.height=unit(.3,"cm"))
}
gi<-hm_indiv(Ea,popv,"empirical: individual genotypes") | hm_indiv(Sa,simpop,"neutral null: individual genotypes")
gi<-gi+plot_annotation(title="Chr14 2-4 Mb -- individual genotypes (aquilonia dosage), columns=SNPs by position, rows=individuals by population (UNSORTED)",
   theme=theme(plot.title=element_text(face="bold",size=10)))
ggsave(file.path(F,"Figures/moduleE_chr14_genotypes.png"),gi,width=14,height=8,dpi=150)
ggsave(file.path(F,"Figures/moduleE_chr14_genotypes.pdf"),gi,width=14,height=8)
gf<-hm_freq(Pe,"empirical: population allele frequencies") | hm_freq(Ps,"neutral null: deme allele frequencies")
gf<-gf+plot_annotation(title="Chr14 2-4 Mb -- population/deme aquilonia allele frequencies, columns=SNPs by position, rows=populations (UNSORTED)",
   theme=theme(plot.title=element_text(face="bold",size=10)))
ggsave(file.path(F,"Figures/moduleE_chr14_popfreq.png"),gf,width=14,height=5,dpi=150)
ggsave(file.path(F,"Figures/moduleE_chr14_popfreq.pdf"),gf,width=14,height=5)
cat("saved genotype + popfreq heatmaps\n")
