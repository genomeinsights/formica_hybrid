## =============================================================================
## Module E -- VARIABLE-ALPHA null: does admixture-proportion variation reproduce
## the empirical F_ST-vs-DI rise (and does the LD dissociation still separate it)?
## =============================================================================
## Compares three data sets on identical markers:
##   - VARIABLE-alpha null (output_varalpha): 20 demes, per-deme alpha matched to real
##     populations (0.35-0.64), bottleneck depth held constant. K=6250 gen 60.
##   - FIXED-alpha null    (output_replicates/rep01): 20 demes, all alpha=0.536. control.
##   - EMPIRICAL.
## STEP 1 (manipulation check): realized per-deme alpha must span the empirical range
##   in the variable run and be ~flat in the fixed control -- verify BEFORE interpreting.
## STEP 2 (test): F_ST vs DI and adjacent-LD vs DI, all downsampled to matched empirical
##   population sizes. Question: does variable-alpha lift diagnostic F_ST toward empirical
##   (a neutral route to the rise), and does it keep HIGH diagnostic LD (so the empirical
##   low-LD dissociation still separates it)?
## Output: Figures/moduleE_varalpha.{pdf,png} + data/moduleE_sim/moduleE_varalpha.rds
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
F<-"/Users/petrikem/gitlab/formica_hybrid"
VADIR <- file.path(F,"moduleE_slim/output_varalpha")
FXDIR <- file.path(F,"moduleE_slim/output_replicates/rep01")     # fixed-alpha control (gen 60)
FVCF  <- file.path(F,"moduleE_slim/founders/maf015_DIstrat1500_burnin")
TAG<-"Naq450_Npol195"; K<-6250; GEN<-60L
b<-readRDS(file.path(F,"moduleE_slim/inputs/empirical_bundle.rds")); uni<-as.data.table(b$universe); setkey(uni,marker)
cfg<-fread(file.path(F,"moduleE_slim/inputs/varalpha_founding.csv"))     # deme,pop,alpha,N_aq,N_pol
fm<-rbindlist(lapply(list.files(FVCF,"founders_ch.*[.]vcf$",full.names=TRUE),
     function(f) fread(f,skip="#CHROM",select=c(1,2,3),header=TRUE,col.names=c("chs","Pos","marker"))))
fm[,`:=`(Chr=sub("ch","Chr",chs),chr_id=as.integer(sub("ch","",chs)))]; setorder(fm,chr_id,Pos)
markers<-fm$marker; DI<-uni[markers,DI]; emp_ns<-b$emp_ns
faq<-b$f_aq_par[markers]; fpol<-b$f_pol_par[markers]
BR<-c(-Inf,-90,-75,-60,-50,-40,-30,-25,-20,-15,Inf); dib<-cut(DI,BR)
mids<-sapply(levels(dib),function(l) mean(DI[dib==l],na.rm=TRUE)); DIAG<-which(DI> -25)
adx<-which(DI> -25 & is.finite(faq)&is.finite(fpol)&abs(faq-fpol)>0.3)      # loci for the hybrid index

read_dose<-function(f){dt<-fread(f,skip="#CHROM",header=TRUE,sep="\t",showProgress=FALSE);mk<-paste(sub("ch","Chr",dt[[1]]),dt[[2]],sep=":");G<-as.matrix(dt[,10:ncol(dt)]);gt<-sub(":.*$","",G);d<-matrix(as.integer(sub("[|/].*$","",gt))+as.integer(sub("^.*[|/]","",gt)),nrow=nrow(G));m<-match(markers,mk);keep<-!is.na(m);D<-matrix(0L,length(markers),ncol(d));D[keep,]<-d[m[keep],,drop=FALSE];t(D)}
fst_bin<-function(Fr){pb<-colMeans(Fr,na.rm=T);Ht<-2*pb*(1-pb);Hs<-colMeans(2*Fr*(1-Fr),na.rm=T);data.table(dib,Ht,Hs)[,.(v=sum(Ht-Hs,na.rm=T)/sum(Ht,na.rm=T)),keyby=dib]$v}
adjld<-function(G){res<-rep(NA_real_,length(markers));for(cid in unique(fm$chr_id)){j<-which(fm$chr_id==cid);if(length(j)<2)next;for(a in seq_along(j)){nb<-j[c(a-1,a+1)];nb<-nb[nb>=min(j)&nb<=max(j)&nb!=j[a]];nb<-nb[!is.na(nb)];if(!length(nb))next;res[j[a]]<-mean(suppressWarnings(cor(G[,j[a]],G[,nb,drop=FALSE]))^2,na.rm=T)}};res}
ld_bin<-function(a){data.table(dib,a)[,.(v=median(a,na.rm=T)),keyby=dib]$v}
alpha_of<-function(freq){v<-(freq[adx]-fpol[adx])/(faq[adx]-fpol[adx]); mean(pmin(pmax(v,0),1),na.rm=TRUE)}

## file path per deme
vfile<-function(d){cf<-cfg[deme==d]; file.path(VADIR,sprintf("hyb_neutral_realfounders_Naq%d_Npol%d_K%d_rep%d_gen%d.vcf",cf$N_aq,cf$N_pol,K,d,GEN))}
ffile<-function(d) file.path(FXDIR,sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf",TAG,K,d,GEN))

## ---- read both sims (full samples) ------------------------------------------
readset<-function(fn){lapply(1:20,function(d){f<-fn(d); if(!file.exists(f)){warning("missing ",f);return(NULL)}; read_dose(f)})}
VA<-readset(vfile); FX<-readset(ffile)

## ---- STEP 1: manipulation check -- realized per-deme alpha -------------------
al_va<-sapply(VA,function(D) if(is.null(D)) NA else alpha_of(colMeans(D)/2))
al_fx<-sapply(FX,function(D) if(is.null(D)) NA else alpha_of(colMeans(D)/2))
al_emp<-readRDS(file.path(F,"data/moduleE_sim/empirical_hybrid_index.rds"))$HI
cat("=== STEP 1: realized alpha spread (manipulation check) ===\n")
cat(sprintf("variable-alpha null : mean=%.3f sd=%.3f range=[%.3f,%.3f]\n",mean(al_va,na.rm=T),sd(al_va,na.rm=T),min(al_va,na.rm=T),max(al_va,na.rm=T)))
cat(sprintf("fixed-alpha control : mean=%.3f sd=%.3f range=[%.3f,%.3f]\n",mean(al_fx,na.rm=T),sd(al_fx,na.rm=T),min(al_fx,na.rm=T),max(al_fx,na.rm=T)))
cat(sprintf("empirical           : mean=%.3f sd=%.3f range=[%.3f,%.3f]\n",mean(al_emp),sd(al_emp),min(al_emp),max(al_emp)))
cat(sprintf("intended vs realized (variable) correlation: r=%.3f\n",cor(cfg[order(deme)]$alpha,al_va,use="complete.obs")))

## ---- STEP 2: F_ST-vs-DI and LD-vs-DI, downsampled to empirical sizes ---------
downcurves<-function(set){
  subs<-lapply(seq_along(set),function(i){D<-set[[i]];if(is.null(D))return(NULL);set.seed(400+i);D[sample(nrow(D),min(emp_ns[i],nrow(D))),,drop=FALSE]})
  Fr<-t(vapply(subs,function(D)colMeans(D)/2,numeric(length(markers)))); Gp<-do.call(rbind,subs); ald<-adjld(Gp)
  list(fst=fst_bin(Fr), ld=ld_bin(ald), diagfst=fst_bin(Fr)[nlevels(dib)], rise=fst_bin(Fr)[nlevels(dib)]-fst_bin(Fr)[1],
       ratio=fst_bin(Fr)[nlevels(dib)]/(median(ald[DIAG],na.rm=T)+0.02))
}
cv_va<-downcurves(VA); cv_fx<-downcurves(FX)
## empirical
Fe<-b$emp_mean[,markers,drop=FALSE]/2000     # 20 x M per-pop; downsample not needed (already pop-level means)? use full genotypes:
e<-new.env();load(file.path(F,"data/hybrids_only_maf005.Rdata"),envir=e);sdh<-as.data.table(e$sample_data);maph<-as.data.table(e$map_hyb_005)
hyb<-sdh[!Population%in%c("aquilonia_parent","polyctena_parent"),Sample_ID]
Ge<-e$GTs_hybrids_005[match(hyb,rownames(e$GTs_hybrids_005)),match(markers,maph$marker),drop=FALSE];Ge[is.na(Ge)]<-0L
popv<-sdh[match(hyb,Sample_ID),Population]; rm(e);gc()
Frp<-t(vapply(split(seq_along(popv),popv),function(ix) colMeans(Ge[ix,,drop=FALSE])/2,numeric(length(markers))))
emp_fst<-fst_bin(Frp); emp_ld<-ld_bin(adjld(Ge)); emp_diagfst<-emp_fst[nlevels(dib)]; emp_rise<-emp_fst[nlevels(dib)]-emp_fst[1]
emp_ratio<-emp_fst[nlevels(dib)]/(median(adjld(Ge)[DIAG],na.rm=T)+0.02)

cat("\n=== STEP 2: diagnostic F_ST rise + ratio ===\n")
cat(sprintf("fixed-alpha   : diag F_ST=%.3f  rise=%.3f  ratio=%.1f\n",cv_fx$diagfst,cv_fx$rise,cv_fx$ratio))
cat(sprintf("variable-alpha: diag F_ST=%.3f  rise=%.3f  ratio=%.1f\n",cv_va$diagfst,cv_va$rise,cv_va$ratio))
cat(sprintf("empirical     : diag F_ST=%.3f  rise=%.3f  ratio=%.1f\n",emp_diagfst,emp_rise,emp_ratio))

saveRDS(list(al_va=al_va,al_fx=al_fx,al_emp=al_emp,mids=mids,cv_va=cv_va,cv_fx=cv_fx,
             emp_fst=emp_fst,emp_ld=emp_ld,emp_diagfst=emp_diagfst,emp_rise=emp_rise,emp_ratio=emp_ratio),
        file.path(F,"data/moduleE_sim/moduleE_varalpha.rds"))

## ---- figure ------------------------------------------------------------------
cols<-c("fixed-alpha (control)"="#7fbf7b","variable-alpha"="#0b6b4f","empirical"="#d95f02")
dA<-rbindlist(list(data.table(a=al_fx,src="fixed-alpha (control)"),data.table(a=al_va,src="variable-alpha"),data.table(a=al_emp,src="empirical")))
pA<-ggplot(dA,aes(src,a,colour=src))+geom_jitter(width=.12,height=0,size=1.6)+
  stat_summary(fun=mean,geom="crossbar",width=.4,colour="grey30")+scale_colour_manual(values=cols,guide="none")+
  labs(x=NULL,y="realized aquilonia ancestry (alpha)",
       title="A  manipulation check FAILS: realized alpha ~0.50 for all demes\n(founding count ratio does not control ancestry; only empirical varies)")+
  theme_bw(base_size=10)+theme(axis.text.x=element_text(angle=20,hjust=1))
dF<-rbindlist(list(data.table(DI=mids,v=cv_fx$fst,src="fixed-alpha (control)"),data.table(DI=mids,v=cv_va$fst,src="variable-alpha"),data.table(DI=mids,v=emp_fst,src="empirical")))
pF<-ggplot(dF,aes(DI,v,colour=src))+geom_line(linewidth=.9)+geom_point(size=1.5)+scale_colour_manual(values=cols,name=NULL)+
  labs(x="DiagnosticIndex",y=expression(F[ST]),title="B  F_ST vs DI")+theme_bw(base_size=10)+theme(legend.position="bottom")
dL<-rbindlist(list(data.table(DI=mids,v=cv_fx$ld,src="fixed-alpha (control)"),data.table(DI=mids,v=cv_va$ld,src="variable-alpha"),data.table(DI=mids,v=emp_ld,src="empirical")))
pL<-ggplot(dL,aes(DI,v,colour=src))+geom_line(linewidth=.9)+geom_point(size=1.5)+scale_colour_manual(values=cols,name=NULL)+
  labs(x="DiagnosticIndex",y="adjacent-SNP LD",title="C  adjacent-LD vs DI (the dissociation)")+theme_bw(base_size=10)+theme(legend.position="bottom")
comp<-pA|pF|pL
comp<-comp+plot_annotation(title="Variable-alpha attempt: founding COUNT ratios do NOT create ancestry variation in this haplodiploid model",
  subtitle="realized alpha ~0.50 for all demes (queen[pol] x male[aq] pins F1 at 50%); test is invalid until founders carry mixed ancestry in both sexes",
  theme=theme(plot.title=element_text(face="bold",size=12)))
ggsave(file.path(F,"Figures/moduleE_varalpha.pdf"),comp,width=14,height=5.5)
ggsave(file.path(F,"Figures/moduleE_varalpha.png"),comp,width=14,height=5.5,dpi=140)
cat("saved figure + rds\n")
