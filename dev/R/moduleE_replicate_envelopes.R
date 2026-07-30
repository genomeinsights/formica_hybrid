## =============================================================================
## Module E -- whole-experiment null ENVELOPES (assessment section 3)
## =============================================================================
## The 20 demes are ONE realization; here we repeat the complete 20-deme experiment
## 20 times (independent seeds + fresh 450/195 founder draws from the diversified
## burn-in pool; run_replicates.sh) and take the distribution ACROSS replicates as
## the null envelope for the experiment-level statistics:
##   - F_ST vs DI (dose-response)      - adjacent-SNP LD vs DI
##   - the diagnostic F_ST/LD ratio    - the neutral->diagnostic F_ST rise (a slope)
## The empirical value is then read against the envelope (z-score / outside-range).
## Output: Figures/moduleE_replicate_envelopes.{pdf,png} + data/moduleE_sim/moduleE_replicate_envelopes.rds
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
F      <- "/Users/petrikem/gitlab/formica_hybrid"
REPDIR <- file.path(F, "moduleE_slim/output_replicates")
FVCF   <- file.path(F, "moduleE_slim/founders/maf015_DIstrat1500_burnin")
TAG<-"Naq450_Npol195"; K<-6250; GEN<-60L
b <- readRDS(file.path(F,"moduleE_slim/inputs/empirical_bundle.rds")); uni<-as.data.table(b$universe); setkey(uni,marker)
fm <- rbindlist(lapply(list.files(FVCF,"founders_ch.*[.]vcf$",full.names=TRUE),
       function(f) fread(f,skip="#CHROM",select=c(1,2,3),header=TRUE,col.names=c("chs","Pos","marker"))))
fm[,`:=`(Chr=sub("ch","Chr",chs),chr_id=as.integer(sub("ch","",chs)))]; setorder(fm,chr_id,Pos)
markers<-fm$marker; DI<-uni[markers,DI]; emp_ns<-b$emp_ns
BR<-c(-Inf,-90,-75,-60,-50,-40,-30,-25,-20,-15,Inf); dib<-cut(DI,BR)
mids<-sapply(levels(dib),function(l) mean(DI[dib==l],na.rm=TRUE)); DIAG<-which(DI> -25); NEU<-which(DI<= -90)

read_dose<-function(f){dt<-fread(f,skip="#CHROM",header=TRUE,sep="\t",showProgress=FALSE);mk<-paste(sub("ch","Chr",dt[[1]]),dt[[2]],sep=":");G<-as.matrix(dt[,10:ncol(dt)]);gt<-sub(":.*$","",G);d<-matrix(as.integer(sub("[|/].*$","",gt))+as.integer(sub("^.*[|/]","",gt)),nrow=nrow(G));m<-match(markers,mk);keep<-!is.na(m);D<-matrix(0L,length(markers),ncol(d));D[keep,]<-d[m[keep],,drop=FALSE];t(D)}  # robust to monomorphic markers dropped from the VCF (filled 0)
fst_bin<-function(Fr){pb<-colMeans(Fr,na.rm=T);Ht<-2*pb*(1-pb);Hs<-colMeans(2*Fr*(1-Fr),na.rm=T);data.table(dib,Ht,Hs)[,.(v=sum(Ht-Hs,na.rm=T)/sum(Ht,na.rm=T)),keyby=dib]$v}
adjld<-function(G){res<-rep(NA_real_,length(markers));for(cid in unique(fm$chr_id)){j<-which(fm$chr_id==cid);if(length(j)<2)next;for(a in seq_along(j)){nb<-j[c(a-1,a+1)];nb<-nb[nb>=min(j)&nb<=max(j)&nb!=j[a]];nb<-nb[!is.na(nb)];if(!length(nb))next;res[j[a]]<-mean(suppressWarnings(cor(G[,j[a]],G[,nb,drop=FALSE]))^2,na.rm=T)}};res}
ld_bin<-function(a){data.table(dib,a)[,.(v=median(a,na.rm=T)),keyby=dib]$v}

## ---- empirical curves (from the compact bundle + 1 GB genotypes) ------------
Fe<-b$emp_mean[,markers,drop=FALSE]/2000
e<-new.env();load(file.path(F,"data/hybrids_only_maf005.Rdata"),envir=e);sdh<-as.data.table(e$sample_data);maph<-as.data.table(e$map_hyb_005)
hyb<-sdh[!Population%in%c("aquilonia_parent","polyctena_parent"),Sample_ID]
Ge<-e$GTs_hybrids_005[match(hyb,rownames(e$GTs_hybrids_005)),match(markers,maph$marker),drop=FALSE];Ge[is.na(Ge)]<-0L;rm(e);gc()
alde<-adjld(Ge)
emp_fst<-fst_bin(Fe); emp_ld<-ld_bin(alde)
emp_ratio<- emp_fst[length(levels(dib))] / ( median(alde[DIAG],na.rm=T)+0.02 )   # diagnostic F_ST / median diagnostic adjacent-LD
emp_rise<- emp_fst[length(levels(dib))] - emp_fst[1]                 # diagnostic - neutral F_ST

## ---- per-replicate null curves ---------------------------------------------
reps<-sort(list.dirs(REPDIR,recursive=FALSE))
cat(sprintf("replicates found: %d\n", length(reps)))
FST<-LD<-matrix(NA_real_,length(reps),nlevels(dib)); RATIO<-RISE<-numeric(length(reps))
for(r in seq_along(reps)){
  subs<-lapply(1:20,function(i){f<-file.path(reps[r],sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf",TAG,K,i,GEN));if(!file.exists(f))return(NULL);D<-read_dose(f);set.seed(400+i);D[sample(nrow(D),min(emp_ns[i],nrow(D))),,drop=FALSE]})
  Fr<-t(vapply(subs,function(D)colMeans(D)/2,numeric(length(markers)))); Gp<-do.call(rbind,subs); ald<-adjld(Gp)
  FST[r,]<-fst_bin(Fr); LD[r,]<-ld_bin(ald)
  RATIO[r]<-fst_bin(Fr)[nlevels(dib)]/(median(ald[DIAG],na.rm=T)+0.02); RISE[r]<-FST[r,nlevels(dib)]-FST[r,1]
  cat(sprintf("  %s done\n", basename(reps[r])))
}
env<-function(M) as.data.table(t(apply(M,2,quantile,c(.025,.5,.975),na.rm=TRUE)))
Fq<-env(FST); Lq<-env(LD); setnames(Fq,c("lo","mid","hi")); setnames(Lq,c("lo","mid","hi"))
Fq[,`:=`(DI=mids,src="neutral null (20 replicates)")]; Lq[,`:=`(DI=mids,src="neutral null (20 replicates)")]
saveRDS(list(mids=mids,FST=FST,LD=LD,RATIO=RATIO,RISE=RISE,emp_fst=emp_fst,emp_ld=emp_ld,emp_ratio=emp_ratio,emp_rise=emp_rise),
        file.path(F,"data/moduleE_sim/moduleE_replicate_envelopes.rds"))

## ---- report ----------------------------------------------------------------
zr<-(emp_ratio-mean(RATIO))/sd(RATIO); zs<-(emp_rise-mean(RISE))/sd(RISE)
cat(sprintf("\ndiagnostic F_ST/LD ratio  : empirical %.2f ; null 20-rep range [%.2f, %.2f] mean %.2f sd %.3f -> z=%.0f\n",
  emp_ratio,min(RATIO),max(RATIO),mean(RATIO),sd(RATIO),zr))
cat(sprintf("neutral->diagnostic F_ST rise: empirical %.3f ; null range [%.3f, %.3f] mean %.3f sd %.4f -> z=%.0f\n",
  emp_rise,min(RISE),max(RISE),mean(RISE),sd(RISE),zs))
cat(sprintf("empirical outside the full null range on both? %s\n", (emp_ratio>max(RATIO) & emp_rise>max(RISE))))

## ---- figure: dose-responses with null envelope + empirical -----------------
pF<-ggplot()+geom_ribbon(data=Fq,aes(DI,ymin=lo,ymax=hi),fill="#1b9e77",alpha=.3)+
  geom_line(data=Fq,aes(DI,mid),colour="#1b9e77")+
  geom_line(data=data.table(DI=mids,v=emp_fst),aes(DI,v),colour="#d95f02",linewidth=.9)+
  geom_point(data=data.table(DI=mids,v=emp_fst),aes(DI,v),colour="#d95f02",size=1.6)+
  labs(x="DiagnosticIndex",y=expression(F[ST]),title="A  F_ST vs DI: empirical vs null 20-replicate envelope")+theme_bw(base_size=10)
pL<-ggplot()+geom_ribbon(data=Lq,aes(DI,ymin=lo,ymax=hi),fill="#1b9e77",alpha=.3)+
  geom_line(data=Lq,aes(DI,mid),colour="#1b9e77")+
  geom_line(data=data.table(DI=mids,v=emp_ld),aes(DI,v),colour="#d95f02",linewidth=.9)+
  geom_point(data=data.table(DI=mids,v=emp_ld),aes(DI,v),colour="#d95f02",size=1.6)+
  labs(x="DiagnosticIndex",y="adjacent-SNP LD",title="B  adjacent-LD vs DI: empirical vs null envelope")+theme_bw(base_size=10)
pR<-ggplot(data.table(RATIO),aes(x="",y=RATIO))+geom_boxplot(width=.4,fill="#1b9e77",alpha=.4,outlier.size=.5)+
  geom_hline(yintercept=emp_ratio,colour="#d95f02",linewidth=.9)+annotate("text",x=1,y=emp_ratio,label=sprintf("empirical %.1f",emp_ratio),vjust=-0.6,colour="#d95f02",size=3)+
  scale_y_log10()+labs(x=NULL,y="diagnostic F_ST/LD ratio",title=sprintf("C  ratio: emp vs null (z=%.0f)",zr))+theme_bw(base_size=10)
comp<-(pF|pL)/pR + plot_layout(heights=c(1.3,1))+plot_annotation(title="Whole-experiment null envelope (20 replicates of the 20-deme design)",theme=theme(plot.title=element_text(face="bold",size=12)))
ggsave(file.path(F,"Figures/moduleE_replicate_envelopes.pdf"),comp,width=11,height=8)
ggsave(file.path(F,"Figures/moduleE_replicate_envelopes.png"),comp,width=11,height=8,dpi=140)
cat("saved figure + rds\n")
