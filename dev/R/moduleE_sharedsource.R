## =============================================================================
## Module E -- SHARED-SOURCE (complex-demography) null: split-time dose-response
## =============================================================================
## Compares the empirical F_ST/LD pattern against neutral nulls that share ancestry
## before splitting into the 20 populations (run_sharedsource.sh). Split time T1
## (total 50 gen) is the lever: T1=0 fully independent, T1=25 half-shared, T1=45
## near-panmictic. For each T1 we read the 20 sub-demes (NREP=1 proof), downsample
## each to its matched empirical population size, and compute the same experiment-
## level curves as the envelope analysis:
##   - F_ST vs DI        - adjacent-SNP LD vs DI       - diagnostic F_ST/LD ratio + rise
## Overlaid against the empirical curve, this shows whether shared ancestry moves the
## diagnostic dissociation toward empirical (needle moves) or only lowers neutral F_ST.
## Output: Figures/moduleE_sharedsource.{pdf,png} + data/moduleE_sim/moduleE_sharedsource.rds
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
F      <- "/Users/petrikem/gitlab/formica_hybrid"
SSDIR  <- file.path(F, "moduleE_slim/output_sharedsource")
FVCF   <- file.path(F, "moduleE_slim/founders/maf015_DIstrat1500_burnin")
TAG<-"Naq450_Npol195"; K<-6250
T1S <- c(0L,25L,45L); T2S <- 50L - T1S; names(T2S)<-T1S      # gen suffix per split time
LAB <- c("0"="T1=0 (independent)","25"="T1=25 (half-shared)","45"="T1=45 (near-panmictic)")

b <- readRDS(file.path(F,"moduleE_slim/inputs/empirical_bundle.rds")); uni<-as.data.table(b$universe); setkey(uni,marker)
fm <- rbindlist(lapply(list.files(FVCF,"founders_ch.*[.]vcf$",full.names=TRUE),
       function(f) fread(f,skip="#CHROM",select=c(1,2,3),header=TRUE,col.names=c("chs","Pos","marker"))))
fm[,`:=`(Chr=sub("ch","Chr",chs),chr_id=as.integer(sub("ch","",chs)))]; setorder(fm,chr_id,Pos)
markers<-fm$marker; DI<-uni[markers,DI]; emp_ns<-b$emp_ns
BR<-c(-Inf,-90,-75,-60,-50,-40,-30,-25,-20,-15,Inf); dib<-cut(DI,BR)
mids<-sapply(levels(dib),function(l) mean(DI[dib==l],na.rm=TRUE)); DIAG<-which(DI> -25)

read_dose<-function(f){dt<-fread(f,skip="#CHROM",header=TRUE,sep="\t",showProgress=FALSE);mk<-paste(sub("ch","Chr",dt[[1]]),dt[[2]],sep=":");G<-as.matrix(dt[,10:ncol(dt)]);gt<-sub(":.*$","",G);d<-matrix(as.integer(sub("[|/].*$","",gt))+as.integer(sub("^.*[|/]","",gt)),nrow=nrow(G));m<-match(markers,mk);keep<-!is.na(m);D<-matrix(0L,length(markers),ncol(d));D[keep,]<-d[m[keep],,drop=FALSE];t(D)}
fst_bin<-function(Fr){pb<-colMeans(Fr,na.rm=T);Ht<-2*pb*(1-pb);Hs<-colMeans(2*Fr*(1-Fr),na.rm=T);data.table(dib,Ht,Hs)[,.(v=sum(Ht-Hs,na.rm=T)/sum(Ht,na.rm=T)),keyby=dib]$v}
adjld<-function(G){res<-rep(NA_real_,length(markers));for(cid in unique(fm$chr_id)){j<-which(fm$chr_id==cid);if(length(j)<2)next;for(a in seq_along(j)){nb<-j[c(a-1,a+1)];nb<-nb[nb>=min(j)&nb<=max(j)&nb!=j[a]];nb<-nb[!is.na(nb)];if(!length(nb))next;res[j[a]]<-mean(suppressWarnings(cor(G[,j[a]],G[,nb,drop=FALSE]))^2,na.rm=T)}};res}
ld_bin<-function(a){data.table(dib,a)[,.(v=median(a,na.rm=T)),keyby=dib]$v}

## ---- empirical curves --------------------------------------------------------
Fe<-b$emp_mean[,markers,drop=FALSE]/2000
e<-new.env();load(file.path(F,"data/hybrids_only_maf005.Rdata"),envir=e);sdh<-as.data.table(e$sample_data);maph<-as.data.table(e$map_hyb_005)
hyb<-sdh[!Population%in%c("aquilonia_parent","polyctena_parent"),Sample_ID]
Ge<-e$GTs_hybrids_005[match(hyb,rownames(e$GTs_hybrids_005)),match(markers,maph$marker),drop=FALSE];Ge[is.na(Ge)]<-0L;rm(e);gc()
alde<-adjld(Ge)
emp_fst<-fst_bin(Fe); emp_ld<-ld_bin(alde)
emp_ratio<- emp_fst[length(levels(dib))]/(median(alde[DIAG],na.rm=T)+0.02); emp_rise<-emp_fst[length(levels(dib))]-emp_fst[1]

## ---- per-split-time null curves (NREP=1 proof) ------------------------------
res <- list(); ratios<-rises<-setNames(numeric(length(T1S)),as.character(T1S))
FSTL<-LDL<-list()
for (t1 in T1S) {
  gen<-T2S[as.character(t1)]; dir<-file.path(SSDIR,sprintf("T1_%02d",t1),"rep01")
  subs<-lapply(1:20,function(i){f<-file.path(dir,sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf",TAG,K,i,gen));if(!file.exists(f)){warning("missing ",f);return(NULL)};D<-read_dose(f);set.seed(400+i);D[sample(nrow(D),min(emp_ns[i],nrow(D))),,drop=FALSE]})
  Fr<-t(vapply(subs,function(D)colMeans(D)/2,numeric(length(markers)))); Gp<-do.call(rbind,subs); ald<-adjld(Gp)
  FSTL[[as.character(t1)]]<-fst_bin(Fr); LDL[[as.character(t1)]]<-ld_bin(ald)
  ratios[as.character(t1)]<-fst_bin(Fr)[nlevels(dib)]/(median(ald[DIAG],na.rm=T)+0.02); rises[as.character(t1)]<-FSTL[[as.character(t1)]][nlevels(dib)]-FSTL[[as.character(t1)]][1]
  cat(sprintf("T1=%2d (gen %2d): diag F_ST=%.3f  neutral F_ST=%.3f  rise=%.3f  ratio=%.1f\n",
    t1,gen,FSTL[[as.character(t1)]][nlevels(dib)],FSTL[[as.character(t1)]][1],rises[as.character(t1)],ratios[as.character(t1)]))
}
saveRDS(list(mids=mids,FSTL=FSTL,LDL=LDL,ratios=ratios,rises=rises,
             emp_fst=emp_fst,emp_ld=emp_ld,emp_ratio=emp_ratio,emp_rise=emp_rise,LAB=LAB),
        file.path(F,"data/moduleE_sim/moduleE_sharedsource.rds"))
cat(sprintf("\nempirical: diag F_ST=%.3f  neutral F_ST=%.3f  rise=%.3f  ratio=%.1f\n",
  emp_fst[nlevels(dib)],emp_fst[1],emp_rise,emp_ratio))

## ---- figure ------------------------------------------------------------------
dtF<-rbindlist(lapply(names(FSTL),function(k) data.table(DI=mids,v=FSTL[[k]],src=LAB[k])))
dtL<-rbindlist(lapply(names(LDL),function(k) data.table(DI=mids,v=LDL[[k]],src=LAB[k])))
cols<-c("T1=0 (independent)"="#7fbf7b","T1=25 (half-shared)"="#1b9e77","T1=45 (near-panmictic)"="#0b6b4f")
pF<-ggplot()+geom_line(data=dtF,aes(DI,v,colour=src),linewidth=.8)+geom_point(data=dtF,aes(DI,v,colour=src),size=1.2)+
  geom_line(data=data.table(DI=mids,v=emp_fst),aes(DI,v),colour="#d95f02",linewidth=1.1)+
  geom_point(data=data.table(DI=mids,v=emp_fst),aes(DI,v),colour="#d95f02",size=1.8)+
  scale_colour_manual(values=cols,name=NULL)+
  labs(x="DiagnosticIndex",y=expression(F[ST]),title="A  F_ST vs DI: empirical (orange) vs shared-source nulls")+
  theme_bw(base_size=10)+theme(legend.position="bottom")
pL<-ggplot()+geom_line(data=dtL,aes(DI,v,colour=src),linewidth=.8)+geom_point(data=dtL,aes(DI,v,colour=src),size=1.2)+
  geom_line(data=data.table(DI=mids,v=emp_ld),aes(DI,v),colour="#d95f02",linewidth=1.1)+
  geom_point(data=data.table(DI=mids,v=emp_ld),aes(DI,v),colour="#d95f02",size=1.8)+
  scale_colour_manual(values=cols,name=NULL)+
  labs(x="DiagnosticIndex",y="adjacent-SNP LD",title="B  adjacent-LD vs DI: empirical (orange) vs shared-source nulls")+
  theme_bw(base_size=10)+theme(legend.position="bottom")
tab<-data.table(`split time`=c(LAB[as.character(T1S)],"EMPIRICAL"),
                `diag F_ST/LD ratio`=round(c(ratios,emp_ratio),1),
                `neutral->diag F_ST rise`=round(c(rises,emp_rise),3))
pT<-gridExtra::tableGrob(tab,rows=NULL,theme=gridExtra::ttheme_minimal(base_size=9))
comp<-(pF|pL)/pT + plot_layout(heights=c(3,1))+
  plot_annotation(title="Shared-source null: does shared ancestry move the diagnostic F_ST/LD dissociation toward empirical?",
                  subtitle="NREP=1 proof; 20 sub-demes downsampled to matched empirical population sizes",
                  theme=theme(plot.title=element_text(face="bold",size=12)))
ggsave(file.path(F,"Figures/moduleE_sharedsource.pdf"),comp,width=12,height=8)
ggsave(file.path(F,"Figures/moduleE_sharedsource.png"),comp,width=12,height=8,dpi=140)
cat("saved figure + rds\n")
