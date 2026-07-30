## =============================================================================
## Module E -- ld_w Manhattan: empirical vs finite-founder null vs BURN-IN null
## =============================================================================
## Does the per-species recombination burn-in (diversified founders) pull the NEUTRAL
## (low-DI) LD background down to empirical while leaving the DIAGNOSTIC admixture LD
## (the sorting signal) intact? Compares three pooled ld_w landscapes (gen 60, K6250,
## 450/195 founding, per-deme N matched to the data):
##   empirical | null (finite founders) | null (burn-in diversified)
## Outputs the Manhattan, a DI-split table, and neutral-locus ld_w by recombination.
##   Figures/moduleE_ldw_manhattan_burnin.{pdf,png}
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2) })
F        <- "/Users/petrikem/gitlab/formica_hybrid"
FVCF_DIR <- file.path(F, "moduleE_slim/founders/maf015_DIstrat1500_burnin")
BUNDLE   <- file.path(F, "moduleE_slim/inputs/empirical_bundle.rds")
EMP_RD   <- file.path(F, "data/hybrids_only_maf005.Rdata")
FIGDIR   <- file.path(F, "Figures")
GEN <- 60L; K <- 6250; TAG <- "Naq450_Npol195"; WIN <- 1e5L; MAXP <- 100L; MAF <- 0.1
SRCS <- list(list(lab="null (finite founders)",   dir=file.path(F,"moduleE_slim/output_foundersweep")),
             list(lab="null (burn-in diversified)",dir=file.path(F,"moduleE_slim/output_burnin")))

local_ld <- function(g, Chr, Pos) { res <- rep(NA_real_, ncol(g))
  for (ch in unique(Chr)) { j <- which(Chr==ch); j <- j[order(Pos[j])]; pos <- Pos[j]; G <- g[,j,drop=FALSE]
    p <- colMeans(G,na.rm=TRUE)/2; poly <- is.finite(p)&pmin(p,1-p)>=MAF
    for (a in seq_along(j)) { if(!poly[a]) next
      b <- which(pos<=pos[a]+WIN & seq_along(pos)>a & poly); if(length(b)>MAXP) b<-b[seq_len(MAXP)]; if(!length(b)) next
      res[j[a]] <- median(suppressWarnings(cor(G[,a],G[,b,drop=FALSE],use="pairwise.complete.obs"))^2,na.rm=TRUE) } }
  pmax(res,0) }
read_dose <- function(f) { dt<-fread(f,skip="#CHROM",header=TRUE,sep="\t",showProgress=FALSE)
  mk<-paste(sub("ch","Chr",dt[[1]]),dt[[2]],sep=":"); G<-as.matrix(dt[,10:ncol(dt)]); gt<-sub(":.*$","",G)
  d<-matrix(suppressWarnings(as.integer(sub("[|/].*$","",gt)))+suppressWarnings(as.integer(sub("^.*[|/]","",gt))),nrow=nrow(G))
  m<-match(markers,mk); keep<-!is.na(m); D<-matrix(0L,length(markers),ncol(d)); D[keep,]<-d[m[keep],,drop=FALSE]; t(D) }

fm <- rbindlist(lapply(list.files(FVCF_DIR,"founders_ch.*[.]vcf$",full.names=TRUE),
  function(f) fread(f,skip="#CHROM",select=c(1,2,3),header=TRUE,col.names=c("chs","Pos","marker"))))
fm[,`:=`(Chr=sub("ch","Chr",chs),chr_id=as.integer(sub("ch","",chs)))]; setorder(fm,chr_id,Pos); markers<-fm$marker
b <- readRDS(BUNDLE); uni<-as.data.table(b$universe); setkey(uni,marker); emp_ns<-b$emp_ns
DI <- uni[markers, DI]
## local recombination rate per marker (cM/Mb)
fm[, rate := NA_real_]
for (cid in unique(fm$chr_id)) { rc<-fread(file.path(F,"moduleE_slim/rec_maps",sprintf("ch_%d.recmap",cid)),header=FALSE,col.names=c("pos","r"))
  rc<-rc[is.finite(pos)&is.finite(r)][order(pos)]; fm[chr_id==cid, rate := approxfun(rc$pos,rc$r,rule=2)(Pos)] }

message("[1] empirical ...")
e<-new.env(); load(EMP_RD,envir=e); GTh<-e$GTs_hybrids_005; sdh<-as.data.table(e$sample_data); maph<-as.data.table(e$map_hyb_005)
hyb<-which(!sdh$Population%in%c("aquilonia_parent","polyctena_parent"))
eg<-GTh[hyb,match(markers,maph$marker),drop=FALSE]; eg[is.na(eg)]<-0L; rm(e); gc()
dt <- data.table(marker=markers, ld_w=local_ld(eg,fm$Chr,fm$Pos), src="empirical")
for (s in SRCS) { message("[2] ",s$lab," ...")
  parts<-lapply(1:20,function(i){f<-file.path(s$dir,sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf",TAG,K,i,GEN))
    if(!file.exists(f))return(NULL); g<-read_dose(f); tgt<-emp_ns[(i-1)%%length(emp_ns)+1]; set.seed(100+i); if(nrow(g)>tgt)g<-g[sample(nrow(g),tgt),,drop=FALSE]; g})
  dt<-rbind(dt, data.table(marker=markers, ld_w=local_ld(do.call(rbind,parts),fm$Chr,fm$Pos), src=s$lab)) }
dt <- fm[,.(marker,Chr,Pos,chr_id,rate)][dt,on="marker"]; dt[, DI := uni[marker,DI]]
saveRDS(dt, file.path(F,"data/moduleE_sim/moduleE_ldw_manhattan_burnin.rds"))

## ---- DI split + neutral-by-recombination summary ---------------------------
dt[, dic := fifelse(DI<=-90,"neutral (control)", fifelse(DI> -25,"diagnostic (signal)","mid"))]
cat("\n=== median ld_w by DI class ===\n")
print(dcast(dt[,.(ldw=round(median(ld_w,na.rm=TRUE),4)),by=.(src,dic)], src~dic, value.var="ldw"))
neu <- dt[DI<=-90 & is.finite(rate)]; qs<-quantile(unique(neu[,.(marker,rate)])$rate,c(0,1/3,2/3,1),na.rm=TRUE)
neu[, recbin:=cut(rate,qs,include.lowest=TRUE,labels=c("low-recomb","mid","high-recomb"))]
cat("\n=== NEUTRAL-locus median ld_w by local recombination ===\n")
print(dcast(neu[,.(ldw=round(median(ld_w,na.rm=TRUE),4)),by=.(recbin,src)], recbin~src, value.var="ldw"))

## ---- Manhattan (loess trend per chromosome) --------------------------------
lv <- c("empirical", sapply(SRCS,`[[`,"lab")); dt[, src:=factor(src,levels=lv)]
tr <- dt[is.finite(ld_w), { if(.N>=10){f<-pmax(0,predict(loess(ld_w~Pos,span=0.2)));.(Pos=Pos,y=f)} else .(Pos=Pos,y=ld_w) }, by=.(chr_id,Chr,src)]
tr[, Chr:=factor(Chr,levels=paste0("Chr",sort(unique(fm$chr_id))))]
cols<-c("empirical"="#d95f02","null (finite founders)"="#7570b3","null (burn-in diversified)"="#1b9e77")
p<-ggplot(tr,aes(Pos/1e6,y,colour=src))+geom_line(linewidth=.5)+
  facet_grid(.~Chr,scales="free_x",space="free_x")+scale_colour_manual(values=cols,name=NULL)+
  labs(x="position (Mbp)",y="local LD (median r^2, 100kb window)",
       title="ld_w landscape: empirical vs finite-founder null vs burn-in-diversified null (K6250, gen60, pooled)",
       subtitle="burn-in relaxes within-species founder LD -> neutral background drops toward empirical; diagnostic admixture LD (signal) retained")+
  theme_bw(base_size=8)+theme(panel.grid.minor=element_blank(),panel.spacing=unit(1,"pt"),
    strip.text.x=element_text(size=5,angle=90),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="top")
ggsave(file.path(FIGDIR,"moduleE_ldw_manhattan_burnin.pdf"),p,width=14,height=4.4)
ggsave(file.path(FIGDIR,"moduleE_ldw_manhattan_burnin.png"),p,width=14,height=4.4,dpi=130)
cat("\nsaved: Figures/moduleE_ldw_manhattan_burnin.{pdf,png}\n")
