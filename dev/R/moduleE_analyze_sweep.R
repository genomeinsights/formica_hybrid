## =============================================================================
## Module E -- analyze the parameter sweep (output_track_length): diversity, F_ST,
## adjacent-LD vs time, per (K, founders) combination. Preliminary Ne + tract-clock.
## =============================================================================
## The sweep = 5 combinations (base K=6250/f150; K in {3125,6250,12500}; total founders
## in {100,150,200}), 20 demes each, sampled every 25 ticks to 250. Markers = 9,902
## high-DI (DI>-25) loci. For each combination x tick we compute:
##   - He_within  : mean within-deme expected heterozygosity (drift clock -> Ne)
##   - Fst_among  : among-deme F_ST at diagnostic loci (rises with drift)
##   - LD_adj     : mean within-deme adjacent-marker r^2 (admixture LD; falls with
##                  recombination = ancestry-tract breaking -> tract clock)
## Ticks -> generations via the model generation time ~3.4 ticks. Empirical values on
## the SAME markers are overlaid as reference lines.
## Output: Figures/moduleE_sweep_summary.{pdf,png} + data/moduleE_sim/moduleE_sweep.rds
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
F<-"/Users/petrikem/gitlab/formica_hybrid"; SW<-file.path(F,"moduleE_slim/output_track_length")
MK<-fread("/Users/petrikem/gitlab/formica_neutral_sim_sweep_handoff/markers.tsv")
MK[,chr_id:=as.integer(sub("Chr","",Chr))]; setorder(MK,chr_id,Pos)
markers<-MK$marker; TGEN<-3.4
combos<-list(K3125_f150=list(K=3125,naq=105,npol=45), K6250_f150=list(K=6250,naq=105,npol=45),
             K12500_f150=list(K=12500,naq=105,npol=45), K6250_f100=list(K=6250,naq=70,npol=30),
             K6250_f200=list(K=6250,naq=140,npol=60))
## fix: filenames carry the exact Naq/Npol; read them from disk instead of hardcoding
gens<-c(25,50,75,100,125,150,175,200,225,250)

read_dose<-function(f){dt<-fread(f,skip="#CHROM",header=TRUE,sep="\t",showProgress=FALSE);mk<-paste(sub("ch","Chr",dt[[1]]),dt[[2]],sep=":");G<-as.matrix(dt[,10:ncol(dt)]);gt<-sub(":.*$","",G);d<-matrix(as.integer(sub("[|/].*$","",gt))+as.integer(sub("^.*[|/]","",gt)),nrow=nrow(G));m<-match(markers,mk);keep<-!is.na(m);D<-matrix(0L,length(markers),ncol(d));D[keep,]<-d[m[keep],,drop=FALSE];t(D)}
nb_idx<-local({ res<-vector("list",nrow(MK)); for(cid in unique(MK$chr_id)){j<-which(MK$chr_id==cid); for(a in seq_along(j)){nb<-j[c(a-1,a+1)]; nb<-nb[nb>=min(j)&nb<=max(j)&nb!=j[a]]; res[[j[a]]]<-nb[!is.na(nb)]}}; res })
adj_r2<-function(G){ v<-rep(NA_real_,ncol(G)); for(a in seq_len(ncol(G))){ nb<-nb_idx[[a]]; if(!length(nb))next; v[a]<-mean(suppressWarnings(cor(G[,a],G[,nb,drop=FALSE]))^2,na.rm=TRUE)}; median(v,na.rm=TRUE) }

res<-list()
for(cl in names(combos)){ dir<-file.path(SW,cl)
  for(g in gens){ files<-list.files(dir,sprintf("_rep[0-9]+_gen%d\\.vcf$",g),full.names=TRUE)
    if(!length(files))next
    Ds<-lapply(files,read_dose)
    Fr<-t(vapply(Ds,function(D)colMeans(D)/2,numeric(length(markers))))               # deme x marker freqs
    He_within<-mean(rowMeans(2*Fr*(1-Fr),na.rm=TRUE))                                  # mean within-deme He
    pb<-colMeans(Fr); Ht<-2*pb*(1-pb); Hs<-colMeans(2*Fr*(1-Fr)); ok<-Ht>0; Fst<-sum((Ht-Hs)[ok])/sum(Ht[ok])
    ld<-mean(vapply(Ds,adj_r2,numeric(1)),na.rm=TRUE)                                  # within-deme adjacent LD
    res[[length(res)+1]]<-data.table(combo=cl,K=combos[[cl]]$K,
      founders=combos[[cl]]$naq+combos[[cl]]$npol, tick=g, gen=g/TGEN,
      He=He_within, Fst=Fst, LDadj=ld, ndeme=length(files))
    cat(sprintf("%-12s gen %3d (%.0f gens): He=%.4f Fst=%.4f LDadj=%.3f (%d demes)\n",cl,g,g/TGEN,He_within,Fst,ld,length(files)))
  }
}
S<-rbindlist(res)
## rough Ne per combination from within-deme He decay: log He ~ a - gen/(2Ne)
Ne<-S[,{fit<-lm(log(He)~gen); .(Ne=round(-1/(2*coef(fit)[2])), r2=round(summary(fit)$r.squared,3))},by=.(combo,K,founders)]
cat("\n=== preliminary Ne (from within-deme He decay, generations = tick/3.4) ===\n"); print(Ne)
saveRDS(list(S=S,Ne=Ne),file.path(F,"data/moduleE_sim/moduleE_sweep.rds"))

## empirical reference on the SAME 9902 markers -------------------------------
b<-readRDS(file.path(F,"moduleE_slim/inputs/empirical_bundle.rds"))
e<-new.env();load(file.path(F,"data/hybrids_only_maf005.Rdata"),envir=e);sdh<-as.data.table(e$sample_data);maph<-as.data.table(e$map_hyb_005)
hyb<-sdh[!Population%in%c("aquilonia_parent","polyctena_parent"),Sample_ID]
mi<-match(markers,maph$marker); have<-!is.na(mi)
Ge<-matrix(0L,length(hyb),length(markers)); Ge[,have]<-e$GTs_hybrids_005[match(hyb,rownames(e$GTs_hybrids_005)),mi[have],drop=FALSE]; Ge[is.na(Ge)]<-0L
popv<-sdh[match(hyb,Sample_ID),Population];rm(e);gc()
emp_mats<-lapply(split(seq_along(popv),popv),function(ix) Ge[ix,,drop=FALSE])
Frp<-t(vapply(emp_mats,function(M)colMeans(M)/2,numeric(length(markers))))
pb<-colMeans(Frp);Ht<-2*pb*(1-pb);Hs<-colMeans(2*Frp*(1-Frp));ok<-Ht>0
emp_fst<-sum((Ht-Hs)[ok])/sum(Ht[ok]); emp_He<-mean(rowMeans(2*Frp*(1-Frp)))
emp_ld<-mean(vapply(emp_mats,function(M) adj_r2(M),numeric(1)),na.rm=TRUE)
cat(sprintf("\nempirical on these %d markers (%d present): He=%.4f Fst=%.4f LDadj=%.3f\n",length(markers),sum(have),emp_He,emp_fst,emp_ld))

## figure ---------------------------------------------------------------------
S[,combo:=factor(combo,levels=names(combos))]
pal<-c(K3125_f150="#66c2a5",K6250_f150="#1b9e77",K12500_f150="#0b6b4f",K6250_f100="#fdae61",K6250_f200="#d95f02")
base<-function(y,lab,eref) ggplot(S,aes(gen,get(y),colour=combo))+geom_line(linewidth=.8)+geom_point(size=1.3)+
  geom_hline(yintercept=eref,linetype=2,colour="grey30")+annotate("text",x=max(S$gen),y=eref,label="empirical",hjust=1,vjust=-0.4,size=3,colour="grey30")+
  scale_colour_manual(values=pal,name=NULL)+labs(x="generations (= ticks / 3.4)",y=lab)+theme_bw(base_size=10)+theme(legend.position="bottom")
pHe<-base("He",  "within-deme He (diversity)", emp_He)+ggtitle("A  diversity decay -> Ne")
pF <-base("Fst", expression(among-deme~F[ST]), emp_fst)+ggtitle("B  among-deme F_ST")
pL <-base("LDadj","within-deme adjacent-marker r^2", emp_ld)+ggtitle("C  adjacent LD (tract clock)")
comp<-(pHe|pF|pL)+plot_layout(guides="collect")&theme(legend.position="bottom")
comp<-comp+plot_annotation(title="Parameter sweep (output_track_length): diversity, differentiation, adjacent-LD vs time",
  subtitle="9,902 high-DI markers; 20 demes/combination; empirical (dashed) on the same markers",theme=theme(plot.title=element_text(face="bold",size=12)))
ggsave(file.path(F,"Figures/moduleE_sweep_summary.pdf"),comp,width=13,height=5)
ggsave(file.path(F,"Figures/moduleE_sweep_summary.png"),comp,width=13,height=5,dpi=140)
cat("saved figure + rds\n")
