## =============================================================================
## Module E -- LD-PERMUTATION FLOOR: the irreducible between-population LD that
## recombination CANNOT remove (tests whether ANY recombinant-source model can work)
## =============================================================================
## Recombination removes within-population haplotype LD but preserves every population's
## per-locus allele frequency -- hence the BETWEEN-population covariance across adjacent
## loci. Pooled LD therefore has a floor = LD_between. We estimate it by permuting
## genotypes among individuals independently at each marker, WITHIN each population:
## this destroys within-population LD while holding every population allele frequency
## (and F_ST) fixed. The resulting pooled adjacent-SNP LD is LD_pooled(inf) ~ LD_between.
##
##   * sim FLOOR > empirical observed LD  -> no amount of recombination rescues that
##       demographic model (adjacent loci distinguish the same populations too consistently)
##   * empirical observed ~ empirical FLOOR -> empirical within-pop LD ~0 (all between-pop)
##   * gap (observed - floor) = the reducible within-population LD recombination could remove
## Compares empirical vs the current fixed-alpha neutral null, matched sample sizes.
## Output: console table + Figures/moduleE_ld_floor.{pdf,png} + data/moduleE_sim/moduleE_ld_floor.rds
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2) })
F<-"/Users/petrikem/gitlab/formica_hybrid"
FXDIR<-file.path(F,"moduleE_slim/output_replicates/rep01")               # current null (fixed-alpha, gen60)
FVCF <-file.path(F,"moduleE_slim/founders/maf015_DIstrat1500_burnin")
TAG<-"Naq450_Npol195"; K<-6250; GEN<-60L; NPERM<-5L
b<-readRDS(file.path(F,"moduleE_slim/inputs/empirical_bundle.rds")); uni<-as.data.table(b$universe); setkey(uni,marker)
fm<-rbindlist(lapply(list.files(FVCF,"founders_ch.*[.]vcf$",full.names=TRUE),
   function(f) fread(f,skip="#CHROM",select=c(1,2,3),header=TRUE,col.names=c("chs","Pos","marker"))))
fm[,`:=`(Chr=sub("ch","Chr",chs),chr_id=as.integer(sub("ch","",chs)))]; setorder(fm,chr_id,Pos)
markers<-fm$marker; DI<-uni[markers,DI]; emp_ns<-b$emp_ns
DIAG<-which(DI> -25)

read_dose<-function(f){dt<-fread(f,skip="#CHROM",header=TRUE,sep="\t",showProgress=FALSE);mk<-paste(sub("ch","Chr",dt[[1]]),dt[[2]],sep=":");G<-as.matrix(dt[,10:ncol(dt)]);gt<-sub(":.*$","",G);d<-matrix(as.integer(sub("[|/].*$","",gt))+as.integer(sub("^.*[|/]","",gt)),nrow=nrow(G));m<-match(markers,mk);keep<-!is.na(m);D<-matrix(0L,length(markers),ncol(d));D[keep,]<-d[m[keep],,drop=FALSE];t(D)}
## per-SNP mean r^2 to immediate chromosome neighbours, pooled genotype matrix G (indiv x markers)
adjld<-function(G){res<-rep(NA_real_,length(markers));for(cid in unique(fm$chr_id)){j<-which(fm$chr_id==cid);if(length(j)<2)next;for(a in seq_along(j)){nb<-j[c(a-1,a+1)];nb<-nb[nb>=min(j)&nb<=max(j)&nb!=j[a]];nb<-nb[!is.na(nb)];if(!length(nb))next;res[j[a]]<-mean(suppressWarnings(cor(G[,j[a]],G[,nb,drop=FALSE]))^2,na.rm=T)}};res}
diag_med<-function(v) median(v[DIAG],na.rm=TRUE)
## permute each marker within each population block -> destroys within-pop LD, keeps pop freqs
perm_within<-function(mats){do.call(rbind,lapply(mats,function(M){apply(M,2,function(col) col[sample.int(length(col))])}))}

## ---- assemble matched-size population matrices -------------------------------
## empirical
e<-new.env();load(file.path(F,"data/hybrids_only_maf005.Rdata"),envir=e);sdh<-as.data.table(e$sample_data);maph<-as.data.table(e$map_hyb_005)
hyb<-sdh[!Population%in%c("aquilonia_parent","polyctena_parent"),Sample_ID]
Ge<-e$GTs_hybrids_005[match(hyb,rownames(e$GTs_hybrids_005)),match(markers,maph$marker),drop=FALSE];Ge[is.na(Ge)]<-0L
popv<-sdh[match(hyb,Sample_ID),Population];rm(e);gc()
emp_mats<-lapply(split(seq_along(popv),popv),function(ix) Ge[ix,,drop=FALSE])
## sim (fixed-alpha), downsample each deme to matched empirical size
sim_mats<-lapply(1:20,function(i){D<-read_dose(file.path(FXDIR,sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf",TAG,K,i,GEN)));set.seed(400+i);D[sample(nrow(D),min(emp_ns[i],nrow(D))),,drop=FALSE]})

## ---- observed pooled adjacent LD + permutation floor ------------------------
floor_of<-function(mats){v<-numeric(NPERM);for(p in 1:NPERM){set.seed(1000+p);v[p]<-diag_med(adjld(perm_within(mats)))};v}
emp_obs<-diag_med(adjld(do.call(rbind,emp_mats))); emp_floor<-floor_of(emp_mats)
sim_obs<-diag_med(adjld(do.call(rbind,sim_mats))); sim_floor<-floor_of(sim_mats)

tab<-data.table(data=c("empirical","fixed-alpha null (current)"),
  observed_pooled_LD=c(emp_obs,sim_obs),
  floor_between_pop=c(mean(emp_floor),mean(sim_floor)),
  floor_sd=c(sd(emp_floor),sd(sim_floor)),
  reducible_within=c(emp_obs-mean(emp_floor),sim_obs-mean(sim_floor)))
cat("=== diagnostic adjacent-SNP LD: observed vs between-population floor ===\n"); print(tab)
cat(sprintf("\nKEY: sim floor (%.3f) vs empirical observed LD (%.3f) -> %s\n",
  mean(sim_floor),emp_obs, ifelse(mean(sim_floor)>emp_obs,
    "sim FLOOR ABOVE empirical LD: recombination CANNOT rescue this demographic model",
    "sim floor below empirical LD: solving for T_rec is meaningful")))
cat(sprintf("empirical observed (%.3f) vs empirical floor (%.3f): within-pop LD ~ %.3f (%.0f%% of observed)\n",
  emp_obs,mean(emp_floor),emp_obs-mean(emp_floor),100*(emp_obs-mean(emp_floor))/emp_obs))
## ---- single genome-wide ancestry axis, SCALED to empirical F_ST -------------
## What between-pop floor would a pure ancestry-axis model produce if it somehow reached
## the empirical diagnostic F_ST? p_jl = clamp(fpol_l + alpha*_j * delta_l); alpha*_j is
## the empirical alpha re-centred and scaled by s. Floor = pooled adjacent LD of genotypes
## drawn (independent loci, no within-pop LD) at matched sizes. If this floor >> empirical
## floor, adjacent loci distinguish the SAME populations too consistently -> a single axis
## cannot reach empirical F_ST with empirical (low) floor.
faq<-b$f_aq_par[markers]; fpol<-b$f_pol_par[markers]; delta<-faq-fpol
HI<-readRDS(file.path(F,"data/moduleE_sim/empirical_hybrid_index.rds"))$HI
dg<-which(DI> -25 & is.finite(faq)&is.finite(fpol)&abs(delta)>0.3)
fst_diag<-function(P){pb<-colMeans(P[,dg]);Ht<-2*pb*(1-pb);Hs<-colMeans(2*P[,dg]*(1-P[,dg]));ok<-Ht>0;sum((Ht-Hs)[ok])/sum(Ht[ok])}
axis_P<-function(s){a<-pmin(pmax(0.5+s*(HI-mean(HI)),0.001),0.999); P<-outer(a,delta); P<-sweep(P,2,fpol,`+`); pmin(pmax(P,0),1)}
## find scale s that hits empirical diagnostic F_ST
target<-fst_diag(do.call(rbind,lapply(emp_mats,function(M)matrix(colMeans(M)/2,1))))  # empirical diag F_ST
s<-uniroot(function(s) fst_diag(axis_P(s))-target, c(1,40))$root
Paxis<-axis_P(s)
## floor of this axis model: draw genotypes at matched sizes, independent loci, pool, adjld
ns<-sapply(emp_mats,nrow)
floor_from_freq<-function(P){v<-numeric(NPERM);for(p in 1:NPERM){set.seed(2000+p)
  mats<-lapply(1:20,function(j) matrix(rbinom(ns[j]*length(markers),2,rep(P[j,],each=ns[j])),nrow=ns[j],byrow=FALSE))
  v[p]<-diag_med(adjld(do.call(rbind,mats)))};v}
axis_floor<-floor_from_freq(Paxis)
cat(sprintf("\n=== single ancestry axis scaled to empirical diagnostic F_ST=%.2f (scale x%.1f) ===\n",target,s))
cat(sprintf("axis-model between-pop floor: %.3f  vs  empirical floor: %.3f  vs  empirical observed LD: %.3f\n",
  mean(axis_floor),mean(emp_floor),emp_obs))
cat(sprintf(" -> %s\n", ifelse(mean(axis_floor)>emp_obs,
  "single-axis floor EXCEEDS empirical observed LD: a genome-wide ancestry axis at empirical F_ST cannot reach the empirical dissociation (floor too high)",
  "single-axis floor below empirical observed LD")))

saveRDS(list(tab=tab,emp_floor=emp_floor,sim_floor=sim_floor,emp_obs=emp_obs,sim_obs=sim_obs,
             axis_floor=axis_floor,axis_scale=s,target_fst=target),
        file.path(F,"data/moduleE_sim/moduleE_ld_floor.rds"))

## ---- figure ------------------------------------------------------------------
dt<-rbindlist(list(
  data.table(x="empirical\n(observed)",     v=emp_obs,          kind="observed pooled LD"),
  data.table(x="empirical\nfloor",           v=mean(emp_floor),  kind="between-pop floor (recomb limit)"),
  data.table(x="fixed-alpha null\nfloor",    v=mean(sim_floor),  kind="between-pop floor (recomb limit)"),
  data.table(x="single ancestry axis\n@ F_ST=0.30 floor", v=mean(axis_floor), kind="between-pop floor (recomb limit)")))
dt[,x:=factor(x,levels=c("empirical\n(observed)","empirical\nfloor","fixed-alpha null\nfloor","single ancestry axis\n@ F_ST=0.30 floor"))]
g<-ggplot(dt,aes(x,v,fill=kind))+geom_col(width=.62)+
  geom_hline(yintercept=emp_obs,linetype=2,colour="#d95f02")+
  annotate("text",x=3.5,y=emp_obs,label="empirical observed LD (target)",vjust=-0.6,hjust=1,size=3,colour="#d95f02")+
  geom_text(aes(label=sprintf("%.3f",v)),vjust=-0.4,size=3)+
  scale_fill_manual(values=c("observed pooled LD"="#8c8c8c","between-pop floor (recomb limit)"="#1b6ca8"),name=NULL)+
  labs(x=NULL,y="diagnostic adjacent-SNP LD (median r^2)",
       title="LD-permutation floor: the between-population LD that recombination CANNOT remove",
       subtitle="A single genome-wide ancestry axis at empirical F_ST has a floor ABOVE the empirical total LD -> unreachable by any recombinant model")+
  theme_bw(base_size=11)+theme(legend.position="bottom")
ggsave(file.path(F,"Figures/moduleE_ld_floor.pdf"),g,width=9,height=5.2)
ggsave(file.path(F,"Figures/moduleE_ld_floor.png"),g,width=9,height=5.2,dpi=140)
cat("saved figure + rds\n")
