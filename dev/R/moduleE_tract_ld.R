## =============================================================================
## Module E -- tract length, recombination time, and the source of pooled LD
## =============================================================================
## Four-panel synthesis (colleague's proposal, adapted for the F_ST confound: no neutral
## sim reaches the empirical F_ST, so the "genome-wide ancestry" comparator is a single
## ancestry axis SCALED to empirical F_ST=0.30, not the raw low-F_ST null).
##   A  ancestry-LD decay vs genetic distance: empirical vs the gen-60 null (tract age).
##   B  the same decay decomposed into within-population LD (tract-driven, recombination-
##      removable) and the between-population floor (covariance, irreducible).
##   C  DECISIVE: distribution of the 20-population allele-frequency-vector correlation of
##      ADJACENT diagnostic loci -- empirical vs the single-axis genome-wide model.
##   D  within- vs between-population contribution to pooled adjacent LD: empirical,
##      single-axis@F_ST, gen-60 null. The recombination-aging limits (t->inf) are the
##      between-pop floors; empirical (0.127) sits between the null floor (0.010, too weak)
##      and the single-axis floor (0.137, too correlated).
## Output: Figures/moduleE_tract_ld.{pdf,png} + data/moduleE_sim/moduleE_tract_ld.rds
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
F<-"/Users/petrikem/gitlab/formica_hybrid"
FVCF<-file.path(F,"moduleE_slim/founders/maf015_DIstrat1500_burnin"); FXDIR<-file.path(F,"moduleE_slim/output_replicates/rep01"); RECDIR<-file.path(F,"moduleE_slim/rec_maps")
TAG<-"Naq450_Npol195"; K<-6250; GEN<-60L
b<-readRDS(file.path(F,"moduleE_slim/inputs/empirical_bundle.rds")); uni<-as.data.table(b$universe); setkey(uni,marker)
fm<-rbindlist(lapply(list.files(FVCF,"founders_ch.*[.]vcf$",full.names=TRUE),
   function(f) fread(f,skip="#CHROM",select=c(1,2,3),header=TRUE,col.names=c("chs","Pos","marker"))))
fm[,`:=`(chr_id=as.integer(sub("ch","",chs)))]; setorder(fm,chr_id,Pos)
markers<-fm$marker; DI<-uni[markers,DI]; faq<-b$f_aq_par[markers]; fpol<-b$f_pol_par[markers]; delta<-faq-fpol; emp_ns<-b$emp_ns
HI<-readRDS(file.path(F,"data/moduleE_sim/empirical_hybrid_index.rds"))$HI
## cM position per marker
cM<-rep(NA_real_,nrow(fm))
for(id in unique(fm$chr_id)){rc<-fread(file.path(RECDIR,sprintf("ch_%d.recmap",id)),header=FALSE,col.names=c("pos","rate"));rc<-rc[is.finite(pos)&is.finite(rate)][order(pos)]
  j<-which(fm$chr_id==id); cM[j]<-approxfun(rc$pos,c(0,cumsum(rc$rate[-nrow(rc)]*diff(rc$pos)/1e6)),rule=2)(fm$Pos[j])}
fm[,cM:=cM]

read_dose<-function(f){dt<-fread(f,skip="#CHROM",header=TRUE,sep="\t",showProgress=FALSE);mk<-paste(sub("ch","Chr",dt[[1]]),dt[[2]],sep=":");G<-as.matrix(dt[,10:ncol(dt)]);gt<-sub(":.*$","",G);d<-matrix(as.integer(sub("[|/].*$","",gt))+as.integer(sub("^.*[|/]","",gt)),nrow=nrow(G));m<-match(markers,mk);keep<-!is.na(m);D<-matrix(0L,length(markers),ncol(d));D[keep,]<-d[m[keep],,drop=FALSE];t(D)}

## ---- per-population genotype matrices (matched sizes) ------------------------
e<-new.env();load(file.path(F,"data/hybrids_only_maf005.Rdata"),envir=e);sdh<-as.data.table(e$sample_data);maph<-as.data.table(e$map_hyb_005)
hyb<-sdh[!Population%in%c("aquilonia_parent","polyctena_parent"),Sample_ID]
Ge<-e$GTs_hybrids_005[match(hyb,rownames(e$GTs_hybrids_005)),match(markers,maph$marker),drop=FALSE];Ge[is.na(Ge)]<-0L
popv<-sdh[match(hyb,Sample_ID),Population];rm(e);gc()
emp_mats<-lapply(split(seq_along(popv),popv),function(ix) Ge[ix,,drop=FALSE])
sim_mats<-lapply(1:20,function(i){D<-read_dose(file.path(FXDIR,sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf",TAG,K,i,GEN)));set.seed(400+i);D[sample(nrow(D),min(emp_ns[i],nrow(D))),,drop=FALSE]})
Gemp<-do.call(rbind,emp_mats); Gsim<-do.call(rbind,sim_mats)
perm_within<-function(mats){do.call(rbind,lapply(mats,function(M) apply(M,2,function(col) col[sample.int(length(col))])))}
set.seed(11); Gemp_f<-perm_within(emp_mats); Gsim_f<-perm_within(sim_mats)

## ---- Panels A/B: LD (r^2) vs genetic distance for diagnostic pairs ----------
finite<-apply(do.call(rbind,lapply(emp_mats,function(M)colMeans(M)/2)),2,function(c)all(is.finite(c)))
dgi<-which(DI> -25 & abs(delta)>0.3 & finite)                        # diagnostic marker indices
maxcM<-1.5; anchors<-dgi[seq(1,length(dgi),by=3)]                    # subsample anchors for speed
r2_decay<-function(G){res<-list();for(a in anchors){cid<-fm$chr_id[a];nb<-dgi[fm$chr_id[dgi]==cid & fm$cM[dgi]>fm$cM[a] & fm$cM[dgi]-fm$cM[a]<=maxcM]
    if(!length(nb))next;r2<-suppressWarnings(cor(G[,a],G[,nb,drop=FALSE]))^2;res[[length(res)+1]]<-data.table(d=fm$cM[nb]-fm$cM[a],r2=as.numeric(r2))};rbindlist(res)}
bincurve<-function(dt){dt[,db:=cut(d,seq(0,maxcM,0.1))][!is.na(db),.(d=mean(d),r2=mean(r2,na.rm=T)),by=db]}
A_emp<-bincurve(r2_decay(Gemp)); A_sim<-bincurve(r2_decay(Gsim))
A_empf<-bincurve(r2_decay(Gemp_f)); A_simf<-bincurve(r2_decay(Gsim_f))   # between-pop floor curves

## ---- Panel C: adjacent-loci 20-pop freq-vector correlation ------------------
Pemp<-do.call(rbind,lapply(emp_mats,function(M)colMeans(M)/2)); Psim<-do.call(rbind,lapply(sim_mats,function(M)colMeans(M)/2))
fl<-which(delta<0 & !is.na(delta)); Pemp[,fl]<-1-Pemp[,fl]; Psim[,fl]<-1-Psim[,fl]     # orient to aq allele
ad<-abs(delta); Paxis<-function(s){a<-pmin(pmax(0.5+s*(HI-mean(HI)),0.001),0.999); pmin(pmax(sweep(outer(a,ad),2,pmin(fpol,faq),`+`),0),1)}
fst<-function(P,ix){pb<-colMeans(P[,ix]);Ht<-2*pb*(1-pb);Hs<-colMeans(2*P[,ix]*(1-P[,ix]));ok<-Ht>0;sum((Ht-Hs)[ok])/sum(Ht[ok])}
s<-uniroot(function(s) fst(Paxis(s),dgi)-fst(Pemp,dgi),c(1,40))$root; Pax<-Paxis(s)
adjpairs<-do.call(rbind,lapply(unique(fm$chr_id),function(cid){j<-which(fm$chr_id==cid & DI> -25 & abs(delta)>0.3 & finite);if(length(j)<2)return(NULL);cbind(j[-length(j)],j[-1])}))
cp<-function(P) apply(adjpairs,1,function(p){v1<-P[,p[1]];v2<-P[,p[2]];if(sd(v1)<1e-9||sd(v2)<1e-9)NA else cor(v1,v2)})
cE<-cp(Pemp); cA<-cp(Pax); cS<-cp(Psim)

## ---- Panel D: within vs between pooled adjacent LD --------------------------
adj<-function(G){res<-rep(NA_real_,length(markers));for(cid in unique(fm$chr_id)){j<-which(fm$chr_id==cid);if(length(j)<2)next;for(a in seq_along(j)){nb<-j[c(a-1,a+1)];nb<-nb[nb>=min(j)&nb<=max(j)&nb!=j[a]];nb<-nb[!is.na(nb)];if(!length(nb))next;res[j[a]]<-mean(suppressWarnings(cor(G[,j[a]],G[,nb,drop=FALSE]))^2,na.rm=T)}};res}
md<-function(v) median(v[dgi],na.rm=TRUE)
axis_geno<-function(P){ns<-sapply(emp_mats,nrow);do.call(rbind,lapply(1:20,function(j) matrix(rbinom(ns[j]*length(markers),2,rep(P[j,],each=ns[j])),nrow=ns[j])))}
set.seed(7); Gax<-axis_geno(Pax); Gax_f<-perm_within(lapply(split(1:nrow(Gax),rep(1:20,sapply(emp_mats,nrow))),function(ix)Gax[ix,,drop=FALSE]))
D<-data.table(model=c("empirical","single axis @F_ST","gen-60 null"),
  observed=c(md(adj(Gemp)),md(adj(Gax)),md(adj(Gsim))),
  between =c(md(adj(Gemp_f)),md(adj(Gax_f)),md(adj(Gsim_f))))
D[,within:=observed-between]
saveRDS(list(A_emp=A_emp,A_sim=A_sim,A_empf=A_empf,A_simf=A_simf,cE=cE,cA=cA,cS=cS,D=D,axis_scale=s),file.path(F,"data/moduleE_sim/moduleE_tract_ld.rds"))
cat("=== panel D (within/between pooled adjacent LD) ===\n"); print(D)
cat(sprintf("\npanel C adjacent freq-vector corr: emp median=%.2f  axis median=%.2f  null median=%.2f\n",median(cE,na.rm=T),median(cA,na.rm=T),median(cS,na.rm=T)))

## ---- figure ------------------------------------------------------------------
cols<-c("empirical"="#d95f02","gen-60 null"="#1b9e77","single axis @F_ST"="#1b6ca8")
pA<-ggplot()+geom_line(data=A_emp,aes(d,r2,colour="empirical"),linewidth=1)+geom_line(data=A_sim,aes(d,r2,colour="gen-60 null"),linewidth=1)+
  scale_colour_manual(values=cols,name=NULL)+labs(x="genetic distance (cM)",y=expression(r^2~"(diagnostic pairs)"),
  title="A  ancestry-LD decay vs distance\n(null tracts decay slower = younger)")+theme_bw(base_size=9)+theme(legend.position="bottom")
pB<-ggplot()+geom_line(data=A_emp,aes(d,r2,colour="empirical: total"),linewidth=.9)+geom_line(data=A_empf,aes(d,r2,colour="empirical: between-pop floor"),linewidth=.9,linetype=2)+
  geom_line(data=A_sim,aes(d,r2,colour="null: total"),linewidth=.9)+geom_line(data=A_simf,aes(d,r2,colour="null: between-pop floor"),linewidth=.9,linetype=2)+
  scale_colour_manual(values=c("empirical: total"="#d95f02","empirical: between-pop floor"="#f0a060","null: total"="#1b9e77","null: between-pop floor"="#8fd0b8"),name=NULL)+
  labs(x="genetic distance (cM)",y=expression(r^2),title="B  total vs between-pop floor\n(within-pop LD = tract-driven, removable)")+theme_bw(base_size=9)+theme(legend.position="bottom",legend.text=element_text(size=6))
dc<-rbindlist(list(data.table(r=cE,src="empirical"),data.table(r=cS,src="gen-60 null"),data.table(r=cA,src="single axis @F_ST")))
pC<-ggplot(dc[is.finite(r)],aes(r,colour=src,fill=src))+geom_density(alpha=.15,linewidth=.9,bw=.05)+
  scale_colour_manual(values=cols,name=NULL)+scale_fill_manual(values=cols,name=NULL)+
  labs(x="20-pop freq-vector correlation of adjacent diagnostic loci",y="density",
  title=sprintf("C  DECISIVE: do adjacent loci sort the SAME pops?\nmedian r: empirical %.2f < null %.2f < single-axis 1.00",median(cE,na.rm=T),median(cS,na.rm=T)))+
  theme_bw(base_size=9)+theme(legend.position="bottom")
Dm<-melt(D[,.(model,within,between)],id.vars="model",variable.name="component",value.name="LD")
Dm[,model:=factor(model,levels=c("gen-60 null","empirical","single axis @F_ST"))]
pD<-ggplot(Dm,aes(model,LD,fill=component))+geom_col(width=.6)+geom_hline(yintercept=D[model=="empirical",observed],linetype=2,colour="#d95f02")+
  scale_fill_manual(values=c("within"="#8c8c8c","between"="#1b6ca8"),name=NULL)+
  labs(x=NULL,y="pooled adjacent LD (median r^2)",title="D  within vs between-pop LD\n(recomb removes within; between is the floor)")+
  theme_bw(base_size=9)+theme(legend.position="bottom",axis.text.x=element_text(angle=15,hjust=1))
comp<-(pA|pB)/(pC|pD)+plot_annotation(title="Tract age is not the missing ingredient: empirical sorts adjacent loci along different population axes",
  subtitle="null too young AND too weak in F_ST; a genome-wide axis at empirical F_ST over-correlates adjacent loci (floor above empirical). Empirical is locus-specific.",
  theme=theme(plot.title=element_text(face="bold",size=11)))
ggsave(file.path(F,"Figures/moduleE_tract_ld.pdf"),comp,width=12,height=9)
ggsave(file.path(F,"Figures/moduleE_tract_ld.png"),comp,width=12,height=9,dpi=140)
cat("saved figure + rds\n")
