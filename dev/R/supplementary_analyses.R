## ==========================================================================
## Supplementary analyses: ancestry sorting and its predictability across
## 20 replicate Formica polyctena x F. aquilonia hybrid populations
## ==========================================================================
## Reproduces all statistics and figures in the supplementary methods.
## Function library: dev/R/parallelism_stats.R (+ dev/R/Ohta.R for prep).
##
## Sections
##   0. Setup and data
##   1. Per-locus sorting statistics (parallelism_stats) + sort_class
##   2. Choice of null for the parallelism test (pooled null is circular)
##   3. Differentiation gate (DI vs parent_diff)
##   4. Tier 1: genomic architecture of differentiation
##   5. Tier 2: sorting vs recombination at the independent-unit level
##   6. LD structure for neutral-null design (within-species vs admixture LD)
## ==========================================================================

suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
source("dev/R/Ohta.R")               # ohta_fast_prepare()
source("dev/R/parallelism_stats.R")  # parallelism_stats(), sort_class, etc.
set.seed(1)
dir.create("Figures", showWarnings = FALSE)

## --------------------------------------------------------------------------
## 0. Setup and data
## --------------------------------------------------------------------------
e <- new.env(); load("./data/hybrids_and_parents_maf005.Rdata", envir = e)
GTs <- e$GTs_with_parents            # 194 samples x 1,114,340 SNPs, dosage 0/1/2
sd  <- e$sample_data_with_parents    # Population, Sample_ID, PC1, PC2
map <- copy(e$map_hyb_005)           # marker, Chr, Pos, Polarity, DiagnosticIndex, ...

aqu <- "aquilonia_parent"; pol <- "polyctena_parent"
hyb <- setdiff(unique(sd$Population), c(aqu, pol))   # 20 hybrid populations

## per-marker recombination rate (cM/Mb), interpolated from the linkage map
rec <- fread("./data/Frufa_DTOL_PR.ref_genome.recmap")
setnames(rec, c("chr","pos","cM","cMMb")); rec[, Chr := sub("chromosome_","Chr",chr)]
map[, recomb := NA_real_]
for (ch in unique(map$Chr)) {
  r <- rec[Chr == ch]; if (nrow(r) < 2) next
  idx <- map[, which(Chr == ch)]
  map[idx, recomb := approx(r$pos, r$cMMb, xout = map$Pos[idx], rule = 2)$y]
}
recomb_by <- setNames(map$recomb, map$marker)
DI_by     <- setNames(map$DiagnosticIndex, map$marker)

## --------------------------------------------------------------------------
## 1. Per-locus sorting statistics + sort_class classification
##    (diagnostic SNP set: DiagnosticIndex > -15)
## --------------------------------------------------------------------------
di15 <- which(map$DiagnosticIndex > -15)                 # 16,077 SNPs
prep15 <- ohta_fast_prepare(GTs[, di15], pops = sd$Population)
r15 <- parallelism_stats(prep15, hyb, aqu, pol, DI = DI_by,
                         min_DI = -15, sort_th = 0.5, fix_th = 0.1, null_prob = 0.5)

## sort_class proportions (unified threshold sort_th = 0.5)
sort_tab <- r15[differentiated == TRUE, .N, by = sort_class][order(-N)]
sort_tab[, pct := round(100 * N / sum(N), 1)]
print(sort_tab)                                          # 89.9% unsorted; aqu > pol

## --------------------------------------------------------------------------
## 2. Choice of null for the parallelism test: the per-locus "pooled" null
##    is circular (its null prob tracks the observed outcome), so we use the
##    symmetric null (justified because global mean ancestry ~ 0.50).
## --------------------------------------------------------------------------
r_sym  <- parallelism_stats(prep15, hyb, aqu, pol, DI = DI_by, min_DI = -15, null_prob = 0.5)
r_pool <- parallelism_stats(prep15, hyb, aqu, pol, DI = DI_by, min_DI = -15, null_prob = "pooled")
te <- !is.na(r_pool$p_binom)
cat("pooled null_p vs observed aqu-fraction, correlation:",
    round(cor(r_pool$null_p[te], (r_pool$n_aqu/r_pool$n_fixed)[te]), 3), "\n")
cat("significant (q<0.05)  symmetric:", sum(r_sym$q_binom < 0.05, na.rm = TRUE),
    " pooled:", sum(r_pool$q_binom < 0.05, na.rm = TRUE), "\n")
g_aqu <- mean(r_pool$f_aqu_pooled[r_pool$differentiated], na.rm = TRUE)
cat("global mean aquilonia ancestry (diagnostic loci):", round(g_aqu, 3), "\n")

## --------------------------------------------------------------------------
## 3. Differentiation gate: DI and empirical parental differentiation agree
## --------------------------------------------------------------------------
pa0 <- colMeans(GTs[sd$Population == aqu, ], na.rm = TRUE) / 2
pp0 <- colMeans(GTs[sd$Population == pol, ], na.rm = TRUE) / 2
map[, parent_diff := abs(pa0 - pp0)]
for (th in c(-15, -20, -25, -30)) {
  s <- map[DiagnosticIndex > th]
  cat(sprintf("DI > %d : frac parent_diff>=0.5 = %.3f (median %.2f)\n",
              th, mean(s$parent_diff >= 0.5, na.rm = TRUE), median(s$parent_diff, na.rm = TRUE)))
}

## --------------------------------------------------------------------------
## 4. Tier 1: genomic architecture of differentiation
##    (per-marker pi within species, dxy absolute, Fst relative, cluster size)
## --------------------------------------------------------------------------
g <- readRDS("data/eMLG_5loci_0025_cM05.rds")$groups
memb <- data.table(marker = unlist(g$members),
                   cluster_size = rep(g$n_loci, lengths(g$members)))
map <- memb[map, on = "marker"]

Hs_a <- 2*pa0*(1-pa0); Hs_p <- 2*pp0*(1-pp0)
map[, pi_within := (Hs_a + Hs_p)/2]                      # within-species diversity
map[, dxy       := pa0*(1-pp0) + pp0*(1-pa0)]            # absolute divergence
pbar <- (pa0+pp0)/2; Ht <- 2*pbar*(1-pbar)
map[, Fst := ifelse(Ht > 0, (Ht - (Hs_a+Hs_p)/2)/Ht, NA_real_)]  # relative differentiation

d1 <- map[is.finite(recomb) & is.finite(DiagnosticIndex) & !is.na(cluster_size)]
sp <- function(a,b) round(cor(a,b,method="spearman",use="complete.obs"),3)
cat("\nTier1 Spearman -- DI~recomb:", sp(d1$DiagnosticIndex,d1$recomb),
    " DI~cluster:", sp(d1$DiagnosticIndex,d1$cluster_size),
    " DI~pi:", sp(d1$DiagnosticIndex,d1$pi_within),
    " DI~dxy:", sp(d1$DiagnosticIndex,d1$dxy), "\n")

d1[, rbin := cut(recomb, quantile(recomb,0:10/10,na.rm=TRUE), include.lowest=TRUE, labels=FALSE)]
t1 <- d1[, .(med_recomb=round(median(recomb),1), DI=round(mean(DiagnosticIndex),1),
             Fst=round(mean(Fst,na.rm=TRUE),3), dxy=round(mean(dxy,na.rm=TRUE),3),
             pi_within=round(mean(pi_within,na.rm=TRUE),3),
             cluster_size=round(mean(cluster_size),0)), by=rbin][order(rbin)]
print(t1)

## Figure S1
z <- function(x)(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
th_bw <- theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),legend.position="bottom")
ch <- melt(t1[,.(med_recomb, DI=z(DI), dxy=z(dxy), pi_within=z(pi_within), Fst=z(Fst))],
           id.vars="med_recomb", variable.name="metric", value.name="zval")
pA <- ggplot(ch, aes(med_recomb, zval, color=metric))+geom_line(linewidth=1)+geom_point()+
  scale_x_log10()+labs(x="recombination (cM/Mb, log)", y="standardised (z)",
  title="(A) relative (DI/Fst) vs absolute (dxy) vs pi vs recombination")+th_bw
pB <- ggplot(t1, aes(med_recomb, cluster_size))+geom_line(linewidth=1,color="grey30")+geom_point()+
  scale_x_log10()+scale_y_log10()+labs(x="recombination (cM/Mb, log)", y="mean cluster size",
  title="(B) cluster size vs recombination")+th_bw
d1[, cbin := cut(cluster_size, unique(quantile(cluster_size,0:12/12)), include.lowest=TRUE)]
dc <- d1[!is.na(cbin), .(med_cs=median(cluster_size), DI=mean(DiagnosticIndex)), by=cbin][order(med_cs)]
pC <- ggplot(dc, aes(med_cs, DI))+geom_line(linewidth=1,color="firebrick")+geom_point()+
  scale_x_log10()+labs(x="cluster size (log)", y="mean DI", title="(C) DI vs cluster size")+th_bw
ggsave("Figures/tier1_architecture.png", (pA|pB)/(pC|plot_spacer()), width=13, height=9, dpi=110)

## --------------------------------------------------------------------------
## 5. Tier 2: sorting vs recombination at the independent-unit (eMLG) level
##    Units = LD-pruning representatives (one tag SNP per cluster).
## --------------------------------------------------------------------------
reps  <- intersect(g$representative, colnames(GTs))
cs_by <- setNames(g$n_loci, g$representative)

run_ps <- function(markers) {
  p <- ohta_fast_prepare(GTs[, markers], pops = sd$Population)
  r <- parallelism_stats(p, hyb, aqu, pol, DI = DI_by,
                         min_parent_diff = 0.5, sort_th = 0.5, fix_th = 0.1)
  r[, recomb := recomb_by[marker]][]
}
r_unit <- run_ps(reps);            r_unit[, n_loci := as.double(cs_by[marker])]
r_snp  <- run_ps(sample(colnames(GTs), 200000))

brk <- quantile(map$recomb, 0:10/10, na.rm=TRUE)
binit <- function(d){ d[, rbin := cut(recomb, brk, include.lowest=TRUE, labels=FALSE)][] }
r_unit <- binit(r_unit); r_snp <- binit(r_snp)
srt <- function(sc) sc %in% c("aquilonia","polyctena","bidirectional")
summ <- function(d, lab) d[differentiated==TRUE, .(level=lab, med_r=median(recomb,na.rm=TRUE),
   frac_sorted=mean(srt(sort_class)), frac_uni=mean(sort_class%in%c("aquilonia","polyctena")),
   frac_bi=mean(sort_class=="bidirectional")), by=rbin][order(rbin)]
su <- summ(r_unit,"unit (eMLG)"); ss <- summ(r_snp,"SNP")
print(su)

## regressions (standardised): sorting ~ recombination + DI (separable), + collider
du <- r_unit[differentiated==TRUE & is.finite(recomb) & is.finite(DI)]
du[, `:=`(zr=scale(log10(recomb+0.1)), zDI=scale(DI), zcs=scale(log10(n_loci)))]
cat("\nTier2 prop_fixed ~ recomb + DI:\n");        print(round(coef(summary(lm(prop_fixed~zr+zDI,       du))),4))
cat("\nTier2 + cluster size (collider):\n");        print(round(coef(summary(lm(prop_fixed~zr+zDI+zcs,   du))),4))

## Figure S2
cmp <- rbindlist(list(su, ss))
qA <- ggplot(cmp, aes(med_r, frac_sorted, color=level))+geom_line(linewidth=1)+geom_point()+
  scale_x_log10()+scale_color_manual(values=c("SNP"="grey55","unit (eMLG)"="firebrick"))+
  labs(x="recombination (cM/Mb, log)", y="fraction sorted", color=NULL,
  title="(A) SNP-level vs independent-unit sorting vs recombination")+th_bw
qB <- ggplot(melt(su[,.(med_r, directional=frac_uni, bidirectional=frac_bi)], id.vars="med_r"),
  aes(med_r, value, color=variable))+geom_line(linewidth=1)+geom_point()+scale_x_log10()+
  labs(x="recombination (cM/Mb, log)", y="fraction of units", color=NULL,
  title="(B) unit-level: directional vs bidirectional")+th_bw
du[, rtert := cut(recomb, quantile(recomb,0:3/3,na.rm=TRUE), include.lowest=TRUE, labels=c("low r","mid r","high r"))]
du[, DIbin := cut(DI, quantile(DI,0:6/6,na.rm=TRUE), include.lowest=TRUE)]
pd <- du[!is.na(DIbin), .(med_DI=median(DI), prop_fixed=mean(prop_fixed,na.rm=TRUE)), by=.(rtert,DIbin)]
qC <- ggplot(pd, aes(med_DI, prop_fixed, color=rtert))+geom_line(linewidth=1)+geom_point()+
  labs(x="DI", y="prop_fixed", color=NULL, title="(C) DI effect within recombination tertiles")+th_bw
ggsave("Figures/tier2_sorting_vs_recomb.png", qA/(qB|qC), width=12, height=10, dpi=110)

## --------------------------------------------------------------------------
## 6. LD structure for neutral-null design: within-species vs admixture LD
##    (composite r^2 decay, matched panels of n = 15 diploids)
## --------------------------------------------------------------------------
DI_TH <- -25; MAXDIST <- 300000; MAXPART <- 60; N_IND <- 15; MAF_MIN <- 0.05
di <- map[DiagnosticIndex > DI_TH, which = TRUE]
rec_split <- median(map$recomb[di], na.rm = TRUE)
rows_of <- function(pops,n){ i<-which(sd$Population%in%pops); if(length(i)>n) i<-sample(i,n); i }
panels <- list(aqu=rows_of(aqu,N_IND), pol=rows_of(pol,N_IND),
               pooled=c(sample(which(sd$Population==aqu),8), sample(which(sd$Population==pol),7)),
               hybrid=rows_of("LangholmenW",N_IND))
dist_breaks <- c(0,1,2,5,10,20,50,100,200,300)*1000
ld_decay_panel <- function(rows){
  G <- GTs[rows, di, drop=FALSE]; info <- map[di, .(Chr,Pos,recomb)]
  poly <- apply(G,2,function(x){v<-x[!is.na(x)]; length(v)>=6 && length(unique(v))>1})
  out <- list()
  for (ch in unique(info$Chr)){
    j <- which(info$Chr==ch & poly); if(length(j)<2) next
    o <- order(info$Pos[j]); j <- j[o]; pos<-info$Pos[j]; rr<-info$recomb[j]; Gc<-G[,j,drop=FALSE]
    for (a in seq_len(length(j)-1)){
      b <- (a+1):min(a+MAXPART,length(j)); dd<-pos[b]-pos[a]; keep<-dd<=MAXDIST
      if(!any(keep)) next; b<-b[keep]; dd<-dd[keep]
      rv <- suppressWarnings(cor(Gc[,a], Gc[,b,drop=FALSE], use="pairwise.complete.obs"))
      out[[length(out)+1]] <- data.table(dist=dd, r2=as.numeric(rv)^2, recomb=rr[a])
    }
  }
  rbindlist(out)[is.finite(r2)]
}
decay <- rbindlist(lapply(names(panels), function(nm){ d<-ld_decay_panel(panels[[nm]]); d[,panel:=nm][] }))
decay[, rec_stratum := fifelse(recomb<=rec_split,"low recomb","high recomb")]
decay[, dbin := cut(dist, dist_breaks, labels=FALSE)]
decay[, dmid := (dist_breaks[dbin]+dist_breaks[dbin+1])/2]
lab <- c(aqu="aquilonia (within-sp)", pol="polyctena (within-sp)",
         pooled="pooled 50/50 (admixture)", hybrid="LangholmenW (hybrid)")
summ2 <- decay[, .(mean_r2=mean(r2), n=.N), by=.(panel, rec_stratum, dmid)]
summ2[, panel_lab := factor(lab[panel], levels=lab)]
pS3 <- ggplot(summ2[n>=20], aes(dmid/1000, mean_r2, color=panel_lab))+geom_line(linewidth=0.9)+geom_point(size=1.3)+
  facet_wrap(~rec_stratum)+scale_x_log10()+
  scale_color_manual(values=c("#1b9e77","#7570b3","#d95f02","grey40"), name=NULL)+
  labs(x="distance (kb, log)", y=expression("mean composite "*r^2),
       title=sprintf("LD decay at DI > %d markers (n = %d per panel)", DI_TH, N_IND))+
  theme_bw(base_size=13)+theme(legend.position="bottom", panel.grid.minor=element_blank())
ggsave("Figures/parent_LD_diagnostic.png", pS3, width=11, height=5, dpi=110)

cat("\nDone. Figures written to Figures/. See manuscript_notes/supplementary_methods.{html,pdf}\n")
