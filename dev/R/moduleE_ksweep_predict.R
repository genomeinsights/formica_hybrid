## =============================================================================
## Module E -- K (Ne) sweep at fixed founding (450/195): does LD scale with 1/K,
## and can we PREDICT the K that reaches the empirical level?
## =============================================================================
## Drift-LD theory: at effective size Ne, E[r^2] ~ 1/(1+4*Ne*c) (+ sampling floor).
## The queen-fecundity skew makes Ne = alpha*K, so at fixed generation LD level and
## neutral F_ST should be ~ LINEAR in 1/K. We run K DOWNWARD (cheap) to trace the
## line, then EXTRAPOLATE upward to the K where each statistic hits its empirical
## target. A good linear-in-1/K fit both confirms the Ne mechanism and predicts K*.
##
## Per K (founding 450/195, gen 60, pooled, per-deme N matched to empirical):
##   - mean local ld_w  (linked LD landscape level; the 2.2x-too-high statistic)
##   - background r^2   (within-deme INTER-chromosome unlinked pairs; pure Ne probe)
##   - neutral-bin F_ST and pi
## Output: data/moduleE_sim/moduleE_ksweep_predict.rds + Figures/moduleE_ksweep_predict.{pdf,png}
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

FORMICA  <- "/Users/petrikem/gitlab/formica_hybrid"
SIMDIR   <- file.path(FORMICA, "moduleE_slim/output_foundersweep")
FVCF_DIR <- file.path(FORMICA, "moduleE_slim/founders/maf015_DIstrat1500_x40")
BUNDLE   <- file.path(FORMICA, "moduleE_slim/inputs/empirical_bundle.rds")
EMP_RD   <- file.path(FORMICA, "data/hybrids_only_maf005.Rdata")
OUTRDS   <- file.path(FORMICA, "data/moduleE_sim/moduleE_ksweep_predict.rds")
FIGDIR   <- file.path(FORMICA, "Figures"); dir.create(FIGDIR, showWarnings=FALSE)
source(file.path(FORMICA, "moduleE_slim/scripts/parallelism_stats.R"))

TAG <- "Naq450_Npol195"; GEN <- 60L
WIN_BP <- 100000L; MAXPART <- 100L; MAF_MIN <- 0.1
DI_BREAKS <- c(-Inf,-90,-75,-60,-50,-40,-30,-25,-20,-15,Inf)
FIX_TH<-0.15; SORT_TH<-0.5; MIN_PARENT_MAF<-0.15; NULL_PROB<-0.5

## ---- markers, empirical reference ------------------------------------------
fm <- rbindlist(lapply(list.files(FVCF_DIR,"founders_ch.*[.]vcf$",full.names=TRUE),
  function(f) fread(f,skip="#CHROM",select=c(1,2,3),header=TRUE,col.names=c("chs","Pos","marker"))))
fm[, `:=`(Chr=sub("ch","Chr",chs), chr_id=as.integer(sub("ch","",chs)))]; setorder(fm,chr_id,Pos)
markers <- fm$marker
b <- readRDS(BUNDLE); uni <- as.data.table(b$universe); setkey(uni,marker)
um <- uni[.(markers)]; DI <- um$DI; dibin <- cut(DI, DI_BREAKS)
f_aq <- b$f_aq_par[markers]; f_pol <- b$f_pol_par[markers]
parent_maf <- setNames(um$parent_maf, markers); emp_ns <- b$emp_ns
neu <- levels(dibin)[1]; hi <- levels(dibin)[nlevels(dibin)]

local_ld <- function(geno, Chr, Pos) {
  res <- rep(NA_real_, ncol(geno))
  for (ch in unique(Chr)) { j <- which(Chr==ch); j <- j[order(Pos[j])]
    pos <- Pos[j]; G <- geno[,j,drop=FALSE]; p <- colMeans(G,na.rm=TRUE)/2
    poly <- is.finite(p) & pmin(p,1-p)>=MAF_MIN
    for (a in seq_along(j)) { if(!poly[a]) next
      bb <- which(pos<=pos[a]+WIN_BP & seq_along(pos)>a & poly); if(length(bb)>MAXPART) bb<-bb[seq_len(MAXPART)]
      if(!length(bb)) next
      r <- suppressWarnings(cor(G[,a], G[,bb,drop=FALSE], use="pairwise.complete.obs"))
      res[j[a]] <- median(as.numeric(r)^2, na.rm=TRUE) } }
  pmax(res,0)
}
bg_r2 <- function(geno, Chr) {                       # inter-chromosome unlinked mean r^2 (400 mk)
  set.seed(7); s <- sample(ncol(geno), min(400,ncol(geno))); G <- geno[,s,drop=FALSE]; ch <- Chr[s]
  v <- apply(G,2,var); G <- G[,v>0,drop=FALSE]; ch <- ch[v>0]
  R <- suppressWarnings(cor(G))^2; ut <- upper.tri(R); un <- outer(ch,ch,"!=")[ut]
  mean(R[ut][un], na.rm=TRUE)
}
read_dose <- function(f) { dt <- fread(f,skip="#CHROM",header=TRUE,sep="\t",showProgress=FALSE)
  mk <- paste(sub("ch","Chr",dt[[1]]),dt[[2]],sep=":"); G <- as.matrix(dt[,10:ncol(dt)]); gt <- sub(":.*$","",G)
  d <- matrix(suppressWarnings(as.integer(sub("[|/].*$","",gt)))+suppressWarnings(as.integer(sub("^.*[|/]","",gt))),nrow=nrow(G))
  m <- match(markers,mk); keep <- !is.na(m); D <- matrix(0L,length(markers),ncol(d)); D[keep,] <- d[m[keep],,drop=FALSE]; t(D) }
by_bin <- function(mm) { pm <- rbind(mm, aquilonia_parent=f_aq*2, polyctena_parent=f_pol*2); colnames(pm) <- markers
  s <- parallelism_stats(list(pop_means=pm), hybrid_pops=rownames(mm), aqu_pops="aquilonia_parent",
        pol_pops="polyctena_parent", fix_th=FIX_TH, sort_th=SORT_TH, null_prob=NULL_PROB,
        parent_maf=parent_maf, min_parent_maf=MIN_PARENT_MAF)
  sc <- s[match(markers,marker),sort_class]; p <- mm/2; pb <- colMeans(p,na.rm=TRUE)
  Hs <- colMeans(2*p*(1-p),na.rm=TRUE); Ht <- 2*pb*(1-pb)
  data.table(dibin=dibin, Hs=Hs, Ht=Ht)[, .(pi=mean(Hs,na.rm=TRUE), fst=sum(Ht-Hs,na.rm=TRUE)/sum(Ht,na.rm=TRUE)), by=dibin] }

## ---- empirical targets ------------------------------------------------------
e <- new.env(); load(EMP_RD, envir=e); GTh <- e$GTs_hybrids_005; sdh <- as.data.table(e$sample_data)
maph <- as.data.table(e$map_hyb_005); hyb <- which(!sdh$Population %in% c("aquilonia_parent","polyctena_parent"))
emp_gt <- GTh[hyb, match(markers,maph$marker), drop=FALSE]; emp_gt[is.na(emp_gt)] <- 0L; rm(e); gc()
emp_ldw <- mean(local_ld(emp_gt, fm$Chr, fm$Pos), na.rm=TRUE)
emp_bg  <- bg_r2(emp_gt, fm$Chr)
emp_M <- (b$emp_mean[, markers, drop=FALSE])/1000; eb <- by_bin(emp_M)
emp_fst <- eb[dibin==neu,fst]; emp_pi <- eb[dibin==neu,pi]
cat(sprintf("EMPIRICAL: mean ld_w=%.3f  bg r^2=%.3f  neutral F_ST=%.3f  pi=%.3f\n", emp_ldw, emp_bg, emp_fst, emp_pi))

## ---- per K ------------------------------------------------------------------
Ks <- sort(as.integer(unique(na.omit(sub(sprintf(".*_%s_K([0-9]+)_rep.*",TAG),"\\1",
        grep(sprintf("_%s_K",TAG), list.files(SIMDIR,"\\.vcf$"), value=TRUE))))))
cat("K values:", paste(Ks,collapse=", "), "\n")
res <- rbindlist(lapply(Ks, function(K) {
  parts <- lapply(1:20, function(i){ f <- file.path(SIMDIR, sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf",TAG,K,i,GEN))
    if(!file.exists(f)) return(NULL); g <- read_dose(f); tgt <- emp_ns[(i-1)%%length(emp_ns)+1]
    set.seed(100+i); if(nrow(g)>tgt) g <- g[sample(nrow(g),tgt),,drop=FALSE]; g })
  gt <- do.call(rbind, parts)
  ldw <- mean(local_ld(gt, fm$Chr, fm$Pos), na.rm=TRUE); bg <- bg_r2(gt, fm$Chr)
  ## per-pop means for F_ST/pi (rows = demes, sample matched)
  pm <- t(vapply(1:20, function(i){ f <- file.path(SIMDIR, sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf",TAG,K,i,GEN))
    if(!file.exists(f)) return(rep(NA_real_,length(markers))); D <- read_dose(f); tgt <- emp_ns[(i-1)%%length(emp_ns)+1]
    set.seed(100+i); if(nrow(D)>tgt) D <- D[sample(nrow(D),tgt),,drop=FALSE]; colMeans(D,na.rm=TRUE) }, numeric(length(markers))))
  rownames(pm) <- paste0("deme",1:20); colnames(pm) <- markers; bb <- by_bin(pm)
  data.table(K=K, ldw=ldw, bg=bg, neu_fst=bb[dibin==neu,fst], neu_pi=bb[dibin==neu,pi])
}))
setorder(res, K); print(res)

## ---- fit vs 1/K and predict K* ---------------------------------------------
res[, invK := 1/K]
predict_K <- function(y, target) { m <- lm(y ~ invK, data=res); a <- coef(m)[1]; slope <- coef(m)[2]
  Kstar <- slope/(target - a); list(a=unname(a), slope=unname(slope), r2=summary(m)$r.squared, Kstar=unname(Kstar)) }
pl <- predict_K(res$ldw, emp_ldw); pf <- predict_K(res$neu_fst, emp_fst); pb2 <- predict_K(res$bg, emp_bg)
saveRDS(list(res=res, emp=list(ldw=emp_ldw,bg=emp_bg,fst=emp_fst,pi=emp_pi),
             fit=list(ldw=pl,fst=pf,bg=pb2)), OUTRDS)

fmtK <- function(p, target, floor_name){ if (p$a >= target)
    sprintf("floor a=%.3f is ABOVE target %.3f -> no finite K reaches it (%s ceiling)", p$a, target, floor_name)
  else if (p$Kstar < 0) sprintf("K*=%.0f (check)", p$Kstar)
  else sprintf("predicted K* = %.0f  (fit: y=%.3f + %.2f/K, R^2=%.3f)", p$Kstar, p$a, p$slope, p$r2) }
cat("\n=== linear-in-1/K fits + extrapolation to empirical ===\n")
cat(sprintf("  mean ld_w  -> target %.3f : %s\n", emp_ldw, fmtK(pl, emp_ldw, "recombination/panel")))
cat(sprintf("  neutral F_ST-> target %.3f : %s\n", emp_fst, fmtK(pf, emp_fst, "sampling")))
cat(sprintf("  background r^2 -> target %.3f : %s\n", emp_bg, fmtK(pb2, emp_bg, "sampling")))

## ---- figure: y vs 1/K with fit line, empirical target, predicted K* ---------
mk_panel <- function(col, yl, emp, p, ttl){
  gx <- data.table(invK=seq(0, max(res$invK)*1.05, length=100)); gx[, yhat := p$a + p$slope*invK]
  ggplot() + geom_line(data=gx, aes(invK, yhat), linetype=2, colour="grey40") +
    geom_point(data=res, aes(invK, get(col)), colour="#1b9e77", size=2) +
    geom_hline(yintercept=emp, colour="#d95f02") +
    annotate("text", x=0, y=emp, label="empirical", hjust=-0.05, vjust=-0.5, colour="#d95f02", size=3) +
    { if (is.finite(p$Kstar) && p$Kstar>0 && p$a<emp) geom_point(aes(x=1/p$Kstar,y=emp), colour="black", shape=4, size=3) } +
    labs(x="1/K", y=yl, title=ttl) + theme_bw(base_size=9)
}
p1 <- mk_panel("ldw", "mean local ld_w", emp_ldw, pl,
        if (pl$a<emp_ldw) sprintf("local LD: predict K*~%.0f", pl$Kstar) else "local LD: floor above empirical")
p2 <- mk_panel("neu_fst","neutral F_ST", emp_fst, pf,
        if (pf$a<emp_fst) sprintf("neutral F_ST: predict K*~%.0f", pf$Kstar) else "neutral F_ST: floor above empirical")
ggsave(file.path(FIGDIR,"moduleE_ksweep_predict.pdf"), p1|p2, width=10, height=4)
ggsave(file.path(FIGDIR,"moduleE_ksweep_predict.png"), p1|p2, width=10, height=4, dpi=130)
cat("\nsaved: ", OUTRDS, " ; Figures/moduleE_ksweep_predict.{pdf,png}\n")
