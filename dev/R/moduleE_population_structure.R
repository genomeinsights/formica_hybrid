## =============================================================================
## Module E -- are the 20 hybrid populations independent? Neutral vs diagnostic
## =============================================================================
## The neutral null assumes 20 INDEPENDENT demes and over-differentiates them
## (neutral F_ST 0.14 vs empirical 0.05). Two explanations for the low empirical
## value: large Ne (little drift) OR non-independence (shared ancestry / gene flow
## among "populations"). These are separable:
##   - NEUTRAL markers (DI<=-75) carry demographic history -> test independence.
##   - DIAGNOSTIC markers (DI>-25) carry the ancestry/sorting signal (Nouhaud).
## If populations cluster the same way on both -> the ancestry correlation is
## shared demography (non-independence). If neutral shows no clustering but
## diagnostic does -> demographically independent + parallel sorting.
##
## Outputs: pairwise neutral F_ST matrix + dendrogram, PCA (neutral & diagnostic),
##          per-pop ancestry, and a Mantel-style test of neutral-distance vs
##          ancestry-distance.
##   Figures/moduleE_structure.{pdf,png} ; data/moduleE_sim/moduleE_structure.rds
## =============================================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

FORMICA <- "/Users/petrikem/gitlab/formica_hybrid"
POOL    <- file.path(FORMICA, "data/moduleE_founders/moduleE_founder_haplotypes.rds")
EMP_RD  <- file.path(FORMICA, "data/hybrids_only_maf005.Rdata")
OUTRDS  <- file.path(FORMICA, "data/moduleE_sim/moduleE_structure.rds")
FIGDIR  <- file.path(FORMICA, "Figures"); dir.create(FIGDIR, showWarnings=FALSE)

## ---- markers by DI (from founder pool = parental reference) -----------------
ph <- readRDS(POOL); pmap <- as.data.table(ph$map)
p_pool <- (colSums(ph$aqu)+colSums(ph$pol))/(nrow(ph$aqu)+nrow(ph$pol))
pmap[, parent_maf := pmin(p_pool,1-p_pool)]
f_aq_par <- colMeans(ph$aqu); f_pol_par <- colMeans(ph$pol)
neutral_mk <- pmap[parent_maf>=0.15 & DI<=-90, marker]   # strictly neutral (F_ST~0.05)
diag_mk    <- pmap[parent_maf>=0.15 & DI> -25, marker]

## ---- empirical hybrids ------------------------------------------------------
e <- new.env(); load(EMP_RD, envir=e)
GTh <- e$GTs_hybrids_005; sdh <- as.data.table(e$sample_data); maph <- as.data.table(e$map_hyb_005)
sdh <- sdh[!Population %in% c("aquilonia_parent","polyctena_parent")]
hyb <- match(sdh$Sample_ID, rownames(GTh)); if (anyNA(hyb)) hyb <- which(!e$sample_data$Population %in% c("aquilonia_parent","polyctena_parent"))
pops <- unique(sdh$Population)
sub <- function(mk){ ci <- match(mk, maph$marker); ci <- ci[!is.na(ci)]
  g <- GTh[hyb, ci, drop=FALSE]; g }

## per-pop allele freq + allele count n (2*diploids)
pop_freq <- function(g){
  do.call(rbind, lapply(pops, function(p){ x <- g[sdh$Population==p,,drop=FALSE]
    colMeans(x, na.rm=TRUE)/2 })) }
Gn <- sub(neutral_mk); Fn <- pop_freq(Gn); rownames(Fn) <- pops
npop <- sapply(pops, function(p) 2*sum(sdh$Population==p))

## ---- Hudson pairwise F_ST (sample-size corrected), neutral markers ---------
hudson <- function(p1,p2,n1,n2){
  num <- (p1-p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
  den <- p1*(1-p2) + p2*(1-p1)
  ok <- is.finite(num) & is.finite(den) & den>0
  sum(num[ok])/sum(den[ok]) }
K <- length(pops); FST <- matrix(0, K, K, dimnames=list(pops,pops))
for (i in 1:(K-1)) for (j in (i+1):K){
  f <- hudson(Fn[i,], Fn[j,], npop[i], npop[j]); FST[i,j] <- FST[j,i] <- max(f,0) }
cat(sprintf("neutral markers: %d ; pairwise F_ST  min=%.3f median=%.3f max=%.3f\n",
            length(neutral_mk), min(FST[upper.tri(FST)]), median(FST[upper.tri(FST)]), max(FST[upper.tri(FST)])))

hc <- hclust(as.dist(FST), method="average")

## ---- PCA (individual) on neutral and diagnostic markers --------------------
pca_of <- function(g){
  g <- g[, apply(g,2,function(x){ p<-mean(x,na.rm=TRUE)/2; is.finite(p) && pmin(p,1-p)>=0.05 }), drop=FALSE]
  for (k in which(colSums(is.na(g))>0)) g[is.na(g[,k]),k] <- mean(g[,k],na.rm=TRUE)
  pr <- prcomp(g, center=TRUE, scale.=FALSE)
  data.table(Population=sdh$Population, PC1=pr$x[,1], PC2=pr$x[,2],
             ve=paste0(round(100*(pr$sdev^2/sum(pr$sdev^2))[1:2],1),"%")[1:2][c(1,2)]) }
pn <- pca_of(Gn); pd <- pca_of(sub(diag_mk))

## per-pop mean ancestry (oriented aq) on diagnostic markers
Gd <- sub(diag_mk); dci <- match(diag_mk, pmap$marker)
sgn <- sign(f_aq_par[diag_mk] - f_pol_par[diag_mk])
anc <- sapply(pops, function(p){ f <- colMeans(Gd[sdh$Population==p,,drop=FALSE], na.rm=TRUE)/2
  mean(ifelse(sgn>0, f, 1-f)[sgn!=0], na.rm=TRUE) })
names(anc) <- pops

saveRDS(list(FST=FST, hc=hc, pca_neutral=pn, pca_diag=pd, ancestry=anc,
             n_neutral=length(neutral_mk), n_diag=length(diag_mk)), OUTRDS)

## ---- Mantel-style: does neutral demographic distance predict ancestry diff? -
anc_d <- as.matrix(dist(anc)); ut <- upper.tri(FST)
mant <- cor(FST[ut], anc_d[ut], method="spearman")
cat(sprintf("Spearman(neutral F_ST, |ancestry diff|) across pop pairs: %.3f\n", mant))
## least-differentiated (candidate non-independent) pairs
pr <- data.table(a=rownames(FST)[row(FST)[ut]], b=colnames(FST)[col(FST)[ut]], fst=FST[ut])[order(fst)]
cat("\n5 LEAST-differentiated pop pairs (neutral F_ST -- candidate non-independence):\n"); print(head(pr,5))
cat("5 MOST-differentiated:\n"); print(tail(pr,5))

## ---- figures ---------------------------------------------------------------
## F_ST heatmap ordered by hierarchical clustering (low-F_ST blocks = candidate
## non-independent groups)
ord <- hc$order; ll <- rownames(FST)[ord]
hm <- as.data.table(as.table(FST)); setnames(hm, c("a","b","fst"))
hm[, `:=`(a=factor(a, levels=ll), b=factor(b, levels=ll))]
p_tree <- ggplot(hm, aes(a, b, fill=fst)) + geom_tile() +
  scale_fill_viridis_c(option="magma", name="F_ST", direction=-1) +
  labs(title="Neutral-marker pairwise F_ST (clustered)", x=NULL, y=NULL) +
  theme_minimal(base_size=7) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=5),
        axis.text.y=element_text(size=5), panel.grid=element_blank(), legend.key.size=unit(6,"pt"))
p_pn <- ggplot(pn, aes(PC1,PC2,colour=Population)) + geom_point(size=1) +
  stat_ellipse(level=0.6, linewidth=0.3) + guides(colour="none") +
  labs(title=sprintf("PCA: NEUTRAL markers (%d)", length(neutral_mk))) + theme_bw(base_size=8)
p_pd <- ggplot(pd, aes(PC1,PC2,colour=Population)) + geom_point(size=1) +
  stat_ellipse(level=0.6, linewidth=0.3) + guides(colour="none") +
  labs(title=sprintf("PCA: DIAGNOSTIC markers (%d)", length(diag_mk))) + theme_bw(base_size=8)
fig <- (p_tree | p_pn | p_pd) + plot_annotation(
  title="Population structure among the 20 hybrid populations: demographic (neutral) vs ancestry (diagnostic)")
ggsave(file.path(FIGDIR,"moduleE_structure.pdf"), fig, width=13, height=4.5)
ggsave(file.path(FIGDIR,"moduleE_structure.png"), fig, width=13, height=4.5, dpi=130)
cat("\nsaved: figures + ", OUTRDS, "\n")
