## =========================================================
## Tier 1: genomic architecture of differentiation (parents only)
## =========================================================
## Q: do high-DI / large-cluster loci sit in low-recombination regions,
##    and is the DI-r association a within-species DIVERSITY artifact
##    (relative differentiation inflated where pi is low) or GENUINE excess
##    between-species divergence (elevated dxy)?  [Cruickshank & Hahn 2014]

suppressMessages({library(data.table); library(ggplot2); library(patchwork)})

e <- new.env(); load("./data/hybrids_and_parents_maf005.Rdata", envir=e)
GTs <- e$GTs_with_parents; sd <- e$sample_data_with_parents; map <- copy(e$map_hyb_005)

## ---- per-marker recombination rate (cM/Mb) ----
rec <- fread("./data/Frufa_DTOL_PR.ref_genome.recmap")
setnames(rec, c("chr","pos","cM","cMMb")); rec[, Chr := sub("chromosome_","Chr",chr)]
map[, recomb := NA_real_]
for (ch in unique(map$Chr)) {
  r <- rec[Chr==ch]; if (nrow(r) < 2) next
  idx <- map[, which(Chr==ch)]
  map[idx, recomb := approx(r$pos, r$cMMb, xout=map$Pos[idx], rule=2)$y]
}

## ---- per-marker cluster size (eMLG complexity reduction) ----
g <- readRDS("data/eMLG_5loci_0025_cM05.rds")$groups
memb <- data.table(marker = unlist(g$members),
                   cluster_size = rep(g$n_loci, lengths(g$members)))
map <- memb[map, on="marker"]

## ---- parental allele frequencies -> pi (within), dxy (absolute), Fst (relative) ----
aqu_rows <- which(sd$Population=="aquilonia_parent")
pol_rows <- which(sd$Population=="polyctena_parent")
pa <- colMeans(GTs[aqu_rows, ], na.rm=TRUE)/2      # aligned to map (column order)
pp <- colMeans(GTs[pol_rows, ], na.rm=TRUE)/2
Hs_a <- 2*pa*(1-pa); Hs_p <- 2*pp*(1-pp)
map[, pi_within := (Hs_a + Hs_p)/2]                 # within-species diversity
map[, dxy       := pa*(1-pp) + pp*(1-pa)]           # absolute divergence
pbar <- (pa+pp)/2; Ht <- 2*pbar*(1-pbar)
map[, Fst := ifelse(Ht > 0, (Ht - (Hs_a+Hs_p)/2)/Ht, NA_real_)]  # relative differentiation

d <- map[is.finite(recomb) & is.finite(DiagnosticIndex) & !is.na(cluster_size)]
cat("markers used:", nrow(d), "\n")

## ---- correlations (Spearman) ----
sp <- function(a,b) round(cor(a,b,method="spearman",use="complete.obs"),3)
cat("\n=== Spearman correlations ===\n")
cat("  DI vs recomb        :", sp(d$DiagnosticIndex, d$recomb), "\n")
cat("  DI vs cluster_size  :", sp(d$DiagnosticIndex, d$cluster_size), "\n")
cat("  DI vs pi_within     :", sp(d$DiagnosticIndex, d$pi_within), "\n")
cat("  DI vs dxy           :", sp(d$DiagnosticIndex, d$dxy), "\n")
cat("  recomb vs cluster   :", sp(d$recomb, d$cluster_size), "\n")
cat("  recomb vs pi_within :", sp(d$recomb, d$pi_within), "\n")
cat("  recomb vs dxy       :", sp(d$recomb, d$dxy), "\n")
cat("  recomb vs Fst       :", sp(d$recomb, d$Fst), "\n")

## ---- binned by recombination quantile ----
d[, rbin := cut(recomb, quantile(recomb, 0:10/10, na.rm=TRUE), include.lowest=TRUE, labels=FALSE)]
tab <- d[, .(med_recomb=round(median(recomb),1), n=.N,
             DI=round(mean(DiagnosticIndex),1),
             Fst=round(mean(Fst,na.rm=TRUE),3),
             dxy=round(mean(dxy,na.rm=TRUE),3),
             pi_within=round(mean(pi_within,na.rm=TRUE),3),
             cluster_size=round(mean(cluster_size),1)), by=rbin][order(rbin)]
cat("\n=== means by recombination decile (rbin 1 = lowest recomb) ===\n"); print(tab)

## ---- plots ----
th <- theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
## CH panel: standardise DI, dxy, pi vs recomb so trends are comparable
z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
ch <- melt(tab[, .(med_recomb, DI=z(DI), dxy=z(dxy), pi_within=z(pi_within), Fst=z(Fst))],
           id.vars="med_recomb", variable.name="metric", value.name="z")
p_ch <- ggplot(ch, aes(med_recomb, z, color=metric)) +
  geom_line(linewidth=1) + geom_point() +
  scale_x_log10() +
  labs(x="recombination (cM/Mb, log)", y="standardised (z)",
       title="(A) Cruickshank-Hahn: relative (DI/Fst) vs absolute (dxy) vs pi vs recomb") + th

p_clust <- ggplot(tab, aes(med_recomb, cluster_size)) + geom_line(linewidth=1,color="grey30") + geom_point() +
  scale_x_log10() + scale_y_log10() +
  labs(x="recombination (cM/Mb, log)", y="mean cluster size", title="(B) cluster size vs recomb") + th

## DI vs cluster size (binned by cluster-size quantile, hexbin too dense)
d[, cbin := cut(cluster_size, unique(quantile(cluster_size, 0:12/12)), include.lowest=TRUE)]
dc <- d[, .(med_cs=median(cluster_size), DI=mean(DiagnosticIndex), dxy=mean(dxy,na.rm=TRUE),
            pi_within=mean(pi_within,na.rm=TRUE)), by=cbin][!is.na(cbin)][order(med_cs)]
p_dic <- ggplot(dc, aes(med_cs, DI)) + geom_line(linewidth=1,color="firebrick") + geom_point() +
  scale_x_log10() + labs(x="cluster size (log)", y="mean DI", title="(C) DI vs cluster size") + th

p <- (p_ch | p_clust) / (p_dic | plot_spacer())
dir.create("Figures", showWarnings=FALSE)
ggsave("Figures/tier1_architecture.png", p, width=13, height=9, dpi=110)
cat("\nsaved Figures/tier1_architecture.png\n")
