## Does the Module-D "structure" module (the low-recombination hub loci LABELLED
## structure by moduleD_network_build.R) separate individuals by LOCATION? If these
## hubs are founding-admixture co-ancestry retained in low-recombination regions (not
## pair-specific incompatibilities), a PCA over just those loci should recover the
## among-population geographic structure -- i.e. cluster individuals by population and
## line up with the genome-wide PCs. Reads data/moduleD_network.rds + the clustering +
## the hybrid sample_data. Writes Figures/moduleD_structure_pca.png. Run from repo root.
suppressPackageStartupMessages({ library(data.table) })
dir.create("Figures", showWarnings = FALSE)

N <- readRDS("data/moduleD_network.rds"); mn <- as.data.table(N$meta_nodes)
clust <- readRDS("data/eMLG_5loci_0025_cM05.rds"); eMLG <- clust$eMLG
e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1)
sd <- as.data.table(e1$sample_data); pop <- sd[match(rownames(eMLG), Sample_ID), Population]

## structure-module loci = member clusters of every structure-labelled meta-node
sm_ids <- intersect(unlist(strsplit(mn[structure == TRUE, members], ";")), colnames(eMLG))
X <- eMLG[, sm_ids, drop = FALSE]
X <- apply(X, 2, function(v) { v[!is.finite(v)] <- mean(v, na.rm = TRUE); v })
pca <- prcomp(X, center = TRUE, scale. = TRUE)
ve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sc <- as.data.table(pca$x[, 1:4]); sc[, pop := pop]

## how much of PC1/PC2 is between-population (does it separate by location?)
r2pop <- function(v) summary(lm(v ~ pop))$r.squared
## agreement with the genome-wide PCs (recapitulates the same structure axis?)
gw1 <- sd[match(rownames(eMLG), Sample_ID), PC1]; gw2 <- sd[match(rownames(eMLG), Sample_ID), PC2]
cat(sprintf("structure module: %d loci x %d individuals; PC1 %.1f%%, PC2 %.1f%% var\n", length(sm_ids), nrow(X), ve[1], ve[2]))
cat(sprintf("between-population R^2: PC1 %.2f, PC2 %.2f\n", r2pop(sc$PC1), r2pop(sc$PC2)))

## baselines: candidate-module loci and random size-matched differentiated loci
cl <- readRDS("data/moduleC_C3_cl.rds"); diff_ids <- intersect(cl[differentiated == TRUE, group_id], colnames(eMLG))
r2_set <- function(ids) { Xi <- eMLG[, ids, drop = FALSE]
  Xi <- apply(Xi, 2, function(v) { v[!is.finite(v)] <- mean(v, na.rm = TRUE); v })
  p <- prcomp(Xi, center = TRUE, scale. = TRUE)$x[, 1:2]; c(r2pop(p[, 1]), r2pop(p[, 2])) }
cand_ids <- intersect(unlist(strsplit(mn[structure == FALSE, members], ";")), colnames(eMLG))
set.seed(1); rnd <- rowMeans(sapply(1:20, function(i) r2_set(sample(diff_ids, length(sm_ids)))))
cat(sprintf("baseline between-pop R^2 (PC1/PC2): candidate module %.2f/%.2f ; random differentiated %.2f/%.2f\n",
    r2_set(cand_ids)[1], r2_set(cand_ids)[2], rnd[1], rnd[2]))
## PC-independent: mean per-locus among-population differentiation (FST-like R^2 of pop)
mean_r2 <- function(ids) mean(sapply(ids, function(g) { v <- eMLG[, g]; ok <- is.finite(v)
  if (sd(v[ok]) == 0) NA_real_ else summary(lm(v[ok] ~ pop[ok]))$r.squared }), na.rm = TRUE)
set.seed(1); rnd_pl <- mean(sapply(1:20, function(i) mean_r2(sample(diff_ids, length(sm_ids)))))
cat(sprintf("mean per-locus among-pop R^2: structure %.2f ; candidate %.2f ; random differentiated %.2f\n",
    mean_r2(sm_ids), mean_r2(cand_ids), rnd_pl))
cat(sprintf("|cor| with genome-wide PCs: PC1 vs gwPC1 %.2f / gwPC2 %.2f ; PC2 vs gwPC1 %.2f / gwPC2 %.2f\n",
    abs(cor(sc$PC1, gw1)), abs(cor(sc$PC1, gw2)), abs(cor(sc$PC2, gw1)), abs(cor(sc$PC2, gw2))))

## ---- plot: PC1-PC2 coloured by population, with population centroids labelled -----
pops <- sort(unique(pop)); pal <- setNames(grDevices::hcl.colors(length(pops), "Spectral"), pops)
cen <- sc[, .(PC1 = median(PC1), PC2 = median(PC2)), by = pop]
png("Figures/moduleD_structure_pca.png", width = 1500, height = 1500, res = 200)
par(mar = c(4, 4, 3, 1))
plot(sc$PC1, sc$PC2, col = pal[sc$pop], pch = 19, cex = 1.1,
     xlab = sprintf("PC1 (%.1f%%)", ve[1]), ylab = sprintf("PC2 (%.1f%%)", ve[2]),
     main = sprintf("Structure-module PCA (%d low-recomb hub loci) -- individuals coloured by population", length(sm_ids)))
text(cen$PC1, cen$PC2, cen$pop, cex = 0.6, font = 2)
legend("topright", bty = "n", cex = 0.7, legend = c(
  sprintf("between-pop R2: PC1=%.2f, PC2=%.2f", r2pop(sc$PC1), r2pop(sc$PC2)),
  sprintf("mean per-locus among-pop R2:"),
  sprintf("  structure %.2f  <  random-diff %.2f  <  candidate %.2f", mean_r2(sm_ids), rnd_pl, mean_r2(cand_ids))))
dev.off()
cat("wrote Figures/moduleD_structure_pca.png\n")
