## module_external_validation / 01b: WorldClim-consistent climate PCA.
##
## 01_ showed the original bioclimatic_variables.csv is NOT plain WorldClim 2.1 (1970-2000)
## at 30 s: the original values are ~1.65 C warmer in bio1, ~2.9 C warmer in bio6 and 7-9%
## wetter (uneven across sites; Sielva differs by only 0.2 C), i.e. a different product/period.
## Mixing the two sources would bias the projection, so here the PCA is rebuilt from WorldClim
## ALONE: bio1-19 at the original 20 sites -> prcomp(scale.=TRUE) -> project the new samples.
## Fidelity check: correlate this WorldClim-based PC1/PC2 with the stored original PC1/PC2.
## LangholmenR falls on a no-data (sea) cell; it is ~at the same position as LangholmenW and
## has identical stored PCs, so LangholmenW's raster values are used for it.
##
## Run from the repo root:  Rscript module_external_validation/R/01b_climate_pca_worldclim_consistent.R [wc_dir]
## Writes: data/new_sample_climate_pca_wc.tsv, results/climate_pca_wc_consistent.txt, results/climate_pca_wc_consistent.png

suppressPackageStartupMessages({ library(data.table); library(terra); library(readxl); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE); WC <- path.expand(if (length(args)) args[1] else "~/data/worldclim_2.1_30s")
OUT <- "module_external_validation"; bio <- paste0("bio", 1:19)
r <- rast(file.path(WC, sprintf("wc2.1_30s_bio_%d.tif", 1:19))); names(r) <- bio
ex <- function(lon, lat) as.data.frame(extract(r, cbind(lon, lat)))[, bio]

bc <- fread("data/bioclimatic_variables.csv")
pcs <- unique(fread("data/Sample_info_outlier_analysis_2026.txt")[, .(Population, PC1, PC2)])
pcs[Population == "Nyrhispera74", Population := "Nyrhispera1"]; pcs[Population == "Nyrhispera2" , Population := "Nyrhispera2"]
pcs[Population == "Nyrhispera75", Population := "Nyrhispera2"]
m <- merge(bc[, .(Location, Longitude, Latitude)], pcs, by.x = "Location", by.y = "Population")
W <- ex(m$Longitude, m$Latitude)
na_site <- which(!complete.cases(W)); if (length(na_site)) { stopifnot(m$Location[na_site] == "LangholmenR"); W[na_site, ] <- W[m$Location == "LangholmenW", ] }
p <- prcomp(W, scale. = TRUE)
sg <- c(sign(cor(p$x[, 1], m$PC1)), sign(cor(p$x[, 2], m$PC2)))      # align arbitrary signs to the stored PCs
m[, `:=`(wcPC1 = sg[1] * p$x[, 1], wcPC2 = sg[2] * p$x[, 2])]

nw <- as.data.table(read_excel("data/NooraNewHybridSamples.xlsx"))
Z <- predict(p, ex(nw$longitude, nw$latitude)); stopifnot(!anyNA(Z))
nw[, `:=`(PC1 = sg[1] * Z[, 1], PC2 = sg[2] * Z[, 2])]
rng <- sapply(c("wcPC1", "wcPC2"), function(v) range(m[[v]]))
nw[, outside_orig_range := PC1 < rng[1, 1] | PC1 > rng[2, 1] | PC2 < rng[1, 2] | PC2 > rng[2, 2]]
## which variables drive the extreme PC2 samples?
ld <- p$rotation[, 2] * sg[2]
fwrite(nw[, .(ID_short, BAM_ID, country, latitude, longitude, Elevation, PC1, PC2, outside_orig_range)], file.path(OUT, "data/new_sample_climate_pca_wc.tsv"), sep = "\t")

sink(file.path(OUT, "results/climate_pca_wc_consistent.txt"), split = TRUE)
cat("WorldClim-only PCA on 20 original sites; variance explained:", round(summary(p)$importance[2, 1:3], 3), "\n")
cat("Fidelity vs stored original PCs: PC1 r =", round(cor(m$wcPC1, m$PC1), 4), " PC2 r =", round(cor(m$wcPC2, m$PC2), 4),
    " | Spearman", round(cor(m$wcPC1, m$PC1, method = "spearman"), 3), round(cor(m$wcPC2, m$PC2, method = "spearman"), 3), "\n")
cat("WorldClim-PC range (orig sites): PC1", round(rng[, 1], 2), " PC2", round(rng[, 2], 2), "\n")
cat("\nTop PC2 loadings:\n"); print(round(sort(ld, decreasing = TRUE)[c(1:4, 16:19)], 3))
cat("\nNew samples projected:\n"); print(nw[, .(ID_short, country, lat = round(latitude, 2), PC1 = round(PC1, 2), PC2 = round(PC2, 2), outside_orig_range)])
cat("\nn outside original range:", sum(nw$outside_orig_range), "of", nrow(nw), "\n")
sink()
g <- ggplot() + geom_point(data = m, aes(wcPC1, wcPC2), colour = "grey40", size = 2) +
  geom_point(data = nw, aes(PC1, PC2, colour = country), size = 2.5) + geom_text(data = nw, aes(PC1, PC2, label = ID_short), size = 2.5, nudge_y = .4) +
  labs(title = "WorldClim-based climate PCA: original sites (grey) and external samples", x = "PC1", y = "PC2") + theme_classic()
ggsave(file.path(OUT, "results/climate_pca_wc_consistent.png"), g, width = 8, height = 6, dpi = 200)
