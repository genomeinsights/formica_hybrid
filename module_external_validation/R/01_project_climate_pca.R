## module_external_validation / 01: project the external hybrid samples into the
## climate PCA used for the original hybrid data set (PC1/PC2 covariates).
##
## The original PC1/PC2 reproduce EXACTLY (|r| = 1) as prcomp() on the 19
## bio1-bio19 values of the 20 sampling sites (data/bioclimatic_variables.csv,
## unscaled bio* with scale.=TRUE == the *_scaled columns; PC2 sign flipped).
## New sites need bio1-19 from the SAME raster source; the resolution is
## identified here by reproducing the 20 original sites' values.
##
## Run from the repo root:  Rscript module_external_validation/R/01_project_climate_pca.R [path/to/wc_dir]
## Reads : data/bioclimatic_variables.csv, data/Sample_info_outlier_analysis_2026.txt,
##         data/NooraNewHybridSamples.xlsx, WorldClim 2.1 bio rasters (wc2.1_<res>_bio_<n>.tif)
## Writes: module_external_validation/data/new_sample_climate_pca.tsv
##         module_external_validation/results/climate_pca_projection_checks.txt

suppressPackageStartupMessages({ library(data.table); library(terra); library(readxl) })
args <- commandArgs(trailingOnly = TRUE)
WC <- if (length(args)) args[1] else "wc"
RES <- if (length(args) > 1) args[2] else "2.5m"
OUT <- "module_external_validation"

bc <- fread("data/bioclimatic_variables.csv")
bio <- paste0("bio", 1:19)
pcs <- unique(fread("data/Sample_info_outlier_analysis_2026.txt")[, .(Population, PC1, PC2)])
pcs[Population == "Nyrhispera74", Population := "Nyrhispera1"]
pcs[Population == "Nyrhispera75", Population := "Nyrhispera2"]
m <- merge(bc, pcs, by.x = "Location", by.y = "Population")
stopifnot(nrow(m) == 20)

## ---- (a) does the raster source reproduce the original site values? ---------
## WC = local directory of wc2.1_<RES>_bio_<n>.tif, or "remote" to read single pixels
## straight out of the WorldClim zip (GDAL /vsizip//vsicurl; avoids the 10 GB 30s download)
lyr <- if (WC == "remote") sprintf("/vsizip//vsicurl/https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_%s_bio.zip/wc2.1_%s_bio_%d.tif", RES, RES, 1:19) else file.path(WC, sprintf("wc2.1_%s_bio_%d.tif", RES, 1:19))
r <- rast(lyr); names(r) <- bio
ext_orig <- as.data.frame(extract(r, as.matrix(m[, .(Longitude, Latitude)])))[, bio]
chk <- sapply(bio, function(b) cor(ext_orig[[b]], m[[b]]))
dif <- sapply(bio, function(b) max(abs(ext_orig[[b]] - m[[b]])))
cat("Reproduction of original site bio values from", RES, "WorldClim 2.1:\n")
print(round(rbind(r = chk, max_abs_diff = dif), 3))

## ---- (b) PCA on the original 20 sites; check against stored PC1/PC2 ---------
p <- prcomp(m[, ..bio], scale. = TRUE)
cat("\ncor(PC1)=", round(cor(p$x[, 1], m$PC1), 5), " cor(PC2)=", round(cor(p$x[, 2], m$PC2), 5), "\n")
sgn <- sign(cor(p$x[, 2], m$PC2))     # stored PC2 has the opposite sign

## ---- (c) project the new samples ---------------------------------------------
nw <- as.data.table(read_excel("data/NooraNewHybridSamples.xlsx"))
ex <- as.data.frame(extract(r, as.matrix(nw[, .(longitude, latitude)])))[, bio]
if (anyNA(ex)) message("NA bio values for: ", paste(nw$ID_short[!complete.cases(ex)], collapse = ", "))
Z <- predict(p, ex)
nw[, `:=`(PC1 = Z[, 1], PC2 = sgn * Z[, 2])]
nw[, c(bio) := as.data.table(ex)]
## extrapolation check: distance in PC space beyond the original range
rng <- sapply(1:2, function(i) range(m[[paste0("PC", i)]]))
nw[, outside_orig_range := PC1 < rng[1, 1] | PC1 > rng[2, 1] | PC2 < rng[1, 2] | PC2 > rng[2, 2]]
fwrite(nw[, c("ID_short", "BAM_ID", "country", "latitude", "longitude", "Elevation", "PC1", "PC2", "outside_orig_range", bio), with = FALSE],
       file.path(OUT, "data/new_sample_climate_pca.tsv"), sep = "\t")
cat("\nOriginal PC range: PC1", round(rng[, 1], 2), " PC2", round(rng[, 2], 2), "\n")
print(nw[, .(ID_short, country, latitude, PC1 = round(PC1, 2), PC2 = round(PC2, 2), outside_orig_range)])
cat("\nvariance explained (orig PCA):", round(summary(p)$importance[2, 1:3], 3), "\n")
