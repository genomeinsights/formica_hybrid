## =========================================================
## module_di25 (high-DI analyses) -- DIEM circos plot, per-SNP (R recreation)
## =========================================================
## Recreates, in R, the DIEM genotype circos originally produced in Python,
## at the per-SNP level, so extra tracks/annotation can be added later.
##
## Input : data/species_diagnostic_markers_DI25_20pops.tsv.gz
##   195 individual genotype columns + `chromosome` + `position`.
##   Genotypes coded relative to the two parental species:
##     0 = missing, 1 = hom F. aquilonia, 2 = het, 3 = hom F. polyctena.
##   (Confirmed: the 15 Faqu parent columns ~all 1, the Fpol parents ~all 3.)
##
## Layout: one sector per chromosome (arc proportional to marker count); each
##   individual is a concentric ring, ordered by genome-wide hybrid index so the
##   innermost rings are the most F. aquilonia and the outermost the most
##   F. polyctena. Rendering is done by render_diem_circos() (diem_circos_core.R).
##
## Run from the repo root:  Rscript module_di25/R/diem_circos.R
## =========================================================

suppressMessages(library(data.table))
source("module_di25/R/diem_circos_core.R")

INFILE   <- "data/species_diagnostic_markers_DI25_20pops.tsv.gz"
OUTPNG   <- "module_di25/Figures/diem_circos.png"
DI_LABEL <- -25

message("[diem-circos SNP] reading ", INFILE)
d   <- fread(INFILE)
pos <- d$position
gt  <- as.matrix(d[, !c("chromosome", "position"), with = FALSE])   # markers x individuals
storage.mode(gt) <- "integer"
chr_num <- as.integer(sub("chromosome_", "", d$chromosome))
message(sprintf("  %d markers x %d individuals", nrow(gt), ncol(gt)))

## order markers by chromosome then position
ord_m <- order(chr_num, pos)
gt <- gt[ord_m, ]; chr_num <- chr_num[ord_m]

## order individuals by hybrid index (0 = aquilonia .. 1 = polyctena)
hi <- apply(gt, 2, function(g) { g <- g[g > 0]; if (!length(g)) NA_real_ else mean((g - 1) / 2) })
gt <- gt[, order(hi)]

n_miss <- sum(gt == 0L)
ttl <- sprintf("Hybrids masked @ DI = %d: %s SNPs, %s missing genotypes (%.1f%%) dithered",
               DI_LABEL, format(nrow(gt), big.mark = ","),
               format(n_miss, big.mark = ","), 100 * n_miss / length(gt))

message("[diem-circos SNP] rendering -> ", OUTPNG)
render_diem_circos(gt, chr_num, OUTPNG, title = ttl)
message("[diem-circos SNP] done")
