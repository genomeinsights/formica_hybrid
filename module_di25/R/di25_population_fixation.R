## =========================================================
## module_di25 (high-DI analyses) -- per-POPULATION near-fixation circos
## =========================================================
## The population analog of the aggregate sorting sweep: instead of one cross-
## population sort class per unit, show each of the 20 hybrid populations' OWN
## near-fixation direction at every unit (phi = 0.85). This exposes repeatability
## directly -- a region near-fixed for the SAME ancestry across all population
## rings is a repeatably sorted region.
##   purple = population near-fixed for F. aquilonia (oriented freq >= phi)
##   yellow = population near-fixed for F. polyctena  (oriented freq <= 1-phi)
##   grey   = not near-fixed in that population
## Rings inner -> outer are ordered by overall aquilonia ancestry. Rendered at
## both levels (SNP and eMLG 5 cM).
##
## Run from the repo root:  Rscript module_di25/R/di25_population_fixation.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")
source("module_di25/R/diem_circos_core.R")

PHI <- 0.85                                          # near-fixation floor
PAL <- c("#F4F4F4", "#21918C", "#D3C93B")            # code 0 not-fixed / 1 aqu / 2 pol

## orient a dosage matrix (units x individuals) so 2 = F. aquilonia allele
orient_aqu <- function(M, faqu, fpol) {
  flip <- which(rowMeans(M[, faqu, drop = FALSE], na.rm = TRUE) <
                rowMeans(M[, fpol, drop = FALSE], na.rm = TRUE))
  M[flip, ] <- 2 - M[flip, ]; M
}

## ---- inputs -------------------------------------------------------------
inp <- readRDS("module_di25/data/di25_inputs.rds")
e2  <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd  <- e2$sample_data_with_parents
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)
pops    <- sd$Population[match(rownames(GTs_all), sd$Sample_ID)]
faqu <- grep("^Faqu", rownames(GTs_all)); fpol <- grep("^Fpol", rownames(GTs_all))
hybrid_pops <- setdiff(unique(pops[!is.na(pops)]), c("aquilonia_parent", "polyctena_parent"))

## per-population near-fixation code (units x 20 pops), from an oriented matrix
popfix <- function(a) {
  F <- vapply(hybrid_pops, function(p) rowMeans(a[, which(pops == p), drop = FALSE], na.rm = TRUE) / 2,
              numeric(nrow(a)))                      # units x pops, oriented aquilonia freq
  ord_p <- order(colMeans(F, na.rm = TRUE), decreasing = TRUE)  # inner ring = most aquilonia
  F <- F[, ord_p]
  code <- matrix(0L, nrow(F), ncol(F))
  code[F >= PHI]      <- 1L                          # near-fixed aquilonia
  code[F <= (1 - PHI)] <- 2L                         # near-fixed polyctena
  list(code = code, pop_order = hybrid_pops[ord_p])
}

render_level <- function(a_units, chr_num, level, n_lab) {
  pf <- popfix(a_units)
  outpng <- sprintf("module_di25/Figures/di25_popfix_%s.png", tolower(level))
  frac_pol <- round(100 * mean(pf$code == 2L), 1); frac_aqu <- round(100 * mean(pf$code == 1L), 1)
  message("[popfix] ", level, ": ", frac_aqu, "% cells aqu-fixed, ", frac_pol, "% pol-fixed -> ", outpng)
  render_circos_raster(
    pf$code, chr_num, palette = PAL, outpng = outpng, ring_sep = FALSE,
    title = sprintf("Per-population near-fixation (phi = %.2f, %s)  |  %s %s x 20 populations",
                    PHI, level, format(nrow(pf$code), big.mark = ","), n_lab),
    ring_labels = pf$pop_order, ring_label_cex = 0.34, open_deg = 30,
    legend_labels = c("near-fixed F. aquilonia", "near-fixed F. polyctena", "not near-fixed"),
    legend_cols = PAL[c(2, 3, 1)])
}

## =========================================================================
## per-SNP
## =========================================================================
Msnp <- t(GTs_all)
chr_snp <- as.integer(sub("Chr", "", sub(":.*", "", rownames(Msnp))))
pos_snp <- as.integer(sub(".*:", "", rownames(Msnp)))
ord_m <- order(chr_snp, pos_snp)
render_level(orient_aqu(Msnp[ord_m, ], faqu, fpol), chr_snp[ord_m], "SNP", "diagnostic SNPs")

## =========================================================================
## per-eMLG (5 cM)
## =========================================================================
g <- readRDS("module_di25/data/di25_clustering_cM5.rds")$groups
is_emlg <- g$n_loci > 2
D <- vapply(seq_len(nrow(g)), function(i) {
  if (is_emlg[i]) consensus_dosage(GTs_all, g$members[[i]]) else GTs_all[, g$representative[i]]
}, numeric(nrow(GTs_all)))
D <- t(D)
rep_chr <- as.integer(sub("Chr", "", sub(":.*", "", g$representative)))
rep_pos <- as.integer(sub(".*:", "", g$representative))
ord_u <- order(rep_chr, rep_pos)
render_level(orient_aqu(D[ord_u, ], faqu, fpol), rep_chr[ord_u], "eMLG", "LD-reduced units")

message("[popfix] done")
