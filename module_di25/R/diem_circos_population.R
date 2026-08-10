## =========================================================
## module_di25 (high-DI analyses) -- per-POPULATION DIEM circos (both levels)
## =========================================================
## The DIEM circos summarised per population: the 195 individual rings collapse
## to 20 hybrid-population rings, each cell = that population's MEAN F. aquilonia-
## allele frequency at the unit (continuous). Rendered at both levels:
##   SNP  -- 51,612 diagnostic markers
##   eMLG -- 11,052 LD-reduced units (5 cM)
##
## Colour (viridis): purple = population near-fixed for F. aquilonia (freq ~1),
## yellow = near-fixed for F. polyctena (freq ~0), teal ~ balanced. Continuous
## frequency is rendered by binning into a 100-step palette via render_circos_raster().
## Rings inner -> outer are ordered by the population's overall aquilonia ancestry.
##
## Run from the repo root:  Rscript module_di25/R/diem_circos_population.R
## =========================================================

suppressMessages({ library(data.table); library(viridisLite) })
devtools::load_all("~/gitlab/LDscnR/")
source("module_di25/R/diem_circos_core.R")

## continuous species ramp: code 0 = missing (grey); codes 1..100 = freq 1..0
## green (F. aquilonia, freq 1) -> teal (balanced) -> yellow (F. polyctena, freq 0)
GRAMP <- colorRampPalette(c("#21918C", "#D3C93B"))(100)
PAL <- c("#D9D9D9", GRAMP)                     # palette[code+1]; code1=green(freq1)..code100=yellow(freq0)
freq_to_code <- function(f) { v <- 1L + as.integer(round((1 - f) * 99)); v[is.na(f)] <- 0L; v }

## orient a dosage matrix (units x individuals) so 2 = F. aquilonia allele
orient_aqu <- function(M, faqu, fpol) {
  flip <- which(rowMeans(M[, faqu, drop = FALSE], na.rm = TRUE) <
                rowMeans(M[, fpol, drop = FALSE], na.rm = TRUE))
  M[flip, ] <- 2 - M[flip, ]
  M
}

## ---- inputs -------------------------------------------------------------
inp <- readRDS("module_di25/data/di25_inputs.rds")
e2  <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd  <- e2$sample_data_with_parents
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)                       # 195 x markers
pops    <- sd$Population[match(rownames(GTs_all), sd$Sample_ID)] # NA for the unmatched hybrid
faqu <- grep("^Faqu", rownames(GTs_all)); fpol <- grep("^Fpol", rownames(GTs_all))
hybrid_pops <- setdiff(unique(pops[!is.na(pops)]), c("aquilonia_parent", "polyctena_parent"))

## population-mean aquilonia frequency (units x 20 pops) from an oriented matrix
pop_freq <- function(a) {
  vapply(hybrid_pops, function(p) rowMeans(a[, which(pops == p), drop = FALSE], na.rm = TRUE) / 2,
         numeric(nrow(a)))                                       # units x pops
}

render_pop <- function(a_units, chr_num, level, n_lab) {
  F <- pop_freq(a_units)                                         # units x 20
  ord_p <- order(colMeans(F, na.rm = TRUE), decreasing = TRUE)  # inner ring = most aquilonia
  F <- F[, ord_p]; pop_order <- hybrid_pops[ord_p]
  code <- matrix(freq_to_code(F), nrow = nrow(F))
  outpng <- sprintf("module_di25/Figures/diem_circos_population_%s.png", tolower(level))
  message("[pop-circos] rendering ", level, " -> ", outpng)
  render_circos_raster(
    code, chr_num, palette = PAL, outpng = outpng, ring_sep = FALSE,
    title = sprintf("Per-population ancestry (%s): mean F. aquilonia-allele freq  |  %s %s",
                    level, format(nrow(code), big.mark = ","), n_lab),
    ring_labels = pop_order, ring_label_cex = 0.34, open_deg = 30,
    legend_labels = c("F. aquilonia (1.0)", "balanced (0.5)", "F. polyctena (0.0)"),
    legend_cols = c(GRAMP[1], GRAMP[50], GRAMP[100]))
}

## =========================================================================
## per-SNP
## =========================================================================
Msnp <- t(GTs_all)                                              # markers x individuals
chr_snp <- as.integer(sub("Chr", "", sub(":.*", "", rownames(Msnp))))
pos_snp <- as.integer(sub(".*:", "", rownames(Msnp)))
ord_m <- order(chr_snp, pos_snp)
a_snp <- orient_aqu(Msnp[ord_m, ], faqu, fpol)
render_pop(a_snp, chr_snp[ord_m], "SNP", "diagnostic SNPs")

## =========================================================================
## per-eMLG (5 cM): consensus for >2-marker clusters, representative otherwise
## =========================================================================
g <- readRDS("module_di25/data/di25_clustering_cM5.rds")$groups
is_emlg <- g$n_loci > 2
D <- vapply(seq_len(nrow(g)), function(i) {
  if (is_emlg[i]) consensus_dosage(GTs_all, g$members[[i]]) else GTs_all[, g$representative[i]]
}, numeric(nrow(GTs_all)))
D <- t(D)                                                       # units x individuals
rep_chr <- as.integer(sub("Chr", "", sub(":.*", "", g$representative)))
rep_pos <- as.integer(sub(".*:", "", g$representative))
ord_u <- order(rep_chr, rep_pos)
a_emlg <- orient_aqu(D[ord_u, ], faqu, fpol)
render_pop(a_emlg, rep_chr[ord_u], "eMLG", "LD-reduced units")

message("[pop-circos] done")
