## =========================================================
## module_manuscript_rho05 -- per-population near-fixation circos, eMLG (rho05)
## =========================================================
## eMLG-level companion to module_di25/R_legacy/di25_population_fixation.R,
## reclustered under min_r2_rho = 0.5. The per-SNP panel does not depend on
## clustering at all and is reused unchanged (module_di25/Figures/di25_popfix_snp.png).
##
## Run from the repo root: Rscript module_manuscript_rho05/R/di25_population_fixation_emlg_rho05.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")
source("module_di25/R/diem_circos_core.R")

PHI <- 0.85
PAL <- c("#F4F4F4", "#21918C", "#D3C93B")
FIGDIR <- "module_manuscript_rho05/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

orient_aqu <- function(M, faqu, fpol) {
  flip <- which(rowMeans(M[, faqu, drop = FALSE], na.rm = TRUE) <
                rowMeans(M[, fpol, drop = FALSE], na.rm = TRUE))
  M[flip, ] <- 2 - M[flip, ]; M
}

inp <- readRDS("module_di25/data/di25_inputs.rds")
e2  <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd  <- e2$sample_data_with_parents
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)
pops    <- sd$Population[match(rownames(GTs_all), sd$Sample_ID)]
faqu <- grep("^Faqu", rownames(GTs_all)); fpol <- grep("^Fpol", rownames(GTs_all))
hybrid_pops <- setdiff(unique(pops[!is.na(pops)]), c("aquilonia_parent", "polyctena_parent"))

popfix <- function(a) {
  F <- vapply(hybrid_pops, function(p) rowMeans(a[, which(pops == p), drop = FALSE], na.rm = TRUE) / 2,
              numeric(nrow(a)))
  ord_p <- order(colMeans(F, na.rm = TRUE), decreasing = TRUE)
  F <- F[, ord_p]
  code <- matrix(0L, nrow(F), ncol(F))
  code[F >= PHI]      <- 1L
  code[F <= (1 - PHI)] <- 2L
  list(code = code, pop_order = hybrid_pops[ord_p])
}

render_level <- function(a_units, chr_num, level, n_lab) {
  pf <- popfix(a_units)
  outpng <- sprintf("%s/di25_popfix_%s.png", FIGDIR, tolower(level))
  frac_pol <- round(100 * mean(pf$code == 2L), 1); frac_aqu <- round(100 * mean(pf$code == 1L), 1)
  message("[popfix-rho05] ", level, ": ", frac_aqu, "% cells aqu-fixed, ", frac_pol, "% pol-fixed -> ", outpng)
  render_circos_raster(
    pf$code, chr_num, palette = PAL, outpng = outpng, ring_sep = FALSE,
    title = sprintf("Per-population near-fixation (phi = %.2f, %s)  |  %s %s x 20 populations",
                    PHI, level, format(nrow(pf$code), big.mark = ","), n_lab),
    ring_labels = pf$pop_order, ring_label_cex = 0.55, open_deg = 30,
    npx = 3800, res = 200, chr_label_cex = 1.5,
    cex_main = 1.8, main_line = 1.0, legend_cex = 1.4,
    legend_labels = c("near-fixed F. aquilonia", "near-fixed F. polyctena", "not near-fixed"),
    legend_cols = PAL[c(2, 3, 1)])
}

## =========================================================================
## per-eMLG (5 cM, min_r2_rho = 0.5)
## =========================================================================
g <- readRDS("module_di25_rho05/data/di25_clustering_cM5_rho05.rds")$groups
is_emlg <- g$n_loci > 2
D <- vapply(seq_len(nrow(g)), function(i) {
  if (is_emlg[i]) consensus_dosage(GTs_all, g$members[[i]]) else GTs_all[, g$representative[i]]
}, numeric(nrow(GTs_all)))
D <- t(D)
rep_chr <- as.integer(sub("Chr", "", sub(":.*", "", g$representative)))
rep_pos <- as.integer(sub(".*:", "", g$representative))
ord_u <- order(rep_chr, rep_pos)
render_level(orient_aqu(D[ord_u, ], faqu, fpol), rep_chr[ord_u], "eMLG_rho05", "LD-reduced units")

message("[popfix-rho05] done (SNP panel unaffected -- see module_di25/Figures/di25_popfix_snp.png)")
