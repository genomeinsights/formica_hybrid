## =========================================================
## module_di25_rho05 -- per-SNP vs LD-reduced (0.2) vs LD-reduced (rho=0.5) DIEM circos
## =========================================================
## Three panels on the SAME DI25 markers and the SAME 195 individuals:
##   (a) per-SNP    -- 51,612 diagnostic SNPs                 (unaffected by min_r2)
##   (b) LD-reduced -- canonical clustering, min_r2 = 0.2      (di25_clustering_cM5.rds)
##   (c) LD-reduced -- min_r2_rho = 0.5                        (di25_clustering_cM5_rho05.rds)
## Individuals are ordered ONCE (by the per-SNP hybrid index) and that order is
## applied to all three panels, so a ring at a given radius is the same
## individual everywhere -- only the LD reduction differs between (b) and (c).
##
## Run from the repo root:  Rscript module_di25_rho05/R/diem_circos_compare_rho05.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")
source("module_di25/R/diem_circos_core.R")

TSV        <- "data/species_diagnostic_markers_DI25_20pops.tsv.gz"
CLUST_02   <- "module_di25/data/di25_clustering_cM5.rds"
CLUST_RHO  <- "module_di25_rho05/data/di25_clustering_cM5_rho05.rds"
INPUTS     <- "module_di25/data/di25_inputs.rds"
OUTPNG     <- "module_di25_rho05/Figures/diem_circos_compare_snp_vs_ldreduced_02_vs_rho05.png"
dir.create(dirname(OUTPNG), showWarnings = FALSE, recursive = TRUE)

## ---- (a) per-SNP: markers x individuals, code 1/2/3 (0 miss) -------------
message("[compare-rho05] per-SNP panel")
d <- fread(TSV)
snp_chr <- as.integer(sub("chromosome_", "", d$chromosome))
snp_gt  <- as.matrix(d[, !c("chromosome", "position"), with = FALSE])   # 1/2/3, 0=miss
storage.mode(snp_gt) <- "integer"
ord_m   <- order(snp_chr, d$position)
snp_gt  <- snp_gt[ord_m, ]; snp_chr <- snp_chr[ord_m]
snp_inds <- colnames(snp_gt)

inp <- readRDS(INPUTS)
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)                              # 195 x markers (012)
faqu <- grep("^Faqu", rownames(GTs_all)); fpol <- grep("^Fpol", rownames(GTs_all))

## shared LD-reduced-panel builder (best-SNP representation, as in
## diem_circos_compare.R / di25_sorting.R)
build_ld_panel <- function(clust_path, tag) {
  message("[compare-rho05] LD-reduced panel (", tag, ")")
  res <- readRDS(clust_path); g <- res$groups
  is_emlg <- g$n_loci > 2
  best <- eMLG_best_snp(res, inp$GTs_hyb, fill = FALSE)
  bm   <- setNames(best$stats$best_marker, best$stats$group_id)
  rep_snp <- g$representative
  rep_snp[is_emlg] <- bm[g$group_id[is_emlg]]
  D <- t(GTs_all[, rep_snp, drop = FALSE])                              # units x individuals
  flip <- which(rowMeans(D[, faqu, drop = FALSE], na.rm = TRUE) <
                rowMeans(D[, fpol, drop = FALSE], na.rm = TRUE))
  D[flip, ] <- 2 - D[flip, ]                                             # 2 = aquilonia
  code <- 3L - round(D); code[is.na(code)] <- 0L
  storage.mode(code) <- "integer"
  colnames(code) <- rownames(GTs_all)
  chr <- as.integer(sub("Chr", "", sub(":.*", "", rep_snp)))
  pos <- as.integer(sub(".*:", "", rep_snp))
  ord <- order(chr, pos)
  list(code = code[ord, ], chr = chr[ord], n_emlg = sum(is_emlg), n_rep = sum(!is_emlg),
       n_units = nrow(g))
}

panel_02  <- build_ld_panel(CLUST_02,  "min_r2=0.2")
panel_rho <- build_ld_panel(CLUST_RHO, "min_r2_rho=0.5")

## ---- shared individual order (by per-SNP hybrid index) ------------------
hi <- apply(snp_gt, 2, function(x) { x <- x[x > 0]; if (!length(x)) NA_real_ else mean((x - 1) / 2) })
ord_names <- snp_inds[order(hi)]                                       # inner = most aquilonia
snp_gt         <- snp_gt[, ord_names]
panel_02$code  <- panel_02$code[, match(ord_names, colnames(panel_02$code))]
panel_rho$code <- panel_rho$code[, match(ord_names, colnames(panel_rho$code))]

## ---- draw all three panels into one wide canvas --------------------------
message("[compare-rho05] rendering -> ", OUTPNG)
png(OUTPNG, width = 7800, height = 3000, res = 300)
par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))
render_diem_circos(snp_gt, snp_chr, new_device = FALSE,
                   title = sprintf("a   Per-SNP: %s diagnostic SNPs", format(nrow(snp_gt), big.mark = ",")),
                   cex_main = 1.0, chr_label_cex = 0.85)
render_diem_circos(panel_02$code, panel_02$chr, new_device = FALSE,
                   title = sprintf("b   LD-reduced (min_r2=0.2): %s units = %s best-SNP + %s rep-SNP",
                                   format(panel_02$n_units, big.mark = ","),
                                   format(panel_02$n_emlg, big.mark = ","),
                                   format(panel_02$n_rep, big.mark = ",")),
                   cex_main = 1.0, chr_label_cex = 0.85)
render_diem_circos(panel_rho$code, panel_rho$chr, new_device = FALSE,
                   title = sprintf("c   LD-reduced (min_r2_rho=0.5): %s units = %s best-SNP + %s rep-SNP",
                                   format(panel_rho$n_units, big.mark = ","),
                                   format(panel_rho$n_emlg, big.mark = ","),
                                   format(panel_rho$n_rep, big.mark = ",")),
                   cex_main = 1.0, chr_label_cex = 0.85)
mtext("DIEM ancestry (purple = F. aquilonia, teal = het, yellow = F. polyctena)  |  rings inner to outer: most aquilonia to most polyctena",
      outer = TRUE, cex = 0.8, line = 0.3)
dev.off()
message("[compare-rho05] done")
