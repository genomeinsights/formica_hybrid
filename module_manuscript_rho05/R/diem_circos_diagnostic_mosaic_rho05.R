## =========================================================
## module_manuscript_rho05 -- Fig [diagnostic_mosaic]: per-SNP vs LD-reduced
## (min_r2_rho = 0.5) DIEM circos, two panels (matches the manuscript figure)
## =========================================================
## Adapted from module_di25/R/diem_circos_compare.R, substituting the
## min_r2_rho = 0.5 clustering (module_di25_rho05/data/di25_clustering_cM5_rho05.rds)
## for the canonical min_r2 = 0.2 one. Best-SNP representation, as in
## di25_sorting.R. Two panels only (a: per-SNP, b: LD-reduced), matching the
## manuscript's Fig [diagnostic_mosaic] structure -- see
## diem_circos_compare_rho05.R (module_di25_rho05/R/) for the 3-panel
## 0.2-vs-rho05 internal QA comparison this was built from.
##
## Run from the repo root:  Rscript module_manuscript_rho05/R/diem_circos_diagnostic_mosaic_rho05.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")
source("module_di25/R/diem_circos_core.R")

TSV      <- "data/species_diagnostic_markers_DI25_20pops.tsv.gz"
CLUST    <- "module_di25_rho05/data/di25_clustering_cM5_rho05.rds"
INPUTS   <- "module_di25/data/di25_inputs.rds"
OUTPNG   <- "module_manuscript_rho05/Figures/diem_circos_diagnostic_mosaic_rho05.png"
dir.create(dirname(OUTPNG), showWarnings = FALSE, recursive = TRUE)

## ---- (a) per-SNP: markers x individuals, code 1/2/3 (0 miss) -------------
message("[diagnostic-mosaic-rho05] per-SNP panel")
d <- fread(TSV)
snp_chr <- as.integer(sub("chromosome_", "", d$chromosome))
snp_gt  <- as.matrix(d[, !c("chromosome", "position"), with = FALSE])
storage.mode(snp_gt) <- "integer"
ord_m   <- order(snp_chr, d$position)
snp_gt  <- snp_gt[ord_m, ]; snp_chr <- snp_chr[ord_m]
snp_inds <- colnames(snp_gt)

## ---- (b) LD-reduced: units x individuals, code 1/2/3 (min_r2_rho = 0.5) --
message("[diagnostic-mosaic-rho05] LD-reduced panel")
res <- readRDS(CLUST); g <- res$groups
inp <- readRDS(INPUTS)
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)
faqu <- grep("^Faqu", rownames(GTs_all)); fpol <- grep("^Fpol", rownames(GTs_all))
is_emlg <- g$n_loci > 2
best <- eMLG_best_snp(res, inp$GTs_hyb, fill = FALSE)
bm   <- setNames(best$stats$best_marker, best$stats$group_id)
rep_snp <- g$representative
rep_snp[is_emlg] <- bm[g$group_id[is_emlg]]
D <- t(GTs_all[, rep_snp, drop = FALSE])
flip <- which(rowMeans(D[, faqu, drop = FALSE], na.rm = TRUE) <
              rowMeans(D[, fpol, drop = FALSE], na.rm = TRUE))
D[flip, ] <- 2 - D[flip, ]
emlg_code <- 3L - round(D); emlg_code[is.na(emlg_code)] <- 0L
storage.mode(emlg_code) <- "integer"
colnames(emlg_code) <- rownames(GTs_all)
rep_chr <- as.integer(sub("Chr", "", sub(":.*", "", rep_snp)))
rep_pos <- as.integer(sub(".*:", "", rep_snp))
ord_u   <- order(rep_chr, rep_pos)
emlg_code <- emlg_code[ord_u, ]; emlg_chr <- rep_chr[ord_u]

## ---- shared individual order (by per-SNP hybrid index) ------------------
hi <- apply(snp_gt, 2, function(x) { x <- x[x > 0]; if (!length(x)) NA_real_ else mean((x - 1) / 2) })
ord_names <- snp_inds[order(hi)]
snp_gt    <- snp_gt[, ord_names]
emlg_code <- emlg_code[, match(ord_names, colnames(emlg_code))]

## ---- draw both panels ------------------------------------------------
message("[diagnostic-mosaic-rho05] rendering -> ", OUTPNG)
png(OUTPNG, width = 5400, height = 3000, res = 300)
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
render_diem_circos(snp_gt, snp_chr, new_device = FALSE,
                   title = sprintf("a   Per-SNP: %s diagnostic SNPs", format(nrow(snp_gt), big.mark = ",")),
                   cex_main = 1.0, chr_label_cex = 0.85)
render_diem_circos(emlg_code, emlg_chr, new_device = FALSE,
                   title = sprintf("b   LD-reduced (min_r2_rho=0.5): %s units = %s best-SNP + %s rep-SNP",
                                   format(nrow(emlg_code), big.mark = ","),
                                   format(sum(is_emlg), big.mark = ","),
                                   format(sum(!is_emlg), big.mark = ",")),
                   cex_main = 1.0, chr_label_cex = 0.85)
mtext("DIEM ancestry (purple = F. aquilonia, teal = het, yellow = F. polyctena)  |  rings inner to outer: most aquilonia to most polyctena",
      outer = TRUE, cex = 0.8, line = 0.3)
dev.off()
message("[diagnostic-mosaic-rho05] done")
