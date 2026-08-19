## =========================================================
## module_di25 (high-DI analyses) -- per-SNP vs LD-reduced DIEM circos (side by side)
## =========================================================
## One figure, two panels on the SAME DI25 markers and the SAME 195 individuals:
##   (a) per-SNP        -- 51,612 diagnostic SNPs        (diem_circos.R)
##   (b) LD-reduced     -- 11,052 high-DI-only units, 5cM (diem_circos_di25_eMLG.R)
## Individuals are ordered ONCE (by the per-SNP hybrid index) and that order is
## applied to both panels, so a ring at a given radius is the same individual in
## both -- the LD reduction is the only thing that changes between panels.
##
## Run from the repo root:  Rscript module_di25/R/diem_circos_compare.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")
source("module_di25/R/diem_circos_core.R")

TSV      <- "data/species_diagnostic_markers_DI25_20pops.tsv.gz"
CLUST    <- "module_di25/data/di25_clustering_cM5.rds"
INPUTS   <- "module_di25/data/di25_inputs.rds"
OUTPNG   <- "module_di25/Figures/diem_circos_compare_snp_vs_ldreduced.png"

## ---- (a) per-SNP: markers x individuals, code 1/2/3 (0 miss) -------------
message("[compare] per-SNP panel")
d <- fread(TSV)
snp_chr <- as.integer(sub("chromosome_", "", d$chromosome))
snp_gt  <- as.matrix(d[, !c("chromosome", "position"), with = FALSE])   # 1/2/3, 0=miss
storage.mode(snp_gt) <- "integer"
ord_m   <- order(snp_chr, d$position)
snp_gt  <- snp_gt[ord_m, ]; snp_chr <- snp_chr[ord_m]
snp_inds <- colnames(snp_gt)

## ---- (b) LD-reduced: units x individuals, code 1/2/3 --------------------
message("[compare] LD-reduced panel")
res <- readRDS(CLUST); g <- res$groups
inp <- readRDS(INPUTS)
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)                              # 195 x markers (012)
faqu <- grep("^Faqu", rownames(GTs_all)); fpol <- grep("^Fpol", rownames(GTs_all))
is_emlg <- g$n_loci > 2
## BEST-SNP representation (consistent with di25_sorting.R): each unit is one real
## SNP -- the consensus-optimal member (best_marker) for >2-marker clusters, the
## representative for 1-2-marker clusters -- not the averaged consensus dosage.
best <- eMLG_best_snp(res, inp$GTs_hyb, fill = FALSE)
bm   <- setNames(best$stats$best_marker, best$stats$group_id)
rep_snp <- g$representative
rep_snp[is_emlg] <- bm[g$group_id[is_emlg]]
D <- t(GTs_all[, rep_snp, drop = FALSE])                              # units x individuals (real best-SNP)
flip <- which(rowMeans(D[, faqu, drop = FALSE], na.rm = TRUE) <
              rowMeans(D[, fpol, drop = FALSE], na.rm = TRUE))
D[flip, ] <- 2 - D[flip, ]                                             # 2 = aquilonia
emlg_code <- 3L - round(D); emlg_code[is.na(emlg_code)] <- 0L
storage.mode(emlg_code) <- "integer"
colnames(emlg_code) <- rownames(GTs_all)
rep_chr <- as.integer(sub("Chr", "", sub(":.*", "", rep_snp)))
rep_pos <- as.integer(sub(".*:", "", rep_snp))
ord_u   <- order(rep_chr, rep_pos)
emlg_code <- emlg_code[ord_u, ]; emlg_chr <- rep_chr[ord_u]

## ---- shared individual order (by per-SNP hybrid index) ------------------
hi <- apply(snp_gt, 2, function(x) { x <- x[x > 0]; if (!length(x)) NA_real_ else mean((x - 1) / 2) })
ord_names <- snp_inds[order(hi)]                                       # inner = most aquilonia
snp_gt    <- snp_gt[, ord_names]
emlg_code <- emlg_code[, match(ord_names, colnames(emlg_code))]

## ---- draw both panels into one wide canvas ------------------------------
message("[compare] rendering -> ", OUTPNG)
png(OUTPNG, width = 5400, height = 3000, res = 300)
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
render_diem_circos(snp_gt, snp_chr, new_device = FALSE,
                   title = sprintf("a   Per-SNP: %s diagnostic SNPs", format(nrow(snp_gt), big.mark = ",")),
                   cex_main = 1.0, chr_label_cex = 0.85)
render_diem_circos(emlg_code, emlg_chr, new_device = FALSE,
                   title = sprintf("b   LD-reduced: %s units (5 cM) = %s best-SNP + %s rep-SNP",
                                   format(nrow(emlg_code), big.mark = ","),
                                   format(sum(is_emlg), big.mark = ","),
                                   format(sum(!is_emlg), big.mark = ",")),
                   cex_main = 1.0, chr_label_cex = 0.85)
mtext("DIEM ancestry (purple = F. aquilonia, teal = het, yellow = F. polyctena)  |  rings inner to outer: most aquilonia to most polyctena",
      outer = TRUE, cex = 0.8, line = 0.3)
dev.off()
message("[compare] done")
