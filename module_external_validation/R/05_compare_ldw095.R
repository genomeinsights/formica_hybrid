## module_external_validation / 05: how well does the local-LD-support vector ld_w_095
## (median r^2 to neighbours within the rho = 0.95 decay window) agree between the
## original hybrid data (165 inds; map_hyb_005$ld_w_095) and the external hybrids
## (21 inds, MAF>0.1 SNPs; computed here from the slide = 100 decay fit)?
##
## Caveats: (i) different marker sets/densities and different n, so absolute levels are not
## expected to match (r^2 from n = 21 carries a ~1/n upward bias) -- agreement is judged by
## correlation, per SNP, per chromosome, and in physical bins (100 kb); (ii) the window length
## d = rho/(a(1-rho)) uses each fit's own a_pred.
##
## Run from the repo root:  Rscript module_external_validation/R/05_compare_ldw095.R
## Writes: data/ext_ldw095.rds, results/ldw095_comparison.txt, results/ldw095_comparison.png

suppressPackageStartupMessages({ library(data.table); library(SNPRelate); library(ggplot2); library(patchwork) })
devtools::load_all("~/gitlab/LDscnR/", quiet = TRUE)
OUT <- "module_external_validation"; SLIDE <- 100L
ldw_f <- file.path(OUT, "data/ext_ldw095.rds")

d <- readRDS(file.path(OUT, "data/ext_geno_maf10.rds")); emap <- d$map
if (!file.exists(ldw_f)) {
  ld_decay <- readRDS(file.path(OUT, sprintf("data/ld_decay_ext_maf10_slide%d.rds", SLIDE)))
  gds <- snpgdsOpen(file.path(OUT, sprintf("data/ext_maf10_slide%d.gds", SLIDE)))
  w <- compute_ld_w(ld_decay, rho = 0.95, cores = 8, gds = gds, slide_win_ld = SLIDE, ld_method = "corr")
  snpgdsClose(gds); saveRDS(w, ldw_f)
}
w <- readRDS(ldw_f); stopifnot(all(names(w) %in% emap$marker))
emap[, ldw_ext := w[marker]]

e <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e)
omap <- as.data.table(e$map_hyb_005)[, .(marker, Chr, Pos, ldw_orig = ld_w_095, maf_orig = maf_hyb)]
m <- merge(emap[, .(marker, ldw_ext, maf21)], omap, by = "marker")
m <- m[is.finite(ldw_ext) & is.finite(ldw_orig)]
m[, chr_n := as.integer(sub("Chr", "", Chr))]

sink(file.path(OUT, "results/ldw095_comparison.txt"), split = TRUE)
cat("ld_w_095 external: n =", sum(is.finite(emap$ldw_ext)), "SNPs (NA for", sum(!is.finite(emap$ldw_ext)), "); median", round(median(emap$ldw_ext, na.rm = TRUE), 3), "\n")
cat("ld_w_095 original (all):  median", round(median(omap$ldw_orig, na.rm = TRUE), 3), "; on shared SNPs", round(median(m$ldw_orig), 3), "(external on shared:", round(median(m$ldw_ext), 3), ")\n")
cat("shared SNPs with finite ld_w in both:", nrow(m), "\n\n")
cc <- function(x, y) c(pearson = cor(x, y), spearman = cor(x, y, method = "spearman"))
cat("Per-SNP correlation (all shared):\n"); print(round(cc(m$ldw_orig, m$ldw_ext), 4))
cat("\nPer-chromosome per-SNP Spearman:\n")
pc <- m[, .(n = .N, spearman = cor(ldw_orig, ldw_ext, method = "spearman"), pearson = cor(ldw_orig, ldw_ext)), by = chr_n][order(chr_n)]
print(pc[, lapply(.SD, function(x) round(x, 3))])
cat("median per-chromosome Spearman:", round(median(pc$spearman), 3), "\n")
for (bw in c(25e3, 100e3, 500e3)) {
  m[, bin := floor(Pos / bw)]
  b <- m[, .(o = mean(ldw_orig), e = mean(ldw_ext), n = .N), by = .(chr_n, bin)][n >= 5]
  cat(sprintf("\nBinned (%d kb, >=5 shared SNPs/bin): %d bins; ", bw / 1e3, nrow(b))); print(round(cc(b$o, b$e), 4))
}
## does the agreement hold when the original is restricted to the same MAF>0.1 marker set?
m10 <- m[maf_orig >= 0.1]; cat("\nShared SNPs also MAF>=0.1 in the original hybrids (n =", nrow(m10), "):\n"); print(round(cc(m10$ldw_orig, m10$ldw_ext), 4))
sink()

m[, bin := floor(Pos / 100e3)]; b <- m[, .(o = mean(ldw_orig), e = mean(ldw_ext), n = .N), by = .(chr_n, bin)][n >= 5]
p1 <- ggplot(m[sample(.N, min(.N, 50000))], aes(ldw_orig, ldw_ext)) + geom_bin2d(bins = 60) + scale_fill_viridis_c(trans = "log10") +
  labs(x = "ld_w_095 original", y = "ld_w_095 external", title = sprintf("per SNP (Spearman %.2f)", cor(m$ldw_orig, m$ldw_ext, method = "spearman"))) + theme_classic()
p2 <- ggplot(b, aes(o, e)) + geom_point(alpha = .3, size = .8) + geom_smooth(method = "lm", se = FALSE, colour = "firebrick") +
  labs(x = "mean ld_w_095 original", y = "mean ld_w_095 external", title = sprintf("100 kb bins (Spearman %.2f)", cor(b$o, b$e, method = "spearman"))) + theme_classic()
ggsave(file.path(OUT, "results/ldw095_comparison.png"), p1 + p2, width = 11, height = 4.8, dpi = 200)
