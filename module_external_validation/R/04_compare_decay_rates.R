## module_external_validation / 04: compare per-chromosome LD-decay rates between the
## original hybrid data (165 inds, 1.1M SNPs) and the external hybrids (21 inds,
## MAF>0.1 SNPs, slide = 100 fit).
##   raw       = per-chromosome fitted decay rate `a` (decay_sum$a)
##   projected = `a_pred`, the rate predicted from chromosome size by the genome-wide
##               robust regression (decay_sum$a_pred) -- carries only the size trend
## Also compares the asymptote c / c_pred and chromosome size for context.
##
## Run from the repo root:  Rscript module_external_validation/R/04_compare_decay_rates.R
## Writes: results/decay_rate_comparison.tsv, results/decay_rate_comparison.txt, results/decay_rate_comparison.png

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
OUT <- "module_external_validation"
o <- as.data.table(readRDS("module0_ld_pruning/data/ld_decay_DIEM_100w.rds")$decay_sum)
e <- as.data.table(readRDS(file.path(OUT, "data/ld_decay_ext_maf10_slide100.rds"))$decay_sum)
stopifnot("Chr" %in% names(o), "Chr" %in% names(e))
m <- merge(o[, .(Chr, size_o = chr_size, a_o = a, apred_o = a_pred, c_o = c)],
           e[, .(Chr, size_e = chr_size, a_e = a, apred_e = a_pred, c_e = c)], by = "Chr")
m[, chr_n := as.integer(sub("Chr", "", Chr))]; setorder(m, chr_n)
fwrite(m, file.path(OUT, "results/decay_rate_comparison.tsv"), sep = "\t")

ct <- function(x, y) { p <- cor.test(x, y); s <- suppressWarnings(cor.test(x, y, method = "spearman"))
  c(pearson = unname(p$estimate), p_pearson = p$p.value, spearman = unname(s$estimate), p_spearman = s$p.value) }
res <- rbind(raw_a = ct(m$a_o, m$a_e), projected_a_pred = ct(m$apred_o, m$apred_e),
             raw_c = ct(m$c_o, m$c_e), chr_size = ct(m$size_o, m$size_e),
             orig_a_vs_logsize = ct(m$a_o, log(m$size_o)), ext_a_vs_logsize = ct(m$a_e, log(m$size_e)))
sink(file.path(OUT, "results/decay_rate_comparison.txt"), split = TRUE)
cat("Per-chromosome decay: original vs external (n =", nrow(m), "chromosomes)\n"); print(round(res, 4))
cat("\nMedian a: original", signif(median(m$a_o), 3), " external", signif(median(m$a_e), 3), " ratio ext/orig", round(median(m$a_e / m$a_o), 3), "\n")
cat("Median a_pred: original", signif(median(m$apred_o), 3), " external", signif(median(m$apred_e), 3), "\n")
sink()

pl <- function(x, y, lab) {
  r <- cor(m[[x]], m[[y]]); ggplot(m, aes(.data[[x]], .data[[y]])) + geom_abline(linetype = 2, colour = "grey60") +
    geom_point() + geom_text(aes(label = chr_n), nudge_y = diff(range(m[[y]])) * 0.04, size = 3) +
    labs(x = paste(lab, "- original"), y = paste(lab, "- external"), title = sprintf("%s (r = %.2f)", lab, r)) + theme_classic()
}
ggsave(file.path(OUT, "results/decay_rate_comparison.png"),
       pl("a_o", "a_e", "raw a") + pl("apred_o", "apred_e", "projected a_pred") + pl("c_o", "c_e", "asymptote c"),
       width = 13, height = 4.5, dpi = 200)
