## =========================================================================
## module_allele_specific_sorting -- 02: F_ST vs LD along a distance axis (A1)
##
## Per distance bin (bp; cM as sensitivity) plus an "unlinked" category
## (cross-chromosome pairs), chromosome-block bootstrap (resampling the 26
## chromosomes; cross-chromosome pairs resampled as pairs of drawn chromosomes):
##   ceiling      mean G_i G_j        among-population LD^2 the two loci would
##                                    show if they partitioned populations identically
##   among LD     mean rST^2          observed among-population LD^2
##   realised     sum rST^2 / sum G_i G_j   fraction of the ceiling realised
##   within LD    mean r_w^2, r_w_adj^2
##   concordance  mean conc, conc_resid (signed)
##
## Expected under heterogeneous, allele-specific sorting: the ceiling is high
## and flat with distance (it only depends on the loci's own F_ST), while the
## realised fraction is substantial only over linkage distances and ~0 for
## unlinked pairs; within-population LD decays like neutral LD.
##
## Inputs: data/01_agg_within.rds, data/01_agg_cross.rds
## Outputs: data/02_decay.rds, Figures/02_fst_vs_ld_decay.{png,pdf},
##          Figures/02_fst_vs_ld_decay_by_fst_class.png, Figures/02_decay_cM.png
## Run from the formica_hybrid repo root.
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
B <- 2000L; SEED <- 1L

aw <- readRDS(file.path(OUT_DATA, "01_agg_within.rds"))
ac <- readRDS(file.path(OUT_DATA, "01_agg_cross.rds"))
chrs <- sort(unique(c(ac$ChrA, ac$ChrB)))
set.seed(SEED)
draws <- replicate(B, tabulate(match(sample(chrs, length(chrs), replace = TRUE), chrs), length(chrs)),
                   simplify = FALSE)

## the quantities: name -> (numerator column, denominator column)
QTY <- list(ceiling    = c("s_ceil2",    "n_ceil2"),
            among_LD   = c("s_rST2",     "n_rST2"),
            realised   = c("s_rST2",     "s_ceil2"),
            within_LD  = c("s_r_w2",     "n_r_w2"),
            within_LD_adj = c("s_r_w_adj2", "n_r_w_adj2"),
            conc       = c("s_conc",     "n_conc"),
            conc_resid = c("s_conc_resid", "n_conc_resid"))

## within-chromosome bins: resample chromosomes (draw counts as weights)
boot_within <- function(agg, bin_col, qty, classes = 1:3) {
  a <- agg[fst_class %in% classes]
  rbindlist(lapply(names(qty), function(q) {
    num <- qty[[q]][1]; den <- qty[[q]][2]
    a[, {
      s <- tapply(get(num), Chr, sum)[chrs]; n <- tapply(get(den), Chr, sum)[chrs]
      s[is.na(s)] <- 0; n[is.na(n)] <- 0
      bt <- vapply(draws, function(w) sum(w * s) / sum(w * n), numeric(1))
      list(qty = q, mean = sum(s) / sum(n), lo = quantile(bt, 0.025, na.rm = TRUE),
           hi = quantile(bt, 0.975, na.rm = TRUE))
    }, by = bin_col] |> setnames(bin_col, "bin")
  }))
}

## cross-chromosome: pair (a,b) weight = w_a * w_b (copies of one chromosome are not cross pairs)
boot_cross <- function(agg, qty, classes = 1:3) {
  a <- agg[fst_class %in% classes]
  rbindlist(lapply(names(qty), function(q) {
    num <- qty[[q]][1]; den <- qty[[q]][2]
    S <- N <- matrix(0, length(chrs), length(chrs), dimnames = list(chrs, chrs))
    s <- a[, .(v = sum(get(num))), by = .(ChrA, ChrB)]; n <- a[, .(v = sum(get(den))), by = .(ChrA, ChrB)]
    S[cbind(s$ChrA, s$ChrB)] <- s$v; N[cbind(n$ChrA, n$ChrB)] <- n$v
    S <- S + t(S); N <- N + t(N)                 # symmetric, zero diagonal
    bt <- vapply(draws, function(w) { W <- outer(w, w); sum(W * S) / sum(W * N) }, numeric(1))
    data.table(bin = "unlinked", qty = q, mean = sum(S) / sum(N),
               lo = quantile(bt, 0.025, na.rm = TRUE), hi = quantile(bt, 0.975, na.rm = TRUE))
  }))
}

run_all <- function(classes, label) {
  d_bp <- rbind(boot_within(aw$bp, "bp_bin", QTY, classes), boot_cross(ac, QTY, classes))
  d_cm <- rbind(boot_within(aw$cm, "cm_bin", QTY, classes), boot_cross(ac, QTY, classes))
  d_bp[, bin := factor(bin, levels = c(BP_LABELS, "unlinked"))]
  d_cm[, bin := factor(bin, levels = c(CM_LABELS, "unlinked"))]
  d_bp[, fst_set := label]; d_cm[, fst_set := label]
  list(bp = d_bp, cm = d_cm)
}

res_all <- run_all(1:3, "all pairs")
res_cls <- lapply(1:3, function(k) run_all(k, c("min F_ST: low", "min F_ST: mid", "min F_ST: high")[k]))
by_class_bp <- rbindlist(c(list(res_all$bp), lapply(res_cls, `[[`, "bp")))

cat("\n[02] all pairs, by physical distance (mean [95% block-bootstrap CI]):\n")
print(dcast(res_all$bp[, .(bin, qty, v = sprintf("%.4f [%.4f,%.4f]", mean, lo, hi))], bin ~ qty, value.var = "v"))

saveRDS(list(bp = res_all$bp, cm = res_all$cm, by_class_bp = by_class_bp, B = B,
             fst_breaks = aw$fst_breaks),
        file.path(OUT_DATA, "02_decay.rds"))

## ---- figures -----------------------------------------------------------------
## r_w^2 and r_w_adj^2 are near-identical at all distances (checked in the printed table), so only
## the hybrid-index-adjusted one is drawn
ld_lab <- c(ceiling = "ceiling: G_i G_j (shared partition)", among_LD = "among-population LD (rST^2)",
            within_LD_adj = "within-population LD (r_w^2, hybrid-index adj.)")
N_POP <- 20L   # chance floor for the realised fraction: E[r^2] of two unrelated 20-population profiles = 1/(N_POP-1)
plot_ld <- function(d) {
  ggplot(d[qty %in% names(ld_lab)], aes(bin, mean, colour = qty, group = qty)) +
    geom_line() + geom_point(size = 1.4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
    scale_y_log10() + scale_colour_brewer(palette = "Dark2", labels = ld_lab, name = NULL) +
    labs(x = NULL, y = "mean per pair (log scale)") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom") +
    guides(colour = guide_legend(ncol = 2))
}
plot_conc <- function(d) {
  dd <- d[qty %in% c("realised", "conc", "conc_resid")]
  dd[, qty := factor(qty, levels = c("realised", "conc", "conc_resid"),
                     labels = c("fraction of ceiling realised", "concordance (signed r)", "concordance, ancestry-residualised"))]
  ggplot(dd, aes(bin, mean, colour = qty, group = qty)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_hline(yintercept = 1 / (N_POP - 1), colour = "#e41a1c", linetype = 3) +
    geom_line() + geom_point(size = 1.4) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
    scale_colour_brewer(palette = "Set1", name = NULL) +
    labs(x = NULL, y = "mean per pair", caption = "dotted: chance level of the realised fraction, 1/19 (20 populations)") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom") +
    guides(colour = guide_legend(ncol = 1))
}
if (requireNamespace("patchwork", quietly = TRUE)) {
  p <- patchwork::wrap_plots(plot_ld(res_all$bp) + ggtitle("a  LD components vs distance"),
                             plot_conc(res_all$bp) + ggtitle("b  shared partitions vs distance"), nrow = 1)
  ggsave(file.path(OUT_FIG, "02_fst_vs_ld_decay.png"), p, width = 11, height = 4.8, dpi = 200)
  ggsave(file.path(OUT_FIG, "02_fst_vs_ld_decay.pdf"), p, width = 11, height = 4.8)
  pc <- patchwork::wrap_plots(plot_ld(res_all$cm) + ggtitle("a  vs genetic distance"),
                              plot_conc(res_all$cm) + ggtitle("b"), nrow = 1)
  ggsave(file.path(OUT_FIG, "02_decay_cM.png"), pc, width = 11, height = 4.8, dpi = 200)
}
pf <- plot_ld(by_class_bp) + facet_wrap(~ fst_set, nrow = 1)
ggsave(file.path(OUT_FIG, "02_fst_vs_ld_decay_by_fst_class.png"), pf, width = 13, height = 4.5, dpi = 200)
cat("[02] done\n")
