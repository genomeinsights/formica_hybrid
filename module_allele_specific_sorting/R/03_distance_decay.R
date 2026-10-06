## =========================================================================
## module_allele_specific_sorting -- 03: differentiation potential vs shared
## sorting along a distance axis (headline analysis)
##
## Per distance bin (genetic distance = headline; physical distance =
## sensitivity) plus an "unlinked" category (all cross-chromosome pairs):
##   ceiling       mean G_i G_j    among-population LD^2 two loci with their
##                                 observed differentiation would show if they
##                                 partitioned the populations identically
##   among LD      mean rST^2      observed among-population LD^2
##   realised      sum rST^2 / sum G_i G_j, with its population-label
##                 permutation baseline b (02) and the excess over it
##   realised_share (realised - b) / (1 - b): the share of the ceiling, beyond
##                 what 20 unrelated population profiles give by chance, that
##                 is actually realised (0 = no shared partition beyond chance,
##                 1 = identical partitions). Interval: bootstrap interval of
##                 the realised fraction transformed with b fixed (b's own
##                 permutation interval is negligible, e.g. +-1e-5 unlinked).
##   within LD     mean r_w^2, hybrid-index adjusted (the only within-population
##                 statistic reported: the unadjusted one is numerically identical)
##   concordance   mean signed c, raw and ancestry-residualised, with baseline
## Uncertainty: chromosome-block bootstrap (2,000 resamples of the 26
## chromosomes; cross-chromosome pairs weighted by the product of the two
## chromosomes' draw counts, so copies of one chromosome are never "unlinked").
##
## Inputs : data/01_agg_within.rds, data/01_agg_cross.rds, data/02_baseline.rds
## Outputs: data/03_decay.rds, Figures/03_decay_cM.{png,pdf} (headline),
##          Figures/03_decay_bp.{png,pdf} (sensitivity)
## Run from the formica_hybrid repo root.
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
B <- 2000L; SEED <- 1L

aw <- readRDS(file.path(OUT_DATA, "01_agg_within.rds"))
ac <- readRDS(file.path(OUT_DATA, "01_agg_cross.rds"))
bl <- readRDS(file.path(OUT_DATA, "02_baseline.rds"))
chrs <- sort(unique(c(ac$ChrA, ac$ChrB)))
draws <- chrom_draws(chrs, B, SEED)

QTY <- list(ceiling       = c("s_ceil2",      "n_ceil2"),
            among_LD      = c("s_rST2",       "n_rST2"),
            realised      = c("s_rST2",       "s_ceil2"),
            within_LD_adj = c("s_r_w_adj2",   "n_r_w_adj2"),
            conc          = c("s_conc",       "n_conc"),
            conc_resid    = c("s_conc_resid", "n_conc_resid"))

boot_within <- function(agg, bin_col) {
  rbindlist(lapply(names(QTY), function(q) {
    num <- QTY[[q]][1]; den <- QTY[[q]][2]
    agg[, {
      s <- tapply(get(num), Chr, sum)[chrs]; n <- tapply(get(den), Chr, sum)[chrs]
      s[is.na(s)] <- 0; n[is.na(n)] <- 0
      bt <- vapply(draws, function(w) sum(w * s) / sum(w * n), numeric(1))
      list(qty = q, mean = sum(s) / sum(n), lo = quantile(bt, 0.025), hi = quantile(bt, 0.975))
    }, by = bin_col] |> setnames(bin_col, "bin")
  }))
}
boot_cross <- function(agg) {
  rbindlist(lapply(names(QTY), function(q) {
    num <- QTY[[q]][1]; den <- QTY[[q]][2]
    S <- N <- matrix(0, length(chrs), length(chrs), dimnames = list(chrs, chrs))
    s <- agg[, .(v = sum(get(num))), by = .(ChrA, ChrB)]; n <- agg[, .(v = sum(get(den))), by = .(ChrA, ChrB)]
    S[cbind(s$ChrA, s$ChrB)] <- s$v; N[cbind(n$ChrA, n$ChrB)] <- n$v
    S <- S + t(S); N <- N + t(N)
    bt <- vapply(draws, function(w) { W <- outer(w, w); sum(W * S) / sum(W * N) }, numeric(1))
    data.table(bin = "unlinked", qty = q, mean = sum(S) / sum(N),
               lo = quantile(bt, 0.025), hi = quantile(bt, 0.975))
  }))
}

decay <- function(agg, bin_col, labels, base) {
  d <- rbind(boot_within(agg, bin_col), boot_cross(ac))
  d[, bin := factor(as.character(bin), levels = c(labels, "unlinked"))]
  ## attach the permutation baseline for realised / conc / conc_resid
  b <- melt(base[, .(bin = as.character(bin),
                     realised_null, realised_lo, realised_hi, conc_null, conc_lo, conc_hi,
                     conc_resid_null, conc_resid_lo, conc_resid_hi)], id.vars = "bin")
  b[, qty := sub("_(null|lo|hi)$", "", variable)][, part := sub(".*_", "", variable)]
  b <- dcast(b, bin + qty ~ part, value.var = "value")
  setnames(b, c("null", "lo", "hi"), c("null_mean", "null_lo", "null_hi"))
  d <- merge(d, b[, bin := factor(bin, levels = c(labels, "unlinked"))], by = c("bin", "qty"), all.x = TRUE)
  d[, excess := mean - null_mean]
  sh <- d[qty == "realised"][, `:=`(qty = "realised_share",
                                    mean = (mean - null_mean) / (1 - null_mean),
                                    lo = (lo - null_mean) / (1 - null_mean),
                                    hi = (hi - null_mean) / (1 - null_mean),
                                    null_mean = NA_real_, null_lo = NA_real_, null_hi = NA_real_, excess = NA_real_)]
  d <- rbind(d, sh)
  setorder(d, qty, bin); d
}
res_cm <- decay(aw$cm, "cm_bin", CM_LABELS, bl$cm)
res_bp <- decay(aw$bp, "bp_bin", BP_LABELS, bl$bp)

cat("\n[03] genetic distance (headline): mean [95% block-bootstrap CI]; permutation baseline; excess\n")
print(res_cm[qty %in% c("realised", "realised_share", "conc_resid", "within_LD_adj"),
             .(qty, bin, obs = sprintf("%.4f [%.4f,%.4f]", mean, lo, hi),
               baseline = ifelse(is.na(null_mean), "", sprintf("%.4f [%.4f,%.4f]", null_mean, null_lo, null_hi)),
               excess = ifelse(is.na(excess), "", sprintf("%.4f", excess)))])
saveRDS(list(cm = res_cm, bp = res_bp, B = B, n_perm = bl$n_perm), file.path(OUT_DATA, "03_decay.rds"))

## ---- figures -------------------------------------------------------------------
ld_lab <- c(ceiling = "ceiling: maximum among-population LD",
            among_LD = "observed among-population LD",
            within_LD_adj = "within-population LD")
plot_ld <- function(d, xlab) {
  ggplot(d[qty %in% names(ld_lab)], aes(bin, mean, colour = qty, group = qty)) +
    geom_line() + geom_point(size = 1.4) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
    scale_y_log10() +
    scale_colour_manual(values = c(ceiling = "#d95f02", among_LD = "#1b9e77", within_LD_adj = "#7570b3"),
                        labels = ld_lab, name = NULL) +
    labs(x = xlab, y = expression("mean per pair ("*r^2*" scale, log)")) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom") +
    guides(colour = guide_legend(ncol = 1))
}
plot_shared <- function(d, xlab) {
  dd <- d[qty %in% c("realised", "conc_resid")]
  dd[, qty := factor(qty, levels = c("realised", "conc_resid"),
                     labels = c("realised fraction of the ceiling", "ancestry-profile similarity (residualised)"))]
  ggplot(dd, aes(bin, group = qty, colour = qty, fill = qty)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_ribbon(aes(ymin = null_lo, ymax = null_hi), colour = NA, alpha = 0.25) +
    geom_line(aes(y = null_mean), linetype = 3) +
    geom_line(aes(y = mean)) + geom_point(aes(y = mean), size = 1.4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
    scale_colour_manual(values = c("#e41a1c", "#377eb8"), aesthetics = c("colour", "fill"), name = NULL) +
    labs(x = xlab, y = "mean per pair",
         caption = sprintf("dotted line + band: population-label permutation baseline (%d permutations; mean, 95%%)", bl$n_perm)) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom") +
    guides(colour = guide_legend(ncol = 1))
}
mk <- function(d, xlab, stem) {
  p <- patchwork::wrap_plots(plot_ld(d, xlab) + ggtitle("a  differentiation vs LD"),
                             plot_shared(d, xlab) + ggtitle("b  shared population partitions"), nrow = 1)
  ggsave(file.path(OUT_FIG, paste0(stem, ".png")), p, width = 11, height = 5, dpi = 200)
  ggsave(file.path(OUT_FIG, paste0(stem, ".pdf")), p, width = 11, height = 5)
}
mk(res_cm, "genetic distance between units", "03_decay_cM")
mk(res_bp, "physical distance between units", "03_decay_bp")
cat("[03] done\n")
