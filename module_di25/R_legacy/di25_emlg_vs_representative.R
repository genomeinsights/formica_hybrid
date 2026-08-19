## =========================================================
## module_di25 -- how well does the eMLG consensus track the representative SNP?
## =========================================================
## Some information in these genotypes is SNP-specific and can differ even
## between strongly correlated SNPs in the same LD cluster. For such signals the
## single representative SNP may be preferable to the block eMLG consensus. This
## script quantifies, per eMLG block, the Pearson correlation across the 165
## hybrids between:
##   * the eMLG CONSENSUS genotype (res$eMLG column), and
##   * its REPRESENTATIVE SNP genotype (res$groups$representative, pulled from
##     the raw hybrid 012 matrix GTs_hyb).
## Only blocks with a stored consensus (has_eMLG = TRUE, i.e. > 2 markers) can
## differ from their representative; 1-2-marker groups are represented by the SNP
## itself and are excluded here by construction.
##
## Run from the repo root:  Rscript module_di25/R/di25_emlg_vs_representative.R
## =========================================================

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

OUTDIR  <- "module_di25/data"
FIGDIR  <- "module_di25/Figures"
CM_SWEEP <- c("cM05", "cM1", "cM2", "cM5", "cM10")   # sweep for robustness
CM_MAIN  <- "cM5"                                     # primary clustering

inp <- readRDS(file.path(OUTDIR, "di25_inputs.rds"))
GTs <- inp$GTs_hyb                                    # hybrids x markers (012, NA = miss)

## ---- per-block consensus vs representative correlation, for one cM cap -----
block_corr <- function(cm_stamp) {
  res    <- readRDS(file.path(OUTDIR, sprintf("di25_clustering_%s.rds", cm_stamp)))
  emlg   <- res$eMLG                                 # hybrids x blocks (consensus)
  groups <- as.data.table(res$groups)

  ## align hybrids (both are the same 165, match by name to be safe)
  stopifnot(all(rownames(emlg) %in% rownames(GTs)))
  G <- GTs[rownames(emlg), , drop = FALSE]

  gi      <- groups[match(colnames(emlg), group_id)]
  rep_mrk <- gi$representative
  stopifnot(all(rep_mrk %in% colnames(G)))

  r <- vapply(seq_len(ncol(emlg)), function(j) {
    cons    <- emlg[, j]
    rep_snp <- G[, rep_mrk[j]]
    ok <- is.finite(cons) & is.finite(rep_snp)
    if (sum(ok) < 3L) return(NA_real_)
    if (sd(cons[ok]) == 0 || sd(rep_snp[ok]) == 0) return(NA_real_)  # monomorphic
    cor(cons[ok], rep_snp[ok])
  }, numeric(1))

  data.table(cM = cm_stamp, group_id = colnames(emlg),
             n_loci = gi$n_loci, score = gi$score, r = r, abs_r = abs(r))
}

## ---- run the sweep --------------------------------------------------------
all_dt <- rbindlist(lapply(CM_SWEEP, block_corr))
saveRDS(all_dt, file.path(OUTDIR, "di25_emlg_vs_representative.rds"))

## ---- summaries ------------------------------------------------------------
summ <- all_dt[, .(
  n_blocks       = .N,
  n_monomorphic  = sum(is.na(r)),
  median_r       = median(r, na.rm = TRUE),
  q25_r          = quantile(r, .25, na.rm = TRUE),
  q05_r          = quantile(r, .05, na.rm = TRUE),
  frac_r_ge_0.9  = mean(r >= 0.9, na.rm = TRUE),
  frac_r_lt_0.7  = mean(r <  0.7, na.rm = TRUE),
  frac_negative  = mean(r <  0,   na.rm = TRUE)
), by = cM][match(CM_SWEEP, cM)]

cat("\n===== eMLG consensus vs representative SNP: correlation across 165 hybrids =====\n")
print(summ)

main <- all_dt[cM == CM_MAIN & !is.na(r)]
cat(sprintf("\n[%s] %d consensus blocks (>2 markers). Pearson r(consensus, representative):\n", CM_MAIN, nrow(main)))
print(round(quantile(main$r, c(0, .05, .10, .25, .5, .75, .90, .95, 1)), 3))

## correlation vs cluster size (the key question: do big blocks diverge?)
main[, size_bin := cut(n_loci, breaks = c(2, 3, 5, 10, 20, 50, Inf),
                       labels = c("3", "4-5", "6-10", "11-20", "21-50", ">50"),
                       right = TRUE)]
cat("\nBy cluster size (n_loci) -- median r and fraction with r < 0.7:\n")
print(main[, .(n_blocks = .N, median_r = round(median(r), 3),
               frac_r_lt_0.7 = round(mean(r < 0.7), 3)), by = size_bin][order(size_bin)])

## ---- figure ---------------------------------------------------------------
p1 <- ggplot(main, aes(r)) +
  geom_histogram(binwidth = 0.02, fill = "#21918C", colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = median(main$r), linetype = 2, colour = "firebrick") +
  labs(x = "Pearson r (eMLG consensus vs representative SNP)", y = "eMLG blocks",
       title = sprintf("Consensus–representative agreement (%s, %d blocks)", CM_MAIN, nrow(main)),
       subtitle = sprintf("median r = %.3f; %.1f%% of blocks have r < 0.7",
                          median(main$r), 100 * mean(main$r < 0.7))) +
  theme_bw(base_size = 11)

p2 <- ggplot(main, aes(n_loci, r)) +
  geom_point(alpha = 0.25, size = 0.7, colour = "#440154") +
  geom_smooth(method = "loess", se = TRUE, colour = "firebrick", linewidth = 0.7) +
  scale_x_log10() +
  labs(x = "cluster size (n markers, log scale)", y = "r (consensus vs representative)",
       title = "Larger blocks: does the consensus diverge from the representative SNP?") +
  theme_bw(base_size = 11)

dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(FIGDIR, "di25_emlg_vs_representative.png"),
       p1 / p2, width = 8, height = 8, dpi = 200)
cat(sprintf("\nWrote %s and %s\n",
            file.path(OUTDIR, "di25_emlg_vs_representative.rds"),
            file.path(FIGDIR, "di25_emlg_vs_representative.png")))
