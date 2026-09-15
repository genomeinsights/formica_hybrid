## =========================================================
## module_localscore_crosscheck -- two threshold-free calibration
## diagnostics, complementing the window-count check
## (local_score_null_check_S1units.R): (1) correlation of each null
## draw's Stage-1-cluster BF/C2 vector with the REAL observed vector,
## histogrammed per covariate, with the real pairwise correlations
## between observed covariates marked for reference; (2) QQ-plots of the
## observed BF/C2 distribution against the pooled null distribution.
## =========================================================
## Both diagnostics answer: "how similar to the observed result does pure
## population structure alone typically produce?" -- (1) at the level of a
## single summary number (genome-wide correlation) per null draw, (2) at
## the level of the full value distribution (does the observed depart from
## null only in the extreme tail, or broadly across the whole distribution).
##
## Uses the SAME 1,000-draw null pools (first 5 persisted batches) as
## local_score_null_check_S1units.R, at Stage-1-cluster resolution
## (18,361 units, MRK order -- NOT the position-sorted order local-score
## needs, so no re-sorting required here).
##
## Reads : module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/
##           {PC1,PC2,bio_winter}_S1units_..._betai_reg.out,
##           mito_C2_S1units_summary_contrast.out,
##           null/bf_matrices/cRegen_bf_b01..05.rds, mitoC2_bf_b01..05.rds
## Writes: module_localscore_crosscheck/Figures/null_correlation_histograms.{png,pdf}
##         module_localscore_crosscheck/Figures/null_qq_plots.{png,pdf}
##         module_localscore_crosscheck/data/null_correlation_summary.tsv
##
## Run from the repo root:
##   Rscript module_localscore_crosscheck/R/null_calibration_correlation_qq.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

UNIT_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded_S1units"
BFDIR    <- file.path(UNIT_DIR, "null", "bf_matrices")
DATA     <- "module_localscore_crosscheck/data"
FIGDIR   <- "module_localscore_crosscheck/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

NBATCH_SCAN <- 5L; BATCH <- 200L; N_NULL <- NBATCH_SCAN * BATCH

## ---- observed Stage-1-cluster vectors (MRK order, matching null matrices) --
obs <- list(
  PC1        = fread(file.path(UNIT_DIR, "PC1_S1units_withOmega_summary_betai_reg.out"))$`BF(dB)`,
  PC2        = fread(file.path(UNIT_DIR, "PC2_S1units_withOmega_summary_betai_reg.out"))$`BF(dB)`,
  bio_winter = fread(file.path(UNIT_DIR, "bio_winter_S1units_withOmega_summary_betai_reg.out"))$`BF(dB)`,
  mitoC2     = fread(file.path(UNIT_DIR, "mito_C2_S1units_summary_contrast.out"))$`log10(1/pval)`
)
NM <- length(obs$PC1)
stopifnot(all(vapply(obs, length, integer(1)) == NM))

## ---- load the two null pools (1000 draws each) -----------------------------
load_pool <- function(pattern) {
  mats <- lapply(seq_len(NBATCH_SCAN), function(b) readRDS(sprintf(pattern, b)))
  do.call(cbind, mats)   # NM x N_NULL
}
cont_pool <- load_pool(file.path(BFDIR, "cRegen_bf_b%02d.rds"))    # shared: PC1/PC2/bio_winter
mito_pool <- load_pool(file.path(BFDIR, "mitoC2_bf_b%02d.rds"))
stopifnot(nrow(cont_pool) == NM, ncol(cont_pool) == N_NULL, nrow(mito_pool) == NM, ncol(mito_pool) == N_NULL)

## ========================================================================
## 1. correlation of each null draw with the observed vector, per covariate
## ========================================================================
POOL_FOR <- list(PC1 = cont_pool, PC2 = cont_pool, bio_winter = cont_pool, mitoC2 = mito_pool)
cor_null <- lapply(names(obs), function(tag) {
  P <- POOL_FOR[[tag]]; o <- obs[[tag]]
  vapply(seq_len(ncol(P)), function(j) cor(P[, j], o), numeric(1))
})
names(cor_null) <- names(obs)

## real pairwise BF-vs-BF correlations between the three continuous
## covariates' OBSERVED vectors (a direct, single-number test of the
## "shared architecture" pattern noted qualitatively for the local-score
## windows). PC1 and PC2 are PCA axes -- their sign is an arbitrary
## artifact of the eigendecomposition, not biologically meaningful -- so
## the DEFENSIBLE summary statistic for "is this relationship stronger
## than chance" is |r| (equivalently, r^2; both give identical percentile
## rankings), not signed r. Signed r is still reported/plotted since
## bio_winter and mitoC2 individually DO have a meaningful, fixed sign,
## but every percentile/significance claim below uses |r|.
real_cor <- combn(c("PC1", "PC2", "bio_winter"), 2, function(p)
  data.table(a = p[1], b = p[2], r = cor(obs[[p[1]]], obs[[p[2]]])), simplify = FALSE)
real_cor <- rbindlist(real_cor)

## empirical percentile of |real correlation| within the corresponding
## covariate's own |null-correlation| distribution
pctl_abs <- function(x, null_dist) mean(abs(null_dist) <= abs(x))
real_cor[, pctl_abs_in_a := mapply(function(aa, rr) pctl_abs(rr, cor_null[[aa]]), a, r)]
real_cor[, pctl_abs_in_b := mapply(function(bb, rr) pctl_abs(rr, cor_null[[bb]]), b, r)]
real_cor[, `:=`(pctl_abs_in_a_label = sprintf("%.1f%%", 100 * pctl_abs_in_a),
                pctl_abs_in_b_label = sprintf("%.1f%%", 100 * pctl_abs_in_b))]
cat("\nReal pairwise BF-vs-BF correlations, with |r| percentile within each covariate's own |null-correlation| distribution:\n")
print(real_cor[, .(a, b, r, pctl_abs_in_a_label, pctl_abs_in_b_label)])

## ---- histogram figure: |r|, since PC1/PC2's sign is arbitrary ---------
## Each panel: |Pearson r| between a null draw's Stage-1-cluster BF/C2
## vector and THAT PANEL'S OBSERVED BF/C2 vector (both axes of every
## correlation here are BF/C2 values -- never raw covariate values).
long <- rbindlist(lapply(names(cor_null), function(tag) data.table(covariate = tag, abs_r = abs(cor_null[[tag]]))))
long[, covariate := factor(covariate, levels = c("PC1", "PC2", "bio_winter", "mitoC2"))]

marks <- rbindlist(list(
  data.table(covariate = "PC1", other = "PC2", abs_r = abs(real_cor[a == "PC1" & b == "PC2", r])),
  data.table(covariate = "PC1", other = "bio_winter", abs_r = abs(real_cor[a == "PC1" & b == "bio_winter", r])),
  data.table(covariate = "PC2", other = "PC1", abs_r = abs(real_cor[a == "PC1" & b == "PC2", r])),
  data.table(covariate = "PC2", other = "bio_winter", abs_r = abs(real_cor[a == "PC2" & b == "bio_winter", r])),
  data.table(covariate = "bio_winter", other = "PC1", abs_r = abs(real_cor[a == "PC1" & b == "bio_winter", r])),
  data.table(covariate = "bio_winter", other = "PC2", abs_r = abs(real_cor[a == "PC2" & b == "bio_winter", r]))
))
marks[, covariate := factor(covariate, levels = c("PC1", "PC2", "bio_winter", "mitoC2"))]

p_hist <- ggplot(long, aes(abs_r)) +
  geom_histogram(bins = 60, fill = "grey75", colour = NA) +
  geom_vline(data = marks, aes(xintercept = abs_r, colour = other), linewidth = 0.8) +
  geom_text(data = marks, aes(x = abs_r, y = Inf, label = other, colour = other),
            angle = 90, vjust = -0.4, hjust = 1.05, size = 2.8, show.legend = FALSE) +
  facet_wrap(~ covariate, scales = "free", ncol = 2) +
  scale_colour_manual(values = c(PC1 = "#0072B2", PC2 = "#D55E00", bio_winter = "#009E73"), name = "real covariate") +
  coord_cartesian(clip = "off") +
  labs(x = expression("|Pearson "*italic(r)*"| between a null draw's BF/C2 and the panel's OBSERVED BF/C2"),
       y = sprintf("count (of %d null draws)", N_NULL),
       title = "|Correlation| of each Omega-structured null draw's BF/C2 with the observed Stage-1-cluster BF/C2",
       subtitle = "|r| used throughout (PC1/PC2 sign is an arbitrary PCA artifact). Coloured lines: |real BF-vs-BF correlation| with another observed covariate") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"))
ggsave(file.path(FIGDIR, "null_correlation_histograms.png"), p_hist, width = 180, height = 140, units = "mm", dpi = 300)
ggsave(file.path(FIGDIR, "null_correlation_histograms.pdf"), p_hist, width = 180, height = 140, units = "mm")

fwrite(rbindlist(lapply(names(cor_null), function(tag) data.table(
  covariate = tag, mean_r = mean(cor_null[[tag]]), sd_r = sd(cor_null[[tag]]),
  min_r = min(cor_null[[tag]]), max_r = max(cor_null[[tag]])))),
  file.path(DATA, "null_correlation_summary.tsv"), sep = "\t")

## ========================================================================
## 2. QQ-plots: observed vs. pooled-null value distribution, per covariate
## ========================================================================
qq_dt <- rbindlist(lapply(names(obs), function(tag) {
  o <- sort(obs[[tag]])
  pool_vals <- as.vector(POOL_FOR[[tag]])
  null_q <- quantile(pool_vals, probs = ppoints(length(o)), names = FALSE, type = 7)
  data.table(covariate = tag, observed = o, null = sort(null_q))
}))
qq_dt[, covariate := factor(covariate, levels = c("PC1", "PC2", "bio_winter", "mitoC2"))]

p_qq <- ggplot(qq_dt, aes(null, observed)) +
  geom_point(size = 0.3, alpha = 0.15, colour = "#0072B2") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "black") +
  facet_wrap(~ covariate, scales = "free", ncol = 2) +
  labs(x = sprintf("null quantile (pooled across %d draws)", N_NULL), y = "observed quantile",
       title = "QQ-plot: observed Stage-1-cluster BF/C2 vs. pooled Omega-structured null",
       subtitle = "Dashed line = 1:1. Points following the line = observed matches null; departure (esp. upper tail) = beyond null expectation") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
ggsave(file.path(FIGDIR, "null_qq_plots.png"), p_qq, width = 180, height = 140, units = "mm", dpi = 300)
ggsave(file.path(FIGDIR, "null_qq_plots.pdf"), p_qq, width = 180, height = 140, units = "mm")

cat("\n[null-calibration-extra] wrote null_correlation_histograms.{png,pdf}, null_qq_plots.{png,pdf}, null_correlation_summary.tsv\n")
