## =========================================================
## module_manuscript_rho05 -- Module C combined calibration figure
## PC1 + PC2 + bio_winter on the same panels
## =========================================================
## Purely a figure-combination step: PC1/PC2 (moduleC_stage1_results.rds)
## and bio_winter (moduleC_stage1_bio_winter_results.rds) were tested
## against the IDENTICAL 10,000-null-covariate distribution at the primary
## cell (min05_tau06) -- verified byte-for-byte equal below -- so all three
## observed covariates can be drawn on one set of null histograms without
## re-running or re-deriving anything.
##
## Reads : module_manuscript_rho05/data/moduleC_stage1_results.rds
##         module_manuscript_rho05/data/moduleC_stage1_bio_winter_results.rds
## Writes: module_manuscript_rho05/Figures/moduleC_stage1_combined_null_calibration.{png,pdf}
##         module_manuscript_rho05/Figures/moduleC_stage1_combined_threshold_sensitivity.{png,pdf}
##
## Run from the repo root:
##   Rscript module_manuscript_rho05/R/moduleC_stage1_combined_calibration_figure.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2) })

BASE <- "module_manuscript_rho05"
FIG  <- file.path(BASE, "Figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

res <- readRDS(file.path(BASE, "data", "moduleC_stage1_results.rds"))          # PC1, PC2
bw  <- readRDS(file.path(BASE, "data", "moduleC_stage1_bio_winter_results.rds"))   # bio_winter

stopifnot(
  "PC1/PC2 and bio_winter were tested at different primary cells" =
    identical(res$meta$primary_cell, bw$meta$primary_cell),
  "PC1/PC2 and bio_winter null distributions differ -- do not combine" =
    isTRUE(all.equal(res$null_stats, bw$null_stats))
)
null <- res$null_stats
obs  <- rbind(res$observed, bio_winter = bw$observed[colnames(res$observed)])
AXES <- rownames(obs)   # PC1, PC2, bio_winter
COLS <- c(PC1 = "#0072B2", PC2 = "#D55E00", bio_winter = "#009E73")

LAB <- c(rho_DI = "DI (Spearman rho)", rho_rec = "recombination (Spearman rho)",
         sort_gap_differentiated = "sorting, differentiated only (BF percentile gap)",
         rho_sort_magnitude = "sorting magnitude / prop_fixed (Spearman rho)")

## ========================================================================
## FIGURE 1 -- structured-null calibration, PC1 + PC2 + bio_winter overlaid
## ========================================================================
fig_stats <- names(LAB)
long <- rbindlist(lapply(fig_stats, function(s)
  data.table(label = LAB[[s]], value = null[, s])))
long[, label := factor(label, levels = LAB[fig_stats])]

mk <- rbindlist(lapply(fig_stats, function(s) rbindlist(lapply(AXES, function(ax)
  data.table(label = LAB[[s]], axis = ax,
             median = median(null[, s]),
             lo = quantile(null[, s], 0.025), hi = quantile(null[, s], 0.975),
             value = obs[ax, s])))))
mk[, label := factor(label, levels = LAB[fig_stats])]
mk[, axis := factor(axis, levels = AXES)]
mk_band <- unique(mk[, .(label, median, lo, hi)])

p1 <- ggplot(long, aes(value)) +
  geom_histogram(bins = 60, fill = "grey80", colour = NA) +
  geom_rect(data = mk_band, inherit.aes = FALSE, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            fill = "#999999", alpha = 0.18) +
  geom_vline(data = mk_band, aes(xintercept = median), linetype = 2, colour = "grey40") +
  geom_vline(data = mk, aes(xintercept = value, colour = axis), linewidth = 0.7) +
  facet_wrap(~ label, scales = "free", ncol = 2) +
  scale_colour_manual(values = COLS, name = NULL, breaks = AXES) +
  labs(x = "genome-wide statistic", y = "count (of 10,000 null covariates)",
       title = "Module C (Stage-1-direct): observed PC1 / PC2 / bio_winter vs 10,000 Omega-structured null covariates",
       subtitle = "shaded = null 95% interval; dashed = null median; PC1/PC2/bio_winter share the identical null (verified)") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8),
        strip.background = element_blank())
ggsave(file.path(FIG, "moduleC_stage1_combined_null_calibration.png"), p1,
       width = 175, height = 140, units = "mm", dpi = 300)
ggsave(file.path(FIG, "moduleC_stage1_combined_null_calibration.pdf"), p1,
       width = 175, height = 140, units = "mm")

## ========================================================================
## FIGURE 2 -- rank-threshold sensitivity, PC1 + PC2 + bio_winter overlaid
## ========================================================================
th_stats <- grep("^top", colnames(null), value = TRUE)
frac_lookup <- c(top0001 = "top 0.1%", top0005 = "top 0.5%", top0010 = "top 1%")
metric_lookup <- c(DI_med = "median DI", rec_med = "median recomb",
                   pdiff = "prop. differentiated", pdir_diff = "prop. directional | diff.")
meta_th <- data.table(stat = th_stats)
meta_th[, threshold_tag := sub("_.*", "", stat)]
meta_th[, metric_tag    := sub("^top[0-9]+_", "", stat)]
meta_th[, frac   := unname(frac_lookup[threshold_tag])]
meta_th[, metric := unname(metric_lookup[metric_tag])]
stopifnot(all(!is.na(meta_th$frac)), all(!is.na(meta_th$metric)))

lev  <- c("top 0.1%", "top 0.5%", "top 1%")
mlev <- unname(metric_lookup)
setfac <- function(d) { d[, frac := factor(frac, levels = lev)]
                        d[, metric := factor(metric, levels = mlev)]; d }

thl <- rbindlist(lapply(th_stats, function(s) data.table(stat = s, value = null[, s])))
thl <- setfac(meta_th[thl, on = "stat"])

thm <- rbindlist(lapply(th_stats, function(s) rbindlist(lapply(AXES, function(ax)
  data.table(stat = s, axis = ax, value = obs[ax, s])))))
thm <- setfac(meta_th[thm, on = "stat"])
thm[, axis := factor(axis, levels = AXES)]
thm_band <- setfac(meta_th[, .(stat, frac, metric)][
  data.table(stat = th_stats, median = apply(null[, th_stats, drop = FALSE], 2, median, na.rm = TRUE),
             lo = apply(null[, th_stats, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE),
             hi = apply(null[, th_stats, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)), on = "stat"])

p2 <- ggplot(thl, aes(value)) +
  geom_histogram(bins = 45, fill = "grey80", colour = NA) +
  geom_rect(data = thm_band, inherit.aes = FALSE, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            fill = "#999999", alpha = 0.18) +
  geom_vline(data = thm_band, aes(xintercept = median), linetype = 2, colour = "grey40") +
  geom_vline(data = thm, aes(xintercept = value, colour = axis), linewidth = 0.6) +
  facet_wrap(vars(metric, frac), scales = "free", ncol = 3,
             labeller = label_wrap_gen(multi_line = FALSE)) +
  scale_colour_manual(values = COLS, name = NULL, breaks = AXES) +
  labs(x = "value within the top-ranked fraction of Stage-1 units", y = "count (null covariates)",
       title = "Module C (Stage-1-direct) threshold sensitivity: PC1 / PC2 / bio_winter vs null",
       subtitle = "prop. directional is computed among differentiated units in the fraction") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom", strip.text = element_text(size = 7),
        strip.background = element_blank())
ggsave(file.path(FIG, "moduleC_stage1_combined_threshold_sensitivity.png"), p2,
       width = 185, height = 165, units = "mm", dpi = 300)
ggsave(file.path(FIG, "moduleC_stage1_combined_threshold_sensitivity.pdf"), p2,
       width = 185, height = 165, units = "mm")

cat("\n[combined-fig] wrote 2 combined figures (calibration + threshold sensitivity), PC1/PC2/bio_winter overlaid\n")
