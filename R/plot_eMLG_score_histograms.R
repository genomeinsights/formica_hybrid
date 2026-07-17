library(data.table)
library(ggplot2)

# ------------------------------------------------------------
# Score histogram diagnostic across ld_w_threshold values.
#
# Illustrates the argument in dev/methods_notes.md: because
# dynamic_cut_eMLG() only accepts a merge while score_eMLG stays above
# score_threshold (0.80, dashed line here), a score distribution piled up
# near 1.0 indicates under-merging relative to the available quality
# budget (merging stopped for a reason other than quality -- ran out of
# eligible candidates within a distance-restricted run), while a
# distribution that spreads out toward the 0.80 floor indicates the cut is
# more fully exploiting that budget. Lowering ld_w_threshold flags more
# Stage-1 clusters into Stage 2's runs, giving the dynamic cut more
# material to merge across.
#
# Input: ./data/results_min_loci5.rds -- a list of ld_prune_and_eMLG()
# outputs (list(eMLG, groups, pruned)) at ld_w_threshold = 0.025, 0.05, 0.1,
# 0.15, 0.2 (min_n_loci_eMLG = 5 throughout, hence the filename). Each run's
# groups$th column records which threshold produced it; the number of runs
# is read from the data, not hardcoded.
# ------------------------------------------------------------

results <- readRDS("./data/results_min_loci5.rds")

flagged <- rbindlist(lapply(results, function(r) {
  g <- r$groups[startsWith(r$groups$group_id, "F")]
  g[, .(score, th, n_loci)]
}))
flagged[, th_label := sprintf(
  "ld_w_threshold = %s   n = %s   mean = %.3f   %% > 0.99 = %.1f%%",
  th, format(.N, big.mark = ","), mean(score), 100 * mean(score > 0.99)
), by = th]
flagged[, th_label := factor(th_label, levels = unique(th_label[order(th)]))]

## A cluster at score ~= 1 can mean two very different things: a singleton
## (n_loci == 1) that was flagged but had no viable merge partner at all
## (trivially perfect score_eMLG, not evidence of anything), or a genuine
## multi-locus merge that reached near-perfect fidelity. Splitting the
## histogram by this distinguishes "ceiling mass is mostly untouched
## singletons" from "ceiling mass includes real, large, successful merges" --
## these have very different implications for whether ld_w_threshold is
## leaving mergeable material uncorrected.
flagged[, cluster_type := fifelse(n_loci == 1, "singleton (n_loci = 1, no merge possible)",
                                   "merged (n_loci > 1)")]

score_threshold <- 0.80
binwidth <- 0.01

## proportions (not raw counts) within each threshold's flagged-cluster set,
## so panel shapes are directly comparable despite very different n across
## thresholds (tens of thousands of flagged clusters at low ld_w_threshold
## down to hundreds at high ld_w_threshold) -- the point is where the mass
## sits relative to the 0.80 floor, not the absolute count
n_th <- flagged[, .N, by = th_label]
hist_dt <- flagged[, {
  breaks <- seq(score_threshold, 1, by = binwidth)
  h <- hist(score, breaks = breaks, plot = FALSE)
  .(bin_mid = h$mids, count = h$counts)
}, by = .(th_label, cluster_type)]
hist_dt <- n_th[hist_dt, on = "th_label"]
hist_dt[, pct := 100 * count / N]

p_score_hist <- ggplot(hist_dt, aes(bin_mid, pct, fill = cluster_type)) +
  geom_col(color = "white", linewidth = 0.1, width = binwidth) +
  geom_vline(xintercept = score_threshold, linetype = "dashed",
             color = "#E41A1C", linewidth = 0.6) +
  scale_fill_manual(values = c("singleton (n_loci = 1, no merge possible)" = "grey70",
                                "merged (n_loci > 1)" = "#377EB8"), name = NULL) +
  facet_wrap(~th_label, ncol = 1) +
  coord_cartesian(xlim = c(score_threshold, 1)) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "top") +
  labs(
    x = expression(score[eMLG] == cor(round(x), x)^2),
    y = "% of that run's flagged clusters",
    title = "Stage 2 (dynamic cut) score distribution by ld_w_threshold",
    subtitle = paste0(
      "Dashed line = score_threshold (", score_threshold, "), the floor every accepted merge must clear.\n",
      "Grey = singleton clusters (flagged but never merged); blue = genuine multi-locus merges."
    )
  )

dir.create("./Figures/", showWarnings = FALSE, recursive = TRUE)
n_panels <- uniqueN(flagged$th_label)
ggsave("./Figures/eMLG_score_histograms_by_threshold.png", p_score_hist,
       width = 7.5, height = 2.4 * n_panels + 1.5, dpi = 150, limitsize = FALSE)
message("Saved: ./Figures/eMLG_score_histograms_by_threshold.png")

p_score_hist
