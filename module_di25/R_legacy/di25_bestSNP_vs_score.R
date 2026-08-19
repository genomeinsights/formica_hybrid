## =========================================================
## module_di25 -- does the best single-SNP-consensus correlation track the eMLG
## fidelity score?
## =========================================================
## Joins the per-block best-member correlation (di25_emlg_best_member.rds) to the
## eMLG fidelity score (res$groups$score, the score_eMLG gate, threshold 0.80) and
## asks whether blocks the pipeline flags as high-fidelity are also the ones where
## a single SNP best reproduces the consensus.
##
## Run from the repo root:  Rscript module_di25/R/di25_bestSNP_vs_score.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2) })

OUTDIR <- "module_di25/data"
FIGDIR <- "module_di25/Figures"
CM     <- "cM5"

dt     <- as.data.table(readRDS(file.path(OUTDIR, "di25_emlg_best_member.rds")))
groups <- as.data.table(readRDS(file.path(OUTDIR, sprintf("di25_clustering_%s.rds", CM)))$groups)

score_lookup <- setNames(groups$score, groups$group_id)
dt[, score := score_lookup[group_id]]
d <- dt[is.finite(best_abs_r) & is.finite(score)]

cat(sprintf("\n===== [%s] best single-SNP |r| vs eMLG fidelity score (%d blocks) =====\n", CM, nrow(d)))
cat(sprintf("score range: [%.3f, %.3f];  best_abs_r range: [%.3f, %.3f]\n",
            min(d$score), max(d$score), min(d$best_abs_r), max(d$best_abs_r)))

cat("\nCorrelation of consensus-agreement metrics with the fidelity score:\n")
cat(sprintf("  best member |r|   vs score :  Pearson %+.3f   Spearman %+.3f\n",
            cor(d$best_abs_r, d$score), cor(d$best_abs_r, d$score, method = "spearman")))
cat(sprintf("  representative|r| vs score :  Pearson %+.3f   Spearman %+.3f\n",
            cor(d$rep_abs_r, d$score, use = "complete.obs"),
            cor(d$rep_abs_r, d$score, use = "complete.obs", method = "spearman")))
cat(sprintf("  gain (best-rep)   vs score :  Pearson %+.3f   Spearman %+.3f\n",
            cor(d$gain, d$score, use = "complete.obs"),
            cor(d$gain, d$score, use = "complete.obs", method = "spearman")))

## partial: does score add anything beyond cluster size?
cat(sprintf("\n  best|r| vs score Pearson, controlling for log10(n_loci): %+.3f\n",
            {
              rs <- resid(lm(best_abs_r ~ log10(n_loci), d))
              rc <- resid(lm(score      ~ log10(n_loci), d))
              cor(rs, rc)
            }))

cat("\nMedian best|r| by fidelity-score bin:\n")
d[, score_bin := cut(score, c(-Inf, .85, .90, .95, .99, Inf),
                     labels = c("<0.85", "0.85-0.90", "0.90-0.95", "0.95-0.99", ">=0.99"))]
print(d[, .(n = .N, median_best_r = round(median(best_abs_r), 3),
            median_rep_r = round(median(rep_abs_r, na.rm = TRUE), 3),
            median_n_loci = as.double(median(n_loci))), by = score_bin][order(score_bin)])

## ---- figure ---------------------------------------------------------------
p <- ggplot(d, aes(score, best_abs_r)) +
  geom_point(aes(colour = log10(n_loci)), alpha = 0.35, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE, colour = "firebrick", linewidth = 0.7) +
  scale_colour_viridis_c(name = expression(log[10]~n~markers)) +
  labs(x = "eMLG fidelity score", y = "best member |r| with consensus",
       title = sprintf("Best single-SNP proxy vs eMLG fidelity score (%s, %d blocks)", CM, nrow(d)),
       subtitle = sprintf("Pearson %+.3f, Spearman %+.3f",
                          cor(d$best_abs_r, d$score),
                          cor(d$best_abs_r, d$score, method = "spearman"))) +
  theme_bw(base_size = 11)

dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(FIGDIR, "di25_bestSNP_vs_score.png"), p, width = 8, height = 5.5, dpi = 200)
cat(sprintf("\nWrote %s\n", file.path(FIGDIR, "di25_bestSNP_vs_score.png")))
