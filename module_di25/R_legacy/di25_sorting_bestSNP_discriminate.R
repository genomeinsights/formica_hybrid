## =========================================================
## module_di25 -- are the best-SNP-only sorted units real signal or thin calls?
## =========================================================
## The filled best-SNP sorts more units than the consensus. This asks whether the
## EXTRA units (sorted by best-SNP, not consensus) are:
##   * real recovered signal  -> many populations near-fixed (high n_obs, high
##     n_fixed), comparable to the units both methods agree on; or
##   * thin near-fixation calls -> clear tau only because few populations had data
##     (low n_obs) or only a handful are fixed (low n_fixed).
## prop_fixed = n_fixed / n_obs over 20 hybrid pops; sorted when prop_fixed >= tau.
##
## Run from the repo root:  Rscript module_di25/R/di25_sorting_bestSNP_discriminate.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()

OUTDIR <- "module_di25/data"; FIGDIR <- "module_di25/Figures"
TAU <- 0.6; SORT_RULE <- "binom"; ALPHA <- 0.05

cmp <- readRDS(file.path(OUTDIR, "di25_sorting_bestSNP_vs_consensus.rds"))
ps_cons <- as.data.table(cmp$ps_cons); ps_best <- as.data.table(cmp$ps_best)
bm  <- as.data.table(readRDS(file.path(OUTDIR, "di25_emlg_best_member.rds")))
score_lookup <- { g <- as.data.table(readRDS(sprintf("module_di25/data/di25_clustering_cM5.rds"))$groups)
                  setNames(g$score, g$group_id) }

## per-unit near-fixation stats + sorted flag at TAU, for each representation
enrich <- function(ps) {
  ps <- copy(ps)
  ps[, n_fixed := n_aqu + n_pol]
  ps[, prop_fixed := ifelse(n_obs > 0, n_fixed / n_obs, NA_real_)]
  cls <- classify_sort(ps$n_aqu, ps$n_pol, ps$n_obs, sort_th = TAU,
                       sort_rule = SORT_RULE, alpha = ALPHA)
  ps[, sorted := differentiated == TRUE & n_obs > 0 &
       cls %in% c("aquilonia","polyctena","unresolved","ambiguous")]
  ps[]
}
C <- enrich(ps_cons); B <- enrich(ps_best)

## join the two representations per block
J <- merge(
  C[, .(marker, c_nobs = n_obs, c_nfix = n_fixed, c_prop = prop_fixed, c_sorted = sorted)],
  B[, .(marker, b_nobs = n_obs, b_nfix = n_fixed, b_prop = prop_fixed, b_sorted = sorted)],
  by = "marker")
J[, set := fifelse(c_sorted & b_sorted, "both",
            fifelse(!c_sorted & b_sorted, "bestSNP_only",
            fifelse(c_sorted & !b_sorted, "consensus_only", "neither")))]
J[, `:=`(n_loci = bm$n_loci[match(marker, bm$group_id)],
         best_abs_r = bm$best_abs_r[match(marker, bm$group_id)],
         score = score_lookup[marker])]

cat(sprintf("\n===== discriminating check @ tau = %.2f (n = %d blocks) =====\n", TAU, nrow(J)))
cat("set sizes:\n"); print(J[, .N, by = set][order(-N)])

## ---- core comparison: best-SNP stats, bestSNP_only vs both -----------------
cat("\n-- best-SNP near-fixation stats by set (median [IQR]) --\n")
sm <- J[set %in% c("both","bestSNP_only"),
        .(n = .N,
          nobs_med    = median(b_nobs),
          nfix_med    = median(b_nfix),
          nfix_q25    = quantile(b_nfix, .25),
          prop_med    = round(median(b_prop), 3),
          thin_nfix_le3 = round(mean(b_nfix <= 3), 3),
          nobs_lt15   = round(mean(b_nobs < 15), 3),
          n_loci_med  = as.double(median(n_loci)),
          bestr_med   = round(median(best_abs_r), 3),
          score_med   = round(median(score), 3)),
        by = set]
print(sm)

## ---- why did the consensus miss these? ------------------------------------
bo <- J[set == "bestSNP_only"]
cat(sprintf("\n-- for the %d best-SNP-only units, what did the consensus look like? --\n", nrow(bo)))
cat(sprintf("  consensus prop_fixed: median %.3f (tau = %.2f)\n", median(bo$c_prop), TAU))
cat(sprintf("  near-miss (consensus prop in [tau-0.1, tau)) : %.1f%%\n",
            100 * mean(bo$c_prop >= TAU - 0.1 & bo$c_prop < TAU)))
cat(sprintf("  genuine gap (consensus prop < tau-0.1)       : %.1f%%\n",
            100 * mean(bo$c_prop < TAU - 0.1)))
cat(sprintf("  filling added populations (b_nobs > c_nobs)  : %.1f%% (median +%.0f pops)\n",
            100 * mean(bo$b_nobs > bo$c_nobs), median(bo$b_nobs - bo$c_nobs)))
cat(sprintf("  best-SNP gained fixed pops (b_nfix > c_nfix)  : %.1f%% (median +%.0f)\n",
            100 * mean(bo$b_nfix > bo$c_nfix), median(bo$b_nfix - bo$c_nfix)))

## ---- verdict heuristic ----------------------------------------------------
both_nfix <- median(J[set=="both"]$b_nfix)
cat(sprintf("\n-- verdict inputs --\n  both-sorted median n_fixed = %d;  bestSNP_only median n_fixed = %d\n",
            both_nfix, median(bo$b_nfix)))
cat(sprintf("  bestSNP_only that are THIN (n_fixed <= 3 OR n_obs < 12): %.1f%%\n",
            100 * mean(bo$b_nfix <= 3 | bo$b_nobs < 12)))
cat(sprintf("  bestSNP_only that are ROBUST (n_obs >= 15 & n_fixed >= %d): %.1f%%\n",
            round(TAU*15), 100 * mean(bo$b_nobs >= 15 & bo$b_nfix >= round(TAU*15))))

## ---- figure ---------------------------------------------------------------
plotd <- J[set %in% c("both","bestSNP_only")]
p1 <- ggplot(plotd, aes(b_nobs, fill = set)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.55) +
  scale_fill_manual(values = c(both = "#440154", bestSNP_only = "#21918C")) +
  labs(x = "n populations with data (n_obs)", y = "units", fill = NULL,
       title = "Data completeness") + theme_bw(base_size = 11)
p2 <- ggplot(plotd, aes(b_nfix, fill = set)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.55) +
  geom_vline(xintercept = 3.5, linetype = 3) +
  scale_fill_manual(values = c(both = "#440154", bestSNP_only = "#21918C")) +
  labs(x = "n populations near-fixed (n_fixed)", y = "units", fill = NULL,
       title = "Breadth of near-fixation (dotted = thin cutoff 3)") + theme_bw(base_size = 11)
ggsave(file.path(FIGDIR, "di25_sorting_bestSNP_discriminate.png"),
       (p1 / p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
       width = 8, height = 7, dpi = 200)
saveRDS(J, file.path(OUTDIR, "di25_sorting_bestSNP_discriminate.rds"))
cat(sprintf("\nWrote %s and the figure.\n", file.path(OUTDIR, "di25_sorting_bestSNP_discriminate.rds")))
