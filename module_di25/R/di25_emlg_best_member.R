## =========================================================
## module_di25 -- best single-SNP proxy for the eMLG consensus
## =========================================================
## Follow-up to di25_emlg_vs_representative.R. The representative SNP is chosen by
## the clustering (centrality / ld_w), NOT to maximise agreement with the block
## consensus. Here, for every eMLG block we correlate EACH member SNP with the
## consensus and take the maximum, so we can ask: how well can the BEST single SNP
## in a block stand in for the consensus, and how far short does the representative
## fall of that best member?
##
## Correlation is on |r| (orientation is free -- you would flip a SNP's coding if
## needed), so "best single-SNP proxy" = member with the highest |r| to consensus.
##
## Run from the repo root:  Rscript module_di25/R/di25_emlg_best_member.R
## =========================================================

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

OUTDIR <- "module_di25/data"
FIGDIR <- "module_di25/Figures"
CM     <- "cM5"

inp <- readRDS(file.path(OUTDIR, "di25_inputs.rds"))
GTs <- inp$GTs_hyb                                    # hybrids x markers (012, NA = miss)

res    <- readRDS(file.path(OUTDIR, sprintf("di25_clustering_%s.rds", CM)))
emlg   <- res$eMLG                                    # hybrids x blocks (consensus)
groups <- as.data.table(res$groups)
G      <- GTs[rownames(emlg), , drop = FALSE]         # align hybrids

gi <- groups[match(colnames(emlg), group_id)]
## members is a list-column of character vectors (one per group)
members_list <- gi$members

## ---- per-block: correlate every member with the consensus -----------------
rows <- lapply(seq_len(ncol(emlg)), function(j) {
  cons  <- emlg[, j]
  mrk   <- members_list[[j]]
  mrk   <- mrk[mrk %in% colnames(G)]
  Gsub  <- G[, mrk, drop = FALSE]
  rr    <- suppressWarnings(cor(Gsub, cons, use = "pairwise.complete.obs"))[, 1]  # one r per member
  ar    <- abs(rr)
  rep_m <- gi$representative[j]
  best_i <- which.max(ar)
  data.table(
    group_id  = colnames(emlg)[j],
    n_loci    = gi$n_loci[j],
    rep_marker  = rep_m,
    rep_abs_r   = if (rep_m %in% mrk) ar[match(rep_m, mrk)] else NA_real_,
    best_marker = if (length(best_i)) mrk[best_i] else NA_character_,
    best_abs_r  = if (length(best_i)) ar[best_i] else NA_real_,
    rep_is_best = length(best_i) > 0 && mrk[best_i] == rep_m
  )
})
dt <- rbindlist(rows)
dt[, gain := best_abs_r - rep_abs_r]     # how much the best member beats the representative
saveRDS(dt, file.path(OUTDIR, "di25_emlg_best_member.rds"))

## ---- summaries ------------------------------------------------------------
main <- dt[is.finite(best_abs_r)]
cat(sprintf("\n===== [%s] best single-SNP proxy for the eMLG consensus (%d blocks) =====\n",
            CM, nrow(main)))

cat("\nBest member |r| with consensus (max over members):\n")
print(round(quantile(main$best_abs_r, c(0, .05, .10, .25, .5, .75, .90, .95, 1)), 3))

cat("\nRepresentative |r| with consensus (for comparison):\n")
print(round(quantile(main$rep_abs_r, c(0, .05, .10, .25, .5, .75, .90, .95, 1), na.rm = TRUE), 3))

cat(sprintf("\nRepresentative IS the best member in %.1f%% of blocks.\n",
            100 * mean(main$rep_is_best)))
cat(sprintf("Median gain (best - representative |r|): %.3f;  90th pct gain: %.3f\n",
            median(main$gain, na.rm = TRUE), quantile(main$gain, .90, na.rm = TRUE)))
cat(sprintf("Blocks where best member reaches |r| >= 0.95: %.1f%%  (representative: %.1f%%)\n",
            100 * mean(main$best_abs_r >= 0.95),
            100 * mean(main$rep_abs_r >= 0.95, na.rm = TRUE)))

cat("\nBy cluster size (n_loci): median best |r|, median representative |r|, median gain:\n")
main[, size_bin := cut(n_loci, breaks = c(2, 3, 5, 10, 20, 50, Inf),
                       labels = c("3", "4-5", "6-10", "11-20", "21-50", ">50"), right = TRUE)]
print(main[, .(n = .N,
               best_r = round(median(best_abs_r), 3),
               rep_r  = round(median(rep_abs_r, na.rm = TRUE), 3),
               gain   = round(median(gain, na.rm = TRUE), 3),
               rep_is_best = round(mean(rep_is_best), 3)),
           by = size_bin][order(size_bin)])

## ---- figure ---------------------------------------------------------------
long <- rbind(
  data.table(which = "representative SNP", abs_r = main$rep_abs_r),
  data.table(which = "best member SNP",    abs_r = main$best_abs_r)
)
p1 <- ggplot(long, aes(abs_r, fill = which)) +
  geom_histogram(binwidth = 0.02, position = "identity", alpha = 0.55, colour = NA) +
  scale_fill_manual(values = c("best member SNP" = "#21918C", "representative SNP" = "#FDE725")) +
  labs(x = "|r| with eMLG consensus", y = "eMLG blocks", fill = NULL,
       title = sprintf("Best single-SNP proxy vs representative (%s, %d blocks)", CM, nrow(main)),
       subtitle = sprintf("median best |r| = %.3f vs representative %.3f; representative is best in %.0f%% of blocks",
                          median(main$best_abs_r), median(main$rep_abs_r, na.rm = TRUE),
                          100 * mean(main$rep_is_best))) +
  theme_bw(base_size = 11) + theme(legend.position = "top")

p2 <- ggplot(main, aes(n_loci, gain)) +
  geom_point(alpha = 0.25, size = 0.7, colour = "#440154") +
  geom_smooth(method = "loess", se = TRUE, colour = "firebrick", linewidth = 0.7) +
  scale_x_log10() +
  labs(x = "cluster size (n markers, log scale)",
       y = "gain in |r| (best member - representative)",
       title = "How much a better-chosen single SNP would improve on the representative") +
  theme_bw(base_size = 11)

dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(FIGDIR, "di25_emlg_best_member.png"), p1 / p2, width = 8, height = 8, dpi = 200)
cat(sprintf("\nWrote %s and %s\n",
            file.path(OUTDIR, "di25_emlg_best_member.rds"),
            file.path(FIGDIR, "di25_emlg_best_member.png")))
