## =========================================================
## MODULE 0 -- consensus fidelity histogram (manuscript Fig S6)
## =========================================================
## Distribution of consensus round-trip fidelity
##   score_eMLG = cor(round(x), x)^2
## across the ld_w flagging thresholds of the >=5-marker eMLG construction sweep,
## split by SINGLETON clusters (n_loci == 1, no merge possible -- trivially perfect
## score) vs genuine MULTI-MARKER merges. The dashed 0.80 line is the fidelity floor
## ld_prune_and_eMLG() enforces when merging. Separating the two clarifies that the
## mass near the upper limit at permissive thresholds reflects isolated markers, not
## merges stopping early.
##
## Ported from the (deprecated) flat R/fig_supplementary.R figS06 block; this is a
## Stage-0 complexity-reduction figure, so it now lives with Module 0. The full
## five-threshold version is R/plot_eMLG_score_histograms.R (kept in /doc); here the
## ld_w = 0.2 panel is dropped so the four panels fit the page.
##
## Reads : module0_ld_pruning/data/results_min_loci5.rds  (threshold sweep of the
##         >=5-marker eMLG construction; produced by module0_ld_pruning_DIEM.R)
## Writes: module0_ld_pruning/Figures/module0_fidelity_hist.{pdf,png}
##         manuscript/figures/figS06_fidelity_hist.{pdf,png}   (the manuscript figure)
## Run from the formica_hybrid repo root.

suppressMessages({ library(data.table); library(ggplot2) })

## manuscript palette (matches R/fig_supplementary.R)
Accent <- "#315B7D"; Pol <- "#D08A45"; Grey <- "#C7CDD2"
DROP_TH <- 0.2                                   # panel dropped to fit the page (full set in /doc)

res <- readRDS("module0_ld_pruning/data/results_min_loci5.rds")
fl <- rbindlist(lapply(res, function(r) {
  g <- r$groups[startsWith(r$groups$group_id, "F")]
  g[, .(score, th, n_loci)]
}))
fl <- fl[th != DROP_TH]
fl[, th_label := sprintf("ld_w threshold = %s   (n = %s, mean = %.3f)",
                         th, format(.N, big.mark = ","), mean(score)), by = th]
fl[, th_label := factor(th_label, levels = unique(th_label[order(th)]))]
fl[, cluster_type := fifelse(n_loci == 1, "singleton (no merge possible)", "merged (n_loci > 1)")]
hb <- fl[, { br <- seq(0.80, 1, by = 0.01); h <- hist(score, breaks = br, plot = FALSE)
  .(bin_mid = h$mids, count = h$counts) }, by = .(th_label, cluster_type)]
hb <- fl[, .(N = .N), by = th_label][hb, on = "th_label"][, pct := 100 * count / N]

p <- ggplot(hb, aes(bin_mid, pct, fill = cluster_type)) +
  geom_col(colour = "white", linewidth = 0.1, width = 0.01) +
  geom_vline(xintercept = 0.80, linetype = 2, colour = Pol, linewidth = 0.5) +
  scale_fill_manual(values = c("singleton (no merge possible)" = Grey, "merged (n_loci > 1)" = Accent),
                    name = NULL) +
  facet_wrap(~ th_label, ncol = 1) + coord_cartesian(xlim = c(0.80, 1)) +
  labs(x = expression(score[eMLG] == cor(round(x), x)^2),
       y = "% of that run's flagged clusters") +
  theme_classic(base_size = 9) +
  theme(legend.position = "top", legend.text = element_text(size = 7),
        strip.text = element_text(size = 7.5, face = "bold"), axis.title = element_text(size = 8))

sv <- function(dir, f) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(dir, paste0(f, ".pdf")), p, width = 120, height = 150, units = "mm")
  ggsave(file.path(dir, paste0(f, ".png")), p, width = 120, height = 150, units = "mm", dpi = 300)
}
sv("module0_ld_pruning/Figures", "module0_fidelity_hist")
sv("manuscript/figures", "figS06_fidelity_hist")
cat("Saved: module0_ld_pruning/Figures/module0_fidelity_hist.{pdf,png},",
    "manuscript/figures/figS06_fidelity_hist.{pdf,png}\n")
