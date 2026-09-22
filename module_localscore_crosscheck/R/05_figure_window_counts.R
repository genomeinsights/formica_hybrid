## =========================================================
## module_localscore_crosscheck -- window-count comparison figure
## (structured vs unstructured, paired by replicate, per mode)
## =========================================================
suppressMessages({ library(data.table); library(ggplot2) })

counts <- readRDS("module_localscore_crosscheck/full_snp_null10/localscore_window_counts.rds")
counts[, mode := factor(mode, levels = c("continuous", "mitoC2"),
                         labels = c("continuous (BF)", "mitoC2 (contrast)"))]
counts[, nulltype := factor(nulltype, levels = c("structured", "unstructured"))]

p <- ggplot(counts, aes(nulltype, n_windows)) +
  geom_line(aes(group = draw), colour = "grey70", linewidth = 0.3) +
  geom_point(aes(colour = nulltype), size = 2.2) +
  facet_wrap(~mode, scales = "free_y") +
  scale_colour_manual(values = c(structured = "#21918C", unstructured = "#B03A2E"), guide = "none") +
  labs(x = NULL, y = "significant local-score windows per replicate") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank())

dir.create("module_localscore_crosscheck/Figures", showWarnings = FALSE, recursive = TRUE)
ggsave("module_localscore_crosscheck/Figures/window_counts_structured_vs_unstructured.pdf", p, width = 7, height = 4)
ggsave("module_localscore_crosscheck/Figures/window_counts_structured_vs_unstructured.png", p, width = 7, height = 4, dpi = 300)
cat("wrote window_counts_structured_vs_unstructured.{pdf,png}\n")
