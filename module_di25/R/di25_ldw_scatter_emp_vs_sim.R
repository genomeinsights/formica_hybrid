## =============================================================================
## module_di25 -- ld_w_095 windowed median: empirical vs simulated-median SCATTER
## =============================================================================
## One point per genomic window: empirical windowed-median ld_w_0.95 vs the
## across-replicate MEDIAN of the simulated windowed medians (1000 DIEM reps).
## Reads the aggregate produced by di25_ldw_manhattan_envelope.R.
##   Rscript module_di25/R/di25_ldw_scatter_emp_vs_sim.R
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2) })

IN     <- "module_di25/data/di25_ldw_envelope.rds"
OUTPNG <- "module_di25/Figures/di25_ldw_scatter_emp_vs_sim.png"
OUTPDF <- sub("\\.png$", ".pdf", OUTPNG)

o  <- readRDS(IN)
pl <- as.data.table(o$pl)[is.finite(emp) & is.finite(sim_med)]

r_p <- cor(pl$emp, pl$sim_med, method = "pearson")
r_s <- cor(pl$emp, pl$sim_med, method = "spearman")
lab <- sprintf("n = %d windows\nPearson r = %.3f\nSpearman rho = %.3f",
               nrow(pl), r_p, r_s)

p <- ggplot(pl, aes(sim_med, emp)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55", linewidth = 0.7) +
  geom_point(alpha = 0.4, size = 1.7, colour = "#1b6f5f", stroke = 0) +
  geom_smooth(method = "lm", se = FALSE, colour = "#d95f02", linewidth = 1.1) +
  annotate("text", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 5, label = lab) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = expression("simulated median (1000 reps)  " * ld[w][0.95]),
       y = expression("empirical  " * ld[w][0.95])) +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank())

ggsave(OUTPDF, p, width = 6.5, height = 6.5)
ggsave(OUTPNG, p, width = 6.5, height = 6.5, dpi = 150)
cat(sprintf("n=%d  Pearson=%.3f  Spearman=%.3f\nsaved: %s\n", nrow(pl), r_p, r_s, OUTPNG))
