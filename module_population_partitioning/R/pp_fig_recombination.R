## =========================================================================
## module_population_partitioning -- 12: recombination-rate figure.
## Three panels: (a) local similarity vs GENETIC (cM) distance, the direct
## analogue of fig2's physical-distance curve; (b) at fixed physical
## distance (100-500kb), local similarity by recombination-rate tertile,
## with block-bootstrap CI on the low-vs-high contrast; (c) per-unit local
## (<=100kb) concordance vs the unit's own local recombination-rate decile.
##
## Run from the formica_hybrid repo root, after pp_recombination.R:
##   Rscript module_population_partitioning/R/pp_fig_recombination.R
## Writes: module_population_partitioning/Figures/fig5_recombination.png
## =========================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
DATADIR <- "module_population_partitioning/data"
FIGDIR  <- "module_population_partitioning/Figures"
rc <- readRDS(file.path(DATADIR, "pp_recombination.rds"))
u <- rc$u; setDT(u)
stopifnot("unit table must have exactly the 20,807 DI25 rho05 units (min_r2_rho=0.5) -- check for a legacy/stale input" = nrow(u) == 20807L)

theme_ms <- theme_bw(base_size = 12) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank())
COL_SIGNED <- "#1b9e77"; COL_LOW <- "#2166ac"; COL_MID <- "grey50"; COL_HIGH <- "#b2182b"

## (a) genetic-distance decay
cmd <- copy(rc$cm_decay)
fig_a <- ggplot(cmd, aes(cmbin, mean_r, group = 1)) +
  geom_line(colour = COL_SIGNED) + geom_point(size = 2.2, colour = COL_SIGNED) +
  labs(x = "genetic distance between unit pair", y = "mean signed r",
      title = "(a) similarity vs genetic distance") +
  theme_ms + theme(axis.text.x = element_text(angle = 40, hjust = 1))

## (b) fixed physical distance (100-500kb), by recombination-rate tertile
st <- rc$strat[dbin == "100-500kb"]
ci <- data.table(recomb_tertile = c("low", "high"),
                 lo = quantile(rc$boot_contrast[, "100-500kb"], 0.025, na.rm = TRUE),
                 hi = quantile(rc$boot_contrast[, "100-500kb"], 0.975, na.rm = TRUE))
## (b) audit fix item 4: subtitle now reports BOTH the original single-bin
## descriptive contrast AND the distance-adjusted (20kb strata within the
## same 100-500kb scope) contrast + its own chromosome-block bootstrap CI,
## so the more precise result sits directly alongside the descriptive one
## rather than replacing it.
adj <- rc$adjusted_contrast
fig_b <- ggplot(st, aes(recomb_tertile, mean_r, fill = recomb_tertile)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c(low = COL_LOW, mid = COL_MID, high = COL_HIGH), guide = "none") +
  labs(x = "local recombination-rate tertile", y = "mean signed r",
      title = "(b) same physical distance (100-500kb): low-recomb\nunits stay more concordant",
      subtitle = sprintf(paste0("single-bin contrast = %.3f, CI [%.3f, %.3f]\n",
                                "distance-adjusted (20kb strata) contrast = %.3f, CI [%.3f, %.3f] -- supports an effect"),
                         rc$obs_contrast["100-500kb"],
                         quantile(rc$boot_contrast[, "100-500kb"], 0.025, na.rm = TRUE),
                         quantile(rc$boot_contrast[, "100-500kb"], 0.975, na.rm = TRUE),
                         adj$obs, adj$ci[1], adj$ci[2])) +
  theme_ms + theme(plot.subtitle = element_text(size = 7.8))

## (c) per-unit local concordance vs local recombination-rate decile
rd <- copy(rc$recomb_dec_tab)
fig_c <- ggplot(rd, aes(mean_recomb, mean_near_r)) +
  geom_point(size = 2.3, colour = COL_SIGNED) + geom_line(colour = COL_SIGNED) +
  labs(x = "local recombination-rate decile mean (cM/Mb)", y = "mean local signed r (<=100kb)",
      title = "(c) higher local recombination -> lower local concordance") +
  theme_ms

fig5 <- fig_a + fig_b + fig_c + plot_layout(ncol = 3)
ggsave(file.path(FIGDIR, "fig5_recombination.png"), fig5, width = 15, height = 5, dpi = 200)
cat("[fig5] saved fig5_recombination.png\n")
