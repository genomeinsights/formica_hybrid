## =============================================================================
## module_manuscript_rho05 -- Fig 3 (main) multipanel, rho05 version
## =============================================================================
## Adapted from module_di25/R/di25_sorting_multipanel.R. Only panel (b)
## (sorting vs recombination across tau) depends on eMLG clustering -- swapped
## to di25_recomb_tau_sweep_rho05.rds (min_r2_rho=0.5 Stage-2, built earlier
## this session). Panels a/c/d are SNP-level (di25_sorting_emp_vs_sim.rds,
## di25_sorting_snp.rds, di25_ld_decay.rds) and are IDENTICAL regardless of
## min_r2 -- reused unchanged, including the cached c/d PNGs (no need to
## recompute the 5000-permutation BDMI-overlap histogram or the circos render).
##
## Run from the repo root: Rscript module_manuscript_rho05/R/di25_sorting_multipanel_rho05.R
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(cowplot); library(magick) })
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()
source("module_di25/R/diem_circos_core.R")          # render_circos_raster()

FIG <- "module_manuscript_rho05/Figures"; DATA_DI25 <- "module_di25/data"
DATA_RHO05 <- "module_manuscript_rho05/data"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
OUTPNG <- file.path(FIG, "di25_sorting_multipanel_rho05.png"); OUTPDF <- sub("png$", "pdf", OUTPNG)
## c/d panels are SNP-level -- reuse the EXISTING cache (module_di25/data/),
## identical content regardless of min_r2, no need to duplicate or recompute.
C_PNG  <- file.path(DATA_DI25, "mp_bdmi_hist_notitle.png")
D_PNG  <- file.path(DATA_DI25, "mp_bdmi_circos_notitle.png")
stopifnot("panel c/d cache missing -- run module_di25/R/di25_sorting_multipanel.R once first" =
            file.exists(C_PNG), file.exists(D_PNG))

## ============================ panel a (ggplot, SNP-level, unchanged) ========
o   <- readRDS(file.path(DATA_DI25, "di25_sorting_emp_vs_sim.rds"))
sim <- as.data.table(o$sim); emp <- as.data.table(o$emp)
pA <- ggplot(sim, aes(factor(tau), pct_sorted)) +
  geom_jitter(width = 0.18, height = 0, colour = "#1b9e77", alpha = 0.25, size = 0.5) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = NA, colour = "#0b3d2e") +
  geom_point(data = emp, aes(factor(tau), pct_sorted), colour = "#d95f02", size = 3.2, shape = 18) +
  labs(x = expression("sorting threshold  " * tau), y = "% SNPs sorted") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.margin = margin(t = 10, r = 5, b = 3, l = 34))

## ============================ panel b (ggplot, rho05 recomb-tau-sweep) ======
rs   <- readRDS(file.path(DATA_RHO05, "di25_recomb_tau_sweep_rho05.rds"))
dec  <- as.data.table(rs$deciles); sweep <- as.data.table(rs$sweep)
ne   <- as.data.table(rs$n_units_per_decile)[order(rbin)]
de   <- dec[level == "eMLG"]; de[, tau := factor(tau)]
med_lab <- de[, .(m = med_recomb[1]), by = rbin][order(rbin)]$m
ymax <- max(de$frac_sorted); bs <- ymax * 0.98 / max(ne$N); N_MIN <- 100L
lo_bins <- ne[N < N_MIN, rbin]
pB1 <- ggplot(de, aes(rbin)) +
  { if (length(lo_bins)) annotate("rect", xmin = min(lo_bins) - 0.5, xmax = max(lo_bins) + 0.5, ymin = -Inf, ymax = Inf, fill = "grey95") } +
  geom_col(data = ne, aes(rbin, N * bs), fill = "grey86", width = 0.85, inherit.aes = FALSE) +
  geom_line(aes(y = frac_sorted, colour = tau), linewidth = 0.9) +
  geom_point(aes(y = frac_sorted, colour = tau), size = 1.8) +
  scale_colour_viridis_d(end = 0.9, name = expression(tau)) +
  scale_x_continuous(breaks = 1:10, labels = med_lab) +
  scale_y_continuous(name = "fraction sorted (eMLG units, rho05)",
                     sec.axis = sec_axis(~ . / bs, name = "n independent units")) +
  labs(x = "recombination decile (median cM/Mb)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.margin = margin(t = 14, r = 5, b = 3, l = 8))
pB2 <- ggplot(sweep, aes(factor(tau), coef_recomb)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_errorbar(aes(ymin = boot_lo, ymax = boot_hi), width = 0.18, colour = "#315B7D") +
  geom_point(size = 2.6, colour = "#315B7D") +
  labs(x = expression(tau), y = "recomb coef (logit / log10 cM/Mb)") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
pB <- plot_grid(pB1, pB2, ncol = 2, rel_widths = c(2, 1))

## ================================ compose ===================================
pC <- ggdraw() + draw_image(C_PNG, y = 0, height = 0.92)
pD <- ggdraw() + draw_image(D_PNG, y = 0, height = 0.96)
row_bc <- plot_grid(pA, pC, ncol = 2, rel_widths = c(1, 2), labels = c("b", "c"), label_size = 24, label_fontface = "bold")
right  <- plot_grid(row_bc, pB, ncol = 1, rel_heights = c(1, 1.8), labels = c("", "d"), label_size = 24, label_fontface = "bold")
full   <- plot_grid(pD, right, ncol = 2, rel_widths = c(1, 1.3), labels = c("a", ""), label_size = 24, label_fontface = "bold")
full   <- full + theme(plot.margin = margin(12, 6, 4, 6))
ggsave(OUTPNG, full, width = 17, height = 6.8, dpi = 200, bg = "white")
ggsave(OUTPDF, full, width = 17, height = 6.8, device = cairo_pdf, bg = "white")
cat("saved:", OUTPNG, "\n       ", OUTPDF, "\n")
