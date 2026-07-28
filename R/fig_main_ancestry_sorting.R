## =========================================================
## MAIN MANUSCRIPT FIGURE: ancestry sorting (4 panels A-D)
## =========================================================
## Assembles the central narrative figure from existing analysis outputs, in the
## manuscript palette (matches manuscript/main_text_*.tex):
##   A  directional sorting: sort-class proportions among ORIENTED gated loci
##   B  direction reverses with DI: P(aquilonia | unidirectional) across DI deciles
##   C  LD reveals pseudoreplication: fraction sorted vs recombination, SNP vs unit
##   D  climate effects not robust: DI-enrichment coefficient, naive vs block-bootstrap CI
##
## Reads : data/moduleA_snp.rds, data/moduleA_di_asymmetry.rds,
##         data/moduleB_architecture.rds, data/sensitivity_block_bootstrap.rds
## Writes: manuscript/figures/main_ancestry_sorting.{pdf,png}
## Run from the repo root. Point estimates untouched; this only re-plots them.

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

## ---- manuscript palette (from main_text_*.tex) --------------------------
Accent <- "#315B7D"; MidBlue <- "#8EB4CE"; Muted <- "#5B6570"
Aqu <- "#3A7D77"; Pol <- "#D08A45"; Grey <- "#C7CDD2"
th <- theme_classic(base_size = 9) +
  theme(plot.tag = element_text(face = "bold", size = 11, colour = Accent),
        plot.title = element_text(size = 10, face = "bold", colour = Accent),
        axis.title = element_text(size = 8), legend.position = "none",
        legend.title = element_blank(), legend.text = element_text(size = 6.5),
        legend.key.size = unit(3, "mm"),
        legend.background = element_rect(fill = scales::alpha("white", 0.7), colour = NA),
        plot.margin = margin(4, 6, 4, 4))

## ---- A: sort-class proportions among oriented gated loci ----------------
lev <- c("unsorted", "polyctena", "aquilonia", "bidirectional")
pcol <- c(unsorted = Grey, polyctena = Pol, aquilonia = Aqu, bidirectional = Muted)
a <- as.data.table(readRDS("data/moduleA_snp.rds"))[
  differentiated == TRUE & !is.na(sort_class), .N, by = sort_class]
a[, pct := 100 * N / sum(N)]
pctv <- a[, setNames(pct, as.character(sort_class))]
a[, cls := factor(as.character(sort_class), levels = lev)]
pA <- ggplot(a, aes(cls, pct, fill = cls)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.45, size = 2.5, colour = Muted) +
  scale_fill_manual(values = pcol[lev]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), limits = c(0, NA)) +
  labs(title = "A", x = NULL, y = "% of oriented parent-polymorphic loci") +
  th + theme(axis.text.x = element_text(size = 7))

## ---- B: P(aquilonia | unidirectional) across DI deciles -----------------
di <- as.data.table(readRDS("data/moduleA_di_asymmetry.rds")$di_decile)
di[, side := ifelse(frac_aqu_of_unidir >= 0.5, "aquilonia", "polyctena")]
pB <- ggplot(di, aes(DI_decile, frac_aqu_of_unidir)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = Muted) +
  geom_line(colour = Muted, linewidth = 0.5) +
  geom_point(aes(colour = side), size = 2) +
  scale_colour_manual(values = c(aquilonia = Aqu, polyctena = Pol)) +
  scale_x_continuous(breaks = 1:10) + ylim(0, 1) +
  annotate("text", x = 8.2, y = 0.90, label = "toward aquilonia", colour = Aqu, size = 2.4) +
  annotate("text", x = 2.8, y = 0.06, label = "toward polyctena", colour = Pol, size = 2.4) +
  labs(title = "B", x = "diagnostic-index decile (low → high)",
       y = "fraction fixing toward aquilonia") + th

## ---- C: fraction sorted vs recombination, SNP vs LD-reduced unit ---------
b <- readRDS("data/moduleB_architecture.rds")
cdt <- rbind(
  as.data.table(b$unit_by_recomb)[, .(med_r, frac_sorted, level = "LD-reduced unit")],
  as.data.table(b$snp_by_recomb)[,  .(med_r, frac_sorted, level = "SNP")])
cdt[, level := factor(level, levels = c("SNP", "LD-reduced unit"))]
pC <- ggplot(cdt, aes(med_r, frac_sorted, colour = level)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.7) + scale_x_log10() +
  scale_colour_manual(values = c("LD-reduced unit" = Accent, "SNP" = Grey)) +
  labs(title = "C", x = "recombination rate (cM/Mb, log)", y = "fraction sorted") +
  th + theme(legend.position = c(0.68, 0.16))

## ---- D: DI-enrichment coefficient, naive vs chromosome block-bootstrap ---
r <- as.data.table(readRDS("data/sensitivity_block_bootstrap.rds")$results)
d <- r[grepl("diagnostic", target)][, pc := ifelse(grepl("PC1", target), "PC1", "PC2")]
dd <- rbind(
  d[, .(pc, method = "model-based", est = estimate, lo = naive_lo, hi = naive_hi)],
  d[, .(pc, method = "chromosome block", est = estimate, lo = chr_lo, hi = chr_hi)])
dd[, pc := factor(pc, levels = c("PC2", "PC1"))]
dd[, method := factor(method, levels = c("model-based", "chromosome block"))]
pD <- ggplot(dd, aes(est, pc, colour = method)) +
  geom_vline(xintercept = 0, linetype = 2, colour = Muted) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18,
                 position = position_dodge(width = 0.55), linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.55), size = 1.8) +
  scale_colour_manual(values = c("model-based" = MidBlue, "chromosome block" = Accent)) +
  labs(title = "D", x = "diagnostic enrichment (pp per +10 pp climate rate)", y = NULL) +
  th + theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
             axis.text.y = element_text(face = "bold", colour = Accent))

fig <- (pA | pB) / (pC | pD) + plot_layout(heights = c(0.9, 1))
dir.create("manuscript/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("manuscript/figures/main_ancestry_sorting.pdf", fig, width = 180, height = 150, units = "mm")
ggsave("manuscript/figures/main_ancestry_sorting.png", fig, width = 180, height = 150, units = "mm", dpi = 300)
cat("Saved manuscript/figures/main_ancestry_sorting.{pdf,png}\n")
