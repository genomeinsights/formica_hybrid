## Figures for the Module D work-in-progress summary (EMMAX arm, paralogy filter).
## Reads data/moduleD_emmax.rds + data/moduleD_ohta.rds. Writes two PNGs to Figures/.
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
dir.create("moduleD_among_region_association/Figures", showWarnings = FALSE)

E <- readRDS("moduleD_among_region_association/data/moduleD_emmax.rds")
chr_all <- sort(unique(E$results[[1]]$Chr))
chr_order <- paste0("Chr", sort(as.integer(sub("Chr", "", chr_all))))

## ---- Fig E: two representative structure-corrected Manhattans -----------
man <- function(gid, main) {
  dt <- copy(E$results[[gid]])[Chr %in% chr_order]
  dt[, Chr := factor(Chr, levels = chr_order)]; setorder(dt, Chr, cM); dt[, x := .I]
  dt[, lp := -log10(pmax(pval, .Machine$double.xmin))]
  col <- setNames(rep(c("grey62", "grey38"), length.out = length(chr_order)), chr_order)
  sig <- 0.05 / nrow(dt)
  plot(dt$x, dt$lp, pch = 20, cex = 0.3, col = col[as.character(dt$Chr)], xaxt = "n",
       xlab = "", ylab = expression(-log[10]~p), main = main, cex.main = 0.95)
  abline(h = -log10(sig), lty = 2, col = "grey70")
  points(dt[partner == gid, x], dt[partner == gid, lp], col = "red", pch = 20, cex = 1.1)  # self
  ph <- E$hits[focal == gid]
  hh <- dt[unlinked == TRUE & pval < sig]
  if (nrow(hh)) {
    hh[, par := partner %in% ph[paralog == TRUE, partner]]
    points(hh[par == FALSE, x], hh[par == FALSE, lp], col = "#E69F00", pch = 1, cex = 1.4, lwd = 1.6)
    points(hh[par == TRUE,  x], hh[par == TRUE,  lp], col = "grey45", pch = 4, cex = 1.2, lwd = 1.3)
  }
  cm <- dt[, .(mid = mean(x)), by = Chr]; axis(1, at = cm$mid, labels = sub("Chr", "", cm$Chr), cex.axis = 0.5, las = 2)
}
png("moduleD_among_region_association/Figures/moduleD_emmax_writeup.png", width = 2100, height = 1350, res = 200)
par(mfrow = c(2, 1), mar = c(2.6, 4, 2.2, 1), mgp = c(2, 0.5, 0))
man("F4250",  "a  Focal F4250 (Chr2): self-peak (red) + unlinked coupling peaks (clean = orange, paralog = grey x)")
man("F14324", "b  Focal F14324, a consensus-bidirectional locus (Chr4): self-peak only — null background, structure-corrected")
dev.off()

## ---- Fig P: the paralogy discriminant (within-pop |r| and excess Ho) -----
eh <- E$hits[, .(within_pop_r, maxHo = pmax(het1, het2), paralog, arm = "EMMAX hits")]
O  <- readRDS("moduleD_among_region_association/data/moduleD_ohta.rds")$pairs
oh <- O[, .(within_pop_r, maxHo = pmax(het1, het2), paralog, arm = "Ohta survivors")]
set.seed(1); ohs <- oh[sample(.N, min(.N, 5000))]
d <- rbind(eh, ohs)
th <- theme_classic(base_size = 9) +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(size = 9.5),
        legend.key.size = unit(3, "mm"))

pa <- ggplot(d, aes(within_pop_r, colour = arm)) +
  geom_vline(xintercept = 0.9, linetype = 2, colour = "grey60") +
  geom_density(na.rm = TRUE) + scale_colour_manual(values = c("EMMAX hits" = "#0072B2", "Ohta survivors" = "#009E73")) +
  annotate("text", x = 0.92, y = 0, label = "flag > 0.9", hjust = 0, vjust = -0.5, size = 2.6, colour = "grey40") +
  labs(x = "within-population |r|", y = "density", title = "a  Duplication discriminant") + th

pb <- ggplot(d, aes(within_pop_r, maxHo, colour = paralog)) +
  geom_vline(xintercept = 0.9, linetype = 2, colour = "grey60") +
  geom_point(size = 0.5, alpha = 0.4, na.rm = TRUE) +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#D55E00"), labels = c("candidate", "duplicate")) +
  labs(x = "within-population |r|", y = "max heterozygosity of the pair",
       title = "b  |r| ~ 1 co-occurs with excess Ho") + th

ggsave("moduleD_among_region_association/Figures/moduleD_paralogy_writeup.png", pa + pb, width = 210, height = 90, units = "mm", dpi = 200)
cat("wrote Figures/moduleD_emmax_writeup.png, Figures/moduleD_paralogy_writeup.png\n")
