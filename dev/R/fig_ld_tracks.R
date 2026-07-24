## =========================================================
## Supplementary figure: ld_w co-localises with the recombination environment
## =========================================================
## Shows, for two example chromosomes (Chr26 and Chr10, each with a pronounced
## low-recombination region), that per-SNP local LD support (ld_w) spikes exactly
## where the LD-decay rate (a) is LOW and where the map-based recombination rate
## is LOW -- i.e. in low-recombination regions. Note the directions: a is a decay
## RATE, so low a means slow decay (high LD, low recombination); a and the
## recombination rate are positive correlates of each other, while ld_w is their
## inverse. A 2x2 panel: rows = chromosome, columns = the two comparisons
## (ld_w vs a, ld_w vs recombination rate).
##
## Pulls only the lean saved inputs written by R/LD_decay_from_DIEM.R, plus the
## recombination map -- it does NOT load the genotype matrices or the 966 MB
## decay object:
##   data/ld_tracks_a_windows.rds   -- windowed decay rate a (Chr, mid, a, regime)
##   data/ld_tracks_ldw_persnp.rds  -- per-SNP ld_w (Chr, Pos, ld_w_095)
##   data/Frufa_DTOL_PR.ref_genome.recmap -- map-based recombination (cM/Mb)
## Run from repo root.

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

CHRS <- c("Chr26", "Chr10")

aw  <- readRDS("data/ld_tracks_a_windows.rds")[Chr %in% CHRS & regime == "structured" & is.finite(a)]
ldw <- readRDS("data/ld_tracks_ldw_persnp.rds")[Chr %in% CHRS & is.finite(ld_w_095)]
rec <- fread("data/Frufa_DTOL_PR.ref_genome.recmap")
rec[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
rec <- rec[Chr %in% CHRS, .(Chr, pos, cMMb = `cM/Mb`)]

## each series is on its own scale; min-max normalise so co-localisation is
## visible on a common 0-1 axis within each panel.
sc01 <- function(v) { v <- as.numeric(v); (v - min(v, na.rm = TRUE)) / diff(range(v, na.rm = TRUE)) }

th <- theme_classic(base_size = 8) +
  theme(plot.tag = element_text(face = "bold", size = 10),
        plot.title = element_text(size = 8.5, hjust = 0),
        axis.title = element_text(size = 8), axis.text = element_text(size = 6.5),
        legend.position = "none")

col_a <- "#D55E00"; col_r <- "#0072B2"; col_ldw <- "grey55"

panel <- function(ch, which = c("a", "rec")) {
  which <- match.arg(which)
  lw <- ldw[Chr == ch]; lw[, `:=`(x = Pos / 1e6, y = sc01(ld_w_095))]
  if (which == "a") {
    tr <- aw[Chr == ch][order(mid)]; tr[, `:=`(x = mid / 1e6, y = sc01(a))]
    col <- col_a; ttl <- sprintf("%s: ld_w and decay rate a", ch)
  } else {
    tr <- rec[Chr == ch][order(pos)]; tr[, `:=`(x = pos / 1e6, y = sc01(cMMb))]
    col <- col_r; ttl <- sprintf("%s: ld_w and recombination rate", ch)
  }
  ggplot() +
    geom_point(data = lw, aes(x, y), colour = col_ldw, size = 0.22, alpha = 0.35) +
    geom_line(data = tr, aes(x, y), colour = col, linewidth = 0.45) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(x = "position (Mbp)", y = "scaled (0-1)", title = ttl) + th
}

## no in-figure caption: this figure carries a LaTeX \caption in the methods doc.
fig <- (panel("Chr26", "a") | panel("Chr26", "rec")) /
       (panel("Chr10", "a") | panel("Chr10", "rec")) +
  plot_annotation(tag_levels = "a")

dir.create("Figures", showWarnings = FALSE)
ggsave("Figures/ld_tracks_chr26_chr10.pdf", fig, width = 180, height = 110, units = "mm")
ggsave("Figures/ld_tracks_chr26_chr10.png", fig, width = 180, height = 110, units = "mm", dpi = 300)
cat("Saved: Figures/ld_tracks_chr26_chr10.pdf/.png\n")
