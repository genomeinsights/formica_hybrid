## =========================================================================
## DI25 among-region -- ALL member-SNP genotypes for the top leads (raw, not consensus).
## =========================================================================
## The lead heatmaps (di25D_lead_heatmaps.R) show round(eMLG CONSENSUS). This shows every
## MEMBER SNP of each lead's two clusters, raw (di25_inputs GTs), oriented to F. aquilonia
## per marker via the parent rows, individuals grouped BY POPULATION. Lets you see the true
## MAF / het structure the consensus hides (these leads are near-fixed one direction + hets).
## Reads di25D_network.rds + di25D_units.rds + di25_inputs.rds. Writes di25D_lead_snp_heatmaps.png.

suppressPackageStartupMessages(library(data.table))
MAX_LEADS <- 6
N <- readRDS("module_di25/among_region_association/data/di25D_network.rds")
u <- readRDS("module_di25/among_region_association/data/di25D_units.rds")
inp <- readRDS("module_di25/data/di25_inputs.rds")
g <- u$groups; setkey(g, group_id)
pops <- as.character(u$pops)
GTh <- inp$GTs_hyb[rownames(u$dosage), , drop = FALSE]     # align hybrids to the unit matrix (raw 012)
GTp <- inp$GTs_par
faqu <- grep("^Faqu", rownames(GTp)); fpol <- grep("^Fpol", rownames(GTp))
gcol <- c("#B2182B", "#F0F0F0", "#2166AC")                 # 0 pol / 1 het / 2 aqu

## per-marker orientation to aqu=2 (raw diagnostic coding is 0=aqu-hom, but polarise via parents)
orient <- function(m) { x <- GTh[, m]
  if (mean(GTp[faqu, m], na.rm = TRUE) <= mean(GTp[fpol, m], na.rm = TRUE)) 2 - x else x }

pop_anc  <- tapply(rowMeans(u$dosage, na.rm = TRUE), pops, mean)
pop_order <- names(sort(pop_anc, decreasing = TRUE))
popcol <- setNames(grDevices::hcl.colors(length(pop_order), "Spectral"), pop_order)
leads <- N$meta_edges[order(-abs(within_pop_r))][seq_len(min(MAX_LEADS, .N))]

panel_snp <- function(a, b, ttl) {
  mA <- intersect(g[a, members][[1]], colnames(GTh)); mB <- intersect(g[b, members][[1]], colnames(GTh))
  nA <- length(mA); nB <- length(mB); nS <- nA + nB
  X <- vapply(c(mA, mB), orient, numeric(nrow(GTh)))        # individuals x SNPs, aqu-oriented
  ord <- order(match(pops, pop_order))
  X <- X[ord, , drop = FALSE]; pv <- factor(pops[ord], levels = pop_order); n <- nrow(X); yy <- n:1
  plot(NA, xlim = c(-1.7, nS + 2), ylim = c(0, n), axes = FALSE, xlab = "", ylab = "", main = ttl, cex.main = 0.8)
  rect(-0.05, yy - 1, 0.6, yy, col = popcol[as.character(pv)], border = NA)   # population bar
  for (j in seq_len(nS)) {
    x0 <- if (j <= nA) j else j + 0.6                       # gap between the two clusters' blocks
    cv <- X[, j]; rect(x0, yy - 1, x0 + 1, yy, col = ifelse(is.na(cv), "grey85", gcol[cv + 1]), border = NA)
  }
  bl <- rle(as.integer(pv)); ends <- cumsum(bl$lengths); mids <- ends - bl$lengths / 2
  text(-0.12, n - mids, pop_order[bl$values], cex = 0.38, adj = 1, xpd = NA)
  abline(h = n - ends[-length(ends)], col = "white", lwd = 0.8)
  text(1 + (nA - 1) / 2, -n * 0.03, sprintf("%s (%d SNP)", a, nA), cex = 0.6, xpd = NA)
  text(nA + 1.6 + (nB - 1) / 2, -n * 0.03, sprintf("%s (%d SNP)", b, nB), cex = 0.6, xpd = NA)
}

nL <- nrow(leads)
png("module_di25/among_region_association/Figures/di25D_lead_snp_heatmaps.png",
    width = 1900, height = 470 * ceiling(nL / 2) + 120, res = 190)
layout(matrix(seq_len(2 * ceiling(nL / 2)), ncol = 2, byrow = TRUE))
par(mar = c(2.2, 1, 2.4, 1))
for (k in seq_len(nL)) {
  ld <- leads[k]; mafA <- mean(u$dosage[, ld$ma], na.rm = TRUE) / 2; mafB <- mean(u$dosage[, ld$mb], na.rm = TRUE) / 2
  panel_snp(ld$ma, ld$mb, sprintf("%s (%s) x %s (%s)  [%s]  within-pop r=%.2f  MAF~%.2f/%.2f",
    ld$ma, g[ld$ma, Chr], ld$mb, g[ld$mb, Chr], ld$coupling, ld$within_pop_r,
    min(mafA, 1 - mafA), min(mafB, 1 - mafB)))
}
par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0)); plot.new()
legend("bottom", horiz = TRUE, bty = "n", cex = 0.75,
       legend = c("aquilonia (2)", "het (1)", "polyctena (0)", "missing"),
       fill = c(gcol[3], gcol[2], gcol[1], "grey85"))
dev.off()
cat("[done] Figures/di25D_lead_snp_heatmaps.png\n")
