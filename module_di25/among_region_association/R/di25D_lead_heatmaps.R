## =========================================================================
## DI25 among-region -- genotype view of the surviving pairwise-DMI LEADS.
## =========================================================================
## The 2 leads (di25D_network_build.R, after FDR+paralogy+merge+leverage) are SINGLE-SNP
## units, so instead of member-marker blocks we show the two-locus JOINT genotype:
##   (left)  individuals x the two loci (aqu-oriented dosage), ordered by locus A, with a
##           POPULATION colour bar -- reveals co-/anti-variation AND whether it is driven by
##           among-population structure (the R_st axis) rather than a within-pop incompatibility.
##   (right) the 3x3 genotype contingency (obs count / obs:exp under independence) -- the DMI
##           read: repulsion depletes the "wrong" ancestry combinations, coupling the diagonal.
## Genotype: 2 = F. aquilonia hom (blue), 1 = het (grey), 0 = F. polyctena hom (red).
##
## Reads di25D_network.rds (leads) + di25D_units.rds (aqu-oriented dosages, pops). Writes
## Figures/di25D_lead_heatmaps.png. NB descriptive only -- the DMI call is Module E's.

suppressPackageStartupMessages(library(data.table))
MAX_LEADS <- 8   # cap the figure to the top leads by |within-pop r| (q<0.1 yields ~100)
N <- readRDS("module_di25/among_region_association/data/di25D_network.rds")
u <- readRDS("module_di25/among_region_association/data/di25D_units.rds")
D <- u$dosage; pops <- as.character(u$pops); chr_of <- u$chr_of
leads <- N$meta_edges[order(-abs(within_pop_r))][seq_len(min(MAX_LEADS, .N))]
gcol <- c("#B2182B", "#F0F0F0", "#2166AC")            # 0 pol / 1 het / 2 aqu
## population order shared by both leads: by overall ancestry (mean dosage over ALL units),
## most-aquilonia populations at the top, so genotype-vs-population alignment is comparable.
pop_anc  <- tapply(rowMeans(D, na.rm = TRUE), pops, mean)
pop_order <- names(sort(pop_anc, decreasing = TRUE))
popcol <- setNames(grDevices::hcl.colors(length(pop_order), "Spectral"), pop_order)

## rows grouped BY POPULATION (not by genotype), populations in the shared ancestry order;
## within a population, individuals ordered by locus A so within-pop variation is visible.
panel_indiv <- function(A, B, ttl) {
  gA <- round(D[, A]); gB <- round(D[, B]); n <- length(gA)
  ord <- order(match(pops, pop_order), -gA, -gB)
  gA <- gA[ord]; gB <- gB[ord]; pv <- factor(pops[ord], levels = pop_order)
  cA <- ifelse(is.na(gA), "grey85", gcol[gA + 1]); cB <- ifelse(is.na(gB), "grey85", gcol[gB + 1])
  plot(NA, xlim = c(-1.7, 3.2), ylim = c(0, n), axes = FALSE, xlab = "", ylab = "", main = ttl, cex.main = 0.85)
  yy <- n:1
  rect(-0.05, yy - 1, 0.6, yy, col = popcol[as.character(pv)], border = NA)   # population bar
  rect(1, yy - 1, 2, yy, col = cA, border = NA)                               # locus A
  rect(2, yy - 1, 3, yy, col = cB, border = NA)                               # locus B
  ## population separators + labels at each block's vertical midpoint
  bl <- rle(as.integer(pv)); ends <- cumsum(bl$lengths); starts <- ends - bl$lengths + 1
  mids <- (starts + ends) / 2
  abline(h = n - ends[-length(ends)] + 0.5, col = "white", lwd = 1.2)
  text(-0.15, n - mids + 1, pop_order[bl$values], cex = 0.5, adj = 1, xpd = NA)
  text(c(0.28, 1.5, 2.5), -n * 0.03, c("pop", A, B), cex = 0.7, xpd = NA)
}

panel_joint <- function(A, B, coup, wpr, wprloo, rst) {
  gA <- round(D[, A]); gB <- round(D[, B]); ok <- is.finite(gA) & is.finite(gB)
  tab <- table(factor(gA[ok], 2:0), factor(gB[ok], 2:0))        # rows/cols: aqu, het, pol
  N0 <- sum(tab); exp <- outer(rowSums(tab), colSums(tab)) / N0
  oe <- tab / exp
  plot(NA, xlim = c(0, 3), ylim = c(0, 3.6), axes = FALSE, xlab = "", ylab = "",
       main = sprintf("%s  within-pop r=%.2f (LOO %.2f)  R_st=%.2f", coup, wpr, wprloo, rst), cex.main = 0.8)
  lab <- c("aqu", "het", "pol")
  for (i in 1:3) for (j in 1:3) {
    cnt <- tab[i, j]; ratio <- oe[i, j]
    shade <- if (cnt == 0) "#000000" else grDevices::adjustcolor(if (ratio >= 1) "#2166AC" else "#B2182B",
                                                                 min(0.85, abs(log2(pmax(ratio, 1e-3))) / 3 + 0.12))
    rect(j - 1, 3 - i, j, 4 - i, col = shade, border = "white")
    text(j - 0.5, 3.5 - i, if (cnt == 0) "0 (absent)" else sprintf("%d\n%.1fx", cnt, ratio),
         cex = 0.62, col = if (cnt == 0) "white" else "black")
  }
  text(-0.15, 3.5 - (1:3), lab, cex = 0.7, xpd = NA); text(0.5 + 0:2, 3.15, lab, cex = 0.7, xpd = NA)
  text(-0.15, 3.35, sprintf("%s\\%s", A, B), cex = 0.6, xpd = NA, adj = 1)
}

nL <- nrow(leads)
png("module_di25/among_region_association/Figures/di25D_lead_heatmaps.png",
    width = 1700, height = 560 * nL + 120, res = 200)
layout(matrix(seq_len(2 * nL), ncol = 2, byrow = TRUE), widths = c(1.5, 1))
par(mar = c(2.5, 1, 2.5, 1))
for (k in seq_len(nL)) {
  A <- leads$ma[k]; B <- leads$mb[k]
  panel_indiv(A, B, sprintf("%s (%s) x %s (%s)  [%s]", A, chr_of[A], B, chr_of[B], leads$coupling[k]))
  panel_joint(A, B, leads$coupling[k], leads$within_pop_r[k], leads$wpr_loo[k], leads$rst[k])
}
## shared legend
par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0)); plot.new()
legend("bottom", horiz = TRUE, bty = "n", cex = 0.7,
       legend = c("aquilonia (2)", "het (1)", "polyctena (0)", "missing"),
       fill = c(gcol[3], gcol[2], gcol[1], "grey85"))
dev.off()
cat("[done] Figures/di25D_lead_heatmaps.png\n")
