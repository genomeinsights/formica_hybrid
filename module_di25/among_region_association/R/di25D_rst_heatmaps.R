## =========================================================================
## DI25 among-region -- genotype heatmaps for the top AMONG-POPULATION (R_st) pairs.
## =========================================================================
## The top |R_st| unlinked trans pairs are the AMONG-population / parallel-co-sorting
## candidates: their per-population mean dosages track each other strongly across the 20
## replicates. They are the equal-fitness-DMI / systematic-among-pop axis (D'2st), but
## structure-CONFOUNDED -- these pairs have within_pop_r ~ 0 and EMMAX Rsq ~ 0 (the
## structure correction strips them). We cannot rule out shared-founding co-differentiation
## (drift), so these are CANDIDATES ONLY; shown descriptively to see what they look like.
## Contrast with di25D_lead_heatmaps.R (the WITHIN-population, EMMAX-surviving leads).
## Rows grouped by population -- the among-pop signal is the aligned population blocks.
## Reads di25D_emmax.rds + di25D_units.rds. Writes Figures/di25D_rst_heatmaps.png.

suppressPackageStartupMessages(library(data.table))
N_EACH <- 3   # top-N positive (coupling) + top-N negative (repulsion) different-chromosome R_st pairs
e <- readRDS("module_di25/among_region_association/data/di25D_emmax.rds")
u <- readRDS("module_di25/among_region_association/data/di25D_units.rds")
D <- u$dosage; pops <- as.character(u$pops); chr_of <- u$chr_of; m <- u$metrics
gcol <- c("#B2182B", "#F0F0F0", "#2166AC")            # 0 pol / 1 het / 2 aqu
pop_anc  <- tapply(rowMeans(D, na.rm = TRUE), pops, mean)
pop_order <- names(sort(pop_anc, decreasing = TRUE))
popcol <- setNames(grDevices::hcl.colors(length(pop_order), "Spectral"), pop_order)

## pick the top different-chromosome pairs by R_st (balanced coupling / repulsion)
p <- e$pairs[unlinked == TRUE & chr_of[focal] != chr_of[partner]]
pr <- rbind(head(p[order(-R_st)], N_EACH), head(p[order(R_st)], N_EACH))
pr[, `:=`(ma = focal, mb = partner)]

panel_indiv <- function(A, B, ttl) {
  gA <- round(D[, A]); gB <- round(D[, B]); n <- length(gA)
  ord <- order(match(pops, pop_order), -gA, -gB)
  gA <- gA[ord]; gB <- gB[ord]; pv <- factor(pops[ord], levels = pop_order)
  cA <- ifelse(is.na(gA), "grey85", gcol[gA + 1]); cB <- ifelse(is.na(gB), "grey85", gcol[gB + 1])
  plot(NA, xlim = c(-1.7, 3.2), ylim = c(0, n), axes = FALSE, xlab = "", ylab = "", main = ttl, cex.main = 0.85)
  yy <- n:1
  rect(-0.05, yy - 1, 0.6, yy, col = popcol[as.character(pv)], border = NA)
  rect(1, yy - 1, 2, yy, col = cA, border = NA); rect(2, yy - 1, 3, yy, col = cB, border = NA)
  bl <- rle(as.integer(pv)); ends <- cumsum(bl$lengths); mids <- ends - bl$lengths / 2
  abline(h = n - ends[-length(ends)] + 0.5, col = "white", lwd = 1.2)
  text(-0.15, n - mids + 1, pop_order[bl$values], cex = 0.5, adj = 1, xpd = NA)
  text(c(0.28, 1.5, 2.5), -n * 0.03, c("pop", A, B), cex = 0.7, xpd = NA)
}
panel_joint <- function(A, B, rst, wpr, mafA, mafB) {
  gA <- round(D[, A]); gB <- round(D[, B]); ok <- is.finite(gA) & is.finite(gB)
  tab <- table(factor(gA[ok], 2:0), factor(gB[ok], 2:0)); N0 <- sum(tab)
  oe <- tab / (outer(rowSums(tab), colSums(tab)) / N0)
  plot(NA, xlim = c(0, 3), ylim = c(0, 3.6), axes = FALSE, xlab = "", ylab = "",
       main = sprintf("%s  R_st=%.2f  within-pop r=%.2f  MAF %.2f/%.2f",
                      ifelse(rst > 0, "coupling", "repulsion"), rst, wpr, mafA, mafB), cex.main = 0.72)
  lab <- c("aqu", "het", "pol")
  for (i in 1:3) for (j in 1:3) {
    cnt <- tab[i, j]; ratio <- oe[i, j]
    shade <- if (cnt == 0) "#000000" else grDevices::adjustcolor(if (ratio >= 1) "#2166AC" else "#B2182B",
                                                                 min(0.85, abs(log2(pmax(ratio, 1e-3))) / 3 + 0.12))
    rect(j - 1, 3 - i, j, 4 - i, col = shade, border = "white")
    text(j - 0.5, 3.5 - i, if (cnt == 0) "0" else sprintf("%d\n%.1fx", cnt, ratio),
         cex = 0.62, col = if (cnt == 0) "white" else "black")
  }
  text(-0.15, 3.5 - (1:3), lab, cex = 0.7, xpd = NA); text(0.5 + 0:2, 3.15, lab, cex = 0.7, xpd = NA)
  text(-0.15, 3.35, sprintf("%s\\%s", A, B), cex = 0.6, xpd = NA, adj = 1)
}

nL <- nrow(pr)
png("module_di25/among_region_association/Figures/di25D_rst_heatmaps.png",
    width = 1700, height = 560 * nL + 120, res = 200)
layout(matrix(seq_len(2 * nL), ncol = 2, byrow = TRUE), widths = c(1.5, 1))
par(mar = c(2.5, 1, 2.5, 1))
for (k in seq_len(nL)) {
  A <- pr$ma[k]; B <- pr$mb[k]
  panel_indiv(A, B, sprintf("%s (%s) x %s (%s)  [top R_st]", A, chr_of[A], B, chr_of[B]))
  panel_joint(A, B, pr$R_st[k], pr$within_pop_r[k], m$maf[match(A, m$group_id)], m$maf[match(B, m$group_id)])
}
par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0)); plot.new()
legend("bottom", horiz = TRUE, bty = "n", cex = 0.75,
       legend = c("aquilonia (2)", "het (1)", "polyctena (0)", "missing"),
       fill = c(gcol[3], gcol[2], gcol[1], "grey85"))
dev.off()
cat(sprintf("[done] Figures/di25D_rst_heatmaps.png | %d top-R_st trans pairs (|R_st| %.2f-%.2f; within_pop_r ~ %.2f)\n",
            nL, min(abs(pr$R_st)), max(abs(pr$R_st)), median(pr$within_pop_r)))
