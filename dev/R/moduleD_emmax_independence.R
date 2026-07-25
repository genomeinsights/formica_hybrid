## LD-cluster independence check: EMMAX Manhattans WITHOUT the 10 cM unlinked limit.
## If the complexity reduction produced independent clusters, a focal locus should give
## a SINGLE local peak (its self-peak) with flat neighbours; several significant peaks in
## the same region (as individual SNPs would give) means the local LD block was split into
## still-correlated clusters (residual redundancy from the conservative 0.5 cM Stage-2 merge).
## Reads data/moduleD_emmax.rds (no re-run). Writes a per-focal zoomed multipage PDF + a summary.
suppressMessages(library(data.table))
dir.create("Figures", showWarnings = FALSE)
E <- readRDS("data/moduleD_emmax.rds"); focal <- E$focal
sig <- 0.05 / E$params$scope_n
fcm <- setNames(focal$cM, focal$group_id); fchr <- setNames(focal$Chr, focal$group_id)
## the Bonferroni p-threshold maps to a single r^2 cutoff (n = 164 fixed, df2 = n-2), so we
## plot the structure-corrected r^2 and draw the equivalent line; highlighting stays p-based.
n_ind <- 164          # hybrid individuals (df1 = 1, df2 = n - 2)
df2 <- n_ind - 2
r2crit <- { Fc <- qf(1 - sig, 1, df2); Fc / (Fc + df2) }

## zoomed view of the focal CHROMOSOME (cM axis), y = r^2 (variance explained). The focal
## locus itself is NOT plotted (its self-r^2~1 stretches the axis); its position is a dashed
## red line. Other same-chr clusters that are genome-wide significant are marked blue (within
## 10 cM = excluded from DMI) or orange (>10 cM = unlinked); grey band = the +/-10 cM zone.
zoom <- function(gid, main = NULL) {
  ch <- fchr[gid]; f0 <- fcm[gid]
  dt <- E$results[[gid]][Chr == ch & is.finite(cM) & partner != gid]    # EXCLUDE self
  if (is.null(main)) main <- sprintf("%s [%s] %s", gid, focal[group_id == gid, set], ch)
  plot(dt$cM, dt$Rsq, pch = 20, cex = 0.5, col = "grey65", xlab = paste(ch, "position (cM)"),
       ylab = expression(r^2), main = main, cex.main = 0.9,
       ylim = c(0, max(0.05, max(dt$Rsq, na.rm = TRUE))))
  if (is.finite(f0)) rect(f0 - 10, -1, f0 + 10, 2, col = "#00000010", border = NA)  # +/-10 cM zone
  abline(h = r2crit, lty = 2, col = "grey70")                            # p-threshold in r^2 units
  abline(v = f0, lty = 2, col = "red")                                   # focal position
  s <- dt[pval < sig]; s[, near := is.finite(cM) & abs(cM - f0) <= 10]
  points(s[near == TRUE, cM],  s[near == TRUE, Rsq],  col = "#0072B2", pch = 1, cex = 1.3, lwd = 1.6)  # <=10cM
  points(s[near == FALSE, cM], s[near == FALSE, Rsq], col = "#E69F00", pch = 1, cex = 1.3, lwd = 1.6)  # >10cM
  legend("topright", bty = "n", cex = 0.7, lty = c(2, NA, NA), pch = c(NA, 1, 1),
         col = c("red", "#0072B2", "#E69F00"),
         legend = c("focal position", "sig <=10 cM (excluded)", "sig >10 cM (unlinked)"))
}

## per-focal zoomed multipage PDF (the "without 10 cM limit" local view for every focal)
pdf("Figures/moduleD_emmax_independence.pdf", width = 8, height = 4.2)
par(mar = c(4, 4, 2.5, 1), mgp = c(2.2, 0.6, 0))
for (g in focal$group_id) zoom(g)
dev.off()

## summary figure: a redundant example, a clean example, and the distribution
n_within10 <- vapply(focal$group_id, function(g) {
  ch <- fchr[g]; f0 <- fcm[g]; d <- E$results[[g]][Chr == ch & partner != g & pval < sig]
  sum(is.finite(d$cM) & abs(d$cM - f0) <= 10) }, integer(1))
clean_g <- focal$group_id[which(n_within10 == 0)][1]
redun_g <- focal$group_id[which.max(n_within10)]

png("Figures/moduleD_emmax_independence_summary.png", width = 2200, height = 780, res = 200)
par(mfrow = c(1, 3), mar = c(4, 4, 2.5, 1), mgp = c(2.2, 0.6, 0))
zoom(redun_g, sprintf("a  %s: local block split (%d sig clusters <=10 cM)", redun_g, max(n_within10)))
zoom(clean_g, sprintf("b  %s: independent (single local peak)", clean_g))
barplot(table(factor(n_within10, levels = 0:max(n_within10))), col = "grey70",
        xlab = "# same-chr significant clusters within 10 cM", ylab = "focal loci",
        main = "c  Local redundancy across focal loci")
dev.off()
cat(sprintf("wrote Figures/moduleD_emmax_independence.pdf (%d pages) + _summary.png | redundant=%s clean=%s\n",
            nrow(focal), redun_g, clean_g))
