## =========================================================================
## DI25 among-region EMMAX -- candidate pairwise-LD CIRCOS.
## =========================================================================
## All candidates (bidirectional, bi_score >= 0.2) placed on the genome; every UNLINKED
## candidate pair with r^2 > EDGE_R2 drawn as a faint chord, and the pooled-FDR-significant
## pairs (r^2 >= r2crit) highlighted by coupling sign (coupling = blue, repulsion = red;
## paralog-flagged = grey dashed). The faint background shows the diffuse residual
## co-ancestry (lambda = 1.34) among the sorted candidates; the few bright chords are the
## only pairs that survive the whole-experiment FDR -- the leads carried to Module E.
##
## Reads di25D_emmax.rds (pairs, edges, candidates, params) + di25D_units.rds (marker
## positions). Writes Figures/di25D_emmax_candidate_circos.png.

suppressPackageStartupMessages({ library(data.table); library(circlize) })
EMMAX  <- "module_di25/among_region_association/data/di25D_emmax.rds"
NET    <- "module_di25/among_region_association/data/di25D_network.rds"
UNITS  <- "module_di25/among_region_association/data/di25D_units.rds"
OUT    <- "module_di25/among_region_association/Figures/di25D_emmax_candidate_circos.png"
EDGE_R2 <- 0.10       # draw a chord for every unlinked candidate pair above this r^2

E <- readRDS(EMMAX); u <- readRDS(UNITS); Nw <- readRDS(NET)
cand <- E$candidates; r2crit <- E$params$r2crit; lambda <- E$params$lambda
## leverage-surviving leads (di25D_network_build.R, within-population LOO)
surv_keys <- paste(pmin(Nw$meta_edges$ma, Nw$meta_edges$mb), pmax(Nw$meta_edges$ma, Nw$meta_edges$mb))
g <- u$groups
pos_of <- setNames(as.integer(sub(".*:", "", g$representative)), g$group_id)   # representative bp
chr_of <- setNames(as.character(g$Chr), g$group_id)
cand[, pos := pos_of[group_id]]

## chromosome order + lengths (max representative position on each chr)
chr_lv <- paste0("Chr", sort(as.integer(sub("Chr", "", unique(chr_of)))))
chr_lv <- chr_lv[chr_lv %in% chr_of]
chrlen <- g[, .(end = max(pos_of[group_id], na.rm = TRUE)), by = .(chr = as.character(Chr))]
chrdf  <- data.frame(chr = chr_lv, start = 0, end = chrlen[match(chr_lv, chr), end])

## edge sets ---------------------------------------------------------------
p <- E$pairs[unlinked == TRUE & Rsq > EDGE_R2]
sig <- E$edges                                   # FDR-significant (may include a paralog)
sig_keys <- paste(pmin(sig$cluster1, sig$cluster2), pmax(sig$cluster1, sig$cluster2))
p[, key2 := paste(pmin(focal, partner), pmax(focal, partner))]
bg <- p[!key2 %in% sig_keys]                     # r^2>EDGE_R2 but not FDR-significant
message(sprintf("[circos] %d candidates | %d chords r^2>%.2f (%d background + %d FDR); lambda=%.2f, r2crit=%.3f",
                nrow(cand), nrow(p), EDGE_R2, nrow(bg), nrow(sig), lambda, r2crit))

## node colour by sorting direction
dircol <- c(aquilonia = "#0072B2", polyctena = "#D55E00", unresolved = "#999999")
cand[, ncol := dircol[as.character(direction)]]; cand[is.na(ncol), ncol := "#999999"]

png(OUT, width = 1800, height = 1800, res = 200)
circos.clear(); circos.par(gap.degree = 1.2, start.degree = 90, cell.padding = c(0, 0, 0, 0))
circos.genomicInitialize(chrdf, plotType = NULL)
## chromosome labels
circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA, panel.fun = function(x, y)
  circos.text(CELL_META$xcenter, 0.35, gsub("Chr", "", CELL_META$sector.index),
              cex = 0.6, facing = "clockwise", niceFacing = TRUE))
## candidate node track: a point per candidate, coloured by direction
circos.genomicTrack(cand[, .(chr = Chr, start = pos, end = pos, y = 0.5, ncol)], ylim = c(0, 1),
  track.height = 0.06, bg.border = "grey85",
  panel.fun = function(region, value, ...)
    circos.genomicPoints(region, value, numeric.column = "y", pch = 16, cex = 0.28, col = value$ncol))
## chords = point-to-point curves (circos.link draws a real line, unlike the ribbon that
## circos.genomicLink makes between two narrow regions). Faint grey background first, then
## the FDR-significant chords on top coloured by coupling (paralog = grey dashed).
draw_links <- function(d, col, lwd, lty = 1) if (nrow(d) > 0) for (i in seq_len(nrow(d)))
  circos.link(chr_of[d[[1]][i]], pos_of[d[[1]][i]], chr_of[d[[2]][i]], pos_of[d[[2]][i]],
              col = col[i], lwd = lwd, lty = lty[i])
if (nrow(bg) > 0)
  draw_links(bg[, .(focal, partner)], col = rep("#88888838", nrow(bg)), lwd = 0.6, lty = rep(1, nrow(bg)))
## FDR edges: leverage-SURVIVING leads bright (by coupling); FDR edges dropped by the
## within-population leverage (or paralog) drawn thin grey dashed.
if (nrow(sig) > 0) {
  sig[, skey := paste(pmin(cluster1, cluster2), pmax(cluster1, cluster2))]
  sig[, survivor := skey %in% surv_keys]
  drp <- sig[survivor == FALSE]; sv <- sig[survivor == TRUE]
  if (nrow(drp) > 0) draw_links(drp[, .(cluster1, cluster2)],
                                col = rep("#AAAAAA", nrow(drp)), lwd = 1.3, lty = rep(2, nrow(drp)))
  if (nrow(sv) > 0) draw_links(sv[, .(cluster1, cluster2)],
                               col = ifelse(sv$coupling == "coupling", "#0072B2", "#D55E00"),
                               lwd = 3.0, lty = rep(1, nrow(sv)))
}
title(sprintf("DI25 candidate pairwise LD (n=%d, bi_score>=%.1f): %d chords r2>%.2f, %d FDR, %d within-pop leads\nfaint = residual co-ancestry (lambda=%.2f); bright = leverage-surviving leads for Module E",
              nrow(cand), E$params$BI_MIN, nrow(p), EDGE_R2, nrow(sig), length(surv_keys), lambda), cex.main = 0.6, line = -1.2)
legend("bottomleft", bty = "n", cex = 0.65, title = "candidate direction",
       fill = dircol[c("aquilonia", "polyctena", "unresolved")], border = NA,
       legend = c("aquilonia", "polyctena", "unresolved"))
legend("bottomright", bty = "n", cex = 0.65, title = "leverage-surviving leads",
       col = c("#0072B2", "#D55E00", "#AAAAAA"), lwd = c(3, 3, 1.3), lty = c(1, 1, 2),
       legend = c("coupling", "repulsion", "FDR dropped by leverage"))
dev.off()
cat(sprintf("\n[done] %s | %d chords (%d FDR, %d within-pop leads) over %d candidates\n",
            OUT, nrow(p), nrow(sig), length(surv_keys), nrow(cand)))
