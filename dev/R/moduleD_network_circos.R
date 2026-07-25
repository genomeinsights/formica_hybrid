## Illustrate the EMMAX trans/cis network and its hubs two ways: (1) a force-directed
## network (hub-and-spoke structure; hubs = high-degree nodes, e.g. F33028), and (2) a
## circular genome plot (chromosomes as sectors, a recombination-rate track, and links
## between significantly associated clusters coloured cis/trans) -- which shows the hubs
## sitting in the low-recombination centromeric dips and their cross-chromosome reach.
## Reads data/moduleD_emmax.rds, the clustering, recmap. Writes two PNGs to Figures/.
suppressPackageStartupMessages({ library(data.table); library(igraph); library(circlize); library(wesanderson) })
dir.create("Figures", showWarnings = FALSE)

E <- readRDS("data/moduleD_emmax.rds"); clust <- readRDS("data/eMLG_5loci_0025_cM05.rds"); groups <- clust$groups
e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1)
map <- as.data.table(e1$map_hyb_005)[, .(marker, Chr = as.character(Chr), Pos, ldw = ld_w_095)]
chr_of <- setNames(as.character(groups$Chr), groups$group_id)
## per-cluster physical position (median member Pos) and recombination rate
rec <- fread("data/Frufa_DTOL_PR.ref_genome.recmap"); setnames(rec, 1:4, c("chr","pos","cM","cMMb")); rec[, Chr := sub("chromosome_","Chr",chr)]
map[, rate := NA_real_]
for (ch in unique(map$Chr)) { r <- rec[Chr == ch]; if (nrow(r) < 2) next
  ix <- map[, which(Chr == ch)]; map[ix, rate := approx(r$pos, r$cMMb, xout = Pos, rule = 2)$y] }
mk2 <- groups[, .(marker = unlist(members)), by = group_id]
mk2 <- map[mk2, on = "marker"]
cpos <- mk2[, .(pos = median(Pos, na.rm = TRUE), start = min(Pos, na.rm = TRUE), end = max(Pos, na.rm = TRUE),
                rate = median(rate, na.rm = TRUE), ldw = median(ldw, na.rm = TRUE)), by = group_id]
pos_of <- setNames(cpos$pos, cpos$group_id)

## edges: clean, unlinked, deduplicated pairs, with cis/trans
h <- E$hits[paralog == FALSE, .(a = focal, b = partner, coupling)]
h[, key := ifelse(a < b, paste(a, b), paste(b, a))]; ed <- h[!duplicated(key)]
cat(sprintf("[net] %d clean unlinked edges (%d trans, %d cis) over %d clusters\n",
            nrow(ed), sum(ed$coupling == "repulsion"), sum(ed$coupling == "coupling"),
            uniqueN(c(ed$a, ed$b))))

## ---- (1) force-directed network -----------------------------------------
g <- graph_from_data_frame(ed[, .(a, b)], directed = FALSE)
vchr <- chr_of[V(g)$name]; deg <- as.integer(igraph::degree(g))
chr_lv <- paste0("Chr", sort(as.integer(sub("Chr", "", unique(vchr)))))
chr_col <- setNames(circlize::rand_color(length(chr_lv), luminosity = "bright"), chr_lv)
E(g)$color <- ifelse(ed$coupling == "repulsion", "#D55E00AA", "#0072B233")
set.seed(1); lay <- layout_with_fr(g)
png("Figures/moduleD_trans_network.png", width = 1900, height = 1500, res = 200)
par(mar = c(0, 0, 2, 0))
plot(g, layout = lay, vertex.color = chr_col[vchr],
     vertex.size = 2 + 2.4 * sqrt(deg), vertex.frame.color = NA,
     vertex.label = ifelse(deg >= 10, V(g)$name, NA), vertex.label.cex = 0.6,
     vertex.label.color = "black", vertex.label.dist = 0.4,
     edge.width = 0.6, main = "EMMAX trans/cis network (nodes = clusters, size ~ degree; red = trans, blue = cis)")
legend("bottomleft", legend = c("trans (repulsion)", "cis (coupling)"), col = c("#D55E00", "#0072B2"), lwd = 3, bty = "n", cex = 0.8)
dev.off()

## ---- (2) circular genome plot -------------------------------------------
chr_len <- map[, .(end = max(Pos, na.rm = TRUE)), by = Chr]
chr_len <- chr_len[Chr %in% chr_lv][order(match(Chr, chr_lv))]
chrdf <- data.frame(chr = chr_len$Chr, start = 0, end = chr_len$end)
## link beds (cluster positions), trans in orange, cis faint blue
lk <- ed[a %in% names(pos_of) & b %in% names(pos_of)]
b1 <- data.frame(chr = chr_of[lk$a], start = pos_of[lk$a], end = pos_of[lk$a] + 1)
b2 <- data.frame(chr = chr_of[lk$b], start = pos_of[lk$b], end = pos_of[lk$b] + 1)
lcol <- ifelse(lk$coupling == "repulsion", "#D55E0066", "#0072B222")

## ld_w as a scatter track in the style of supplementary Fig S4: per member marker of the
## has_eMLG clusters, y = local LD support (ld_w), coloured by LD-cluster (rolling palette).
he <- groups[has_eMLG == TRUE, group_id]
md <- mk2[group_id %in% he & is.finite(Pos) & is.finite(ldw), .(group_id, chr = Chr, pos = Pos, ldw)]
md <- md[chr %in% chr_lv]
rollpal <- grDevices::hcl.colors(12, "Dark 3")
cl_ord <- md[, .(chr = chr[1], mp = median(pos)), by = group_id][order(match(chr, chr_lv), mp)]
cl_ord[, roll := rollpal[(seq_len(.N) - 1) %% length(rollpal) + 1]]
md[cl_ord, on = "group_id", roll := i.roll]
setorder(md, chr, pos)

png("Figures/moduleD_trans_circos.png", width = 1700, height = 1700, res = 200)
circos.clear(); circos.par(gap.degree = 1.4, start.degree = 90, cell.padding = c(0, 0, 0, 0))
circos.genomicInitialize(chrdf, plotType = NULL)
## chromosome-number labels (radial, non-overlapping)
circos.track(ylim = c(0, 1), track.height = 0.04, bg.border = NA, panel.fun = function(x, y)
  circos.text(CELL_META$xcenter, 0.3, gsub("Chr", "", CELL_META$sector.index),
              cex = 0.6, facing = "clockwise", niceFacing = TRUE))
## ld_w scatter track (per marker, coloured by LD-cluster) -- as supplementary Fig S4
circos.genomicTrack(md[, .(chr, start = pos, end = pos, ldw, roll)], ylim = c(0, 1), track.height = 0.24, bg.border = "grey85",
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, numeric.column = 1, col = value$roll, pch = 16, cex = 0.08)
    if (CELL_META$sector.index == chr_lv[1])
      circos.yaxis(side = "left", at = c(0, 0.5, 1), labels.cex = 0.4, tick.length = 0.3) })
circos.genomicLink(b1, b2, col = lcol, lwd = 0.6)
title("Genome-wide trans (orange) / cis (blue) associations; track = ld_w per marker, coloured by LD-cluster",
      cex.main = 0.72, line = -0.5)
legend("bottomright", legend = c("trans", "cis"), col = c("#D55E00", "#0072B2"), lwd = 3, bty = "n", cex = 0.7)
dev.off()
cat("wrote Figures/moduleD_trans_network.png, Figures/moduleD_trans_circos.png\n")
