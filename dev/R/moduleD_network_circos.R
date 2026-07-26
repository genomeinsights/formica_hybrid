## Illustrate the Module D trans/cis association network built by moduleD_network_build.R
## (global FDR q<0.01 -> paralogy filter -> third-level within-region merge -> leverage
## filter, with low-recombination hubs LABELLED 'structure' rather than removed). Two views:
##  (1) a force-directed network (meta-nodes sized by degree, structure hubs greyed), and
##  (2) a circular genome plot where each association is drawn as a REGION-SPANNING BAND
##      (the full extent an LD-cluster/meta-node maps to) with opacity scaled by degree, so
##      large multi-connected hubs pop out and singletons recede. Tracks: ld_w>0.2 per
##      marker (local LD support) and the genome-wide sorting bar. The two heatmap hub
##      pairs are drawn as thick black bands.
## Reads data/moduleD_network.rds + the clustering + recmap. Writes two PNGs to Figures/.
suppressPackageStartupMessages({ library(data.table); library(igraph); library(circlize); library(wesanderson) })
dir.create("Figures", showWarnings = FALSE)

N  <- readRDS("data/moduleD_network.rds"); mn <- as.data.table(N$meta_nodes); me <- as.data.table(N$meta_edges)
clust <- readRDS("data/eMLG_5loci_0025_cM05.rds"); groups <- clust$groups
e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1)
map <- as.data.table(e1$map_hyb_005)[, .(marker, Chr = as.character(Chr), Pos, ldw = ld_w_095)]
chr_of <- setNames(as.character(groups$Chr), groups$group_id)
mk2  <- map[groups[, .(marker = unlist(members)), by = group_id], on = "marker"]
cpos <- mk2[, .(start = min(Pos, na.rm = TRUE), end = max(Pos, na.rm = TRUE)), by = group_id]
chr_lv <- paste0("Chr", sort(as.integer(sub("Chr", "", unique(mn$chr)))))

## ---- meta-network edges with degree-scaled opacity -----------------------
deg_of <- setNames(mn$degree, mn$meta)
me[, w := pmax(deg_of[ma], deg_of[mb])]                    # opacity driver = the hub end's degree
al <- function(x) { a <- pmin(0.95, pmax(0.15, 0.15 + 0.8 * (x - 1) / (max(me$w) - 1)))
  sprintf("%02X", as.integer(round(255 * a))) }
me[, base := ifelse(coupling == "trans", "#D55E00", "#0072B2")]
me[, col := paste0(base, al(w))]
## heatmap hub pairs -> their meta-nodes (drawn thick black on top)
meta_of_cluster <- function(h) mn$meta[grepl(paste0("(^|;)", h, "(;|$)"), mn$members)][1]
hub_pairs <- list(c("F27700","F33028"), c("F51717","F57289"))
hp <- rbindlist(lapply(hub_pairs, function(p) data.table(a = meta_of_cluster(p[1]), b = meta_of_cluster(p[2]))))
hp <- hp[!is.na(a) & !is.na(b)]
hp[, key := ifelse(a < b, paste(a, b), paste(b, a))]
me[, key := ifelse(ma < mb, paste(ma, mb), paste(mb, ma))]
hp <- hp[key %in% me$key]                                   # only pairs that survived the filters

cat(sprintf("[net] %d meta-edges (%d trans / %d cis) over %d meta-nodes; %d structure-labelled; %d heatmap hub pairs drawn\n",
            nrow(me), me[coupling=="trans",.N], me[coupling=="cis",.N], nrow(mn), mn[structure==TRUE,.N], nrow(hp)))

## ---- (1) force-directed network ------------------------------------------
g <- graph_from_data_frame(me[, .(ma, mb)], directed = FALSE, vertices = mn[, .(meta, chr, degree, structure)])
E(g)$color <- me$col
node_fill <- ifelse(V(g)$structure, "#BBBBBB", "#E69F00")   # grey = structure, gold = candidate region
png("Figures/moduleD_trans_network.png", width = 1900, height = 1500, res = 200)
par(mar = c(0, 0, 2, 0)); set.seed(1); lay <- layout_with_fr(g)
plot(g, layout = lay, vertex.color = node_fill, vertex.frame.color = NA,
     vertex.size = 3 + 2.6 * sqrt(as.integer(V(g)$degree)),
     vertex.label = ifelse(as.integer(V(g)$degree) >= 8, V(g)$name, NA), vertex.label.cex = 0.6,
     vertex.label.color = "black", vertex.label.dist = 0.4, edge.width = 1.1,
     main = "Module D meta-network (size ~ degree; grey = structure, gold = candidate)")
legend("bottomleft", legend = c("trans (repulsion)", "cis (coupling)", "structure region", "candidate region"),
       col = c("#D55E00", "#0072B2", "#BBBBBB", "#E69F00"), pch = c(NA, NA, 19, 19), lwd = c(3, 3, NA, NA), bty = "n", cex = 0.8)
dev.off()

## ---- (2) circular genome plots: region-spanning bands --------------------
chr_len <- map[Chr %in% chr_lv, .(end = max(Pos, na.rm = TRUE)), by = Chr][order(match(Chr, chr_lv))]
chrdf <- data.frame(chr = chr_len$Chr, start = 0, end = chr_len$end)
sp <- mn[, .(meta, chr, start, end, structure)]              # each meta-node's full mapped extent
hbA <- sp[hp, on = c(meta = "a")][, .(chr, start, end)]      # heatmap hub pairs (both structure)
hbB <- sp[hp, on = c(meta = "b")][, .(chr, start, end)]

## shared context tracks (identical on every plot): genome-wide sorting bar + ld_w>0.2 scatter
he <- groups[has_eMLG == TRUE, group_id]
cl <- readRDS("data/moduleC_C3_cl.rds"); sortcat <- setNames(as.character(cl$sort_class), cl$group_id)
bi_reps <- tryCatch(readRDS("data/moduleD_bidirectional.rds")$reps$group_id, error = function(e) character(0))
scat <- cpos[group_id %in% he]; scat[, `:=`(chr = chr_of[group_id], cat = sortcat[group_id])]
scat[group_id %in% bi_reps, cat := "bidirectional"]
scat[!cat %in% c("aquilonia","polyctena","bidirectional"), cat := "unsorted"]
scat <- scat[chr %in% chr_lv & is.finite(start) & is.finite(end)]
sortcol <- c(aquilonia = "#0072B2", polyctena = "#D55E00", bidirectional = "#CC79A7", unsorted = "#E6E6E6")
scat[, `:=`(col = sortcol[cat], uord = as.integer(cat != "unsorted"))]; setorder(scat, uord)
md <- mk2[group_id %in% he & is.finite(Pos) & is.finite(ldw) & ldw > 0.2, .(group_id, chr = Chr, pos = Pos, ldw)][chr %in% chr_lv]
rollpal <- grDevices::hcl.colors(12, "Dark 3")
cl_ord <- md[, .(chr = chr[1], mp = median(pos)), by = group_id][order(match(chr, chr_lv), mp)]
cl_ord[, roll := rollpal[(seq_len(.N) - 1) %% length(rollpal) + 1]]; md[cl_ord, on = "group_id", roll := i.roll]; setorder(md, chr, pos)

## draw one circos for an edge subset; band opacity ~ degree WITHIN that subset
## widen any band end narrower than min_span (bp) to that width, centred, so thin
## bands stay visible (band thickness in circos.genomicLink is the region's genomic span)
widen <- function(bed, min_span) { bed <- as.data.frame(bed); if (min_span <= 0) return(bed)
  s <- which((bed$end - bed$start) < min_span); mid <- (bed$start[s] + bed$end[s]) / 2
  bed$start[s] <- pmax(0, mid - min_span / 2); bed$end[s] <- mid + min_span / 2; bed }

draw_circos <- function(esub, file, subtitle, hub_beds = NULL, flat_alpha = FALSE, min_span = 0,
                        edge_col = NULL, node_col = NULL) {
  esub <- copy(esub); d2 <- table(c(esub$ma, esub$mb))
  wv <- pmax(d2[esub$ma], d2[esub$mb]); mx <- max(wv)
  alpha <- if (flat_alpha) rep(1, nrow(esub)) else pmin(0.95, pmax(0.18, 0.18 + 0.77 * (wv - 1) / max(1, mx - 1)))
  esub[, col := paste0(ifelse(coupling == "trans", "#D55E00", "#0072B2"),
                       sprintf("%02X", as.integer(round(255 * alpha))))]
  if (!is.null(edge_col)) esub[, col := edge_col]            # override cis/trans colour (e.g. by module)
  bA <- widen(sp[esub, on = c(meta = "ma")][, .(chr, start, end)], min_span)
  bB <- widen(sp[esub, on = c(meta = "mb")][, .(chr, start, end)], min_span)
  ns <- sp[meta %in% unique(c(esub$ma, esub$mb)) & chr %in% chr_lv]
  nsub <- data.table(chr = ns$chr, start = ns$start, end = ns$end,
    col = if (!is.null(node_col)) node_col[ns$meta] else ifelse(ns$structure, "#BBBBBB", "#E69F00"))
  png(file, width = 1700, height = 1700, res = 200)
  circos.clear(); circos.par(gap.degree = 1.4, start.degree = 90, cell.padding = c(0, 0, 0, 0))
  circos.genomicInitialize(chrdf, plotType = NULL)
  circos.track(ylim = c(0, 1), track.height = 0.04, bg.border = NA, panel.fun = function(x, y)
    circos.text(CELL_META$xcenter, 0.3, gsub("Chr", "", CELL_META$sector.index), cex = 0.6, facing = "clockwise", niceFacing = TRUE))
  circos.genomicTrack(md[, .(chr, start = pos, end = pos, ldw, roll)], ylim = c(0.2, 1), track.height = 0.08, bg.border = "grey85",
    panel.fun = function(region, value, ...) { circos.genomicPoints(region, value, numeric.column = 1, col = value$roll, pch = 16, cex = 0.1)
      if (CELL_META$sector.index == chr_lv[1]) circos.yaxis(side = "left", at = c(0.2, 0.6, 1), labels.cex = 0.4, tick.length = 0.3) })
  circos.genomicTrack(scat[, .(chr, start, end, col)], ylim = c(0, 1), track.height = 0.05, bg.border = "grey85",
    panel.fun = function(region, value, ...) circos.genomicRect(region, value, ytop = 1, ybottom = 0, col = value$col, border = NA))
  circos.genomicTrack(nsub, ylim = c(0, 1), track.height = 0.035, bg.border = "grey85",
    panel.fun = function(region, value, ...) circos.genomicRect(region, value, ytop = 1, ybottom = 0, col = value$col, border = NA))
  ordv <- order(wv)                                          # faint (low-degree) first, hubs on top
  circos.genomicLink(bA[ordv, ], bB[ordv, ], col = esub$col[ordv], border = NA)
  if (!is.null(hub_beds)) circos.genomicLink(hub_beds$a, hub_beds$b, col = "black", lwd = 2.4, border = NA)
  title(subtitle, cex.main = 0.64, line = -0.6)
  if (is.null(edge_col)) {                                   # cis/trans + structure/candidate legends
    legend("bottomright", legend = c("trans", "cis", if (!is.null(hub_beds)) "heatmap hub"),
           col = c("#D55E00", "#0072B2", if (!is.null(hub_beds)) "black"), lwd = c(3, 3, if (!is.null(hub_beds)) 2.4), bty = "n", cex = 0.7)
    legend("topleft", legend = c("structure (low-recomb)", "candidate"), fill = c("#BBBBBB", "#E69F00"), border = NA, bty = "n", cex = 0.7, title = "meta-node")
  } else {
    legend("bottomright", legend = "within-module link (colour = module)", col = "grey40", lwd = 3, bty = "n", cex = 0.7)
  }
  legend("bottomleft", legend = c("aquilonia", "polyctena", "bidirectional", "unsorted"), title = "sorting",
         fill = sortcol[c("aquilonia","polyctena","bidirectional","unsorted")], border = NA, bty = "n", cex = 0.7)
  dev.off()
}

str_of <- setNames(mn$structure, mn$meta); me[, `:=`(a_str = str_of[ma], b_str = str_of[mb])]
hub_beds <- if (nrow(hp)) list(a = hbA, b = hbB) else NULL
draw_circos(me, "Figures/moduleD_trans_circos.png",
  "Module D meta-network: all association bands, opacity ~ degree (grey = structure, gold = candidate)", hub_beds)
draw_circos(me[a_str == FALSE & b_str == FALSE], "Figures/moduleD_candidate_circos.png",
  "Module D candidate loci: F33454 aquilonia co-sorting module (Chr10 anchor) + F11431-F49480 trans (Chr3-16)",
  flat_alpha = TRUE, min_span = 3e6)
draw_circos(me[a_str == TRUE & b_str == TRUE], "Figures/moduleD_structure_circos.png",
  "Module D structure / co-ancestry module: low-recomb admixture-block association bands", hub_beds)

## within-module circos: only links whose two meta-nodes are in the SAME correlation module,
## coloured by module (so each co-ancestry module reads as its own coherent block of links)
mod_of <- setNames(mn$module, mn$meta); me[, `:=`(moda = mod_of[ma], modb = mod_of[mb])]
wm <- me[moda == modb]
modpal <- grDevices::hcl.colors(12, "Dark 3"); mods_u <- sort(unique(wm$moda))
modcol <- setNames(modpal[(seq_along(mods_u) - 1) %% length(modpal) + 1], mods_u)
node_modcol <- setNames(modcol[as.character(mn$module)], mn$meta)
draw_circos(wm, "Figures/moduleD_module_circos.png",
  sprintf("Module D within-module links only (|r|>%.1f modules): %d of %d links; node & link colour = module",
          N$params$MODULE_R, nrow(wm), nrow(me)),
  flat_alpha = TRUE, min_span = 2e6,
  edge_col = paste0(modcol[as.character(wm$moda)], "CC"), node_col = node_modcol)
cat(sprintf("wrote network + 4 circos (combined, candidate, structure, within-module: %d/%d links)\n", nrow(wm), nrow(me)))
