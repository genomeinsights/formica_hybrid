## =========================================================
## MODULE B -- eMLG-level climate-association Manhattan + eMLG-only findings
## =========================================================
## One point per LD-cluster eMLG (primary config: Aland-excluded, fixed
## LD-pruned Omega), placed at its representative marker's position, coloured by
## how much INDIVIDUAL member-SNP support the cluster has among the eMLGs that
## themselves reach BF(dB) >= EMLG_THR:
##   "eMLG-only"  -- eMLG significant but ZERO member SNPs reach BF(dB) >= SIG_THR
##   "1-4 sig SNPs"
##   ">=5 sig SNPs (candidate)" -- the member-count candidate set (moduleB_climate_candidates.R)
## The eMLG-only class is the point of the plot: coherent cluster-level climate
## association that the member-count criterion (>=5 significant members) cannot
## detect, because a small cluster can never reach 5 significant members. These
## are genuine cluster-level signals recovered by aggregating to the consensus.
##
## If the sim-FDR results (moduleB_BayPass/data/moduleB_eMLG_null.rds, produced by
## moduleB_eMLG_null.R -- runs AFTER this script in the pipeline) are present,
## two extra things are drawn: (a) in the eMLG Manhattan the FDR floor-survivors
## (k = 0 of 10,000 structure nulls) are marked as LARGER TRIANGLES; (b) a
## separate SNP-level Manhattan over ALL SNPs where only the member SNPs of the
## FDR-passing eMLG clusters are coloured (a unique colour per cluster, as in the
## Fig. 2 SNP Manhattans). Both are skipped with a note if the null is absent, so
## re-run this script after moduleB_eMLG_null.R to regenerate them.
##
## LD-cluster / eMLG association approach: Li Z, Kemppainen P, Rastas P, Merila J.
## 2018. Linkage disequilibrium clustering-based approach for association mapping
## with tightly linked genomewide data. Mol Ecol Resour 18:809-824.
##
## Reads : data/eMLG_5loci_0025_cM05.rds, data/hybrids_only_maf005.Rdata (map),
##         aland_excluded/PC{1,2}_DIEM_withOmega_summary_betai_reg.out (member SNP BF),
##         aland_excluded_eMLG/PC{1,2}_eMLG_withOmega_summary_betai_reg.out (eMLG BF),
##         aland_excluded_eMLG/eMLG_group_order.txt,
##         moduleB_BayPass/data/moduleB_eMLG_null.rds  (OPTIONAL -- FDR floor flags for the overlays)
## Writes: moduleB_BayPass/data/moduleB_eMLG_association.rds, moduleB_BayPass/data/moduleB_eMLG_outliers.csv,
##         moduleB_BayPass/Figures/moduleB_eMLG_manhattan.png       (Fig 3; FDR survivors = triangles),
##         moduleB_BayPass/Figures/moduleB_fdr_snp_manhattan.png    (all SNPs; FDR clusters coloured)
## Run from the formica_hybrid repo root.

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
SIG_THR <- 15; EMLG_THR <- 15; MIN5 <- 5
PRIM_POP <- "aland_excluded"; PRIM_OM <- "withOmega"; D <- "aland_excluded_eMLG"
dir.create("moduleB_BayPass/Figures", showWarnings = FALSE)

## ---- clustering + representative positions ------------------------------
g  <- readRDS("data/eMLG_5loci_0025_cM05.rds")$groups
he <- g[has_eMLG == TRUE]
e <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e)
map <- as.data.table(e$map_hyb_005)[, .(marker, Chr, Pos)]; mk <- map$marker; rm(e); invisible(gc())
pos <- map[he[, .(group_id, marker = representative)], on = "marker"][, .(group_id, Chr, Pos)]

## ---- member-SNP significance per cluster (vectorised: lookup OUTSIDE by) --
m2g <- he[, .(marker = unlist(members)), by = group_id]
imp <- function(pc) { r <- fread(sprintf("./%s/%s_DIEM_%s_summary_betai_reg.out", PRIM_POP, pc, PRIM_OM))
  stopifnot(nrow(r) == length(mk), identical(r$MRK, seq_len(nrow(r)))); setNames(r$`BF(dB)`, mk) }
b1 <- imp("PC1"); b2 <- imp("PC2")
m2g[, `:=`(s1 = b1[marker] >= SIG_THR, s2 = b2[marker] >= SIG_THR)]
memb <- m2g[, .(size = .N, nsig1 = sum(s1, na.rm = TRUE), nsig2 = sum(s2, na.rm = TRUE)), by = group_id]

## ---- eMLG-level BF, mapped to group_id ---------------------------------
grp <- readLines(file.path(D, "eMLG_group_order.txt"))
rd <- function(pc) { r <- fread(sprintf("%s/%s_eMLG_%s_summary_betai_reg.out", D, pc, PRIM_OM))
  stopifnot(nrow(r) == length(grp), identical(r$MRK, seq_len(nrow(r)))); setNames(r$`BF(dB)`, grp) }
eBF1 <- rd("PC1"); eBF2 <- rd("PC2")

dt <- pos[memb, on = "group_id"]
dt[, `:=`(eBF1 = eBF1[group_id], eBF2 = eBF2[group_id])]
cat_axis <- function(eBF, nsig) fifelse(eBF < EMLG_THR, "ns",
                          fifelse(nsig == 0, "eMLG-only (0 sig SNPs)",
                          fifelse(nsig < MIN5, "1-4 sig SNPs", ">=5 sig SNPs (candidate)")))
dt[, `:=`(cat1 = cat_axis(eBF1, nsig1), cat2 = cat_axis(eBF2, nsig2))]

## ---- OPTIONAL: sim-FDR floor-survivor flags (from the downstream null) --
NULLF <- "moduleB_BayPass/data/moduleB_eMLG_null.rds"; have_fdr <- file.exists(NULLF)
if (have_fdr) {
  nul <- readRDS(NULLF)[, .(group_id, floor1, floor2)]
  dt <- nul[dt, on = "group_id"]
  dt[is.na(floor1), floor1 := FALSE][is.na(floor2), floor2 := FALSE]
  cat(sprintf("sim-FDR floor-survivors: PC1 %d, PC2 %d\n", dt[floor1 == TRUE, .N], dt[floor2 == TRUE, .N]))
} else {
  dt[, `:=`(floor1 = FALSE, floor2 = FALSE)]
  message("note: ", NULLF, " not found -- FDR overlays/figure skipped (run moduleB_eMLG_null.R, then re-run this).")
}

## ---- save the association table + the outlier-eMLG summary --------------
saveRDS(dt, "moduleB_BayPass/data/moduleB_eMLG_association.rds")
mkout <- function(bfcol, nsigcol, catcol, ax)
  dt[get(bfcol) >= EMLG_THR, .(group_id, Chr, Pos, size, axis = ax,
      eMLG_BF = get(bfcol), n_sig = get(nsigcol),
      pct_sig = round(100 * get(nsigcol) / size, 1), support = get(catcol))]
outl <- rbind(mkout("eBF1", "nsig1", "cat1", "PC1"),
              mkout("eBF2", "nsig2", "cat2", "PC2"))[order(axis, -eMLG_BF)]
fwrite(outl, "moduleB_BayPass/data/moduleB_eMLG_outliers.csv")
cat(sprintf("outlier eMLGs (BF>=%d): PC1 %d, PC2 %d  |  of which eMLG-only (0 sig SNPs): PC1 %d, PC2 %d\n",
            EMLG_THR, outl[axis == "PC1", .N], outl[axis == "PC2", .N],
            outl[axis == "PC1" & n_sig == 0, .N], outl[axis == "PC2" & n_sig == 0, .N]))

## ---- shared genome coordinate transform (from ALL markers) -------------
mapx <- copy(map)
mapx[, chr_num := suppressWarnings(as.integer(gsub("[^0-9]", "", Chr)))]
clen <- mapx[, .(len = max(Pos)), by = .(Chr, chr_num)][order(chr_num)]
spg  <- 0.01 * sum(clen$len)
clen[, start := c(0, head(cumsum(len + spg), -1))]
mapx <- clen[, .(Chr, start)][mapx, on = "Chr"]; mapx[, x := Pos + start]
xkey <- setNames(mapx$x, mapx$marker)          # marker -> genome x
cmid <- clen[, .(Chr, mid = start + len / 2)]
shade <- clen[chr_num %% 2 == 0, .(xmin = start, xmax = start + len)]
chr_lab <- function() scale_x_continuous(breaks = cmid$mid, labels = gsub("^Chr", "", cmid$Chr), expand = c(0.01, 0.01))
dt <- clen[, .(Chr, start)][dt, on = "Chr"]; dt[, x := Pos + start]  # eMLG x = representative position

## ---- Fig 3: eMLG Manhattan (FDR floor-survivors = larger triangles) -----
pal <- c("ns" = "grey78", "eMLG-only (0 sig SNPs)" = "#D81B60",
         "1-4 sig SNPs" = "#1E88E5", ">=5 sig SNPs (candidate)" = "#000000")
lev <- c("ns", "1-4 sig SNPs", "eMLG-only (0 sig SNPs)", ">=5 sig SNPs (candidate)")
panel <- function(bfcol, catcol, floorcol, lab) {
  d <- dt[, .(x, BF = get(bfcol), cat = factor(get(catcol), levels = lev), fl = get(floorcol))]
  setorder(d, cat)
  ggplot() +
    geom_rect(data = shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey93") +
    geom_point(data = d[cat == "ns"], aes(x, BF), color = "grey78", size = 0.5) +
    geom_point(data = d[cat != "ns" & !fl], aes(x, BF, color = cat), size = 1.7) +
    ## FDR floor-survivors: larger triangles with a dark outline, on top
    ## (kept out of the colour legend -- the caption explains the shape)
    geom_point(data = d[fl == TRUE], aes(x, BF, color = cat), shape = 17, size = 3.6, show.legend = FALSE) +
    geom_point(data = d[fl == TRUE], aes(x, BF), shape = 2, size = 3.6, color = "grey15", stroke = 0.5, show.legend = FALSE) +
    geom_hline(yintercept = EMLG_THR, linetype = 2, color = "red", linewidth = 0.4) +
    scale_color_manual(values = pal, name = NULL, drop = FALSE) +
    chr_lab() +
    labs(x = if (lab == "PC2") "Chromosome" else NULL, y = sprintf("eMLG BF(dB) -- %s", lab)) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "top")
}
p <- panel("eBF1", "cat1", "floor1", "PC1") / panel("eBF2", "cat2", "floor2", "PC2") +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
if (have_fdr) p <- p + plot_annotation(
  caption = "Triangles = FDR floor-survivors (real eMLG BF exceeds all 10,000 population-structure nulls).",
  theme = theme(plot.caption = element_text(size = 9, hjust = 0)))
ggsave("moduleB_BayPass/Figures/moduleB_eMLG_manhattan.png", p, width = 15, height = 7, dpi = 150)

## ---- NEW: all-SNP Manhattan, only FDR-passing eMLG clusters coloured ----
if (have_fdr) {
  fdr_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628",
               "#F781BF", "#1B9E77", "#D95F02", "#7570B3")
  snp_panel <- function(bf, fdr_groups, lab) {
    bg <- data.table(x = mapx$x, BF = bf[mapx$marker])                       # ALL SNPs
    fg <- m2g[group_id %in% fdr_groups, .(marker, group_id)]
    fg[, `:=`(x = xkey[marker], BF = bf[marker], grp = factor(group_id, levels = fdr_groups))]
    cols <- setNames(fdr_pal[seq_along(fdr_groups)], fdr_groups)
    ggplot() +
      geom_rect(data = shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey93") +
      geom_point(data = bg, aes(x, BF), color = "grey80", size = 0.3, alpha = 0.5) +
      geom_point(data = fg, aes(x, BF, color = grp), size = 1.3) +
      geom_hline(yintercept = SIG_THR, linetype = 2, color = "red", linewidth = 0.4) +
      scale_color_manual(values = cols, name = "FDR eMLG", drop = FALSE) +
      chr_lab() +
      labs(x = if (lab == "PC2") "Chromosome" else NULL, y = sprintf("SNP BF(dB) -- %s", lab)) +
      theme_bw(base_size = 11) +
      theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "right")
  }
  fdr1 <- dt[floor1 == TRUE, group_id]; fdr2 <- dt[floor2 == TRUE, group_id]
  ps <- snp_panel(b1, fdr1, "PC1") / snp_panel(b2, fdr2, "PC2") +
    plot_annotation(title = "All SNPs (grey); member SNPs of FDR-passing eMLG clusters coloured (one colour per cluster)",
                    theme = theme(plot.title = element_text(size = 11)))
  ggsave("moduleB_BayPass/Figures/moduleB_fdr_snp_manhattan.png", ps, width = 15, height = 7, dpi = 150)
}

cat("\nWrote moduleB_BayPass/data/moduleB_eMLG_association.rds, moduleB_BayPass/data/moduleB_eMLG_outliers.csv, moduleB_BayPass/Figures/moduleB_eMLG_manhattan.png",
    if (have_fdr) ", moduleB_BayPass/Figures/moduleB_fdr_snp_manhattan.png" else "", "\n", sep = "")
