## One genotype heatmap PER cross-chromosome correlation MODULE (moduleD_network_build.R,
## meta_nodes$module: average-linkage on meta-node consensus, cut at |r|>MODULE_R). Each
## module is a coherent co-ancestry block, so its heatmap shows related SNPs only -- unlike
## the all-structure heatmap, which mixes several unrelated modules. Columns = the module's
## member-cluster SNPs (ward.D2, slice dendrogram); colour bars bottom->top: cluster, meta,
## sorting, DI; rows = individuals along the module's co-ancestry PC1 (each cluster polarised
## to its consensus then oriented to that PC1). Writes Figures/moduleD_module_<id>_heatmap.png
## for every module with >= MIN_META meta-nodes. Run from repo root.
suppressPackageStartupMessages({ library(data.table); library(ComplexHeatmap); library(circlize) })
dir.create("Figures", showWarnings = FALSE)
MIN_META <- 3                                              # only modules this size get a heatmap

N <- readRDS("data/moduleD_network.rds"); mn <- as.data.table(N$meta_nodes)
clust <- readRDS("data/eMLG_5loci_0025_cM05.rds"); groups <- clust$groups
cl <- readRDS("data/moduleC_C3_cl.rds")
e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1)
GT <- e1$GTs_hybrids_005; sdt <- as.data.table(e1$sample_data); pop <- setNames(sdt$Population, sdt$Sample_ID)
map <- as.data.table(e1$map_hyb_005)[, .(marker, Chr = as.character(Chr), Pos)]
DI_of <- setNames(cl$DI, cl$group_id); sc_of <- setNames(as.character(cl$sort_class), cl$group_id)
chrn <- function(c) as.integer(sub("Chr", "", c))
roll <- function(ids, pal) setNames(pal[(match(ids, unique(ids)) - 1) %% length(pal) + 1], ids)
sortcol <- c(aquilonia = "#0072B2", polyctena = "#D55E00", bidirectional = "#CC79A7", unsorted = "#E6E6E6")
di_col <- colorRamp2(c(-120, -40, -5), c("#D55E00", "#F7F7F7", "#0072B2"))

draw_module <- function(mod) {
  mm <- mn[module == mod][order(chrn(chr), cM)]
  memb <- rbindlist(lapply(seq_len(nrow(mm)), function(i)
    data.table(meta = mm$meta[i], cluster = strsplit(mm$members[i], ";")[[1]])))
  mk <- groups[.(memb$cluster), on = "group_id", .(marker = unlist(members)), by = group_id]
  setnames(mk, "group_id", "cluster"); mk <- map[mk, on = "marker"]
  mk[memb, on = "cluster", meta := i.meta]
  mk[, `:=`(sorting = sc_of[cluster], DI = DI_of[cluster], meta = factor(meta, levels = mm$meta))]
  mk[, clpos := median(Pos, na.rm = TRUE), by = cluster]; setorder(mk, meta, clpos, Pos)
  mk <- mk[marker %in% colnames(GT)]
  if (uniqueN(mk$cluster) < 2 || nrow(mk) < 4) { cat(sprintf("  module %d skipped (too small)\n", mod)); return(invisible()) }
  X <- GT[, mk$marker, drop = FALSE]; storage.mode(X) <- "double"
  consm <- sapply(split(mk$marker, mk$cluster), function(cols) rowMeans(X[, cols, drop = FALSE], na.rm = TRUE))
  Cim <- apply(consm, 2, function(v) { v[!is.finite(v)] <- mean(v, na.rm = TRUE); v })
  pc1 <- prcomp(Cim, center = TRUE, scale. = TRUE)$x[, 1]
  for (cc in unique(mk$cluster)) {
    cols <- mk[cluster == cc, marker]; cv <- rowMeans(X[, cols, drop = FALSE], na.rm = TRUE)
    r <- suppressWarnings(as.vector(cor(X[, cols, drop = FALSE], cv, use = "pairwise.complete.obs")))
    fl <- cols[which(r < 0)]; if (length(fl)) X[, fl] <- 2 - X[, fl]
    if (cor(cv, pc1, use = "pairwise.complete.obs") < 0) X[, cols] <- 2 - X[, cols]
  }
  row_ord <- order(pc1)
  top_ann <- HeatmapAnnotation(DI = mk$DI, sorting = mk$sorting, meta = as.character(mk$meta), cluster = mk$cluster,
    col = list(DI = di_col, sorting = sortcol, meta = roll(as.character(mk$meta), grDevices::hcl.colors(9, "Set 3")),
               cluster = roll(mk$cluster, grDevices::hcl.colors(8, "Dark 3"))),
    annotation_name_gp = gpar(fontsize = 8), show_legend = c(DI = TRUE, sorting = TRUE, meta = FALSE, cluster = FALSE),
    simple_anno_size = unit(4, "mm"), gap = unit(0.5, "mm"))
  pop_v <- pop[rownames(X)]; pop_col <- setNames(grDevices::hcl.colors(uniqueN(pop_v), "Spectral"), sort(unique(pop_v)))
  left_ann <- rowAnnotation(population = pop_v, col = list(population = pop_col),
    annotation_name_gp = gpar(fontsize = 8), simple_anno_size = unit(3, "mm"))
  f <- sprintf("Figures/moduleD_module_%02d_heatmap.png", mod)
  png(f, width = 14, height = 8.5, units = "in", res = 200)
  ht <- Heatmap(X[row_ord, ], name = "genotype", col = colorRamp2(c(0, 1, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
    na_col = "grey90", cluster_rows = FALSE, cluster_columns = TRUE, show_column_dend = TRUE,
    clustering_method_columns = "ward.D2", column_split = factor(mk$cluster, levels = unique(mk$cluster)),
    cluster_column_slices = TRUE, column_dend_height = unit(14, "mm"), column_gap = unit(0, "mm"), column_title = NULL,
    show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE, raster_quality = 2,
    top_annotation = top_ann, left_annotation = left_ann[row_ord, ], heatmap_legend_param = list(at = c(0, 1, 2)))
  draw(ht, merge_legend = TRUE, column_title_gp = gpar(fontsize = 11),
    column_title = sprintf("Module %d (|r|>%.1f): %d SNPs x %d individuals; %d clusters / %d meta-nodes; %s (chr %s)",
      mod, N$params$MODULE_R, nrow(mk), nrow(X), uniqueN(mk$cluster), nrow(mm),
      if (all(mm$structure)) "structure" else if (!any(mm$structure)) "candidate" else "mixed",
      paste(sort(unique(chrn(mm$chr))), collapse = ",")))
  dev.off(); cat(sprintf("  module %02d -> %s (%d SNPs, %d clusters)\n", mod, f, nrow(mk), uniqueN(mk$cluster)))
}
mods <- mn[, .N, by = module][N >= MIN_META, module][order(-mn[, .N, by = module][N >= MIN_META, N])]
cat(sprintf("modules with >= %d meta-nodes: %s\n", MIN_META, paste(sort(mods), collapse = ",")))
for (m in mods) draw_module(m)
cat("done\n")
