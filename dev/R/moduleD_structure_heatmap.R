## Genotype heatmap of the Module-D "structure" / co-ancestry network: all member SNPs of
## every structure-labelled meta-node (individuals x markers), to see whether these
## low-recombination hub clusters really are one co-inherited admixture-block axis across
## individuals. Column colour bars (bottom->top, nearest the heatmap first): original
## LD-cluster, region-merged meta-node, ancestry-sorting class, and DI. Row order and the
## whole-cluster orientation follow the shared co-ancestry axis (PC1 of the cluster-consensus
## matrix); each cluster is first polarised to its own consensus so within-cluster markers are
## coherent. Reads data/moduleD_network.rds + the clustering + moduleC gate + genotypes.
## Writes Figures/moduleD_structure_heatmap.png. Run from repo root.
suppressPackageStartupMessages({ library(data.table); library(ComplexHeatmap); library(circlize) })
dir.create("Figures", showWarnings = FALSE)

N <- readRDS("data/moduleD_network.rds"); mn <- as.data.table(N$meta_nodes)
clust <- readRDS("data/eMLG_5loci_0025_cM05.rds"); groups <- clust$groups
cl <- readRDS("data/moduleC_C3_cl.rds")
e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1)
GT <- e1$GTs_hybrids_005                                     # individuals x markers, 0/1/2
sdt <- as.data.table(e1$sample_data); pop <- setNames(sdt$Population, sdt$Sample_ID)
map <- as.data.table(e1$map_hyb_005)[, .(marker, Chr = as.character(Chr), Pos)]
DI_of <- setNames(cl$DI, cl$group_id); sc_of <- setNames(as.character(cl$sort_class), cl$group_id)

## ---- ordered marker list: meta (chr,cM) -> cluster (cM) -> marker (Pos) --------
sm <- mn[structure == TRUE][order(match(chr, paste0("Chr", sort(as.integer(sub("Chr","",unique(chr)))))), cM)]
memb <- rbindlist(lapply(seq_len(nrow(sm)), function(i)
  data.table(meta = sm$meta[i], cluster = strsplit(sm$members[i], ";")[[1]])))
cmv <- setNames(mn$cM, mn$meta)                             # cluster cM only exists for reps; use per-cluster below
mk <- groups[.(memb$cluster), on = "group_id", .(marker = unlist(members)), by = group_id]
setnames(mk, "group_id", "cluster"); mk <- map[mk, on = "marker"]
mk[memb, on = "cluster", meta := i.meta]
mk[, `:=`(sorting = sc_of[cluster], DI = DI_of[cluster])]
## order: meta order (as in sm) -> cluster median Pos -> marker Pos
mk[, meta := factor(meta, levels = sm$meta)]
mk[, clpos := median(Pos, na.rm = TRUE), by = cluster]
setorder(mk, meta, clpos, Pos)
mk <- mk[marker %in% colnames(GT)]
cat(sprintf("structure heatmap: %d individuals x %d markers over %d clusters / %d meta-nodes\n",
            nrow(GT), nrow(mk), uniqueN(mk$cluster), uniqueN(mk$meta)))

## ---- genotype block, polarised per cluster then oriented to the shared co-ancestry axis ----
X <- GT[, mk$marker, drop = FALSE]; storage.mode(X) <- "double"
cons <- sapply(split(mk$marker, mk$cluster), function(cols) rowMeans(X[, cols, drop = FALSE], na.rm = TRUE))
## PC1 of the 164 x nCluster consensus matrix = the dominant co-ancestry axis
Ci <- apply(cons, 2, function(v) { v[!is.finite(v)] <- mean(v, na.rm = TRUE); v })
pc1 <- prcomp(Ci, center = TRUE, scale. = TRUE)$x[, 1]
for (cc in unique(mk$cluster)) {                            # within-cluster polarise, then flip cluster to pc1
  cols <- mk[cluster == cc, marker]; cv <- rowMeans(X[, cols, drop = FALSE], na.rm = TRUE)
  r <- suppressWarnings(as.vector(cor(X[, cols, drop = FALSE], cv, use = "pairwise.complete.obs")))
  fl <- cols[which(r < 0)]; if (length(fl)) X[, fl] <- 2 - X[, fl]
  if (cor(cv, pc1, use = "pairwise.complete.obs") < 0) X[, cols] <- 2 - X[, cols]
}
row_ord <- order(pc1)                                       # individuals along the co-ancestry axis

## ---- column annotations (bottom = cluster, then meta, then sorting, top = DI) ---
roll <- function(ids, pal) setNames(pal[(match(ids, unique(ids)) - 1) %% length(pal) + 1], ids)
cl_lv <- unique(mk$cluster); mt_lv <- levels(droplevels(mk$meta))
cl_col <- roll(mk$cluster, grDevices::hcl.colors(8, "Dark 3"))
mt_col <- roll(as.character(mk$meta), grDevices::hcl.colors(9, "Set 3"))
sortcol <- c(aquilonia = "#0072B2", polyctena = "#D55E00", bidirectional = "#CC79A7", unsorted = "#E6E6E6")
mk[is.na(sorting) | !sorting %in% names(sortcol), sorting := "unsorted"]
di_col <- colorRamp2(c(-120, -40, -5), c("#D55E00", "#F7F7F7", "#0072B2"))   # low DI polyctena / high aquilonia
## HeatmapAnnotation stacks top->bottom in list order; put DI top, cluster nearest heatmap
top_ann <- HeatmapAnnotation(
  DI       = mk$DI,
  sorting  = mk$sorting,
  meta     = as.character(mk$meta),
  cluster  = mk$cluster,
  col = list(DI = di_col, sorting = sortcol, meta = mt_col, cluster = cl_col),
  annotation_name_gp = gpar(fontsize = 8),
  show_legend = c(DI = TRUE, sorting = TRUE, meta = FALSE, cluster = FALSE),
  simple_anno_size = unit(4, "mm"), gap = unit(0.5, "mm"))
pop_v <- pop[rownames(X)]; pop_col <- setNames(grDevices::hcl.colors(uniqueN(pop_v), "Spectral"), sort(unique(pop_v)))
left_ann <- rowAnnotation(population = pop_v, col = list(population = pop_col),
  annotation_name_gp = gpar(fontsize = 8), show_legend = TRUE, simple_anno_size = unit(3, "mm"))

## cluster the columns at the CLUSTER level: 87 blocks reordered by centroid similarity
## (slice dendrogram), SNPs kept contiguous & genomic within each block
col_split <- factor(mk$cluster, levels = unique(mk$cluster))
png("Figures/moduleD_structure_heatmap.png", width = 16, height = 9.5, units = "in", res = 200)
ht <- Heatmap(X[row_ord, ], name = "genotype", col = colorRamp2(c(0, 1, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
  na_col = "grey90", cluster_rows = FALSE, cluster_columns = TRUE, show_column_dend = TRUE,
  clustering_method_columns = "ward.D2",
  column_split = col_split, cluster_column_slices = TRUE,
  column_dend_height = unit(16, "mm"), column_gap = unit(0, "mm"), column_title = NULL,
  show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE, raster_quality = 2,
  top_annotation = top_ann, left_annotation = left_ann[row_ord, ],
  heatmap_legend_param = list(at = c(0, 1, 2)))
draw(ht, merge_legend = TRUE,
  column_title = sprintf("Module D structure / co-ancestry network: %d SNPs x %d individuals (rows ~ co-ancestry PC1; %d clusters ordered by correlation, slice dendrogram)", nrow(mk), nrow(X), uniqueN(mk$cluster)),
  column_title_gp = gpar(fontsize = 11))
dev.off()
cat("wrote Figures/moduleD_structure_heatmap.png\n")
