## Visualize the strongest EMMAX trans-pairs (unlinked clusters that sort in OPPOSITE
## ancestry directions -- repulsion, negative structure-corrected sign) as genotype
## heatmaps via LDscnR::plot_genotype_heatmap(). For a trans pair the two clusters
## must render as MIRROR images, so polarize = FALSE and cluster_columns = FALSE (the
## default polarisation would flip the anti-correlated cluster and hide the trans
## signal); individuals are ordered by the focal cluster's consensus dosage so the
## opposite gradients are obvious. Reads data/moduleD_emmax.rds. Writes PNGs to Figures/.
suppressPackageStartupMessages({ library(data.table); devtools::load_all("~/gitlab/LDscnR/", quiet = TRUE) })

E     <- readRDS("data/moduleD_emmax.rds")
clust <- readRDS("data/eMLG_5loci_0025_cM05.rds"); groups <- clust$groups; eMLG <- clust$eMLG
load("data/hybrids_only_maf005.Rdata")                 # GTs_hybrids_005, sample_data, ...
sd <- as.data.table(sample_data)
pop_by_ind <- setNames(sd$Population, sd$Sample_ID)
chr_of <- setNames(as.character(groups$Chr), groups$group_id)
mem_of <- function(g) groups[.(g), on = "group_id", members][[1]]

pop_colors <- c(
  Aland = "#E41A1C", Bunkkeri = "#377EB8", Grundsund = "#4DAF4A",
  Heinamaki = "#984EA3", Hiivola = "#FF7F00", Jarvenpaa = "#FFFF33",
  Karsikas = "#A65628", Katiskoski = "#F781BF", Kummunmaki = "#999999",
  LangholmenR = "#66C2A5", LangholmenW = "#FC8D62", Nyrhispera74 = "#8DA0CB",
  Nyrhispera75 = "#E78AC3", Parikkala = "#A6D854", Pikkala = "#FFD92F",
  Sielva = "#E5C494", Svanvik1 = "#B3B3B3", Svanvik2 = "#1B9E77",
  Tvarminne = "#D95F02", Vuosaari = "#7570B3")

## strongest clean trans pairs (deduplicated on the unordered pair)
tr <- E$hits[paralog == FALSE & coupling == "repulsion"][order(-F)]
tr[, key := ifelse(focal < partner, paste(focal, partner), paste(partner, focal))]
tr <- tr[!duplicated(key)]
top <- head(tr, 4)

## orient each marker to POSITIVELY correlate with its cluster's consensus, so each
## block is an internally coherent stripe while the A-vs-B (trans) relationship is
## preserved -- unlike plot_genotype_heatmap's global polarisation, which would flip
## the whole B block to match A and erase the trans signal.
polarize_to_ref <- function(M, ref) {
  for (j in seq_len(ncol(M))) {
    cc <- suppressWarnings(cor(M[, j], ref, use = "pairwise.complete.obs"))
    if (is.finite(cc) && cc < 0) M[, j] <- 2 - M[, j]
  }
  M
}

for (i in seq_len(nrow(top))) {
  A <- top$focal[i]; B <- top$partner[i]
  mkA <- intersect(mem_of(A), colnames(GTs_hybrids_005))
  mkB <- intersect(mem_of(B), colnames(GTs_hybrids_005))
  mk  <- c(mkA, mkB)
  GTa <- polarize_to_ref(GTs_hybrids_005[, mkA, drop = FALSE], eMLG[, A])
  GTb <- polarize_to_ref(GTs_hybrids_005[, mkB, drop = FALSE], eMLG[, B])
  GT  <- cbind(GTa, GTb)
  cl_lab <- setNames(c(rep(sprintf("%s (%s)", A, chr_of[A]), length(mkA)),
                       rep(sprintf("%s (%s)", B, chr_of[B]), length(mkB))), mk)
  cl_cols <- setNames(c("#1B9E77", "#D95F02"),
                      c(sprintf("%s (%s)", A, chr_of[A]), sprintf("%s (%s)", B, chr_of[B])))
  ## order individuals by the focal cluster's RAW member-mean dosage (matches the
  ## displayed genotypes' polarity, unlike the eMLG consensus) so cluster A renders as
  ## a clean gradient and a trans partner B as its mirror image.
  aMean <- rowMeans(GTa, na.rm = TRUE)
  ord <- rownames(GT)[order(aMean, na.last = TRUE)]
  r   <- cor(eMLG[, A], eMLG[, B], use = "pairwise.complete.obs")
  out <- sprintf("Figures/moduleD_trans_%02d_%s_%s.png", i, A, B)
  plot_genotype_heatmap(
    GT, polarize = FALSE, cluster_columns = FALSE, row_order = ord,
    col_annotation = cl_lab[mk], col_annotation_name = "cluster (chr)",
    col_annotation_colors = cl_cols, col_annotation_legend = TRUE,
    row_annotation = pop_by_ind[ord], row_annotation_name = "Population",
    row_annotation_colors = pop_colors, row_annotation_legend = TRUE,
    legend_name = "dosage (cluster-oriented)",
    title = sprintf("trans pair: %s (%s) vs %s (%s)   r = %.2f, F = %.0f",
                    A, chr_of[A], B, chr_of[B], r, top$F[i]),
    out_file = out, width = 11, height = 7)
  message(sprintf("[%d] %s (%s, %d loci) vs %s (%s, %d loci) | r=%.2f -> %s",
                  i, A, chr_of[A], length(mkA), B, chr_of[B], length(mkB), r, out))
}
cat("done: strongest trans-pair heatmaps in Figures/\n")
