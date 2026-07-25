## Isolate and visualize the PAIR-SPECIFIC (non-hub) trans candidates from the
## PC-conditioned EMMAX scan. Most trans signal runs through hubs -- one focal locus
## anti-sorting with a whole co-varying module (a secondary structure axis / candidate
## chromosomal rearrangement, e.g. F33028). A genuine pair-specific incompatibility is
## instead a trans edge whose BOTH endpoints have low connectivity: they associate
## with each other but not with a module. This script scores each cluster's degree in
## the trans network, keeps the low-degree ("non-hub") trans pairs, and heatmaps them.
## Reads data/moduleD_emmax.rds. Writes a ranked table + PNGs to Figures/.
suppressPackageStartupMessages({ library(data.table); devtools::load_all("~/gitlab/LDscnR/", quiet = TRUE) })

MAX_DEG <- 3        # a pair is "non-hub" if BOTH endpoints have <= MAX_DEG trans/cis partners
MIN_MAF <- 0.15     # require appreciable within-hybrid MAF at BOTH loci -- a near-fixed locus
                    # (MAF-LD bound) cannot carry real LD, so its correlation is noise/leverage
N_PLOT  <- 3        # heatmap the strongest ROBUST non-hub pairs

E <- readRDS("data/moduleD_emmax.rds")
clust <- readRDS("data/eMLG_5loci_0025_cM05.rds"); eMLG <- clust$eMLG; groups <- clust$groups
load("data/hybrids_only_maf005.Rdata"); sd <- as.data.table(sample_data)
pop_by_ind <- setNames(sd$Population, sd$Sample_ID)
chr_of <- setNames(as.character(groups$Chr), groups$group_id)
mem_of <- function(g) groups[.(g), on = "group_id", members][[1]]
pop_colors <- c(
  Aland = "#E41A1C", Bunkkeri = "#377EB8", Grundsund = "#4DAF4A", Heinamaki = "#984EA3",
  Hiivola = "#FF7F00", Jarvenpaa = "#FFFF33", Karsikas = "#A65628", Katiskoski = "#F781BF",
  Kummunmaki = "#999999", LangholmenR = "#66C2A5", LangholmenW = "#FC8D62", Nyrhispera74 = "#8DA0CB",
  Nyrhispera75 = "#E78AC3", Parikkala = "#A6D854", Pikkala = "#FFD92F", Sielva = "#E5C494",
  Svanvik1 = "#B3B3B3", Svanvik2 = "#1B9E77", Tvarminne = "#D95F02", Vuosaari = "#7570B3")

## ---- score connectivity and isolate the non-hub trans pairs -------------
h  <- E$hits[paralog == FALSE]
fd <- setNames(h[, .N, by = focal]$N, h[, .N, by = focal]$focal)          # focal out-degree
pin<- setNames(h[, uniqueN(focal), by = partner]$V1, h[, uniqueN(focal), by = partner]$partner)  # partner in-degree
deg_of <- function(g) (if (!is.na(fd[g])) fd[g] else 0L) + (if (!is.na(pin[g])) pin[g] else 0L)

tr <- h[coupling == "repulsion"]
tr[, key := ifelse(focal < partner, paste(focal, partner), paste(partner, focal))]
tr <- tr[order(-F)][!duplicated(key)]
tr[, `:=`(deg_focal = sapply(focal, deg_of), deg_partner = sapply(partner, deg_of))]
tr[, r := mapply(function(a, b) cor(eMLG[, a], eMLG[, b], use = "pairwise.complete.obs"), focal, partner)]
## within-hybrid MAF (near-fixed => noise) and leverage (r after dropping the 5 most
## influential individuals => a robust edge keeps most of its correlation)
maf <- function(g) { p <- mean(eMLG[, g], na.rm = TRUE) / 2; min(p, 1 - p) }
lev <- function(a, b) { x <- eMLG[, a]; y <- eMLG[, b]; ok <- is.finite(x) & is.finite(y)
  infl <- order(-abs((x - mean(x[ok])) * (y - mean(y[ok]))))[1:5]; k <- ok; k[infl] <- FALSE; cor(x[k], y[k]) }
tr[, `:=`(maf_focal = sapply(focal, maf), maf_partner = sapply(partner, maf),
          r_drop5 = mapply(lev, focal, partner))]
nonhub <- tr[deg_focal <= MAX_DEG & deg_partner <= MAX_DEG]
## a genuine pair-specific candidate: appreciable MAF at both loci AND leverage-robust
nonhub[, robust := maf_focal >= MIN_MAF & maf_partner >= MIN_MAF & abs(r_drop5) >= 0.7 * abs(r)]
nonhub <- nonhub[order(-robust, r)]                                       # robust first, then most negative
nonhub[, `:=`(Chr_focal = chr_of[focal], Chr_partner = chr_of[partner])]
## per-cluster DI + parental allele frequencies (aqu-oriented allele freq in each parent),
## for the summary table: DI = max across members (cl); parental freqs = mean over members.
cl <- readRDS("data/moduleC_C3_cl.rds"); snp <- readRDS("data/moduleA_snp.rds")
sm <- snp[groups[, .(marker = unlist(members)), by = group_id], on = "marker"]
cs <- sm[!is.na(f_aqu_parent), .(faqu = mean(f_aqu_parent, na.rm = TRUE), fpol = mean(f_pol_parent, na.rm = TRUE)), by = group_id]
DI <- setNames(cl$DI, cl$group_id); nl <- setNames(groups$n_loci, groups$group_id)
faqu <- setNames(cs$faqu, cs$group_id); fpol <- setNames(cs$fpol, cs$group_id)
nonhub[, `:=`(n_focal = nl[focal], n_partner = nl[partner], DI_focal = DI[focal], DI_partner = DI[partner],
              faqu_focal = faqu[focal], fpol_focal = fpol[focal], faqu_partner = faqu[partner], fpol_partner = fpol[partner],
              disposition = fifelse(robust, "candidate", fifelse(pmin(maf_focal, maf_partner) < MIN_MAF, "near-fixed", "leverage")))]
saveRDS(nonhub, "data/moduleD_nonhub_trans.rds")
message(sprintf("[nonhub] %d non-hub trans pairs of %d clean trans pairs; %d ROBUST (MAF>=%.2f both, leverage-stable)",
                nrow(nonhub), nrow(tr), sum(nonhub$robust), MIN_MAF))
print(nonhub[, .(focal, Chr_focal, partner, Chr_partner, deg = sprintf("%d/%d", deg_focal, deg_partner),
                 F = round(F), r = round(r, 2), r_drop5 = round(r_drop5, 2),
                 maf = sprintf("%.2f/%.2f", maf_focal, maf_partner), robust)])
nonhub <- nonhub[robust == TRUE]                                          # heatmap only genuine candidates

## ---- heatmap the strongest non-hub pairs (trans mirror; see moduleD_trans_heatmaps) --
dir.create("Figures", showWarnings = FALSE)
polarize_to_ref <- function(M, ref) {
  for (j in seq_len(ncol(M))) { cc <- suppressWarnings(cor(M[, j], ref, use = "pairwise.complete.obs"))
    if (is.finite(cc) && cc < 0) M[, j] <- 2 - M[, j] }; M
}
for (i in seq_len(min(N_PLOT, nrow(nonhub)))) {
  A <- nonhub$focal[i]; B <- nonhub$partner[i]
  mkA <- intersect(mem_of(A), colnames(GTs_hybrids_005)); mkB <- intersect(mem_of(B), colnames(GTs_hybrids_005))
  GTa <- polarize_to_ref(GTs_hybrids_005[, mkA, drop = FALSE], eMLG[, A])
  GTb <- polarize_to_ref(GTs_hybrids_005[, mkB, drop = FALSE], eMLG[, B])
  GT <- cbind(GTa, GTb); mk <- c(mkA, mkB)
  lab <- setNames(c(rep(sprintf("%s (%s)", A, chr_of[A]), length(mkA)),
                    rep(sprintf("%s (%s)", B, chr_of[B]), length(mkB))), mk)
  lc <- setNames(c("#1B9E77", "#D95F02"), c(sprintf("%s (%s)", A, chr_of[A]), sprintf("%s (%s)", B, chr_of[B])))
  ord <- rownames(GT)[order(rowMeans(GTa, na.rm = TRUE), na.last = TRUE)]
  plot_genotype_heatmap(
    GT, polarize = FALSE, cluster_columns = FALSE, row_order = ord,
    col_annotation = lab[mk], col_annotation_name = "cluster (chr)", col_annotation_colors = lc,
    row_annotation = pop_by_ind[ord], row_annotation_name = "Population", row_annotation_colors = pop_colors,
    legend_name = "dosage (cluster-oriented)",
    title = sprintf("non-hub trans: %s (%s, deg %d) vs %s (%s, deg %d)   r = %.2f, F = %.0f",
                    A, chr_of[A], nonhub$deg_focal[i], B, chr_of[B], nonhub$deg_partner[i], nonhub$r[i], nonhub$F[i]),
    out_file = sprintf("Figures/moduleD_nonhub_trans_%02d_%s_%s.png", i, A, B), width = 11, height = 7)
  message(sprintf("[%d] %s vs %s | r=%.2f -> Figures/moduleD_nonhub_trans_%02d_%s_%s.png", i, A, B, nonhub$r[i], i, A, B))
}
cat("done\n")
