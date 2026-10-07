## =========================================================================
## module_allele_specific_sorting -- EXPLORATORY: do directionally sorted units
## cluster along the genome, and do nearby sorted units sort in the same populations?
##
## (1) Clustering. For each distance window (cM, and kb for map cold spots), count pairs
##     of sorted units on the same chromosome. Null: permute the sorted labels
##       (a) within chromosome;
##       (b) within chromosome x F_ST decile x recombination tercile, so that clustering
##           is not merely "sorted units sit in differentiated, low-recombination regions".
##     Reported as observed / mean null pairs, with permutation p (one-sided, excess).
## (2) Concordance. For sorted pairs by window: share sorted in the same direction
##     (aquilonia vs polyctena), and the correlation of their population ancestry
##     profiles (oriented allele frequency across the 20 populations, Fmat), against
##     sorted pairs on different chromosomes.
## Output: data/explore_sorted_clustering.rds, printed tables.
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/exploratory/explore_sorted_clustering.R [N_PERM=1000]
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
args <- commandArgs(trailingOnly = TRUE)
N_PERM <- if (length(args) >= 1) as.integer(args[1]) else 1000L
set.seed(1)
u <- readRDS(file.path(OUT_DATA, "01_units.rds"))[, .(group_id, Chr, ChrNum, Pos, cM_pos, FST, recomb_rate, sort_class, sorted)]
U <- load_units(); Fm <- U$Fmat[, u$group_id]
u[, idx := .I]
cat(sprintf("[clust] %d units, %d sorted (%d aquilonia, %d polyctena)\n", nrow(u), sum(u$sorted),
            sum(u$sort_class == "aquilonia"), sum(u$sort_class == "polyctena")))

## strata for null (b): F_ST decile and recombination tercile (genome-wide), within chromosome
u[, fst_dec := cut(FST, quantile(FST, 0:10 / 10, na.rm = TRUE), include.lowest = TRUE, labels = FALSE)]
u[, rec_ter := cut(recomb_rate, quantile(recomb_rate, 0:3 / 3, na.rm = TRUE), include.lowest = TRUE, labels = FALSE)]
u[is.na(fst_dec), fst_dec := 0L][is.na(rec_ter), rec_ter := 0L]
u[, stratum := paste(Chr, fst_dec, rec_ter)]

## same-chromosome pairs within 5 Mb, with cM and kb windows
CM_W <- c(0, 0.01, 0.2, 1, 5); KB_W <- c(0, 10, 100, 1000, 5000)
pr <- u[, {
  n <- .N; i <- rep(seq_len(n), each = n); j <- rep(seq_len(n), n); k <- i < j & abs(Pos[j] - Pos[i]) <= 5e6
  .(a = idx[i[k]], b = idx[j[k]], kb = abs(Pos[j[k]] - Pos[i[k]]) / 1e3, cm = abs(cM_pos[j[k]] - cM_pos[i[k]]))
}, by = Chr]
pr[, cm_w := cut(cm, CM_W, right = FALSE, labels = FALSE)][, kb_w := cut(kb, KB_W, right = FALSE, labels = FALSE)]
cat(sprintf("[clust] %d same-chromosome pairs within 5 Mb\n", nrow(pr)))
lab_cm <- c("<0.01 cM", "0.01-0.2 cM", "0.2-1 cM", "1-5 cM"); lab_kb <- c("<10 kb", "10-100 kb", "0.1-1 Mb", "1-5 Mb")

count_pairs <- function(s) {
  both <- s[pr$a] & s[pr$b]
  c(tabulate(pr$cm_w[both & !is.na(pr$cm_w)], 4), tabulate(pr$kb_w[both], 4))
}
obs <- count_pairs(u$sorted)
perm_within <- function(groups) {
  s <- u$sorted
  for (g in groups) if (length(g) > 1) s[g] <- s[g][sample.int(length(g))]
  s
}
g_chr <- split(u$idx, u$Chr); g_str <- split(u$idx, u$stratum)
null_a <- replicate(N_PERM, count_pairs(perm_within(g_chr)))
null_b <- replicate(N_PERM, count_pairs(perm_within(g_str)))
clus <- data.table(window = c(lab_cm, lab_kb), scale = rep(c("cM", "kb"), each = 4), observed = obs,
                   exp_chr = rowMeans(null_a), ratio_chr = obs / rowMeans(null_a), p_chr = (1 + rowSums(null_a >= obs)) / (N_PERM + 1),
                   exp_strat = rowMeans(null_b), ratio_strat = obs / rowMeans(null_b), p_strat = (1 + rowSums(null_b >= obs)) / (N_PERM + 1))
cat("\n[clust] (1) sorted-sorted pairs: observed vs permutation nulls (a: within chromosome; b: chromosome x F_ST decile x recomb tercile)\n")
print(clus, digits = 3)

## (2) concordance of nearby sorted pairs
sp <- pr[u$sorted[a] & u$sorted[b]]
dirc <- function(a, b) mean(u$sort_class[a] == u$sort_class[b])
pc <- function(a, b) mapply(function(x, y) suppressWarnings(cor(Fm[, x], Fm[, y])), a, b)
sp[, `:=`(same_dir = u$sort_class[a] == u$sort_class[b], prof_cor = pc(a, b))]
conc <- rbind(sp[!is.na(cm_w), .(n_pairs = .N, same_direction = mean(same_dir), profile_cor = mean(prof_cor, na.rm = TRUE)), by = .(window = lab_cm[cm_w])],
              sp[, .(n_pairs = .N, same_direction = mean(same_dir), profile_cor = mean(prof_cor, na.rm = TRUE)), by = .(window = lab_kb[kb_w])])
## reference: sorted pairs on different chromosomes (random sample)
si <- which(u$sorted); ra <- sample(si, 2e5, TRUE); rb <- sample(si, 2e5, TRUE); k <- u$Chr[ra] != u$Chr[rb]
conc <- rbind(conc, data.table(window = "different chromosomes", n_pairs = sum(k), same_direction = dirc(ra[k], rb[k]),
                               profile_cor = mean(pc(ra[k][1:20000], rb[k][1:20000]), na.rm = TRUE)))
## same, for unsorted pairs matched on nothing (context)
cat("\n[clust] (2) concordance of sorted-sorted pairs by distance\n"); print(conc, digits = 3)
saveRDS(list(clustering = clus, concordance = conc, n_perm = N_PERM), file.path(OUT_DATA, "explore_sorted_clustering.rds"))
cat("[clust] done\n")
