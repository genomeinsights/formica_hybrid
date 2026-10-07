## =========================================================================
## module_allele_specific_sorting -- EXPLORATORY: size distribution of sorted clusters
##
## A cluster = consecutive directionally sorted units on one chromosome with no gap
## between neighbouring sorted units larger than GAP_KB (single linkage; unsorted
## units may lie in between). For each cluster: number of sorted units, physical
## and genetic span, number of unsorted units inside the span (purity = sorted /
## all units in span), direction (share aquilonia), F_ST.
## Reference distributions, same procedure:
##   null  : sorted labels permuted within chromosome x F_ST decile x recombination
##           tercile (N_PERM permutations), i.e. the clustering expected from where
##           differentiated, low-recombination units lie;
##   topFST: the same number of the most differentiated unsorted units.
## Output: data/explore_sorted_cluster_sizes.rds, Figures/explore_sorted_cluster_sizes.png
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/exploratory/explore_sorted_cluster_sizes.R [N_PERM=200]
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
args <- commandArgs(trailingOnly = TRUE)
N_PERM <- if (length(args) >= 1) as.integer(args[1]) else 200L
GAPS <- c(20, 50, 100)
set.seed(1)
u <- readRDS(file.path(OUT_DATA, "01_units.rds"))[, .(group_id, Chr, ChrNum, Pos, cM_pos, FST, recomb_rate, sort_class, sorted, n_loci)]
setorder(u, ChrNum, Pos); u[, idx := .I]
u[, fst_dec := cut(FST, quantile(FST, 0:10 / 10, na.rm = TRUE), include.lowest = TRUE, labels = FALSE)]
u[, rec_ter := cut(recomb_rate, quantile(recomb_rate, 0:3 / 3, na.rm = TRUE), include.lowest = TRUE, labels = FALSE)]
u[is.na(fst_dec), fst_dec := 0L][is.na(rec_ter), rec_ter := 0L]
g_str <- split(u$idx, paste(u$Chr, u$fst_dec, u$rec_ter))

clusters <- function(s, gap_kb) {
  x <- u[s]
  x[, cl := cumsum(c(TRUE, diff(Pos) > gap_kb * 1e3 | diff(ChrNum) != 0))]
  x[, .(n_sorted = .N, Chr = Chr[1], start = min(Pos), end = max(Pos), span_kb = (max(Pos) - min(Pos)) / 1e3,
        span_cM = diff(range(cM_pos, na.rm = TRUE)), first = min(idx), last = max(idx),
        share_aqu = mean(sort_class == "aquilonia"), FST = mean(FST, na.rm = TRUE), snps = sum(n_loci)), by = cl][
    , n_units_in_span := last - first + 1L][, purity := n_sorted / n_units_in_span][]
}
size_class <- function(n) cut(n, c(0, 1, 2, 3, 5, 10, Inf), labels = c("1", "2", "3", "4-5", "6-10", ">10"))
summ <- function(cl) {
  t <- table(size_class(cl$n_sorted))
  c(n_clusters = nrow(cl), as.list(setNames(as.numeric(t), paste0("n_", names(t)))),
    units_in_multi = sum(cl$n_sorted[cl$n_sorted > 1]) / sum(cl$n_sorted), max_size = max(cl$n_sorted))
}
top <- rep(FALSE, nrow(u)); top[u[sort_class == "unsorted" & !is.na(FST)][order(-FST)][1:sum(u$sorted), idx]] <- TRUE
out <- list(); tabs <- list()
for (G in GAPS) {
  obs <- clusters(u$sorted, G)
  nul <- rbindlist(lapply(seq_len(N_PERM), function(k) {
    s <- u$sorted; for (g in g_str) if (length(g) > 1) s[g] <- s[g][sample.int(length(g))]
    as.data.table(summ(clusters(s, G))) }))
  tp <- clusters(top, G)
  tabs[[as.character(G)]] <- rbind(data.table(gap_kb = G, set = "sorted (observed)", as.data.table(summ(obs))),
    data.table(gap_kb = G, set = "null (stratified permutation), mean", nul[, lapply(.SD, mean)]),
    data.table(gap_kb = G, set = "null, 97.5% quantile", nul[, lapply(.SD, quantile, 0.975)]),
    data.table(gap_kb = G, set = "top-F_ST unsorted", as.data.table(summ(tp))))
  out[[as.character(G)]] <- list(observed = obs, null = nul, topFST = tp)
}
tab <- rbindlist(tabs)
options(width = 220)
cat(sprintf("[sizes] %d sorted units; clusters = sorted units linked by gaps <= G kb\n", sum(u$sorted)))
print(tab, digits = 3)

obs <- out[["50"]]$observed
cat("\n[sizes] multi-unit clusters (gap <= 50 kb): span, purity, direction\n")
print(obs[n_sorted > 1, .(n = .N, median_span_kb = median(span_kb), q90_span_kb = quantile(span_kb, 0.9),
                         median_span_cM = median(span_cM, na.rm = TRUE), median_purity = median(purity),
                         share_same_direction = mean(share_aqu %in% c(0, 1)), median_snps = as.numeric(median(snps))),
          by = .(size = size_class(n_sorted))][order(size)], digits = 3)
cat("\n[sizes] largest clusters (gap <= 50 kb)\n")
print(obs[order(-n_sorted)][1:15, .(Chr, start, end, n_sorted, n_units_in_span, purity, span_kb, span_cM, share_aqu, FST, snps)], digits = 3)
saveRDS(list(table = tab, clusters = out, gaps = GAPS, n_perm = N_PERM), file.path(OUT_DATA, "explore_sorted_cluster_sizes.rds"))

## figure: cluster-size distribution (units per size class), observed vs null vs top-F_ST, per gap
fd <- rbindlist(lapply(as.character(GAPS), function(G) {
  o <- out[[G]]
  share <- function(cl) { s <- size_class(cl$n_sorted); tapply(cl$n_sorted, s, sum) / sum(cl$n_sorted) }
  nl <- rbindlist(lapply(seq_len(min(N_PERM, 100)), function(k) {
    s <- u$sorted; for (g in g_str) if (length(g) > 1) s[g] <- s[g][sample.int(length(g))]
    data.table(size = levels(size_class(1)), v = as.numeric(share(clusters(s, as.numeric(G))))) }))
  rbind(data.table(gap = paste0("gap <= ", G, " kb"), set = "sorted", size = levels(size_class(1)), v = as.numeric(share(o$observed))),
        nl[, .(v = mean(v, na.rm = TRUE)), by = size][, `:=`(gap = paste0("gap <= ", G, " kb"), set = "null (stratified permutation)")],
        data.table(gap = paste0("gap <= ", G, " kb"), set = "top-F_ST unsorted", size = levels(size_class(1)), v = as.numeric(share(o$topFST))))
}))
fd[is.na(v), v := 0][, size := factor(size, levels = levels(size_class(1)))][, gap := factor(gap, levels = paste0("gap <= ", GAPS, " kb"))]
p <- ggplot(fd, aes(size, v, fill = set)) + geom_col(position = position_dodge(0.8), width = 0.75) + facet_wrap(~ gap, nrow = 1) +
  scale_fill_manual(values = c("grey60", "black", "#3182bd"), name = NULL) +
  labs(x = "cluster size (sorted units)", y = "share of units in clusters of this size") + theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG, "explore_sorted_cluster_sizes.png"), p, width = 11, height = 3.8, dpi = 200)
cat("[sizes] done\n")
