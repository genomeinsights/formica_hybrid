## =========================================================================
## module_population_partitioning -- 02: local partition concordance between
## neighbouring DI25 LD-reduced units.
##
## For every within-chromosome pair of units: signed Pearson r and |r| between
## their (population x unit) oriented aquilonia-frequency profiles, Euclidean
## distance between standardized profiles, and physical distance (bp).
## Then: adjacent-unit pairs, distance-binned similarity, and FST vs local
## similarity (adjacent-neighbour and <=100kb window), stratified by FST
## quartile, sort_class, and current_map_DI decile.
##
## Run from the formica_hybrid repo root, after pp_prep_units.R:
##   Rscript module_population_partitioning/R/pp_local_concordance.R
## Reads : module_population_partitioning/data/pp_units_Fmat.rds
## Writes: module_population_partitioning/data/pp_concordance_results.rds
##           list(u = <unit table + nxt_r/nxt_absr/near_r/near_absr/FST_q/DI_decile>,
##                pairs_summary = <similarity by physical-distance bin>)
##         module_population_partitioning/data/pp_all_pairs.csv.gz
##           (~2.7M within-chromosome unit pairs: Chr, i, j, dist_bp, r, eucl, absr)
## =========================================================================
suppressMessages({ library(data.table); library(ggplot2) })
OUTDIR <- "module_population_partitioning/data"
obj <- readRDS(file.path(OUTDIR, "pp_units_Fmat.rds"))
u <- obj$u; Fmat <- obj$Fmat; hybrid_pops <- obj$hybrid_pops
setDT(u); setorder(u, ChrNum, Pos)
stopifnot(identical(colnames(Fmat), u$group_id))

## ---------------------------------------------------------------------
## A. per-chromosome ALL-PAIRS: signed r, |r|, Euclidean on standardized
##    profiles, physical distance (bp). Accumulated as flat vectors
##    (not a big object) -- ~2-3M pairs genome-wide.
## ---------------------------------------------------------------------
chrs <- unique(u$Chr)
acc <- vector("list", length(chrs))
Rm_list <- vector("list", length(chrs)); names(Rm_list) <- chrs
idx_list <- vector("list", length(chrs)); names(idx_list) <- chrs
for (k in seq_along(chrs)) {
  ch <- chrs[k]
  idx <- which(u$Chr == ch)                      # already position-ordered
  idx_list[[ch]] <- idx
  if (length(idx) < 2) { acc[[k]] <- NULL; next }
  sub <- Fmat[, idx, drop = FALSE]                # pops x n_units_chr
  pos <- u$Pos[idx]
  Rm  <- suppressWarnings(cor(sub, use = "pairwise.complete.obs"))
  Rm_list[[ch]] <- Rm
  Zm  <- scale(sub, center = TRUE, scale = TRUE)
  Zm[is.na(Zm)] <- 0                              # for Euclidean only: missing pop contributes 0 deviation
  n <- ncol(sub)
  ii <- rep(seq_len(n), times = n); jj <- rep(seq_len(n), each = n)
  keep <- ii < jj
  ii <- ii[keep]; jj <- jj[keep]
  ## BUGFIX (caught by AUDIT.md): this must be rowSums (one value per PAIR,
  ## i.e. per row of the ii/jj-indexed difference matrix), not colSums (which
  ## silently returned one value per POPULATION, length 20, recycled into the
  ## per-pair column -- the source of the "recycling" warnings previously
  ## seen and of the previous, meaningless ~100-scale eucl means). Note: for
  ## standardized profiles with negligible missingness (this dataset),
  ## eucl = sqrt(2*(n_pop-1)*(1-r)) is a deterministic function of the signed
  ## correlation r, so it carries no information beyond r/absr -- kept only
  ## for completeness against the original brief, not used in any conclusion.
  eucl <- sqrt(rowSums((t(Zm)[ii, , drop = FALSE] - t(Zm)[jj, , drop = FALSE])^2)) / sqrt(nrow(Zm))
  acc[[k]] <- data.table(Chr = ch, i = idx[ii], j = idx[jj],
                         dist_bp = abs(pos[ii] - pos[jj]),
                         r = Rm[cbind(ii, jj)], eucl = eucl)
}
pairs <- rbindlist(acc)
pairs[, absr := abs(r)]
cat(sprintf("[concordance] %d within-chromosome unit pairs (all physical distances)\n", nrow(pairs)))
cat(sprintf("[concordance] genome-wide: mean r = %.3f, mean |r| = %.3f, median dist = %.0f bp\n",
            mean(pairs$r, na.rm = TRUE), mean(pairs$absr, na.rm = TRUE), median(pairs$dist_bp)))

## ---------------------------------------------------------------------
## B. adjacent-unit pairs only (immediate physical neighbour, same chr)
## ---------------------------------------------------------------------
u[, nxt_r := NA_real_][, nxt_absr := NA_real_][, nxt_eucl := NA_real_][, nxt_dist := NA_real_]
for (ch in chrs) {
  idx <- which(u$Chr == ch)
  if (length(idx) < 2) next
  for (k in seq_len(length(idx) - 1)) {
    i <- idx[k]; j <- idx[k + 1]
    r <- suppressWarnings(cor(Fmat[, i], Fmat[, j], use = "pairwise.complete.obs"))
    zi <- scale(Fmat[, i]); zj <- scale(Fmat[, j])
    zi[is.na(zi)] <- 0; zj[is.na(zj)] <- 0
    e <- sqrt(sum((zi - zj)^2)) / sqrt(length(zi))
    u$nxt_r[i] <- r; u$nxt_absr[i] <- abs(r); u$nxt_eucl[i] <- e
    u$nxt_dist[i] <- u$Pos[j] - u$Pos[i]
  }
}
cat(sprintf("\n[concordance] adjacent-unit pairs (n=%d): mean r = %.3f, mean |r| = %.3f\n",
            sum(!is.na(u$nxt_r)), mean(u$nxt_r, na.rm = TRUE), mean(u$nxt_absr, na.rm = TRUE)))
cat("[concordance] adjacent |r| by sort_class (this unit):\n")
print(u[!is.na(nxt_r), .(n = .N, mean_absr = round(mean(nxt_absr), 3),
                        mean_r = round(mean(nxt_r), 3)), by = sort_class][order(-n)])

## ---------------------------------------------------------------------
## C. distance-binned local similarity (all-pairs, pooled genome-wide)
## ---------------------------------------------------------------------
BRK <- c(0, 5e3, 2e4, 1e5, 5e5, 2e6, 1e7, Inf)
LAB <- c("0-5kb","5-20kb","20-100kb","100-500kb","0.5-2Mb","2-10Mb",">10Mb")
pairs[, dbin := cut(dist_bp, BRK, labels = LAB)]
dbin_summary <- pairs[, .(n = .N, mean_r = mean(r, na.rm = TRUE), se_r = sd(r, na.rm = TRUE) / sqrt(.N),
                          mean_absr = mean(absr, na.rm = TRUE), se_absr = sd(absr, na.rm = TRUE) / sqrt(.N),
                          mean_eucl = mean(eucl, na.rm = TRUE)), by = dbin][order(dbin)]
cat("\n[concordance] similarity vs physical distance bin (all within-chromosome pairs):\n"); print(dbin_summary)

## ---------------------------------------------------------------------
## D. FST vs local similarity
## ---------------------------------------------------------------------
u[, near_absr := NA_real_][, near_r := NA_real_][, n_near := 0L]
NEAR_BP <- 1e5   # <=100kb window, excluding self -- reuses the Rm matrices from section A (no recomputation)
for (ch in chrs) {
  idx <- idx_list[[ch]]; if (length(idx) < 2) next
  pos <- u$Pos[idx]; Rm <- Rm_list[[ch]]
  for (k in seq_along(idx)) {
    d <- abs(pos - pos[k]); near <- which(d > 0 & d <= NEAR_BP)
    if (!length(near)) next
    rs <- Rm[k, near]
    u$near_r[idx[k]]    <- mean(rs, na.rm = TRUE)
    u$near_absr[idx[k]] <- mean(abs(rs), na.rm = TRUE)
    u$n_near[idx[k]]    <- length(near)
  }
}
cat(sprintf("\n[concordance] units with >=1 neighbour within %.0f kb: %d / %d\n", NEAR_BP/1e3, sum(u$n_near > 0), nrow(u)))
cat(sprintf("[concordance] Spearman FST vs adjacent |r|: rho = %.3f (n=%d)\n",
            cor(u$FST, u$nxt_absr, use = "pairwise.complete.obs", method = "spearman"),
            sum(!is.na(u$FST) & !is.na(u$nxt_absr))))
cat(sprintf("[concordance] Spearman FST vs near(<=100kb) |r|: rho = %.3f (n=%d)\n",
            cor(u$FST, u$near_absr, use = "pairwise.complete.obs", method = "spearman"),
            sum(!is.na(u$FST) & !is.na(u$near_absr))))

u[, FST_q := cut(FST, quantile(FST, seq(0, 1, 0.25), na.rm = TRUE), include.lowest = TRUE,
                 labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"))]
cat("\n[concordance] local similarity by FST quartile:\n")
print(u[!is.na(FST_q), .(n = .N, mean_nxt_absr = round(mean(nxt_absr, na.rm = TRUE), 3),
                         mean_near_absr = round(mean(near_absr, na.rm = TRUE), 3)), by = FST_q][order(FST_q)])

cat("\n[concordance] local similarity by sort_class x FST(sorted units only):\n")
print(u[sort_class != "unsorted", .(n = .N, median_FST = round(median(FST, na.rm = TRUE), 3),
                                    mean_near_absr = round(mean(near_absr, na.rm = TRUE), 3)), by = sort_class])

u[!is.na(current_map_DI), DI_decile := cut(current_map_DI, quantile(current_map_DI, seq(0, 1, 0.1), na.rm = TRUE),
                                           include.lowest = TRUE, labels = FALSE)]
cat("\n[concordance] local similarity by current_map_DI decile (1=most negative/least diagnostic):\n")
print(u[!is.na(DI_decile), .(n = .N, mean_FST = round(mean(FST, na.rm = TRUE), 3),
                             mean_near_absr = round(mean(near_absr, na.rm = TRUE), 3)), by = DI_decile][order(DI_decile)])

saveRDS(list(u = u, pairs_summary = dbin_summary), file.path(OUTDIR, "pp_concordance_results.rds"))
fwrite(pairs, file.path(OUTDIR, "pp_all_pairs.csv.gz"))
cat("\n[concordance] saved -> pp_concordance_results.rds, pp_all_pairs.csv.gz\n")
