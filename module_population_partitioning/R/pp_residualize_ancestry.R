## =========================================================================
## module_population_partitioning -- 08: residualize each unit's population
## profile against genome-wide ancestry (AUDIT.md "most informative next
## analysis").
##
## The raw oriented-frequency profiles mix two things: (1) genome-wide
## differences in ancestry among populations (each population's overall
## hybrid index), and (2) locus-specific departures from that background.
## Local concordance and the FST relationship computed on RAW profiles
## (pp_local_concordance.R) cannot distinguish:
##   (i)  high-FST units repeatedly just reflecting the same genome-wide
##        ancestry gradient (populations that are more aquilonia-admixed
##        overall tend to be more aquilonia-fixed at MOST high-FST loci), vs
##  (ii) high-FST units being driven by different populations departing from
##        their OWN genome-wide ancestry at different loci (the proposed
##        population-partitioning mechanism).
## Only (ii) is genuine locus-specific multilocus partitioning.
##
## For each population p and unit u on chromosome c(u):
##   H_p^{-c(u)} = leave-one-chromosome-out genome-wide ancestry, i.e. the
##     mean oriented-aquilonia frequency of population p across all 11,052
##     DI25 units NOT on chromosome c(u) (excludes the focal unit's own
##     chromosome so within-chromosome LD cannot leak into its own covariate).
## Per unit: f_{u,p} = alpha_u + beta_u * H_p^{-c(u)} + residual_{u,p}
##   (simple OLS across the 20 hybrid populations, one fit per unit).
## Then: R^2_u (fraction of the unit's among-population variance explained by
## genome-wide ancestry) vs FST_u; and the FULL local-concordance analysis
## (adjacent / <=100kb / distance-binned / FST-decile) is re-run on the
## RESIDUAL profiles, mirroring pp_local_concordance.R exactly, for direct
## comparison against the raw-profile results.
##
## Run from the formica_hybrid repo root, after pp_prep_units.R:
##   Rscript module_population_partitioning/R/pp_residualize_ancestry.R
## Reads : module_population_partitioning/data/pp_units_Fmat.rds
##         module_population_partitioning/data/pp_concordance_results.rds  (for FST, sort_class)
## Writes: module_population_partitioning/data/pp_residual_ancestry.rds
## =========================================================================
suppressMessages(library(data.table))
OUTDIR <- "module_population_partitioning/data"
obj <- readRDS(file.path(OUTDIR, "pp_units_Fmat.rds"))
res <- readRDS(file.path(OUTDIR, "pp_concordance_results.rds"))
u <- res$u; Fmat <- obj$Fmat; hybrid_pops <- obj$hybrid_pops
setDT(u); setorder(u, ChrNum, Pos)
stopifnot(identical(colnames(Fmat), u$group_id))

## ---------------------------------------------------------------------
## 1. leave-one-chromosome-out genome-wide ancestry H_p^{-c}, per population
## ---------------------------------------------------------------------
chrs <- unique(u$Chr)
S_p <- rowSums(Fmat, na.rm = TRUE); N_p <- rowSums(!is.na(Fmat))         # totals, all units, all pops
H_loco <- matrix(NA_real_, nrow(Fmat), length(chrs), dimnames = list(rownames(Fmat), chrs))
for (ch in chrs) {
  idx <- which(u$Chr == ch)
  sub <- Fmat[, idx, drop = FALSE]
  sum_c <- rowSums(sub, na.rm = TRUE); n_c <- rowSums(!is.na(sub))
  H_loco[, ch] <- (S_p - sum_c) / (N_p - n_c)
}
cat("[residualize] leave-one-chromosome-out genome-wide ancestry (H_p), by population:\n")
print(round(sort(rowMeans(H_loco), decreasing = TRUE), 3))

## ---------------------------------------------------------------------
## 2. per-unit regression f_{u,p} ~ H_p^{-c(u)}; keep alpha, beta, R^2, residuals
## ---------------------------------------------------------------------
n_units <- ncol(Fmat); n_pops <- nrow(Fmat)
Resid <- matrix(NA_real_, n_pops, n_units, dimnames = dimnames(Fmat))
alpha_u <- beta_u <- R2_u <- rep(NA_real_, n_units); names(alpha_u) <- names(beta_u) <- names(R2_u) <- colnames(Fmat)

chr_of_unit <- u$Chr[match(colnames(Fmat), u$group_id)]
for (ch in chrs) {
  idx <- which(chr_of_unit == ch)
  Hc <- H_loco[, ch]                                     # population-named LOCO ancestry for units on this chromosome
  for (k in idx) {
    y <- Fmat[, k]; ok <- !is.na(y) & !is.na(Hc)
    if (sum(ok) < 4) next
    fit <- lm(y[ok] ~ Hc[ok])
    alpha_u[k] <- coef(fit)[1]; beta_u[k] <- coef(fit)[2]
    R2_u[k] <- summary(fit)$r.squared
    Resid[ok, k] <- resid(fit); Resid[!ok, k] <- NA_real_
  }
}
u[, R2_ancestry := R2_u[group_id]]
u[, beta_ancestry := beta_u[group_id]]
cat(sprintf("\n[residualize] per-unit R^2 (variance explained by genome-wide ancestry): median %.3f, mean %.3f\n",
            median(u$R2_ancestry, na.rm = TRUE), mean(u$R2_ancestry, na.rm = TRUE)))
cat(sprintf("[residualize] Spearman FST vs R2_ancestry: rho = %.3f (n=%d)\n",
            cor(u$FST, u$R2_ancestry, use = "pairwise.complete.obs", method = "spearman"),
            sum(!is.na(u$FST) & !is.na(u$R2_ancestry))))
cat("[residualize] R2_ancestry by sort_class:\n")
print(u[, .(n = .N, median_R2 = round(median(R2_ancestry, na.rm = TRUE), 3)), by = sort_class][order(-n)])

## ---------------------------------------------------------------------
## 3. re-run the FULL local-concordance analysis on the RESIDUAL profiles
##    (mirrors pp_local_concordance.R sections A/C/D exactly, on Resid
##    instead of Fmat)
## ---------------------------------------------------------------------
acc <- vector("list", length(chrs)); Rm_list <- vector("list", length(chrs)); names(Rm_list) <- chrs
idx_list <- vector("list", length(chrs)); names(idx_list) <- chrs
for (k in seq_along(chrs)) {
  ch <- chrs[k]
  idx <- which(u$Chr == ch); idx_list[[ch]] <- idx
  if (length(idx) < 2) { acc[[k]] <- NULL; next }
  sub <- Resid[, idx, drop = FALSE]
  pos <- u$Pos[idx]
  Rm <- suppressWarnings(cor(sub, use = "pairwise.complete.obs"))
  Rm_list[[ch]] <- Rm
  n <- ncol(sub); ii <- rep(seq_len(n), times = n); jj <- rep(seq_len(n), each = n)
  keep <- ii < jj; ii <- ii[keep]; jj <- jj[keep]
  acc[[k]] <- data.table(Chr = ch, i = idx[ii], j = idx[jj], dist_bp = abs(pos[ii] - pos[jj]), r = Rm[cbind(ii, jj)])
}
pairs_resid <- rbindlist(acc); pairs_resid[, absr := abs(r)]
cat(sprintf("\n[residualize] residual-profile all-pairs: mean r = %.3f, mean |r| = %.3f (raw-profile values were 0.031 / 0.193)\n",
            mean(pairs_resid$r, na.rm = TRUE), mean(pairs_resid$absr, na.rm = TRUE)))

BRK <- c(0, 5e3, 2e4, 1e5, 5e5, 2e6, 1e7, Inf)
LAB <- c("0-5kb","5-20kb","20-100kb","100-500kb","0.5-2Mb","2-10Mb",">10Mb")
pairs_resid[, dbin := cut(dist_bp, BRK, labels = LAB)]
dbin_resid <- pairs_resid[!is.na(dbin), .(n = .N, mean_r = mean(r, na.rm = TRUE), mean_absr = mean(absr, na.rm = TRUE)), by = dbin][order(dbin)]
cat("\n[residualize] residual-profile similarity vs physical-distance bin:\n"); print(dbin_resid)

u[, near_r_resid := NA_real_][, near_absr_resid := NA_real_]
NEAR_BP <- 1e5
for (ch in chrs) {
  idx <- idx_list[[ch]]; if (length(idx) < 2) next
  pos <- u$Pos[idx]; Rm <- Rm_list[[ch]]
  for (k in seq_along(idx)) {
    d <- abs(pos - pos[k]); near <- which(d > 0 & d <= NEAR_BP)
    if (!length(near)) next
    rs <- Rm[k, near]
    u$near_r_resid[idx[k]] <- mean(rs, na.rm = TRUE)
    u$near_absr_resid[idx[k]] <- mean(abs(rs), na.rm = TRUE)
  }
}
cat(sprintf("\n[residualize] Spearman FST vs residual near(<=100kb) |r|: rho = %.3f (n=%d)  [raw-profile value was 0.081]\n",
            cor(u$FST, u$near_absr_resid, use = "pairwise.complete.obs", method = "spearman"),
            sum(!is.na(u$FST) & !is.na(u$near_absr_resid))))
cat(sprintf("[residualize] Spearman FST vs residual near(<=100kb) signed r: rho = %.3f (n=%d)  [raw-profile value was 0.114]\n",
            cor(u$FST, u$near_r_resid, use = "pairwise.complete.obs", method = "spearman"),
            sum(!is.na(u$FST) & !is.na(u$near_r_resid))))

u[, FST_dec := cut(FST, quantile(FST, seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)]
cat("\n[residualize] residual local similarity by FST decile:\n")
print(u[!is.na(FST_dec), .(n = .N, mean_near_r_resid = round(mean(near_r_resid, na.rm = TRUE), 3),
                          mean_near_absr_resid = round(mean(near_absr_resid, na.rm = TRUE), 3)), by = FST_dec][order(FST_dec)])

saveRDS(list(H_loco = H_loco, alpha_u = alpha_u, beta_u = beta_u, R2_u = R2_u, Resid = Resid,
            u = u, pairs_resid_summary = dbin_resid),
        file.path(OUTDIR, "pp_residual_ancestry.rds"))
cat("\n[residualize] saved -> pp_residual_ancestry.rds\n")
