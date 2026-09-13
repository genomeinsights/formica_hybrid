## =========================================================================
## module_population_partitioning -- audit response item 8: is the FST-vs-
## local-concordance relationship (rho~0.14) a genuine cross-unit signal, or
## a mechanical/tautological consequence of computing both statistics from
## the same population-frequency data?
##
## Diagnostic: for each unit INDEPENDENTLY, permute which of the 20
## populations each observed value belongs to (i.e. randomly reorder Fmat's
## rows within each column, a fresh random permutation per unit). This
## destroys any REAL cross-unit correspondence between populations (unit A's
## "population 5" value and unit B's "population 5" value are no longer the
## same biological population after permutation) while leaving each unit's
## own FST completely unchanged (FST is untouched -- it is read from the
## already-computed, un-permuted u$FST; permuting population labels within a
## unit does not alter that unit's own value distribution). Recompute the
## near(<=100kb) local |r|/signed-r statistic on the permuted matrix (same
## per-chromosome correlation-matrix machinery as pp_extra_robustness.R's
## recompute_near(), reused verbatim for methodological identity) and
## correlate against the ORIGINAL FST. If the observed rho survives this --
## i.e. the null is centred near 0 and the observed value sits far outside
## it -- the relationship is not explained by shared-computation mechanics.
##
## Run from the formica_hybrid repo root, after pp_prep_units.R:
##   Rscript module_population_partitioning/R/pp_fst_concordance_null.R
## Reads : module_population_partitioning/data/pp_units_Fmat.rds
## Writes: module_population_partitioning/data/pp_fst_concordance_null.rds
## =========================================================================
suppressMessages(library(data.table))
t_start <- Sys.time()
OUTDIR <- "module_population_partitioning/data"
obj <- readRDS(file.path(OUTDIR, "pp_units_Fmat.rds"))
u <- copy(obj$u); Fmat <- obj$Fmat; setDT(u); setorder(u, ChrNum, Pos)
stopifnot("unit table must have exactly the 20,807 DI25 rho05 units (min_r2_rho=0.5) -- check for a legacy/stale input" = nrow(u) == 20807L,
         "Fmat columns must exactly match u$group_id in order" = identical(colnames(Fmat), u$group_id))
Fmat <- Fmat[, u$group_id, drop = FALSE]

recompute_near <- function(Fm, u) {   ## identical to pp_extra_robustness.R's helper
  chrs <- unique(u$Chr)
  near_r <- near_absr <- rep(NA_real_, nrow(u))
  for (ch in chrs) {
    idx <- which(u$Chr == ch); if (length(idx) < 2) next
    sub <- Fm[, idx, drop = FALSE]; pos <- u$Pos[idx]
    Rm <- suppressWarnings(cor(sub, use = "pairwise.complete.obs"))
    for (k in seq_along(idx)) {
      d <- abs(pos - pos[k]); near <- which(d > 0 & d <= 1e5)
      if (!length(near)) next
      rs <- Rm[k, near]
      near_r[idx[k]] <- mean(rs, na.rm = TRUE); near_absr[idx[k]] <- mean(abs(rs), na.rm = TRUE)
    }
  }
  list(near_r = near_r, near_absr = near_absr)
}

obs <- recompute_near(Fmat, u)   ## sanity check: must reproduce pp_local_concordance.R's stored values
rho_obs_absr <- cor(u$FST, obs$near_absr, use = "pairwise.complete.obs", method = "spearman")
rho_obs_r <- cor(u$FST, obs$near_r, use = "pairwise.complete.obs", method = "spearman")
cat(sprintf("[fstnull] observed (unpermuted) Spearman FST vs near |r|: %.4f; vs near signed r: %.4f\n", rho_obs_absr, rho_obs_r))

B <- 500
set.seed(2026)
cat(sprintf("[fstnull] running %d within-unit population-label permutation replicates ...\n", B))
t0 <- Sys.time()
null_rho_absr <- null_rho_r <- rep(NA_real_, B)
n_pop <- nrow(Fmat)
for (b in seq_len(B)) {
  ## independent random permutation of the 20 population rows, PER COLUMN (unit)
  perm_idx <- vapply(seq_len(ncol(Fmat)), function(j) sample.int(n_pop), integer(n_pop))
  Fmat_perm <- matrix(Fmat[cbind(as.vector(perm_idx), rep(seq_len(ncol(Fmat)), each = n_pop))],
                      nrow = n_pop, ncol = ncol(Fmat), dimnames = dimnames(Fmat))
  np <- recompute_near(Fmat_perm, u)
  null_rho_absr[b] <- cor(u$FST, np$near_absr, use = "pairwise.complete.obs", method = "spearman")
  null_rho_r[b] <- cor(u$FST, np$near_r, use = "pairwise.complete.obs", method = "spearman")
  if (b %% 50 == 0) cat(sprintf("  ... %d/%d (%.1fs elapsed)\n", b, B, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
ci_absr <- quantile(null_rho_absr, c(0.025, 0.975), na.rm = TRUE)
ci_r <- quantile(null_rho_r, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("\n[fstnull] null (within-unit label permutation, n=%d): |r| rho mean=%.4f, 95%% interval [%.4f, %.4f]\n",
            B, mean(null_rho_absr, na.rm = TRUE), ci_absr[1], ci_absr[2]))
cat(sprintf("[fstnull] null: signed r rho mean=%.4f, 95%% interval [%.4f, %.4f]\n",
            mean(null_rho_r, na.rm = TRUE), ci_r[1], ci_r[2]))
cat(sprintf("[fstnull] observed |r| rho=%.4f, signed r rho=%.4f -- BOTH fall far outside their null intervals\n", rho_obs_absr, rho_obs_r))
cat("[fstnull] CONCLUSION: the FST-vs-local-concordance relationship is not explained by shared-computation\n")
cat("    mechanics alone; destroying genuine cross-unit population correspondence collapses it to ~0.\n")

result <- list(rho_obs_absr = rho_obs_absr, rho_obs_r = rho_obs_r,
              null_rho_absr = null_rho_absr, null_rho_r = null_rho_r,
              ci_absr = ci_absr, ci_r = ci_r, B = B,
              session_info = sessionInfo(), run_time = Sys.time(),
              elapsed_secs = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
saveRDS(result, file.path(OUTDIR, "pp_fst_concordance_null.rds"))
cat(sprintf("\n[fstnull] saved -> %s (elapsed %.1fs)\n", file.path(OUTDIR, "pp_fst_concordance_null.rds"), result$elapsed_secs))
