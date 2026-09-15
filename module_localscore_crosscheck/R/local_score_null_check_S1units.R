## =========================================================
## module_manuscript_rho05 -- does the local-score method's analytic
## threshold survive the structured-Omega null? (Stage-1-cluster resolution)
## =========================================================
## moduleB_stage1_local_score_regions_S1units.R found 1/0/1/8 "significant"
## local-score windows for PC1/PC2/bio_winter/mitoC2 using BayPass's own
## ANALYTIC per-chromosome threshold (Fariello et al. 2017 / Bonhomme et
## al. 2019) -- not the 10,000-draw Omega-structured null used everywhere
## else in this pipeline. This script asks the obvious calibration
## question directly: if a NULL covariate/contrast (no causal link to
## climate or mitotype, only shares the populations' Omega covariance
## structure) is run through the SAME local-score procedure, how many
## "significant" windows does it produce purely by chance?
##
## Uses the SAME persisted null BF/log10p matrices already built for the
## Module C genome-wide calibration:
##  - continuous-covariate null (shared by PC1/PC2/bio_winter): pulled back
##    from mini2 (baypass_stage1_S1units/null/bf_matrices/cRegen_bf_b##.rds,
##    5 batches = 1000 draws pulled for this check).
##  - mitoC2's own dedicated contrast-mode null (already local from the
##    earlier mitoC2 Module C enrichment work, mitoC2_bf_b##.rds) -- same
##    1000-draw subset (first 5 batches) for a fair, equal-N comparison.
## Both were built at Stage-1-cluster resolution (18,361 units) -- there is
## no full-SNP-resolution null (would require ~9-10h of additional BayPass
## compute per batch), so this check is necessarily at the coarser
## resolution, matching moduleB_stage1_local_score_regions_S1units.R.
##
## Reports the distribution of "significant window" counts across 1000
## null draws per covariate-type, and saves the single WORST-CASE null
## draw (most windows) for each, so it can be visualised the same way as
## the real covariates (moduleB_stage1_local_score_manhattan_null_S1units.R).
##
## Reads : module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/
##           null/bf_matrices/cRegen_bf_b01..05.rds (continuous null)
##           null/bf_matrices/mitoC2_bf_b01..05.rds (contrast null)
##           PC1_S1units_..._pi_xtx.out, mito_C2_S1units_..._pi_xtx.out
## Writes: module_manuscript_rho05/data/moduleB_stage1_localscore_nullcheck_S1units.rds
##
## Run from the repo root:
##   Rscript module_manuscript_rho05/R/moduleB_stage1_local_score_null_check_S1units.R
## =========================================================

suppressMessages(library(data.table))
source("~/gitlab/baypass_public-master/utils/baypass_utils.R")

UNIT_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded_S1units"
BFDIR    <- file.path(UNIT_DIR, "null", "bf_matrices")
DATA     <- "module_manuscript_rho05/data"
dir.create(DATA, showWarnings = FALSE, recursive = TRUE)

stage1 <- readRDS("module0_ld_pruning/data/pruned_stage1.rds")
cl <- as.data.table(stage1$clusters)
cl5 <- cl[n_snps >= 5]
group_order <- readLines(file.path(UNIT_DIR, "S1units_group_order.txt"))
stopifnot(identical(group_order, paste0("S1_", cl5$CL_id)))

base_pos <- data.table(row = seq_len(nrow(cl5)), Chr = cl5$Chr,
                       Pos = as.integer(sub(".*:", "", cl5$core_snp)))
setorder(base_pos, Chr, Pos)
ord <- base_pos$row
pos_df <- as.data.frame(base_pos[, .(Chr, Pos)])
MIN_NSNP <- 100L

pi_pc1 <- fread(file.path(UNIT_DIR, "PC1_S1units_withOmega_summary_pi_xtx.out"), select = c("MRK", "M_P"))$M_P
pi_c2  <- fread(file.path(UNIT_DIR, "mito_C2_S1units_summary_pi_xtx.out"), select = c("MRK", "M_P"))$M_P

NBATCH_SCAN <- 5L   # 5 x 200 = 1000 null draws per covariate-type
BATCH <- 200L

scan_nulls <- function(label, file_pattern, pi_vec, is_bf) {
  message("\n=== scanning ", label, " null draws (", NBATCH_SCAN * BATCH, " total) ===")
  n_win <- integer(NBATCH_SCAN * BATCH)
  best_win <- NULL; best_n <- -1L; best_idx <- NA_integer_
  k <- 0L
  for (b in seq_len(NBATCH_SCAN)) {
    f <- sprintf(file_pattern, b)
    M <- readRDS(f)
    stopifnot(nrow(M) == nrow(base_pos), ncol(M) == BATCH)
    t0 <- Sys.time()
    for (j in seq_len(BATCH)) {
      k <- k + 1L
      out <- if (is_bf) {
        invisible(capture.output(res <- compute.local.scores(
          snp.position = pos_df, snp.pi = pi_vec[ord], snp.bf = M[ord, j],
          xi = 1, pval.local.score.thres = 0.01, min.maf = 0.2, min.nsnp = MIN_NSNP)))
      } else {
        invisible(capture.output(res <- compute.local.scores(
          snp.position = pos_df, snp.pi = pi_vec[ord], snp.pvalue = M[ord, j],
          xi = 1, pval.local.score.thres = 0.01, min.maf = 0.2, min.nsnp = MIN_NSNP)))
      }
      nw <- if (is.null(res$significant.windows)) 0L else nrow(res$significant.windows)
      n_win[k] <- nw
      if (nw > best_n) { best_n <- nw; best_win <- res$significant.windows; best_idx <- k }
    }
    message(sprintf("  batch %d/%d done in %.1f s (cumulative max so far: %d windows)",
                    b, NBATCH_SCAN, as.numeric(difftime(Sys.time(), t0, units = "secs")), best_n))
    rm(M); invisible(gc())
  }
  list(n_win = n_win, best_win = best_win, best_idx = best_idx, best_n = best_n)
}

cont_null <- scan_nulls("continuous-covariate (PC1/PC2/bio_winter)",
                        file.path(BFDIR, "cRegen_bf_b%02d.rds"), pi_pc1, is_bf = TRUE)
mito_null <- scan_nulls("mitoC2 contrast", file.path(BFDIR, "mitoC2_bf_b%02d.rds"), pi_c2, is_bf = FALSE)

cat("\n=== continuous-covariate null: window-count distribution (n=", length(cont_null$n_win), ") ===\n")
print(table(cont_null$n_win))
cat("mean:", mean(cont_null$n_win), " median:", median(cont_null$n_win),
    " max:", cont_null$best_n, " (draw #", cont_null$best_idx, ")\n")
cat("fraction of null draws with >=1 window:", mean(cont_null$n_win >= 1), "\n")
cat("fraction of null draws with >=1 window (PC1's real count):", mean(cont_null$n_win >= 1), "\n")

cat("\n=== mitoC2 null: window-count distribution (n=", length(mito_null$n_win), ") ===\n")
print(table(mito_null$n_win))
cat("mean:", mean(mito_null$n_win), " median:", median(mito_null$n_win),
    " max:", mito_null$best_n, " (draw #", mito_null$best_idx, ")\n")
cat("fraction of null draws with >=8 windows (mitoC2's real count):", mean(mito_null$n_win >= 8), "\n")

saveRDS(list(cont_null = cont_null, mito_null = mito_null,
            real_counts = c(PC1 = 1L, PC2 = 0L, bio_winter = 1L, mitoC2 = 8L)),
       file.path(DATA, "moduleB_stage1_localscore_nullcheck_S1units.rds"))
cat("\n[null-check] wrote module_manuscript_rho05/data/moduleB_stage1_localscore_nullcheck_S1units.rds\n")
