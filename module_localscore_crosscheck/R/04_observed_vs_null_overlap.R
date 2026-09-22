## =========================================================
## module_localscore_crosscheck -- observed local-score windows vs. the
## 40 null-replicate windows: physical overlap check.
## =========================================================
## Fresh implementation. Computes compute.local.scores() on the authoritative
## OBSERVED full-SNP BayPass scans in module_manuscript_rho05/baypass_stage1/
## aland_excluded/ (PC1, PC2, bio6, bio11, bio_winter as continuous BF;
## mito_C2 as p-value contrast), using the identical call convention as the
## null sweep (R/03_compute_local_scores_fullsnp.R), then checks whether any
## null-replicate window (from all_significant_windows.rds, already
## committed) physically overlaps an observed window on the same chromosome.
##
## Does NOT reuse the exploratory module's nullcont_draw74/nullmito_draw212
## runs or any of its window tables -- only the 40 freshly generated null
## replicates in this module and the authoritative observed scans.
##
## Reads:
##   module_manuscript_rho05/baypass_stage1/aland_excluded/
##     {PC1,PC2,bio6,bio11,bio_winter}_fullSNP_stage1Omega_summary_betai_reg.out
##     {same}_summary_pi_xtx.out
##     mito_C2_fullSNP_stage1Omega_summary_contrast.out (+ its pi_xtx)
##   module_localscore_crosscheck/data/fullsnp_chr_pos.rds
##   module_localscore_crosscheck/full_snp_null10/all_significant_windows.rds
## Writes:
##   module_localscore_crosscheck/full_snp_null10/observed_windows.rds/.tsv
##   module_localscore_crosscheck/full_snp_null10/observed_vs_null_overlap.tsv
##
## Run from the repo root:
##   Rscript module_localscore_crosscheck/R/04_observed_vs_null_overlap.R
## =========================================================

suppressMessages(library(data.table))

source("~/gitlab/baypass_public-master/utils/baypass_utils.R")

BP_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded"
OUT_DIR <- "module_localscore_crosscheck/full_snp_null10"
pos <- readRDS(file.path("module_localscore_crosscheck/data/fullsnp_chr_pos.rds"))
stopifnot(nrow(pos) == 1114423)

run_bf <- function(tag, bf_file, pi_file) {
  bf <- fread(bf_file, select = c("MRK", "BF(dB)"))
  pi <- fread(pi_file, select = c("MRK", "M_P"))
  stopifnot(all(bf$MRK == seq_len(nrow(bf))), all(pi$MRK == seq_len(nrow(pi))),
            nrow(bf) == nrow(pos), nrow(pi) == nrow(pos))
  res <- compute.local.scores(snp.position = pos, snp.pi = pi$M_P, snp.bf = bf$`BF(dB)`,
                               xi = 1, pval.local.score.thres = 0.01, min.maf = 0.2, min.nsnp = 1e4)
  w <- res$significant.windows
  if (is.null(w) || nrow(w) == 0) return(data.table())
  dt <- as.data.table(w); dt[, covariate := tag]; dt
}

run_pval <- function(tag, pval_file, pi_file) {
  pv <- fread(pval_file, select = c("MRK", "log10(1/pval)"))
  pi <- fread(pi_file, select = c("MRK", "M_P"))
  stopifnot(all(pv$MRK == seq_len(nrow(pv))), all(pi$MRK == seq_len(nrow(pi))),
            nrow(pv) == nrow(pos), nrow(pi) == nrow(pos))
  res <- compute.local.scores(snp.position = pos, snp.pi = pi$M_P, snp.pvalue = pv$`log10(1/pval)`,
                               xi = 1, pval.local.score.thres = 0.01, min.maf = 0.2, min.nsnp = 1e4)
  w <- res$significant.windows
  if (is.null(w) || nrow(w) == 0) return(data.table())
  dt <- as.data.table(w); dt[, covariate := tag]; dt
}

observed <- rbindlist(list(
  run_bf("PC1",        file.path(BP_DIR, "PC1_fullSNP_stage1Omega_summary_betai_reg.out"),
                        file.path(BP_DIR, "PC1_fullSNP_stage1Omega_summary_pi_xtx.out")),
  run_bf("PC2",        file.path(BP_DIR, "PC2_fullSNP_stage1Omega_summary_betai_reg.out"),
                        file.path(BP_DIR, "PC2_fullSNP_stage1Omega_summary_pi_xtx.out")),
  run_bf("bio6",       file.path(BP_DIR, "bio6_fullSNP_stage1Omega_summary_betai_reg.out"),
                        file.path(BP_DIR, "bio6_fullSNP_stage1Omega_summary_pi_xtx.out")),
  run_bf("bio11",      file.path(BP_DIR, "bio11_fullSNP_stage1Omega_summary_betai_reg.out"),
                        file.path(BP_DIR, "bio11_fullSNP_stage1Omega_summary_pi_xtx.out")),
  run_bf("bio_winter", file.path(BP_DIR, "bio_winter_fullSNP_stage1Omega_summary_betai_reg.out"),
                        file.path(BP_DIR, "bio_winter_fullSNP_stage1Omega_summary_pi_xtx.out")),
  run_pval("mitoC2",   file.path(BP_DIR, "mito_C2_fullSNP_stage1Omega_summary_contrast.out"),
                        file.path(BP_DIR, "mito_C2_fullSNP_stage1Omega_summary_pi_xtx.out"))
), fill = TRUE)

cat("[observed] windows per covariate:\n")
print(observed[, .N, by = covariate])
fwrite(observed, file.path(OUT_DIR, "observed_windows.tsv"), sep = "\t")
saveRDS(observed, file.path(OUT_DIR, "observed_windows.rds"))

## ---- overlap check against all 40 null-replicate windows -------------------
null_windows <- readRDS(file.path(OUT_DIR, "all_significant_windows.rds"))

physical_overlap <- function(chr1, beg1, end1, chr2, beg2, end2)
  chr1 == chr2 & beg1 <= end2 & beg2 <= end1

overlap_rows <- list()
for (i in seq_len(nrow(observed))) {
  o <- observed[i]
  hits <- null_windows[physical_overlap(chr, beg, end, o$chr, o$beg, o$end)]
  if (nrow(hits) > 0) {
    overlap_rows[[length(overlap_rows) + 1]] <- cbind(
      observed_covariate = o$covariate, observed_chr = o$chr, observed_beg = o$beg, observed_end = o$end,
      hits[, .(null_mode = mode, null_nulltype = nulltype, null_draw = draw, null_chr = chr, null_beg = beg, null_end = end)]
    )
  }
}
overlap_dt <- if (length(overlap_rows)) rbindlist(overlap_rows) else data.table()

cat("\n[overlap] observed windows with at least one overlapping null-replicate window:", length(overlap_rows), "of", nrow(observed), "\n")
if (nrow(overlap_dt) > 0) print(overlap_dt)
fwrite(overlap_dt, file.path(OUT_DIR, "observed_vs_null_overlap.tsv"), sep = "\t")

message("[overlap] done")
