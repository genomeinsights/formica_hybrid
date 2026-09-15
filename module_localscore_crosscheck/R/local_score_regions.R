## =========================================================
## module_localscore_crosscheck -- BayPass's own local-score outlier-region
## method (Fariello et al. 2017; Bonhomme et al. 2019), applied to the
## full-SNP Stage-1-direct scan for all four covariates
## =========================================================
## BayPass ships compute.local.scores() (baypass_public/utils/baypass_utils.R)
## as its recommended way to turn per-SNP evidence into genomic OUTLIER
## REGIONS: each SNP gets a score S_i = stat_i - xi (stat = BF(dB)/10, or
## -log10(p) for a p-value-based statistic; xi is a fixed per-step penalty,
## default 1); a per-chromosome Lindley process L_i = max(0, L_{i-1}+S_i) is
## then accumulated along the chromosome (floored at 0, so it resets
## whenever cumulative evidence goes negative) -- this naturally links
## neighbouring markers with even individually-weak signal into one
## contiguous peak, using the data's own LD/autocorrelation rather than an
## externally-imposed LD-clustering rule. An analytic significance
## threshold (Fariello et al. 2017 for xi in {1,2}; Bonhomme et al. 2019 for
## xi in {3,4}) is derived per chromosome from its length and
## autocorrelation, and every excursion of the Lindley process above that
## threshold is reported as one candidate window (chr/start/end/size/nsnps
## + the peak position of both the Lindley score and the ORIGINAL
## statistic within the window).
##
## This is a DIFFERENT region-definition method from the rest of this
## pipeline (Stage-2 rho05 LD-clustering + null-calibrated floor-survivor
## test): it builds windows directly from the statistic's own spatial
## profile and gets its significance threshold analytically rather than
## from the 10,000-Omega-structured-null simulation. Run here as an
## independent cross-check against the existing floor-survivor regions,
## using BayPass's own recommended defaults (xi=1, pval.local.score.thres=
## 0.01, min.maf=0.2, min.nsnp=1e4) -- not yet folded into the manuscript.
##
## mitoC2 uses -log10(p) (the manual's recommended input for C2, avoiding
## the negative-BF-to-uniform-draw step); PC1/PC2/bio_winter use BF(dB)
## directly. bio_winter's full-SNP run has no persisted pi_xtx output (only
## BF), so its Pi is taken from PC1's contemporaneous pi_xtx (same
## genotypes/Omega; Pi differs from PC1's own posterior mean only by MCMC
## noise, ~3% mean relative difference between independent BayPass runs on
## the same data) -- documented here as an approximation, not exact.
##
## Reads : baypass_public-master/utils/baypass_utils.R (compute.local.scores)
##         data/hybrids_only_maf005.Rdata (map_hyb_005, full-SNP order)
##         module_manuscript_rho05/baypass_stage1/aland_excluded/
##           {PC1,PC2,mito_C2}_fullSNP_stage1Omega_summary_{betai_reg,pi_xtx,contrast}.out
##           bio_winter_fullSNP_stage1Omega_summary_betai_reg.out
## Writes: module_localscore_crosscheck/data/localscore_<tag>.rds
##         (list(res.local.scores, significant.windows) per covariate)
##         module_localscore_crosscheck/data/localscore_summary.tsv
##
## Run from the repo root:
##   Rscript module_localscore_crosscheck/R/local_score_regions.R
## =========================================================

suppressMessages(library(data.table))
source("~/gitlab/baypass_public-master/utils/baypass_utils.R")

## compute.local.scores() draws a RANDOM p-value (-log10(runif())) for every
## NEGATIVE BF value (BayPass's own code, not this script's), so re-running
## this script can change the exact window count/boundaries for the
## BF-based covariates from one run to the next. Seeded here so the numbers
## reported in doc/ are exactly reproducible; see doc/ for a worked example
## of the resulting run-to-run instability itself.
set.seed(1)

BP_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded"
DATA   <- "module_localscore_crosscheck/data"
dir.create(DATA, showWarnings = FALSE, recursive = TRUE)

load("data/hybrids_only_maf005.Rdata")   # map_hyb_005 (Chr, Pos, marker), full-SNP MRK order
pos <- data.frame(Chr = map_hyb_005$Chr, Pos = map_hyb_005$Pos)

run_bf <- function(tag, bf_file, pi_file) {
  message("\n=== ", tag, " (BF-based) ===")
  bf <- fread(bf_file, select = c("MRK", "BF(dB)"))
  pi <- fread(pi_file, select = c("MRK", "M_P"))
  stopifnot(all(bf$MRK == seq_len(nrow(bf))), all(pi$MRK == seq_len(nrow(pi))),
            nrow(bf) == nrow(pos), nrow(pi) == nrow(pos))
  compute.local.scores(snp.position = pos, snp.pi = pi$M_P, snp.bf = bf$`BF(dB)`,
                       xi = 1, pval.local.score.thres = 0.01, min.maf = 0.2, min.nsnp = 1e4)
}

run_pval <- function(tag, pval_file, pi_file) {
  message("\n=== ", tag, " (p-value-based) ===")
  pv <- fread(pval_file, select = c("MRK", "log10(1/pval)"))
  pi <- fread(pi_file, select = c("MRK", "M_P"))
  stopifnot(all(pv$MRK == seq_len(nrow(pv))), all(pi$MRK == seq_len(nrow(pi))),
            nrow(pv) == nrow(pos), nrow(pi) == nrow(pos))
  compute.local.scores(snp.position = pos, snp.pi = pi$M_P, snp.pvalue = pv$`log10(1/pval)`,
                       xi = 1, pval.local.score.thres = 0.01, min.maf = 0.2, min.nsnp = 1e4)
}

res_pc1 <- run_bf("PC1", file.path(BP_DIR, "PC1_fullSNP_stage1Omega_summary_betai_reg.out"),
                  file.path(BP_DIR, "PC1_fullSNP_stage1Omega_summary_pi_xtx.out"))
res_pc2 <- run_bf("PC2", file.path(BP_DIR, "PC2_fullSNP_stage1Omega_summary_betai_reg.out"),
                  file.path(BP_DIR, "PC2_fullSNP_stage1Omega_summary_pi_xtx.out"))
res_bw  <- run_bf("bio_winter", file.path(BP_DIR, "bio_winter_fullSNP_stage1Omega_summary_betai_reg.out"),
                  file.path(BP_DIR, "PC1_fullSNP_stage1Omega_summary_pi_xtx.out"))  # Pi reused from PC1 (see header)
res_c2  <- run_pval("mitoC2", file.path(BP_DIR, "mito_C2_fullSNP_stage1Omega_summary_contrast.out"),
                    file.path(BP_DIR, "mito_C2_fullSNP_stage1Omega_summary_pi_xtx.out"))

all_res <- list(PC1 = res_pc1, PC2 = res_pc2, bio_winter = res_bw, mitoC2 = res_c2)
for (tag in names(all_res)) saveRDS(all_res[[tag]], file.path(DATA, sprintf("localscore_%s.rds", tag)))

## ---- summary across covariates --------------------------------------------
summ <- rbindlist(lapply(names(all_res), function(tag) {
  w <- all_res[[tag]]$significant.windows
  if (is.null(w) || nrow(w) == 0) return(data.table(covariate = tag, n_windows = 0L))
  dt <- as.data.table(w)
  dt[, covariate := tag]
  dt
}), fill = TRUE)
fwrite(summ, file.path(DATA, "localscore_summary.tsv"), sep = "\t")

cat("\n=== Local-score significant windows, all covariates ===\n")
print(summ)
cat("\n[local-score] wrote per-covariate RDS + summary TSV to ", DATA, "\n")
