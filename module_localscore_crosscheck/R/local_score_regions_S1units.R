## =========================================================
## module_manuscript_rho05 -- BayPass's local-score outlier-region method
## applied to the Stage-1-DIRECT CLUSTER scan (18,361 units, n_snps>=5),
## not the full-SNP scan
## =========================================================
## Same method as moduleB_stage1_local_score_regions.R (compute.local.scores(),
## Fariello et al. 2017 / Bonhomme et al. 2019), but run on the coarser
## Stage-1-unit resolution used by Module C (one best/representative SNP
## per LD-pruned cluster, core_snp) instead of every genome-wide SNP.
##
## PARAMETER DEVIATION FROM BAYPASS DEFAULTS: min.nsnp defaults to 1e4
## (chromosomes with fewer markers after filtering are excluded), which
## would drop EVERY chromosome here -- the largest chromosome has only
## 1,078 Stage-1 units (vs >=21,000 full SNPs/chromosome in the full-SNP
## scan). Set to min.nsnp=100 instead (every chromosome has >=310 tested
## units) so the method can run at all at this resolution; this is a
## necessary adaptation, not a tuning choice made to change the result.
## xi, pval.local.score.thres, and min.maf are kept at BayPass's defaults
## (1, 0.01, 0.2).
##
## Neighbouring Stage-1 units are, by construction, much LESS correlated
## than neighbouring full SNPs (the whole point of the rho05 LD-pruning
## step) -- the local-score method adapts to this automatically (its
## analytic threshold is derived per-chromosome from the data's own
## marker-to-marker p-value autocorrelation), but the resulting windows
## are expected to look different (typically narrower / fewer merged
## clusters) than the full-SNP local-score windows for the same covariate.
##
## Reads : baypass_public-master/utils/baypass_utils.R (compute.local.scores)
##         module0_ld_pruning/data/pruned_stage1.rds (cl5, core_snp positions)
##         module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/
##           {PC1,PC2,bio_winter}_S1units_withOmega_summary_{betai_reg,pi_xtx}.out
##           mito_C2_S1units_summary_{contrast,pi_xtx}.out
## Writes: module_manuscript_rho05/data/moduleB_stage1_localscore_S1units_<tag>.rds
##         module_manuscript_rho05/data/moduleB_stage1_localscore_S1units_summary.tsv
##
## Run from the repo root:
##   Rscript module_manuscript_rho05/R/moduleB_stage1_local_score_regions_S1units.R
## =========================================================

suppressMessages(library(data.table))
source("~/gitlab/baypass_public-master/utils/baypass_utils.R")

UNIT_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded_S1units"
DATA     <- "module_manuscript_rho05/data"
dir.create(DATA, showWarnings = FALSE, recursive = TRUE)

stage1 <- readRDS("module0_ld_pruning/data/pruned_stage1.rds")
cl <- as.data.table(stage1$clusters)
cl5 <- cl[n_snps >= 5]
group_order <- readLines(file.path(UNIT_DIR, "S1units_group_order.txt"))
stopifnot(identical(group_order, paste0("S1_", cl5$CL_id)))   # confirms MRK row order == cl5 row order

## compute.local.scores() requires snp.position sorted by chromosome+position
## (within chromosome) -- cl5 is in CL_id/discovery order, NOT position
## order, so every per-covariate vector must be re-ordered together before
## calling it (not just the position data frame on its own).
base_pos <- data.table(row = seq_len(nrow(cl5)), Chr = cl5$Chr,
                       Pos = as.integer(sub(".*:", "", cl5$core_snp)))
setorder(base_pos, Chr, Pos)

MIN_NSNP <- 100L   # see header -- default 1e4 would drop every chromosome at this resolution

run_bf <- function(tag, bf_file, pi_file) {
  message("\n=== ", tag, " S1units (BF-based) ===")
  bf <- fread(bf_file, select = c("MRK", "BF(dB)"))
  pi <- fread(pi_file, select = c("MRK", "M_P"))
  stopifnot(all(bf$MRK == seq_len(nrow(bf))), all(pi$MRK == seq_len(nrow(pi))),
            nrow(bf) == nrow(base_pos), nrow(pi) == nrow(base_pos))
  ord <- base_pos$row
  compute.local.scores(snp.position = as.data.frame(base_pos[, .(Chr, Pos)]),
                       snp.pi = pi$M_P[ord], snp.bf = bf$`BF(dB)`[ord],
                       xi = 1, pval.local.score.thres = 0.01, min.maf = 0.2, min.nsnp = MIN_NSNP)
}

run_pval <- function(tag, pval_file, pi_file) {
  message("\n=== ", tag, " S1units (p-value-based) ===")
  pv <- fread(pval_file, select = c("MRK", "log10(1/pval)"))
  pi <- fread(pi_file, select = c("MRK", "M_P"))
  stopifnot(all(pv$MRK == seq_len(nrow(pv))), all(pi$MRK == seq_len(nrow(pi))),
            nrow(pv) == nrow(base_pos), nrow(pi) == nrow(base_pos))
  ord <- base_pos$row
  compute.local.scores(snp.position = as.data.frame(base_pos[, .(Chr, Pos)]),
                       snp.pi = pi$M_P[ord], snp.pvalue = pv$`log10(1/pval)`[ord],
                       xi = 1, pval.local.score.thres = 0.01, min.maf = 0.2, min.nsnp = MIN_NSNP)
}

res_pc1 <- run_bf("PC1", file.path(UNIT_DIR, "PC1_S1units_withOmega_summary_betai_reg.out"),
                  file.path(UNIT_DIR, "PC1_S1units_withOmega_summary_pi_xtx.out"))
res_pc2 <- run_bf("PC2", file.path(UNIT_DIR, "PC2_S1units_withOmega_summary_betai_reg.out"),
                  file.path(UNIT_DIR, "PC2_S1units_withOmega_summary_pi_xtx.out"))
res_bw  <- run_bf("bio_winter", file.path(UNIT_DIR, "bio_winter_S1units_withOmega_summary_betai_reg.out"),
                  file.path(UNIT_DIR, "bio_winter_S1units_withOmega_summary_pi_xtx.out"))
res_c2  <- run_pval("mitoC2", file.path(UNIT_DIR, "mito_C2_S1units_summary_contrast.out"),
                    file.path(UNIT_DIR, "mito_C2_S1units_summary_pi_xtx.out"))

all_res <- list(PC1 = res_pc1, PC2 = res_pc2, bio_winter = res_bw, mitoC2 = res_c2)
for (tag in names(all_res)) saveRDS(all_res[[tag]], file.path(DATA, sprintf("moduleB_stage1_localscore_S1units_%s.rds", tag)))

summ <- rbindlist(lapply(names(all_res), function(tag) {
  w <- all_res[[tag]]$significant.windows
  if (is.null(w) || nrow(w) == 0) return(data.table(covariate = tag, n_windows = 0L))
  dt <- as.data.table(w)
  dt[, covariate := tag]
  dt
}), fill = TRUE)
fwrite(summ, file.path(DATA, "moduleB_stage1_localscore_S1units_summary.tsv"), sep = "\t")

cat("\n=== Local-score significant windows, Stage-1-unit resolution, all covariates ===\n")
print(summ[, .N, by = covariate])
print(summ)
cat("\n[local-score-S1units] wrote per-covariate RDS + summary TSV to ", DATA, "\n")
