## =========================================================
## MODULE C (revised) -- statistic-stability probe
## =========================================================
## The determinism probe showed BayPass BF are NOT bit-reproducible across runs
## (MCMC noise, ~+/-1.4 dB typical). That rules out a strict max|dk|==0 gate. The
## question that actually matters: are the GENOME-WIDE STATISTICS (aggregates over
## 32,840 eMLGs) reproducible across independent MCMC realizations of the SAME
## covariates? If the within-covariate MCMC scatter of a statistic is small
## relative to its BETWEEN-covariate null spread, the null distribution is governed
## by real covariate variation, not MCMC noise, and Module C is valid despite
## per-BF non-determinism.
##
## Runs BayPass TWICE on the same K-covariate efile and, per covariate, computes the
## Module C statistics under each realization. Reports, per statistic:
##   within-cov MCMC SD  = sd(statA - statB) / sqrt(2)   (per-realization noise)
##   between-cov SD      = sd(statA)                       (real null spread)
##   ratio (noise / spread)  -- small => statistics reproducible / null valid
##
## Reads : aland_excluded_eMLG/null/{u_eMLG.geno,omega_mat_omega.out,u_DIEM.size,null_b01.env}
##         moduleC_climate_vs_sorting/data/moduleC_annotations.rds
## Writes: moduleC_climate_vs_sorting/data/stat_stability_probe.rds  (+ prints table)
##         aland_excluded_eMLG/null/sprobe_{a,b}*  (deleted on success)
## Run from the formica_hybrid repo root. ~15-25 min.

suppressMessages(library(data.table))
source("moduleC_climate_vs_sorting/R/moduleC_stat_functions.R")
ND      <- "aland_excluded_eMLG/null"
BAYPASS <- "/Users/petrikem/gitlab/baypass_public-master/sources/g_baypass"
OPT     <- "-nthreads 6 -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74"
K       <- 20L                                 # covariates to test
stopifnot(file.exists(BAYPASS))

ann <- readRDS("moduleC_climate_vs_sorting/data/moduleC_annotations.rds")
A   <- prepare_annotation_ranks(ann); NM <- A$N

Y  <- as.matrix(fread(file.path(ND, "null_b01.env")))
ef <- file.path(ND, "sprobe.env")
write.table(Y[seq_len(K), , drop = FALSE], ef, quote = FALSE, row.names = FALSE, col.names = FALSE)

run <- function(tag) {
  pref <- file.path(ND, paste0("sprobe_", tag))
  st <- system(paste0(BAYPASS, " -countdatafile ", ND, "/u_eMLG.geno -omegafile ", ND,
    "/omega_mat_omega.out -efile ", ef, " -poolsizefile ", ND, "/u_DIEM.size ",
    OPT, " -outprefix ", pref, " > ", pref, "_stdout.log 2>&1"))
  stopifnot("stability probe BayPass failed" = st == 0)
  nb <- fread(paste0(pref, "_summary_betai_reg.out"), select = c("COVARIABLE", "MRK", "BF(dB)"))
  matrix(nb$`BF(dB)`, nrow = NM, ncol = K)
}
message("[sprobe] run A ..."); MA <- run("a")
message("[sprobe] run B ..."); MB <- run("b")

SA <- t(apply(MA, 2, compute_covariate_stats, A = A))   # K x NSTAT
SB <- t(apply(MB, 2, compute_covariate_stats, A = A))
stats <- colnames(SA)
tab <- data.table(
  statistic   = stats,
  within_mcmc_sd = apply(SA - SB, 2, sd) / sqrt(2),
  between_cov_sd = apply(SA, 2, sd),
  max_abs_AB     = apply(abs(SA - SB), 2, max))
tab[, ratio := within_mcmc_sd / between_cov_sd]
saveRDS(list(SA = SA, SB = SB, tab = tab, K = K,
             bf_max_abs_diff = max(abs(MA - MB))),
        "moduleC_climate_vs_sorting/data/stat_stability_probe.rds")

cat(sprintf("\n[sprobe] K=%d covariates, 2 MCMC realizations. per-BF max|A-B| = %.2f dB\n",
            K, max(abs(MA - MB))))
cat("[sprobe] statistic reproducibility (within-covariate MCMC noise vs between-covariate spread):\n")
print(tab[statistic %in% c("rho_DI", "rho_rec", "sort_gap_differentiated", "sort_gap_all",
                           "rho_sort_magnitude", "rho_sort_orientation"),
          .(statistic, within_mcmc_sd = signif(within_mcmc_sd, 3),
            between_cov_sd = signif(between_cov_sd, 3),
            ratio = signif(ratio, 3), max_abs_AB = signif(max_abs_AB, 3))])
cat("\n[sprobe] interpretation: ratio << 1 => the null distribution of the statistic is\n",
    "governed by real between-covariate variation, not MCMC noise -> statistics reproducible.\n")

file.remove(Sys.glob(file.path(ND, "sprobe_*"))); file.remove(ef)
