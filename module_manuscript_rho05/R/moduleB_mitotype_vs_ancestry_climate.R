## =========================================================
## module_manuscript_rho05 -- mitotype vs ancestry/PC1/PC2/bio6/bio11
## =========================================================
## Does mitotype (aquilonia-like vs polyctena-like maternal lineage, uniform
## within every aland_excluded population -- see u.mito_contrast) predict
## genome-wide nuclear ancestry or any of the climate variables? Point-
## biserial correlation (= Pearson r with a 0/1 mitotype indicator) against
## the SAME two non-independence-aware tests used for
## moduleB_ancestry_vs_winter_climate.R: the Omega-structured null (reuses
## the Stage-1 Omega) and the shared-origin block permutation (17
## independent lineage units).
##
## Run from the repo root: Rscript module_manuscript_rho05/R/moduleB_mitotype_vs_ancestry_climate.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2) })

OMEGA_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded"
FIGDIR <- "module_manuscript_rho05/Figures"
NSIM <- 10000; NPERM <- 10000
SEED_SIM <- 2026; SEED_PERM <- 2027

obj <- readRDS("module_manuscript_rho05/data/moduleB_ancestry_vs_winter_climate.rds")
dt19 <- copy(obj$dt19)

## ---- mitotype per population (uniform within pop; Aland already excluded) -
load("data/hybrids_only_maf005.Rdata")   # sample_data
mito_by_pop <- unique(sample_data[Population %in% dt19$Population & !is.na(Mitotype),
                                  .(Population, Mitotype)])
stopifnot(uniqueN(mito_by_pop$Population) == nrow(dt19))   # every pop resolved, no duplicates
dt19[mito_by_pop, on = "Population", Mitotype := i.Mitotype]
dt19[, mito01 := as.integer(Mitotype == "Faquilonia")]     # 1 = aquilonia-like, 0 = polyctena-like
cat("mitotype counts (aland_excluded, n=19): ",
    sum(dt19$mito01 == 1), " aquilonia-like, ", sum(dt19$mito01 == 0), " polyctena-like\n\n", sep = "")
print(dt19[, .(Population, Mitotype, ancestry, PC1, PC2, bio6, bio11)])

VARS <- c("ancestry", "PC1", "PC2", "bio6", "bio11")

cat("\n=== plain point-biserial correlations with mitotype (1=aquilonia-like) ===\n")
for (v in VARS) cat(sprintf("  cor(mito01, %-8s) = %+.3f\n", v, cor(dt19$mito01, dt19[[v]])))

## ---- Omega-structured null ------------------------------------------------
Omega <- as.matrix(fread(file.path(OMEGA_DIR, "omega_mat_omega.out"))); Omega <- (Omega + t(Omega)) / 2
eg <- eigen(Omega, symmetric = TRUE); vals <- pmax(eg$values, 0); P <- nrow(Omega)
stopifnot(P == nrow(dt19))
set.seed(SEED_SIM)
NullMat <- sapply(seq_len(NSIM), function(k) as.numeric(scale(eg$vectors %*% (sqrt(vals) * rnorm(P)))))
mito_z <- as.numeric(scale(dt19$mito01))

omega_test <- rbindlist(lapply(VARS, function(v) {
  obs <- cor(dt19$mito01, dt19[[v]])
  ## null distribution for the mito-vs-v correlation: correlate each
  ## structure-preserving null draw against the OBSERVED mitotype vector
  ## (mitotype held fixed; the null draw stands in for "a structure-
  ## consistent v" with no real link to mitotype)
  nc_mito <- as.numeric(cor(NullMat, mito_z))
  k <- sum(abs(nc_mito) >= abs(obs))
  data.table(variable = v, observed_r = round(obs, 3),
             null_median = round(median(nc_mito), 3),
             null_95 = sprintf("[%.3f,%.3f]", quantile(nc_mito, .025), quantile(nc_mito, .975)),
             k = k, p_omega_null = round((1 + k) / (NSIM + 1), 4))
}))
cat("\n=== Omega-structured null (n=", NSIM, " draws) ===\n", sep = "")
print(omega_test)

## ---- shared-origin block permutation --------------------------------------
unit <- dt19$Population
unit[dt19$Population %in% c("LangholmenW", "LangholmenR")] <- "unit_Lang"
unit[dt19$Population %in% c("Bunkkeri", "Grundsund")]      <- "unit_BunGru"
units <- unique(unit)
cat("\nindependent lineage units:", length(units), "\n")

set.seed(SEED_PERM)
block_perm_test <- rbindlist(lapply(VARS, function(v) {
  obs <- cor(dt19$mito01, dt19[[v]])
  unit_val <- setNames(dt19[[v]][match(units, unit)], units)
  null_r <- vapply(seq_len(NPERM), function(i) {
    perm_units <- sample(units)
    perm_val_by_unit <- setNames(unit_val[perm_units], units)
    cor(dt19$mito01, perm_val_by_unit[unit])
  }, numeric(1))
  k <- sum(abs(null_r) >= abs(obs))
  data.table(variable = v, observed_r = round(obs, 3),
             perm_median = round(median(null_r), 3),
             perm_95 = sprintf("[%.3f,%.3f]", quantile(null_r, .025), quantile(null_r, .975)),
             k = k, p_block_perm = round((1 + k) / (NPERM + 1), 4))
}))
cat("\n=== shared-origin block permutation (n=", NPERM, " perms, ", length(units), " units) ===\n", sep = "")
print(block_perm_test)

summ <- omega_test[block_perm_test[, .(variable, perm_median, perm_95, p_block_perm)], on = "variable"]
cat("\n=== combined summary ===\n"); print(summ)
saveRDS(list(dt19 = dt19, omega_test = omega_test, block_perm_test = block_perm_test, summ = summ),
        "module_manuscript_rho05/data/moduleB_mitotype_vs_ancestry_climate.rds")

## ---- figure: each variable split by mitotype ------------------------------
m <- melt(dt19, id.vars = c("Population", "Mitotype"), measure.vars = VARS,
          variable.name = "variable", value.name = "value")
m[summ, on = "variable", `:=`(r = i.observed_r, p_om = i.p_omega_null, p_bp = i.p_block_perm)]
m[, lab := sprintf("%s: r=%+.2f, p_Om=%.3f, p_blk=%.3f", variable, r, p_om, p_bp)]
p <- ggplot(m, aes(Mitotype, value)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, colour = "grey60") +
  geom_jitter(width = 0.08, size = 1.8, colour = "#315B7D") +
  facet_wrap(~ lab, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL, title = "Mitotype vs ancestry/climate variables (aland_excluded, n=19)") +
  theme_classic(base_size = 9)
ggsave(file.path(FIGDIR, "moduleB_mitotype_vs_ancestry_climate.png"), p, width = 13, height = 4, dpi = 200)
cat("\nSaved figure + moduleB_mitotype_vs_ancestry_climate.rds\n")
