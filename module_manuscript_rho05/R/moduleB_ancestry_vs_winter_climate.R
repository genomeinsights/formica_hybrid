## =========================================================
## module_manuscript_rho05 -- genome-wide ancestry vs winter-temperature
## climate (bio6, bio11), non-independence-aware significance
## =========================================================
## Extends moduleB_climate_GEA/R/moduleB_ancestry_confound.R (which
## established cor(PC1,ancestry)=-0.23/-0.50-no-Aland, cor(PC2,ancestry)=
## +0.57/+0.38-no-Aland) to the new winter-temperature covariates
## (bio6 = min temp of coldest month, bio11 = mean temp of coldest quarter),
## and adds TWO non-independence-aware significance tests the original plain
## Pearson r doesn't have:
##
##   1. Omega-structured null (reuses the Stage-1 Omega already estimated for
##      the BayPass scan, aland_excluded/omega_mat_omega.out): draw 10,000
##      structure-preserving random vectors from Omega's eigendecomposition
##      (same mechanism as moduleB_stage1_S1units_null.R), correlate each
##      against the OBSERVED ancestry vector, empirical p = (1+k)/(nsim+1)
##      where k = # null |r| >= observed |r|. Answers: is this correlation
##      stronger than expected given the populations' allele-frequency
##      covariance alone?
##   2. Shared-origin block permutation: LangholmenW+LangholmenR and
##      Bunkkeri+Grundsund have a documented shared origin (phylogenetic
##      network nesting, Nouhaud et al. 2022) -- treating them as 19
##      independent populations overstates n. Permutes climate values across
##      17 INDEPENDENT lineage units (the two shared pairs collapsed to one
##      unit each), recomputing the correlation each time. Answers: does the
##      correlation survive treating known non-independent pairs as one data
##      point?
##
## Run on the aland_excluded population set (19 pops, matching the Stage-1
## Omega and the BayPass covariate scans) -- ancestry recomputed exactly as
## moduleB_ancestry_confound.R (diagnostic markers, DI>-25 & |dp|>=0.5).
##
## Run from the repo root: Rscript module_manuscript_rho05/R/moduleB_ancestry_vs_winter_climate.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2) })

OMEGA_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded"
FIGDIR <- "module_manuscript_rho05/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
NSIM <- 10000
SEED_SIM <- 2026
SEED_PERM <- 2027

## ---- 1. ancestry, exactly as moduleB_ancestry_confound.R -----------------
e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
sd  <- as.data.table(e$sample_data_with_parents)
GT  <- e$GTs_with_parents
map <- copy(e$map_hyb_005)

aqu <- "aquilonia_parent"; pol <- "polyctena_parent"
h   <- sd[!grepl("_parent$", Population)]
hyb_pops <- unique(h$Population)

pa_v <- colMeans(GT[sd$Population == aqu, ], na.rm = TRUE) / 2
pp_v <- colMeans(GT[sd$Population == pol, ], na.rm = TRUE) / 2
map[, `:=`(pa = pa_v[marker], pp = pp_v[marker])]
map[, sign_aqu := sign(pa - pp)]
diag_mk <- map[DiagnosticIndex > -25 & abs(pa - pp) >= 0.5 & sign_aqu != 0, marker]
diag_mk <- intersect(diag_mk, colnames(GT))
s_aqu   <- setNames(map$sign_aqu, map$marker)[diag_mk]

anc <- vapply(hyb_pops, function(P) {
  f <- colMeans(GT[h[Population == P, Sample_ID], diag_mk, drop = FALSE], na.rm = TRUE) / 2
  mean(ifelse(s_aqu > 0, f, 1 - f), na.rm = TRUE)
}, numeric(1))

dt <- h[, .(PC1 = PC1[1], PC2 = PC2[1]), by = Population]
dt[, ancestry := anc[Population]]

## ---- 2. winter-temperature covariates -------------------------------------
bc <- fread("data/bioclimatic_variables.csv")
bc[Location == "Nyrhispera1", Location := "Nyrhispera74"]
bc[Location == "Nyrhispera2", Location := "Nyrhispera75"]
dt[bc, on = .(Population = Location), `:=`(bio6 = i.bio6, bio11 = i.bio11)]
stopifnot(!anyNA(dt$bio6), !anyNA(dt$bio11))

## ---- 3. restrict to the aland_excluded population set (matches the Omega
## and the BayPass covariate scans) + confirm order matches u_DIEM.size -----
pop_order <- unique(h[Population != "Aland", Population])
src_size <- scan(file.path(OMEGA_DIR, "u_DIEM.size"), what = integer(), quiet = TRUE)
stopifnot(length(pop_order) == length(src_size))
dt19 <- dt[match(pop_order, Population)]
stopifnot(!anyNA(dt19$Population))

cat("=== plain Pearson correlations with ancestry ===\n")
cat("  all 20 populations:\n")
for (v in c("PC1", "PC2", "bio6", "bio11")) cat(sprintf("    cor(%-5s, ancestry) = %+.3f\n", v, cor(dt[[v]], dt$ancestry)))
cat("  aland_excluded (n=19, matches the Omega/BayPass scan set):\n")
for (v in c("PC1", "PC2", "bio6", "bio11")) cat(sprintf("    cor(%-5s, ancestry) = %+.3f\n", v, cor(dt19[[v]], dt19$ancestry)))

## ---- 4. Omega-structured null (reuses the Stage-1 Omega) ------------------
Omega <- as.matrix(fread(file.path(OMEGA_DIR, "omega_mat_omega.out"))); Omega <- (Omega + t(Omega)) / 2
eg <- eigen(Omega, symmetric = TRUE); vals <- pmax(eg$values, 0); P <- nrow(Omega)
stopifnot(P == nrow(dt19))
set.seed(SEED_SIM)
NullMat <- sapply(seq_len(NSIM), function(k) as.numeric(scale(eg$vectors %*% (sqrt(vals) * rnorm(P)))))
null_cor <- function(x) { xz <- as.numeric(scale(x)); as.numeric(cor(NullMat, xz)) }

omega_test <- rbindlist(lapply(c("PC1", "PC2", "bio6", "bio11"), function(v) {
  obs <- cor(dt19[[v]], dt19$ancestry)
  nc <- null_cor(dt19[[v]])
  k <- sum(abs(nc) >= abs(obs))
  data.table(variable = v, n_pop = nrow(dt19), observed_r = round(obs, 3),
             null_median = round(median(nc), 3),
             null_95 = sprintf("[%.3f,%.3f]", quantile(nc, .025), quantile(nc, .975)),
             k = k, p_omega_null = round((1 + k) / (NSIM + 1), 4))
}))
cat("\n=== Omega-structured null (n=", NSIM, " draws, aland_excluded) ===\n", sep = "")
print(omega_test)

## ---- 5. shared-origin block permutation -----------------------------------
## LangholmenW+LangholmenR, Bunkkeri+Grundsund share origins (phylogenetic
## nesting, Nouhaud et al. 2022) -- collapse each pair to ONE independent unit.
unit <- dt19$Population
unit[dt19$Population %in% c("LangholmenW", "LangholmenR")] <- "unit_Lang"
unit[dt19$Population %in% c("Bunkkeri", "Grundsund")]      <- "unit_BunGru"
units <- unique(unit)
cat("\nindependent lineage units:", length(units), "(of", nrow(dt19), "populations)\n")

set.seed(SEED_PERM)
NPERM <- 10000
block_perm_test <- rbindlist(lapply(c("PC1", "PC2", "bio6", "bio11"), function(v) {
  obs <- cor(dt19[[v]], dt19$ancestry)
  ## per-unit climate value (constant within a shared-origin pair already)
  unit_val <- setNames(dt19[[v]][match(units, unit)], units)
  null_r <- vapply(seq_len(NPERM), function(i) {
    perm_units <- sample(units)                       # permute unit->value assignment
    perm_val_by_unit <- setNames(unit_val[perm_units], units)
    perm_climate <- perm_val_by_unit[unit]             # broadcast back to populations
    cor(perm_climate, dt19$ancestry)
  }, numeric(1))
  k <- sum(abs(null_r) >= abs(obs))
  data.table(variable = v, n_units = length(units), observed_r = round(obs, 3),
             perm_median = round(median(null_r), 3),
             perm_95 = sprintf("[%.3f,%.3f]", quantile(null_r, .025), quantile(null_r, .975)),
             k = k, p_block_perm = round((1 + k) / (NPERM + 1), 4))
}))
cat("\n=== shared-origin block permutation (n=", NPERM, " perms, ", length(units), " independent units) ===\n", sep = "")
print(block_perm_test)

## ---- 6. combined summary + figure -----------------------------------------
summ <- omega_test[block_perm_test[, .(variable, perm_median, perm_95, p_block_perm)], on = "variable"]
cat("\n=== combined summary ===\n"); print(summ)
saveRDS(list(dt = dt, dt19 = dt19, omega_test = omega_test, block_perm_test = block_perm_test, summ = summ,
             units = units, NSIM = NSIM, NPERM = NPERM),
        "module_manuscript_rho05/data/moduleB_ancestry_vs_winter_climate.rds")

m <- melt(dt19, id.vars = c("Population", "ancestry"),
          measure.vars = c("bio6", "bio11"), variable.name = "variable", value.name = "value")
m[summ, on = "variable", `:=`(r = i.observed_r, p_om = i.p_omega_null, p_bp = i.p_block_perm)]
m[, lab := sprintf("%s: r=%+.2f, p_Omega=%.3f, p_block=%.3f", variable, r, p_om, p_bp)]
p <- ggplot(m, aes(ancestry, value)) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey70", linewidth = 0.5) +
  geom_point(size = 1.8, colour = "#315B7D") +
  geom_text(aes(label = Population), size = 2.2, vjust = -0.7, colour = "grey30") +
  facet_wrap(~ lab, scales = "free_y") +
  labs(x = "genome-wide aquilonia ancestry (per population, aland_excluded)",
       y = "winter-temperature bioclim value",
       title = "Ancestry vs winter temperature, aland_excluded (n=19)") +
  theme_classic(base_size = 9)
ggsave(file.path(FIGDIR, "moduleB_ancestry_vs_winter_climate.png"), p, width = 8, height = 4.2, dpi = 200)
cat("\nSaved figure + module_manuscript_rho05/data/moduleB_ancestry_vs_winter_climate.rds\n")
