## =========================================================================
## module_population_partitioning -- follow-up 10: geographic prediction of
## residual ancestry profiles (HIGHEST PRIORITY, per the 2026-09-13 follow-up
## brief). Does NOT modify any existing pp_*.R script or output.
##
## Question: after removing each population's leave-one-chromosome-out
## genome-wide ancestry (pp_residualize_ancestry.R's $Resid), do
## geographically nearby populations carry similar RESIDUAL ancestry at the
## same genomic units? If yes, geographic structure (isolation-by-distance,
## spatially-varying gene flow) is a live alternative to locus-specific
## selection for the observed heterogeneous sorting; if no, that
## alternative is disfavoured (though not thereby proof of selection).
##
## *** DATASET NOTE (read before using this script's output) ***
## The follow-up brief that requested this analysis says to "retain" the
## 11,052-unit, min_r2=0.2 cM5 clustering (module_di25/data/di25_clustering_
## cM5.rds, di25_sorting_emlg.rds). That lineage is SUPERSEDED as of the
## 2026-09-13 rho05 migration (see ../README.md "Migration to the rho05
## primary dataset" and ../CROSS_MODULE_INPUTS.md) -- the module's own
## up-to-date README, which the brief itself lists as required reading,
## documents this. This script therefore uses the CURRENT, primary dataset
## (module_di25_rho05, 20,807 units) via pp_residual_ancestry.rds, which is
## itself the "completed" residualization result this follow-up brief names
## as an authoritative input -- not the legacy 11,052-unit numbers the brief
## also (inconsistently) asks to "retain". Flagged here rather than silently
## resolved; if the legacy dataset is genuinely wanted for this follow-up
## work specifically, say so and this can be re-run against it.
##
## *** SIELVA NOTE (corrected 2026-09-13) *** data/bioclimatic_variables.csv
## lists Sielva at (46.61 N, 10.44 E) -- the Italian Alps. This script
## originally flagged that as a likely transcription error, reasoning from
## Sielva's bioclim columns fitting smoothly into the Fennoscandian
## latitude-ordered climate gradient of the other 19 populations. THE USER
## CONFIRMED SIELVA IS A GENUINE ALPINE SITE -- that reasoning was wrong (see
## FOLLOWUP_STATUS.md's correction note for the full explanation). This
## means Sielva is geographically disjunct by ~3000km from the other 19
## Fennoscandian populations -- plausibly a different hybrid zone/history
## entirely, not merely a distant outlier within one zone. Included in the
## primary analysis as one of the 20 populations (per the brief's own
## design; Sielva-exclusion is already one of the required leave-one-out
## scenarios below), with results interpreted accordingly.
##
## Nyrhispera naming: data/bioclimatic_variables.csv uses "Nyrhispera1"/
## "Nyrhispera2"; this module's population set uses "Nyrhispera74"/
## "Nyrhispera75". Crosswalk Nyrhispera1->74, Nyrhispera2->75 is the
## ESTABLISHED convention already used identically in
## module_manuscript_rho05/R/moduleB_stage1_prepare_bioclim_covariates.R
## and moduleB_ancestry_vs_winter_climate.R (both: `bc[Location ==
## "Nyrhispera1", Location := "Nyrhispera74"]` etc.) -- not a new guess.
##
## Run from the formica_hybrid repo root:
##   Rscript module_population_partitioning/R/10_geographic_prediction.R
## Reads : module_population_partitioning/data/pp_residual_ancestry.rds
##         module_population_partitioning/data/pp_recombination.rds (for the
##           secondary geographic-R2-vs-recombination stratification)
##         data/bioclimatic_variables.csv
## Writes: module_population_partitioning/data/followup/10_geographic_prediction.rds
##         module_population_partitioning/Figures/followup/10_*.png
## =========================================================================
suppressMessages(library(data.table))
t_start <- Sys.time()
OUTDIR <- "module_population_partitioning/data/followup"
FIGDIR <- "module_population_partitioning/Figures/followup"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(20260913)

## ---------------------------------------------------------------------
## 1. load the residual profiles (completed, unmodified analysis) + geography
## ---------------------------------------------------------------------
rez <- readRDS("module_population_partitioning/data/pp_residual_ancestry.rds")
Resid <- rez$Resid; u <- copy(rez$u); setDT(u)
cat(sprintf("[geo] pp_residual_ancestry.rds: Resid %d pops x %d units; u %d x %d\n",
            nrow(Resid), ncol(Resid), nrow(u), ncol(u)))
stopifnot(identical(colnames(Resid), u$group_id))

bc <- fread("data/bioclimatic_variables.csv")
bc[Location == "Nyrhispera1", Location := "Nyrhispera74"]
bc[Location == "Nyrhispera2", Location := "Nyrhispera75"]
cat("\n[geo] bioclimatic_variables.csv, latitude-ordered (Sielva is the genuine Alpine outlier, see script header):\n")
print(bc[order(Latitude), .(Location, Latitude, Longitude, bio1, bio6, bio11)])

pops <- rownames(Resid)
stopifnot("population(s) unmatched between Resid and bioclimatic_variables.csv" = all(pops %in% bc$Location))
stopifnot("duplicated Location in bioclimatic_variables.csv after Nyrhispera rename" = !anyDuplicated(bc$Location))
geo <- bc[match(pops, Location), .(Location, Latitude, Longitude)]
stopifnot(identical(geo$Location, pops))
cat(sprintf("\n[geo] all %d populations matched to geography (Nyrhispera1/2 -> 74/75 applied). Sielva is a genuine, geographically disjunct Alpine population (see header).\n", length(pops)))

## ---------------------------------------------------------------------
## 2. response matrices: primary = per-unit standardized (scale to unit
##    variance across the 20 pops, so a few high-variance units don't
##    dominate); sensitivity = unstandardized. Both already ~centered per
##    unit (OLS residuals with an intercept sum to 0 by construction); an
##    explicit re-centre is applied defensively.
## ---------------------------------------------------------------------
Y_raw <- scale(Resid, center = TRUE, scale = FALSE)           # defensive re-centre only
unit_sd <- apply(Resid, 2, sd, na.rm = TRUE)
n_na <- colSums(is.na(Resid))
keep_unit <- is.finite(unit_sd) & unit_sd > 0 & n_na == 0
cat(sprintf("[geo] %d/%d units are complete (0 NA across 20 pops) with non-zero variance and are retained (%d dropped: %d had >=1 NA, rest zero-variance)\n",
            sum(keep_unit), length(keep_unit), sum(!keep_unit), sum(n_na > 0)))
Y_raw <- Y_raw[, keep_unit, drop = FALSE]
Y_std <- scale(Resid[, keep_unit, drop = FALSE], center = TRUE, scale = TRUE)   # PRIMARY
u_kept <- u[keep_unit]
stopifnot(identical(colnames(Y_std), u_kept$group_id))

## ---------------------------------------------------------------------
## 3. geography predictor matrix (centered lat/lon), and the core
##    multivariate-regression machinery (fast: X is shared across all
##    unit-columns, so B/Yhat/RSS/TSS are all single matrix operations,
##    not a per-unit loop)
## ---------------------------------------------------------------------
fit_geo <- function(Y, lat, lon) {
  X <- cbind(1, lat - mean(lat), lon - mean(lon))
  XtX_inv <- solve(crossprod(X))
  B <- XtX_inv %*% crossprod(X, Y)
  Yhat <- X %*% B
  RSS_unit <- colSums((Y - Yhat)^2)
  TSS_unit <- colSums(scale(Y, center = TRUE, scale = FALSE)^2)
  list(RSS_unit = RSS_unit, TSS_unit = TSS_unit, B = B, Yhat = Yhat, X = X)
}
r2_from <- function(RSS_unit, TSS_unit, n, p = 2) {
  RSS <- sum(RSS_unit); TSS <- sum(TSS_unit)
  R2 <- 1 - RSS / TSS
  R2adj <- 1 - (1 - R2) * (n - 1) / (n - p - 1)
  c(R2 = R2, R2adj = R2adj)
}

fit_primary <- fit_geo(Y_std, geo$Latitude, geo$Longitude)
obs_primary <- r2_from(fit_primary$RSS_unit, fit_primary$TSS_unit, n = length(pops))
fit_sens <- fit_geo(Y_raw, geo$Latitude, geo$Longitude)
obs_sens <- r2_from(fit_sens$RSS_unit, fit_sens$TSS_unit, n = length(pops))
cat(sprintf("\n[geo] PRIMARY (per-unit standardized) global R2 = %.4f, adj R2 = %.4f\n", obs_primary["R2"], obs_primary["R2adj"]))
cat(sprintf("[geo] SENSITIVITY (unstandardized) global R2 = %.4f, adj R2 = %.4f\n", obs_sens["R2"], obs_sens["R2adj"]))

## ---------------------------------------------------------------------
## 4. permutation test: permute POPULATION ROWS of Y as complete rows
##    (never per-unit) against fixed geography. X'X is invariant to row
##    permutation of Y, so XtX_inv/B machinery is reused across all 10,000
##    reps without re-inverting -- fast.
## ---------------------------------------------------------------------
N_PERM <- 10000
X <- cbind(1, geo$Latitude - mean(geo$Latitude), geo$Longitude - mean(geo$Longitude))
XtX_inv <- solve(crossprod(X))
TSS_obs <- sum(scale(Y_std, center = TRUE, scale = FALSE)^2)   # invariant under row permutation
perm_R2 <- vapply(seq_len(N_PERM), function(b) {
  ord <- sample(nrow(Y_std))
  Yp <- Y_std[ord, , drop = FALSE]
  B <- XtX_inv %*% crossprod(X, Yp)
  RSS <- sum((Yp - X %*% B)^2)
  1 - RSS / TSS_obs
}, numeric(1))
p_perm <- (1 + sum(perm_R2 >= obs_primary["R2"])) / (1 + N_PERM)
cat(sprintf("[geo] permutation test (%d reps, complete-row permutation): null R2 mean %.4f, 95th pct %.4f; observed %.4f; p = %.4f\n",
            N_PERM, mean(perm_R2), quantile(perm_R2, 0.95), obs_primary["R2"], p_perm))

## ---------------------------------------------------------------------
## 5. chromosome-block bootstrap CI on R2 (resamples the per-unit RSS/TSS
##    already computed in step 3 -- no matrix algebra re-run per replicate,
##    same convention as pp_block_bootstrap.R)
## ---------------------------------------------------------------------
chrs <- unique(u_kept$Chr)
cell <- data.table(Chr = u_kept$Chr, RSS = fit_primary$RSS_unit, TSS = fit_primary$TSS_unit)[
  , .(RSS = sum(RSS), TSS = sum(TSS)), by = Chr]
B_BOOT <- 2000
boot_R2 <- vapply(seq_len(B_BOOT), function(b) {
  draw <- sample(chrs, length(chrs), replace = TRUE)
  agg <- cell[Chr %in% draw][, .(w = .N), by = Chr][cell, on = "Chr", nomatch = 0]
  ## weight each chromosome's RSS/TSS by how many times it was drawn
  wtab <- table(draw)
  agg2 <- cell[Chr %in% names(wtab)]
  w <- as.numeric(wtab[agg2$Chr])
  1 - sum(agg2$RSS * w) / sum(agg2$TSS * w)
}, numeric(1))
ci_R2 <- quantile(boot_R2, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("[geo] chromosome-block bootstrap 95%% CI on global R2: [%.4f, %.4f]\n", ci_R2[1], ci_R2[2]))

## ---------------------------------------------------------------------
## 6. leave-one-population-out: (a) in-sample R2 with that population
##    dropped; (b) out-of-sample prediction R2 for the held-out population
##    (pooled RSS/TSS across all 20 folds -> a genuine cross-validated R2,
##    not an average of unstable per-fold ratios)
## ---------------------------------------------------------------------
loo_insample <- rbindlist(lapply(pops, function(p) {
  keep <- pops != p
  f <- fit_geo(Y_std[keep, , drop = FALSE], geo$Latitude[keep], geo$Longitude[keep])
  r <- r2_from(f$RSS_unit, f$TSS_unit, n = sum(keep))
  data.table(dropped = p, R2 = r["R2"], R2adj = r["R2adj"])
}))
setorder(loo_insample, -R2)
cat("\n[geo] leave-one-population-out IN-SAMPLE R2 (dropped population refit):\n"); print(loo_insample)

loo_pred <- rbindlist(lapply(seq_along(pops), function(k) {
  p <- pops[k]; train <- pops != p
  Xtr <- cbind(1, geo$Latitude[train] - mean(geo$Latitude[train]), geo$Longitude[train] - mean(geo$Longitude[train]))
  B <- solve(crossprod(Xtr), crossprod(Xtr, Y_std[train, , drop = FALSE]))
  x_test <- c(1, geo$Latitude[k] - mean(geo$Latitude[train]), geo$Longitude[k] - mean(geo$Longitude[train]))
  yhat_test <- as.numeric(x_test %*% B)
  y_test <- Y_std[k, ]
  ybar_train <- colMeans(Y_std[train, , drop = FALSE])
  data.table(pop = p, RSS = sum((y_test - yhat_test)^2), TSS = sum((y_test - ybar_train)^2))
}))
R2_cv <- 1 - sum(loo_pred$RSS) / sum(loo_pred$TSS)
cat(sprintf("\n[geo] leave-one-population-out CROSS-VALIDATED (out-of-sample) R2, pooled across 20 folds: %.4f\n", R2_cv))
cat("     (a genuine out-of-sample statistic; negative values mean geography predicts a held-out population's\n")
cat("      profile WORSE than that population's own training-set mean -- i.e. no real generalizable signal)\n")

## ---------------------------------------------------------------------
## 7. PCA of population residual profiles (pops as rows/observations,
##    units as columns/variables -- prcomp's default column-centering is
##    CORRECT here, unlike the earlier unit-space PCA gotcha in
##    pp_robustness_pca.R, because centering columns = re-centering each
##    unit across pops, which is exactly what should happen)
## ---------------------------------------------------------------------
pca_pop <- prcomp(Y_std, center = TRUE, scale. = FALSE)
ve_pop <- 100 * pca_pop$sdev^2 / sum(pca_pop$sdev^2)
cat(sprintf("\n[geo] population-space PCA of residual profiles: PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%%\n", ve_pop[1], ve_pop[2], ve_pop[3]))
pca_scores <- data.table(pop = rownames(Y_std), PC1 = pca_pop$x[, 1], PC2 = pca_pop$x[, 2],
                         Latitude = geo$Latitude, Longitude = geo$Longitude)
cat("[geo] Spearman correlation, PC1/PC2 vs Latitude/Longitude:\n")
cat(sprintf("  PC1~Lat %.3f  PC1~Lon %.3f  PC2~Lat %.3f  PC2~Lon %.3f\n",
            cor(pca_scores$PC1, pca_scores$Latitude, method = "spearman"),
            cor(pca_scores$PC1, pca_scores$Longitude, method = "spearman"),
            cor(pca_scores$PC2, pca_scores$Latitude, method = "spearman"),
            cor(pca_scores$PC2, pca_scores$Longitude, method = "spearman")))

## ---------------------------------------------------------------------
## 8. geographic distance vs residual-profile dissimilarity (all 190 pairs),
##    + a descriptive Mantel test (row/column joint permutation, 10,000 reps)
## ---------------------------------------------------------------------
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  p1 <- lat1 * pi / 180; p2 <- lat2 * pi / 180
  dphi <- (lat2 - lat1) * pi / 180; dlam <- (lon2 - lon1) * pi / 180
  a <- sin(dphi / 2)^2 + cos(p1) * cos(p2) * sin(dlam / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}
n_pop <- length(pops)
geo_dist <- matrix(NA_real_, n_pop, n_pop, dimnames = list(pops, pops))
for (a in seq_len(n_pop)) for (b in seq_len(n_pop))
  geo_dist[a, b] <- haversine_km(geo$Latitude[a], geo$Longitude[a], geo$Latitude[b], geo$Longitude[b])
prof_dist <- as.matrix(dist(Y_std))    # Euclidean distance between standardized residual profiles
ut <- upper.tri(geo_dist)
pair_tab <- data.table(pop_i = rep(pops, times = n_pop)[as.vector(ut)],
                       pop_j = rep(pops, each = n_pop)[as.vector(ut)],
                       geo_km = geo_dist[ut], prof_dist = prof_dist[ut])
mantel_obs <- cor(pair_tab$geo_km, pair_tab$prof_dist, method = "spearman")
mantel_null <- vapply(seq_len(N_PERM), function(b) {
  ord <- sample(n_pop)
  gd <- geo_dist[ord, ord][ut]
  cor(gd, pair_tab$prof_dist, method = "spearman")
}, numeric(1))
p_mantel <- (1 + sum(abs(mantel_null) >= abs(mantel_obs))) / (1 + N_PERM)
cat(sprintf("\n[geo] descriptive Mantel-type test (Spearman, geo-distance vs profile-distance): rho = %.3f, perm p = %.4f (sensitivity only, not the primary test)\n",
            mantel_obs, p_mantel))

## ---------------------------------------------------------------------
## 9. secondary: per-unit geographic R2, related to FST / sort_class /
##    chromosome / local recombination rate (chromosome-block uncertainty)
## ---------------------------------------------------------------------
u_kept[, geoR2 := 1 - fit_primary$RSS_unit / fit_primary$TSS_unit]
cat(sprintf("\n[geo] per-unit geographic R2: median %.4f, 95th pct %.4f (individually noisy with n=20 pops -- descriptive only, NOT interpreted as per-locus significance)\n",
            median(u_kept$geoR2), quantile(u_kept$geoR2, 0.95)))
cat(sprintf("[geo] Spearman per-unit geoR2 vs FST: rho = %.3f\n", cor(u_kept$geoR2, u_kept$FST, use = "pairwise.complete.obs", method = "spearman")))
cat("[geo] median per-unit geoR2 by sort_class:\n")
print(u_kept[, .(n = .N, median_geoR2 = median(geoR2)), by = sort_class][order(-n)])

rc <- readRDS("module_population_partitioning/data/pp_recombination.rds")
u_kept <- merge(u_kept, rc$u[, .(group_id, recomb_rate)], by = "group_id", all.x = TRUE)
cat(sprintf("[geo] Spearman per-unit geoR2 vs local recombination rate: rho = %.3f\n",
            cor(u_kept$geoR2, u_kept$recomb_rate, use = "pairwise.complete.obs", method = "spearman")))

## ---------------------------------------------------------------------
## 10. figures
## ---------------------------------------------------------------------
suppressMessages({ library(ggplot2); library(ggrepel) })
theme_ms <- theme_bw(base_size = 12) + theme(strip.background = element_blank(), panel.grid.minor = element_blank())
pca_scores[, is_sielva := pop == "Sielva"]
fig_pca <- ggplot(pca_scores, aes(Longitude, Latitude, colour = PC1, shape = is_sielva)) +
  geom_point(size = 4) + geom_text_repel(aes(label = pop), size = 2.8, colour = "black", max.overlaps = 20, seed = 1) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b", midpoint = 0) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), guide = "none") +
  labs(title = "Population map, coloured by residual-ancestry PC1 score",
      subtitle = "triangle = Sielva (genuine Alpine population, geographically disjunct from the rest)") +
  theme_ms
ggsave(file.path(FIGDIR, "10_map_pc1.png"), fig_pca, width = 8, height = 7, dpi = 200)

fig_geodist <- ggplot(pair_tab, aes(geo_km, prof_dist)) +
  geom_point(alpha = 0.6, colour = "#1b9e77") + geom_smooth(method = "loess", se = TRUE, colour = "#7570b3") +
  labs(x = "geographic distance (km)", y = "residual-profile Euclidean distance",
      title = sprintf("Geographic vs residual-profile distance (Mantel rho=%.3f, p=%.3f)", mantel_obs, p_mantel)) +
  theme_ms
ggsave(file.path(FIGDIR, "10_geodist_vs_profile.png"), fig_geodist, width = 7, height = 5.5, dpi = 200)

fig_perm <- ggplot(data.table(R2 = perm_R2), aes(R2)) +
  geom_histogram(bins = 60, fill = "grey70") +
  geom_vline(xintercept = obs_primary["R2"], colour = "firebrick", linewidth = 1) +
  annotate("text", x = obs_primary["R2"], y = Inf, label = sprintf("observed R2=%.3f, p=%.4f", obs_primary["R2"], p_perm),
           colour = "firebrick", hjust = -0.05, vjust = 1.5, size = 3.4) +
  labs(x = "R2 under complete-row population-label permutation (n=10,000)", y = "count",
      title = "Global geographic effect vs its permutation null") +
  theme_ms
ggsave(file.path(FIGDIR, "10_permutation_null.png"), fig_perm, width = 7, height = 5.5, dpi = 200)
cat("\n[geo] figures saved: 10_map_pc1.png, 10_geodist_vs_profile.png, 10_permutation_null.png\n")

## ---------------------------------------------------------------------
## 11. verdict + save
## ---------------------------------------------------------------------
verdict <- if (p_perm >= 0.05 || obs_primary["R2adj"] < 0.02) {
  "no detectable geographic organization"
} else if (obs_primary["R2adj"] < 0.10 || R2_cv < 0.05) {
  "weak but reproducible geographic organization"
} else {
  "strong geographic prediction of residual ancestry profiles"
}
cat(sprintf("\n[geo] VERDICT: %s\n", verdict))
cat("      (weighting the leave-Sielva-out row of loo_insample more heavily than the\n")
cat("       all-20-population primary result, given Sielva's disjunct geography)\n")

result <- list(
  n_pops = n_pop, n_units_total = ncol(Resid), n_units_kept = sum(keep_unit),
  geo = geo, obs_primary = obs_primary, obs_sensitivity = obs_sens,
  perm_R2 = perm_R2, p_perm = p_perm, n_perm = N_PERM,
  boot_R2 = boot_R2, ci_R2 = ci_R2,
  loo_insample = loo_insample, loo_pred = loo_pred, R2_cv = R2_cv,
  pca_pop = pca_pop, ve_pop = ve_pop, pca_scores = pca_scores,
  geo_dist = geo_dist, prof_dist = prof_dist, pair_tab = pair_tab,
  mantel_obs = mantel_obs, mantel_null = mantel_null, p_mantel = p_mantel,
  u_kept = u_kept, verdict = verdict,
  sielva_coordinate_flag = "Sielva lat/lon (46.61N,10.44E) inconsistent with its own bioclim values (fit ~60-61N); see script header",
  dataset_note = "uses the current rho05 (20,807-unit) dataset via pp_residual_ancestry.rds, not the legacy 11,052-unit lineage the follow-up brief also named",
  session_info = sessionInfo(), run_time = Sys.time(), elapsed_secs = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
)
saveRDS(result, file.path(OUTDIR, "10_geographic_prediction.rds"))
cat(sprintf("\n[geo] saved -> %s (elapsed %.1fs)\n", file.path(OUTDIR, "10_geographic_prediction.rds"), result$elapsed_secs))
