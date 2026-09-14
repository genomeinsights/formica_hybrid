## =========================================================================
## module_population_partitioning -- candidate-locus follow-up 21: closing
## instruction of item 5 -- do NOT interpret elevated distant/cross-
## chromosome candidate profile similarity (scripts 18/19/20) as evidence
## of epistasis by itself. Evaluate whether it is (a) explained by a shared
## climate axis (loci that both track the same population-level climate
## gradient will correlate with each other for a reason having nothing to
## do with physical interaction), or (b) concentrated among BDMI-flagged
## loci (which WOULD be consistent with, though not proof of, incompatibility).
##
## Climate covariate: the population-level "bio_winter" axis, computed with
## the IDENTICAL formula module_manuscript_rho05 uses (scale(rowMeans(
## cbind(scale(bio6), scale(bio11)))), documented in that module's own
## moduleB_stage1_prepare_bio_winter_covariate.R) from this repo's
## data/bioclimatic_variables.csv (Nyrhispera1/2 -> 74/75 crosswalk, this
## module's established convention). NOT re-derived differently.
##
## For each target's raw candidate set: (i) each locus's own correlation
## with the population climate axis; (ii) whether cross-chr/same-chr-far
## locus PAIRS with high similarity are disproportionately pairs where BOTH
## loci are themselves strongly climate-correlated; (iii) BDMI-BDMI vs
## BDMI-other vs other-other mean |r| by distance category.
##
## Run from the formica_hybrid repo root, after 20_matched_null.R:
##   Rscript module_population_partitioning/R/21_climate_bdmi_context.R
## Reads : module_population_partitioning/data/followup/{15,18,19,20}_*.rds
##         data/bioclimatic_variables.csv
## Writes: module_population_partitioning/data/followup/21_climate_bdmi_context.rds
##         module_population_partitioning/Figures/followup/21_climate_context_<target>.png
## =========================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
t_start <- Sys.time()
OUTDIR <- "module_population_partitioning/data/followup"
FIGDIR <- "module_population_partitioning/Figures/followup"
r15 <- readRDS(file.path(OUTDIR, "15_candidate_data.rds"))
r18 <- readRDS(file.path(OUTDIR, "18_similarity_raw.rds"))
u <- r15$u; Fmat_all <- r15$Fmat_all
theme_ms <- theme_bw(base_size = 11) + theme(strip.background = element_blank(), panel.grid.minor = element_blank())

## ---------------------------------------------------------------------
## 1. population-level bio_winter climate axis (identical formula to
##    module_manuscript_rho05, Nyrhispera crosswalk per this module's
##    established convention)
## ---------------------------------------------------------------------
bc <- fread("data/bioclimatic_variables.csv")
bc[Location == "Nyrhispera1", Location := "Nyrhispera74"]
bc[Location == "Nyrhispera2", Location := "Nyrhispera75"]
bc <- bc[Location %in% r15$POPS_19]
stopifnot("all 19 populations must be present in bioclimatic_variables.csv" = nrow(bc) == 19L)
bc[, bio_winter := as.numeric(scale(rowMeans(cbind(scale(bio6), scale(bio11)))))]
clim <- setNames(bc$bio_winter, bc$Location)[r15$POPS_19]
cat("[climate] population-level bio_winter climate axis:\n"); print(round(sort(clim, decreasing = TRUE), 3))

## ---------------------------------------------------------------------
## 2. per-locus correlation with the climate axis, per target's raw set
## ---------------------------------------------------------------------
locus_climate_cor <- function(ids) {
  vapply(ids, function(id) suppressWarnings(cor(Fmat_all[, id], clim, use = "pairwise.complete.obs")), numeric(1))
}

## ---------------------------------------------------------------------
## 3. for each target: is high locus-pair similarity concentrated among
##    pairs where BOTH loci are strongly climate-correlated? and among
##    BDMI-BDMI pairs?
## ---------------------------------------------------------------------
run_context <- function(tgt) {
  sr <- r18$sim_results[[paste0(tgt, "_raw")]]
  if (is.null(sr) || is.null(sr$pairs)) { cat(sprintf("[context] %-10s: no similarity result available, skipped\n", tgt)); return(NULL) }
  ids <- sr$ids
  lc <- locus_climate_cor(ids)
  pairs <- copy(sr$pairs)
  pairs[, clim_i := lc[group_id_i]][, clim_j := lc[group_id_j]]
  pairs[, both_climate_strong := abs(clim_i) >= 0.5 & abs(clim_j) >= 0.5]
  bdmi_flag <- setNames(u$bdmi_overlap[match(ids, u$group_id)], ids)
  pairs[, bdmi_i := bdmi_flag[group_id_i]][, bdmi_j := bdmi_flag[group_id_j]]
  pairs[, bdmi_pair_class := fifelse(bdmi_i & bdmi_j, "BDMI-BDMI", fifelse(bdmi_i | bdmi_j, "BDMI-other", "other-other"))]

  long_range <- pairs[dcat != "same-chr near (<=100kb)"]
  clim_strat <- long_range[, .(n = .N, mean_absr = mean(absr, na.rm = TRUE)), by = both_climate_strong][order(-both_climate_strong)]
  bdmi_strat <- long_range[, .(n = .N, mean_absr = mean(absr, na.rm = TRUE)), by = bdmi_pair_class][order(-mean_absr)]
  cat(sprintf("\n[context] %-10s (long-range pairs, n=%d): similarity by joint climate-correlation status:\n", tgt, nrow(long_range)))
  print(clim_strat)
  cat(sprintf("[context] %-10s: similarity by BDMI-pair class:\n", tgt)); print(bdmi_strat)

  df <- data.table(locus_climate_r = lc, group_id = ids)
  p1 <- ggplot(df, aes(locus_climate_r)) + geom_histogram(bins = 20, fill = "grey60") +
    geom_vline(xintercept = c(-0.5, 0.5), colour = "firebrick", linetype = "dashed") +
    labs(title = sprintf("%s raw candidates: per-locus correlation with\npopulation climate axis (bio_winter)", tgt),
        x = "locus-climate Pearson r", y = "n loci") + theme_ms
  p2 <- ggplot(long_range, aes(bdmi_pair_class, absr)) + geom_boxplot() +
    labs(title = "long-range pair similarity by BDMI status", x = NULL, y = "|r| (population-profile correlation)") + theme_ms
  fig <- p1 + p2
  ggsave(file.path(FIGDIR, sprintf("21_climate_context_%s.png", tgt)), fig, width = 10, height = 4.5, dpi = 200)

  list(target = tgt, locus_climate_cor = lc, clim_strat = clim_strat, bdmi_strat = bdmi_strat, n_long_range_pairs = nrow(long_range))
}

context_results <- lapply(names(r15$TARGETS), run_context)
names(context_results) <- names(r15$TARGETS)
context_results <- context_results[!vapply(context_results, is.null, logical(1))]

result <- list(clim = clim, context_results = context_results,
              session_info = sessionInfo(), run_time = Sys.time(),
              elapsed_secs = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
saveRDS(result, file.path(OUTDIR, "21_climate_bdmi_context.rds"))
cat(sprintf("\n[context] saved -> %s (elapsed %.1fs)\n", file.path(OUTDIR, "21_climate_bdmi_context.rds"), result$elapsed_secs))
