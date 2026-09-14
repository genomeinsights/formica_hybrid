## =========================================================================
## module_population_partitioning -- candidate-locus follow-up 16: genome-
## wide population-structure reference, estimated from the CANONICAL
## genome-wide LD-reduced reference panel -- NEVER from candidate loci.
##
## Reference panel: the 17,509-unit LARGE-CLUSTER subset of
## module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05_bestsnp.rds
## ($stats/$geno; the full $rep_snp_all has 661,386 units but only this
## large-cluster subset has a saved genotype matrix -- using it is standard
## practice (well-estimated clusters) and avoids reconstructing dosage for
## 661k markers from scratch). Design choice, documented per the plan.
##
## Population set: 19 hybrid populations, Aland excluded (matches script 15
## and the candidate-side BayPass scan population set).
##
## "Principal population-structure pattern" = PC1 of a row-centered PCA
## (populations as observations, units as features, each unit's own mean
## across the 19 populations subtracted -- i.e. prcomp(Fmat_ref, center=TRUE)
## with populations as ROWS). PC1 is primary; PC1+PC2 saved as a named
## sensitivity variant.
##
## Leave-one-chromosome-out: for each of the 26 chromosomes, PCA is REFIT
## (not resummed) on the reference panel with that chromosome's units
## excluded, so a candidate's own chromosome never contributes to its
## structure covariate. This deliberately ADAPTS, not copies,
## pp_residualize_ancestry.R's LOCO mechanism: PCA is not an additive
## statistic, so a resummation (as LOCO ancestry mean is Deep-computed there)
## is not valid here -- a full refit per excluded chromosome is used instead
## (26 refits, cheap: ~17,500 x 19 each).
##
## Run from the formica_hybrid repo root, after 15_candidate_import.R:
##   Rscript module_population_partitioning/R/16_structure_reference.R
## Reads : module_population_partitioning/data/followup/15_candidate_data.rds
## Writes: module_population_partitioning/data/followup/16_structure_reference.rds
## =========================================================================
suppressMessages(library(data.table))
t_start <- Sys.time()
OUTDIR <- "module_population_partitioning/data/followup"
r15 <- readRDS(file.path(OUTDIR, "15_candidate_data.rds"))
stopifnot("script 15's candidate unit table must have exactly 18,361 rows -- rerun 15_candidate_import.R" = nrow(r15$u) == 18361L,
         "POPS_19 must have exactly 19 populations" = length(r15$POPS_19) == 19L)

## ---------------------------------------------------------------------
## 1. reference panel genotypes -> oriented population x unit frequency
##    matrix (identical formula to script 15 / pp_prep_units.R:120-132)
## ---------------------------------------------------------------------
rp <- readRDS("module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05_bestsnp.rds")
stopifnot("reference panel large-cluster subset must have exactly 17,509 units -- check for a stale input" = nrow(rp$stats) == 17509L,
         "reference geno columns must match stats$group_id" = identical(colnames(rp$geno), rp$stats$group_id))

e2 <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd <- e2$sample_data_with_parents; setDT(sd)
GTs_with_parents <- e2$GTs_with_parents
aqu_ids <- sd[Population == "aquilonia_parent", Sample_ID]
pol_ids <- sd[Population == "polyctena_parent", Sample_ID]

hyb_ids <- rownames(rp$geno)
mm <- match(hyb_ids, sd$Sample_ID)
stopifnot("every reference-panel individual must resolve to sample metadata" = !anyNA(mm))
hyb_pop <- sd$Population[mm]
keep_ind <- hyb_pop %in% r15$POPS_19

E <- rp$geno[keep_ind, , drop = FALSE]
pops19 <- hyb_pop[keep_ind]
pop_mean_dosage <- function(G, pop) {
  levs <- unique(pop)
  P <- matrix(NA_real_, length(levs), ncol(G), dimnames = list(levs, colnames(G)))
  for (k in seq_along(levs)) P[k, ] <- colMeans(G[pop == levs[k], , drop = FALSE], na.rm = TRUE)
  P
}
P_hyb <- pop_mean_dosage(E, pops19) / 2

par_geno <- GTs_with_parents[c(aqu_ids, pol_ids), rp$stats$best_marker, drop = FALSE]
f_aqu_par <- colMeans(par_geno[aqu_ids, , drop = FALSE], na.rm = TRUE) / 2
f_pol_par <- colMeans(par_geno[pol_ids, , drop = FALSE], na.rm = TRUE) / 2
sign_aqu  <- sign(f_aqu_par - f_pol_par)

Fmat_ref <- P_hyb[r15$POPS_19, , drop = FALSE]
flip <- which(sign_aqu < 0); undef <- which(is.na(sign_aqu) | sign_aqu == 0)
if (length(flip))  Fmat_ref[, flip]  <- 1 - Fmat_ref[, flip]
if (length(undef)) Fmat_ref[, undef] <- NA_real_
stopifnot("Fmat_ref columns must match rp$stats$group_id in order" = identical(colnames(Fmat_ref), rp$stats$group_id))

ref_chr <- sub(":.*", "", rp$stats$best_marker)
complete_unit <- colSums(is.na(Fmat_ref)) == 0
cat(sprintf("[structure] reference panel: %d populations x %d units; %d units dropped for missingness (%d complete)\n",
            nrow(Fmat_ref), ncol(Fmat_ref), sum(!complete_unit), sum(complete_unit)))
Fmat_ref <- Fmat_ref[, complete_unit, drop = FALSE]
ref_chr  <- ref_chr[complete_unit]
ref_chrs <- sort(unique(ref_chr))
cat(sprintf("[structure] reference panel spans %d chromosomes (%s)\n", length(ref_chrs), paste(ref_chrs, collapse = ",")))
stopifnot("reference panel must span exactly the same 26 chromosomes as the candidate universe" = identical(ref_chrs, r15$chrs))

## ---------------------------------------------------------------------
## 2. full-panel PCA (populations as observations, units as features,
##    column-centered = each unit's own mean across the 19 pops removed)
## ---------------------------------------------------------------------
pc_full <- prcomp(Fmat_ref, center = TRUE, scale. = FALSE)
ve <- 100 * pc_full$sdev^2 / sum(pc_full$sdev^2)
cat("\n[structure] full reference-panel PCA, variance explained (%):\n")
print(round(setNames(ve[1:5], paste0("PC", 1:5)), 2))
cat(sprintf("[structure] PC1 alone explains %.1f%% of among-population structure variance in the reference panel\n", ve[1]))

## ---------------------------------------------------------------------
## 3. leave-one-chromosome-out PC1(+PC2) population scores: refit PCA once
##    per excluded chromosome (not a resummation -- PCA isn't additive)
## ---------------------------------------------------------------------
H_loco_pc1 <- matrix(NA_real_, nrow(Fmat_ref), length(ref_chrs), dimnames = list(rownames(Fmat_ref), ref_chrs))
H_loco_pc2 <- matrix(NA_real_, nrow(Fmat_ref), length(ref_chrs), dimnames = list(rownames(Fmat_ref), ref_chrs))
for (ch in ref_chrs) {
  keep <- ref_chr != ch
  pc <- prcomp(Fmat_ref[, keep, drop = FALSE], center = TRUE, scale. = FALSE)
  ## sign of prcomp's PCs is arbitrary -- anchor to the full-panel PC so LOCO
  ## scores stay comparable across chromosomes (correlate the LOCO loadings
  ## against the full-panel loadings on the shared unit subset, flip if negative)
  shared_full_load <- pc_full$rotation[keep, 1]
  if (cor(pc$rotation[, 1], shared_full_load) < 0) pc$x[, 1] <- -pc$x[, 1]
  H_loco_pc1[, ch] <- pc$x[, 1]
  if (ncol(pc$x) >= 2) {
    shared_full_load2 <- pc_full$rotation[keep, 2]
    if (cor(pc$rotation[, 2], shared_full_load2) < 0) pc$x[, 2] <- -pc$x[, 2]
    H_loco_pc2[, ch] <- pc$x[, 2]
  }
}
cat(sprintf("\n[structure] leave-one-chromosome-out PC1 refit complete for %d chromosomes\n", length(ref_chrs)))
cat("[structure] full-panel PC1 population scores:\n"); print(round(sort(pc_full$x[, 1], decreasing = TRUE), 3))

## ---------------------------------------------------------------------
## 4. save
## ---------------------------------------------------------------------
result <- list(Fmat_ref = Fmat_ref, ref_chr = ref_chr, ref_chrs = ref_chrs,
              pc_full = pc_full, ve = ve,
              H_loco_pc1 = H_loco_pc1, H_loco_pc2 = H_loco_pc2,
              POPS_19 = r15$POPS_19, n_units_dropped = sum(!complete_unit),
              session_info = sessionInfo(), run_time = Sys.time(),
              elapsed_secs = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
saveRDS(result, file.path(OUTDIR, "16_structure_reference.rds"))
cat(sprintf("\n[structure] saved -> %s (elapsed %.1fs)\n", file.path(OUTDIR, "16_structure_reference.rds"), result$elapsed_secs))
