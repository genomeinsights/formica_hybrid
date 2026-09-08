## =========================================================
## module_manuscript_rho05 -- BayPass inputs, Stage-1-direct methodology
## =========================================================
## Follows LDscnR-paper/module_3sp's finding that Stage-2 over-merging costs
## power for outlier-style scans (dilution, documented in
## LDscnR_manuscript/Supplementary.tex sec:goldilocks-units) -- but adapted
## for BayPass, which (unlike EMMAX+Simes) needs one representative genotype
## PER UNIT FED IN, not per-marker p-values combined after testing. So Stage
## 2 is bypassed for BayPass entirely here, not merely re-parameterised:
##
##   * GRM/Omega basis: ALL Stage-1 clusters (698,251, incl. singletons),
##     represented by core_snp -- mirrors module_3sp's
##     `grm_markers <- unique(na.omit(stage1$pruned))`. Reuses the existing
##     write_baypass_inputs() "pruned genotype file -> BayPass estimates its
##     own Omega" mechanism (aland_excluded/run_baypass.sh step 1); only the
##     marker SET changes, not the mechanism.
##   * Scan units: Stage-1 clusters with >= 5 loci only (18,361 of 698,251;
##     matches the min_n_loci_eMLG = 5 convention already used everywhere
##     else in this pipeline, e.g. eMLG_5loci_0025_cM05). Represented by
##     BEST-SNP (not averaged consensus, module_3sp's own reported choice) --
##     required here because downstream needs exact per-marker metadata (DI,
##     position) attached, which only a real SNP carries.
##
## Stage-1/decay reused UNCHANGED (unseeded, Aug 13) -- per-session decision:
## a fresh seeded regen was checked against this run and found to differ by
## <=1.5% in b, with ld_w and a effectively unchanged (see this session's
## seed-equivalence verification), so re-deriving from scratch was judged not
## worth the multi-hour edge-list rebuild.
##
## Run from the repo root:  Rscript module_manuscript_rho05/R/moduleB_stage1_prepare_baypass_inputs.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")
source("moduleB_climate_GEA/R/moduleB_write_baypass_inputs.R")   # write_baypass_inputs()

MIN_N_LOCI <- 5L
OUTROOT <- "module_manuscript_rho05/baypass_stage1"
OMEGA_DIR <- file.path(OUTROOT, "aland_excluded")
UNIT_DIR  <- file.path(OUTROOT, "aland_excluded_S1units")
dir.create(OMEGA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(UNIT_DIR,  showWarnings = FALSE, recursive = TRUE)

## ---- inputs ----------------------------------------------------------------
load("data/hybrids_only_maf005.Rdata")   # GTs_hybrids_005, map_hyb_005, sample_data, ld_decay
stage1 <- readRDS("module0_ld_pruning/data/pruned_stage1.rds")
cl <- as.data.table(stage1$clusters)
cat("Stage-1 clusters total:", nrow(cl), " (singletons:", sum(cl$n_snps == 1), ")\n")

## =============================================================================
## 1. GRM/Omega basis: ALL Stage-1 clusters, one representative (core_snp) each
## =============================================================================
grm_markers <- cl$core_snp
stopifnot(!anyNA(grm_markers), length(grm_markers) == nrow(cl))
cat("GRM/Omega basis (all Stage-1 core_snp):", length(grm_markers), "markers\n")

write_baypass_inputs(
  GTs = GTs_hybrids_005, map = map_hyb_005, sample_data = sample_data,
  pruned_markers = grm_markers, out_folder = paste0(OMEGA_DIR, "/"),
  exclude_population = "Aland"
)

## =============================================================================
## 2. Scan units: Stage-1 clusters with >= 5 loci, best-SNP representation
## =============================================================================
cl5 <- cl[n_snps >= MIN_N_LOCI]
cat("\nStage-1 clusters with >=", MIN_N_LOCI, "loci:", nrow(cl5), "\n")

message("Computing consensus dosage for ", nrow(cl5), " clusters...")
cons_list <- vector("list", nrow(cl5))
for (i in seq_len(nrow(cl5))) cons_list[[i]] <- consensus_dosage(GTs_hybrids_005, cl5$members[[i]])
eMLG_S1 <- do.call(cbind, cons_list)
colnames(eMLG_S1) <- paste0("S1_", cl5$CL_id)
rownames(eMLG_S1) <- rownames(GTs_hybrids_005)

groups_S1 <- data.table(
  group_id = colnames(eMLG_S1), Chr = cl5$Chr, representative = cl5$core_snp,
  n_loci = cl5$n_snps, score = vapply(cons_list, score_eMLG, numeric(1)),
  has_eMLG = TRUE, members = cl5$members
)

message("Selecting best-SNP per unit (eMLG_best_snp)...")
best_S1 <- eMLG_best_snp(list(eMLG = eMLG_S1, groups = groups_S1),
                         GTs_hybrids_005, fill = TRUE, round_fill = TRUE)
saveRDS(list(groups = groups_S1, best = best_S1),
        "module_manuscript_rho05/data/moduleB_stage1_units_bestsnp.rds")
cat("rep_is_best (best-SNP == centrality representative): ",
    round(100 * mean(best_S1$stats$rep_is_best), 1), "%\n", sep = "")

## ---- write the Stage-1-unit BayPass geno file, reusing Omega/PC/poolsize --
sd_ex <- sample_data[Population != "Aland"]
E <- best_S1$geno[sd_ex$Sample_ID, , drop = FALSE]
stopifnot(all(E %in% c(0L, 1L, 2L, NA)))
pop <- sd_ex$Population
pop_order <- unique(pop)

geno <- do.call(cbind, lapply(pop_order, function(y) {
  t(apply(E[pop == y, , drop = FALSE], 2, function(x) {
    c(sum(x == 0, na.rm = TRUE) * 2 + sum(x == 1, na.rm = TRUE),
      sum(x == 2, na.rm = TRUE) * 2 + sum(x == 1, na.rm = TRUE))
  }))
}))
stopifnot(nrow(geno) == ncol(E), ncol(geno) == 2 * length(pop_order))

poolsize <- as.integer(table(pop)[pop_order])
src_size <- scan(file.path(OMEGA_DIR, "u_DIEM.size"), what = integer(), quiet = TRUE)
stopifnot(
  "rebuilt pool sizes disagree with the omega run's u_DIEM.size -- population order mismatch" =
    identical(poolsize, src_size)
)

write.table(geno, file.path(UNIT_DIR, "u_S1units.geno"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
writeLines(colnames(E), file.path(UNIT_DIR, "S1units_group_order.txt"))
for (f in c("u.PC1", "u.PC2", "u_DIEM.size"))
  file.copy(file.path(OMEGA_DIR, f), file.path(UNIT_DIR, f), overwrite = TRUE)

cat("\nWrote ", UNIT_DIR, "/u_S1units.geno (", nrow(geno), " units x ", ncol(geno),
    " cols)\n", sep = "")
cat("Omega will be estimated in ", OMEGA_DIR, " (step 1 of run_baypass.sh);",
    " copy omega_mat_omega.out into ", UNIT_DIR, " before the unit scan.\n", sep = "")
message("[moduleB-stage1-prepare] done")
