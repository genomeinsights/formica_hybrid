## =========================================================
## module_manuscript_rho05 -- BayPass C2 contrast file: mitotype
## =========================================================
## Mitonuclear-incompatibility test (manuscript Methods: "This analysis was
## performed with BayPass assessing genome-wide associations with
## mitochondria as a binary variable") -- BayPass's C2 statistic
## (Olazcuaga et al. 2020), contrasting standardized allele frequencies
## between two population groups, computed under the core model
## (-omegafile <existing Omega> -contrastfile <this file>, no -efile).
##
## sample_data$Mitotype is per-INDIVIDUAL but uniform within every
## population except Aland (7 Faquilonia + 3 Fpolyctena) -- already excluded
## from the aland_excluded population set, so every remaining population has
## an unambiguous group assignment. Verified directly (not assumed) before
## writing: stops if any surviving population turns out non-unanimous.
##
## Contrast coding: 1 = Faquilonia-like mitotype, -1 = Fpolyctena-like.
## One row (single contrast), population columns in the SAME order as
## u_DIEM.size / u.PC1 / u.PC2 in aland_excluded/ (and staged into the
## Stage-1-unit dir too), so it drops straight into either scan.
##
## Run from the repo root:  Rscript module_manuscript_rho05/R/moduleB_stage1_prepare_mito_contrast.R
## =========================================================

suppressMessages(library(data.table))

OMEGA_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded"
UNIT_DIR  <- "module_manuscript_rho05/baypass_stage1/aland_excluded_S1units"

load("data/hybrids_only_maf005.Rdata")   # sample_data (has Mitotype)
sd_ex <- sample_data[Population != "Aland"]
pop_order <- unique(sd_ex$Population)    # must match write_baypass_inputs()'s pop_order exactly

## ---- verify unanimous mitotype per surviving population -------------------
mito_by_pop <- sd_ex[!is.na(Mitotype), .(n_types = uniqueN(Mitotype)), by = Population]
stopifnot(
  "a surviving population has >1 mitotype -- contrast coding below assumes unanimity" =
    all(mito_by_pop$n_types == 1)
)

pop_mito <- sd_ex[!is.na(Mitotype), .(Mitotype = unique(Mitotype)), by = Population]
mito_vec <- setNames(pop_mito$Mitotype, pop_mito$Population)[pop_order]
stopifnot("every aland_excluded population must have a resolved mitotype" = !anyNA(mito_vec))

contrast <- ifelse(mito_vec == "Faquilonia", 1L, -1L)
cat("Mitotype contrast (1 = Faquilonia-like, -1 = Fpolyctena-like):\n")
print(data.table(Population = pop_order, Mitotype = mito_vec, contrast = contrast))
cat("\nGroup sizes: Faquilonia-like =", sum(contrast == 1),
    " populations, Fpolyctena-like =", sum(contrast == -1), " populations\n")

## ---- verify pop_order matches the already-written u_DIEM.size -------------
src_size <- scan(file.path(OMEGA_DIR, "u_DIEM.size"), what = integer(), quiet = TRUE)
stopifnot(
  "pop_order here disagrees with aland_excluded/u_DIEM.size -- population order mismatch" =
    length(src_size) == length(pop_order)
)

## ---- write + stage --------------------------------------------------------
write.table(t(as.matrix(contrast)), file.path(OMEGA_DIR, "u.mito_contrast"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
if (dir.exists(UNIT_DIR)) {
  file.copy(file.path(OMEGA_DIR, "u.mito_contrast"), file.path(UNIT_DIR, "u.mito_contrast"),
            overwrite = TRUE)
  cat("\nStaged into", UNIT_DIR, "too.\n")
} else {
  cat("\n", UNIT_DIR, " doesn't exist yet (Stage-1 unit prep still running) -- ",
      "re-run this script (or just copy u.mito_contrast) once it's done.\n", sep = "")
}
cat("Wrote", file.path(OMEGA_DIR, "u.mito_contrast"), "\n")
message("[moduleB-stage1-mito-contrast] done")
