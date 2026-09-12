## =========================================================
## module_manuscript_rho05 -- AUDIT FIX (Issue 1): standardize BayPass
## observed covariates to match the structured-null convention
## =========================================================
## The observed u.PC1/u.PC2 files are raw (unstandardized: SD 3.267/2.236,
## mean != 0) while every Omega-structured null covariate is drawn via
## scale(...) (mean 0, SD 1), and BayPass is run with -nocovscaling for BOTH
## -- so observed and null Bayes factors are computed under a mismatched
## effective covariate scale (BayPass's default beta-prior grid, disabled
## from auto-rescaling by -nocovscaling, is calibrated for a unit-SD
## covariate). Verified real (not assumed) before writing this fix:
##   u.PC1: n=19, mean=0.190, SD=3.267 ;  u.PC2: n=19, mean=0.283, SD=2.236
##
## bio6/bio11 (u.bio6/u.bio11) are written from the source CSV's
## "_scaled" columns, which turn out to already be close to standardized
## (mean~-0.08, SD~0.96) but NOT computed over exactly this 19-population,
## Aland-excluded subset -- they were evidently scaled over a different/
## larger reference set. Because z-scoring is invariant to any prior
## invertible affine transform (scale(a*x+b) == scale(x) for any x), it does
## not matter whether we start from raw bio6/bio11 or from bio6_scaled/
## bio11_scaled: re-applying scale() to the CURRENT u.bio6/u.bio11 file
## content (already correctly ordered to the 19-pop BayPass order by the
## existing prep script) gives EXACTLY the same z-scores as redoing the
## standardization from the raw temperatures. This lets one code path fix
## all four covariates uniformly and correctly.
##
## RHO05-SPECIFIC FIX, not a change to the shared legacy
## moduleB_climate_GEA/R/moduleB_write_baypass_inputs.R (which also writes
## unstandardized u.PC1/u.PC2 and is used by other, unaudited modules --
## touching it would silently change those modules too; out of scope here,
## flagged in AUDIT_FIXES.md instead). This script POST-PROCESSES the
## already-written covariate files in place.
##
## Omega is NOT touched (frozen, checksum recorded in AUDIT_FIXES.md).
##
## Reads : module_manuscript_rho05/baypass_stage1/aland_excluded/u.{PC1,PC2,bio6,bio11}
##         (source of truth; aland_excluded_S1units/ carries identical copies)
## Writes: the same 4 files in BOTH aland_excluded/ and aland_excluded_S1units/,
##         now standardized (mean 0, SD 1);
##         module_manuscript_rho05/data/moduleB_stage1_covariate_standardization.csv
##         (provenance: Population, variable, raw, standardized, mean, sd, pop_index)
## Pre-fix files are moved (not overwritten in place without a copy) to
## module_manuscript_rho05/stale_pre_fix_20260912/baypass_stage1/<dir>/u.<var>
##
## Run from the repo root:
##   Rscript module_manuscript_rho05/R/moduleB_stage1_standardize_covariates.R
## =========================================================

suppressMessages(library(data.table))

OMEGA_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded"
UNIT_DIR  <- "module_manuscript_rho05/baypass_stage1/aland_excluded_S1units"
STALE     <- "module_manuscript_rho05/stale_pre_fix_20260912"
VARS      <- c("PC1", "PC2", "bio6", "bio11")

## ---- population order (must match u_DIEM.size / every other covariate) ---
load("data/hybrids_only_maf005.Rdata")   # sample_data
sd_ex <- sample_data[Population != "Aland"]
pop_order <- unique(sd_ex$Population)
src_size <- scan(file.path(OMEGA_DIR, "u_DIEM.size"), what = integer(), quiet = TRUE)
stopifnot("pop_order length != u_DIEM.size length" = length(pop_order) == length(src_size),
          "expected exactly 19 Aland-excluded populations" = length(pop_order) == 19)
message("[standardize] population order (n=", length(pop_order), "): ", paste(pop_order, collapse = ", "))

## ---- quarantine pre-fix files (copy, not silently overwrite) -------------
for (d in c(OMEGA_DIR, UNIT_DIR)) {
  dest <- file.path(STALE, d)
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  for (v in VARS) {
    f <- file.path(d, paste0("u.", v))
    if (file.exists(f)) file.copy(f, file.path(dest, paste0("u.", v)), overwrite = TRUE)
  }
}
message("[standardize] quarantined pre-fix covariate files under ", STALE)

## ---- standardize each covariate, build provenance, write back -----------
prov <- list()
for (v in VARS) {
  raw <- scan(file.path(OMEGA_DIR, paste0("u.", v)), quiet = TRUE)
  if (length(raw) != length(pop_order)) stop(sprintf("u.%s: wrong length", v))
  z <- as.numeric(scale(raw))
  if (abs(mean(z)) >= 1e-8)   stop(sprintf("u.%s: standardized mean not ~0", v))
  if (abs(sd(z) - 1) >= 1e-8) stop(sprintf("u.%s: standardized SD not ~1", v))
  if (length(z) != 19)        stop(sprintf("u.%s: length != 19", v))
  prov[[v]] <- data.table(Population = pop_order, pop_index = seq_along(pop_order),
                          variable = v, raw = raw, standardized = z,
                          raw_mean = mean(raw), raw_sd = sd(raw))
  for (d in c(OMEGA_DIR, UNIT_DIR))
    write.table(t(as.matrix(z)), file.path(d, paste0("u.", v)),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  message(sprintf("[standardize] u.%s: raw mean=%.4f sd=%.4f -> standardized mean=%.2e sd=%.6f (n=%d)",
                  v, mean(raw), sd(raw), mean(z), sd(z), length(z)))
}
prov <- rbindlist(prov)
dir.create("module_manuscript_rho05/data", showWarnings = FALSE, recursive = TRUE)
fwrite(prov, "module_manuscript_rho05/data/moduleB_stage1_covariate_standardization.csv")

## ---- Omega: verify frozen (checksum only, not touched) -------------------
om_md5 <- unname(tools::md5sum(file.path(OMEGA_DIR, "omega_mat_omega.out")))
cat(sprintf("\n[standardize] Omega checksum (FROZEN, unchanged by this fix): %s\n", om_md5))
cat("[standardize] wrote module_manuscript_rho05/data/moduleB_stage1_covariate_standardization.csv\n")
cat("[standardize] done -- rerun the 8 observed BayPass scans next (moduleB_stage1_rerun_observed.sh)\n")
