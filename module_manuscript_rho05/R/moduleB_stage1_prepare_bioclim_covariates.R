## =========================================================
## module_manuscript_rho05 -- BayPass covariate files: bio6/bio11 (winter temp)
## =========================================================
## Winter-temperature GEA, suggested by a colleague: bio6 (min temperature of
## the coldest month) and bio11 (mean temperature of the coldest quarter),
## same structure as the existing PC1/PC2 climate scan. Scaled variants used
## (bio6_scaled/bio11_scaled) WITH -nocovscaling, matching how PC1/PC2 were
## run (pre-scaled values, BayPass's own scaling disabled) -- avoids mixing
## scaling conventions across covariates in the same analysis.
##
## Population order must match aland_excluded/u_DIEM.size exactly. Name
## reconciliation: bioclimatic_variables.csv uses "Nyrhispera1"/"Nyrhispera2",
## genotype data uses "Nyrhispera74"/"Nyrhispera75" -- confirmed same order
## (1=74, 2=75).
##
## Run from the repo root:  Rscript module_manuscript_rho05/R/moduleB_stage1_prepare_bioclim_covariates.R
## =========================================================

suppressMessages(library(data.table))

OMEGA_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded"
UNIT_DIR  <- "module_manuscript_rho05/baypass_stage1/aland_excluded_S1units"

bc <- fread("data/bioclimatic_variables.csv")
bc[Location == "Nyrhispera1", Location := "Nyrhispera74"]
bc[Location == "Nyrhispera2", Location := "Nyrhispera75"]

load("data/hybrids_only_maf005.Rdata")   # sample_data
sd_ex <- sample_data[Population != "Aland"]
pop_order <- unique(sd_ex$Population)

stopifnot("bioclimatic_variables.csv is missing a population present in aland_excluded" =
            all(pop_order %in% bc$Location))

src_size <- scan(file.path(OMEGA_DIR, "u_DIEM.size"), what = integer(), quiet = TRUE)
stopifnot(length(src_size) == length(pop_order))

bc_ord <- bc[match(pop_order, Location)]
cat("Population order + winter-temperature covariates (aland_excluded order):\n")
print(bc_ord[, .(Location, C2_category, bio6, bio6_scaled, bio11, bio11_scaled)])

for (v in c("bio6_scaled", "bio11_scaled")) {
  vals <- bc_ord[[v]]
  stopifnot(!anyNA(vals))
  fn <- sub("_scaled$", "", v)
  write.table(t(as.matrix(vals)), file.path(OMEGA_DIR, paste0("u.", fn)),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  file.copy(file.path(OMEGA_DIR, paste0("u.", fn)), file.path(UNIT_DIR, paste0("u.", fn)), overwrite = TRUE)
  cat("wrote u.", fn, " (", v, ", n=", length(vals), ")\n", sep = "")
}
message("[moduleB-stage1-bioclim-prepare] done")
