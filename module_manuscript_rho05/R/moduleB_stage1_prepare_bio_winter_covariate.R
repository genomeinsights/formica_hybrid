## =========================================================
## AUDIT FIX (Issue 7): reduce bio6/bio11 to ONE winter-temperature axis
## =========================================================
## bio6 and bio11 are essentially the same signal (raw covariate r=0.96,
## observed BF r=0.94 -- checked before this fix was written). DECISION
## (pipeline author, 2026-09-12): reduce to one axis BEFORE null calibration
## -- their mean (equivalent to PC1 up to scale, since both are already
## standardized to mean 0/SD 1 with near-equal variance) -- rather than
## calibrating both as if they were independent tests, which risks reporting
## whichever of two near-duplicate variables looks stronger after the fact.
##
## Built from the ALREADY-STANDARDIZED u.bio6/u.bio11 (Issue 1 fix), then
## re-standardized itself so it carries the same mean-0/SD-1 convention as
## every other covariate (observed and null alike).
##
## The separate bio6/bio11 full-SNP scans (queued before this decision was
## made) are kept for the Manhattan-plot visualization only, captioned as two
## highly correlated components of one axis -- NOT as two independent
## calibrated tests. Only bio_winter gets its own floor-survivor calibration.
##
## Run from the repo root:
##   Rscript module_manuscript_rho05/R/moduleB_stage1_prepare_bio_winter_covariate.R
## =========================================================

suppressMessages(library(data.table))
OMEGA_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded"
UNIT_DIR  <- "module_manuscript_rho05/baypass_stage1/aland_excluded_S1units"

z6  <- scan(file.path(OMEGA_DIR, "u.bio6"),  quiet = TRUE)
z11 <- scan(file.path(OMEGA_DIR, "u.bio11"), quiet = TRUE)
stopifnot(length(z6) == 19, length(z11) == 19,
          abs(mean(z6)) < 1e-6, abs(sd(z6) - 1) < 1e-6,
          abs(mean(z11)) < 1e-6, abs(sd(z11) - 1) < 1e-6)
r <- cor(z6, z11)
message(sprintf("[bio_winter] bio6 vs bio11 correlation (standardized): r=%.4f", r))

combined <- rowMeans(cbind(z6, z11))
z_winter <- as.numeric(scale(combined))
stopifnot(abs(mean(z_winter)) < 1e-8, abs(sd(z_winter) - 1) < 1e-8, length(z_winter) == 19)

for (d in c(OMEGA_DIR, UNIT_DIR))
  write.table(t(as.matrix(z_winter)), file.path(d, "u.bio_winter"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)

message(sprintf("[bio_winter] wrote u.bio_winter (mean=%.2e sd=%.6f, n=%d) to both dirs",
                mean(z_winter), sd(z_winter), length(z_winter)))
cat("[bio_winter] done -- rerun the observed Stage-1-unit BayPass scan for bio_winter next\n")
