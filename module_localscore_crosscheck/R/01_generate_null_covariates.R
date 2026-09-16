## =========================================================
## module_localscore_crosscheck -- structured vs. paired-unstructured
## null covariate generation (design: config/null_covariate_design.md)
## =========================================================
## Fresh implementation (not sourced/copied from the exploratory module or
## from module_manuscript_rho05). Draws:
##   10 structured continuous null covariates   (Omega-eigenvector MVN, z-scored)
##   10 paired unstructured continuous nulls    (within-pair permutation of the
##                                                structured draw's 19 values)
##   10 structured mitoC2 null contrasts        (rank-threshold of the
##                                                structured continuous draws
##                                                to the real 7/12 split)
##   10 paired unstructured mitoC2 null contrasts (same threshold, applied to
##                                                the unstructured draws)
##
## Reads (authoritative inputs, from module_manuscript_rho05):
##   module_manuscript_rho05/baypass_stage1/aland_excluded/omega_mat_omega.out
##   module_manuscript_rho05/baypass_stage1/aland_excluded/u_DIEM.size
##   module_manuscript_rho05/baypass_stage1/aland_excluded/u.mito_contrast
##   data/hybrids_only_maf005.Rdata   (sample_data, for pop_order)
##
## Writes:
##   module_localscore_crosscheck/config/null_covariates/
##     null_structured_continuous.env       (10 x 19, BayPass -efile format)
##     null_unstructured_continuous.env     (10 x 19)
##     null_structured_mitoC2.contrast      (10 x 19, BayPass -contrastfile format)
##     null_unstructured_mitoC2.contrast    (10 x 19)
##   module_localscore_crosscheck/manifests/null_covariates_manifest.rds
##   module_localscore_crosscheck/manifests/null_covariates_manifest.md
##
## Run from the repo root:
##   Rscript module_localscore_crosscheck/R/01_generate_null_covariates.R
## =========================================================

suppressMessages(library(data.table))

OMEGA_DIR  <- "module_manuscript_rho05/baypass_stage1/aland_excluded"
OUT_DIR    <- "module_localscore_crosscheck/config/null_covariates"
MANIF_DIR  <- "module_localscore_crosscheck/manifests"
dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(MANIF_DIR, showWarnings = FALSE, recursive = TRUE)

N_REP <- 10L
STRUCTURED_SEEDS <- 20100L + seq_len(N_REP)   # 20101..20110, prespecified
PERMUTE_SEEDS    <- 20200L + seq_len(N_REP)   # 20201..20210, prespecified
stopifnot(!any(STRUCTURED_SEEDS %in% PERMUTE_SEEDS))

## ---- 1. authoritative population order -------------------------------------
load("data/hybrids_only_maf005.Rdata")   # -> sample_data
sd_ex <- sample_data[Population != "Aland"]
pop_order <- unique(sd_ex$Population)
P <- length(pop_order)
stopifnot("expected 19 aland_excluded populations" = P == 19)

## ---- 2. cross-checks against independently-written authoritative files -----
src_size <- scan(file.path(OMEGA_DIR, "u_DIEM.size"), what = integer(), quiet = TRUE)
stopifnot("u_DIEM.size length must match pop_order" = length(src_size) == P)

real_mito <- scan(file.path(OMEGA_DIR, "u.mito_contrast"), what = integer(), quiet = TRUE)
stopifnot("u.mito_contrast length must match pop_order" = length(real_mito) == P)
n_pos <- sum(real_mito ==  1L)
n_neg <- sum(real_mito == -1L)
stopifnot("real mitoC2 contrast must be a +1/-1 partition" = n_pos + n_neg == P)
message(sprintf("[nullcov] real mitoC2 split: %d populations +1, %d populations -1", n_pos, n_neg))
K_POS <- n_pos   # group size to match when rank-thresholding null contrasts (7)

## ---- 3. load and validate Omega ---------------------------------------------
Omega_raw <- as.matrix(fread(file.path(OMEGA_DIR, "omega_mat_omega.out"), header = FALSE))
stopifnot("Omega must be P x P" = all(dim(Omega_raw) == P))
asym <- max(abs(Omega_raw - t(Omega_raw)))
message(sprintf("[nullcov] Omega max asymmetry before symmetrizing: %.3e", asym))
Omega <- (Omega_raw + t(Omega_raw)) / 2

eg <- eigen(Omega, symmetric = TRUE)
n_neg_eig <- sum(eg$values < 0)
message(sprintf("[nullcov] Omega eigenvalues: %d negative (clamped to 0) of %d", n_neg_eig, P))
vals <- pmax(eg$values, 0)

## ---- 4. draw structured + paired unstructured continuous nulls -------------
structured   <- matrix(NA_real_, nrow = N_REP, ncol = P)
unstructured <- matrix(NA_real_, nrow = N_REP, ncol = P)
perm_record  <- matrix(NA_integer_, nrow = N_REP, ncol = P)

for (i in seq_len(N_REP)) {
  set.seed(STRUCTURED_SEEDS[i])
  raw_i <- eg$vectors %*% (sqrt(vals) * rnorm(P))
  structured[i, ] <- as.numeric(scale(raw_i))

  set.seed(PERMUTE_SEEDS[i])
  perm_i <- sample.int(P)
  perm_record[i, ] <- perm_i
  unstructured[i, ] <- structured[i, perm_i]
}
colnames(structured) <- colnames(unstructured) <- pop_order

## ---- 5. validation: draws ----------------------------------------------------
for (i in seq_len(N_REP)) {
  stopifnot(
    "structured draw not mean-0" = abs(mean(structured[i, ])) < 1e-8,
    "structured draw not sd-1"   = abs(sd(structured[i, ]) - 1) < 1e-8,
    "unstructured draw must be a permutation (same sorted values)" =
      isTRUE(all.equal(sort(unname(structured[i, ])), sort(unname(unstructured[i, ])))),
    "unstructured draw must not equal structured draw in original order" =
      !isTRUE(all.equal(unname(structured[i, ]), unname(unstructured[i, ])))
  )
}
message("[nullcov] continuous draw validation: OK (mean/sd, permutation integrity)")

## ---- 6. derive mitoC2 (contrast) nulls by rank-thresholding -----------------
top_k <- function(x, k) {
  p <- rep(-1L, length(x))
  p[order(x, decreasing = TRUE)[seq_len(k)]] <- 1L
  p
}
mito_structured   <- t(apply(structured,   1, top_k, k = K_POS))
mito_unstructured <- t(apply(unstructured, 1, top_k, k = K_POS))
colnames(mito_structured) <- colnames(mito_unstructured) <- pop_order

for (i in seq_len(N_REP)) {
  stopifnot(
    "structured mitoC2 null must match real group sizes" =
      sum(mito_structured[i, ] == 1L) == K_POS && sum(mito_structured[i, ] == -1L) == (P - K_POS),
    "unstructured mitoC2 null must match real group sizes" =
      sum(mito_unstructured[i, ] == 1L) == K_POS && sum(mito_unstructured[i, ] == -1L) == (P - K_POS)
  )
}
message(sprintf("[nullcov] mitoC2 null validation: OK (%d/%d split, both structured and unstructured)",
                 K_POS, P - K_POS))

## ---- 7. write BayPass-format files -------------------------------------------
write_baypass <- function(mat, path) {
  write.table(mat, path, quote = FALSE, row.names = FALSE, col.names = FALSE)
}
write_baypass(structured,        file.path(OUT_DIR, "null_structured_continuous.env"))
write_baypass(unstructured,      file.path(OUT_DIR, "null_unstructured_continuous.env"))
write_baypass(mito_structured,   file.path(OUT_DIR, "null_structured_mitoC2.contrast"))
write_baypass(mito_unstructured, file.path(OUT_DIR, "null_unstructured_mitoC2.contrast"))

## ---- 8. manifest ---------------------------------------------------------------
manifest <- list(
  generated_at        = Sys.time(),
  omega_source        = file.path(OMEGA_DIR, "omega_mat_omega.out"),
  omega_md5           = tools::md5sum(file.path(OMEGA_DIR, "omega_mat_omega.out")),
  omega_max_asymmetry = asym,
  omega_n_neg_eigenvalues = n_neg_eig,
  pop_order           = pop_order,
  P                   = P,
  real_mito_split     = c(pos = n_pos, neg = n_neg),
  n_replicates        = N_REP,
  structured_seeds    = STRUCTURED_SEEDS,
  permute_seeds       = PERMUTE_SEEDS,
  permutations        = perm_record,
  structured_continuous   = structured,
  unstructured_continuous = unstructured,
  structured_mitoC2       = mito_structured,
  unstructured_mitoC2     = mito_unstructured
)
saveRDS(manifest, file.path(MANIF_DIR, "null_covariates_manifest.rds"))

md_lines <- c(
  "# Null covariate generation manifest",
  sprintf("Generated: %s", format(manifest$generated_at)),
  "",
  sprintf("Omega source: `%s`", manifest$omega_source),
  sprintf("Omega MD5: `%s`", manifest$omega_md5),
  sprintf("Omega max asymmetry (pre-symmetrizing): %.3e", asym),
  sprintf("Omega negative eigenvalues clamped: %d / %d", n_neg_eig, P),
  "",
  sprintf("Population order (P=%d, source: data/hybrids_only_maf005.Rdata, Aland excluded):", P),
  paste0("  ", paste(pop_order, collapse = ", ")),
  "",
  sprintf("Real mitoC2 split (source: u.mito_contrast): %d = +1 (Faquilonia-like), %d = -1 (Fpolyctena-like)",
          n_pos, n_neg),
  "",
  sprintf("Replicates: %d", N_REP),
  sprintf("Structured seeds: %s", paste(STRUCTURED_SEEDS, collapse = ", ")),
  sprintf("Permute seeds: %s", paste(PERMUTE_SEEDS, collapse = ", ")),
  "",
  "All validation checks in R/01_generate_null_covariates.R passed (mean/sd of",
  "structured draws; permutation integrity of unstructured draws; exact",
  "7/12 group sizes for both mitoC2 null variants)."
)
writeLines(md_lines, file.path(MANIF_DIR, "null_covariates_manifest.md"))

message("[nullcov] wrote BayPass-format null covariate/contrast files to ", OUT_DIR)
message("[nullcov] wrote manifest to ", MANIF_DIR)
