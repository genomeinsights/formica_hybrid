## =========================================================
## module_manuscript_rho05 -- CORRECTED ancestry/climate/mitotype analysis
## =========================================================
## Supersedes moduleB_ancestry_vs_winter_climate.R,
## moduleB_mitotype_vs_ancestry_climate.R, and
## moduleB_PC1_ancestry_partial_mitotype.R, per a colleague's audit that
## found two real bugs:
##
##  BUG 1 (block permutation): `unit_val <- ...[match(units, unit)]` took the
##  FIRST population's value per shared-origin pair and broadcast it to both
##  members during permutation -- harmless for LangholmenW/R (identical
##  climate values) but WRONG for Bunkkeri/Grundsund (PC1 -3.409 vs -3.899;
##  Grundsund's value was silently discarded). Fixed here by actually
##  averaging ancestry/climate/mitotype within each shared-origin pair to
##  build a genuine 17-observation dataset, and permuting/correlating on
##  THAT reduced dataset -- not permuting-then-broadcasting-back-to-19.
##
##  BUG 2 (partial-correlation Omega-null): null draws replaced PC1 (fixed
##  ancestry+mitotype) instead of replacing ancestry (fixed PC1+mitotype).
##  Omega is a NUCLEAR allele-frequency covariance matrix -- ancestry is the
##  variable mechanistically of that kind (computed directly from allele
##  frequencies at diagnostic markers); PC1/PC2/bio6/bio11 (climate) and
##  mitotype (uniparental, non-recombining) are not. So Omega-null draws
##  should always replace ANCESTRY specifically, holding the other variable
##  fixed. Fixed throughout. For mitotype-vs-climate tests, NEITHER variable
##  is the Omega-modelled quantity, so the Omega-null is dropped there as
##  not well-motivated -- block permutation only.
##
## Omega itself is reduced to match the 17-unit dataset via M %*% Omega %*%
## t(M), M = the same averaging operator used on the phenotypes, rather than
## redrawing null vectors at dimension 17 from an arbitrary matrix -- this
## propagates the ACTUAL estimated relatedness structure through the same
## linear reduction the phenotypes underwent.
##
## Adds Benjamini-Hochberg correction across the 4-variable climate family
## and the 5-variable mitotype family, using the corrected block-permutation
## p-values as primary (consistent across all tests, unlike Omega-null which
## isn't well-motivated for every comparison).
##
## AUDIT FIX (2026-09-12, Issue 3): Omega17 was built in `units`
## (first-appearance) order but dt17 is sorted alphabetically
## (setorder(dt17, unit_id)) -- genuinely different orderings, verified. The
## Omega-null (NullMat17, from Omega17's own eigendecomposition) was compared
## directly against dt17's columns with no name-based alignment, silently
## misaligning which null draw stood in for which lineage unit. Fixed by
## explicitly reordering Omega17 to dt17$unit_id via match(), with a hard
## rownames/colnames assertion. Raw correlations and block-permutation
## results do NOT depend on Omega17's row order and are unaffected (asserted
## below); only the Omega-null p-values change.
##
## Run from the repo root: Rscript module_manuscript_rho05/R/moduleB_ancestry_climate_mitotype_CORRECTED.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2) })

OMEGA_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded"
FIGDIR <- "module_manuscript_rho05/Figures"
NSIM <- 10000; NPERM <- 10000
SEED_SIM <- 2026; SEED_PERM <- 2027

## ---- reuse the already-verified ancestry/PC1/PC2/bio6/bio11/mitotype -----
obj <- readRDS("module_manuscript_rho05/data/moduleB_mitotype_vs_ancestry_climate.rds")
dt19 <- copy(obj$dt19)   # Population, PC1, PC2, ancestry, bio6, bio11, Mitotype, mito01
VARS <- c("PC1", "PC2", "bio6", "bio11")

## =============================================================================
## 1. build a GENUINE 17-lineage-unit dataset (average paired populations)
## =============================================================================
unit <- dt19$Population
unit[dt19$Population %in% c("LangholmenW", "LangholmenR")] <- "unit_Lang"
unit[dt19$Population %in% c("Bunkkeri", "Grundsund")]      <- "unit_BunGru"
units <- unique(unit)
n17 <- length(units)
cat("independent lineage units:", n17, "(collapsed from", nrow(dt19), "populations)\n")

num_cols <- c("PC1", "PC2", "ancestry", "bio6", "bio11", "mito01")
dt17 <- dt19[, lapply(.SD, mean), by = .(unit_id = unit), .SDcols = num_cols]
setorder(dt17, unit_id)
cat("\nBunkkeri/Grundsund pair -- before (2 rows) vs after (1 averaged row):\n")
print(dt19[Population %in% c("Bunkkeri", "Grundsund"), c("Population", num_cols), with = FALSE])
print(dt17[unit_id == "unit_BunGru"])
cat("\nfull 17-unit dataset:\n"); print(dt17)

## ---- reduce Omega to match, via the SAME averaging operator ---------------
Omega19 <- as.matrix(fread(file.path(OMEGA_DIR, "omega_mat_omega.out")))
Omega19 <- (Omega19 + t(Omega19)) / 2
pop_order19 <- dt19$Population
stopifnot(nrow(Omega19) == length(pop_order19))
M <- matrix(0, nrow = n17, ncol = length(pop_order19), dimnames = list(units, pop_order19))
for (u in units) {
  members <- pop_order19[unit == u]
  M[u, members] <- 1 / length(members)
}
Omega17 <- M %*% Omega19 %*% t(M)
Omega17 <- (Omega17 + t(Omega17)) / 2
dimnames(Omega17) <- list(units, units)   # explicit, don't rely on %*% propagation

## AUDIT FIX (Issue 3): Omega17's row/col order here is `units` -- unique()'s
## FIRST-APPEARANCE order in dt19$Population -- while dt17 was built via
## setorder(dt17, unit_id), i.e. ALPHABETICAL order. These are genuinely
## different orderings (verified: identical(units, sort(units)) is FALSE),
## so every downstream Omega-null draw (NullMat17, built from Omega17's own
## eigendecomposition) was being compared directly against dt17's columns
## with NO name-based alignment -- silently misaligning which null value
## stood in for which lineage unit's ancestry. Fixed by explicitly reordering
## Omega17 to dt17's exact row order before building the null, with a hard
## assertion that the names now agree exactly.
ord17 <- match(dt17$unit_id, rownames(Omega17))
stopifnot("AUDIT FIX (Issue 3): some dt17$unit_id not found in Omega17 rownames" = !anyNA(ord17))
Omega17 <- Omega17[ord17, ord17]
stopifnot("AUDIT FIX (Issue 3): Omega17 rownames != dt17$unit_id after reorder" =
            identical(rownames(Omega17), dt17$unit_id),
          "AUDIT FIX (Issue 3): Omega17 colnames != dt17$unit_id after reorder" =
            identical(colnames(Omega17), dt17$unit_id))
message("[AUDIT FIX Issue 3] Omega17 reordered to match dt17$unit_id exactly (was: ",
        paste(units, collapse = ", "), ")")

eg17 <- eigen(Omega17, symmetric = TRUE); vals17 <- pmax(eg17$values, 0); P17 <- nrow(Omega17)
stopifnot(P17 == nrow(dt17))
set.seed(SEED_SIM)
NullMat17 <- sapply(seq_len(NSIM), function(k) as.numeric(scale(eg17$vectors %*% (sqrt(vals17) * rnorm(P17)))))

## =============================================================================
## 2. ancestry vs climate, on the 17-unit dataset: raw r, Omega-null
##    (replaces ANCESTRY), block permutation (on the 17 units directly)
## =============================================================================
res_climate <- rbindlist(lapply(VARS, function(v) {
  obs <- cor(dt17[[v]], dt17$ancestry)
  ## Omega-null: draws replace ancestry, climate variable v fixed
  nc <- as.numeric(cor(NullMat17, as.numeric(scale(dt17[[v]]))))
  k_om <- sum(abs(nc) >= abs(obs))
  p_om <- (1 + k_om) / (NSIM + 1)
  ## block permutation: shuffle v across the 17 units directly (no broadcast bug)
  set.seed(SEED_PERM)
  null_perm <- vapply(seq_len(NPERM), function(i) cor(sample(dt17[[v]]), dt17$ancestry), numeric(1))
  k_bp <- sum(abs(null_perm) >= abs(obs))
  p_bp <- (1 + k_bp) / (NPERM + 1)
  data.table(variable = v, n_units = n17, r_17unit = round(obs, 3),
             p_omega_null = round(p_om, 4), p_block_perm = round(p_bp, 4))
}))
res_climate[, p_block_perm_BH := round(p.adjust(p_block_perm, method = "BH"), 4)]
cat("\n=== CORRECTED: ancestry vs climate, 17 lineage units ===\n"); print(res_climate)

## =============================================================================
## 3. mitotype vs {ancestry, climate}: Omega-null only where ancestry is
##    involved (replaces ancestry); block permutation for all
## =============================================================================
res_mito <- rbindlist(lapply(c("ancestry", VARS), function(v) {
  obs <- cor(dt17$mito01, dt17[[v]])
  if (v == "ancestry") {
    ## ancestry is being tested here; Omega-null draw replaces ANCESTRY,
    ## correlated against FIXED mitotype
    nc <- as.numeric(cor(NullMat17, as.numeric(scale(dt17$mito01))))
    k_om <- sum(abs(nc) >= abs(obs)); p_om <- round((1 + k_om) / (NSIM + 1), 4)
  } else {
    p_om <- NA_real_   # neither mitotype nor a climate variable is Omega's target quantity
  }
  set.seed(SEED_PERM)
  null_perm <- vapply(seq_len(NPERM), function(i) cor(sample(dt17[[v]]), dt17$mito01), numeric(1))
  k_bp <- sum(abs(null_perm) >= abs(obs))
  data.table(variable = v, n_units = n17, r_17unit = round(obs, 3),
             p_omega_null = p_om, p_block_perm = round((1 + k_bp) / (NPERM + 1), 4))
}))
res_mito[, p_block_perm_BH := round(p.adjust(p_block_perm, method = "BH"), 4)]
cat("\n=== CORRECTED: mitotype vs ancestry/climate, 17 lineage units ===\n"); print(res_mito)

## =============================================================================
## 4. PC1-ancestry partial correlation controlling for mitotype, 17 units,
##    Omega-null correctly replacing ANCESTRY (fixed PC1 + mitotype)
## =============================================================================
partial_cor <- function(x, y, z) {
  rxy <- cor(x, y); rxz <- cor(x, z); ryz <- cor(y, z)
  (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
}
raw_r17     <- cor(dt17$PC1, dt17$ancestry)
partial_r17 <- partial_cor(dt17$PC1, dt17$ancestry, dt17$mito01)
cat(sprintf("\nraw cor(PC1, ancestry), 17 units          = %+.3f\n", raw_r17))
cat(sprintf("partial cor(PC1, ancestry | mito), 17 units = %+.3f\n", partial_r17))

## Omega-null: draws replace ANCESTRY (fixed PC1, fixed mitotype)
null_partial_om <- apply(NullMat17, 2, function(a_null) partial_cor(dt17$PC1, a_null, dt17$mito01))
k_om_p <- sum(abs(null_partial_om) >= abs(partial_r17))
p_om_partial <- (1 + k_om_p) / (NSIM + 1)

## block permutation on the partial statistic, 17 units (permute PC1 directly)
set.seed(SEED_PERM)
null_partial_perm <- vapply(seq_len(NPERM), function(i) {
  partial_cor(sample(dt17$PC1), dt17$ancestry, dt17$mito01)
}, numeric(1))
k_bp_p <- sum(abs(null_partial_perm) >= abs(partial_r17))
p_perm_partial <- (1 + k_bp_p) / (NPERM + 1)

cat(sprintf("Omega-null (replaces ancestry) on partial r : p=%.4f\n", p_om_partial))
cat(sprintf("block-permutation on partial r, 17 units    : p=%.4f\n", p_perm_partial))

## =============================================================================
## 5. save + combined printout
## =============================================================================
saveRDS(list(dt17 = dt17, Omega17 = Omega17, res_climate = res_climate, res_mito = res_mito,
             raw_r17 = raw_r17, partial_r17 = partial_r17,
             p_om_partial = p_om_partial, p_perm_partial = p_perm_partial,
             n17 = n17), "module_manuscript_rho05/data/moduleB_ancestry_climate_mitotype_CORRECTED.rds")

cat("\n\n================ FULL CORRECTED SUMMARY ================\n")
cat("-- ancestry vs climate (17 units, BH across 4 tests) --\n"); print(res_climate)
cat("\n-- mitotype vs ancestry/climate (17 units, BH across 5 tests) --\n"); print(res_mito)
cat(sprintf("\n-- PC1-ancestry partial | mitotype (17 units) --\nraw r=%.3f -> partial r=%.3f\nOmega-null p=%.4f  block-perm p=%.4f\n",
            raw_r17, partial_r17, p_om_partial, p_perm_partial))
message("\n[CORRECTED] done")
