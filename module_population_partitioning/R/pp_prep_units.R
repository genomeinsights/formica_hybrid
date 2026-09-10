## =========================================================================
## module_population_partitioning -- 01: build the population x unit
## ancestry-oriented allele-frequency matrix + per-unit FST for the DI25
## LD-reduced units.
##
## Question this module addresses: do different high-DI LD-reduced units
## partition the 20 hybrid populations differently? I.e. can two units both
## be strongly differentiated (high FST) while distinguishing different
## subsets of populations -- which would mean high marginal differentiation
## does not by itself imply a shared multilocus population partition.
##
## Uses ONLY authoritative module_di25 objects (verified against the current
## module_di25/R/ scripts, not the module README, which lags):
##   module_di25/data/di25_inputs.rds          -- map, GTs_hyb, GTs_par
##                                                 (the frozen DI>-25 panel, 51,612 markers)
##   module_di25/data/di25_clustering_cM5.rds  -- $groups: the from-scratch,
##                                                 DI25-only LD-reduced units (11,052;
##                                                 the clustering module_di25/R/di25_sorting.R
##                                                 itself uses, CM_STAMP = "cM5")
##   module_di25/data/di25_sorting_emlg.rds    -- per-unit sorting stats;
##                                                 unit_marker = rep_snp = a REAL SNP
##                                                 (best-SNP for >2-marker clusters via
##                                                 eMLG_best_snp(fill=FALSE) i.e. strictly
##                                                 OBSERVED calls, no eMLG consensus; the
##                                                 cluster's own representative for 1-2
##                                                 marker clusters); current_map_DI is the
##                                                 LATER full-data map's DI (NOT the
##                                                 ascertainment DI), NA for ~1.2% of units.
##
## Orientation: reconstructs the SAME aquilonia-allele-frequency orientation
## parallelism_stats() computes internally (sign from the parental populations'
## mean frequency difference) but never saves -- see
## moduleA_sorting/R/parallelism_stats.R lines ~254-286.
##
## sort_class here is re-derived from the stored (n_aqu, n_pol, n_obs) at the
## LOCKED Module A convention (tau = 0.6, phi = 0.85 i.e. fix_th = 0.15 already
## baked into n_aqu/n_pol, sort_rule = "binom", alpha = 0.05) via the pipeline's
## own classify_sort() -- no new threshold invented.
##
## Run from the formica_hybrid repo root:
##   Rscript module_population_partitioning/R/pp_prep_units.R
## Writes: module_population_partitioning/data/pp_units_Fmat.rds
##   list(u = <data.table, one row per unit>, Fmat = <20 pops x 11,052 units
##        oriented aquilonia-allele frequency>, hybrid_pops = <character>)
## =========================================================================
suppressMessages({ library(data.table); library(ggplot2) })
devtools::load_all("~/gitlab/LDscnR/", quiet = TRUE)      # ohta_fast_prepare()
source("moduleA_sorting/R/parallelism_stats.R")            # classify_sort()

OUTDIR <- "module_population_partitioning/data"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------------
## 1. reconstruct individuals x GTs_all, pops, hybrid_pops, parent rows
##    (mirrors module_di25/R/di25_sorting.R lines 51-70 exactly)
## ---------------------------------------------------------------------
inp <- readRDS("module_di25/data/di25_inputs.rds"); map <- inp$map
e2  <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd  <- e2$sample_data_with_parents

GTs_all  <- rbind(inp$GTs_hyb, inp$GTs_par)
keep_ind <- rownames(GTs_all) %in% sd$Sample_ID
GTs_all  <- GTs_all[keep_ind, ]
pops     <- sd$Population[match(rownames(GTs_all), sd$Sample_ID)]
aqu_pops <- "aquilonia_parent"; pol_pops <- "polyctena_parent"
hybrid_pops <- setdiff(unique(pops), c(aqu_pops, pol_pops))
parent_rows <- grepl("_parent$", pops)
cat(sprintf("[audit] %d individuals (%d hybrids in %d pops + %d parents), %d DI25 markers\n",
            nrow(GTs_all), sum(!parent_rows), length(hybrid_pops), sum(parent_rows), ncol(GTs_all)))

## ---------------------------------------------------------------------
## 2. authoritative units: cM5 from-scratch clustering + sorting output
## ---------------------------------------------------------------------
g  <- readRDS("module_di25/data/di25_clustering_cM5.rds")$groups
ps <- readRDS("module_di25/data/di25_sorting_emlg.rds")
stopifnot(all(ps$unit_marker %in% colnames(GTs_all)))
u <- g[, .(group_id, Chr, n_loci_g = n_loci, has_eMLG, score)][ps, on = "group_id"]
u[, Pos := as.integer(sub(".*:", "", unit_marker))]
u[, ChrNum := as.integer(sub("Chr", "", Chr))]
setorder(u, ChrNum, Pos)

u[, sort_class := classify_sort(n_aqu, n_pol, n_obs, sort_th = 0.6, sort_rule = "binom", alpha = 0.05)]

cat(sprintf("\n[audit] %d LD-reduced units (%d best-SNP eMLG >2-marker, %d representative 1-2-marker)\n",
            nrow(u), sum(u$is_emlg), sum(!u$is_emlg)))
cat("[audit] sort_class @ tau=0.6 (locked convention):\n"); print(u[, .N, by = sort_class][order(-N)])
cat(sprintf("[audit] current_map_DI missing for %d/%d units (later full-data map lacks a matching marker)\n",
            sum(is.na(u$current_map_DI)), nrow(u)))
cat("[audit] n_loci (SNPs collapsed per unit) distribution:\n"); print(summary(u$n_loci_g))

## ---------------------------------------------------------------------
## 3. population x unit ancestry-oriented allele-frequency matrix
## ---------------------------------------------------------------------
E <- GTs_all[, u$unit_marker, drop = FALSE]; colnames(E) <- u$group_id   # obs-only, real best-SNP genotypes
prep <- ohta_fast_prepare(E, pops = pops)
P <- prep$pop_means / 2                                    # pop x unit ALLELE FREQ (dosage/2), 0/1/2-coded input

f_aqu_par <- P[aqu_pops, ]; f_pol_par <- P[pol_pops, ]      # single parental "population" rows
sign_aqu  <- sign(f_aqu_par - f_pol_par)
Ph   <- P[hybrid_pops, , drop = FALSE]
flip <- which(sign_aqu < 0); undef <- which(is.na(sign_aqu) | sign_aqu == 0)
Fmat <- Ph
if (length(flip)) Fmat[, flip]  <- 1 - Ph[, flip]
if (length(undef)) Fmat[, undef] <- NA_real_
## Fmat: 20 hybrid populations x 11,052 units, ORIENTED F. aquilonia-allele frequency in [0,1]
stopifnot(identical(colnames(Fmat), u$group_id))

miss_pop <- rowMeans(is.na(Fmat)); miss_unit <- colMeans(is.na(Fmat))
cat(sprintf("\n[audit] population x unit matrix: %d pops x %d units\n", nrow(Fmat), ncol(Fmat)))
cat("[audit] missingness by population (fraction of units with NA pop-mean):\n")
print(round(sort(miss_pop, decreasing = TRUE), 3))
cat(sprintf("[audit] missingness by unit: median %.3f, mean %.3f, max %.3f (n units 100%% missing: %d)\n",
            median(miss_unit), mean(miss_unit), max(miss_unit), sum(miss_unit == 1)))
cat(sprintf("[audit] units with undefined orientation (parents tied/NA): %d\n", length(undef)))

## ---------------------------------------------------------------------
## 4. per-unit FST (Weir & Cockerham 1984), hybrid populations only
##    (identical estimator to module_di25/R/di25_fst_vs_di.R::wc_ac)
## ---------------------------------------------------------------------
wc_ac <- function(G, pop) {
  levs <- unique(pop); M <- ncol(G)
  N <- P2 <- H <- matrix(0, length(levs), M)
  for (k in seq_along(levs)) {
    gk <- G[pop == levs[k], , drop = FALSE]
    n <- colSums(!is.na(gk)); s <- colSums(gk, na.rm = TRUE); het <- colSums(gk == 1, na.rm = TRUE)
    N[k, ] <- n; P2[k, ] <- ifelse(n > 0, s / (2 * n), 0); H[k, ] <- ifelse(n > 0, het / n, 0)
  }
  C <- colSums(N); sumN2 <- colSums(N^2); r <- colSums(N > 0)
  nbar <- C / r; nc <- (C - sumN2 / C) / (r - 1)
  pbar <- colSums(N * P2) / C; hbar <- colSums(N * H) / C
  s2  <- colSums(N * sweep(P2, 2, pbar)^2) / ((r - 1) * nbar)
  msp <- pbar * (1 - pbar)
  a  <- (nbar / nc) * (s2 - (1 / (nbar - 1)) * (msp - ((r - 1) / r) * s2 - 0.25 * hbar))
  b  <- (nbar / (nbar - 1)) * (msp - ((r - 1) / r) * s2 - ((2 * nbar - 1) / (4 * nbar)) * hbar)
  cc <- 0.5 * hbar
  list(a = a, abc = a + b + cc)
}
ac <- wc_ac(E[!parent_rows, , drop = FALSE], pops[!parent_rows])
u[, FST := ifelse(ac$abc > 0, ac$a / ac$abc, NA_real_)]
cat(sprintf("\n[audit] per-unit FST: median %.3f, range [%.3f, %.3f], NA: %d\n",
            median(u$FST, na.rm = TRUE), min(u$FST, na.rm = TRUE), max(u$FST, na.rm = TRUE), sum(is.na(u$FST))))

saveRDS(list(u = u, Fmat = Fmat, hybrid_pops = hybrid_pops), file.path(OUTDIR, "pp_units_Fmat.rds"))
cat("\n[audit] saved unit table + population x unit matrix ->", file.path(OUTDIR, "pp_units_Fmat.rds"), "\n")
