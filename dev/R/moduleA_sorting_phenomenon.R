## =========================================================
## MODULE A -- Sorting phenomenon (descriptive)
## =========================================================
## Merged Formica hybrid pipeline, Module A. Establishes the descriptive
## ancestry-sorting signal at two resolutions, on ONE frozen clustering:
##
##   A1  per-SNP parallelism_stats()  -> the sorted markers   (no eMLG needed)
##   A2  build_sorted_eMLG(): a companion eMLG file covering every cluster the
##       sorted markers touch -- REUSING the frozen clustering's consensus
##       columns where they exist, BUILDING the rest fresh (identical method).
##   A3  add a matched parent side, run parallelism_stats() ONCE at cluster
##       level (one BH-FDR over all units), + a score_eMLG() dilution check.
##
## Two eMLG files, by design:
##   data/eMLG_5loci_0025_cM05.rds  -- canonical, FROZEN, geared to the Ohta
##       analyses; the keystone Modules B/C/D join on. Never regenerated here.
##   data/eMLG_sorted_cM05.rds      -- the A2 companion: same $eMLG/$groups
##       shape, extended down to the sorted clusters below min_n_loci_eMLG.
##       Written by THIS script; depends on the sorting parameters below (they
##       are recorded inside it).
##
## Conventions locked for the merged pipeline:
##   * null_prob = 0.5   -- symmetric null ("pooled" is degenerate/circular)
##   * di_agg    = "max" -- cluster DI = best DI across all its members
##   * gate on pooled-parental MAF (min_parent_maf), NOT on DI -- DI is kept
##     ungated as a covariate so DI-vs-{recomb, cluster size, PC} stay unbiased
##     and "does DI predict sorting?" is non-circular (Module B, MAF-controlled)
##   * sort_class is threshold-based and reported SEPARATELY from a binomial
##     significance flag `sig` (q_binom < SIG_Q). Never conflated: sort_class =
##     descriptive magnitude/direction class; sig = departure from the
##     random-direction null.
##
## Transparency: every step is in this file. It sources three reviewed stat
## files -- Ohta.R (ohta_fast_prepare), parallelism_stats.R (parallelism_stats),
## eMLG_parallelism.R (build_sorted_eMLG, build_group_consensus, cluster_DI) --
## and defines score_eMLG() inline. Run top-to-bottom from the repo root
## (~/gitlab/formica_hybrid): `Rscript dev/R/moduleA_sorting_phenomenon.R` or
## interactively. Writes data/eMLG_sorted_cM05.rds and data/moduleA_*.rds.

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

## ---- PARAMETERS (edit here) --------------------------------------------
## LOCKED for the merged pipeline -- change only by deliberate decision:
NULL_PROB   <- 0.5      # symmetric null; "pooled" correlates ~0.96 with outcome
DI_AGG      <- "max"    # cluster DI = max over members (not representative only)

## Open knobs -- defaults chosen, review/override as needed:
SORT_TH     <- 0.5      # unified sort_class threshold (a population fraction)
MIN_PARENT_MAF <- 0.15  # PRIMARY sorting gate: pooled-parental MAF floor. Keeps loci
                        #   polymorphic in the parents, so a high sorting index isn't just
                        #   a founding near-monomorphism. 0.15 retains ~59% of loci (chosen
                        #   from the A1 diagnostic); revisit that table if you change it.
MIN_DI      <- NULL     # DI is a COVARIATE, not a gate -- keep its full variation
                        #   (gating DI would truncate range & make DI->sorting circular)
FIX_TH      <- 0.15     # per-pop near-fixation tolerance (parallelism_stats default
                        #   is 0.1; small pops of 3-20 diploids are sensitive -- open Q6)
MIN_FIXED   <- 5        # min near-fixed pops before a binomial p is computed
DROP_SIELVA <- FALSE    # exclude Sielva from hybrid_pops (some earlier analyses do;
                        #   there is also a Sielva-excluded BayPass config)
SIG_Q       <- 0.05     # q_binom (BH-FDR) cutoff for the SEPARATE `sig` flag
CORES       <- 8        # mclapply cores for consensus building
SCORE_TH    <- 0.80     # score_eMLG floor for the A3 dilution check == the same
                        #   hard floor ld_prune_and_eMLG() enforced when merging, so
                        #   reused clusters (all >= 0.80) never flag and the check
                        #   focuses on the un-vetted built (small) clusters.

CLUSTERING  <- "data/eMLG_5loci_0025_cM05.rds"   # frozen canonical (input)
SORTED_FILE <- "data/eMLG_sorted_cM05.rds"       # A2 companion (output)

## ---- inputs -------------------------------------------------------------
source("dev/R/Ohta.R")                # ohta_fast_prepare()
source("dev/R/parallelism_stats.R")   # parallelism_stats(), .binom_two_sided_p()
source("dev/R/eMLG_parallelism.R")    # build_sorted_eMLG(), build_group_consensus(), cluster_DI()

## hybrids + parents: parents define the aquilonia/polyctena allele orientation
e2 <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
GTs_wp      <- e2$GTs_with_parents            # samples x markers, dosage 0/1/2
sample_data <- e2$sample_data_with_parents    # Sample_ID, Population
map         <- e2$map_hyb_005                 # marker, DiagnosticIndex, ...

## hybrids only: the exact matrix the eMLG consensus was built from
e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1)
GTs_hyb <- e1$GTs_hybrids_005

## canonical clustering (frozen; we only READ it here)
clust  <- readRDS(CLUSTERING)
groups <- clust$groups                        # group_id, Chr, representative, n_loci, members
eMLG   <- clust$eMLG                           # hybrid individuals x group_id (>=5-loci clusters)
stopifnot(!is.null(groups), !is.null(eMLG), !is.null(clust$params))
cat("Clustering:", CLUSTERING, "|", nrow(groups), "groups |",
    ncol(eMLG), "eMLG clusters | cM_threshold =", clust$params$cM_threshold, "\n")

## parents-only matrix (rows = parent individuals) for the matched-flip consensus
parent_ids  <- sample_data[grepl("_parent$", Population), Sample_ID]
GTs_parents <- GTs_wp[parent_ids, , drop = FALSE]

DI_vec   <- setNames(map$DiagnosticIndex, map$marker)
aqu_pops <- "aquilonia_parent"
pol_pops <- "polyctena_parent"
hybrid_pops <- setdiff(unique(sample_data$Population), c(aqu_pops, pol_pops))
if (DROP_SIELVA) hybrid_pops <- setdiff(hybrid_pops, "Sielva")
cat("Hybrid populations tested:", length(hybrid_pops),
    if (DROP_SIELVA) "(Sielva excluded)" else "", "\n")

t_start <- Sys.time()
elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

## ========================================================================
## A1 -- per-SNP parallelism (stable; no eMLG) -> the sorted markers
## ========================================================================
cat("\n[A1] ohta_fast_prepare over", ncol(GTs_wp), "markers x",
    nrow(GTs_wp), "individuals ...\n")
t0 <- Sys.time()
prep_snp <- ohta_fast_prepare(GTs_wp, pops = sample_data$Population)
cat(sprintf("      prep done | %.0fs\n", elapsed(t0)))

## pooled-parental MAF: allele frequency over ALL parent individuals pooled
## (i.e. pooled parental allele counts), folded to the minor allele. This is
## the PRIMARY sorting gate -- it keeps loci polymorphic in the parents, so
## near-fixation in hybrids is real sorting rather than an allele already
## near-monomorphic in the founding pool. DI is left ungated (a covariate).
parent_freq    <- colMeans(GTs_parents, na.rm = TRUE) / 2
parent_maf_vec <- pmin(parent_freq, 1 - parent_freq)          # marker-named

## ---- A1 DIAGNOSTIC: sorting inflation vs parent_maf (to choose MIN_PARENT_MAF) ----
## Run UNGATED (min_parent_maf = 0) so we can see how the sorted fraction behaves
## at low parent_maf and pick the floor where trivial fixation stops inflating it.
cat("[A1] diagnostic pass (ungated) for the parent_maf threshold ...\n")
par_ung <- parallelism_stats(
  prep_snp, hybrid_pops = hybrid_pops, aqu_pops = aqu_pops, pol_pops = pol_pops,
  DI = DI_vec, min_DI = NULL, parent_maf = parent_maf_vec, min_parent_maf = 0,
  min_fixed = MIN_FIXED, sort_th = SORT_TH, fix_th = FIX_TH, null_prob = NULL_PROB
)
diag_dt <- par_ung[!is.na(sort_class) & !is.na(parent_maf)]
diag_dt[, maf_bin := cut(parent_maf, breaks = c(0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5),
                         include.lowest = TRUE)]
maf_diag <- diag_dt[, .(
  n          = .N,
  pct_sorted = round(100 * mean(sort_class != "unsorted"), 1),
  pct_unidir = round(100 * mean(sort_class %in% c("aquilonia", "polyctena")), 1),
  pct_bidir  = round(100 * mean(sort_class == "bidirectional"), 1),
  median_DI  = round(median(DI, na.rm = TRUE), 1)
), by = maf_bin][order(maf_bin)]
cat("\n--- A1 DIAGNOSTIC: sorting vs pooled-parental MAF bin (ungated) ---\n")
cat("    (pick the MAF floor above which pct_sorted stops being inflated)\n")
print(maf_diag)
cat("\nDI distribution (kept ungated -> retains full variation for Module B/C):\n")
print(round(quantile(par_ung$DI, c(0, .1, .25, .5, .75, .9, 1), na.rm = TRUE), 1))
cat("Spearman cor(DI, parent_maf) over classified loci:",
    round(suppressWarnings(cor(diag_dt$DI, diag_dt$parent_maf,
                               method = "spearman", use = "complete.obs")), 3), "\n")

## ---- real gated pass (primary gate = MIN_PARENT_MAF; DI ungated covariate) ----
cat("\n[A1] gated pass (min_parent_maf =", MIN_PARENT_MAF, ") ...\n")
t0 <- Sys.time()
par_res_snp <- parallelism_stats(
  prep_snp,
  hybrid_pops = hybrid_pops, aqu_pops = aqu_pops, pol_pops = pol_pops,
  DI = DI_vec, min_DI = MIN_DI, parent_maf = parent_maf_vec, min_parent_maf = MIN_PARENT_MAF,
  min_fixed = MIN_FIXED, sort_th = SORT_TH, fix_th = FIX_TH, null_prob = NULL_PROB
)
par_res_snp[, sig := !is.na(q_binom) & q_binom < SIG_Q]   # SEPARATE from sort_class
cat(sprintf("      parallelism_stats done | %.0fs\n", elapsed(t0)))

snp_tab <- par_res_snp[differentiated == TRUE, .N, by = sort_class][order(-N)]
snp_tab[, pct := round(100 * N / sum(N), 1)]
cat("\n--- A1: per-SNP sort_class (parent_maf >=", MIN_PARENT_MAF, ") ---\n")
print(snp_tab)
cat("\nSNP-level sig (q_binom <", SIG_Q, ") by sort_class:\n")
print(par_res_snp[differentiated == TRUE, .(n = .N, n_sig = sum(sig)), by = sort_class][order(-n)])

saveRDS(par_res_snp, "data/moduleA_snp.rds")

sorted_markers <- par_res_snp[differentiated == TRUE & sort_class != "unsorted", marker]
cat("\nSorted markers (differentiated, sort_class != unsorted):", length(sorted_markers), "\n")

## ========================================================================
## A2 -- companion eMLG file covering the clusters behind the sorted markers
## ========================================================================
cat("\n[A2] assembling sorted-cluster consensus (reuse frozen columns + build the rest) ...\n")
t0 <- Sys.time()
sorted_eMLG <- build_sorted_eMLG(
  groups = groups, eMLG = eMLG, GTs_hybrids = GTs_hyb,
  markers = sorted_markers, cores = CORES, progress = TRUE
)
cat(sprintf("      A2 done | %.0fs\n", elapsed(t0)))
## self-documenting provenance: this file depends on the sorting parameters
sorted_eMLG$derived_from    <- CLUSTERING
sorted_eMLG$params          <- clust$params
sorted_eMLG$sorting_params  <- list(
  min_parent_maf = MIN_PARENT_MAF, min_DI = MIN_DI, fix_th = FIX_TH,
  sort_th = SORT_TH, null_prob = NULL_PROB, min_fixed = MIN_FIXED,
  drop_sielva = DROP_SIELVA, di_agg = DI_AGG, n_sorted_markers = length(sorted_markers)
)
saveRDS(sorted_eMLG, SORTED_FILE)
cat("  ", ncol(sorted_eMLG$eMLG), "units:",
    sum(sorted_eMLG$source == "reused"), "reused +",
    sum(sorted_eMLG$source == "built"), "built  ->", SORTED_FILE, "\n")

## ========================================================================
## A3 -- one cluster-level parallelism_stats() on the companion file
## ========================================================================
## Build the matched parent side for exactly the companion file's units, stack
## it under the hybrid-side consensus, and test once. Single call => q_binom is
## a coherent BH-FDR across all units with no pooled re-adjustment needed.
units <- colnames(sorted_eMLG$eMLG)
umem  <- setNames(groups[.(units), on = "group_id", members], units)
cat("\n[A3] building matched parent-side consensus for", length(units), "units ...\n")
t0 <- Sys.time()
par_cons  <- build_group_consensus(umem, GTs_hyb, GTs_parents, cores = CORES,
                                   progress = TRUE, label = "A3 parent")
cat(sprintf("      parent consensus done | %.0fs\n", elapsed(t0)))
GTs_units <- rbind(sorted_eMLG$eMLG, par_cons)          # (hybrids + parents) x units
pops_units <- sample_data[match(rownames(GTs_units), Sample_ID), Population]

cat("[A3] cluster-level parallelism_stats over", length(units), "units ...\n")
t0 <- Sys.time()
DI_units         <- cluster_DI(groups, units, DI_vec, di_agg = DI_AGG)
parent_maf_units <- { pf <- colMeans(par_cons, na.rm = TRUE) / 2; pmin(pf, 1 - pf) }
prep_units <- ohta_fast_prepare(GTs_units, pops = pops_units)
cl_res <- parallelism_stats(
  prep_units,
  hybrid_pops = hybrid_pops, aqu_pops = aqu_pops, pol_pops = pol_pops,
  DI = DI_units, min_DI = MIN_DI,
  parent_maf = parent_maf_units, min_parent_maf = MIN_PARENT_MAF,
  min_fixed = MIN_FIXED, sort_th = SORT_TH, fix_th = FIX_TH, null_prob = NULL_PROB
)
cat(sprintf("      cluster-level test done | %.0fs\n", elapsed(t0)))
setnames(cl_res, "marker", "unit_id")
cl_res[, sig := !is.na(q_binom) & q_binom < SIG_Q]

## attach unit metadata (size + reused/built) for reporting
meta <- sorted_eMLG$groups[, .(unit_id = group_id, n_loci)]
meta[, source := sorted_eMLG$source[unit_id]]
cl_res <- meta[cl_res, on = "unit_id"]

cat("\n--- A3: cluster-level sort_class (one row per independent unit) ---\n")
print(cl_res[differentiated == TRUE, .N, by = sort_class][order(-N)])
cat("\nby unit source:\n")
print(dcast(cl_res[differentiated == TRUE], source ~ sort_class,
            value.var = "unit_id", fun.aggregate = length))
cat("\ncluster-level sig (q_binom <", SIG_Q, ") by sort_class:\n")
print(cl_res[differentiated == TRUE, .(n = .N, n_sig = sum(sig)), by = sort_class][order(-n)])

saveRDS(cl_res, "data/moduleA_clusters.rds")

## ========================================================================
## A3b -- score_eMLG() dilution check on the "washed-out" clusters
## ========================================================================
## Multi-marker clusters that individually contained a sorted SNP but test
## "unsorted" once aggregated: genuinely signal-free, or a real signal DILUTED
## by averaging in less-informative members? score_eMLG(x) = cor(round(x),x)^2
## is the round-trip fidelity of a consensus dosage against its hard call:
##   high => coherent block => a trustworthy "no sorting" call
##   low  => heterogeneous members => dilution plausible, worth revisiting
## Defined inline (one expression) for full transparency of the statistic.
score_eMLG <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2L) return(NA_real_)
  suppressWarnings(stats::cor(round(x), x)^2)
}

hc <- sorted_eMLG$eMLG                          # hybrid-side consensus, all units
scores <- data.table(
  unit_id = colnames(hc),
  score   = vapply(seq_len(ncol(hc)), function(j) score_eMLG(hc[, j]), numeric(1))
)
dil <- scores[cl_res, on = "unit_id"]
dil[, diluted := n_loci > 1 & sort_class == "unsorted" & score < SCORE_TH]

cat("\n--- A3b: score_eMLG for multi-marker clusters, by sort_class ---\n")
print(dil[n_loci > 1, .(n = .N,
                        median_score = round(median(score, na.rm = TRUE), 3),
                        n_low_score  = sum(score < SCORE_TH, na.rm = TRUE)),
          by = sort_class][order(-n)])
cat("\nWashed-out (unsorted) multi-marker clusters that are ALSO low-fidelity",
    "(score <", SCORE_TH, ") => possible diluted signal, worth a second look:\n")
print(dil[diluted == TRUE][order(score),
          .(unit_id, n_loci, source, n_aqu, n_fixed, score)])

saveRDS(dil, "data/moduleA_dilution.rds")

cat(sprintf("\nModule A complete in %.0fs total. Outputs written:\n", elapsed(t_start)),
    " ", SORTED_FILE, "         (A2 companion eMLG: reused + built units)\n",
    "  data/moduleA_snp.rds       (per-SNP parallelism_stats + sig)\n",
    "  data/moduleA_clusters.rds  (one row per independent unit + sig)\n",
    "  data/moduleA_dilution.rds  (score_eMLG dilution check)\n")
