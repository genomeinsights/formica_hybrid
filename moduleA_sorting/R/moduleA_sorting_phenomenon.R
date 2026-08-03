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
##   * di_agg    = "max" -- cluster DI = best DI across all its members
##   * gate on pooled-parental MAF (min_parent_maf), NOT on DI -- DI is kept
##     ungated as a covariate so DI-vs-{recomb, cluster size, PC} stay unbiased
##     and "does DI predict sorting?" is non-circular (Module B, MAF-controlled)
##
## Transparency: every step is in this file. It loads the LDscnR package
## (ohta_fast_prepare) and sources two reviewed stat files -- parallelism_stats.R
## (parallelism_stats), eMLG_parallelism.R (build_sorted_eMLG,
## build_group_consensus, cluster_DI) -- and defines score_eMLG() inline. Run top-to-bottom from the repo root
## (~/gitlab/formica_hybrid): `Rscript R/moduleA_sorting_phenomenon.R` or
## interactively. Writes RDS_data/eMLG_sorted_cM05.rds and moduleA_sorting/data/moduleA_*.rds.

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

## ---- PARAMETERS (edit here) --------------------------------------------
## LOCKED for the merged pipeline -- change only by deliberate decision:
DI_AGG      <- "max"    # cluster DI = max over members (not representative only)

## Open knobs -- defaults chosen, review/override as needed:
SORT_TH     <- 0.5      # unified sort_class threshold (a population fraction)
SORT_RULE   <- "prop_fixed"  # magnitude gate = TOTAL near-fixation (prop_fixed >= sort_th),
                        #   not the larger of |uni_score|/bi_score. Fixes the "valley" that
                        #   demoted both-directions loci to unsorted (bidirectional was
                        #   ~6x under-called under the original "component" rule). See
                        #   parallelism_stats() sort_rule.
MIN_PARENT_MAF <- 0.15  # PRIMARY sorting gate: pooled-parental MAF floor. Keeps loci
                        #   polymorphic in the parents, so a high sorting index isn't just
                        #   a founding near-monomorphism. 0.15 retains ~59% of loci (chosen
                        #   from the A1 diagnostic); revisit that table if you change it.
MIN_DI      <- NULL     # DI is a COVARIATE, not a gate -- keep its full variation
                        #   (gating DI would truncate range & make DI->sorting circular)
FIX_TH      <- 0.1      # per-pop near-fixation tolerance (parallelism_stats default
                        #   is 0.1; small pops of 3-20 diploids are sensitive -- open Q6)
DROP_SIELVA <- FALSE    # exclude Sielva from hybrid_pops (some earlier analyses do;
                        #   there is also a Sielva-excluded BayPass config)
CORES       <- 4        # mclapply cores for consensus building
SCORE_TH    <- 0.80     # score_eMLG floor for the A3 dilution check == the same
                        #   hard floor ld_prune_and_eMLG() enforced when merging, so
                        #   reused clusters (all >= 0.80) never flag and the check
                        #   focuses on the un-vetted built (small) clusters.

CLUSTERING  <- "module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds"  # frozen canonical (input)
SORTED_FILE <- "moduleA_sorting/data/eMLG_sorted_cM05.rds"         # A2 companion (output)

## ---- inputs -------------------------------------------------------------
devtools::load_all("~/gitlab/LDscnR/")             # ohta_fast_prepare() etc. (LDscnR package)
source("moduleA_sorting/R/parallelism_stats.R")    # parallelism_stats(), .binom_two_sided_p()
source("moduleA_sorting/R/eMLG_parallelism.R")     # build_sorted_eMLG(), build_group_consensus(), cluster_DI()

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

## ---- caching helper -----------------------------------------------------
## Read a saved intermediate if it exists (and passes an optional validity
## check), otherwise build it and write it. `valid(obj)` guards against reusing
## a cache produced under different parameters: if it returns non-TRUE the file
## is rebuilt. Set RECOMPUTE <- TRUE to force every cache to rebuild.
RECOMPUTE <- FALSE
cache_rds <- function(path, build, label = basename(path), valid = NULL) {
  if (!RECOMPUTE && file.exists(path)) {
    obj <- readRDS(path)
    if (is.null(valid) || isTRUE(valid(obj))) {
      cat(sprintf("      [cache] read  %s\n", path)); return(obj)
    }
    cat(sprintf("      [cache] %s present but stale -> rebuilding\n", path))
  }
  cat(sprintf("      [build] %s ...\n", label))
  obj <- build(); saveRDS(obj, path)
  cat(sprintf("      [cache] wrote %s\n", path)); obj
}

## ========================================================================
## A1 -- per-SNP parallelism (stable; no eMLG) -> the sorted markers
## ========================================================================
cat("\n[A1] ohta_fast_prepare over", ncol(GTs_wp), "markers x",
    nrow(GTs_wp), "individuals ...\n")
t0 <- Sys.time()
prep_snp <- cache_rds(
  "moduleA_sorting/data/moduleA_prep_snp.rds",
  function() ohta_fast_prepare(GTs_wp, pops = sample_data$Population),
  label = "ohta_fast_prepare (prep_snp)")
cat(sprintf("      prep done | %.0fs\n", elapsed(t0)))

## pooled-parental MAF: allele frequency over ALL parent individuals pooled
## (i.e. pooled parental allele counts), folded to the minor allele. This is
## the PRIMARY sorting gate -- it keeps loci polymorphic in the parents, so
## near-fixation in hybrids is real sorting rather than an allele already
## near-monomorphic in the founding pool. DI is left ungated (a covariate).
parent_freq    <- colMeans(GTs_parents, na.rm = TRUE) / 2
parent_maf_vec <- pmin(parent_freq, 1 - parent_freq)          # marker-named

## ========================================================================
## A1b -- sensitivity of the sorting signal to (sort_th, fix_major)
## ========================================================================
## Diagnostic sweep showing how much of the descriptive per-SNP sorting signal
## depends on the two classification knobs:
##   sort_th -- THE cross-unit comparison threshold (a population fraction).
##   fix_major -- the per-population near-fixation convention, expressed as the
##              MAJOR-allele floor (a pop is near-fixed when the allele it is
##              fixed for reaches >= fix_major). Normally held constant across
##              analyses (see parallelism_stats() docs); swept here ONLY to show
##              the signal is robust to it.
## The ancestry gate (min_parent_maf) is untouched, so the denominator -- the
## number of `differentiated` SNPs -- is CONSTANT across the whole grid and the
## percentages are directly comparable cell-to-cell. prep_snp is reused; only
## the cheap classification is recomputed per (sort_th, fix_major) cell.
suppressPackageStartupMessages({ library(ggplot2); library(wesanderson) })

SORT_TH_GRID   <- c(0.5, 0.6, 0.7, 0.8)
## fixation expressed as the MAJOR-allele floor (more intuitive): a population
## is near-fixed when the allele it is fixed for reaches >= fix_major. This maps
## onto parallelism_stats()'s fix_th (the away-allele tolerance) as
## fix_th = 1 - fix_major, converted at the call site below so the shared,
## reviewed function keeps its locked fix_th convention untouched.
FIX_MAJOR_GRID <- seq(0.80, 0.95, by = 0.05)

sweep_res <- rbindlist(lapply(FIX_MAJOR_GRID, function(fm) {
  rbindlist(lapply(SORT_TH_GRID, function(st) {
    r <- parallelism_stats(
      prep = prep_snp,
      hybrid_pops = hybrid_pops, aqu_pops = aqu_pops, pol_pops = pol_pops,
      DI = DI_vec, min_DI = MIN_DI,
      parent_maf = parent_maf_vec, min_parent_maf = MIN_PARENT_MAF,
      sort_th = st, fix_th = 1 - fm, sort_rule = SORT_RULE   # major-allele floor -> away-allele tol
    )
    d <- r[differentiated == TRUE]
    data.table(
      sort_th = st, fix_major = fm, n_diff = nrow(d),
      aquilonia     = sum(d$sort_class == "aquilonia",na.rm = TRUE),
      polyctena     = sum(d$sort_class == "polyctena",na.rm = TRUE),
      bidirectional = sum(d$sort_class == "bidirectional",na.rm = TRUE),
      sorted        = sum(d$sort_class != "unsorted",na.rm = TRUE)   # any of the three
    )
  }))
}))
sweep_res[, `:=`(
  pct_sorted        = 100 * sorted        / n_diff,
  pct_aquilonia     = 100 * aquilonia     / n_diff,
  pct_polyctena     = 100 * polyctena     / n_diff,
  pct_bidirectional = 100 * bidirectional / n_diff
)]

cat("\n--- A1b: sorting sensitivity sweep -- % of", sweep_res$n_diff[1],
    "differentiated SNPs sorted (rows = fix_major, cols = sort_th) ---\n")
print(dcast(sweep_res, fix_major ~ sort_th, value.var = "pct_sorted"))
saveRDS(sweep_res, "moduleA_sorting/data/moduleA_sortth_fixth_sweep.rds")

 ## long form: one panel for the "any sorted" total + one per sort_class
sweep_long <- melt(
  sweep_res, id.vars = c("sort_th", "fix_major"),
  measure.vars = c("pct_sorted", "pct_aquilonia", "pct_polyctena", "pct_bidirectional"),
  variable.name = "metric", value.name = "pct"
)
sweep_long[, metric := factor(metric,
  levels = c("pct_sorted", "pct_aquilonia", "pct_polyctena", "pct_bidirectional"),
  labels = c("any sorted", "aquilonia", "polyctena", "bidirectional"))]

p_sweep <- ggplot(sweep_long,
                  aes(sort_th, pct, colour = factor(fix_major), group = fix_major)) +
  geom_line() + geom_point() +
  facet_wrap(~ metric, nrow = 1) +
  scale_x_continuous(breaks = SORT_TH_GRID) +
  scale_colour_manual(values = wes_palette("Zissou1", length(FIX_MAJOR_GRID), type = "continuous")) +
  labs(x = "sorting threshold", y = "% of differentiated SNPs",
       colour = "major-allele\nfixation floor",
       title = "Module A -- per-SNP sorting sensitivity to sorting threshold and major-allele fixation floor",
       subtitle = sprintf("min_parent_maf = %.2f; %d differentiated SNPs",
                          MIN_PARENT_MAF, sweep_res$n_diff[1])) +
  theme_bw() +
  theme(strip.background = element_blank())

ggsave("moduleA_sorting/Figures/moduleA_sorting_sweep.png", p_sweep,
       width = 12, height = 3.5, dpi = 300)
cat("Saved: moduleA_sorting/Figures/moduleA_sorting_sweep.png\n")


## ---- real gated pass (primary gate = MIN_PARENT_MAF; DI ungated covariate) ----
cat("\n[A1] gated pass (min_parent_maf =", MIN_PARENT_MAF, ") ...\n")
t0 <- Sys.time()
par_res_snp <- parallelism_stats(
  prep=prep_snp,
  hybrid_pops = hybrid_pops, aqu_pops = aqu_pops, pol_pops = pol_pops,
  DI = DI_vec, min_DI = MIN_DI, parent_maf = parent_maf_vec, min_parent_maf = MIN_PARENT_MAF, sort_th = SORT_TH, fix_th = FIX_TH, sort_rule = SORT_RULE
)

snp_tab <- par_res_snp[differentiated == TRUE, .N, by = sort_class][order(-N)]
snp_tab[, pct := round(100 * N / sum(N), 3)]
cat("\n--- A1: per-SNP sort_class (parent_maf >=", MIN_PARENT_MAF, ") ---\n")
print(snp_tab)

saveRDS(par_res_snp, "moduleA_sorting/data/moduleA_snp.rds")

sorted_markers <- par_res_snp[differentiated == TRUE & sort_class != "unsorted", marker]
cat("\nSorted markers (differentiated, sort_class != unsorted):", length(sorted_markers), "\n")

## ========================================================================
## A1c -- main-figure panel B (DI-decile direction) sensitivity to (sort_th, fix_major)
## ========================================================================
## Main-figure panel B (fig_main_ancestry_sorting.R) shows frac_aqu_of_unidir --
## among UNIDIRECTIONAL loci, the fraction fixing toward aquilonia -- across DI
## deciles, at ONE operating point (sort_th = 0.5, fix_major = 0.90). Here we
## test that the DI-direction reversal is robust to the two classification
## thresholds, on the SAME (sort_th, fix_major) grid used in A1b.
##
## DI deciles are defined ONCE, on the differentiated set at the operating point
## (par_res_snp). That set -- and hence the decile edges -- is threshold-
## independent: `differentiated` depends only on min_parent_maf, and a locus has
## a defined sort_class iff n_obs > 0, which depends only on parental orientation
## (not on sort_th/fix_major). So the bins are identical across every cell and
## the aquilonia fractions are directly comparable. Only the cheap classification
## is recomputed per cell; prep_snp is reused.
uni_cls  <- c("aquilonia", "polyctena")
di_ref   <- par_res_snp[differentiated == TRUE & !is.na(sort_class) & !is.na(DI)]
DI_BREAKS <- quantile(di_ref$DI, 0:10 / 10, na.rm = TRUE)

panelB_sweep <- rbindlist(lapply(FIX_MAJOR_GRID, function(fm) {
  rbindlist(lapply(SORT_TH_GRID, function(st) {
    r <- parallelism_stats(
      prep = prep_snp,
      hybrid_pops = hybrid_pops, aqu_pops = aqu_pops, pol_pops = pol_pops,
      DI = DI_vec, min_DI = MIN_DI,
      parent_maf = parent_maf_vec, min_parent_maf = MIN_PARENT_MAF,
      sort_th = st, fix_th = 1 - fm, sort_rule = SORT_RULE   # major-allele floor -> away-allele tol
    )
    d <- r[differentiated == TRUE & !is.na(sort_class) & !is.na(DI)]
    d[, DI_decile := cut(DI, breaks = DI_BREAKS, include.lowest = TRUE, labels = FALSE)]
    d[, .(
      sort_th = st, fix_major = fm,
      n_unidir           = sum(sort_class %in% uni_cls),
      frac_aqu_of_unidir = sum(sort_class == "aquilonia") /
                             max(1L, sum(sort_class %in% uni_cls))
    ), by = DI_decile][order(DI_decile)]
  }))
}))
saveRDS(panelB_sweep, "moduleA_sorting/data/moduleA_panelB_sweep.rds")

## faceted view: one panel per fix_major, DI decile on x, one line per sort_th.
p_panelB <- ggplot(panelB_sweep,
                   aes(DI_decile, frac_aqu_of_unidir,
                       colour = factor(sort_th), group = sort_th)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey55") +
  geom_line() + geom_point(size = 1.1) +
  facet_wrap(~ fix_major, nrow = 1,
             labeller = as_labeller(
               function(x) paste0("phi == '", formatC(as.numeric(x), format = "f", digits = 2), "'"),
               default = label_parsed)) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = wes_palette("Zissou1", length(SORT_TH_GRID), type = "continuous")) +
  labs(x = "diagnostic-index decile (low -> high)",
       y = "fraction fixing toward aquilonia", colour = expression(tau),
       title = "Panel B sensitivity: DI-decile sorting direction across sorting threshold and major-allele fixation floor",
       subtitle = sprintf("min_parent_maf = %.2f; dashed line = no directional bias (0.5)",
                          MIN_PARENT_MAF)) +
  theme_bw() +
  theme(strip.background = element_blank())

ggsave("moduleA_sorting/Figures/moduleA_panelB_sweep.png", p_panelB,
       width = 12, height = 3.5, dpi = 300)
cat("Saved: moduleA_sorting/Figures/moduleA_panelB_sweep.png\n")


## ---- ancestry-diagnostic gate ----
## Near-fixation is only evidence of SORTING if the locus was segregating
## in the founding admixture, i.e. the parents differ. Loci that are near-
## monomorphic everywhere (low DI / low parent_diff) would otherwise let
## every hybrid pop trivially "near-fix" the same allele -> spurious 20/20
## parallelism. Primary gate is min_DI (matches the DIEM-based pipeline and
## the DI exploration axis); min_parent_diff is an optional mechanistic
## safeguard (default off). Both are applied if both are set.

## ========================================================================
## A2 -- companion eMLG file covering the clusters behind the sorted markers
## ========================================================================
cat("\n[A2] assembling sorted-cluster consensus (reuse frozen columns + build the rest) ...\n")
t0 <- Sys.time()
## cache key = SORTED_FILE; validation rebuilds if the cached file was made from
## a different sorted-marker set (i.e. the sorting thresholds changed).
# sorted_eMLG <- readRDS(SORTED_FILE)
sorted_eMLG <- cache_rds(SORTED_FILE, label = "build_sorted_eMLG (A2)",
  valid = function(se) identical(se$sorting_params$n_sorted_markers, length(sorted_markers)),
  build = function() {
    se <- build_sorted_eMLG(
      groups = groups, eMLG = eMLG, GTs_hybrids = GTs_hyb,
      markers = sorted_markers, cores = CORES, progress = TRUE)
    ## self-documenting provenance: this file depends on the sorting parameters
    se$derived_from   <- CLUSTERING
    se$params         <- clust$params
    se$sorting_params <- list(
      min_parent_maf = MIN_PARENT_MAF, min_DI = MIN_DI, fix_th = FIX_TH,
      sort_th = SORT_TH, drop_sielva = DROP_SIELVA, di_agg = DI_AGG,
      n_sorted_markers = length(sorted_markers))
    se
  })
cat(sprintf("      A2 done | %.0fs\n", elapsed(t0)))
cat("  ", ncol(sorted_eMLG$eMLG), "units:",
    sum(sorted_eMLG$source == "reused"), "reused +",
    sum(sorted_eMLG$source == "built"), "built  ->", SORTED_FILE, "\n")

## ========================================================================
## A3 -- one cluster-level parallelism_stats() on the companion file
## ========================================================================
## Build the matched parent side for exactly the companion file's units, stack
## it under the hybrid-side consensus, and classify once, so every unit is scored
## on the same footing (sort_class + uni_score/bi_score) in a single pass.

units <- colnames(sorted_eMLG$eMLG)
umem  <- setNames(groups[.(units), on = "group_id", members], units)
cat("\n[A3] building matched parent-side consensus for", length(units), "units ...\n")
t0 <- Sys.time()
## cache validation rebuilds if the cached parent side no longer matches the
## current unit set (e.g. sorted_eMLG was rebuilt with different units).
#par_cons <-readRDS("moduleA_sorting/data/eMLG_sorted_cM05_parents.rds")
par_cons <- cache_rds("moduleA_sorting/data/eMLG_sorted_cM05_parents.rds",
  label = "build_group_consensus (A3 parent side)",
  valid = function(pc) identical(colnames(pc), units),
  build = function() build_group_consensus(umem, GTs_hyb, GTs_parents,
                                            cores = CORES, progress = TRUE, label = "A3 parent"))
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
  sort_th = SORT_TH, fix_th = FIX_TH, sort_rule = SORT_RULE
)
cat(sprintf("      cluster-level test done | %.0fs\n", elapsed(t0)))
setnames(cl_res, "marker", "unit_id")
## The cluster-level signal is carried by sort_class and the population-fraction
## scores (uni_score/bi_score). The earlier binomial parallelism test was dropped
## as obsolete -- it added nothing to the story -- so there is no `sig` column.

## attach unit metadata (size + reused/built) for reporting
meta <- sorted_eMLG$groups[, .(unit_id = group_id, n_loci)]
meta[, source := sorted_eMLG$source[unit_id]]
cl_res <- meta[cl_res, on = "unit_id"]

cat("\n--- A3: cluster-level sort_class (one row per independent unit) ---\n")
print(cl_res[differentiated == TRUE, .N, by = sort_class][order(-N)])
cat("\nby unit source:\n")
print(dcast(cl_res[differentiated == TRUE], source ~ sort_class,
            value.var = "unit_id", fun.aggregate = length))

saveRDS(cl_res, "moduleA_sorting/data/moduleA_clusters.rds")

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
# SCORE_TH = 0.8
cat("\n--- A3b: score_eMLG for multi-marker clusters, by sort_class ---\n")
print(dil[n_loci > 1, .(n = .N,
                        median_score = round(median(score, na.rm = TRUE), 3),
                        n_low_score  = sum(score < SCORE_TH, na.rm = TRUE)),
          by = sort_class][order(-n)])
cat("\nWashed-out (unsorted) multi-marker clusters that are ALSO low-fidelity",
    "(score <", SCORE_TH, ") => possible diluted signal, worth a second look:\n")
print(dil[diluted == TRUE][order(score),
          .(unit_id, n_loci, source, n_aqu, n_fixed, score)])

saveRDS(dil, "moduleA_sorting/data/moduleA_dilution.rds")

cat(sprintf("\nModule A complete in %.0fs total. Outputs written:\n", elapsed(t_start)),
    " ", SORTED_FILE, "              (A2 companion eMLG: reused + built units)\n",
    "  moduleA_sorting/data/moduleA_snp.rds       (per-SNP parallelism_stats)\n",
    "  moduleA_sorting/data/moduleA_clusters.rds  (one row per independent unit)\n",
    "  moduleA_sorting/data/moduleA_dilution.rds  (score_eMLG dilution check)\n")
