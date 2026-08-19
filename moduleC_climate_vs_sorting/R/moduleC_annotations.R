## =========================================================
## MODULE C (revised) -- build the aligned per-eMLG annotation table
## =========================================================
## Assembles, for the complete 32,840-eMLG universe (has_eMLG == TRUE; primary
## config = Aland excluded, 19 pops, fixed LD-pruned Omega), the annotations the
## genome-wide climate calibration is tested against:
##   * DI          -- per-eMLG (consensus) Diagnostic Index      [moduleA_cluster_sorting.rds]
##   * directional_tau05/06/08 -- Module A unidirectional-sorting status at each
##                    fixation-threshold tau in the sensitivity series; `directional`
##                    is an alias for the PRIMARY (tau06). Only this column varies
##                    across tau -- everything below is tau-independent.    [stamped files]
##   * prop_fixed  -- continuous sorting MAGNITUDE (degree of fixation) [moduleA_cluster_sorting.rds]
##   * uni_score   -- signed sorting ORIENTATION (supplementary)  [moduleA_cluster_sorting.rds]
##   * sort_class (primary tau06), differentiated                 [moduleA_cluster_sorting.rds]
##   * n_loci      -- cluster size (cross-checked clustering vs cluster_sorting) [clustering]
##   * recomb      -- map-based recombination rate (cM/Mb) at the cluster's
##                    representative marker (moduleB_architecture convention)
##   * eBF1, eBF2  -- observed eMLG BayPass BF on PC1 / PC2       [assoc / null obj]
##
## EVERY join is an explicit group_id key join with a one-to-one assertion; the
## table is finally ordered to the BayPass row order (eMLG_group_order.txt) so the
## regenerated null BF rows align positionally to it (that alignment is re-asserted
## at run time in moduleC_null_regen.R). Positional matching is NEVER assumed for
## the sorting/DI object (its row order differs from the association object).
##
## Reads : moduleB_climate_GEA/data/moduleB_eMLG_association.rds
##         moduleA_sorting/data/moduleA_cluster_sorting_tau{05,06,08}.rds  (tau series)
##         module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds  (groups: representative, n_loci)
##         data/Frufa_DTOL_PR.ref_genome.recmap
##         aland_excluded_eMLG/eMLG_group_order.txt
## Writes: moduleC_climate_vs_sorting/data/moduleC_annotations.rds
## Run from the formica_hybrid repo root.

suppressMessages(library(data.table))
source("moduleC_climate_vs_sorting/R/moduleC_stat_functions.R")   # MODULEC_TAU_SERIES, tauC_stamp

OUT <- "moduleC_climate_vs_sorting/data/moduleC_annotations.rds"
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

TAUS   <- MODULEC_TAU_SERIES          # c(0.5, 0.6, 0.8)
TSTAMP <- tauC_stamp(TAUS)            # "tau05" "tau06" "tau08"
PRIMARY_STAMP <- tauC_stamp(MODULEC_TAU_PRIMARY)   # "tau06"

## ---- inputs -------------------------------------------------------------
assoc <- readRDS("moduleB_climate_GEA/data/moduleB_eMLG_association.rds")   # group_id, eBF1, eBF2, ...
grp   <- readLines("aland_excluded_eMLG/eMLG_group_order.txt")         # BayPass row order
groups <- readRDS("module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds")$groups
he     <- groups[has_eMLG == TRUE, .(group_id, representative, n_loci)]
## canonical best-SNP: read recomb at each cluster's best_marker (its real
## representative SNP), consistent with DI now coming from best_marker
## (moduleA_cluster_sorting.rds), rather than the centrality representative.
rep_snp <- as.data.table(readRDS("module0_ld_pruning/data/eMLG_5loci_0025_cM05_bestsnp.rds")$rep_snp_all)
he[rep_snp, on = "group_id", best_marker := i.rep_snp]
stopifnot("best_marker missing for some has_eMLG cluster" = all(!is.na(he$best_marker)))

## ---- Module A sorting, tau series ---------------------------------------
## Read the stamped cluster-sorting objects; the DIRECTIONAL classification is the
## only tau-dependent column (verified). Build `cl` = tau-independent columns from
## the primary object + one directional_<tau> column per threshold, asserting the
## shared columns are byte-identical across the series (guards against a stale or
## mismatched stamped file silently entering the join).
cl_tau <- lapply(TSTAMP, function(ts)
  readRDS(sprintf("moduleA_sorting/data/moduleA_cluster_sorting_%s.rds", ts)))
names(cl_tau) <- TSTAMP
SHARED <- c("group_id", "DI", "prop_fixed", "uni_score", "differentiated", "n_loci")
for (ts in TSTAMP[-1]) for (col in SHARED)
  if (!identical(cl_tau[[ts]][[col]], cl_tau[[PRIMARY_STAMP]][[col]]))
    stop(sprintf("tau-independent column '%s' differs between %s and primary %s", col, ts, PRIMARY_STAMP))
cl <- copy(cl_tau[[PRIMARY_STAMP]][, c(SHARED, "sort_class"), with = FALSE])
for (i in seq_along(TSTAMP)) {
  d <- cl_tau[[TSTAMP[i]]]
  stopifnot("directional not 0/1 in stamped file" = all(d$directional %in% c(0L, 1L)),
            "stamped file group_id order differs from primary" = identical(d$group_id, cl$group_id))
  set(cl, j = paste0("directional_", TSTAMP[i]), value = d$directional)
}
set(cl, j = "directional", value = cl[[paste0("directional_", PRIMARY_STAMP)]])   # alias == primary tau

## cross-check cluster size between the clustering object and cluster_sorting (must agree)
nl_chk <- he[cl[, .(group_id, n_loci_cl = n_loci)], on = "group_id"]
stopifnot("n_loci disagrees between clustering and cluster_sorting" =
            all(nl_chk$n_loci == nl_chk$n_loci_cl))

N <- length(grp)
message("[ann] BayPass row order: ", N, " eMLGs")

## ---- ordering / membership assertions (fail loudly) ---------------------
stopifnot(
  "association order != BayPass order"     = identical(as.character(assoc$group_id), grp),
  "duplicate group_id in association"      = !any(duplicated(assoc$group_id)),
  "duplicate group_id in cluster_sorting"            = !any(duplicated(cl$group_id)),
  "duplicate group_id in clustering"       = !any(duplicated(he$group_id)),
  "cluster_sorting missing some eMLGs"               = all(grp %in% cl$group_id),
  "clustering missing some eMLGs"          = all(grp %in% he$group_id),
  "cluster_sorting has extra ids beyond universe"    = setequal(cl$group_id, grp),
  "clustering has_eMLG != universe"        = setequal(he$group_id, grp)
)

## ---- recombination rate at the representative marker (cM/Mb) -------------
rec <- fread("data/Frufa_DTOL_PR.ref_genome.recmap")
stopifnot(ncol(rec) >= 4)
setnames(rec, 1:4, c("chr", "pos", "cM", "cMMb"))
rec[, Chr := sub("chromosome_", "Chr", chr)]
he[, `:=`(rep_chr = sub(":.*", "", best_marker),
          rep_pos = as.integer(sub(".*:", "", best_marker)),
          recomb  = NA_real_)]
for (ch in unique(he$rep_chr)) {
  r <- rec[Chr == ch]
  if (nrow(r) > 1) {
    idx <- which(he$rep_chr == ch)
    he[idx, recomb := approx(r$pos, r$cMMb, xout = he$rep_pos[idx], rule = 2)$y]
  }
}
stopifnot("recomb has NAs" = all(is.finite(he$recomb)))

## ---- assemble by explicit group_id joins, then order to BayPass rows -----
DIR_COLS <- paste0("directional_", TSTAMP)                          # directional_tau05/06/08
ann <- data.table(group_id = grp)                                   # canonical order
ann <- assoc[, .(group_id, eBF1, eBF2, Chr, Pos, size)][ann, on = "group_id"]
ann <- cl[, c("group_id", "DI", "directional", DIR_COLS, "prop_fixed", "uni_score",
              "sort_class", "differentiated", "n_loci"), with = FALSE][ann, on = "group_id"]
ann <- he[, .(group_id, recomb)][ann, on = "group_id"]
setcolorder(ann, "group_id")
ann <- ann[match(grp, group_id)]                                    # enforce BayPass order

## final integrity: order + completeness + one-to-one
## prop_fixed (sorting magnitude) and uni_score (orientation) share an established
## set of 66 missing values (the same eMLGs); assert that, and that all finite
## prop_fixed lie in [0, 1]. Correlations use these fixed complete-case subsets.
pf_na <- is.na(ann$prop_fixed); us_na <- is.na(ann$uni_score)
stopifnot(
  "final order != BayPass order"       = identical(ann$group_id, grp),
  "row count changed"                  = nrow(ann) == N,
  "NA in DI"                           = all(is.finite(ann$DI)),
  "NA in recomb"                       = all(is.finite(ann$recomb)),
  "NA in directional"                  = all(!is.na(ann$directional)),
  "directional not 0/1"                = all(ann$directional %in% c(0L, 1L)),
  "NA in a directional_tau column"     = all(vapply(DIR_COLS, function(c) all(!is.na(ann[[c]])), logical(1))),
  "directional_tau not 0/1"            = all(vapply(DIR_COLS, function(c) all(ann[[c]] %in% c(0L, 1L)), logical(1))),
  "directional alias != primary tau"   = identical(ann$directional, ann[[paste0("directional_", PRIMARY_STAMP)]]),
  "NA in differentiated"               = all(!is.na(ann$differentiated)),
  "NA in eBF1/eBF2"                    = all(is.finite(ann$eBF1) & is.finite(ann$eBF2)),
  "prop_fixed/uni_score NA sets differ"= identical(pf_na, us_na),
  "finite prop_fixed out of [0,1]"     = all(ann$prop_fixed[!pf_na] >= 0 & ann$prop_fixed[!pf_na] <= 1)
)
message(sprintf("[ann] prop_fixed (magnitude): %d NA (established), finite range %.2f..%.2f",
                sum(pf_na), min(ann$prop_fixed, na.rm = TRUE), max(ann$prop_fixed, na.rm = TRUE)))

n_dir_by_tau <- setNames(vapply(DIR_COLS, function(c) sum(ann[[c]] == 1L), integer(1)), TSTAMP)
attr(ann, "meta") <- list(
  config = "aland_excluded / withOmega / 19 pops / fixed LD-pruned Omega",
  N_eMLG = N,
  tau_series = TAUS, tau_primary = MODULEC_TAU_PRIMARY, tau_stamp = TSTAMP,
  n_directional = sum(ann$directional == 1L),           # primary (tau06)
  n_directional_by_tau = n_dir_by_tau,
  recomb_source = "map-interpolated cM/Mb at best-SNP marker (canonical best-SNP representation)",
  DI_source = "per-eMLG best-SNP DI (moduleA_cluster_sorting.rds)",
  magnitude_source = "prop_fixed = degree of fixation (moduleA_cluster_sorting.rds); NOT uni_score",
  orientation_source = "uni_score = signed direction (moduleA_cluster_sorting.rds); supplementary only",
  n_prop_fixed_NA = sum(is.na(ann$prop_fixed)),
  built = as.character(Sys.time())
)
saveRDS(ann, OUT)

cat(sprintf("\n[ann] wrote %s\n  N=%d  DI %.1f..%.1f  recomb %.2f..%.2f cM/Mb\n",
            OUT, N, min(ann$DI), max(ann$DI), min(ann$recomb), max(ann$recomb)))
cat("  directional by tau: ",
    paste(sprintf("%s=%d (%.1f%%)", TSTAMP, n_dir_by_tau, 100 * n_dir_by_tau / N), collapse = "  "),
    sprintf("  [primary %s]\n", PRIMARY_STAMP))
