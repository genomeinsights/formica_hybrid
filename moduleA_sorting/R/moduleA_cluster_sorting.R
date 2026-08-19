## =========================================================
## MODULE A -- cluster-level sorting annotation (has_eMLG universe)
## =========================================================
## Produces the canonical per-eMLG sorting call over the FULL has_eMLG universe
## (every >=5-loci cluster with a stored consensus, i.e. colnames(clust$eMLG)),
## in the schema the climate modules join on. This is the home of what used to be
## the flat R/moduleC_sorting_climate.R checkpoint ("C3_cl"): it is a SORTING
## product, so it belongs in Module A. moduleC_climate_vs_sorting (and any other
## climate-side reader) now reads it from here, so the climate arm no longer
## depends on deprecated Module-C code.
##
## NB this is DELIBERATELY the full-universe companion to A3 in
## moduleA_sorting_phenomenon.R: A3's moduleA_clusters.rds is the sorted-focused
## companion set (build_sorted_eMLG reuses only the frozen clusters TOUCHED by a
## sorted marker, plus freshly-built small ones), whereas the climate modules need
## EVERY has_eMLG cluster aligned to the BayPass eMLG row order. So this script
## builds its own matched parent-side consensus over all has_eMLG clusters and
## classifies once.
##
## Runs parallelism_stats() ONCE at cluster level over the eMLG consensus + a
## matched parent-side consensus, with Module A's LOCKED conventions:
##   min_parent_maf 0.15, sort_th 0.5, fix_th 0.1, di_agg "max", DI ungated.
## (The old flat checkpoint used fix_th 0.15 by accident -- corrected here to 0.1
## to match the rest of Module A.)
##
## Reads : data/hybrids_and_parents_maf005.Rdata, data/hybrids_only_maf005.Rdata,
##         module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds
## Writes: moduleA_sorting/data/moduleA_cluster_sorting_{tau05,tau06,tau08}.rds
##         (group_id, n_loci, differentiated, sort_class, DI, prop_fixed,
##          uni_score, directional, sorted) -- one per threshold in the sensitivity
##         series; the primary (tau06) is ALSO written to the unstamped name
##         moduleA_cluster_sorting.rds (what Module C reads). Plus
##         moduleA_cluster_sorting_counts.rds (the tau-independent per-eMLG counts,
##         for cheap re-classification at any tau).
## Run from the formica_hybrid repo root.

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")           # ohta_fast_prepare() etc. (LDscnR package)
source("moduleA_sorting/R/parallelism_stats.R")  # parallelism_stats()
source("moduleA_sorting/R/eMLG_parallelism.R")   # build_group_consensus(), cluster_DI()

## ---- LOCKED parameters (match moduleA_sorting_phenomenon.R) --------------
MIN_PARENT_MAF <- 0.15    # PRIMARY sorting gate (pooled-parental MAF floor)
SORT_TH        <- 0.6     # unified sort_class threshold
SORT_RULE      <- "binom"       # prop_fixed magnitude gate + binomial direction test (see parallelism_stats)
ALPHA          <- 0.05          # binom direction-test significance
FIX_TH         <- 0.15    # per-pop near-fixation tolerance (near-fixation floor 0.85)
MIN_DI         <- NULL    # DI is a covariate, never a gate
## DI_AGG (max-over-members) is superseded: under the best-SNP representation a
## cluster's DI is simply its representative SNP's own DI -- no averaging.
CORES          <- 4
CLUSTERING     <- "module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds"
BESTSNP        <- "module0_ld_pruning/data/eMLG_5loci_0025_cM05_bestsnp.rds"
OUT            <- "moduleA_sorting/data/moduleA_cluster_sorting.rds"
uni_cls        <- c("aquilonia", "polyctena")

## ---- inputs -------------------------------------------------------------
e2 <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
GTs_wp <- e2$GTs_with_parents; sample_data <- e2$sample_data_with_parents; map <- e2$map_hyb_005
e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1); GTs_hyb <- e1$GTs_hybrids_005

clust  <- readRDS(CLUSTERING); groups <- clust$groups
best   <- readRDS(BESTSNP)                    # best-SNP companion: $geno, $stats$best_marker
DI_vec <- setNames(map$DiagnosticIndex, map$marker)
aqu_pops <- "aquilonia_parent"; pol_pops <- "polyctena_parent"
hybrid_pops <- setdiff(unique(sample_data$Population), c(aqu_pops, pol_pops))
parent_ids  <- sample_data[grepl("_parent$", Population), Sample_ID]
GTs_parents <- GTs_wp[parent_ids, , drop = FALSE]

## ---- cluster-level parallelism over the FULL has_eMLG universe -----------
## Best-SNP representation: each has_eMLG cluster is ONE real SNP (best_marker),
## not the averaged consensus. Hybrid side = best$geno (that SNP's calls, missing
## filled from consensus); parent side = the SAME SNP's real parental calls;
## cluster DI = that SNP's own DI (no di_agg averaging). One real SNP throughout.
has_ids <- colnames(best$geno)                             # every >=5-loci cluster
bm      <- setNames(best$stats$best_marker, best$stats$group_id)[has_ids]
stopifnot(!anyNA(bm), all(bm %in% colnames(GTs_parents)),
          identical(rownames(best$geno), rownames(GTs_hyb)) |
            all(rownames(best$geno) %in% rownames(GTs_wp)))
message("[A cluster-sorting] classifying ", length(has_ids), " has_eMLG clusters (best-SNP) ...")
t0 <- Sys.time()
hyb_units <- best$geno[, has_ids, drop = FALSE]                       # hybrids x clusters
par_units <- GTs_parents[, bm, drop = FALSE]; colnames(par_units) <- has_ids  # parents x clusters
GTs_units  <- rbind(hyb_units, par_units)
pops_units <- sample_data[match(rownames(GTs_units), Sample_ID), Population]
DI_units   <- setNames(DI_vec[bm], has_ids)                # single-SNP DI, not max over members
maf_units  <- { pf <- colMeans(par_units, na.rm = TRUE) / 2; pmin(pf, 1 - pf) }
prep_units <- ohta_fast_prepare(GTs_units, pops = pops_units)
ps <- parallelism_stats(prep_units, hybrid_pops = hybrid_pops, aqu_pops = aqu_pops,
                        pol_pops = pol_pops, DI = DI_units, min_DI = MIN_DI,
                        parent_maf = maf_units, min_parent_maf = MIN_PARENT_MAF,
                        sort_th = SORT_TH, fix_th = FIX_TH, sort_rule = SORT_RULE)
setnames(ps, "marker", "group_id")
message(sprintf("      done | %.0fs", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

## The per-eMLG counts are tau-INDEPENDENT: keep them once, then classify at each
## tau in the predefined series (MODULEA_TAU_SERIES). Saving the counts lets any
## future tau be re-classified cheaply (classify_sort) without this expensive build.
base <- groups[.(has_ids), on = "group_id", .(group_id, n_loci)][
  ps[, .(group_id, differentiated, n_aqu, n_pol, n_obs, prop_fixed, uni_score, p_binom, DI)],
  on = "group_id"]
saveRDS(base, "moduleA_sorting/data/moduleA_cluster_sorting_counts.rds")

emit <- function(tau) {
  cl <- copy(base)[, sort_class := NA_character_]
  ok <- cl$differentiated & cl$n_obs > 0 & !is.na(cl$uni_score)
  cl[ok, sort_class := classify_sort(n_aqu, n_pol, n_obs, sort_th = tau,
                                     sort_rule = SORT_RULE, alpha = ALPHA)]
  cl[, `:=`(directional = as.integer(sort_class %in% uni_cls),
            sorted      = as.integer(sort_class %in% c(uni_cls, "unresolved")))]
  out <- cl[, .(group_id, n_loci, differentiated, sort_class, DI, prop_fixed, uni_score,
                directional, sorted)]
  saveRDS(out, sprintf("moduleA_sorting/data/moduleA_cluster_sorting_%s.rds", tau_stamp(tau)))
  if (isTRUE(all.equal(tau, MODULEA_TAU_PRIMARY))) saveRDS(out, OUT)  # primary -> unstamped (Module C)
  cat(sprintf("  %s | %5d directional | %4d unresolved | %5d differentiated\n",
              tau_stamp(tau), sum(out$directional == 1L, na.rm = TRUE),
              sum(out$sort_class == "unresolved", na.rm = TRUE),
              sum(out$differentiated, na.rm = TRUE)))
  invisible(out)
}
cat(sprintf("\n[A cluster-sorting] has_eMLG universe (%d clusters), tau series:\n", nrow(base)))
invisible(lapply(MODULEA_TAU_SERIES, emit))
cat(sprintf("  primary %s also written to %s\n", tau_stamp(MODULEA_TAU_PRIMARY), OUT))
