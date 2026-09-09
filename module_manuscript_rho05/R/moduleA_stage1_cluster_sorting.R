## =========================================================
## module_manuscript_rho05 -- Module A equivalent for Stage-1 units
## =========================================================
## Adapted from moduleA_sorting/R/moduleA_cluster_sorting.R for the
## Stage-1-direct universe (18,361 units, n_snps>=5, best-SNP represented --
## module_manuscript_rho05/data/moduleB_stage1_units_bestsnp.rds, already
## built by moduleB_stage1_prepare_baypass_inputs.R). Needed because
## moduleC_annotations.R (Fig [climate_ass]) joins the BayPass association
## results against Module A's per-eMLG sorting classification (DI,
## prop_fixed, uni_score, sort_class, directional at tau in
## MODULEA_TAU_SERIES) -- a transitive dependency, not something Fig
## [climate_ass] can skip.
##
## LOCKED parameters, unchanged from the canonical Module A build:
## min_parent_maf 0.15, sort_th 0.6 (parallelism_stats' own internal call;
## superseded downstream by the per-tau classify_sort reclassification),
## fix_th 0.15, sort_rule "binom", alpha 0.05, DI ungated. Best-SNP
## representation throughout: hybrid genotypes from best$geno (real SNP
## calls, missing filled from consensus), parent genotypes at the SAME
## best_marker, DI from that SNP directly (no di_agg averaging).
##
## Reads : data/hybrids_and_parents_maf005.Rdata, data/hybrids_only_maf005.Rdata,
##         module_manuscript_rho05/data/moduleB_stage1_units_bestsnp.rds
## Writes: module_manuscript_rho05/data/moduleA_stage1_cluster_sorting_{tau05,tau06,tau08}.rds
##         (group_id, n_loci, differentiated, sort_class, DI, prop_fixed,
##          uni_score, directional, sorted); primary (tau06) also to the
##          unstamped moduleA_stage1_cluster_sorting.rds (what moduleC reads);
##          plus moduleA_stage1_cluster_sorting_counts.rds (tau-independent).
##
## Run from the repo root: Rscript module_manuscript_rho05/R/moduleA_stage1_cluster_sorting.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")
source("moduleA_sorting/R/parallelism_stats.R")   # parallelism_stats(), classify_sort(), tau_stamp(),
                                                   # MODULEA_TAU_SERIES, MODULEA_TAU_PRIMARY

MIN_PARENT_MAF <- 0.15
SORT_TH        <- 0.6
SORT_RULE      <- "binom"
ALPHA          <- 0.05
FIX_TH         <- 0.15
MIN_DI         <- NULL
CORES          <- 4
BESTSNP_OBJ    <- "module_manuscript_rho05/data/moduleB_stage1_units_bestsnp.rds"
OUTDIR         <- "module_manuscript_rho05/data"
OUT            <- file.path(OUTDIR, "moduleA_stage1_cluster_sorting.rds")
uni_cls        <- c("aquilonia", "polyctena")

## ---- inputs ---------------------------------------------------------------
e2 <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
GTs_wp <- e2$GTs_with_parents; sample_data <- e2$sample_data_with_parents; map <- e2$map_hyb_005
e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1); GTs_hyb <- e1$GTs_hybrids_005

obj <- readRDS(BESTSNP_OBJ)
groups <- obj$groups   # group_id, Chr, representative, n_loci, score, has_eMLG, members
best   <- obj$best     # stats (group_id, best_marker, ...), geno (individuals x units)
DI_vec <- setNames(map$DiagnosticIndex, map$marker)
aqu_pops <- "aquilonia_parent"; pol_pops <- "polyctena_parent"
hybrid_pops <- setdiff(unique(sample_data$Population), c(aqu_pops, pol_pops))
parent_ids  <- sample_data[grepl("_parent$", Population), Sample_ID]
GTs_parents <- GTs_wp[parent_ids, , drop = FALSE]

## ---- cluster-level parallelism over the Stage-1 (n_snps>=5) universe ------
has_ids <- colnames(best$geno)                              # all 18,361 Stage-1 units
bm      <- setNames(best$stats$best_marker, best$stats$group_id)[has_ids]
stopifnot(!anyNA(bm), all(bm %in% colnames(GTs_parents)),
          identical(rownames(best$geno), rownames(GTs_hyb)) |
            all(rownames(best$geno) %in% rownames(GTs_wp)))
message("[A-stage1 cluster-sorting] classifying ", length(has_ids), " Stage-1 units (best-SNP) ...")
t0 <- Sys.time()
hyb_units <- best$geno[, has_ids, drop = FALSE]
par_units <- GTs_parents[, bm, drop = FALSE]; colnames(par_units) <- has_ids
GTs_units  <- rbind(hyb_units, par_units)
pops_units <- sample_data[match(rownames(GTs_units), Sample_ID), Population]
DI_units   <- setNames(DI_vec[bm], has_ids)
maf_units  <- { pf <- colMeans(par_units, na.rm = TRUE) / 2; pmin(pf, 1 - pf) }
prep_units <- ohta_fast_prepare(GTs_units, pops = pops_units)
ps <- parallelism_stats(prep_units, hybrid_pops = hybrid_pops, aqu_pops = aqu_pops,
                        pol_pops = pol_pops, DI = DI_units, min_DI = MIN_DI,
                        parent_maf = maf_units, min_parent_maf = MIN_PARENT_MAF,
                        sort_th = SORT_TH, fix_th = FIX_TH, sort_rule = SORT_RULE)
setnames(ps, "marker", "group_id")
message(sprintf("      done | %.0fs", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

base <- groups[group_id %in% has_ids, .(group_id, n_loci)][
  ps[, .(group_id, differentiated, n_aqu, n_pol, n_obs, prop_fixed, uni_score, p_binom, DI)],
  on = "group_id"]
saveRDS(base, file.path(OUTDIR, "moduleA_stage1_cluster_sorting_counts.rds"))

emit <- function(tau) {
  cl <- copy(base)[, sort_class := NA_character_]
  ok <- cl$differentiated & cl$n_obs > 0 & !is.na(cl$uni_score)
  cl[ok, sort_class := classify_sort(n_aqu, n_pol, n_obs, sort_th = tau,
                                     sort_rule = SORT_RULE, alpha = ALPHA)]
  cl[, `:=`(directional = as.integer(sort_class %in% uni_cls),
            sorted      = as.integer(sort_class %in% c(uni_cls, "unresolved")))]
  out <- cl[, .(group_id, n_loci, differentiated, sort_class, DI, prop_fixed, uni_score,
                directional, sorted)]
  saveRDS(out, file.path(OUTDIR, sprintf("moduleA_stage1_cluster_sorting_%s.rds", tau_stamp(tau))))
  if (isTRUE(all.equal(tau, MODULEA_TAU_PRIMARY))) saveRDS(out, OUT)
  cat(sprintf("  %s | %5d directional | %4d unresolved | %5d differentiated\n",
              tau_stamp(tau), sum(out$directional == 1L, na.rm = TRUE),
              sum(out$sort_class == "unresolved", na.rm = TRUE),
              sum(out$differentiated, na.rm = TRUE)))
  invisible(out)
}
cat(sprintf("\n[A-stage1 cluster-sorting] Stage-1 (n_snps>=5) universe (%d clusters), tau series:\n", nrow(base)))
invisible(lapply(MODULEA_TAU_SERIES, emit))
cat(sprintf("  primary %s also written to %s\n", tau_stamp(MODULEA_TAU_PRIMARY), OUT))
