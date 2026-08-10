## =========================================================
## module_di25 (high-DI analyses) -- ancestry sorting on the DI25 data (per SNP + per eMLG)
## =========================================================
## Estimates ancestry sorting exactly as Module A does, but on the high-DI-only
## data set: per SNP (51,612 diagnostic markers) and per eMLG unit (the from-
## scratch 5 cM clustering; eMLG consensus for >2-marker clusters, representative
## SNP for 1-2). Near-fixation floor phi = 0.85 (fix_th = 0.15) held fixed; the
## sorting threshold tau is swept {0.5, 0.6, 0.7, 0.8} as in Module A (Fig S1).
##
## Method (Module A conventions, sort_rule = "binom", alpha = 0.05):
##   * pooled-parental MAF gate >= 0.15 (min_parent_maf); DI left ungated.
##   * parallelism_stats() run ONCE per level -> tau-independent per-unit counts
##     (n_aqu, n_pol, n_obs); classify_sort() re-classifies at each tau.
##   * A unit is `sorted` when prop_fixed >= tau; among sorted, direction is
##     resolved by a two-sided binomial test of n_aqu among n_fixed (null 0.5) --
##     else `unresolved`, or `ambiguous` when too few populations are fixed to
##     ever reach significance.
##
## Individuals: the 164 hybrids + 30 parents present in sample_data (the one
## extra tsv hybrid, Hei159_h, has no population label and is dropped). 20 hybrid
## populations tested; parents supply the aquilonia/polyctena orientation.
##
## Run from the repo root:  Rscript module_di25/R/di25_sorting.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")             # ohta_fast_prepare(), consensus_dosage()
source("moduleA_sorting/R/parallelism_stats.R")    # parallelism_stats(), classify_sort()

## ---- parameters ---------------------------------------------------------
CM_STAMP  <- "cM5"
CLUST     <- sprintf("module_di25/data/di25_clustering_%s.rds", CM_STAMP)
INPUTS    <- "module_di25/data/di25_inputs.rds"
OUTDIR    <- "module_di25/data"
TAU_GRID  <- c(0.5, 0.6, 0.7, 0.8)
FIX_TH    <- 0.15      # phi = 0.85 near-fixation floor
MIN_PMAF  <- 0.15      # pooled-parental MAF gate
SORT_RULE <- "binom"
ALPHA     <- 0.05

## ---- inputs -------------------------------------------------------------
inp <- readRDS(INPUTS); map <- inp$map
e2  <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd  <- e2$sample_data_with_parents
DI_vec <- setNames(e2$map_hyb_005$DiagnosticIndex, e2$map_hyb_005$marker)   # marker-named DI

## matched 194-individual matrix (drop the one hybrid absent from sample_data)
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)                    # 195 x markers (012)
keep_ind <- rownames(GTs_all) %in% sd$Sample_ID
GTs_all  <- GTs_all[keep_ind, ]
pops     <- sd$Population[match(rownames(GTs_all), sd$Sample_ID)]
aqu_pops <- "aquilonia_parent"; pol_pops <- "polyctena_parent"
hybrid_pops <- setdiff(unique(pops), c(aqu_pops, pol_pops))
parent_rows <- grepl("_parent$", pops)
cat(sprintf("Individuals: %d (%d hybrids in %d pops + %d parents)\n",
            nrow(GTs_all), sum(!parent_rows), length(hybrid_pops), sum(parent_rows)))

## pooled-parental MAF per SNP (folded pooled parental allele frequency)
par_freq_snp <- colMeans(GTs_all[parent_rows, , drop = FALSE], na.rm = TRUE) / 2
pmaf_snp     <- pmin(par_freq_snp, 1 - par_freq_snp)          # marker-named

## =========================================================================
## per-SNP sorting
## =========================================================================
cat("\n[per-SNP] ohta_fast_prepare + parallelism_stats over", ncol(GTs_all), "markers ...\n")
prep_snp <- ohta_fast_prepare(GTs_all, pops = pops)
ps_snp <- parallelism_stats(prep_snp, hybrid_pops = hybrid_pops,
                            aqu_pops = aqu_pops, pol_pops = pol_pops,
                            fix_th = FIX_TH, DI = DI_vec, min_DI = NULL,
                            parent_maf = pmaf_snp, min_parent_maf = MIN_PMAF,
                            sort_rule = SORT_RULE, alpha = ALPHA)
saveRDS(ps_snp, file.path(OUTDIR, "di25_sorting_snp.rds"))

## =========================================================================
## per-eMLG sorting (5 cM clustering; consensus for >2 markers, else rep SNP)
## =========================================================================
res <- readRDS(CLUST); g <- res$groups
is_emlg <- g$n_loci > 2
cat(sprintf("\n[per-eMLG] %d units (%d eMLG + %d rep-SNP); building unit matrix ...\n",
            nrow(g), sum(is_emlg), sum(!is_emlg)))
E <- vapply(seq_len(nrow(g)), function(i) {
  if (is_emlg[i]) consensus_dosage(GTs_all, g$members[[i]]) else GTs_all[, g$representative[i]]
}, numeric(nrow(GTs_all)))                                    # individuals x units
colnames(E) <- g$group_id

## per-unit pooled-parental MAF and DI (max over members, ungated covariate)
par_freq_u <- colMeans(E[parent_rows, , drop = FALSE], na.rm = TRUE) / 2
pmaf_u     <- setNames(pmin(par_freq_u, 1 - par_freq_u), g$group_id)
DI_u <- setNames(vapply(g$members, function(mk) {
  v <- DI_vec[mk]; if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
}, numeric(1)), g$group_id)

prep_emlg <- ohta_fast_prepare(E, pops = pops)
ps_emlg <- parallelism_stats(prep_emlg, hybrid_pops = hybrid_pops,
                             aqu_pops = aqu_pops, pol_pops = pol_pops,
                             fix_th = FIX_TH, DI = DI_u, min_DI = NULL,
                             parent_maf = pmaf_u, min_parent_maf = MIN_PMAF,
                             sort_rule = SORT_RULE, alpha = ALPHA)
ps_emlg <- g[, .(group_id, n_loci, is_emlg)][ps_emlg, on = c(group_id = "marker")]
saveRDS(ps_emlg, file.path(OUTDIR, "di25_sorting_emlg.rds"))

## =========================================================================
## tau sweep (phi = 0.85 fixed): classify + tally at each level
## =========================================================================
tally_level <- function(ps, level) {
  base <- ps[differentiated == TRUE & n_obs > 0]
  rbindlist(lapply(TAU_GRID, function(tau) {
    cls <- classify_sort(base$n_aqu, base$n_pol, base$n_obs,
                         sort_th = tau, sort_rule = SORT_RULE, alpha = ALPHA)
    n_diff <- nrow(base)
    n_aqu <- sum(cls == "aquilonia"); n_pol <- sum(cls == "polyctena")
    n_unres <- sum(cls == "unresolved"); n_amb <- sum(cls == "ambiguous")
    n_sorted <- n_aqu + n_pol + n_unres + n_amb
    data.table(level = level, tau = tau, n_differentiated = n_diff,
               n_sorted = n_sorted, pct_sorted = 100 * n_sorted / n_diff,
               toward_aqu = n_aqu, toward_pol = n_pol,
               dir_unresolved = n_unres, ambiguous = n_amb,
               pct_aqu_of_resolved = 100 * n_aqu / (n_aqu + n_pol))
  }))
}
sweep <- rbind(tally_level(ps_snp, "SNP"), tally_level(ps_emlg, "eMLG"))
saveRDS(sweep, file.path(OUTDIR, "di25_sorting_sweep.rds"))

cat("\n===== ancestry sorting, phi = 0.85, tau sweep =====\n")
print(sweep[, .(level, tau, n_differentiated, n_sorted, pct_sorted = round(pct_sorted, 1),
                toward_aqu, toward_pol, dir_unresolved, ambiguous,
                pct_aqu_of_resolved = round(pct_aqu_of_resolved, 1))])
cat("\n[di25-sorting] done\n")
