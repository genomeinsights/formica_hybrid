## =========================================================
## module_di25_rho05 -- ancestry sorting on the DI25 data, per-eMLG side only
## =========================================================
## Mirrors module_di25/R/di25_sorting.R exactly, but consumes the
## min_r2_rho = 0.5 clustering (di25_clustering_cM5_rho05.rds) instead of
## the canonical min_r2 = 0.2 one. The per-SNP sorting (ps_snp) does not
## depend on the eMLG clustering AT ALL -- it's computed directly from the
## 51,612-marker genotype matrix -- so it is IDENTICAL between arms and is
## reused unchanged from module_di25/data/di25_sorting_snp.rds rather than
## recomputed (paired comparison, same rationale as the clustering step).
##
## Same conventions as the canonical script: best-SNP (fill = FALSE),
## no in-function differentiation/MAF gate (already applied upstream via the
## DI > -25 input panel), phi = 0.85, tau in {0.5,0.6,0.7,0.8}, sort_rule="binom".
##
## Run from the repo root:  Rscript module_di25_rho05/R/di25_sorting_rho05.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")             # ohta_fast_prepare(), eMLG_best_snp()
source("moduleA_sorting/R/parallelism_stats.R")    # parallelism_stats(), classify_sort()

## ---- parameters ---------------------------------------------------------
CLUST_RHO <- "module_di25_rho05/data/di25_clustering_cM5_rho05.rds"
CANON_SNP <- "module_di25/data/di25_sorting_snp.rds"    # reused unchanged (see header)
CANON_SWEEP <- "module_di25/data/di25_sorting_sweep.rds"
INPUTS    <- "module_di25/data/di25_inputs.rds"
OUTDIR    <- "module_di25_rho05/data"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
TAU_GRID  <- c(0.5, 0.6, 0.7, 0.8)
FIX_TH    <- 0.15
SORT_RULE <- "binom"
ALPHA     <- 0.05

## ---- inputs -------------------------------------------------------------
inp <- readRDS(INPUTS); map <- inp$map
e2  <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd  <- e2$sample_data_with_parents
DI_vec <- setNames(e2$map_hyb_005$DiagnosticIndex, e2$map_hyb_005$marker)

GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)
keep_ind <- rownames(GTs_all) %in% sd$Sample_ID
GTs_all  <- GTs_all[keep_ind, ]
pops     <- sd$Population[match(rownames(GTs_all), sd$Sample_ID)]
aqu_pops <- "aquilonia_parent"; pol_pops <- "polyctena_parent"
hybrid_pops <- setdiff(unique(pops), c(aqu_pops, pol_pops))
parent_rows <- grepl("_parent$", pops)
cat(sprintf("Individuals: %d (%d hybrids in %d pops + %d parents)\n",
            nrow(GTs_all), sum(!parent_rows), length(hybrid_pops), sum(parent_rows)))

## ---- per-SNP: reused unchanged from the canonical (0.2) run -------------
stopifnot(file.exists(CANON_SNP))
ps_snp <- readRDS(CANON_SNP)
cat("[per-SNP] reused unchanged from ", CANON_SNP, " (", nrow(ps_snp), " markers)\n", sep = "")

## =========================================================================
## per-eMLG sorting (5 cM clustering, min_r2_rho = 0.5) -- BEST-SNP representation
## =========================================================================
res <- readRDS(CLUST_RHO); g <- res$groups
is_emlg <- g$n_loci > 2
best <- eMLG_best_snp(res, inp$GTs_hyb, fill = FALSE)
bm   <- setNames(best$stats$best_marker, best$stats$group_id)
rep_snp <- g$representative
rep_snp[is_emlg] <- bm[g$group_id[is_emlg]]
stopifnot(!anyNA(rep_snp), all(rep_snp %in% colnames(GTs_all)))
cat(sprintf("\n[per-eMLG, rho05] %d units (%d eMLG best-SNP + %d rep-SNP); building unit matrix ...\n",
            nrow(g), sum(is_emlg), sum(!is_emlg)))
E <- GTs_all[, rep_snp, drop = FALSE]
colnames(E) <- g$group_id

par_freq_u <- colMeans(E[parent_rows, , drop = FALSE], na.rm = TRUE) / 2
pmaf_u     <- setNames(pmin(par_freq_u, 1 - par_freq_u), g$group_id)
DI_u       <- setNames(DI_vec[rep_snp], g$group_id)

prep_emlg <- ohta_fast_prepare(E, pops = pops)
ps_emlg <- parallelism_stats(prep_emlg, hybrid_pops = hybrid_pops,
                             aqu_pops = aqu_pops, pol_pops = pol_pops,
                             fix_th = FIX_TH, DI = DI_u, min_DI = NULL,
                             parent_maf = pmaf_u, min_parent_maf = NULL,
                             sort_rule = SORT_RULE, alpha = ALPHA)
ps_emlg <- g[, .(group_id, n_loci, is_emlg)][ps_emlg, on = c(group_id = "marker")]
ps_emlg[, unit_marker := setNames(rep_snp, g$group_id)[group_id]]
setnames(ps_emlg, "DI", "current_map_DI")
saveRDS(ps_emlg, file.path(OUTDIR, "di25_sorting_emlg_rho05.rds"))

## =========================================================================
## tau sweep (phi = 0.85 fixed) -- rho05 eMLG, alongside the reused SNP level
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
sweep_rho05 <- rbind(tally_level(ps_snp, "SNP"), tally_level(ps_emlg, "eMLG_rho05"))
saveRDS(sweep_rho05, file.path(OUTDIR, "di25_sorting_sweep_rho05.rds"))

cat("\n===== ancestry sorting, phi = 0.85, tau sweep -- min_r2_rho = 0.5 =====\n")
print(sweep_rho05[, .(level, tau, n_differentiated, n_sorted, pct_sorted = round(pct_sorted, 1),
                       toward_aqu, toward_pol, dir_unresolved, ambiguous,
                       pct_aqu_of_resolved = round(pct_aqu_of_resolved, 1))])

## ---- side-by-side vs the canonical (min_r2 = 0.2) eMLG sweep -------------
if (file.exists(CANON_SWEEP)) {
  sweep_canon <- readRDS(CANON_SWEEP)[level == "eMLG"]
  cmp <- rbind(
    sweep_canon[, .(level = "eMLG_0.20", tau, pct_sorted = round(pct_sorted, 1),
                     toward_aqu, toward_pol, dir_unresolved, ambiguous,
                     pct_aqu_of_resolved = round(pct_aqu_of_resolved, 1))],
    sweep_rho05[level == "eMLG_rho05", .(level, tau, pct_sorted = round(pct_sorted, 1),
                     toward_aqu, toward_pol, dir_unresolved, ambiguous,
                     pct_aqu_of_resolved = round(pct_aqu_of_resolved, 1))]
  )
  setorder(cmp, tau, level)
  cat("\n===== eMLG-level comparison: min_r2 = 0.2 vs min_r2_rho = 0.5 =====\n")
  print(cmp)
}
cat("\n[di25-sorting-rho05] done\n")
