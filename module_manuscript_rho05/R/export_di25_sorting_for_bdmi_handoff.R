## =========================================================
## module_manuscript_rho05 -- BED+ export of di25 sorted-SNP calls for handoff
## =========================================================
## Prepares the input a colleague needs to take over the BDMI-overlap /
## direction / recombination-control analyses (currently
## module_di25/R/bdmi_sorting_overlap.R, bdmi_sorting_direction.R,
## bdmi_sorting_recomb_controlled.R) -- these are di25-only (the 51,612-marker
## DI25 diagnostic panel; nothing genome-wide uses BDMI overlap).
##
## Matches the coordinate convention of her existing BDMI candidate region
## files (data/liftoff_Frufa_DTOL_PR/bdmi_candidates.cutoff_*.bed:
## "chromosome_N", 0-based start), so she can bedtools-intersect directly with
## no renaming on her end.
##
## Columns beyond the first 3 (standard BED) are extra ("BED+") -- fine for
## bedtools/pandas/R, but not for strict-BED-only tools (e.g. the UCSC
## browser).
##
## Per-population near-fixation columns (popfix_<population>, phi = 0.85,
## code 0 = not near-fixed, 1 = near-fixed F. aquilonia, 2 = near-fixed
## F. polyctena) let her build her own null designs (e.g. a leave-one-
## population-out robustness check) rather than being limited to the
## aggregate n_aqu/n_pol/n_obs counts this repo's own overlap test used.
##
## sort_class_tauXX columns are recomputed at all four Module-A tau levels
## (0.5, 0.6, 0.7, 0.8), phi = 0.85 fixed, sort_rule = "binom", alpha = 0.05
## -- identical conventions to di25_sorting.R / bdmi_sorting_*.R.
##
## Run from the repo root:  Rscript module_manuscript_rho05/R/export_di25_sorting_for_bdmi_handoff.R
## =========================================================

suppressMessages(library(data.table))
source("moduleA_sorting/R/parallelism_stats.R")   # classify_sort()

OUTDIR <- "module_manuscript_rho05/data"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
OUT_BED <- file.path(OUTDIR, "di25_sorted_snps_for_bdmi_handoff.bed")
TAU_GRID <- c(0.5, 0.6, 0.7, 0.8)
PHI <- 0.85

## ---- aggregate sorting table (n_aqu/n_pol/n_obs, differentiated, DI) -----
ps <- readRDS("module_di25/data/di25_sorting_snp.rds")

## ---- per-population near-fixation codes (phi = 0.85), SNP level ----------
inp <- readRDS("module_di25/data/di25_inputs.rds")
e2  <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd  <- e2$sample_data_with_parents
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)
pops    <- sd$Population[match(rownames(GTs_all), sd$Sample_ID)]
faqu <- grep("^Faqu", rownames(GTs_all)); fpol <- grep("^Fpol", rownames(GTs_all))
hybrid_pops <- sort(setdiff(unique(pops[!is.na(pops)]), c("aquilonia_parent", "polyctena_parent")))

M <- t(GTs_all)   # markers x individuals
flip <- which(rowMeans(M[, faqu, drop = FALSE], na.rm = TRUE) <
              rowMeans(M[, fpol, drop = FALSE], na.rm = TRUE))
M[flip, ] <- 2 - M[flip, ]   # 2 = F. aquilonia allele, oriented

F <- vapply(hybrid_pops, function(p) rowMeans(M[, which(pops == p), drop = FALSE], na.rm = TRUE) / 2,
            numeric(nrow(M)))   # markers x pops, oriented aquilonia freq
popfix_code <- matrix(0L, nrow(F), ncol(F), dimnames = list(rownames(M), hybrid_pops))
popfix_code[F >= PHI]       <- 1L   # near-fixed F. aquilonia
popfix_code[F <= (1 - PHI)] <- 2L   # near-fixed F. polyctena
popfix_dt <- as.data.table(popfix_code, keep.rownames = "marker")
setnames(popfix_dt, hybrid_pops, paste0("popfix_", hybrid_pops))

## ---- sort_class at all four tau levels ------------------------------------
for (tau in TAU_GRID) {
  cn <- sprintf("sort_class_tau%s", sub("\\.", "", format(tau, nsmall = 1)))
  ps[[cn]] <- classify_sort(ps$n_aqu, ps$n_pol, ps$n_obs, sort_th = tau, sort_rule = "binom", alpha = 0.05)
}

## ---- assemble BED+ table --------------------------------------------------
map <- inp$map[, .(marker, Chr, Pos)]
tab <- map[ps, on = "marker"][popfix_dt, on = "marker"]
tab[, `:=`(chrom = paste0("chromosome_", sub("Chr", "", Chr)), start = Pos - 1L, end = Pos)]

sort_cols <- grep("^sort_class_tau", names(tab), value = TRUE)
popfix_cols <- paste0("popfix_", hybrid_pops)
tab <- tab[, c("chrom", "start", "end", "marker", "current_map_DI", "differentiated",
               "parent_maf", "n_aqu", "n_pol", "n_obs", sort_cols, popfix_cols), with = FALSE]
setorder(tab, chrom, start)

fwrite(tab, OUT_BED, sep = "\t", quote = FALSE)
cat(sprintf("wrote %s (%d markers x %d columns)\n", OUT_BED, nrow(tab), ncol(tab)))

## ---- companion README ------------------------------------------------------
readme <- c(
  "di25_sorted_snps_for_bdmi_handoff.bed -- column guide",
  "========================================================",
  "One row per DI25 diagnostic SNP (51,612 total; DI > -25 ascertainment panel).",
  "",
  "Columns 1-3 (standard BED): chrom, start (0-based), end -- matches the",
  "  coordinate convention of data/liftoff_Frufa_DTOL_PR/bdmi_candidates.*.bed",
  "  ('chromosome_N', Frufa_DTOL_PR reference), so bedtools intersect works with",
  "  no renaming.",
  "Columns beyond 3 are extra ('BED+') -- fine for bedtools/pandas/R, not for",
  "  strict-BED-only tools (e.g. the UCSC browser).",
  "",
  "marker           Chr:Pos id, matches this repo's own marker naming.",
  "current_map_DI   Diagnostic index from the LATER full-data map -- an ungated",
  "                 covariate, NOT the frozen ascertainment DI that defined the",
  "                 DI25 panel (DIEM was run unseeded for this subset, so its DI",
  "                 values don't exactly match; see di25_sorting.R header).",
  "differentiated   TRUE if the locus carries informative parental frequency",
  "                 difference (prerequisite for a meaningful sort call).",
  "parent_maf       Pooled-parental minor allele frequency (ungated covariate).",
  "n_aqu, n_pol     Number of the (up to 20) hybrid populations near-fixed",
  "                 (phi = 0.85) toward F. aquilonia / F. polyctena.",
  "n_obs            Number of populations with an observed, oriented frequency",
  "                 at this locus (denominator for prop_fixed = (n_aqu+n_pol)/n_obs).",
  "",
  "sort_class_tau05/06/07/08",
  "  Direction call at sorting threshold tau in {0.5,0.6,0.7,0.8}, phi = 0.85",
  "  fixed, two-sided exact binomial direction test (alpha = 0.05):",
  "    'aquilonia' / 'polyctena'  -- sorted, direction significant",
  "    'unresolved'               -- sorted, direction not significant",
  "    'ambiguous'                -- sorted, too few near-fixed populations",
  "                                  to ever reach significance",
  "    'unsorted'                 -- prop_fixed below tau",
  "    NA                         -- n_obs == 0",
  "  tau = 0.6 is the manuscript's primary reported operating point.",
  "",
  sprintf("popfix_<population>  (%d columns, one per hybrid population)", length(hybrid_pops)),
  "  Per-population near-fixation code at phi = 0.85, BEFORE aggregation across",
  "  populations -- lets you build your own null designs (e.g. leave-one-",
  "  population-out) rather than being limited to the aggregate n_aqu/n_pol/n_obs.",
  "    0 = not near-fixed in that population",
  "    1 = near-fixed toward F. aquilonia (oriented freq >= 0.85)",
  "    2 = near-fixed toward F. polyctena (oriented freq <= 0.15)",
  sprintf("  Populations: %s", paste(hybrid_pops, collapse = ", ")),
  "",
  "Scope: this table covers the DI25 diagnostic panel only -- the BDMI-overlap",
  "analyses in this repo (module_di25/R/bdmi_sorting_*.R) are applied exclusively",
  "to this 51,612-marker panel, never to the full genome-wide (~1.1M SNP) dataset."
)
writeLines(readme, file.path(OUTDIR, "di25_sorted_snps_for_bdmi_handoff_README.txt"))
cat("wrote", file.path(OUTDIR, "di25_sorted_snps_for_bdmi_handoff_README.txt"), "\n")
