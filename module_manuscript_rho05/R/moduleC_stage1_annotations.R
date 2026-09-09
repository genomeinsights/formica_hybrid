## =========================================================
## module_manuscript_rho05 -- Module C annotations, Stage-1-direct universe
## =========================================================
## Adapted from moduleC_climate_vs_sorting/R/moduleC_annotations.R for the
## 18,361 Stage-1 units (n_snps>=5, best-SNP represented). Single universe --
## unlike the canonical (min=5/min=10) grid, Stage-1 has no second min level,
## so MINS is fixed to 5 throughout (moduleC_stage1_null_regen.R reuses the
## shared min x tau grid code with a one-element MINS for that reason, not
## because the grid concept applies here).
##
## Reads : module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/
##           {PC1,PC2}_S1units_withOmega_summary_betai_reg.out, S1units_group_order.txt
##         module_manuscript_rho05/data/moduleB_stage1_units_bestsnp.rds  (groups, best$stats)
##         module_manuscript_rho05/data/moduleA_stage1_cluster_sorting_{tau05,tau06,tau08}.rds
##         data/Frufa_DTOL_PR.ref_genome.recmap
## Writes: module_manuscript_rho05/data/moduleC_stage1_annotations.rds
##
## Run from the repo root: Rscript module_manuscript_rho05/R/moduleC_stage1_annotations.R
## =========================================================

suppressMessages(library(data.table))
source("moduleC_climate_vs_sorting/R/moduleC_stat_functions.R")   # MODULEC_TAU_SERIES, tauC_stamp

BP_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded_S1units"
DATA   <- "module_manuscript_rho05/data"
OUT    <- file.path(DATA, "moduleC_stage1_annotations.rds")

TAUS   <- MODULEC_TAU_SERIES          # c(0.5, 0.6, 0.8)
TSTAMP <- tauC_stamp(TAUS)            # "tau05" "tau06" "tau08"
PRIMARY_STAMP <- tauC_stamp(MODULEC_TAU_PRIMARY)   # "tau06"

## ---- inputs ---------------------------------------------------------------
grp <- readLines(file.path(BP_DIR, "S1units_group_order.txt"))     # BayPass row order
N <- length(grp)
message("[ann-S1] BayPass row order: ", N, " Stage-1 units")

obj <- readRDS(file.path(DATA, "moduleB_stage1_units_bestsnp.rds"))
groups <- obj$groups                                                # group_id, Chr, representative, n_loci, has_eMLG
bstats <- obj$best$stats                                            # group_id, best_marker, n_loci

he <- groups[, .(group_id, n_loci)]
he[bstats, on = "group_id", best_marker := i.best_marker]
stopifnot("best_marker missing for some Stage-1 unit" = all(!is.na(he$best_marker)),
          "n_loci disagrees between groups and best$stats" =
            all(he$n_loci == bstats$n_loci[match(he$group_id, bstats$group_id)]))

## ---- Module A sorting, tau series (Stage-1 equivalent, already built) ----
cl_tau <- lapply(TSTAMP, function(ts)
  readRDS(file.path(DATA, sprintf("moduleA_stage1_cluster_sorting_%s.rds", ts))))
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

## cross-check cluster size between the bestsnp object and cluster_sorting (must agree)
nl_chk <- he[cl[, .(group_id, n_loci_cl = n_loci)], on = "group_id"]
stopifnot("n_loci disagrees between bestsnp groups and cluster_sorting" =
            all(nl_chk$n_loci == nl_chk$n_loci_cl))

## ---- ordering / membership assertions (fail loudly) ----------------------
stopifnot(
  "duplicate group_id in cluster_sorting"  = !any(duplicated(cl$group_id)),
  "duplicate group_id in bestsnp groups"   = !any(duplicated(he$group_id)),
  "cluster_sorting missing some Stage-1 units" = all(grp %in% cl$group_id),
  "bestsnp groups missing some Stage-1 units"  = all(grp %in% he$group_id),
  "cluster_sorting has extra ids beyond universe" = setequal(cl$group_id, grp),
  "bestsnp groups set != universe"         = setequal(he$group_id, grp)
)

## ---- observed BayPass PC1/PC2 BF, in BayPass MRK order --------------------
b1 <- fread(file.path(BP_DIR, "PC1_S1units_withOmega_summary_betai_reg.out"))
b2 <- fread(file.path(BP_DIR, "PC2_S1units_withOmega_summary_betai_reg.out"))
stopifnot(nrow(b1) == N, nrow(b2) == N,
          "PC1 MRK not 1..N in order" = all(b1$MRK == seq_len(N)),
          "PC2 MRK not 1..N in order" = all(b2$MRK == seq_len(N)))
eBF <- data.table(group_id = grp, eBF1 = b1$`BF(dB)`, eBF2 = b2$`BF(dB)`)

## ---- recombination rate (cM/Mb) at each unit's best_marker ----------------
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
DIR_COLS <- paste0("directional_", TSTAMP)
ann <- data.table(group_id = grp)                                   # canonical order
ann <- eBF[, .(group_id, eBF1, eBF2)][ann, on = "group_id"]
ann <- cl[, c("group_id", "DI", "directional", DIR_COLS, "prop_fixed", "uni_score",
              "sort_class", "differentiated", "n_loci"), with = FALSE][ann, on = "group_id"]
ann <- he[, .(group_id, recomb, Chr = rep_chr, Pos = rep_pos)][ann, on = "group_id"]
setcolorder(ann, "group_id")
ann <- ann[match(grp, group_id)]                                    # enforce BayPass order

## final integrity: order + completeness + one-to-one
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
message(sprintf("[ann-S1] prop_fixed (magnitude): %d NA (established), finite range %.2f..%.2f",
                sum(pf_na), min(ann$prop_fixed, na.rm = TRUE), max(ann$prop_fixed, na.rm = TRUE)))

n_dir_by_tau <- setNames(vapply(DIR_COLS, function(c) sum(ann[[c]] == 1L), integer(1)), TSTAMP)
attr(ann, "meta") <- list(
  config = "Stage-1-direct (n_snps>=5) / aland_excluded / withOmega / 19 pops",
  N_units = N,
  tau_series = TAUS, tau_primary = MODULEC_TAU_PRIMARY, tau_stamp = TSTAMP,
  n_directional = sum(ann$directional == 1L),
  n_directional_by_tau = n_dir_by_tau,
  recomb_source = "map-interpolated cM/Mb at best-SNP marker",
  DI_source = "per-Stage-1-unit best-SNP DI (moduleA_stage1_cluster_sorting.rds)",
  magnitude_source = "prop_fixed = degree of fixation; NOT uni_score",
  orientation_source = "uni_score = signed direction; supplementary only",
  n_prop_fixed_NA = sum(is.na(ann$prop_fixed)),
  built = as.character(Sys.time())
)
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
saveRDS(ann, OUT)

cat(sprintf("\n[ann-S1] wrote %s\n  N=%d  DI %.1f..%.1f  recomb %.2f..%.2f cM/Mb\n",
            OUT, N, min(ann$DI), max(ann$DI), min(ann$recomb), max(ann$recomb)))
cat("  directional by tau: ",
    paste(sprintf("%s=%d (%.1f%%)", TSTAMP, n_dir_by_tau, 100 * n_dir_by_tau / N), collapse = "  "),
    sprintf("  [primary %s]\n", PRIMARY_STAMP))
