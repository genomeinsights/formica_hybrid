## =========================================================
## module_localscore_crosscheck -- Stage-1-unit genomic Position table,
## using best_marker positions from the outset (not core_snp/representative)
## =========================================================
## Fresh implementation. Reads the already-authoritative Stage-1 best-marker
## object (module_manuscript_rho05/data/moduleB_stage1_units_bestsnp.rds --
## a deterministic function of authoritative upstream inputs, not an
## exploratory-module artifact) and assigns each of the 18,361 Stage-1
## clusters (n_snps >= 5) its best_marker's real genomic Chr/Pos, validated
## against the marker map rather than string-parsed from the marker ID.
##
## This table is the position backbone for constructing Stage-1-resolution
## local-score inputs later -- it replaces the exploratory module's use of
## core_snp/representative for Stage-1 unit placement.
##
## Reads:
##   module_manuscript_rho05/data/moduleB_stage1_units_bestsnp.rds
##   data/hybrids_only_maf005.Rdata   (map_hyb_005, for Chr/Pos validation)
##   module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/S1units_group_order.txt
## Writes:
##   module_localscore_crosscheck/data/stage1_unit_positions.tsv
##   module_localscore_crosscheck/manifests/stage1_unit_positions_manifest.md
##
## Run from the repo root:
##   Rscript module_localscore_crosscheck/R/02_stage1_best_marker_positions.R
## =========================================================

suppressMessages(library(data.table))

BEST_RDS  <- "module_manuscript_rho05/data/moduleB_stage1_units_bestsnp.rds"
GROUP_ORDER_FILE <- "module_manuscript_rho05/baypass_stage1/aland_excluded_S1units/S1units_group_order.txt"
OUT_DIR   <- "module_localscore_crosscheck/data"
MANIF_DIR <- "module_localscore_crosscheck/manifests"
dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(MANIF_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- 1. load authoritative Stage-1 best-marker object -----------------------
obj <- readRDS(BEST_RDS)
stopifnot(all(c("groups", "best") %in% names(obj)))
stopifnot(all(c("stats", "geno") %in% names(obj$best)))

groups <- as.data.table(obj$groups)
stats  <- as.data.table(obj$best$stats)
stopifnot(nrow(groups) == nrow(stats))
stopifnot(identical(groups$group_id, stats$group_id))

## ---- 2. load marker map for independent Chr/Pos validation ------------------
load("data/hybrids_only_maf005.Rdata")   # -> map_hyb_005
stopifnot("map_hyb_005" %in% ls())
map_key <- setNames(map_hyb_005$Pos, map_hyb_005$marker)
map_chr <- setNames(map_hyb_005$Chr, map_hyb_005$marker)

stopifnot(
  "every best_marker must exist in the marker map" =
    all(stats$best_marker %in% map_hyb_005$marker)
)

pos <- as.integer(map_key[stats$best_marker])
chr <- as.character(map_chr[stats$best_marker])

## ---- 3. cross-check against best_marker's own Chr:Pos string encoding -------
parsed_chr <- sub(":.*$", "", stats$best_marker)
parsed_pos <- as.integer(sub("^.*:", "", stats$best_marker))
stopifnot(
  "map-derived Chr disagrees with best_marker's own Chr:Pos encoding" =
    identical(chr, parsed_chr),
  "map-derived Pos disagrees with best_marker's own Chr:Pos encoding" =
    identical(pos, parsed_pos)
)

## ---- 4. cross-check against the group_id order used to build BayPass inputs -
group_order <- readLines(GROUP_ORDER_FILE)
stopifnot(
  "group_id set must match the BayPass Stage-1-unit run order exactly" =
    identical(sort(group_order), sort(groups$group_id))
)

## ---- 5. assemble output table -----------------------------------------------
out <- data.table(
  group_id       = groups$group_id,
  Chr            = chr,
  Pos            = pos,
  best_marker    = stats$best_marker,
  representative = groups$representative,
  best_r         = stats$best_r,
  rep_is_best    = stats$rep_is_best,
  n_loci         = groups$n_loci
)
setorder(out, Chr, Pos)

## ---- 6. validation ------------------------------------------------------------
stopifnot(
  "no duplicate group_id" = !anyDuplicated(out$group_id),
  "no duplicate (Chr,Pos) -- each best_marker must map to a unique physical site" =
    !anyDuplicated(out[, .(Chr, Pos)]),
  "no missing Pos" = !anyNA(out$Pos)
)
n_rep_is_best <- sum(out$rep_is_best)
message(sprintf(
  "[stage1pos] %d / %d Stage-1 units: representative == best_marker (%.1f%%)",
  n_rep_is_best, nrow(out), 100 * n_rep_is_best / nrow(out)))

## ---- 7. write -------------------------------------------------------------------
fwrite(out, file.path(OUT_DIR, "stage1_unit_positions.tsv"), sep = "\t")

md <- c(
  "# Stage-1 unit position table -- manifest",
  sprintf("Source: `%s`", BEST_RDS),
  sprintf("N Stage-1 units: %d", nrow(out)),
  sprintf("Units where representative == best_marker: %d (%.1f%%)",
          n_rep_is_best, 100 * n_rep_is_best / nrow(out)),
  "",
  "Validation performed:",
  "- every best_marker exists in map_hyb_005 (the authoritative marker map)",
  "- map-derived Chr/Pos agrees exactly with best_marker's own Chr:Pos string",
  "- Stage-1 unit group_id set matches the BayPass run's S1units_group_order.txt exactly",
  "- no duplicate group_id, no duplicate (Chr,Pos), no missing Pos",
  "",
  "Output: module_localscore_crosscheck/data/stage1_unit_positions.tsv",
  "(group_id, Chr, Pos, best_marker, representative, best_r, rep_is_best, n_loci)"
)
writeLines(md, file.path(MANIF_DIR, "stage1_unit_positions_manifest.md"))

message("[stage1pos] wrote ", file.path(OUT_DIR, "stage1_unit_positions.tsv"))
message("[stage1pos] wrote ", file.path(MANIF_DIR, "stage1_unit_positions_manifest.md"))
