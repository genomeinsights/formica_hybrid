## =========================================================================
## module_population_partitioning -- candidate-locus follow-up 15: import
## and validate module_manuscript_rho05's frozen Stage-1-direct BayPass
## candidate outputs (PC1, PC2, bio_winter, mitoC2).
##
## module_manuscript_rho05 is treated as READ-ONLY, authoritative upstream:
## Stage-1-direct best-SNP candidates, BayPass BF/C2 statistics,
## floor-survivor status, LD-cluster identities, DI, recombination and
## sorting are all IMPORTED, never recomputed or redefined. This script
## does NOT rerun LD reduction and does NOT alter candidate membership.
##
## Population set: 19 hybrid populations, ALAND EXCLUDED throughout --
## matches exactly the population set module_manuscript_rho05's BayPass
## scans were computed on (confirmed with the user 2026-09-14; this
## module's own primary pipeline elsewhere uses all 20 populations, but
## mixing population sets here would decouple "what made a locus a
## candidate" from "what we analyse").
##
## No pre-saved ORIENTED (toward F. aquilonia) allele-frequency matrix
## exists anywhere in module_manuscript_rho05 -- moduleA_stage1_cluster_
## sorting.rds computes orientation internally via parallelism_stats() but
## never saves it. This script reconstructs it using the IDENTICAL formula
## pp_prep_units.R:120-132 uses (sign of the parental population mean
## difference), then validates agreement against moduleA_stage1_cluster_
## sorting.rds$uni_score (built under the same orient="parents" convention)
## -- this is the "check that allele orientations agree between modules"
## requirement. NOTE: best$stats$flipped is NOT ancestry orientation (it is
## an internal LD-polarity bookkeeping flag from the clustering step) and
## must not be used for this check.
##
## Chromosome set for this candidate universe is 26, NOT 27: Chr23 has zero
## Stage-1-direct units (and zero large-cluster full-genome reference-panel
## units, checked separately in R/16_structure_reference.R). Always derive
## chrs from the data, never hardcode a chromosome count.
##
## Run from the formica_hybrid repo root:
##   Rscript module_population_partitioning/R/15_candidate_import.R
## Writes: module_population_partitioning/data/followup/15_manifest.rds
##         module_population_partitioning/data/followup/15_candidate_data.rds
## =========================================================================
suppressMessages(library(data.table))
t_start <- Sys.time()
MSR   <- "module_manuscript_rho05"                 # upstream module root (read-only)
OUTDIR <- "module_population_partitioning/data/followup"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## adjustable, documented design constants (flagged per the plan, not hidden)
BDMI_CUTOFF_IDX   <- 13L   # "good coverage" middle-ground cutoff, matching module_di25's own single-cutoff precedent
N_TOP             <- 10L   # equal-sized top-ranked sensitivity set, common across all 4 targets
ORIENT_MISMATCH_TOL <- 0.05  # max tolerated fraction of strongly-sorted units where our orientation disagrees with uni_score

## ---------------------------------------------------------------------
## 0. manifest scaffold -- filled in as each object is read below, so the
##    provenance record can never drift from what was actually read.
## ---------------------------------------------------------------------
manifest <- data.table(path = character(0), role = character(0), mtime_recorded = as.POSIXct(character(0)),
                       sha256_first12 = character(0), n_rows_read = integer(0), status = character(0))
record <- function(path, role, n_rows, status = "authoritative, currently read") {
  full <- if (file.exists(path)) path else file.path(MSR, path)
  sha <- tryCatch(substr(system(sprintf("shasum -a 256 '%s' | cut -d' ' -f1", full), intern = TRUE), 1, 12),
                  error = function(e) NA_character_)
  manifest <<- rbind(manifest, data.table(path = path, role = role,
                                          mtime_recorded = file.info(full)$mtime,
                                          sha256_first12 = if (length(sha)) sha else NA_character_,
                                          n_rows_read = n_rows, status = status))
}

## ---------------------------------------------------------------------
## 1. Stage-1-direct unit universe + best-SNP genotype dosage
## ---------------------------------------------------------------------
bs <- readRDS(file.path(MSR, "data/moduleB_stage1_units_bestsnp.rds"))
record("module_manuscript_rho05/data/moduleB_stage1_units_bestsnp.rds",
       "Stage-1-direct unit universe (group_id/Chr/n_loci) + best-SNP genotype dosage (raw, unoriented, hybrids-only)",
       nrow(bs$groups))

stopifnot(
  "Stage-1-direct unit universe must have exactly 18,361 rows" = nrow(bs$groups) == 18361L,
  "groups$group_id must be unique" = !anyDuplicated(bs$groups$group_id),
  "best$stats$group_id must be unique" = !anyDuplicated(bs$best$stats$group_id),
  "best$stats$best_marker must be unique -- one Stage-1 representative per candidate" = !anyDuplicated(bs$best$stats$best_marker),
  "groups and best$stats must be row-aligned" = identical(bs$groups$group_id, bs$best$stats$group_id),
  "geno columns must exactly match group_id order" = identical(colnames(bs$best$geno), bs$best$stats$group_id)
)
GROUP_ORDER_TXT <- file.path(MSR, "baypass_stage1/aland_excluded_S1units/S1units_group_order.txt")
stopifnot("S1units_group_order.txt must exactly match groups$group_id" =
            identical(readLines(GROUP_ORDER_TXT), bs$groups$group_id))
cat(sprintf("[import] Stage-1-direct unit universe: %d units, all uniqueness/alignment checks passed\n", nrow(bs$groups)))

## ---------------------------------------------------------------------
## 2. per-target raw BayPass stat + floor-survivor status (row order
##    IDENTICAL to bs$groups$group_id, verified below)
## ---------------------------------------------------------------------
nul_pc  <- readRDS(file.path(MSR, "data/moduleB_stage1_S1units_null.rds"))
nul_bw  <- readRDS(file.path(MSR, "data/moduleB_stage1_bio_winter_null.rds"))
nul_mc  <- readRDS(file.path(MSR, "data/moduleB_stage1_mitoC2_null.rds"))
record("module_manuscript_rho05/data/moduleB_stage1_S1units_null.rds", "PC1/PC2 raw BF(dB) + null-calibrated floor-survivor status", nrow(nul_pc))
record("module_manuscript_rho05/data/moduleB_stage1_bio_winter_null.rds", "bio_winter raw BF(dB) + null-calibrated floor-survivor status", nrow(nul_bw))
record("module_manuscript_rho05/data/moduleB_stage1_mitoC2_null.rds", "mitoC2 raw C2/log10p + null-calibrated floor-survivor status (floor3 = log10p>=3 & k3==0)", nrow(nul_mc))
stopifnot("PC1/PC2 null-calibration row order must match groups$group_id" = identical(nul_pc$group_id, bs$groups$group_id),
         "bio_winter null-calibration row order must match groups$group_id" = identical(nul_bw$group_id, bs$groups$group_id),
         "mitoC2 null-calibration row order must match groups$group_id" = identical(nul_mc$group_id, bs$groups$group_id))

## explicitly excluded, stale object -- recorded so the exclusion itself is auditable
record("module_manuscript_rho05/data/moduleB_stage1_mitoC2_regions.rds",
       "STALE -- pre-null-calibration mitoC2 crossings (11 raw, mtime predates the current floor3 calibration which correctly finds 0 floor survivors)",
       0L, status = "excluded -- confirmed stale, superseded by moduleB_stage1_mitoC2_null.rds$floor3")

## ---------------------------------------------------------------------
## 3. consolidated per-unit annotation table (recomb, Chr, Pos, DI, eBF1/2)
## ---------------------------------------------------------------------
an <- readRDS(file.path(MSR, "data/moduleC_stage1_annotations.rds"))
record("module_manuscript_rho05/data/moduleC_stage1_annotations.rds",
       "recomb (cM/Mb)/Chr/Pos/DI(full-genome ungated vintage)/prop_fixed/uni_score/eBF1/eBF2 per unit", nrow(an))
stopifnot("annotation table row order must match groups$group_id" = identical(an$group_id, bs$groups$group_id),
         "recomb must have no missing values across all Stage-1-direct units" = sum(is.na(an$recomb)) == 0L)

## ---------------------------------------------------------------------
## 4. sort_class + uni_score (authoritative, tau sweep) -- NOT row-aligned
##    with bs$groups$group_id, must join explicitly by group_id
## ---------------------------------------------------------------------
sc06 <- readRDS(file.path(MSR, "data/moduleA_stage1_cluster_sorting.rds"))              # tau = 0.6, primary
sc05 <- readRDS(file.path(MSR, "data/moduleA_stage1_cluster_sorting_tau05.rds"))
sc08 <- readRDS(file.path(MSR, "data/moduleA_stage1_cluster_sorting_tau08.rds"))
record("module_manuscript_rho05/data/moduleA_stage1_cluster_sorting.rds",
       "sort_class @ tau=0.6 (primary) + uni_score (orientation ground truth, orient=\"parents\" convention)", nrow(sc06))
record("module_manuscript_rho05/data/moduleA_stage1_cluster_sorting_tau05.rds", "sort_class @ tau=0.5 (sensitivity)", nrow(sc05))
record("module_manuscript_rho05/data/moduleA_stage1_cluster_sorting_tau08.rds", "sort_class @ tau=0.8 (sensitivity)", nrow(sc08))
stopifnot("sort_class objects must cover exactly the same 18,361 group_ids (any order)" =
            setequal(sc06$group_id, bs$groups$group_id) && !anyDuplicated(sc06$group_id) &&
            setequal(sc05$group_id, bs$groups$group_id) && setequal(sc08$group_id, bs$groups$group_id))
m06 <- match(bs$groups$group_id, sc06$group_id)
m05 <- match(bs$groups$group_id, sc05$group_id)
m08 <- match(bs$groups$group_id, sc08$group_id)
stopifnot("sort_class(tau06) join must be complete -- no unmatched group_id" = !anyNA(m06))

## ---------------------------------------------------------------------
## 5. sample metadata + parental genotypes (repo-root, read-only)
## ---------------------------------------------------------------------
e2 <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd <- e2$sample_data_with_parents; setDT(sd)
GTs_with_parents <- e2$GTs_with_parents
record("data/hybrids_and_parents_maf005.Rdata", "sample metadata (Population/Sample_ID) + parental genotype dosage for orientation + parental MAF", nrow(sd))

aqu_ids <- sd[Population == "aquilonia_parent", Sample_ID]
pol_ids <- sd[Population == "polyctena_parent", Sample_ID]
stopifnot("expected 15 aquilonia_parent and 15 polyctena_parent samples" = length(aqu_ids) == 15L, length(pol_ids) == 15L)

hyb_ids <- rownames(bs$best$geno)
mm <- match(hyb_ids, sd$Sample_ID)
stopifnot("every hybrid individual in best$geno must resolve to sample metadata" = !anyNA(mm))
hyb_pop <- sd$Population[mm]
POPS_19 <- sort(setdiff(unique(sd$Population), c("Aland", "aquilonia_parent", "polyctena_parent")))
stopifnot("population set must have exactly 19 hybrid populations (Aland excluded)" = length(POPS_19) == 19L,
         "Aland must be present in metadata (confirms it was deliberately dropped, not silently absent)" = "Aland" %in% unique(sd$Population))
keep_ind <- hyb_pop %in% POPS_19
cat(sprintf("[import] %d/%d hybrid individuals retained across %d populations (Aland's %d individuals excluded)\n",
            sum(keep_ind), length(hyb_ids), length(POPS_19), sum(hyb_pop == "Aland")))

## ---------------------------------------------------------------------
## 6. reconstruct oriented population x unit allele-frequency matrix
##    (identical formula to pp_prep_units.R:120-132), for ALL 18,361 units
## ---------------------------------------------------------------------
E <- bs$best$geno[keep_ind, , drop = FALSE]                       # 155 hybrids (19 pops) x 18,361 units, dosage 0/1/2
pops19 <- hyb_pop[keep_ind]
pop_mean_dosage <- function(G, pop) {
  levs <- unique(pop)
  P <- matrix(NA_real_, length(levs), ncol(G), dimnames = list(levs, colnames(G)))
  for (k in seq_along(levs)) P[k, ] <- colMeans(G[pop == levs[k], , drop = FALSE], na.rm = TRUE)
  P
}
P_hyb <- pop_mean_dosage(E, pops19) / 2                            # 19 pops x 18,361, allele freq [0,1]

par_geno <- GTs_with_parents[c(aqu_ids, pol_ids), bs$best$stats$best_marker, drop = FALSE]
f_aqu_par <- colMeans(par_geno[aqu_ids, , drop = FALSE], na.rm = TRUE) / 2
f_pol_par <- colMeans(par_geno[pol_ids, , drop = FALSE], na.rm = TRUE) / 2
sign_aqu  <- sign(f_aqu_par - f_pol_par)

Fmat_all <- P_hyb[POPS_19, , drop = FALSE]
flip <- which(sign_aqu < 0); undef <- which(is.na(sign_aqu) | sign_aqu == 0)
undef_ids_all <- colnames(Fmat_all)[undef]   ## captured NOW, before any later reordering of Fmat_all's columns
if (length(flip))  Fmat_all[, flip]  <- 1 - Fmat_all[, flip]
if (length(undef)) Fmat_all[, undef] <- NA_real_
stopifnot("Fmat_all columns must exactly match bs$groups$group_id in order" = identical(colnames(Fmat_all), bs$groups$group_id))
cat(sprintf("[import] oriented Fmat_all: %d populations x %d units; %d units with undefined orientation (parents tied/NA)\n",
            nrow(Fmat_all), ncol(Fmat_all), length(undef)))

## ---------------------------------------------------------------------
## 7. orientation-agreement validation against uni_score
## ---------------------------------------------------------------------
uni_score_06 <- sc06$uni_score[m06]
strong <- which(abs(uni_score_06) >= 0.6 & !is.na(sign_aqu) & sign_aqu != 0)
our_direction   <- sign(colMeans(Fmat_all[, strong, drop = FALSE], na.rm = TRUE) - 0.5)
their_direction <- sign(uni_score_06[strong])
mismatch_rate <- mean(our_direction != their_direction, na.rm = TRUE)
cat(sprintf("[import] orientation-agreement check: %d strongly-sorted units (|uni_score|>=0.6), mismatch rate = %.4f (tolerance %.2f)\n",
            length(strong), mismatch_rate, ORIENT_MISMATCH_TOL))
stopifnot("orientation reconstruction disagrees with moduleA_stage1_cluster_sorting.rds$uni_score beyond tolerance -- stop, do not proceed with a possibly-flipped Fmat" =
            mismatch_rate <= ORIENT_MISMATCH_TOL)

## ---------------------------------------------------------------------
## 8. parental MAF (folded), all 18,361 units
## ---------------------------------------------------------------------
pf <- colMeans(par_geno[c(aqu_ids, pol_ids), , drop = FALSE], na.rm = TRUE) / 2
pmaf <- pmin(pf, 1 - pf)

## ---------------------------------------------------------------------
## 9. BDMI overlap flag -- lifted verbatim from module_di25/R/bdmi_sorting_circos.R
##    (merge_iv() + in_intervals(), pure data.table/base R, no bedtools/GRanges)
## ---------------------------------------------------------------------
merge_iv <- function(s, e) {
  o <- order(s); s <- s[o]; e <- e[o]
  cs <- s[1L]; ce <- e[1L]; outS <- numeric(0); outE <- numeric(0)
  for (i in seq_along(s)[-1L]) {
    if (s[i] <= ce) ce <- max(ce, e[i])
    else { outS <- c(outS, cs); outE <- c(outE, ce); cs <- s[i]; ce <- e[i] }
  }
  list(s = c(outS, cs), e = c(outE, ce))
}
in_intervals <- function(qpos, iv) {
  if (!length(iv$s)) return(logical(length(qpos)))
  brk <- as.vector(rbind(iv$s, iv$e))
  (findInterval(qpos, brk) %% 2L) == 1L
}
bed_files <- list.files("data/liftoff_Frufa_DTOL_PR", pattern = sprintf("^bdmi_candidates\\.cutoff_%d_.*\\.bed$", BDMI_CUTOFF_IDX), full.names = TRUE)
stopifnot("expected exactly one BDMI bed file for the chosen cutoff index" = length(bed_files) == 1L)
bed <- fread(bed_files[1], header = FALSE, col.names = c("chr", "start", "end"))
bed[, chr := sub("chromosome_", "Chr", chr)]
record(bed_files[1], sprintf("BDMI candidate-region overlap, cutoff index %d (documented default, adjustable -- module_di25's own \"good coverage\" middle-ground single-cutoff precedent)", BDMI_CUTOFF_IDX), nrow(bed))

bdmi_flag <- logical(nrow(an))
for (cc in unique(an$Chr)) {
  bedc <- bed[chr == cc]
  if (!nrow(bedc)) next
  iv <- merge_iv(bedc$start, bedc$end)
  idx <- which(an$Chr == cc)
  bdmi_flag[idx] <- in_intervals(an$Pos[idx], iv)
}
cat(sprintf("[import] BDMI overlap (cutoff_%d): %d/%d Stage-1-direct units overlap a candidate region\n", BDMI_CUTOFF_IDX, sum(bdmi_flag), length(bdmi_flag)))

## ---------------------------------------------------------------------
## 10. assemble the unified unit table (u), 18,361 rows, ChrNum/Pos ordered
## ---------------------------------------------------------------------
u <- copy(an)
u[, ChrNum := as.integer(sub("Chr", "", Chr))]
u[, sort_class := sc06$sort_class[m06]]
u[, sort_class_tau05 := sc05$sort_class[match(group_id, sc05$group_id)]]
u[, sort_class_tau08 := sc08$sort_class[match(group_id, sc08$group_id)]]
u[, uni_score := uni_score_06]
u[, pmaf := pmaf]
u[, bdmi_overlap := bdmi_flag]
u[, BF1 := nul_pc$BF1][, floor1 := nul_pc$floor1]
u[, BF2 := nul_pc$BF2][, floor2 := nul_pc$floor2]
u[, BF_bw := nul_bw$BF][, floor_bw := nul_bw$floor]
u[, C2_mc := nul_mc$C2][, log10p_mc := nul_mc$log10p][, floor_mc := nul_mc$floor3]

setorder(u, ChrNum, Pos)
ord <- match(u$group_id, colnames(Fmat_all))
Fmat_all <- Fmat_all[, ord, drop = FALSE]
stopifnot("final unit table must have exactly the 18,361 Stage-1-direct units -- check for a legacy/stale input" = nrow(u) == 18361L,
         "Fmat_all columns must exactly match u$group_id in order after final sort" = identical(colnames(Fmat_all), u$group_id))
chrs26 <- sort(unique(u$Chr))
cat(sprintf("[import] chromosome set for this candidate universe: %d chromosomes (%s)\n", length(chrs26), paste(chrs26, collapse = ",")))
stopifnot("expected exactly 26 chromosomes in the Stage-1-direct candidate universe (Chr23 absent)" = length(chrs26) == 26L)

## ---------------------------------------------------------------------
## 11. per-target raw-candidate / floor-survivor / top-N sets, and
##     completeness validation of every candidate -> Stage-1 mapping
## ---------------------------------------------------------------------
TARGETS <- list(
  PC1        = list(stat_col = "BF1",     floor_col = "floor1",   raw_rule = quote(BF1 >= 15),
                    is_discovery = FALSE, no_discovery_set = FALSE, rank_desc = TRUE),
  PC2        = list(stat_col = "BF2",     floor_col = "floor2",   raw_rule = quote(BF2 >= 15),
                    is_discovery = FALSE, no_discovery_set = FALSE, rank_desc = TRUE),
  bio_winter = list(stat_col = "BF_bw",   floor_col = "floor_bw", raw_rule = quote(BF_bw >= 15),
                    is_discovery = TRUE,  no_discovery_set = FALSE, rank_desc = TRUE),
  mitoC2     = list(stat_col = "log10p_mc", floor_col = "floor_mc", raw_rule = quote(log10p_mc >= 3),
                    is_discovery = FALSE, no_discovery_set = TRUE,  rank_desc = TRUE)
)

cand_ids <- list(); floor_ids <- list(); topN_ids <- list()
for (tgt in names(TARGETS)) {
  cfg <- TARGETS[[tgt]]
  raw <- u[eval(cfg$raw_rule), group_id]
  flr <- u[[cfg$floor_col]]; flr_ids <- u[flr == TRUE, group_id]
  ord_stat <- order(-u[[cfg$stat_col]], na.last = NA)
  top_ids <- u$group_id[ord_stat][seq_len(min(N_TOP, length(ord_stat)))]

  hit <- match(raw, bs$groups$group_id)
  if (anyNA(hit)) stop(sprintf("[%s] %d/%d raw-candidate group_ids have no Stage-1 representative -- incomplete mapping, stopping", tgt, sum(is.na(hit)), length(raw)))
  if (anyDuplicated(raw)) stop(sprintf("[%s] duplicated raw-candidate group_ids -- stopping", tgt))

  cand_ids[[tgt]] <- raw; floor_ids[[tgt]] <- flr_ids; topN_ids[[tgt]] <- top_ids
  cat(sprintf("[import] %-10s: %4d raw candidates (%s >= threshold), %2d floor survivors, top-%d by %s\n",
              tgt, length(raw), cfg$stat_col, length(flr_ids), N_TOP, cfg$stat_col))
}
cat(sprintf("\n[import] bio_winter floor survivors: %s\n", paste(floor_ids$bio_winter, collapse = ", ")))

## flag units with UNDEFINED orientation (fixed for the same allele in BOTH
## parental species, or missing parental genotype) within each target's sets
## -- their oriented aquilonia frequency is NA by construction (real biology,
## not missing data); this must be surfaced, not silently absorbed as NA
## (undef_ids_all was captured earlier, BEFORE Fmat_all's columns were
## reordered to match u$group_id -- using colnames(Fmat_all) here would be
## silently wrong since positional `undef` indices no longer line up)
undef_counts <- data.table(target = character(0), status = character(0), n_undef = integer(0), n_total = integer(0))
for (tgt in names(TARGETS)) {
  for (status in c("raw", "floor", "topN")) {
    ids <- switch(status, raw = cand_ids[[tgt]], floor = floor_ids[[tgt]], topN = topN_ids[[tgt]])
    undef_counts <- rbind(undef_counts, data.table(target = tgt, status = status,
                                                    n_undef = sum(ids %in% undef_ids_all), n_total = length(ids)))
  }
}
cat("\n[import] units with UNDEFINED ancestry orientation (fixed for the same allele in both parental\n")
cat("    species, or missing parental genotype -- oriented frequency is NA by construction, not missing data):\n")
print(undef_counts[n_total > 0])
stopifnot("bio_winter must have exactly 10 floor survivors" = length(floor_ids$bio_winter) == 10L,
         "PC1 must have exactly 1 floor survivor" = length(floor_ids$PC1) == 1L,
         "PC2 must have exactly 2 floor survivors" = length(floor_ids$PC2) == 2L,
         "mitoC2 must have exactly 0 floor survivors (11 raw crossings)" = length(floor_ids$mitoC2) == 0L,
         "mitoC2 must have exactly 11 raw candidates" = length(cand_ids$mitoC2) == 11L)

## ---------------------------------------------------------------------
## 12. save
## ---------------------------------------------------------------------
saveRDS(manifest, file.path(OUTDIR, "15_manifest.rds"))
cat("\n[import] provenance manifest:\n"); print(manifest[, .(path, role, n_rows_read, status)], width = 200)

result <- list(u = u, Fmat_all = Fmat_all, POPS_19 = POPS_19, TARGETS = TARGETS,
              cand_ids = cand_ids, floor_ids = floor_ids, topN_ids = topN_ids,
              chrs = chrs26, N_TOP = N_TOP, BDMI_CUTOFF_IDX = BDMI_CUTOFF_IDX,
              orient_mismatch_rate = mismatch_rate, undef_orientation_units = length(undef),
              undef_ids_all = undef_ids_all, undef_counts = undef_counts,
              session_info = sessionInfo(), run_time = Sys.time(),
              elapsed_secs = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
saveRDS(result, file.path(OUTDIR, "15_candidate_data.rds"))
cat(sprintf("\n[import] saved -> %s, %s (elapsed %.1fs)\n",
            file.path(OUTDIR, "15_manifest.rds"), file.path(OUTDIR, "15_candidate_data.rds"), result$elapsed_secs))
