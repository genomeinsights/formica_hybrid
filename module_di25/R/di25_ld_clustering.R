## =========================================================
## module_di25 (high-DI analyses) -- LD complexity reduction for the DI25 markers
## =========================================================
## Fresh, from-scratch LD-decay -> two-stage complexity reduction on the
## high-DI (DI > -25) diagnostic markers ONLY, so the resulting eMLG clusters
## summarise high-DI variation exclusively (unlike the canonical module0
## clustering, whose all-marker clusters also average in linked lower-DI SNPs).
##
## Input : data/species_diagnostic_markers_DI25_20pops.tsv.gz
##   195 individual genotype columns + `chromosome` + `position`, coded
##     0 = missing, 1 = hom F. aquilonia, 2 = het, 3 = hom F. polyctena.
##   Translated here to 012 dosage (1->0, 2->1, 3->2, 0->NA). r^2 is
##   orientation-invariant, so the 012 polarity is irrelevant to clustering.
##
## Conventions (following module0, with the deliberate changes requested):
##   * HYBRIDS ONLY for all LD estimation. Parents are near-fixed at diagnostic
##     markers and would manufacture huge artificial LD -- they are dropped here
##     and kept aside (parent 012 matrix saved) only for later orientation.
##   * No loci-count limit entering the clustering: ld_w_threshold = 0 AND
##     min_n_loci_flag = 1, so EVERY Stage-1 cluster is evaluated for Stage-2
##     merging (the marker set is small, ~51k, so the module0 compute gates that
##     restrict the flagged set are unnecessary and would only over-split).
##   * eMLG consensus stored for clusters with > 2 markers (min_n_loci_eMLG = 3);
##     clusters of 1-2 markers are represented by their representative SNP.
##   * Stage-2 quality gates kept canonical: score_threshold 0.80, min_r2 0.2.
##   * cM distance cap on merging swept for sensitivity: 0.5, 1, 2, 5, 10 cM
##     (fewer SNPs -> a tight cap risks splitting low-density regions).
##
## Run from the repo root:  Rscript module_di25/R/di25_ld_clustering.R
## =========================================================

suppressMessages({
  library(data.table)
  library(igraph)
  library(SNPRelate)
})
devtools::load_all("~/gitlab/LDscnR/")

## ---- parameters ---------------------------------------------------------
INFILE    <- "data/species_diagnostic_markers_DI25_20pops.tsv.gz"
RECMAP    <- "data/Frufa_DTOL_PR.ref_genome.recmap"
OUTDIR    <- "module_di25/data"
GDS_PATH  <- file.path(OUTDIR, "di25_hybrids.gds")
CM_SWEEP  <- c(0.5, 1, 2, 5, 10)          # Stage-2 merge distance caps (cM)
CORES     <- 4
RECOMPUTE <- FALSE                         # TRUE forces LD-decay/Stage1 rebuild

## Stage-2 flagging / eMLG rules (see header)
LD_W_THRESHOLD  <- 0      # with min_n_loci_flag = 1 this flags every cluster
MIN_N_LOCI_FLAG <- 1      # -> no limit on what enters Stage-2 merging
MIN_N_LOCI_EMLG <- 3      # eMLG only for clusters with > 2 markers
SCORE_THRESHOLD <- 0.80
MIN_R2          <- 0.2
RHO_STAGE1      <- 0.5    # Stage-1 relative LD threshold
RHO_LDW         <- 0.95   # ld_w neighbourhood

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
cm_stamp <- function(cm) paste0("cM", sub("\\.", "", format(cm, trim = TRUE)))

## =========================================================================
## 1. read + translate to 012, split hybrids / parents
## =========================================================================
message("[di25-clust] reading ", INFILE)
d   <- fread(INFILE)
pos <- d$position
chr <- sub("chromosome_", "", d$chromosome)
G34 <- as.matrix(d[, !c("chromosome", "position"), with = FALSE])   # markers x individuals (1/2/3, 0=miss)
inds <- colnames(G34)

## 012 dosage: 1->0, 2->1, 3->2, 0(missing)->NA
G012 <- G34 - 1L
G012[G34 == 0L] <- NA_integer_

map <- data.table(Chr = paste0("Chr", chr), Pos = pos,
                  marker = paste0("Chr", chr, ":", pos))
setorder(map, Chr, Pos)
ord <- match(map$marker, paste0("Chr", chr, ":", pos))   # reorder G012 rows to map order
G012 <- G012[ord, ]
rownames(G012) <- map$marker                             # markers on rows, individuals on cols
colnames(G012) <- inds

is_parent <- grepl("^Faqu|^Fpol", inds)
GTs_hyb <- t(G012[, !is_parent, drop = FALSE])           # hybrids x markers
GTs_par <- t(G012[,  is_parent, drop = FALSE])           # parents x markers
message(sprintf("  %d markers | %d hybrids (LD) + %d parents (aside)",
                nrow(map), nrow(GTs_hyb), nrow(GTs_par)))

## save prepared inputs for the plotting step (avoids re-reading/translating)
saveRDS(list(map = map, GTs_hyb = GTs_hyb, GTs_par = GTs_par),
        file.path(OUTDIR, "di25_inputs.rds"))

## =========================================================================
## 2. LD decay (hybrids only) + ld_w  -- computed once, cached
## =========================================================================
DECAY_RDS <- file.path(OUTDIR, "di25_ld_decay.rds")
STAGE1_RDS <- file.path(OUTDIR, "di25_stage1.rds")

if (RECOMPUTE || !file.exists(STAGE1_RDS)) {
  message("[di25-clust] building gds + LD decay (hybrids only)")
  gds <- create_gds_from_geno(geno = GTs_hyb, map = map, gds_path = GDS_PATH)
  ld_decay <- compute_LD_decay(gds, keep_el = TRUE, slide = 100,
                               ld_method = "corr", cores = CORES)
  snpgdsClose(gds)
  saveRDS(ld_decay, DECAY_RDS)

  ## ld_w, matched to map$marker by name (compute_ld_w drops names)
  lw  <- compute_ld_w(ld_decay, rho = RHO_LDW, cores = CORES)
  ids <- unlist(lapply(ld_decay$by_chr, function(o) o$snp_ids), use.names = FALSE)
  map[, ld_w_095 := lw[match(marker, ids)]]

  message("[di25-clust] Stage 1 (rho = ", RHO_STAGE1, ")")
  stage1 <- ld_complexity_reduction(map = map, LD_decay = ld_decay,
                                    rho = RHO_STAGE1, cores = 1)
  saveRDS(list(stage1 = stage1, map = map), STAGE1_RDS)
} else {
  message("[di25-clust] loading cached Stage 1")
  tmp <- readRDS(STAGE1_RDS); stage1 <- tmp$stage1; map <- tmp$map
}

## =========================================================================
## 3. genetic map for the cM distance cap
## =========================================================================
rec_map <- fread(RECMAP)
rec_map[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
genetic_map <- rec_map[, .(Chr, Pos = pos, cM)]

## =========================================================================
## 4. Stage-2 merge + eMLG, swept over the cM cap
## =========================================================================
summ <- rbindlist(lapply(CM_SWEEP, function(cm) {
  message(sprintf("[di25-clust] Stage 2 @ %.1f cM", cm))
  res <- ld_prune_and_eMLG(
    GTs = GTs_hyb, stage1 = stage1, ld_w_col = "ld_w_095",
    ld_w_threshold = LD_W_THRESHOLD, min_n_loci_flag = MIN_N_LOCI_FLAG,
    score_threshold = SCORE_THRESHOLD, min_r2 = MIN_R2,
    genetic_map = genetic_map, cM_threshold = cm,
    compute_unflagged_eMLG = TRUE, min_n_loci_eMLG = MIN_N_LOCI_EMLG,
    cores = CORES
  )
  saveRDS(res, file.path(OUTDIR, sprintf("di25_clustering_%s.rds", cm_stamp(cm))))
  g <- res$groups
  data.table(cM = cm, n_clusters = nrow(g),
             n_eMLG = sum(g$has_eMLG), n_singleton = sum(g$n_loci == 1),
             n_2mark = sum(g$n_loci == 2), max_n_loci = max(g$n_loci),
             markers_in = sum(g$n_loci), reduction_pct = 100 * (1 - nrow(g) / sum(g$n_loci)))
}))
saveRDS(summ, file.path(OUTDIR, "di25_clustering_sweep_summary.rds"))
message("[di25-clust] sweep summary:")
print(summ)
message("[di25-clust] done")
