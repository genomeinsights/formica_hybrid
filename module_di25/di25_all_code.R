## ############################################################################
## module_di25 (high-DI analyses) -- ALL R code in one file (for sharing)
## Concatenated 2026-08-07. The ORIGINALS remain the source of truth in R/.
## Each block is preceded by a banner naming its original file.
## Order = logical reading order (shared core first, then pipeline).
## ############################################################################



## ============================================================================
## FILE: module_di25/R/diem_circos_core.R
## ============================================================================

## =========================================================
## module_di25 (high-DI analyses) -- DIEM circos core renderer (shared)
## =========================================================
## render_diem_circos(): draws a DIEM genotype circos by direct polar
## rasterisation. Used by both the per-SNP (diem_circos.R) and per-eMLG
## (diem_circos_eMLG.R) plots so they stay pixel-comparable.
##
## gt      : integer matrix, UNITS (markers/clusters) x INDIVIDUALS, already
##           ordered as it should appear -- rows in the drawing order around the
##           circle, columns from innermost ring to outermost. Codes:
##             0 = missing, 1 = hom F. aquilonia, 2 = het, 3 = hom F. polyctena
## chr_num : integer chromosome of each row of gt (same order as gt rows)
## Everything else is cosmetic (see defaults).

DIEM_COL    <- c("aqu" = "#440154",   # 1 hom F. aquilonia  (viridis purple)
                 "het" = "#21918C",   # 2 heterozygous      (viridis teal)
                 "pol" = "#FDE725")   # 3 hom F. polyctena  (viridis yellow)
DIEM_NA_COL <- "#D9D9D9"
DIEM_BG_COL <- "#FFFFFF"

## ---------------------------------------------------------------------------
## render_circos_raster(): general categorical circos by polar rasterisation.
## Same chromosome geometry as render_diem_circos(), but draws an arbitrary
## integer `code` matrix (units x rings) through an arbitrary `palette`
## (palette[v + 1] colours code value v). No dithering, no aqu/pol semantics.
## Rings are radial (column 1 = innermost); use it for per-unit track panels
## such as the tau sweep. `ring_labels` are drawn in the opening gap at each
## ring's radius; `legend_labels`/`legend_cols` add a corner legend.
## ---------------------------------------------------------------------------
render_circos_raster <- function(code, chr_num, palette, outpng = NULL, title = "",
                                 npx = 2600, r_in = 0.30, r_out = 0.92,
                                 start_deg = 90, clockwise = TRUE,
                                 open_deg = 24, chr_gap_deg = 0.5,
                                 chr_label_prefix = "Chr ", cex_main = 0.9,
                                 ring_labels = NULL, ring_sep = TRUE, ring_label_cex = 0.62,
                                 legend_labels = NULL, legend_cols = NULL,
                                 bg_col = "#FFFFFF", new_device = TRUE,
                                 add = FALSE, draw_chr_labels = TRUE) {
  storage.mode(code) <- "integer"
  nrings <- ncol(code)
  na_bg  <- length(bg_col) == 1L && is.na(bg_col)   # transparent background (overlay)

  chr_levels <- sort(unique(chr_num)); n_chr <- length(chr_levels)
  blk_start  <- vapply(chr_levels, function(k) min(which(chr_num == k)), integer(1))
  blk_n      <- vapply(chr_levels, function(k) sum(chr_num == k), integer(1))
  usable_deg <- 360 - open_deg - (n_chr - 1) * chr_gap_deg
  chr_deg    <- usable_deg * blk_n / sum(blk_n)
  chr_start  <- numeric(n_chr); chr_end <- numeric(n_chr); cur <- 0
  for (i in seq_len(n_chr)) {
    chr_start[i] <- cur; chr_end[i] <- cur + chr_deg[i]; cur <- chr_end[i] + chr_gap_deg
  }

  half <- npx / 2
  gx <- (col(matrix(0, npx, npx)) - 0.5 - half) / half
  gy <- (half - row(matrix(0, npx, npx)) + 0.5) / half
  rr <- sqrt(gx * gx + gy * gy)
  ang_raw <- atan2(gy, gx) * 180 / pi
  theta <- if (clockwise) (start_deg - ang_raw) %% 360 else (ang_raw - (start_deg - 90)) %% 360
  chr_of <- rep(NA_integer_, length(theta))
  for (i in seq_len(n_chr)) chr_of[theta >= chr_start[i] & theta < chr_end[i]] <- i

  cols <- rep(if (na_bg) NA_character_ else bg_col, npx * npx)
  inband <- rr >= r_in & rr <= r_out & !is.na(chr_of)
  rf   <- (rr[inband] - r_in) / (r_out - r_in) * nrings
  ring <- pmin(nrings - 1L, as.integer(floor(rf)))
  ci   <- chr_of[inband]
  loc  <- (theta[inband] - chr_start[ci]) / (chr_end[ci] - chr_start[ci])
  mrow <- blk_start[ci] + pmin(blk_n[ci] - 1L, as.integer(floor(loc * blk_n[ci])))
  val  <- code[mrow + ring * nrow(code)]; val[is.na(val)] <- 0L
  cell_cols <- palette[val + 1L]
  if (ring_sep) {                                   # thin bg gaps between rings
    within <- rf - floor(rf)
    cell_cols[within < 0.02 | within > 0.98] <- if (na_bg) NA_character_ else bg_col
  }
  cols[which(inband)] <- cell_cols
  ras <- as.raster(matrix(cols, npx, npx, byrow = FALSE))

  if (!add) {
    if (new_device) png(outpng, width = npx, height = npx, res = 300)
    op <- par(mar = c(0, 0, 2, 0), xpd = NA)
    on.exit({ par(op); if (new_device) dev.off() }, add = TRUE)
    plot.new(); plot.window(xlim = c(-1, 1), ylim = c(-1, 1), asp = 1)
  } else {                                          # overlay into the caller's plot
    op <- par(xpd = NA); on.exit(par(op), add = TRUE)
  }
  rasterImage(ras, -1, -1, 1, 1, interpolate = FALSE)

  ## chromosome labels outside the rings
  lab_r <- r_out + 0.03
  if (draw_chr_labels) for (i in seq_len(n_chr)) {
    mid <- (chr_start[i] + chr_end[i]) / 2
    a <- (start_deg - mid) * pi / 180
    x <- lab_r * cos(a); y <- lab_r * sin(a)
    deg <- atan2(y, x) * 180 / pi; flip <- deg > 90 || deg < -90
    text(x, y, paste0(chr_label_prefix, chr_levels[i]), cex = 0.5, col = "grey20",
         srt = if (flip) deg + 180 else deg, adj = if (flip) c(1, 0.5) else c(0, 0.5))
  }
  ## ring labels in the opening gap, one at each ring's mid-radius
  if (!is.null(ring_labels)) {
    theta_lab <- 360 - open_deg / 2
    a <- (start_deg - theta_lab) * pi / 180
    for (i in seq_len(nrings)) {
      rm <- r_in + (i - 0.5) / nrings * (r_out - r_in)
      text(rm * cos(a), rm * sin(a), ring_labels[i], cex = ring_label_cex, col = "grey15", adj = c(0.5, 0.5))
    }
  }
  if (!is.null(legend_labels))
    legend("bottomleft", legend = legend_labels, fill = legend_cols, bty = "n",
           cex = 0.75, border = NA, inset = c(0.02, 0.02))
  if (nzchar(title)) title(main = title, cex.main = cex_main, line = 0.2)
  invisible(list(n_units = nrow(code), n_rings = nrings))
}

render_diem_circos <- function(gt, chr_num, outpng = NULL, title = "",
                               col = DIEM_COL, na_col = DIEM_NA_COL,
                               bg_col = DIEM_BG_COL, dither = TRUE,
                               npx = 2600, r_in = 0.30, r_out = 0.92,
                               start_deg = 90, clockwise = TRUE,
                               open_deg = 12, chr_gap_deg = 0.5,
                               seed = 1, chr_label_prefix = "Chr ",
                               cex_main = 0.85, new_device = TRUE,
                               draw_chr_labels = TRUE) {
  set.seed(seed)
  storage.mode(gt) <- "integer"

  ## dither missing -> random genotype state (cosmetic; keeps the speckle)
  n_miss <- sum(gt == 0L, na.rm = TRUE)
  if (dither) gt[which(gt == 0L)] <- sample.int(3L, n_miss, replace = TRUE)

  ## chromosome blocks in the (already sorted) row order
  chr_levels <- sort(unique(chr_num))
  n_chr      <- length(chr_levels)
  blk_start  <- vapply(chr_levels, function(k) min(which(chr_num == k)), integer(1))
  blk_n      <- vapply(chr_levels, function(k) sum(chr_num == k), integer(1))

  ## angular layout: arc per chromosome proportional to its unit count
  usable_deg <- 360 - open_deg - (n_chr - 1) * chr_gap_deg
  chr_deg    <- usable_deg * blk_n / sum(blk_n)
  chr_start  <- numeric(n_chr); chr_end <- numeric(n_chr); cur <- 0
  for (i in seq_len(n_chr)) {
    chr_start[i] <- cur; chr_end[i] <- cur + chr_deg[i]; cur <- chr_end[i] + chr_gap_deg
  }

  ## pixel grid in fraction-of-half-canvas coordinates
  half <- npx / 2
  gx <- (col(matrix(0, npx, npx)) - 0.5 - half) / half
  gy <- (half - row(matrix(0, npx, npx)) + 0.5) / half
  rr <- sqrt(gx * gx + gy * gy)
  ang_raw <- atan2(gy, gx) * 180 / pi
  theta <- if (clockwise) (start_deg - ang_raw) %% 360 else (ang_raw - (start_deg - 90)) %% 360

  ## which chromosome each pixel falls in (NA in gaps / opening)
  chr_of <- rep(NA_integer_, length(theta))
  for (i in seq_len(n_chr)) chr_of[theta >= chr_start[i] & theta < chr_end[i]] <- i

  ## fill genotype annulus
  idx    <- matrix(0L, npx, npx)
  inband <- rr >= r_in & rr <= r_out & !is.na(chr_of)
  ni     <- ncol(gt)
  ring   <- pmin(ni - 1L, as.integer(floor((rr[inband] - r_in) / (r_out - r_in) * ni)))
  ci     <- chr_of[inband]
  loc    <- (theta[inband] - chr_start[ci]) / (chr_end[ci] - chr_start[ci])
  mrow   <- blk_start[ci] + pmin(blk_n[ci] - 1L, as.integer(floor(loc * blk_n[ci])))
  val    <- gt[mrow + ring * nrow(gt)]
  val[is.na(val)] <- 0L
  idx[inband] <- val

  pal  <- c(bg_col, col[["aqu"]], col[["het"]], col[["pol"]])
  ras  <- as.raster(matrix(pal[idx + 1L], npx, npx, byrow = FALSE))

  ## render (own PNG device, or draw into the caller's current device/panel)
  if (new_device) png(outpng, width = npx, height = npx, res = 300)
  op <- par(mar = c(0, 0, 2, 0), xpd = NA)
  on.exit({ par(op); if (new_device) dev.off() }, add = TRUE)
  plot.new(); plot.window(xlim = c(-1, 1), ylim = c(-1, 1), asp = 1)
  rasterImage(ras, -1, -1, 1, 1, interpolate = FALSE)

  lab_r <- r_out + 0.03
  if (draw_chr_labels) for (i in seq_len(n_chr)) {
    mid <- (chr_start[i] + chr_end[i]) / 2
    a   <- (start_deg - mid) * pi / 180
    x <- lab_r * cos(a); y <- lab_r * sin(a)
    deg <- atan2(y, x) * 180 / pi
    flip <- deg > 90 || deg < -90
    text(x, y, paste0(chr_label_prefix, chr_levels[i]), cex = 0.5, col = "grey20",
         srt = if (flip) deg + 180 else deg, adj = if (flip) c(1, 0.5) else c(0, 0.5))
  }
  if (nzchar(title)) title(main = title, cex.main = cex_main, line = 0.2)
  invisible(list(n_miss = n_miss, n_units = nrow(gt), n_ind = ncol(gt)))
}


## ============================================================================
## FILE: module_di25/R/diem_circos.R
## ============================================================================

## =========================================================
## module_di25 (high-DI analyses) -- DIEM circos plot, per-SNP (R recreation)
## =========================================================
## Recreates, in R, the DIEM genotype circos originally produced in Python,
## at the per-SNP level, so extra tracks/annotation can be added later.
##
## Input : data/species_diagnostic_markers_DI25_20pops.tsv.gz
##   195 individual genotype columns + `chromosome` + `position`.
##   Genotypes coded relative to the two parental species:
##     0 = missing, 1 = hom F. aquilonia, 2 = het, 3 = hom F. polyctena.
##   (Confirmed: the 15 Faqu parent columns ~all 1, the Fpol parents ~all 3.)
##
## Layout: one sector per chromosome (arc proportional to marker count); each
##   individual is a concentric ring, ordered by genome-wide hybrid index so the
##   innermost rings are the most F. aquilonia and the outermost the most
##   F. polyctena. Rendering is done by render_diem_circos() (diem_circos_core.R).
##
## Run from the repo root:  Rscript module_di25/R/diem_circos.R
## =========================================================

suppressMessages(library(data.table))
source("module_di25/R/diem_circos_core.R")

INFILE   <- "data/species_diagnostic_markers_DI25_20pops.tsv.gz"
OUTPNG   <- "module_di25/Figures/diem_circos.png"
DI_LABEL <- -25

message("[diem-circos SNP] reading ", INFILE)
d   <- fread(INFILE)
pos <- d$position
gt  <- as.matrix(d[, !c("chromosome", "position"), with = FALSE])   # markers x individuals
storage.mode(gt) <- "integer"
chr_num <- as.integer(sub("chromosome_", "", d$chromosome))
message(sprintf("  %d markers x %d individuals", nrow(gt), ncol(gt)))

## order markers by chromosome then position
ord_m <- order(chr_num, pos)
gt <- gt[ord_m, ]; chr_num <- chr_num[ord_m]

## order individuals by hybrid index (0 = aquilonia .. 1 = polyctena)
hi <- apply(gt, 2, function(g) { g <- g[g > 0]; if (!length(g)) NA_real_ else mean((g - 1) / 2) })
gt <- gt[, order(hi)]

n_miss <- sum(gt == 0L)
ttl <- sprintf("Hybrids masked @ DI = %d: %s SNPs, %s missing genotypes (%.1f%%) dithered",
               DI_LABEL, format(nrow(gt), big.mark = ","),
               format(n_miss, big.mark = ","), 100 * n_miss / length(gt))

message("[diem-circos SNP] rendering -> ", OUTPNG)
render_diem_circos(gt, chr_num, OUTPNG, title = ttl)
message("[diem-circos SNP] done")


## ============================================================================
## FILE: module_di25/R/diem_circos_eMLG.R
## ============================================================================

## =========================================================
## module_di25 (high-DI analyses) -- DIEM circos plot, per-eMLG (LD-reduced)
## =========================================================
## The LD-reduced counterpart of diem_circos.R: instead of one ring-cell per
## diagnostic SNP, each cell is an LD-cluster eMLG consensus. This collapses
## linked diagnostic markers into their cluster and shows the same picture at
## the level of (approximately) independent units.
##
## Mapping LD-clusters to the DI25 data:
##   * eMLG_sorted_cM05.rds        -- $eMLG (164 hybrids x 96,461 clusters, dosage
##                                    0..2) + $groups (representative, members).
##   * eMLG_sorted_cM05_parents.rds -- 30 parent individuals x the SAME 96,461
##                                    clusters. rbind -> 194 individuals.
##   * A cluster is kept if any of its member markers is a DI25 diagnostic marker
##     (paste0("Chr", chr, ":", pos) from the SNP file). -> ~5,178 clusters.
##
## Orientation: the eMLG consensus is oriented to an arbitrary within-cluster
##   reference allele, so each cluster is re-oriented to F. aquilonia using the
##   parent side (mean Faqu dosage >= mean Fpol dosage => keep, else 2 - x).
##   Oriented dosage a in [0,2] with 2 = aquilonia; rounded and recoded to the
##   shared DIEM scheme 1=aqu, 2=het, 3=pol via 3 - round(a).
##
## Run from the repo root:  Rscript module_di25/R/diem_circos_eMLG.R
## =========================================================

suppressMessages(library(data.table))
source("module_di25/R/diem_circos_core.R")

SORTED_H <- "moduleA_sorting/data/eMLG_sorted_cM05.rds"
SORTED_P <- "moduleA_sorting/data/eMLG_sorted_cM05_parents.rds"
DI25     <- "data/species_diagnostic_markers_DI25_20pops.tsv.gz"
OUTPNG   <- "module_di25/Figures/diem_circos_eMLG.png"
DI_LABEL <- -25

## ---- load eMLG consensus (hybrid + parent) ------------------------------
message("[diem-circos eMLG] loading eMLG consensus")
s    <- readRDS(SORTED_H)
H    <- s$eMLG                                   # 164 hybrids x clusters
g    <- s$groups
P    <- readRDS(SORTED_P)                        # 30 parents  x clusters
stopifnot(identical(colnames(H), colnames(P)))
E    <- rbind(H, P)                              # 194 individuals x clusters

## ---- keep clusters that carry a DI25 diagnostic marker ------------------
d      <- fread(DI25, select = c("chromosome", "position"))
di_ids <- paste0("Chr", sub("chromosome_", "", d$chromosome), ":", d$position)
mm     <- data.table(group_id = rep(g$group_id, lengths(g$members)),
                     marker   = unlist(g$members))
keep   <- unique(mm[marker %in% di_ids, group_id])
gk     <- g[match(keep, group_id)]
E      <- E[, keep, drop = FALSE]
message(sprintf("  %d clusters carry a DI25 marker (from %d hybrids + %d parents)",
                length(keep), nrow(H), nrow(P)))

## ---- orient each cluster to F. aquilonia via the parent side ------------
faqu <- grep("^Faqu", rownames(E)); fpol <- grep("^Fpol", rownames(E))
maqu <- colMeans(E[faqu, , drop = FALSE], na.rm = TRUE)
mpol <- colMeans(E[fpol, , drop = FALSE], na.rm = TRUE)
flip <- which(maqu < mpol)
E[, flip] <- 2 - E[, flip]                       # now 2 = aquilonia allele everywhere

## ---- recode to DIEM scheme (1 aqu, 2 het, 3 pol; 0 missing) --------------
code <- 3L - round(E)                            # 2->1 aqu, 1->2 het, 0->3 pol
code[is.na(code)] <- 0L
storage.mode(code) <- "integer"

## ---- positions: representative marker of each kept cluster --------------
rep_chr <- as.integer(sub("Chr", "", sub(":.*", "", gk$representative)))
rep_pos <- as.integer(sub(".*:", "", gk$representative))

## units x individuals, ordered by chr/pos (rows) and hybrid index (cols)
gt <- t(code)                                    # clusters x individuals
ord_u <- order(rep_chr, rep_pos)
gt <- gt[ord_u, ]; chr_num <- rep_chr[ord_u]

## hybrid index on the shared 1=aqu..3=pol code: 0 = aquilonia, 1 = polyctena
hi <- apply(gt, 2, function(x) { x <- x[x > 0]; if (!length(x)) NA_real_ else mean((x - 1) / 2) })
gt <- gt[, order(hi)]                            # inner = most aquilonia (matches SNP plot)

n_miss <- sum(gt == 0L)
ttl <- sprintf("LD-reduced (eMLG) @ DI = %d: %s clusters, %s missing (%.1f%%) dithered",
               DI_LABEL, format(nrow(gt), big.mark = ","),
               format(n_miss, big.mark = ","), 100 * n_miss / length(gt))

message("[diem-circos eMLG] rendering -> ", OUTPNG)
render_diem_circos(gt, chr_num, OUTPNG, title = ttl)
message("[diem-circos eMLG] done")


## ============================================================================
## FILE: module_di25/R/di25_ld_clustering.R
## ============================================================================

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


## ============================================================================
## FILE: module_di25/R/diem_circos_di25_eMLG.R
## ============================================================================

## =========================================================
## module_di25 (high-DI analyses) -- DIEM circos on the from-scratch DI25 clustering
## =========================================================
## LD-reduced DIEM circos built on the high-DI-ONLY clustering produced by
## di25_ld_clustering.R (5 cM cap), plotting ALL units:
##   * clusters with > 2 markers  -> eMLG consensus dosage
##   * clusters with 1-2 markers  -> representative SNP genotype
##
## Unlike diem_circos_eMLG.R (which used the canonical all-marker eMLGs), every
## unit here summarises high-DI variation only.
##
## Consensus is computed over the COMBINED 165 hybrids + 30 parents in one pass
## (one polarization per unit), so hybrid and parent sides share a sign
## convention; each unit is then oriented to F. aquilonia via the parent rows.
##
## Run from the repo root:  Rscript module_di25/R/diem_circos_di25_eMLG.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")                 # consensus_dosage(), expected_gt_dosage()
source("module_di25/R/diem_circos_core.R")

CM_STAMP <- "cM5"
CLUST    <- sprintf("module_di25/data/di25_clustering_%s.rds", CM_STAMP)
INPUTS   <- "module_di25/data/di25_inputs.rds"
OUTPNG   <- "module_di25/Figures/diem_circos_di25_eMLG_cM5.png"
DI_LABEL <- -25

## ---- load clustering + prepared 012 genotypes ---------------------------
message("[di25-circos] loading ", CLUST)
res <- readRDS(CLUST); g <- res$groups
inp <- readRDS(INPUTS)
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)             # 195 individuals x markers (012, marker colnames)
faqu <- grep("^Faqu", rownames(GTs_all)); fpol <- grep("^Fpol", rownames(GTs_all))
message(sprintf("  %d units | %d individuals (%d hybrids + %d parents)",
                nrow(g), nrow(GTs_all), nrow(inp$GTs_hyb), nrow(inp$GTs_par)))

## ---- one dosage vector per unit (eMLG consensus or representative SNP) ---
is_emlg <- g$n_loci > 2
message(sprintf("  %d eMLG units (>2 markers) + %d representative-SNP units",
                sum(is_emlg), sum(!is_emlg)))
D <- vapply(seq_len(nrow(g)), function(i) {
  if (is_emlg[i]) consensus_dosage(GTs_all, g$members[[i]])   # 0..2, polarized
  else            GTs_all[, g$representative[i]]               # 0/1/2 representative SNP
}, numeric(nrow(GTs_all)))                                     # individuals x units
D <- t(D)                                                      # units x individuals

## ---- orient each unit to F. aquilonia via the parent side ---------------
maqu <- rowMeans(D[, faqu, drop = FALSE], na.rm = TRUE)
mpol <- rowMeans(D[, fpol, drop = FALSE], na.rm = TRUE)
flip <- which(maqu < mpol)
D[flip, ] <- 2 - D[flip, ]                                     # 2 = aquilonia everywhere

## ---- recode to shared DIEM scheme (1 aqu, 2 het, 3 pol; 0 missing) -------
code <- 3L - round(D)                                          # 2->1 aqu, 1->2 het, 0->3 pol
code[is.na(code)] <- 0L
storage.mode(code) <- "integer"

## ---- position by representative marker; order units & individuals -------
rep_chr <- as.integer(sub("Chr", "", sub(":.*", "", g$representative)))
rep_pos <- as.integer(sub(".*:", "", g$representative))
ord_u   <- order(rep_chr, rep_pos)
code    <- code[ord_u, ]; chr_num <- rep_chr[ord_u]

hi <- apply(code, 2, function(x) { x <- x[x > 0]; if (!length(x)) NA_real_ else mean((x - 1) / 2) })
code <- code[, order(hi)]                                      # inner = most aquilonia

## ---- render -------------------------------------------------------------
n_miss <- sum(code == 0L)
ttl <- sprintf("LD-reduced (DI = %d, %s): %s units [%s eMLG + %s rep-SNP], %s missing (%.1f%%)",
               DI_LABEL, sub("cM", "", CM_STAMP) |> (\(x) paste0(x, " cM"))(),
               format(nrow(code), big.mark = ","),
               format(sum(is_emlg), big.mark = ","), format(sum(!is_emlg), big.mark = ","),
               format(n_miss, big.mark = ","), 100 * n_miss / length(code))

message("[di25-circos] rendering -> ", OUTPNG)
render_diem_circos(code, chr_num, OUTPNG, title = ttl, cex_main = 0.8)
message("[di25-circos] done")


## ============================================================================
## FILE: module_di25/R/diem_circos_compare.R
## ============================================================================

## =========================================================
## module_di25 (high-DI analyses) -- per-SNP vs LD-reduced DIEM circos (side by side)
## =========================================================
## One figure, two panels on the SAME DI25 markers and the SAME 195 individuals:
##   (a) per-SNP        -- 51,612 diagnostic SNPs        (diem_circos.R)
##   (b) LD-reduced     -- 11,052 high-DI-only units, 5cM (diem_circos_di25_eMLG.R)
## Individuals are ordered ONCE (by the per-SNP hybrid index) and that order is
## applied to both panels, so a ring at a given radius is the same individual in
## both -- the LD reduction is the only thing that changes between panels.
##
## Run from the repo root:  Rscript module_di25/R/diem_circos_compare.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")
source("module_di25/R/diem_circos_core.R")

TSV      <- "data/species_diagnostic_markers_DI25_20pops.tsv.gz"
CLUST    <- "module_di25/data/di25_clustering_cM5.rds"
INPUTS   <- "module_di25/data/di25_inputs.rds"
OUTPNG   <- "module_di25/Figures/diem_circos_compare_snp_vs_ldreduced.png"

## ---- (a) per-SNP: markers x individuals, code 1/2/3 (0 miss) -------------
message("[compare] per-SNP panel")
d <- fread(TSV)
snp_chr <- as.integer(sub("chromosome_", "", d$chromosome))
snp_gt  <- as.matrix(d[, !c("chromosome", "position"), with = FALSE])   # 1/2/3, 0=miss
storage.mode(snp_gt) <- "integer"
ord_m   <- order(snp_chr, d$position)
snp_gt  <- snp_gt[ord_m, ]; snp_chr <- snp_chr[ord_m]
snp_inds <- colnames(snp_gt)

## ---- (b) LD-reduced: units x individuals, code 1/2/3 --------------------
message("[compare] LD-reduced panel")
res <- readRDS(CLUST); g <- res$groups
inp <- readRDS(INPUTS)
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)                              # 195 x markers (012)
faqu <- grep("^Faqu", rownames(GTs_all)); fpol <- grep("^Fpol", rownames(GTs_all))
is_emlg <- g$n_loci > 2
D <- vapply(seq_len(nrow(g)), function(i) {
  if (is_emlg[i]) consensus_dosage(GTs_all, g$members[[i]]) else GTs_all[, g$representative[i]]
}, numeric(nrow(GTs_all)))
D <- t(D)                                                              # units x individuals
flip <- which(rowMeans(D[, faqu, drop = FALSE], na.rm = TRUE) <
              rowMeans(D[, fpol, drop = FALSE], na.rm = TRUE))
D[flip, ] <- 2 - D[flip, ]                                             # 2 = aquilonia
emlg_code <- 3L - round(D); emlg_code[is.na(emlg_code)] <- 0L
storage.mode(emlg_code) <- "integer"
colnames(emlg_code) <- rownames(GTs_all)
rep_chr <- as.integer(sub("Chr", "", sub(":.*", "", g$representative)))
rep_pos <- as.integer(sub(".*:", "", g$representative))
ord_u   <- order(rep_chr, rep_pos)
emlg_code <- emlg_code[ord_u, ]; emlg_chr <- rep_chr[ord_u]

## ---- shared individual order (by per-SNP hybrid index) ------------------
hi <- apply(snp_gt, 2, function(x) { x <- x[x > 0]; if (!length(x)) NA_real_ else mean((x - 1) / 2) })
ord_names <- snp_inds[order(hi)]                                       # inner = most aquilonia
snp_gt    <- snp_gt[, ord_names]
emlg_code <- emlg_code[, match(ord_names, colnames(emlg_code))]

## ---- draw both panels into one wide canvas ------------------------------
message("[compare] rendering -> ", OUTPNG)
png(OUTPNG, width = 5200, height = 2750, res = 300)
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
render_diem_circos(snp_gt, snp_chr, new_device = FALSE,
                   title = sprintf("a   Per-SNP: %s diagnostic SNPs", format(nrow(snp_gt), big.mark = ",")),
                   cex_main = 1.0)
render_diem_circos(emlg_code, emlg_chr, new_device = FALSE,
                   title = sprintf("b   LD-reduced: %s units (5 cM) = %s eMLG + %s rep-SNP",
                                   format(nrow(emlg_code), big.mark = ","),
                                   format(sum(is_emlg), big.mark = ","),
                                   format(sum(!is_emlg), big.mark = ",")),
                   cex_main = 1.0)
mtext("DIEM ancestry (purple = F. aquilonia, teal = het, yellow = F. polyctena)  |  rings inner to outer: most aquilonia to most polyctena",
      outer = TRUE, cex = 0.8, line = 0.3)
dev.off()
message("[compare] done")


## ============================================================================
## FILE: module_di25/R/di25_sorting.R
## ============================================================================

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


## ============================================================================
## FILE: module_di25/R/diem_circos_sorting_sweep.R
## ============================================================================

## =========================================================
## module_di25 (high-DI analyses) -- ancestry-sorting tau-sweep circos
## =========================================================
## The tau sweep drawn AROUND the genome, in the same circular chromosome layout
## as the DIEM circos -- but the concentric rings are the four sorting thresholds
## tau = 0.5, 0.6, 0.7, 0.8 (inner -> outer), and each unit is coloured by its
## sort class at that tau. phi = 0.85 fixed. No DIEM genotype data here.
##
##   colour: toward F. aquilonia (purple) / toward F. polyctena (yellow) /
##           direction unresolved (teal) / unsorted (near-white)
##
## As tau increases outward, colour thins out (fewer units sorted); direction
## colour shows where each ancestry prevails and how the balance shifts with tau.
##
## LEVEL = "SNP" (51,612 markers, mirrors the flagship DIEM circos) or
##         "eMLG" (11,052 LD-reduced units).
##
## Run from the repo root:  Rscript module_di25/R/diem_circos_sorting_sweep.R
## =========================================================

suppressMessages(library(data.table))
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()
source("module_di25/R/diem_circos_core.R")  # render_circos_raster()

LEVEL    <- Sys.getenv("SWEEP_LEVEL", "SNP")         # "SNP" or "eMLG" (env-overridable)
TAU_GRID <- c(0.5, 0.6, 0.7, 0.8)
OUTPNG   <- sprintf("module_di25/Figures/diem_circos_sorting_sweep_%s.png", tolower(LEVEL))

## sort-class palette (code 0 unsorted .. 3 unresolved)
SORT_PAL <- c("#F4F4F4",   # 0 unsorted / not differentiated
              "#440154",   # 1 toward F. aquilonia (purple)
              "#FDE725",   # 2 toward F. polyctena (yellow)
              "#35B7AE")   # 3 direction unresolved (teal)
CLS_CODE <- c(unsorted = 0L, aquilonia = 1L, polyctena = 2L, unresolved = 3L, ambiguous = 3L)

## ---- load the per-unit sorting counts + genomic positions ---------------
if (LEVEL == "SNP") {
  ps <- readRDS("module_di25/data/di25_sorting_snp.rds")
  chr_num <- as.integer(sub("Chr", "", sub(":.*", "", ps$marker)))
  pos     <- as.integer(sub(".*:", "", ps$marker))
} else {
  ps <- readRDS("module_di25/data/di25_sorting_emlg.rds")
  g  <- readRDS("module_di25/data/di25_clustering_cM5.rds")$groups
  rep_mk  <- g$representative[match(ps$group_id, g$group_id)]
  chr_num <- as.integer(sub("Chr", "", sub(":.*", "", rep_mk)))
  pos     <- as.integer(sub(".*:", "", rep_mk))
}

## order units by chromosome then position
ord <- order(chr_num, pos); ps <- ps[ord]; chr_num <- chr_num[ord]

## ---- classify every unit at each tau -> code matrix (units x rings) ------
ok <- ps$differentiated & ps$n_obs > 0 & is.finite(ps$uni_score)
code <- vapply(TAU_GRID, function(tau) {
  v <- integer(nrow(ps))                                   # 0 = unsorted / not differentiated
  cl <- classify_sort(ps$n_aqu[ok], ps$n_pol[ok], ps$n_obs[ok],
                      sort_th = tau, sort_rule = "binom", alpha = 0.05)
  v[which(ok)] <- CLS_CODE[cl]
  v
}, integer(nrow(ps)))                                      # units x length(TAU_GRID)

## ---- render -------------------------------------------------------------
n_lab <- if (LEVEL == "SNP") "diagnostic SNPs" else "LD-reduced units"
ttl <- sprintf("Ancestry-sorting tau sweep (phi = 0.85) -- %s %s",
               format(nrow(ps), big.mark = ","), n_lab)
message("[sorting-sweep] rendering -> ", OUTPNG)
render_circos_raster(
  code, chr_num, palette = SORT_PAL, outpng = OUTPNG, title = ttl,
  ring_labels = sprintf("tau=%.1f", TAU_GRID),
  legend_labels = c("toward F. aquilonia", "toward F. polyctena",
                    "direction unresolved", "unsorted"),
  legend_cols = SORT_PAL[c(2, 3, 4, 1)]
)
message("[sorting-sweep] done (inner->outer rings: tau ",
        paste(TAU_GRID, collapse = ", "), ")")


## ============================================================================
## FILE: module_di25/R/diem_circos_di25_sortclass.R
## ============================================================================

## =========================================================
## module_di25 (high-DI analyses) -- LD-reduced DIEM circos + sort-class outer track
## =========================================================
## The 5 cM LD-reduced DIEM circos (per-eMLG genotype rings), with an OUTER
## annotation ring giving each unit's ancestry-sorting class at tau = 0.6
## (phi = 0.85). Genotype rings and the sort track share the SAME unit order
## (representative chr/pos), so the outer class aligns radially with the unit
## below it.
##
## Run from the repo root:  Rscript module_di25/R/diem_circos_di25_sortclass.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()
source("module_di25/R/diem_circos_core.R")

TAU_TRACK <- 0.6
OUTPNG    <- "module_di25/Figures/diem_circos_di25_eMLG_sortclass_cM5.png"

## sort-class palette (matches the tau-sweep figure)
SORT_PAL <- c("#F4F4F4", "#440154", "#FDE725", "#35B7AE")   # unsorted/aqu/pol/unresolved
CLS_CODE <- c(unsorted = 0L, aquilonia = 1L, polyctena = 2L, unresolved = 3L, ambiguous = 3L)

## ---- inputs (clustering, prepared genotypes, sorting counts) -------------
g   <- readRDS("module_di25/data/di25_clustering_cM5.rds")$groups
inp <- readRDS("module_di25/data/di25_inputs.rds")
ps  <- readRDS("module_di25/data/di25_sorting_emlg.rds")
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)
faqu <- grep("^Faqu", rownames(GTs_all)); fpol <- grep("^Fpol", rownames(GTs_all))

## shared unit order: representative chr/pos
rep_chr <- as.integer(sub("Chr", "", sub(":.*", "", g$representative)))
rep_pos <- as.integer(sub(".*:", "", g$representative))
ord_u   <- order(rep_chr, rep_pos)
chr_num <- rep_chr[ord_u]

## ---- genotype rings (per-eMLG, oriented to aquilonia) -------------------
is_emlg <- g$n_loci > 2
D <- vapply(seq_len(nrow(g)), function(i) {
  if (is_emlg[i]) consensus_dosage(GTs_all, g$members[[i]]) else GTs_all[, g$representative[i]]
}, numeric(nrow(GTs_all)))
D <- t(D)                                                    # units x individuals
flip <- which(rowMeans(D[, faqu, drop = FALSE], na.rm = TRUE) <
              rowMeans(D[, fpol, drop = FALSE], na.rm = TRUE))
D[flip, ] <- 2 - D[flip, ]
code_geno <- 3L - round(D); code_geno[is.na(code_geno)] <- 0L
storage.mode(code_geno) <- "integer"
colnames(code_geno) <- rownames(GTs_all)
code_geno <- code_geno[ord_u, ]                              # unit order
hi <- apply(code_geno, 2, function(x) { x <- x[x > 0]; if (!length(x)) NA_real_ else mean((x - 1) / 2) })
code_geno <- code_geno[, order(hi)]                          # inner ring = most aquilonia

## ---- sort-class outer track (same unit order) --------------------------
psm <- ps[match(g$group_id, group_id)]                       # align ps to g order
ok  <- psm$differentiated & psm$n_obs > 0 & is.finite(psm$uni_score)
track <- integer(nrow(g))
cl <- classify_sort(psm$n_aqu[ok], psm$n_pol[ok], psm$n_obs[ok],
                    sort_th = TAU_TRACK, sort_rule = "binom", alpha = 0.05)
track[which(ok)] <- CLS_CODE[cl]
track <- matrix(track[ord_u], ncol = 1)                      # units x 1, unit order
storage.mode(track) <- "integer"

## ---- draw: genotype inside, sort-class ring outside --------------------
message("[sortclass] rendering -> ", OUTPNG)
png(OUTPNG, width = 2600, height = 2600, res = 300)
render_diem_circos(code_geno, chr_num, new_device = FALSE, draw_chr_labels = FALSE,
                   r_in = 0.30, r_out = 0.80, cex_main = 0.85,
                   title = sprintf("LD-reduced DIEM (5 cM) + sort class @ tau = %.1f (phi = 0.85)  |  %s units",
                                   TAU_TRACK, format(nrow(g), big.mark = ",")))
render_circos_raster(track, chr_num, palette = SORT_PAL, add = TRUE, bg_col = NA,
                     r_in = 0.825, r_out = 0.90, ring_sep = FALSE, draw_chr_labels = TRUE,
                     ring_labels = sprintf("sort tau=%.1f", TAU_TRACK),
                     legend_labels = c("toward F. aquilonia", "toward F. polyctena",
                                       "direction unresolved", "unsorted"),
                     legend_cols = SORT_PAL[c(2, 3, 4, 1)])
dev.off()
message("[sortclass] done")


## ============================================================================
## FILE: module_di25/R/diem_circos_population.R
## ============================================================================

## =========================================================
## module_di25 (high-DI analyses) -- per-POPULATION DIEM circos (both levels)
## =========================================================
## The DIEM circos summarised per population: the 195 individual rings collapse
## to 20 hybrid-population rings, each cell = that population's MEAN F. aquilonia-
## allele frequency at the unit (continuous). Rendered at both levels:
##   SNP  -- 51,612 diagnostic markers
##   eMLG -- 11,052 LD-reduced units (5 cM)
##
## Colour (viridis): purple = population near-fixed for F. aquilonia (freq ~1),
## yellow = near-fixed for F. polyctena (freq ~0), teal ~ balanced. Continuous
## frequency is rendered by binning into a 100-step palette via render_circos_raster().
## Rings inner -> outer are ordered by the population's overall aquilonia ancestry.
##
## Run from the repo root:  Rscript module_di25/R/diem_circos_population.R
## =========================================================

suppressMessages({ library(data.table); library(viridisLite) })
devtools::load_all("~/gitlab/LDscnR/")
source("module_di25/R/diem_circos_core.R")

## continuous viridis palette: code 0 = missing (grey); codes 1..100 = freq 1..0
PAL <- c("#D9D9D9", viridis(100))              # palette[code+1]; code1=purple(freq1)..code100=yellow(freq0)
freq_to_code <- function(f) { v <- 1L + as.integer(round((1 - f) * 99)); v[is.na(f)] <- 0L; v }

## orient a dosage matrix (units x individuals) so 2 = F. aquilonia allele
orient_aqu <- function(M, faqu, fpol) {
  flip <- which(rowMeans(M[, faqu, drop = FALSE], na.rm = TRUE) <
                rowMeans(M[, fpol, drop = FALSE], na.rm = TRUE))
  M[flip, ] <- 2 - M[flip, ]
  M
}

## ---- inputs -------------------------------------------------------------
inp <- readRDS("module_di25/data/di25_inputs.rds")
e2  <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd  <- e2$sample_data_with_parents
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)                       # 195 x markers
pops    <- sd$Population[match(rownames(GTs_all), sd$Sample_ID)] # NA for the unmatched hybrid
faqu <- grep("^Faqu", rownames(GTs_all)); fpol <- grep("^Fpol", rownames(GTs_all))
hybrid_pops <- setdiff(unique(pops[!is.na(pops)]), c("aquilonia_parent", "polyctena_parent"))

## population-mean aquilonia frequency (units x 20 pops) from an oriented matrix
pop_freq <- function(a) {
  vapply(hybrid_pops, function(p) rowMeans(a[, which(pops == p), drop = FALSE], na.rm = TRUE) / 2,
         numeric(nrow(a)))                                       # units x pops
}

render_pop <- function(a_units, chr_num, level, n_lab) {
  F <- pop_freq(a_units)                                         # units x 20
  ord_p <- order(colMeans(F, na.rm = TRUE), decreasing = TRUE)  # inner ring = most aquilonia
  F <- F[, ord_p]; pop_order <- hybrid_pops[ord_p]
  code <- matrix(freq_to_code(F), nrow = nrow(F))
  outpng <- sprintf("module_di25/Figures/diem_circos_population_%s.png", tolower(level))
  message("[pop-circos] rendering ", level, " -> ", outpng)
  render_circos_raster(
    code, chr_num, palette = PAL, outpng = outpng, ring_sep = FALSE,
    title = sprintf("Per-population ancestry (%s): mean F. aquilonia-allele freq  |  %s %s",
                    level, format(nrow(code), big.mark = ","), n_lab),
    ring_labels = pop_order, ring_label_cex = 0.34, open_deg = 30,
    legend_labels = c("F. aquilonia (1.0)", "balanced (0.5)", "F. polyctena (0.0)"),
    legend_cols = c(viridis(100)[1], viridis(100)[50], viridis(100)[100]))
}

## =========================================================================
## per-SNP
## =========================================================================
Msnp <- t(GTs_all)                                              # markers x individuals
chr_snp <- as.integer(sub("Chr", "", sub(":.*", "", rownames(Msnp))))
pos_snp <- as.integer(sub(".*:", "", rownames(Msnp)))
ord_m <- order(chr_snp, pos_snp)
a_snp <- orient_aqu(Msnp[ord_m, ], faqu, fpol)
render_pop(a_snp, chr_snp[ord_m], "SNP", "diagnostic SNPs")

## =========================================================================
## per-eMLG (5 cM): consensus for >2-marker clusters, representative otherwise
## =========================================================================
g <- readRDS("module_di25/data/di25_clustering_cM5.rds")$groups
is_emlg <- g$n_loci > 2
D <- vapply(seq_len(nrow(g)), function(i) {
  if (is_emlg[i]) consensus_dosage(GTs_all, g$members[[i]]) else GTs_all[, g$representative[i]]
}, numeric(nrow(GTs_all)))
D <- t(D)                                                       # units x individuals
rep_chr <- as.integer(sub("Chr", "", sub(":.*", "", g$representative)))
rep_pos <- as.integer(sub(".*:", "", g$representative))
ord_u <- order(rep_chr, rep_pos)
a_emlg <- orient_aqu(D[ord_u, ], faqu, fpol)
render_pop(a_emlg, rep_chr[ord_u], "eMLG", "LD-reduced units")

message("[pop-circos] done")


## ============================================================================
## FILE: module_di25/R/di25_architecture.R
## ============================================================================

## =========================================================
## module_di25 (high-DI analyses) -- is sorting stronger in low-recombination regions?
## =========================================================
## Module A's architecture analysis (Table 3 / Fig S4b / Table 4), restricted to
## the high-DI (DI > -25) loci. Full-data finding to compare against: at the UNIT
## level sorting was LOWEST in low recombination (increased with recombination),
## while the per-SNP level looked inflated there from LD pseudoreplication.
##
## Recombination (cM/Mb) is interpolated per chromosome from the recmap exactly as
## Module A does (approx(..., rule = 2)). Sorting is taken at the primary tau = 0.6
## (phi = 0.85): sort_class is re-derived from the tau-independent counts.
##
## Run from the repo root:  Rscript module_di25/R/di25_architecture.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()

RECMAP <- "data/Frufa_DTOL_PR.ref_genome.recmap"
TAU    <- 0.6
OUTDIR <- "module_di25/data"; FIGDIR <- "module_di25/Figures"

## ---- recombination interpolation (cM/Mb), Module A convention ------------
rec <- fread(RECMAP); setnames(rec, 1:4, c("chr", "pos", "cM", "cMMb"))
rec[, Chr := sub("chromosome_", "Chr", chr)]
assign_recomb <- function(chr_int, pos) {
  out <- rep(NA_real_, length(pos))
  for (ch in unique(chr_int)) {
    r <- rec[Chr == paste0("Chr", ch)]; if (nrow(r) < 2) next
    idx <- which(chr_int == ch)
    out[idx] <- approx(r$pos, r$cMMb, xout = pos[idx], rule = 2)$y
  }
  out
}

## ---- load sorting counts + positions, add recomb & tau=0.6 sort class ----
add_sort <- function(dt) {
  ok <- dt$differentiated & dt$n_obs > 0 & is.finite(dt$uni_score)
  cl <- rep(NA_character_, nrow(dt))
  cl[ok] <- classify_sort(dt$n_aqu[ok], dt$n_pol[ok], dt$n_obs[ok],
                          sort_th = TAU, sort_rule = "binom", alpha = 0.05)
  dt[, `:=`(cls = cl,
            sorted = cl %in% c("aquilonia", "polyctena", "unresolved", "ambiguous"),
            is_aqu = cl == "aquilonia", directional = cl %in% c("aquilonia", "polyctena"))]
  dt
}

snp <- as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_snp.rds")))
snp[, `:=`(chr = as.integer(sub("Chr", "", sub(":.*", "", marker))),
           pos = as.integer(sub(".*:", "", marker)))]
snp[, recomb := assign_recomb(chr, pos)]
add_sort(snp)

emlg <- as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_emlg.rds")))
g <- readRDS(file.path(OUTDIR, "di25_clustering_cM5.rds"))$groups
emlg[, rep_mk := g$representative[match(group_id, g$group_id)]]
emlg[, `:=`(chr = as.integer(sub("Chr", "", sub(":.*", "", rep_mk))),
            pos = as.integer(sub(".*:", "", rep_mk)))]
emlg[, recomb := assign_recomb(chr, pos)]
add_sort(emlg)

## ---- fraction sorted by recombination decile (breaks from the SNP set) ---
base_s <- snp[differentiated == TRUE & n_obs > 0 & is.finite(recomb)]
base_e <- emlg[differentiated == TRUE & n_obs > 0 & is.finite(recomb)]
brk <- quantile(base_s$recomb, 0:10 / 10, na.rm = TRUE)
wilson <- function(k, n) {                          # Wilson 95% CI for a proportion
  z <- 1.959964; p <- k / n; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  cbind(lo = pmax(0, ctr - hw), hi = pmin(1, ctr + hw))
}
decile_tab <- function(x, lab) {
  x[, rbin := cut(recomb, brk, include.lowest = TRUE, labels = FALSE)]
  t <- x[!is.na(rbin), .(level = lab, med_recomb = round(median(recomb), 2), n = .N,
                         n_sorted = sum(sorted), frac_sorted = mean(sorted),
                         frac_aqu = sum(is_aqu) / sum(directional)), by = rbin][order(rbin)]
  ci <- wilson(t$n_sorted, t$n); t[, `:=`(lo = ci[, 1], hi = ci[, 2])][]
}
dec <- rbind(decile_tab(base_s, "SNP"), decile_tab(base_e, "eMLG"))
saveRDS(dec, file.path(OUTDIR, "di25_architecture_deciles.rds"))

## ---- STATISTICAL TEST: does P(sorted) depend on recombination? ----------
## eMLG units are ~independent (that is what the LD reduction buys), so one
## logistic observation per unit is a valid, non-pseudoreplicated test.
eu <- base_e[is.finite(recomb)][, `:=`(y = as.integer(sorted), lr = log10(recomb + 0.1))]
N_MIN   <- 100L                                     # "reasonable n" per decile
good_bins <- dec[level == "eMLG" & n >= N_MIN, rbin]
r_lo <- min(base_e[rbin %in% good_bins, recomb])    # recomb floor of the reliable range
glm_all  <- glm(y ~ lr,             data = eu,                family = binomial())
glm_allDI<- glm(y ~ lr + scale(DI),  data = eu[is.finite(DI)], family = binomial())
glm_rr   <- glm(y ~ lr,             data = eu[recomb >= r_lo], family = binomial())
mcf <- function(m) { m0 <- update(m, . ~ 1); round(1 - as.numeric(logLik(m) / logLik(m0)), 4) }
sm  <- function(m) { c <- summary(m)$coef["lr", ]; sprintf("%+.3f (p=%.1e)", c[1], c[4]) }
## effect size in probability terms: predicted P(sorted) across the reliable range
pr_at <- function(x) predict(glm_rr, data.frame(lr = log10(x + 0.1)), type = "response")
p_lo <- pr_at(r_lo); p_hi <- pr_at(quantile(eu[recomb >= r_lo]$recomb, 0.95))
test_txt <- sprintf("P(sorted eMLG)~log10(recomb): all %s (R2=%.3f) ; reliable(recomb>=%.1f,n=%d) %s (R2=%.3f) | P(sorted) %.1f->%.1f%% across the reliable range",
                    sm(glm_all), mcf(glm_all), r_lo, nrow(eu[recomb >= r_lo]), sm(glm_rr), mcf(glm_rr),
                    100 * p_lo, 100 * p_hi)

## ---- magnitude & direction regressions (unit level, Module A form) -------
du <- base_e[is.finite(DI)]
du[, `:=`(zr = as.numeric(scale(log10(recomb + 0.1))), zDI = as.numeric(scale(DI)))]
lm_mag <- lm(prop_fixed ~ zr + zDI, data = du)

cu <- base_e[directional == TRUE & is.finite(DI) & is.finite(parent_maf)]
cu[, `:=`(zr = as.numeric(scale(log10(recomb + 0.1))), zDI = as.numeric(scale(DI)),
          zmaf = as.numeric(scale(parent_maf)), zcs = as.numeric(scale(log10(n_loci))))]
glm_dir <- glm(is_aqu ~ zDI + zr + zmaf + zcs, data = cu, family = binomial())

## also the SNP-level magnitude slope, for the pseudoreplication contrast
ds <- base_s[is.finite(DI)]
ds[, `:=`(zr = as.numeric(scale(log10(recomb + 0.1))), zDI = as.numeric(scale(DI)))]
lm_mag_snp <- lm(prop_fixed ~ zr + zDI, data = ds)

sp <- function(a, b) round(cor(a, b, method = "spearman", use = "pairwise.complete.obs"), 3)

## ---- report -------------------------------------------------------------
cat("\n===== DI25: sorting vs recombination (tau = 0.6, phi = 0.85) =====\n")
cat(sprintf("Spearman DI vs recomb (SNP): %s  |  DI vs recomb (eMLG): %s\n",
            sp(base_s$DI, base_s$recomb), sp(base_e$DI, base_e$recomb)))
cat("\nFraction sorted by recombination decile (1 = lowest recomb):\n")
print(dcast(dec, rbin + med_recomb ~ level, value.var = "frac_sorted")[order(rbin)])
cat("\nFraction toward F. aquilonia by recombination decile:\n")
print(dcast(dec, rbin ~ level, value.var = "frac_aqu")[order(rbin)])
cat(sprintf("\nMagnitude  prop_fixed ~ log10(recomb)+DI :\n  eMLG  recomb %+.3f (t=%.1f), DI %+.3f (t=%.1f)\n  SNP   recomb %+.3f (t=%.1f), DI %+.3f (t=%.1f)\n",
            coef(lm_mag)["zr"], summary(lm_mag)$coef["zr","t value"],
            coef(lm_mag)["zDI"], summary(lm_mag)$coef["zDI","t value"],
            coef(lm_mag_snp)["zr"], summary(lm_mag_snp)$coef["zr","t value"],
            coef(lm_mag_snp)["zDI"], summary(lm_mag_snp)$coef["zDI","t value"]))
cat(sprintf("\nDirection  P(aquilonia) ~ DI+recomb+MAF+logNloci (eMLG, directional):\n  DI %+.3f (z=%.1f) ; recomb %+.3f (z=%.1f)\n",
            coef(glm_dir)["zDI"], summary(glm_dir)$coef["zDI","z value"],
            coef(glm_dir)["zr"], summary(glm_dir)$coef["zr","z value"]))
cat("\nSTATISTICAL TEST -- ", test_txt, "\n", sep = "")
cat("n independent eMLG units per recomb decile: ",
    paste(dec[level == "eMLG", n], collapse = ", "), "\n")

saveRDS(list(deciles = dec, lm_mag = coef(summary(lm_mag)),
             lm_mag_snp = coef(summary(lm_mag_snp)), glm_dir = coef(summary(glm_dir)),
             glm_sorted_all = coef(summary(glm_all)), glm_sorted_allDI = coef(summary(glm_allDI)),
             glm_sorted_reliable = coef(summary(glm_rr)), n_min = N_MIN, recomb_floor = r_lo,
             di_recomb_rho = c(SNP = sp(base_s$DI, base_s$recomb), eMLG = sp(base_e$DI, base_e$recomb))),
        file.path(OUTDIR, "di25_architecture.rds"))

## ---- figure: eMLG fraction sorted vs recombination decile, with n bars ---
## x = recombination decile (median cM/Mb labelled). Background grey bars = the
## number of INDEPENDENT eMLG units in each decile (scaled to the right axis) --
## they collapse to a handful in low recombination. eMLG points carry Wilson 95%
## CIs and are sized by n; SNP is a faint reference. Deciles below the reliable-n
## floor are shaded.
e <- dec[level == "eMLG"]; s <- dec[level == "SNP"]
ymax <- max(e$hi, s$frac_sorted, na.rm = TRUE) * 1.02
bs   <- ymax / max(e$n)                              # bar scale to the fraction axis
lo_bins <- e[n < N_MIN, rbin]
pA <- ggplot(e, aes(rbin)) +
  { if (length(lo_bins)) annotate("rect", xmin = min(lo_bins) - 0.5, xmax = max(lo_bins) + 0.5,
                                  ymin = -Inf, ymax = Inf, fill = "grey95") } +
  geom_col(aes(y = n * bs), fill = "grey85", width = 0.85) +
  geom_line(data = s, aes(y = frac_sorted), colour = "grey65", linewidth = 0.6) +
  geom_point(data = s, aes(y = frac_sorted), colour = "grey65", size = 1.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, colour = "#315B7D") +
  geom_line(aes(y = frac_sorted), colour = "#315B7D", linewidth = 0.8) +
  geom_point(aes(y = frac_sorted, size = n), colour = "#315B7D") +
  scale_size_area(max_size = 5, guide = "none") +
  scale_x_continuous(breaks = 1:10, labels = round(e$med_recomb, 1)) +
  scale_y_continuous(name = "fraction sorted (eMLG; SNP faint)",
                     sec.axis = sec_axis(~ . / bs, name = "n independent eMLG units (bars)")) +
  labs(x = "recombination decile (median cM/Mb)",
       title = "Sorting is weakly stronger at lower recombination (high-DI loci)",
       subtitle = sprintf("real & highly significant but SMALL: reliable range %s, McFadden R2=%.3f; P(sorted) %.0f%%->%.0f%% across the range",
                          sm(glm_rr), mcf(glm_rr), 100 * p_lo, 100 * p_hi)) +
  theme_bw(base_size = 10)
ggsave(file.path(FIGDIR, "di25_sorting_vs_recomb.png"), pA, width = 9, height = 5.2, dpi = 200)
cat("\n[di25-architecture] wrote", file.path(FIGDIR, "di25_sorting_vs_recomb.png"), "\n")


## ============================================================================
## FILE: module_di25/R/di25_recomb_tau_sweep.R
## ============================================================================

## =========================================================
## module_di25 (high-DI analyses) -- sorting vs recombination, swept over tau
## =========================================================
## Is high-DI sorting stronger at low recombination, and does that hold across
## the sorting threshold tau? phi = 0.85 fixed; tau in {0.5,0.6,0.7,0.8}.
##
## Inference note: LD-reduced units are NOT independent -- they are much more
## independent than SNPs (residual LD is reduced, not removed). So the per-unit
## logistic model SE is anticonservative. We therefore lead with EFFECT SIZE
## (McFadden R2; predicted P(sorted) across the range) and take CIs on the
## recombination coefficient from a GENOMIC BLOCK BOOTSTRAP that resamples whole
## chromosomes (preserving within-chromosome residual dependence). The per-SNP
## level is shown only as the pseudoreplicated foil.
##
## Run from the repo root:  Rscript module_di25/R/di25_recomb_tau_sweep.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()

TAU_GRID <- c(0.5, 0.6, 0.7, 0.8)
N_BOOT   <- 400L
N_MIN    <- 100L
SEED     <- 1
OUTDIR   <- "module_di25/data"; FIGDIR <- "module_di25/Figures"
set.seed(SEED)

## ---- recombination ------------------------------------------------------
rec <- fread("data/Frufa_DTOL_PR.ref_genome.recmap"); setnames(rec, 1:4, c("chr","p","cM","cMMb"))
rec[, Chr := sub("chromosome_", "Chr", chr)]
rc <- function(ch, pos) { o <- rep(NA_real_, length(pos))
  for (cc in unique(ch)) { r <- rec[Chr == paste0("Chr", cc)]; if (nrow(r) < 2) next
    i <- which(ch == cc); o[i] <- approx(r$p, r$cMMb, xout = pos[i], rule = 2)$y }; o }

## ---- load counts + positions (tau-independent) --------------------------
e <- as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_emlg.rds")))
g <- readRDS(file.path(OUTDIR, "di25_clustering_cM5.rds"))$groups
e[, rmk := g$representative[match(group_id, g$group_id)]]
e[, `:=`(chr = as.integer(sub("Chr","",sub(":.*","",rmk))), pos = as.integer(sub(".*:","",rmk)))]
e[, recomb := rc(chr, pos)]; e <- e[differentiated == TRUE & n_obs > 0 & is.finite(uni_score) & is.finite(recomb)]
e[, lr := log10(recomb + 0.1)]

s <- as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_snp.rds")))
s[, `:=`(chr = as.integer(sub("Chr","",sub(":.*","",marker))), pos = as.integer(sub(".*:","",marker)))]
s[, recomb := rc(chr, pos)]; s <- s[differentiated == TRUE & n_obs > 0 & is.finite(uni_score) & is.finite(recomb)]

## decile breaks from the SNP recomb distribution (shared across levels & tau)
brk <- quantile(s$recomb, 0:10/10, na.rm = TRUE)
e[, rbin := cut(recomb, brk, include.lowest = TRUE, labels = FALSE)]
s[, rbin := cut(recomb, brk, include.lowest = TRUE, labels = FALSE)]
good <- e[!is.na(rbin), .N, by = rbin][N >= N_MIN, rbin]; r_lo <- min(e[rbin %in% good, recomb])
chr_idx <- split(seq_len(nrow(e[recomb >= r_lo])), e[recomb >= r_lo]$chr)   # for block bootstrap

sorted_at <- function(dt, tau) {
  cl <- classify_sort(dt$n_aqu, dt$n_pol, dt$n_obs, sort_th = tau, sort_rule = "binom", alpha = 0.05)
  as.integer(cl %in% c("aquilonia","polyctena","unresolved","ambiguous"))
}
mcf <- function(m) { m0 <- update(m, . ~ 1); as.numeric(1 - logLik(m)/logLik(m0)) }

## ---- sweep --------------------------------------------------------------
dec_all <- list(); summ <- list()
for (tau in TAU_GRID) {
  e[, y := sorted_at(e, tau)]; s[, y := sorted_at(s, tau)]
  de <- e[!is.na(rbin), .(level="eMLG", tau=tau, n=.N, med_recomb=round(median(recomb),1),
                          frac_sorted=mean(y)), by=rbin][order(rbin)]
  ds <- s[!is.na(rbin), .(level="SNP",  tau=tau, n=.N, med_recomb=round(median(recomb),1),
                          frac_sorted=mean(y)), by=rbin][order(rbin)]
  dec_all[[length(dec_all)+1]] <- rbind(de, ds)

  er <- e[recomb >= r_lo]
  m  <- glm(y ~ lr, data = er, family = binomial())
  ## genomic block bootstrap: resample whole chromosomes (residual LD is within-chr)
  chs <- names(chr_idx)
  bc <- vapply(seq_len(N_BOOT), function(b) {
    idx <- unlist(chr_idx[sample(chs, length(chs), replace = TRUE)], use.names = FALSE)
    coef(glm.fit(cbind(1, er$lr[idx]), er$y[idx], family = binomial()))[2]
  }, numeric(1))
  pr <- predict(m, data.frame(lr = log10(c(r_lo, quantile(er$recomb, 0.95)) + 0.1)), type = "response")
  summ[[length(summ)+1]] <- data.table(
    tau = tau, frac_sorted_overall = round(mean(e$y), 3),
    coef_recomb = round(coef(m)["lr"], 3),
    boot_lo = round(quantile(bc, 0.025), 3), boot_hi = round(quantile(bc, 0.975), 3),
    McFadden_R2 = round(mcf(m), 4),
    P_sorted_loRecomb = round(100*pr[1], 1), P_sorted_hiRecomb = round(100*pr[2], 1))
}
dec <- rbindlist(dec_all); sweep <- rbindlist(summ)
saveRDS(list(deciles = dec, sweep = sweep, r_lo = r_lo,
             n_units_per_decile = e[!is.na(rbin), .N, by=rbin][order(rbin)]),
        file.path(OUTDIR, "di25_recomb_tau_sweep.rds"))

cat("\n===== sorting vs recombination, tau sweep (eMLG, reliable range recomb>=", round(r_lo,1),
    "cM/Mb; block-bootstrap CI over chromosomes) =====\n", sep="")
print(sweep)
cat("\nn independent eMLG units per recomb decile (tau-independent): ",
    paste(e[!is.na(rbin), .N, by=rbin][order(rbin)]$N, collapse=", "), "\n")

## ---- figure -------------------------------------------------------------
ne <- e[!is.na(rbin), .N, by=rbin][order(rbin)]
med_lab <- e[!is.na(rbin), round(median(recomb), 1), by=rbin][order(rbin)]$V1
de <- dec[level == "eMLG"]; de[, tau := factor(tau)]
ymax <- max(de$frac_sorted); bs <- ymax*0.98 / max(ne$N)
lo_bins <- ne[N < N_MIN, rbin]
pA <- ggplot(de, aes(rbin)) +
  { if (length(lo_bins)) annotate("rect", xmin=min(lo_bins)-0.5, xmax=max(lo_bins)+0.5, ymin=-Inf, ymax=Inf, fill="grey95") } +
  geom_col(data = ne, aes(rbin, N*bs), fill = "grey86", width = 0.85, inherit.aes = FALSE) +
  geom_line(aes(y = frac_sorted, colour = tau), linewidth = 0.8) +
  geom_point(aes(y = frac_sorted, colour = tau), size = 1.6) +
  scale_colour_viridis_d(end = 0.9, name = expression(tau)) +
  scale_x_continuous(breaks = 1:10, labels = med_lab) +
  scale_y_continuous(name = "fraction sorted (eMLG units)",
                     sec.axis = sec_axis(~ ./bs, name = "n independent units (bars)")) +
  labs(x = "recombination decile (median cM/Mb)",
       title = "Sorting vs recombination across tau (high-DI eMLG units)",
       subtitle = "sorting is weakly higher at low recombination at every tau; grey = low-n deciles (not reliably estimable)") +
  theme_bw(base_size = 10)
pB <- ggplot(sweep, aes(factor(tau), coef_recomb)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_errorbar(aes(ymin = boot_lo, ymax = boot_hi), width = 0.18, colour = "#315B7D") +
  geom_point(size = 2.4, colour = "#315B7D") +
  labs(x = expression(tau), y = "recomb coef (logit / log10 cM/Mb)",
       title = "b  Effect + block-bootstrap 95% CI",
       subtitle = "negative = sorting stronger at low recomb") +
  theme_bw(base_size = 10)
ggsave(file.path(FIGDIR, "di25_recomb_tau_sweep.png"), pA + pB + plot_layout(widths = c(2, 1)),
       width = 13, height = 5.2, dpi = 200)
cat("\n[di25-recomb-sweep] wrote", file.path(FIGDIR, "di25_recomb_tau_sweep.png"), "\n")
