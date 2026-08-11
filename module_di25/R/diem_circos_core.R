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

DIEM_COL    <- c("aqu" = "#21918C",   # 1 hom F. aquilonia  (F. aquilonia teal)
                 "het" = "#440154",   # 2 heterozygous      (heterozygous dark viridis)
                 "pol" = "#D3C93B")   # 3 hom F. polyctena  (F. polyctena darker yellow)
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
                                 chr_label_cex = 0.5, legend_cex = 0.75, main_line = 0.2,
                                 legend_labels = NULL, legend_cols = NULL,
                                 bg_col = "#FFFFFF", new_device = TRUE, res = 300,
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
    if (new_device) png(outpng, width = npx, height = npx, res = res, type = "cairo")
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
    text(x, y, paste0(chr_label_prefix, chr_levels[i]), cex = chr_label_cex, col = "grey20",
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
           cex = legend_cex, border = NA, inset = c(0.02, 0.02))
  if (nzchar(title)) title(main = title, cex.main = cex_main, line = main_line)
  invisible(list(n_units = nrow(code), n_rings = nrings))
}

render_diem_circos <- function(gt, chr_num, outpng = NULL, title = "",
                               col = DIEM_COL, na_col = DIEM_NA_COL,
                               bg_col = DIEM_BG_COL, dither = TRUE,
                               npx = 2600, r_in = 0.30, r_out = 0.92,
                               start_deg = 90, clockwise = TRUE,
                               open_deg = 12, chr_gap_deg = 0.5,
                               seed = 1, chr_label_prefix = "Chr ", chr_label_cex = 0.5,
                               cex_main = 0.85, new_device = TRUE, res = 300,
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
  if (new_device) png(outpng, width = npx, height = npx, res = res, type = "cairo")
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
    text(x, y, paste0(chr_label_prefix, chr_levels[i]), cex = chr_label_cex, col = "grey20",
         srt = if (flip) deg + 180 else deg, adj = if (flip) c(1, 0.5) else c(0, 0.5))
  }
  if (nzchar(title)) title(main = title, cex.main = cex_main, line = 0.2)
  invisible(list(n_miss = n_miss, n_units = nrow(gt), n_ind = ncol(gt)))
}
