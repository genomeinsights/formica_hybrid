plot_ld_outlier_circos <- function(
    el,
    map,
    map_manh,
    r2_th = 0.5,
    only_inter = TRUE,
    core_snps = NULL,
    pad_bp = 5e4,
    highlight_linked_or = TRUE,
    show_legend = TRUE,
    show_chr_labels = TRUE,
    show_chr_ticks = TRUE,
    max_links = NULL,
    bg_col = "white",
    outlier_col = "grey60",
    linked_outlier_col = "grey20",
    chr_line_col = "black",
    palette_name = "Zissou1",
    link_alpha_range = c(0.25, 0.9),
    link_lwd_range = c(0.6, 2.5),
    start_degree = 90,
    gap_degree = 2,
    legend_x_offset = 0.12,
    legend_y_bottom = 0.35,
    legend_y_top = 0.75,
    add_manhattan = FALSE,
    manhattan_stat_col = NULL,
    manhattan_transform = NULL,
    manhattan_type = c("points", "segments"),
    manhattan_track_height = 0.10,
    manhattan_cex = 0.25,
    manhattan_col = "#444444",
    manhattan_outlier_col = "#D55E00",
    manhattan_track_margin = c(0.01, 0.01),
    manhattan_threshold = NULL,
    manhattan_show_labels = FALSE,
    manhattan_label_cex = 0.45,
    DIEM = NULL,
    add_diem = FALSE,
    diem_summary = c("mean", "sd", "both"),
    diem_track_height = 0.08,
    diem_track_margin = c(0.01, 0.01),
    diem_col = "#2C7FB8",
    diem_sd_col = "#7FCDBB",
    diem_outlier_col = "#D95F0E",
    diem_type = c("points", "segments", "line"),
    diem_cex = 0.2,
    diem_lwd = 1,
    diem_ylim = NULL,
    diem_show_labels = TRUE,
    diem_label_cex = 0.45
) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stopifnot(requireNamespace("circlize", quietly = TRUE))
  stopifnot(requireNamespace("wesanderson", quietly = TRUE))

  manhattan_type <- match.arg(manhattan_type)
  diem_summary <- match.arg(diem_summary)
  diem_type <- match.arg(diem_type)

  old_par <- par(no.readonly = TRUE)
  on.exit({
    try(circos.clear(), silent = TRUE)
    try(par(old_par), silent = TRUE)
  }, add = TRUE)

  el <- as.data.table(copy(el))
  map <- as.data.table(copy(map))
  map_manh <- as.data.table(copy(map_manh))

  if (!all(c("Chr", "Pos", "marker") %in% names(map))) {
    stop("map must contain columns: Chr, Pos, marker")
  }
  if (!all(c("marker", "Chr", "Pos", "OR_id") %in% names(map_manh))) {
    stop("map_manh must contain columns: marker, Chr, Pos, OR_id")
  }
  if (!all(c("SNP1", "SNP2", "Chr1", "Chr2", "pos1", "pos2", "r2") %in% names(el))) {
    stop("el must contain columns: SNP1, SNP2, Chr1, Chr2, pos1, pos2, r2")
  }

  genome <- map[, .(chr = as.character(Chr), pos = as.numeric(Pos), snp = as.character(marker))]
  chr_df <- genome[, .(chr_len = max(pos, na.rm = TRUE)), by = chr][order(chr)]

  if (nrow(chr_df) == 0) stop("No chromosome information found in map.")

  outlier_regions <- map_manh[OR_id != "ns" & !is.na(OR_id),
                              .(
                                start = min(Pos, na.rm = TRUE),
                                end   = max(Pos, na.rm = TRUE)
                              ),
                              by = .(chr = as.character(Chr), OR_id = as.character(OR_id))
  ][order(chr, start)]

  if (nrow(outlier_regions) > 0 && !is.null(pad_bp) && pad_bp > 0) {
    outlier_regions[, `:=`(
      start = pmax(0, start - pad_bp),
      end   = end + pad_bp
    )]
  }

  # ---------- Manhattan prep ----------
  manh_list <- list()

  if (isTRUE(add_manhattan)) {
    if (is.null(manhattan_stat_col)) {
      if ("p_BF" %in% names(map_manh)) {
        manhattan_stat_col <- "p_BF"
      } else if ("F_BF" %in% names(map_manh)) {
        manhattan_stat_col <- "F_BF"
      } else {
        stop("add_manhattan=TRUE but no manhattan_stat_col supplied, and neither p_BF nor F_BF found in map_manh.")
      }
    }

    manhattan_stat_col <- as.character(manhattan_stat_col)

    missing_stats <- setdiff(manhattan_stat_col, names(map_manh))
    if (length(missing_stats) > 0) {
      stop("These manhattan_stat_col values are not in map_manh: ", paste(missing_stats, collapse = ", "))
    }

    # normalize transforms to a named list
    if (is.null(manhattan_transform)) {
      manhattan_transform <- lapply(manhattan_stat_col, function(stat) {
        if (identical(stat, "p_BF")) {
          function(x) -log10(pmax(x, 1e-300))
        } else {
          identity
        }
      })
      names(manhattan_transform) <- manhattan_stat_col
    } else if (is.function(manhattan_transform)) {
      tmp <- replicate(length(manhattan_stat_col), manhattan_transform, simplify = FALSE)
      names(tmp) <- manhattan_stat_col
      manhattan_transform <- tmp
    } else if (is.list(manhattan_transform)) {
      if (is.null(names(manhattan_transform))) {
        if (length(manhattan_transform) != length(manhattan_stat_col)) {
          stop("If manhattan_transform is an unnamed list, it must have the same length as manhattan_stat_col.")
        }
        names(manhattan_transform) <- manhattan_stat_col
      }
    } else {
      stop("manhattan_transform must be NULL, a function, or a list of functions.")
    }

    # normalize thresholds
    if (is.null(manhattan_threshold)) {
      threshold_vec <- setNames(rep(NA_real_, length(manhattan_stat_col)), manhattan_stat_col)
    } else if (length(manhattan_threshold) == 1) {
      threshold_vec <- setNames(rep(as.numeric(manhattan_threshold), length(manhattan_stat_col)), manhattan_stat_col)
    } else {
      threshold_vec <- as.numeric(manhattan_threshold)
      if (length(threshold_vec) != length(manhattan_stat_col)) {
        stop("manhattan_threshold must have length 1 or the same length as manhattan_stat_col.")
      }
      names(threshold_vec) <- manhattan_stat_col
    }

    # normalize track heights
    if (length(manhattan_track_height) == 1) {
      track_height_vec <- setNames(rep(as.numeric(manhattan_track_height), length(manhattan_stat_col)), manhattan_stat_col)
    } else {
      track_height_vec <- as.numeric(manhattan_track_height)
      if (length(track_height_vec) != length(manhattan_stat_col)) {
        stop("manhattan_track_height must have length 1 or the same length as manhattan_stat_col.")
      }
      names(track_height_vec) <- manhattan_stat_col
    }

    # build one manhattan table per stat
    for (stat in manhattan_stat_col) {
      tf <- manhattan_transform[[stat]]
      if (!is.function(tf)) {
        stop(sprintf("manhattan_transform for '%s' is not a function.", stat))
      }

      dd <- map_manh[, .(
        chr = as.character(Chr),
        pos = as.numeric(Pos),
        snp = as.character(marker),
        OR_id = as.character(OR_id),
        raw_stat = get(stat)
      )]

      dd[, stat := tf(raw_stat)]
      dd <- dd[is.finite(stat) & !is.na(chr) & !is.na(pos)]

      thr <- threshold_vec[[stat]]
      # keep all points; threshold is only for line, not filtering

      if (nrow(dd) == 0) next

      ylim_here <- c(0, max(dd$stat, na.rm = TRUE))
      if (!is.finite(ylim_here[2]) || ylim_here[2] <= 0) ylim_here <- c(0, 1)

      manh_list[[stat]] <- list(
        data = dd,
        ylim = ylim_here,
        threshold = thr,
        height = track_height_vec[[stat]],
        label = stat
      )
    }

    if (length(manh_list) == 0) {
      warning("No Manhattan tracks could be built; Manhattan plotting will be skipped.")
      add_manhattan <- FALSE
    }
  }

  # ---------- DIEM prep ----------
  diem_tracks <- list()

  if (isTRUE(add_diem)) {
    if (is.null(DIEM)) {
      stop("add_diem=TRUE but DIEM is NULL.")
    }

    if (is.null(colnames(DIEM))) {
      stop("DIEM must have locus names in colnames(DIEM) matching map$marker.")
    }

    diem_markers <- intersect(colnames(DIEM), map$marker)
    if (length(diem_markers) == 0) {
      stop("No overlap between colnames(DIEM) and map$marker.")
    }

    diem_map <- merge(
      data.table(marker = diem_markers),
      map[, .(marker, Chr, Pos)],
      by = "marker",
      all.x = TRUE
    )

    diem_map <- diem_map[!is.na(Chr) & !is.na(Pos)]
    setorder(diem_map, Chr, Pos)

    DIEM_sub <- DIEM[, diem_map$marker, drop = FALSE]

    if (diem_summary %in% c("mean", "both")) {
      dd_mean <- data.table(
        chr = as.character(diem_map$Chr),
        pos = as.numeric(diem_map$Pos),
        marker = diem_map$marker,
        OR_id = map_manh[match(diem_map$marker, marker), OR_id],
        stat = colMeans(DIEM_sub, na.rm = TRUE)
      )

      if (is.null(diem_ylim)) {
        ylim_mean <- c(0, 1)
      } else {
        ylim_mean <- diem_ylim
      }

      diem_tracks[["DIEM mean"]] <- list(
        data = dd_mean,
        ylim = ylim_mean,
        col = diem_col,
        outlier_col = diem_outlier_col,
        height = diem_track_height,
        label = "DIEM mean"
      )
    }

    if (diem_summary %in% c("sd", "both")) {
      dd_sd <- data.table(
        chr = as.character(diem_map$Chr),
        pos = as.numeric(diem_map$Pos),
        marker = diem_map$marker,
        OR_id = map_manh[match(diem_map$marker, marker), OR_id],
        stat = apply(DIEM_sub, 2, sd, na.rm = TRUE)
      )

      ylim_sd <- c(0, max(dd_sd$stat, na.rm = TRUE))
      if (!is.finite(ylim_sd[2]) || ylim_sd[2] <= 0) ylim_sd <- c(0, 1)

      diem_tracks[["DIEM sd"]] <- list(
        data = dd_sd,
        ylim = ylim_sd,
        col = diem_sd_col,
        outlier_col = diem_outlier_col,
        height = diem_track_height,
        label = "DIEM sd"
      )
    }
  }

  # Additional DIEM tracks
  if (isTRUE(add_diem) && length(diem_tracks) > 0) {
    for (nm in names(diem_tracks)) {
      tr <- diem_tracks[[nm]]

      circos.trackPlotRegion(
        ylim = tr$ylim,
        track.height = tr$height,
        track.margin = diem_track_margin,
        bg.border = NA,
        panel.fun = function(x, y) {
          chr_now <- CELL_META$sector.index
          dd <- tr$data[chr == chr_now]

          if (nrow(dd) > 0) {
            cols <- ifelse(dd$OR_id != "ns" & !is.na(dd$OR_id), tr$outlier_col, tr$col)

            if (diem_type == "points") {
              circos.points(
                x = dd$pos,
                y = dd$stat,
                pch = 16,
                cex = diem_cex,
                col = cols
              )
            } else if (diem_type == "segments") {
              circos.segments(
                x0 = dd$pos,
                y0 = tr$ylim[1],
                x1 = dd$pos,
                y1 = dd$stat,
                col = cols,
                lwd = 0.4
              )
            } else if (diem_type == "line") {
              circos.lines(
                x = dd$pos,
                y = dd$stat,
                col = tr$col,
                lwd = diem_lwd
              )

              # optionally overlay outlier points
              dd_or <- dd[OR_id != "ns" & !is.na(OR_id)]
              if (nrow(dd_or) > 0) {
                circos.points(
                  x = dd_or$pos,
                  y = dd_or$stat,
                  pch = 16,
                  cex = diem_cex,
                  col = tr$outlier_col
                )
              }
            }
          }

          if (isTRUE(diem_show_labels)) {
            circos.text(
              CELL_META$xcenter,
              tr$ylim[2] * 0.98,
              labels = tr$label,
              facing = "inside",
              niceFacing = TRUE,
              cex = diem_label_cex
            )
          }
        }
      )
    }
  }
  # ---------- edges ----------
  edges <- el[!is.na(r2) & r2 >= r2_th]

  if (!is.null(core_snps)) {
    core_snps <- as.character(core_snps)
    edges <- edges[SNP1 %in% core_snps & SNP2 %in% core_snps]
  }

  if (only_inter) {
    edges <- edges[Chr1 != Chr2]
  }

  edges <- edges[SNP1 != SNP2]

  if (nrow(edges) == 0) {
    warning("No edges passed the filters.")
  }

  if (nrow(edges) > 0) {
    edges[, pair_id := ifelse(
      paste(Chr1, pos1, SNP1) < paste(Chr2, pos2, SNP2),
      paste(Chr1, pos1, SNP1, Chr2, pos2, SNP2, sep = "|"),
      paste(Chr2, pos2, SNP2, Chr1, pos1, SNP1, sep = "|")
    )]
    edges <- edges[, .SD[which.max(r2)], by = pair_id]
  }

  if (!is.null(max_links) && nrow(edges) > max_links) {
    setorder(edges, -r2)
    edges <- edges[seq_len(max_links)]
  }

  or_map <- unique(map_manh[OR_id != "ns", .(marker = as.character(marker), OR_id = as.character(OR_id))])

  if (nrow(edges) > 0) {
    edges <- merge(edges, or_map, by.x = "SNP1", by.y = "marker", all.x = TRUE)
    setnames(edges, "OR_id", "OR1")
    edges <- merge(edges, or_map, by.x = "SNP2", by.y = "marker", all.x = TRUE)
    setnames(edges, "OR_id", "OR2")
    edges[is.na(OR1), OR1 := "ns"]
    edges[is.na(OR2), OR2 := "ns"]
  } else {
    edges[, `:=`(OR1 = character(), OR2 = character())]
  }

  if (highlight_linked_or && nrow(outlier_regions) > 0 && nrow(edges) > 0) {
    linked_or <- unique(c(edges[OR1 != "ns", OR1], edges[OR2 != "ns", OR2]))
    outlier_regions[, linked := OR_id %in% linked_or]
  } else {
    outlier_regions[, linked := FALSE]
  }

  # ---------- link colors ----------
  if (nrow(edges) > 0) {
    z5 <- wes_palette(palette_name, 5, type = "discrete")
    r2_min <- min(edges$r2, na.rm = TRUE)
    r2_mid <- stats::median(edges$r2, na.rm = TRUE)
    r2_max <- max(edges$r2, na.rm = TRUE)

    if (isTRUE(all.equal(r2_min, r2_max))) {
      col_fun <- function(x) rep(z5[5], length(x))
      s <- rep(1, nrow(edges))
    } else {
      col_fun <- circlize::colorRamp2(
        c(r2_min, r2_mid, r2_max),
        c(z5[1], z5[3], z5[5])
      )
      s <- (edges$r2 - r2_min) / (r2_max - r2_min)
    }

    link_alpha <- link_alpha_range[1] + (link_alpha_range[2] - link_alpha_range[1]) * s
    link_lwd   <- link_lwd_range[1] + (link_lwd_range[2] - link_lwd_range[1]) * s
    base_cols <- col_fun(edges$r2)
    link_cols <- mapply(
      function(cl, a) grDevices::adjustcolor(cl, alpha.f = a),
      base_cols,
      link_alpha,
      USE.NAMES = FALSE
    )
  } else {
    col_fun <- NULL
    link_cols <- character(0)
    link_lwd <- numeric(0)
  }

  # ---------- layout ----------
  par(bg = bg_col, mar = c(1, 1, 1, if (show_legend) 5 else 1), xpd = TRUE)

  circos.clear()
  circos.par(
    start.degree = start_degree,
    gap.degree = gap_degree,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.002, 0.002),
    points.overflow.warning = FALSE
  )

  circos.initialize(
    factors = chr_df$chr,
    xlim = cbind(rep(0, nrow(chr_df)), chr_df$chr_len)
  )

  # Track 1: chromosome line + labels
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.055,
    bg.border = NA,
    panel.fun = function(x, y) {
      xlim <- CELL_META$xlim
      chr_now <- CELL_META$sector.index

      circos.lines(xlim, c(0.35, 0.35), lwd = 1.6, col = chr_line_col)

      if (show_chr_ticks) {
        circos.segments(xlim[1], 0.25, xlim[1], 0.45, lwd = 1.2, col = chr_line_col)
        circos.segments(xlim[2], 0.25, xlim[2], 0.45, lwd = 1.2, col = chr_line_col)
      }

      if (show_chr_labels) {
        circos.text(
          CELL_META$xcenter, 0.82, chr_now,
          facing = "bending.inside",
          niceFacing = TRUE,
          cex = 0.7
        )
      }
    }
  )

  # Track 2: shaded outlier regions
  # circos.trackPlotRegion(
  #   ylim = c(0, 1),
  #   track.height = 0.045,
  #   bg.border = NA,
  #   panel.fun = function(x, y) {
  #     chr_now <- CELL_META$sector.index
  #     rr <- outlier_regions[chr == chr_now]
  #     if (nrow(rr) > 0) {
  #       for (j in seq_len(nrow(rr))) {
  #         fill_col <- if (isTRUE(highlight_linked_or) && isTRUE(rr$linked[j])) {
  #           grDevices::adjustcolor(linked_outlier_col, alpha.f = 1)
  #         } else {
  #           grDevices::adjustcolor(outlier_col, alpha.f = 1)
  #         }
  #
  #         circos.rect(
  #           xleft = rr$start[j],
  #           ybottom = 0.15,
  #           xright = rr$end[j],
  #           ytop = 0.85,
  #           col = fill_col,
  #           border = NA
  #         )
  #       }
  #     }
  #   }
  # )

  # Additional Manhattan tracks
  if (isTRUE(add_manhattan) && length(manh_list) > 0) {
    for (stat in names(manh_list)) {
      tr <- manh_list[[stat]]

      circos.trackPlotRegion(
        ylim = tr$ylim,
        track.height = tr$height,
        track.margin = manhattan_track_margin,
        bg.border = NA,
        panel.fun = function(x, y) {
          chr_now <- CELL_META$sector.index
          dd <- tr$data[chr == chr_now]

          if (nrow(dd) > 0) {
            cols <- ifelse(dd$OR_id != "ns" & !is.na(dd$OR_id), manhattan_outlier_col, manhattan_col)

            if (manhattan_type == "points") {
              circos.points(
                x = dd$pos,
                y = dd$stat,
                pch = 16,
                cex = manhattan_cex,
                col = cols
              )
            } else {
              circos.segments(
                x0 = dd$pos,
                y0 = 0,
                x1 = dd$pos,
                y1 = dd$stat,
                col = cols,
                cex = manhattan_cex
              )
            }
          }

          if (!is.na(tr$threshold)) {
            circos.lines(
              CELL_META$xlim,
              c(tr$threshold, tr$threshold),
              col = "grey30",
              lwd = 0.7,
              lty = 2
            )
          }

          if (isTRUE(manhattan_show_labels)) {
            circos.text(
              CELL_META$xcenter,
              tr$ylim[2] * 0.98,
              labels = tr$label,
              facing = "inside",
              niceFacing = TRUE,
              cex = manhattan_label_cex
            )
          }
        }
      )
    }
  }

  # Links
  if (nrow(edges) > 0) {
    for (i in seq_len(nrow(edges))) {
      circos.link(
        sector.index1 = edges$Chr1[i],
        point1 = edges$pos1[i],
        sector.index2 = edges$Chr2[i],
        point2 = edges$pos2[i],
        col = link_cols[i],
        lwd = link_lwd[i],
        border = NA
      )
    }
  }

  # Legend
  if (show_legend && nrow(edges) > 0) {
    usr <- par("usr")

    x0 <- usr[2] + legend_x_offset * diff(usr[1:2])
    y0 <- usr[3] + legend_y_bottom * diff(usr[3:4])
    y1 <- usr[3] + legend_y_top * diff(usr[3:4])
    w  <- 0.035 * diff(usr[1:2])

    rvals <- seq(min(edges$r2), max(edges$r2), length.out = 100)
    cols <- vapply(
      col_fun(rvals),
      function(cl) grDevices::adjustcolor(cl, alpha.f = 0.9),
      character(1)
    )
    yy <- seq(y0, y1, length.out = length(cols) + 1)

    for (k in seq_along(cols)) {
      rect(x0, yy[k], x0 + w, yy[k + 1], col = cols[k], border = NA, xpd = TRUE)
    }

    text(x0 + w / 2, y1 + 0.03 * diff(usr[3:4]), labels = expression(r^2), xpd = TRUE)
    ticks <- pretty(range(edges$r2), n = 4)
    ticks <- ticks[ticks >= min(edges$r2) & ticks <= max(edges$r2)]

    for (tv in ticks) {
      ty <- y0 + (tv - min(edges$r2)) / (max(edges$r2) - min(edges$r2) + 1e-12) * (y1 - y0)
      segments(x0 + w, ty, x0 + w + 0.01 * diff(usr[1:2]), ty, xpd = TRUE)
      text(
        x0 + w + 0.02 * diff(usr[1:2]), ty,
        labels = format(round(tv, 2), nsmall = 2),
        adj = c(0, 0.5), cex = 0.8, xpd = TRUE
      )
    }

    # legend(
    #   x = x0,
    #   y = y0 - 0.08 * diff(usr[3:4]),
    #   legend = c("Outlier region", "Linked outlier region"),
    #   pch = 15,
    #   pt.cex = 1.4,
    #   col = c(
    #     grDevices::adjustcolor(outlier_col, alpha.f = 1),
    #     grDevices::adjustcolor(linked_outlier_col, alpha.f = 1)
    #   ),
    #   bty = "n",
    #   xpd = TRUE,
    #   cex = 0.8
    # )
  }

  invisible(list(
    edges = edges,
    chr_df = chr_df,
    outlier_regions = outlier_regions,
    manhattan = manh_list
  ))
}

find_hapl_blocks <- function(ld_decay,map,SNPs,rho_ld=0.95,rho_d=0.95,ld_th=NULL,d_th=NULL,col_vector=NULL){
  #ch = "Chr10"
  if(is.null(col_vector)){
    col_vector <- c("#B2DF8A", "#FFD92F", "firebrick", "#33A02C", "#7FC97F", "#CAB2D6",
                    "#FB8072", "grey30", "#E6AB02", "#FDC086", "steelblue", "#1F78B4",
                    "#FB9A99", "#1B9E77", "#BC80BD", "#E31A1C", "#7570B3", "#A6761D",
                    "#A6CEE3", "salmon", "#FFFF33", "forestgreen", "#FDCDAC", "#BF5B17",
                    "#A6761D", "#FBB4AE", "#4DAF4A", "#B3E2CD", "#FDDAEC", "#BEBADA",
                    "#FFF2AE", "#1F78B4", "#66C2A5", "#F0027F", "#E6AB02", "#E78AC3",
                    "#FF7F00", "#8DA0CB", "#6A3D9A", "#B15928", "#E41A1C")
  }

  map_with_hb <- rbindlist(lapply(names(ld_decay$by_chr), function(ch){
    message(ch)
    mp <- copy(map[Chr==ch])

    chr_obj <- ld_decay$by_chr[[ch]]
    chr_obj$el <- fread(chr_obj$el,showProgress = FALSE)

    if(is.null(ld_th) & is.null(d_th)){
      a_chr <- ld_decay$decay_sum[Chr==ch,a_pred]
      b_chr <- ld_decay$decay_sum[Chr==ch,b]
      c_chr <- ld_decay$decay_sum[Chr==ch,c_pred]

      d_th  <- d_from_rho(a_chr, rho = rho_d)
      ld_th <- ld_from_rho(b_chr, c_chr, rho = rho_ld)
      ed <- chr_obj$el[r2>ld_th & d<d_th & (SNP1 %in% SNPs | SNP2 %in% SNPs),.(SNP1,SNP2)]
    }else{
      ed <- chr_obj$el[r2>ld_th & d<d_th & (SNP1 %in% SNPs | SNP2 %in% SNPs),.(SNP1,SNP2)]
    }



    #print(dim(ed))
    g     <- igraph::graph_from_data_frame(ed, directed = FALSE)
    comps <- igraph::components(g)

    ors_chr <- split(names(comps$membership), comps$membership)
    ors_chr <- ors_chr[vapply(ors_chr, length, integer(1)) >= 10]

    hap_bl <- data.table(HB_id=rep(paste0(ch,"_",1:length(ors_chr)),sapply(ors_chr,length)),marker=unlist(ors_chr),col=rep(rep(sample(col_vector),100)[1:length(ors_chr)],sapply(ors_chr,length)))
    mp[,HB_id:=hap_bl$HB_id[match(marker,hap_bl$marker)]]
    mp[,HB_col:=hap_bl$col[match(marker,hap_bl$marker)]]
    mp[is.na(HB_id),HB_col:="grey"]

    return(mp)
  }))
  return(map_with_hb)
}
