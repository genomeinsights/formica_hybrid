# Helpers -----------------------------------------------------------------

default_cluster_colours <- function() {
  c(
    "#B2DF8A", "#FFD92F", "firebrick", "#33A02C", "#7FC97F", "#CAB2D6",
    "#FB8072", "grey30", "#E6AB02", "#FDC086", "steelblue", "#1F78B4",
    "#FB9A99", "#1B9E77", "#BC80BD", "#E31A1C", "#7570B3", "#A6761D",
    "#A6CEE3", "salmon", "#FFFF33", "forestgreen", "#FDCDAC", "#BF5B17",
    "#A6761D", "#FBB4AE", "#4DAF4A", "#B3E2CD", "#FDDAEC", "#BEBADA",
    "#FFF2AE", "#1F78B4", "#66C2A5", "#F0027F", "#E6AB02", "#E78AC3",
    "#FF7F00", "#8DA0CB", "#6A3D9A", "#B15928", "#E41A1C"
  )
}


detect_gt_input <- function(x, tol = 1e-8) {
  vals <- as.numeric(x)
  vals <- vals[is.finite(vals)]

  if (!length(vals)) return("hard")

  is_hard <- all(abs(vals - round(vals)) < tol) &&
    all(vals %in% c(0, 1, 2))

  if (is_hard) "hard" else "dosage"
}


polarize_genotypes <- function(x, cor_threshold = 0) {
  x <- as.matrix(x)

  if (ncol(x) <= 1L) return(x)

  ref_idx <- which.max(colSums(!is.na(x)))
  ref <- x[, ref_idx]

  flip <- vapply(seq_len(ncol(x)), function(j) {
    r <- suppressWarnings(
      stats::cor(ref, x[, j], use = "pairwise.complete.obs")
    )
    is.finite(r) && r < cor_threshold
  }, logical(1))

  if (any(flip)) {
    x[, flip] <- 2 - x[, flip]
  }

  x
}


expected_gt_hard <- function(x) {
  x <- polarize_genotypes(x)

  n <- rowSums(!is.na(x))
  c1 <- rowSums(x == 1, na.rm = TRUE)
  c2 <- rowSums(x == 2, na.rm = TRUE)

  E <- (c1 + 2 * c2) / n
  E[n == 0] <- NA_real_

  E
}


expected_gt_dosage <- function(x) {
  x <- polarize_genotypes(x)

  E <- rowMeans(x, na.rm = TRUE)
  E[is.nan(E)] <- NA_real_

  E
}


score_eMLG <- function(x) {
  r2 <- suppressWarnings(
    stats::cor(round(x), x, use = "pairwise.complete.obs")^2
  )

  if (!is.finite(r2)) NA_real_ else r2
}


weighted_row_mean <- function(x, w) {
  x <- as.matrix(x)
  w <- as.numeric(w)

  ok <- !is.na(x)
  num <- rowSums(sweep(x, 2, w, `*`), na.rm = TRUE)
  den <- rowSums(sweep(ok, 2, w, `*`), na.rm = TRUE)

  y <- num / den
  y[den == 0] <- NA_real_

  y
}


cluster_level_map <- function(map_snp, id_col = "CL_id") {
  data.table::setDT(map_snp)

  map_snp[
    !is.na(get(id_col)),
    {
      pos <- as.numeric(Pos)

      .(
        Chr = Chr[1],
        Pos_min = min(pos, na.rm = TRUE),
        Pos_max = max(pos, na.rm = TRUE),
        Pos_mid = median(pos, na.rm = TRUE),
        span_bp = max(pos, na.rm = TRUE) - min(pos, na.rm = TRUE),
        n_loci = .N,
        r2_eMLG = if ("r2_eMLG" %in% names(.SD)) r2_eMLG[1] else NA_real_,
        col = if ("CL_col" %in% names(.SD)) CL_col[1] else NA_character_
      )
    },
    by = id_col
  ]
}

add_snp_positions <- function(map_snp, map_ref) {
  map_ref <- data.table::as.data.table(map_ref)
  map_snp <- data.table::as.data.table(map_snp)

  pos_cols <- intersect(c("marker", "Chr", "Pos"), names(map_ref))

  merge(
    unique(map_ref[, ..pos_cols]),
    map_snp,
    by = "marker",
    all.y = TRUE,
    sort = FALSE
  )
}

# Main functions -----------------------------------------------------------------
make_eMLGs <- function(GTs,
                       map_cl,
                       input = c("auto", "hard", "dosage"),
                       cor_th = 0.8,
                       l_min = 10,
                       ncores = 1,
                       col_vector = NULL) {
  input <- match.arg(input)

  if (is.null(col_vector)) {
    col_vector <- default_cluster_colours()
  }

  map_cl <- data.table::copy(map_cl)
  data.table::setDT(map_cl)

  map_cl <- map_cl[!is.na(CL_id) & n_loci >= l_min]

  cls <- split(map_cl$marker, map_cl$CL_id)
  cls <- cls[lengths(cls) >= l_min]

  if (!length(cls)) {
    stop("No LD clusters passed l_min.")
  }

  all_markers <- unique(unlist(cls, use.names = FALSE))

  if (input == "auto") {
    input <- detect_gt_input(GTs[, all_markers, drop = FALSE])
    message("Detected genotype input type: ", input)
  }

  expected_fun <- switch(
    input,
    hard = expected_gt_hard,
    dosage = expected_gt_dosage
  )

  process_cluster <- function(markers) {
    expected_fun(GTs[, markers, drop = FALSE])
  }

  res <- parallel::mclapply(cls, process_cluster, mc.cores = ncores)

  eMLG <- do.call(cbind, res)
  colnames(eMLG) <- names(cls)
  rownames(eMLG) <- rownames(GTs)

  r2_eMLG <- apply(eMLG, 2, score_eMLG)

  cl_cols <- rep(col_vector, length.out = length(cls))

  map_snp <- data.table::data.table(
    CL_id = rep(names(cls), lengths(cls)),
    marker = unlist(cls, use.names = FALSE),
    CL_col = rep(cl_cols, lengths(cls)),
    n_loci = rep(lengths(cls), lengths(cls)),
    r2_eMLG = rep(r2_eMLG, lengths(cls))
  )

  map_snp <- add_snp_positions(map_snp, map_cl)

  map_eMLG <- cluster_level_map(
    map_snp,
    id_col = "CL_id"
  )

  map_eMLG <- map_eMLG[match(colnames(eMLG), CL_id)]


  list(
    eMLG = eMLG,
    map_eMLG = map_eMLG,
    map_SNP = map_snp,
    clusters = cls,
    input = input,
    cor_th = cor_th
  )
}

refine_eMLGs_complete <- function(GTs,
                                  MLGs,
                                  input = "dosage",
                                  cor_th = 0.8,
                                  l_min = 10,
                                  hclust_method = "complete",
                                  ncores = 1,
                                  col_vector = NULL) {
  if (is.null(col_vector)) {
    col_vector <- default_cluster_colours()
  }

  cls <- MLGs$clusters
  eMLG <- MLGs$eMLG
  map_snp <- MLGs$map_SNP

  r2_eMLG <- apply(eMLG, 2, score_eMLG)

  good_ids <- names(r2_eMLG)[!is.na(r2_eMLG) & r2_eMLG >= cor_th]
  bad_ids  <- names(r2_eMLG)[is.na(r2_eMLG) | r2_eMLG < cor_th]

  split_cluster <- function(id) {
    markers <- cls[[id]]

    if (length(markers) < 2L) {
      return(list())
    }

    x <- polarize_genotypes(GTs[, markers, drop = FALSE])

    cc <- suppressWarnings(
      stats::cor(x, use = "pairwise.complete.obs")^2
    )

    cc[!is.finite(cc)] <- 0
    diag(cc) <- 1

    hc <- stats::hclust(
      stats::as.dist(1 - cc),
      method = hclust_method
    )

    gr <- stats::cutree(hc, h = 1 - cor_th)

    out <- split(markers, gr)
    out <- out[lengths(out) >= l_min]

    if (!length(out)) {
      return(setNames(list(markers), id))
    }
    names(out) <- paste0(id, ".", seq_along(out))
    out
  }

  split_bad <- unlist(
    lapply(bad_ids, split_cluster),
    recursive = FALSE
  )

  cls_final <- c(cls[good_ids], split_bad)
  cls_final <- cls_final[lengths(cls_final) >= l_min]

  expected_fun <- switch(
    input,
    hard = expected_gt_hard,
    dosage = expected_gt_dosage
  )

  process_cluster <- function(markers) {
    expected_fun(GTs[, markers, drop = FALSE])
  }

  res <- parallel::mclapply(cls_final, process_cluster, mc.cores = ncores)

  eMLG_final <- do.call(cbind, res)
  colnames(eMLG_final) <- names(cls_final)
  rownames(eMLG_final) <- rownames(GTs)

  r2_final <- apply(eMLG_final, 2, score_eMLG)

  keep <- !is.na(r2_final) & r2_final >= cor_th

  eMLG_final <- eMLG_final[, keep, drop = FALSE]
  cls_final <- cls_final[keep]
  r2_final <- r2_final[keep]

  cl_cols <- rep(col_vector, length.out = length(cls_final))

  map_snp <- data.table::data.table(
    CL_id = rep(names(cls_final), lengths(cls_final)),
    marker = unlist(cls_final, use.names = FALSE),
    CL_col = rep(cl_cols, lengths(cls_final)),
    n_loci = rep(lengths(cls_final), lengths(cls_final)),
    r2_eMLG = rep(r2_final, lengths(cls_final))
  )

  map_snp <- add_snp_positions(map_snp, MLGs$map_SNP)

  map_eMLG <- cluster_level_map(
    map_snp,
    id_col = "CL_id"
  )

  list(
    eMLG = eMLG_final,
    map_SNP = map_snp,
    map_eMLG = map_eMLG,
    clusters = cls_final,
    r2_eMLG = r2_final,
    input = input,
    cor_th = cor_th
  )
}

collapse_eMLGs_by_chr_window <- function(MLGs,
                                         distance_threshold = 5e5,
                                         r2_threshold = 0.8,
                                         method = "complete",
                                         prefix = "C_eMLG") {
  eMLG <- as.matrix(MLGs$eMLG)

  map <- data.table::copy(MLGs$map_eMLG)
  data.table::setDT(map)

  data.table::setorderv(map, c("Chr", "Pos_min"))

  map[, gap_bp := Pos_min - data.table::shift(Pos_max), by = Chr]
  map[, new_run := is.na(gap_bp) | gap_bp > distance_threshold, by = Chr]
  map[, run := cumsum(new_run), by = Chr]

  collapse_run <- function(gmap) {
    ids <- intersect(gmap$CL_id, colnames(eMLG))

    if (length(ids) == 1L) {
      y <- eMLG[, ids]

      return(list(
        C_eMLG = matrix(y, ncol = 1, dimnames = list(rownames(eMLG), ids)),
        map = data.table::data.table(
          old_C_eMLG_id = ids,
          n_eMLGs = 1L,
          mean_r2_to_members = 1,
          min_r2_to_members = 1,
          r2_eMLG = score_eMLG(y),
          members = list(ids)
        )
      ))
    }

    x <- eMLG[, ids, drop = FALSE]

    r2 <- suppressWarnings(
      stats::cor(x, use = "pairwise.complete.obs")^2
    )

    r2[!is.finite(r2)] <- 0
    diag(r2) <- 1

    hc <- stats::hclust(
      stats::as.dist(1 - r2),
      method = method
    )

    groups <- split(ids, stats::cutree(hc, h = 1 - r2_threshold))

    score_group <- function(members) {
      z <- x[, members, drop = FALSE]

      if (length(members) == 1L) {
        y <- as.numeric(z[, 1])
        r2_members <- 1
      } else {
        w <- map[match(members, CL_id), n_loci]
        w[is.na(w)] <- 1

        y <- weighted_row_mean(
          polarize_genotypes(z),
          w
        )

        r2_members <- suppressWarnings(
          stats::cor(z, use = "pairwise.complete.obs")^2
        )
      }

      list(
        y = y,
        n_loci = sum(map[match(members, CL_id), n_loci], na.rm = TRUE),
        mean_r2_to_members = mean(r2_members, na.rm = TRUE),
        min_r2_to_members = min(r2_members, na.rm = TRUE),
        r2_eMLG = score_eMLG(y)
      )
    }

    scores <- lapply(groups, score_group)

    C_eMLG <- do.call(cbind, lapply(scores, `[[`, "y"))
    run_prefix <- paste0(gmap$Chr[1], "_run", gmap$run[1], "_")
    colnames(C_eMLG) <- paste0(run_prefix, seq_along(groups))
    rownames(C_eMLG) <- rownames(eMLG)

    map_out <- data.table::data.table(
      old_C_eMLG_id = colnames(C_eMLG),
      n_eMLGs = lengths(groups),
      n_loci = vapply(scores, `[[`, numeric(1), "n_loci"),
      mean_r2_to_members = vapply(scores, `[[`, numeric(1), "mean_r2_to_members"),
      min_r2_to_members = vapply(scores, `[[`, numeric(1), "min_r2_to_members"),
      r2_eMLG = vapply(scores, `[[`, numeric(1), "r2_eMLG"),
      members = I(groups)
    )
    list(C_eMLG = C_eMLG, map = map_out)
  }

  runs <- split(map, by = c("Chr", "run"), keep.by = TRUE)

  res <- lapply(runs, collapse_run)

  C_eMLG <- do.call(cbind, lapply(res, `[[`, "C_eMLG"))
  map_out <- data.table::rbindlist(lapply(res, `[[`, "map"), fill = TRUE)

  new_ids <- paste0(prefix, "_", seq_len(ncol(C_eMLG)))

  colnames(C_eMLG) <- new_ids
  map_out[, C_eMLG_id := new_ids]

  lookup <- data.table::data.table(
    C_eMLG_id = rep(map_out$C_eMLG_id, map_out$n_eMLGs),
    CL_id = unlist(map_out$members, use.names = FALSE)
  )

  map_snp_out <- merge(
    MLGs$map_SNP,
    lookup,
    by = "CL_id",
    all.x = FALSE,
    all.y = TRUE,
    sort = FALSE
  )

  list(
    eMLG = C_eMLG,
    map_eMLG = map_out,
    map_SNP = map_snp_out,
    lookup = lookup,
    clusters = split(lookup$CL_id, lookup$C_eMLG_id),
    params = list(
      stage = "collapse_eMLGs_by_chr_window",
      distance_threshold = distance_threshold,
      r2_threshold = r2_threshold,
      method = method
    )
  )
}

collapse_eMLGs_global <- function(MLGs,
                                  r2_threshold = 0.8,
                                  method = "complete",
                                  prefix = "G_eMLG") {
  eMLG <- as.matrix(MLGs$eMLG)

  r2 <- suppressWarnings(
    stats::cor(eMLG, use = "pairwise.complete.obs")^2
  )

  r2[!is.finite(r2)] <- 0
  diag(r2) <- 1

  hc <- stats::hclust(
    stats::as.dist(1 - r2),
    method = method
  )

  groups <- split(
    colnames(eMLG),
    stats::cutree(hc, h = 1 - r2_threshold)
  )

  score_group <- function(members) {
    x <- eMLG[, members, drop = FALSE]

    if (length(members) == 1L) {
      y <- as.numeric(x[, 1])
      r2_members <- 1
    } else {
      w <- MLGs$map_eMLG[match(members, C_eMLG_id), n_loci]
      w[is.na(w)] <- 1

      y <- weighted_row_mean(
        polarize_genotypes(x),
        w
      )

      r2_members <- suppressWarnings(
        stats::cor(x, use = "pairwise.complete.obs")^2
      )
    }

    list(
      y = y,
      mean_r2_to_members = mean(r2_members, na.rm = TRUE),
      min_r2_to_members = min(r2_members, na.rm = TRUE),
      r2_eMLG = score_eMLG(y)
    )
  }

  scores <- lapply(groups, score_group)

  G_eMLG <- do.call(cbind, lapply(scores, `[[`, "y"))
  G_ids <- paste0(prefix, "_", seq_along(groups))

  colnames(G_eMLG) <- G_ids
  rownames(G_eMLG) <- rownames(eMLG)

  lookup <- data.table::data.table(
    G_eMLG_id = rep(G_ids, lengths(groups)),
    C_eMLG_id = unlist(groups, use.names = FALSE)
  )

  map_eMLG <- data.table::data.table(
    G_eMLG_id = G_ids,
    n_C_eMLGs = lengths(groups),
    mean_r2_to_members = vapply(scores, `[[`, numeric(1), "mean_r2_to_members"),
    min_r2_to_members = vapply(scores, `[[`, numeric(1), "min_r2_to_members"),
    r2_eMLG = vapply(scores, `[[`, numeric(1), "r2_eMLG"),
    members = I(groups)
  )

  map_snp <- merge(
    MLGs$map_SNP,
    lookup,
    by = "C_eMLG_id",
    all.x = FALSE,
    all.y = TRUE,
    sort = FALSE
  )

  data.table::setcolorder(
    map_snp,
    intersect(
      c(
        "marker", "Chr", "Pos",
        "CL_id", "C_eMLG_id", "G_eMLG_id",
        "CL_col", "n_loci", "r2_eMLG"
      ),
      names(map_snp)
    )
  )

  list(
    eMLG = G_eMLG,
    map_eMLG = map_eMLG,
    map_SNP = map_snp,
    lookup = lookup,
    clusters = split(lookup$C_eMLG_id, lookup$G_eMLG_id),
    params = list(
      stage = "collapse_eMLGs_global",
      r2_threshold = r2_threshold,
      method = method
    )
  )
}

# Example use -----------------------------------------------------------------

# map_CL <- LD_clustering(ld_decay,
#                         map=map_DIEM,
#                         ld_th = 0.8,
#                         d_th=1e6,
#                         l_min=10,
#                         cores = 4)
#
# eMLGs <- make_eMLGs(GTs=GTs_DIEM,
#                     map_cl = map_CL[n_loci>=10,],
#                     cor_th = 0.8,
#                     input = "hard",ncores = 8)
#
# eMLGs_ref <- refine_eMLGs_complete(GTs,
#                                    MLGs = eMLGs,
#                                    cor_th = 0.8,
#                                    l_min = 10,
#                                    ncores = 8)
#
# C_eMLGs <- collapse_eMLGs_by_chr_window(MLGs=eMLGs_ref,
#                              distance_threshold = 5e5,
#                              r2_threshold = 0.8,
#                              method = "complete",
#                              prefix = "C_eMLG")
#
#
# GMLGs <- collapse_eMLGs_global(
#   MLGs = C_eMLGs,
#   r2_threshold = 0.8,
#   method = "complete"
# )
