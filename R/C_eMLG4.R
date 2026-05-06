#' Single linkage LD clustering
#'
#' Partitions SNPs into LD clusters using an LD edge list and single-linkage
#' connected components. Edges are retained when pairwise LD and physical distance
#' pass user-supplied thresholds, or thresholds predicted from an LD-decay object.
#'
#' @param ld_decay LD-decay object containing `by_chr` and `decay_sum`. Each
#'   chromosome entry in `ld_decay$by_chr` must contain an edge-list file path in
#'   `el` with columns `SNP1`, `SNP2`, `r2`, and `d`.
#' @param map A `data.table` with at least `Chr`, `Pos`, and `marker` columns.
#' @param SNPs Optional character vector of focal SNPs. If supplied, only LD edges
#'   where `SNP1` or `SNP2` is in `SNPs` are considered.
#' @param rho_ld,rho_d Quantiles/proportions used by `ld_from_rho()` and
#'   `d_from_rho()` when `ld_th` and/or `d_th` are not supplied.
#' @param ld_th Optional fixed LD threshold.
#' @param d_th Optional fixed distance threshold in bp.
#' @param col_vector Optional vector of colours used to label clusters.
#' @param l_min Minimum number of SNPs required for a connected component to be
#'   retained as an LD cluster.
#' @param cores Number of chromosomes to process in parallel.
#'
#' @return A copy of `map` with added columns `CL_id`, `CL_col`, and `n_loci`.
#'
#' @export
LD_clustering <- function(ld_decay,
                          map,
                          SNPs = NULL,
                          rho_ld = 0.95,
                          rho_d = 0.95,
                          ld_th = NULL,
                          d_th = NULL,
                          col_vector = NULL,
                          l_min = 10,
                          cores = 1) {
  data.table::setDT(map)

  if (is.null(col_vector)) {
    col_vector <- c(
      "#B2DF8A", "#FFD92F", "firebrick", "#33A02C", "#7FC97F", "#CAB2D6",
      "#FB8072", "grey30", "#E6AB02", "#FDC086", "steelblue", "#1F78B4",
      "#FB9A99", "#1B9E77", "#BC80BD", "#E31A1C", "#7570B3", "#A6761D",
      "#A6CEE3", "salmon", "#FFFF33", "forestgreen", "#FDCDAC", "#BF5B17",
      "#A6761D", "#FBB4AE", "#4DAF4A", "#B3E2CD", "#FDDAEC", "#BEBADA",
      "#FFF2AE", "#1F78B4", "#66C2A5", "#F0027F", "#E6AB02", "#E78AC3",
      "#FF7F00", "#8DA0CB", "#6A3D9A", "#B15928", "#E41A1C"
    )
  }

  res <- parallel::mclapply(names(ld_decay$by_chr), function(ch) {
    message("Clustering ", ch)

    chr_map <- data.table::copy(map[Chr == ch])
    chr_map[, `:=`(CL_id = NA_character_, CL_col = "grey", n_loci = NA_integer_)]

    chr_obj <- ld_decay$by_chr[[ch]]
    el <- data.table::fread(chr_obj$el, showProgress = FALSE)

    ld_th_chr <- ld_th
    d_th_chr <- d_th

    if (is.null(ld_th_chr)) {
      b_chr <- ld_decay$decay_sum[Chr == ch, b]
      c_chr <- ld_decay$decay_sum[Chr == ch, c_pred]
      ld_th_chr <- ld_from_rho(b_chr, c_chr, rho = rho_ld)
    }

    if (is.null(d_th_chr)) {
      a_chr <- ld_decay$decay_sum[Chr == ch, a_pred]
      d_th_chr <- d_from_rho(a_chr, rho = rho_d)
    }

    ed <- el[r2 > ld_th_chr & d < d_th_chr]

    if (!is.null(SNPs)) {
      ed <- ed[SNP1 %in% SNPs | SNP2 %in% SNPs]
    }

    ed <- ed[, .(SNP1, SNP2)]

    if (!nrow(ed)) {
      return(chr_map)
    }

    g <- igraph::graph_from_data_frame(ed, directed = FALSE)
    comps <- igraph::components(g)
    cls <- split(names(comps$membership), comps$membership)
    cls <- cls[vapply(cls, length, integer(1)) >= l_min]

    if (!length(cls)) {
      return(chr_map)
    }

    cl_cols <- rep(col_vector, length.out = length(cls))

    cl_dt <- data.table::data.table(
      CL_id = rep(paste0(ch, "_", seq_along(cls)), lengths(cls)),
      marker = unlist(cls, use.names = FALSE),
      CL_col = rep(cl_cols, lengths(cls)),
      n_loci = rep(lengths(cls), lengths(cls))
    )

    chr_map[cl_dt, `:=`(
      CL_id = i.CL_id,
      CL_col = i.CL_col,
      n_loci = i.n_loci
    ), on = "marker"]

    chr_map[]
  }, mc.cores = cores)

  data.table::rbindlist(res, fill = TRUE)
}


#' Compute GL-like LD-cluster expected genotypes
#'
#' Collapses SNP genotypes within each LD cluster into an expected genotype dosage.
#' Within each cluster, SNPs are polarized so that negatively correlated SNPs are
#' flipped relative to the most complete reference SNP. For each individual, the
#' empirical frequencies of genotypes 0, 1, and 2 across SNPs in the cluster are
#' treated as GL-like state probabilities. The expected multi-locus genotype (eMLG) is then computed
#' as `P1 + 2 * P2`. An eMLG represents the expected genotype dosage of an LD-cluster
#' derived from empirical genotype-state probabilities across loci.
#'
#' @param GTs Genotype matrix with individuals in rows and SNPs in columns. Values
#'   should be coded as 0, 1, 2, with missing values allowed.
#' @param map A `data.table` containing at least `marker`, `CL_id`, `Chr`, and
#'   `Pos`. Usually this is the output from `LD_clustering()`.
#' @param ncores Number of LD clusters to process in parallel.
#' @param return_probs Logical. If `TRUE`, return per-cluster GL-like probability
#'   arrays in addition to expected genotypes.
#'
#' @return A list with `eMLG`, a matrix of expected cluster genotypes, and `map`,
#'   a cluster-level map. If `return_probs = TRUE`, the list also contains
#'   `GL_like`, a named list of matrices with columns `P0`, `P1`, and `P2`.
#'
#' @export
compute_eMLGs <- function(GTs,
                          map,
                          ncores = 1,
                          return_probs = FALSE,
                          input = c("hard", "dosage")) {
  input <- match.arg(input)

  data.table::setDT(map)
  map <- map[!is.na(CL_id)]

  cls <- split(map$marker, map$CL_id)
  cls <- cls[lengths(cls) > 0L]

  if (input == "auto") {
    input <- detect_gt_input(GTs[, unique(unlist(cls, use.names = FALSE)), drop = FALSE])
    message("Detected genotype input type: ", input)
  }

  cluster_expected_dosage_fast <- function(x) {
    x <- as.matrix(x)

    E <- rowMeans(x, na.rm = TRUE)
    E[is.nan(E)] <- NA_real_

    cbind(
      P0 = NA_real_,
      P1 = NA_real_,
      P2 = NA_real_,
      E = E
    )
  }

  cluster_expected_gt_fast <- function(x) {
    x <- as.matrix(x)

    n <- rowSums(!is.na(x))

    c0 <- rowSums(x == 0, na.rm = TRUE)
    c1 <- rowSums(x == 1, na.rm = TRUE)
    c2 <- rowSums(x == 2, na.rm = TRUE)

    P0 <- c0 / n
    P1 <- c1 / n
    P2 <- c2 / n
    E  <- (c1 + 2 * c2) / n

    P0[n == 0] <- NA_real_
    P1[n == 0] <- NA_real_
    P2[n == 0] <- NA_real_
    E[n == 0]  <- NA_real_

    cbind(P0 = P0, P1 = P1, P2 = P2, E = E)
  }

  process_block <- function(markers) {
    x <- GTs[, markers, drop = FALSE]
    x <- polarize_genotypes(x)

    if (input == "hard") {
      cluster_expected_gt_fast(x)
    } else {
      cluster_expected_dosage_fast(x)
    }
  }

  if (ncores > 1L) {
    res_list <- parallel::mclapply(cls, process_block, mc.cores = ncores)
  } else {
    res_list <- vector("list", length(cls))

    for (i in seq_along(cls)) {
      if (i %% max(1, floor(length(cls) / 20)) == 0L || i == length(cls)) {
        message(sprintf(
          "Processing cluster %d / %d (%.1f%%)",
          i, length(cls), 100 * i / length(cls)
        ))
      }
      res_list[[i]] <- process_block(cls[[i]])
    }
  }

  names(res_list) <- names(cls)

  eMLG <- do.call(cbind, lapply(res_list, function(z) z[, "E"]))
  colnames(eMLG) <- names(res_list)
  rownames(eMLG) <- rownames(GTs)

  CL_map <- map[
    ,
    {
      pos <- as.numeric(Pos)

      .(
        Chr = Chr[1],
        Pos_min = min(pos, na.rm = TRUE),
        Pos_max = max(pos, na.rm = TRUE),
        Pos_mid = median(pos, na.rm = TRUE),
        span_bp = max(pos, na.rm = TRUE) - min(pos, na.rm = TRUE),
        N = .N,
        col = if ("CL_col" %in% names(.SD)) CL_col[1] else NA_character_
      )
    },
    by = CL_id
  ]


  out <- list(
    eMLG = eMLG,
    map = CL_map[match(colnames(eMLG), CL_id)]
  )

  if (isTRUE(return_probs)) {
    if (input == "hard") {
      out$GL_like <- lapply(
        res_list,
        function(z) z[, c("P0", "P1", "P2"), drop = FALSE]
      )
    } else {
      warning("`return_probs = TRUE` is only meaningful for `input = 'hard'`.")
    }
  }

  out
}

detect_gt_input <- function(x, tol = 1e-8) {
  vals <- as.numeric(x)
  vals <- vals[is.finite(vals)]

  if (!length(vals)) {
    return("hard")
  }

  is_hard <- all(abs(vals - round(vals)) < tol) &&
    all(vals %in% c(0, 1, 2))

  if (is_hard) "hard" else "dosage"
}


#' Recursively collapse expected multi-locus genotypes (eMLGs) into composite eMLGs
#'
#' Recursively partitions expected LD-cluster genotypes into composite eMLGs
#' (`C-eMLGs`). Candidate groups are polarized, collapsed into GL-like expected
#' genotypes, and accepted if the mean maximum GL-like state probability exceeds `min_r2_threshold`.
#' If not accepted, the candidate group is split by hierarchical clustering
#' on `1 - abs(correlation)`.
#'
#' @param eMLG Matrix of expected LD-cluster genotypes with individuals in rows
#'   and LD clusters in columns.
#' @param cor_threshold Optional initial absolute-correlation threshold used to
#'   split MLGs into coarse connected components before recursive refinement. Set
#'   to `NULL` to start from one cluster containing all supplied columns.
#' @param min_r2_threshold Minimum accepted average of the per-individual
#'   maximum GL-like state probability.
#' @param min_size Minimum number of eMLGs below which recursive splitting stops.
#' @param hclust_method Hierarchical clustering method. Default is `"single"`.
#' @param prefix Prefix used for composite MLG IDs.
#'
#' @return A list with `C-eMLG`, a collapsed composite expected-genotype matrix, and
#'   `map` containing `C_MLG_id`, `n_eMLGs`,
#'    `min_r2_to_members`, `passed`, and `members`.
#'
#' @export
recursive_collapse_eMLGs <- function(eMLG,
                                     cor_threshold = 0.9,
                                     min_r2_threshold = 0.95,
                                     min_size = 2,
                                     hclust_method = "single",
                                     prefix = "C_eMLG") {
  eMLG <- as.matrix(eMLG)

  if (is.null(rownames(eMLG))) {
    rownames(eMLG) <- paste0("ind_", seq_len(nrow(eMLG)))
  }
  if (is.null(colnames(eMLG))) {
    colnames(eMLG) <- paste0("MLG_", seq_len(ncol(eMLG)))
  }

  score_cluster <- function(cols) {
    x <- eMLG[, cols, drop = FALSE]
    x <- polarize_genotypes(x)

    if (length(cols) == 1L) {
      y <- as.numeric(x[, 1])
      return(list(
        C_eMLG = y,
        mean_r2_to_members = 1,
        min_r2_to_members = 1,
        ok = TRUE
      ))
    }

    y <- rowMeans(x, na.rm = TRUE)
    y[is.nan(y)] <- NA_real_

    r2_to_members <- suppressWarnings(
      stats::cor(y, x, use = "pairwise.complete.obs")^2
    )

    list(
      C_eMLG = y,
      mean_r2_to_members = mean(r2_to_members, na.rm = TRUE),
      min_r2_to_members = min(r2_to_members, na.rm = TRUE),
      ok = min(r2_to_members, na.rm = TRUE) >= min_r2_threshold
    )
  }

  split_cluster <- function(cols) {
    sc <- score_cluster(cols)

    if (sc$ok || length(cols) <= min_size) {
      return(list(list(cols = cols, score = sc)))
    }

    x <- polarize_genotypes(eMLG[, cols, drop = FALSE])
    cc <- suppressWarnings(stats::cor(x, use = "pairwise.complete.obs"))
    cc[!is.finite(cc)] <- 0
    diag(cc) <- 1

    hc <- stats::hclust(stats::as.dist(1 - abs(cc)), method = hclust_method)
    cl <- stats::cutree(hc, k = 2)

    child1 <- cols[cl == 1]
    child2 <- cols[cl == 2]

    if (!length(child1) || !length(child2)) {
      return(list(list(cols = cols, score = sc)))
    }

    c(split_cluster(child1), split_cluster(child2))
  }

  initial_clusters <- initial_eMLG_components(eMLG, cor_threshold = cor_threshold)
  final_blocks <- unlist(lapply(initial_clusters, split_cluster), recursive = FALSE)

  block_scores <- lapply(final_blocks, `[[`, "score")
  block_ids <- paste0(prefix, "_", seq_along(final_blocks))

  C_eMLG <- do.call(cbind, lapply(block_scores, function(sc) sc$C_eMLG))

  colnames(C_eMLG) <- block_ids
  rownames(C_eMLG) <- rownames(eMLG)

  cmlg_map <- data.table::data.table(
    C_eMLG_id = block_ids,
    n_MLGs = vapply(final_blocks, function(z) length(z$cols), integer(1)),
    min_r2_to_members = vapply(block_scores, function(sc) sc$min_r2_to_members, numeric(1)),
    mean_r2_to_members = vapply(block_scores, function(sc) sc$mean_r2_to_members, numeric(1)),
    passed = vapply(block_scores, function(sc) sc$ok, logical(1)),
    members = I(lapply(final_blocks, `[[`, "cols"))
  )

  list(C_eMLG = C_eMLG, map = cmlg_map[])
}


#' Collapse expected MLGs chromosome-wise into composite expected MLGs
#'
#' Applies `recursive_collapse_eMLGs()` separately within chromosomes and,
#' optionally, within distance-delimited genomic windows/components. This prevents
#' composite MLGs from spanning chromosomes and can prevent very long-range
#' composite MLGs within a chromosome.
#'
#' @param MLGs Object with matrix `eMLG` and metadata table `map`. The columns of
#'   `eMLG` must match the IDs in `map[[id_col]]`.
#' @param map_cl Optional original SNP-level cluster map produced by
#'   `LD_clustering()`. If supplied, output ranges and marker lists are collected
#'   from the original SNP markers.
#' @param id_col Name of the ID column in `MLGs$map`. Default is `"CL_id"`.
#' @param chr_col,pos_col Column names for chromosome and position.
#' @param distance_threshold Optional maximum bp distance used to split ordered
#'   MLGs into independent runs. A new run starts whenever the gap between adjacent
#'   MLG positions exceeds this threshold.
#' @param cores Number of chromosome/window groups to process in parallel.
#' @param ... Additional arguments passed to `recursive_collapse_eMLGs()`.
#'
#' @return A list with chromosome-wise collapsed `C_eMLG` and combined `map`.
#'
#' @export
collapse_eMLGs_by_chr <- function(MLGs,
                                  map_cl = NULL,
                                  id_col = "CL_id",
                                  chr_col = "Chr",
                                  pos_col = "Pos",
                                  distance_threshold = NULL,
                                  min_r2_threshold=0.75,
                                  min_size = 2,
                                  hclust_method = "single",
                                  cores = 1,
                                  ...) {
  eMLG <- as.matrix(MLGs$eMLG)
  map <- data.table::copy(MLGs$map)

  required <- c(id_col, chr_col, pos_col)
  missing <- setdiff(required, names(map))
  if (length(missing)) {
    stop("Missing required columns in `map`: ", paste(missing, collapse = ", "))
  }

  map <- data.table::copy(map[get(id_col) %in% colnames(eMLG)])
  data.table::setorderv(map, c(chr_col, pos_col))

  map[, .run := {
    pos <- as.numeric(get(pos_col))
    if (is.null(distance_threshold)) {
      rep(1L, .N)
    } else {
      cumsum(c(TRUE, diff(pos) > distance_threshold))
    }
  }, by = chr_col]

  groups <- split(map, by = c(chr_col, ".run"), keep.by = TRUE)
  #length(groups)
  #i <- 1
  res <- parallel::mclapply(seq_along(groups), function(i) {
    gmap <- groups[[i]]
    ids <- gmap[[id_col]]
    chr <- as.character(gmap[[chr_col]][1])

    message("Collapsing ", chr, ", run ", gmap$.run[1])

    out <- recursive_collapse_eMLGs(
      eMLG = eMLG[, ids, drop = FALSE],
      prefix = paste0("C_eMLG_", chr),
      min_r2_threshold,
      cor_threshold,
      min_size,
      hclust_method
    )

    if (!is.null(map_cl)) {
      range_map <- out$map[
        ,
        {
          member_ids <- unlist(members, use.names = FALSE)
          m <- map_cl[get(id_col) %in% member_ids]

          .(
            Chr = chr,
            Pos_min = min(m[[pos_col]], na.rm = TRUE),
            Pos_max = max(m[[pos_col]], na.rm = TRUE),
            markers = list(m$marker),
            tot_loci = nrow(m),
            group_id = paste0(chr, "_", gmap$.run[1])
          )
        },
        by = C_eMLG_id
      ]
    } else {
      range_map <- out$map[
        ,
        {
          member_ids <- unlist(members, use.names = FALSE)
          m <- gmap[get(id_col) %in% member_ids]

          .(
            Chr = chr,
            Pos_min = min(m[[pos_col]], na.rm = TRUE),
            Pos_max = max(m[[pos_col]], na.rm = TRUE),
            markers = list(m[[id_col]]),
            tot_loci = nrow(m),
            group_id = paste0(chr, "_", gmap$.run[1])
          )
        },
        by = C_eMLG
      ]
    }

    out$map[range_map, `:=`(
      Chr = i.Chr,
      Pos_min = i.Pos_min,
      Pos_max = i.Pos_max,
      markers = i.markers,
      tot_loci = i.tot_loci,
      group_id = i.group_id
    ), on = "C_eMLG_id"]

    out
  }, mc.cores = cores)


  collapsed <- do.call(cbind, lapply(res, `[[`, "C_eMLG"))
  combined_map <- data.table::rbindlist(lapply(res, `[[`, "map"), fill = TRUE)

  list(C_eMLG = collapsed, map = combined_map[])
}


#' Compute squared correlations among expected genotypes
#'
#' Computes an ngsLD-style LD matrix as squared Pearson correlations among
#' expected genotypes.
#'
#' @param eMLG Matrix of expected genotypes with individuals in rows and loci or
#'   clusters in columns.
#'
#' @return A numeric matrix of pairwise `r^2` values.
#'
#' @export
r2_expected_genotypes <- function(eMLG) {
  stats::cor(eMLG, use = "pairwise.complete.obs")^2
}


# Internal helpers ---------------------------------------------------------

#' GL-like expected genotype for rows of a genotype matrix
#'
#' For each row, computes empirical state probabilities for genotypes 0, 1, and
#' 2, then returns the expected genotype `P1 + 2 * P2`.
#'
#' @param x Matrix with individuals in rows and loci or clusters in columns.
#'
#' @return A numeric matrix with columns `P0`, `P1`, `P2`, and `E`.
#'
#' @keywords internal
cluster_expected_gt <- function(x) {
  x <- as.matrix(x)
  x <- polarize_genotypes(x)

  n <- rowSums(!is.na(x))

  c0 <- rowSums(x == 0, na.rm = TRUE)
  c1 <- rowSums(x == 1, na.rm = TRUE)
  c2 <- rowSums(x == 2, na.rm = TRUE)

  P0 <- c0 / n
  P1 <- c1 / n
  P2 <- c2 / n
  E  <- (c1 + 2 * c2) / n

  P0[n == 0] <- NA_real_
  P1[n == 0] <- NA_real_
  P2[n == 0] <- NA_real_
  E[n == 0]  <- NA_real_

  cbind(P0 = P0, P1 = P1, P2 = P2, E = E)
}


#' Polarize genotype columns to a common direction
#'
#' Flips genotype columns coded as 0, 1, 2 to `2 - genotype` when they are
#' negatively correlated with the most complete reference column. This is used so
#' that expected-genotype collapsing is not cancelled out by opposite LD phase.
#'
#' @param x Genotype matrix with individuals in rows and loci in columns.
#'
#' @return A polarized genotype matrix.
#'
#' @keywords internal
polarize_genotypes <- function(x,
                               cor_threshold = 0,
                               dosage_flip_threshold = 1.5) {

  x <- as.matrix(x)

  if (ncol(x) <= 1L) {
    return(x)
  }

  ref_idx <- which.max(colSums(!is.na(x)))
  ref <- x[, ref_idx]

  ## first pass: correlation-based flipping
  flip <- vapply(seq_len(ncol(x)), function(j) {
    r <- suppressWarnings(
      stats::cor(ref, x[, j], use = "pairwise.complete.obs")
    )

    is.finite(r) && r < cor_threshold
  }, logical(1))

  if (any(flip)) {
    x[, flip] <- 2 - x[, flip]
  }

  ## second pass: mean-dosage harmonization
  mu <- colMeans(x, na.rm = TRUE)

  global_mean <- mean(mu, na.rm = TRUE)

  flip2 <- mu > (2 - global_mean)

  if (any(flip2)) {
    x[, flip2] <- 2 - x[, flip2]
  }

  x
}


#' Initial connected components among expected MLGs
#'
#' Builds initial components from expected MLGs using squared Pearson correlation.
#'
#' @param eMLG Matrix of expected MLGs with individuals in rows and clusters in
#'   columns.
#' @param cor_threshold Optional `r^2` threshold for connecting two columns.
#'
#' @return A list of character vectors containing column names.
#' @keywords internal
initial_eMLG_components <- function(eMLG, cor_threshold = 0.9) {
  cols <- colnames(eMLG)

  if (is.null(cor_threshold) || ncol(eMLG) <= 1L) {
    return(list(cols))
  }

  x <- polarize_genotypes(eMLG)

  cc <- suppressWarnings(stats::cor(x, use = "pairwise.complete.obs")^2)
  cc[!is.finite(cc)] <- 0
  diag(cc) <- 1

  edge_idx <- which(cc >= cor_threshold & upper.tri(cc), arr.ind = TRUE)

  if (!nrow(edge_idx)) {
    return(as.list(cols))
  }

  ed <- data.frame(
    from = cols[edge_idx[, 1]],
    to = cols[edge_idx[, 2]],
    stringsAsFactors = FALSE
  )

  g <- igraph::graph_from_data_frame(ed, directed = FALSE, vertices = cols)
  comps <- igraph::components(g)
  split(names(comps$membership), comps$membership)
}

collapse_eMLGs_hclust2 <- function(eMLG,
                                  r2_threshold = 0.75,
                                  prefix = "C_eMLG") {
  eMLG <- as.matrix(eMLG)

  r2 <- suppressWarnings(stats::cor(eMLG, use = "pairwise.complete.obs")^2)
  r2[!is.finite(r2)] <- 0
  diag(r2) <- 1

  d <- stats::as.dist(1 - r2)
  hc <- stats::hclust(d, method = "complete")

  groups <- split(colnames(eMLG), stats::cutree(hc, h = 1 - r2_threshold))


  scores <- lapply(groups, function(cols) {
    x <- polarize_genotypes(eMLG[, cols, drop = FALSE])

    if (length(cols) == 1L) {
      y <- as.numeric(x[, 1])
      r2_members <- 1
    } else {
      y <- rowMeans(x, na.rm = TRUE)
      y[is.nan(y)] <- NA_real_
      r2_members <- suppressWarnings(stats::cor(y, x, use = "pairwise.complete.obs")^2)
    }

    list(
      C_eMLG = y,
      mean_r2_to_members = mean(r2_members, na.rm = TRUE),
      min_r2_to_members = min(r2_members, na.rm = TRUE)
    )
  })



  C_eMLG <- do.call(cbind, lapply(scores, `[[`, "C_eMLG"))
  colnames(C_eMLG) <- paste0(prefix, "_", seq_along(groups))
  rownames(C_eMLG) <- rownames(eMLG)

  map <- data.table::data.table(
    C_eMLG_id = colnames(C_eMLG),
    n_eMLGs = lengths(groups),
    mean_r2_to_members = vapply(scores, `[[`, numeric(1), "mean_r2_to_members"),
    min_r2_to_members = vapply(scores, `[[`, numeric(1), "min_r2_to_members"),
    members = I(groups)
  )

  list(C_eMLG = C_eMLG, map = map)
}
collapse_eMLGs_hclust <- function(eMLG,
                                  r2_threshold = 0.75,
                                  prefix = "C_eMLG",
                                  method = "complete") {
  eMLG <- as.matrix(eMLG)

  r2 <- suppressWarnings(stats::cor(eMLG, use = "pairwise.complete.obs")^2)
  r2[!is.finite(r2)] <- 0
  diag(r2) <- 1

  d <- stats::as.dist(1 - r2)
  hc <- stats::hclust(d, method = "complete")

  groups <- split(colnames(eMLG), stats::cutree(hc, h = 1 - r2_threshold))


  scores <- lapply(groups, function(cols) {
    x <- polarize_genotypes(eMLG[, cols, drop = FALSE])

    if (length(cols) == 1L) {
      y <- as.numeric(x[, 1])
      r2_members <- 1
    } else {
      y <- rowMeans(x, na.rm = TRUE)
      y[is.nan(y)] <- NA_real_
      r2_members <- suppressWarnings(stats::cor(y, x, use = "pairwise.complete.obs")^2)
    }

    list(
      C_eMLG = y,
      mean_r2_to_members = mean(r2_members, na.rm = TRUE),
      min_r2_to_members = min(r2_members, na.rm = TRUE)
    )
  })



  C_eMLG <- do.call(cbind, lapply(scores, `[[`, "C_eMLG"))
  colnames(C_eMLG) <- paste0(prefix, "_", seq_along(groups))
  rownames(C_eMLG) <- rownames(eMLG)

  map <- data.table::data.table(
    C_eMLG_id = colnames(C_eMLG),
    n_eMLGs = lengths(groups),
    mean_r2_to_members = vapply(scores, `[[`, numeric(1), "mean_r2_to_members"),
    min_r2_to_members = vapply(scores, `[[`, numeric(1), "min_r2_to_members"),
    members = I(groups)
  )

  keep <- map[,min_r2_to_members>=r2_threshold]
  map <- map[keep]
  C_eMLG <- C_eMLG[,keep]

  list(C_eMLG = C_eMLG, map = map)
}

collapse_eMLGs_by_chr_window <- function(MLGs,
                                         distance_threshold = 5e5,
                                         r2_threshold = 0.75,
                                         method = "complete") {
  eMLG <- as.matrix(MLGs$eMLG)
  map <- data.table::copy(MLGs$map)

  data.table::setorderv(map, c("Chr", "Pos_min"))

  map[, gap_bp := data.table::shift(Pos_min, type = "lead") - Pos_max]
  map[, new_run := c(TRUE, head(gap_bp, -1) > distance_threshold)]
  map[, run := cumsum(new_run), by = Chr]

  groups <- split(map, by = c("Chr", "run"), keep.by = TRUE)


  res <- lapply(groups, function(gmap) {

    ids <- gmap$CL_id

    if (length(ids) == 1L) {

      x <- eMLG[, ids, drop = FALSE]

      return(list(
        C_eMLG = x,
        map = data.table::data.table(
          C_eMLG_id = ids,
          n_eMLGs = 1L,
          mean_r2_to_members = 1,
          min_r2_to_members = 1,
          members = list(ids)
        )
      ))
    }

    collapse_eMLGs_hclust(
      eMLG = eMLG[, ids, drop = FALSE],
      r2_threshold = r2_threshold,
      method = method
    )
  })

  C_eMLG <- do.call(
    cbind,
    lapply(res, `[[`, "C_eMLG")
  )

  map <- data.table::rbindlist(
    lapply(res, `[[`, "map"),
    fill = TRUE
  )


  new_ids <- paste0("C_eMLG_", seq_len(ncol(C_eMLG)))

  old_ids <- colnames(C_eMLG)

  colnames(C_eMLG) <- new_ids
  map[, old_C_eMLG_id := C_eMLG_id]
  map[, C_eMLG_id := new_ids]

  return(list(map=map,C_eMLG=C_eMLG))
}

### example use ###
map_CL <- LD_clustering(ld_decay,map=map_DIEM,ld_th = 0.8,d_th=1e6,l_min=3,cores = 4)

eMLGs <- compute_eMLGs(GTs=GTs_DIEM,map = map_CL[n_loci>=10,],ncores = 8,return_probs = TRUE,input = "hard")

ce <- collapse_eMLGs_by_chr_window(MLGs=eMLGs,distance_threshold = 5e5,r2_threshold=0.75,method = "complete")
ce$map[n_eMLGs>1 ,hist(min_r2_to_members)]

global_C <- collapse_eMLGs_hclust(eMLG = ce$C_eMLG, r2_threshold = 0.75,method = "complete")

build_C_eMLG_marker_lookup <- function(global_map,
                                       local_map,
                                       map_CL,
                                       global_id_col = "C_eMLG_id",
                                       local_id_col = "C_eMLG_id",
                                       cl_id_col = "CL_id") {

  g2local <- global_map[
    ,
    .(local_C_eMLG_id = unlist(members, use.names = FALSE)),
    by = global_id_col
  ]

  l2cl <- local_map[
    ,
    .(CL_id = unlist(members, use.names = FALSE)),
    by = local_id_col
  ]

  data.table::setnames(l2cl, local_id_col, "local_C_eMLG_id")

  out <- g2local[l2cl, on = "local_C_eMLG_id", allow.cartesian = TRUE]
  out <- out[map_CL, on = c(CL_id = cl_id_col), nomatch = 0]

  out[]
}

run_ohta_for_global_C_eMLG <- function(global_id,
                                       marker_lookup,
                                       GTs,
                                       sample_info,
                                       keep_inds,
                                       pop_col = "Population",
                                       n_pairs = 10000,
                                       cores = 1) {

  message("Processing global C-eMLG ", global_id)

  mrk_dt <- marker_lookup[C_eMLG_id == global_id]
  mrks <- intersect(mrk_dt$marker, colnames(GTs))

  if (length(mrks) < 2L) {
    return(NULL)
  }

  gts <- GTs[keep_inds, mrks, drop = FALSE]
  gts <- polarize_genotypes(gts)

  rownames(gts) <- rownames(GTs)[keep_inds]

  prep <- ohta_fast_prepare(
    data_set = gts,
    pops = sample_info[keep_inds, get(pop_col)]
  )

  pairs <- t(utils::combn(ncol(gts), 2))

  local_ids <- mrk_dt[match(colnames(gts), marker), local_C_eMLG_id]

  pairs <- pairs[local_ids[pairs[, 1]] != local_ids[pairs[, 2]], , drop = FALSE]

  if (!nrow(pairs)) {
    return(NULL)
  }

  if (nrow(pairs) > n_pairs) {
    set.seed(abs(as.integer(factor(global_id))))
    pairs <- pairs[sort(sample(seq_len(nrow(pairs)), n_pairs)), , drop = FALSE]
  }

  ohta <- dstat_unphased_scan(
    pairs = pairs,
    prep = prep,
    tot_maf = 0.01,
    pop_maf = 0.01,
    cores = cores
  )

  ohta <- data.table::data.table(
    global_C_eMLG_id = global_id,
    from = colnames(gts)[pairs[, 1]],
    to = colnames(gts)[pairs[, 2]],
    ohta
  )

  annot <- unique(marker_lookup[, .(
    marker,
    Chr,
    Pos,
    local_C_eMLG_id,
    CL_id
  )])

  ohta[annot, `:=`(
    Chr1 = i.Chr,
    Pos1 = i.Pos,
    local_C_eMLG_id1 = i.local_C_eMLG_id,
    CL_id1 = i.CL_id
  ), on = c(from = "marker")]

  ohta[annot, `:=`(
    Chr2 = i.Chr,
    Pos2 = i.Pos,
    local_C_eMLG_id2 = i.local_C_eMLG_id,
    CL_id2 = i.CL_id
  ), on = c(to = "marker")]

  ohta[, Same_chr := Chr1 == Chr2]
  ohta[, local_C_eMLG_pair := paste(local_C_eMLG_id1, local_C_eMLG_id2, sep = " | ")]

  ohta[]
}

marker_lookup <- build_C_eMLG_marker_lookup(
  global_map = global_C$map,
  local_map = ce$map,
  map_CL = map_CL
)

#image(t(GTs[,marker_lookup[C_eMLG_id=="C_eMLG_2049",marker]]))

marker_lookup
ohta_one <- run_ohta_for_global_C_eMLG(
  global_id = global_C$map$C_eMLG_id[1],
  marker_lookup = marker_lookup,
  GTs = GTs_DIEM,
  sample_info = sample_info,
  keep_inds = keep_inds,
  n_pairs = 10000,
  cores = 8
)
