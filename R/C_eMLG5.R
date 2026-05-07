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
                          cor_th = 0.8,
                          return_probs = FALSE,
                          input = c("hard", "dosage")) {
  input <- match.arg(input)
  #map <- data.table::copy(map_CL[n_loci>=10])
  data.table::setDT(map)
  map <- map[!is.na(CL_id)]

  cls <- split(map$marker, map$CL_id)
  cls <- cls[lengths(cls) > 0L]


  if (input == "auto") {
    input <- detect_gt_input(GTs[, unique(unlist(cls, use.names = FALSE)), drop = FALSE])
    message("Detected genotype input type: ", input)
  }

  cluster_expected_dosage <- function(x) {
    x <- as.matrix(x)

    E <- rowMeans(x, na.rm = TRUE)
    E[is.nan(E)] <- NA_real_

    # cbind(
    #   P0 = NA_real_,
    #   P1 = NA_real_,
    #   P2 = NA_real_,
    #   E = E
    # )
    E
  }

  cluster_expected_gt <- function(x) {
    #x <- as.matrix(x)

    n <- rowSums(!is.na(x))

   # c0 <- rowSums(x == 0, na.rm = TRUE)
    c1 <- rowSums(x == 1, na.rm = TRUE)
    c2 <- rowSums(x == 2, na.rm = TRUE)

    # P0 <- c0 / n
    # P1 <- c1 / n
    # P2 <- c2 / n
    E  <- (c1 + 2 * c2) / n

    # P0[n == 0] <- NA_real_
    # P1[n == 0] <- NA_real_
    # P2[n == 0] <- NA_real_
    E[n == 0]  <- NA_real_

    # cbind(P0 = P0, P1 = P1, P2 = P2, E = E)
    E
  }
  #markers <- cls["Chr18_655"][[1]]
  process_block <- function(markers) {
    x <- GTs[, markers, drop = FALSE]

    if (input == "hard") {
      cluster_expected_gt(x)
    } else {
      cluster_expected_dosage(x)
    }
  }
  #input ="hard"
  #markers <- cls[[1]]
  #ncores = 4
  if (ncores > 1L) {
    eMLG <- do.call(cbind,parallel::mclapply(cls, process_block, mc.cores = ncores))

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
    eMLG <- do.call(cbind,res_list)
  }

  colnames(eMLG) <- names(cls)

  split_cl <- function(markers,th){
    x <- polarize_genotypes(GTs[, markers, drop = FALSE])


    cc <- suppressWarnings(stats::cor(
      x,
      use = "pairwise.complete.obs"
    )^2)


    cc[!is.finite(cc)] <- 0
    diag(cc) <- 1


    hc <- stats::hclust(
      stats::as.dist(1 - cc),
      method = "complete"
    )

    split(markers, stats::cutree(hc, h = 1 - th))
  }

  ## one round of splitting clusters with too weak correlation between expected genotypes and GL
  cors <- apply(eMLG,2,function(x)  suppressWarnings(cor(round(x),x,use= "pairwise.complete.obs")^2))

  redo <- which(cors<cor_th & !is.na(cors))

  cls_split <- unlist(lapply(cls[redo],split_cl,th=cor_th),recursive = FALSE)
  eMLG_split <- do.call(cbind,parallel::mclapply(cls_split, process_block, mc.cores = ncores))
  colnames(eMLG_split) <- names(cls_split)

  cors <- apply(eMLG_split,2,function(x)  suppressWarnings(cor(round(x),x,use= "pairwise.complete.obs")^2))
  keep <- which(cors>cor_th & !is.na(cors) & lengths(cls_split)>10)

  ## combina with good clusters
  eMLG_final <- cbind(eMLG[,-redo],eMLG_split[,keep])
  cls_final <- c(cls[-redo],cls_split[keep])

  ## redo the map
  cl_cols <- rep(col_vector, length.out = length(cls_final))

  cors <- apply(eMLG_final,2,function(x)  suppressWarnings(cor(round(x),x,use= "pairwise.complete.obs")^2))

  cl_dt <- data.table::data.table(
    CL_id    = rep(names(cls_final), lengths(cls_final)),
    marker   = unlist(cls_final, use.names = FALSE),
    CL_col   = rep(cl_cols, lengths(cls_final)),
    r2_eMLG  = rep(cors, lengths(cls_final)),
    n_loci   = rep(lengths(cls_final), lengths(cls_final))
  )

  map[,CL_id := NULL]
  map[cl_dt, `:=`(
    CL_id   = i.CL_id,
    CL_col  = i.CL_col,
    n_loci  = i.n_loci,
    r2_eMLG = i.r2_eMLG
  ), on = "marker"]


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
        r2_eMLG = r2_eMLG[1],
        col = if ("CL_col" %in% names(.SD)) CL_col[1] else NA_character_
      )
    },
    by = CL_id
  ]

  out <- list(
    eMLG = eMLG_final,
    map_eMLG = CL_map[match(colnames(eMLG_final), CL_id)],
    map_SNP = cl_dt
  )

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
                               cor_threshold = 0#,
                               #dosage_flip_threshold = 1.5
                               ) {

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
  # mu <- colMeans(x, na.rm = TRUE)
  #
  # global_mean <- mean(mu, na.rm = TRUE)
  #
  # flip2 <- mu > (2 - global_mean)
  #
  # if (any(flip2)) {
  #   x[, flip2] <- 2 - x[, flip2]
  # }

  x
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
#method = "complete"
#collapse_eMLGs_hclust

#Chr11_344,Chr11_1300
#map <- map_CL
collapse_eMLGs_hclust <- function(MLGs,
                                  map_cl,
                                  r2_threshold = 0.75,
                                  prefix = "C_eMLG",
                                  method = "complete") {
  eMLG <- as.matrix(MLGs$eMLG)

  map <- MLGs$map

  r2 <- suppressWarnings(stats::cor(eMLG, use = "pairwise.complete.obs")^2)
  r2[!is.finite(r2)] <- 0
  diag(r2) <- 1

  d <- stats::as.dist(1 - r2)
  hc <- stats::hclust(d, method = method)

  groups <- split(colnames(eMLG), stats::cutree(hc, h = 1 - r2_threshold))
  #which(sapply(groups, function(x) any(x %in% "Chr11_344")))

  #markers <- groups[[1]]
  score_group <- function(markers) {
    x <- eMLG[, markers, drop = FALSE]


    if (length(markers) == 1L) {
      y <- as.numeric(x[, 1])
      r2_members <- 1
    } else {

      w <- map_cl[match(markers, CL_id), n_loci]
      w[is.na(w)] <- 1
      y <- weighted_row_mean(polarize_genotypes(x), w)

      r2_members <- suppressWarnings(
        stats::cor(x, use = "pairwise.complete.obs")^2
      )
    }
    r2_eMLG <- suppressWarnings(
      stats::cor(round(y),y, use = "pairwise.complete.obs")^2
    )
    list(
      C_eMLG = y,
      mean_r2_to_members = mean(r2_members, na.rm = TRUE),
      min_r2_to_members = min(r2_members, na.rm = TRUE),
      r2_eMLG = r2_eMLG
    )
  }

  scores <- lapply(groups, score_group)

  C_eMLG <- do.call(cbind, lapply(scores, `[[`, "C_eMLG"))
  colnames(C_eMLG) <- paste0(prefix, "_", seq_along(groups))
  rownames(C_eMLG) <- rownames(eMLG)

  map_out <- data.table::data.table(
    C_eMLG_id = colnames(C_eMLG),
    n_eMLGs = lengths(groups),
    mean_r2_to_members = vapply(scores, `[[`, numeric(1), "mean_r2_to_members"),
    min_r2_to_members = vapply(scores, `[[`, numeric(1), "min_r2_to_members"),
    r2_eMLG = vapply(scores, `[[`, numeric(1), "r2_eMLG"),
    members = I(groups)
  )


  list(C_eMLG = C_eMLG, map = map_out)
}

collapse_eMLGs_by_chr_window <- function(MLGs,
                                         map_cl,
                                         distance_threshold = 5e5,
                                         r2_threshold = 0.75,
                                         method = "complete") {
  eMLG <- as.matrix(MLGs$eMLG)
  map <- data.table::copy(MLGs$map_eMLG)
  map_cl <- data.table::copy(MLGs$map_SNP)

  data.table::setorderv(map, c("Chr", "Pos_min"))

  map[, gap_bp := data.table::shift(Pos_min, type = "lead") - Pos_max]
  map[, new_run := c(TRUE, head(gap_bp, -1) > distance_threshold)]
  map[, run := cumsum(new_run), by = Chr]

  groups <- split(map, by = c("Chr", "run"), keep.by = TRUE)
  gmap <- groups$Chr11.0
  #"Chr11_344"
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

    out <- collapse_eMLGs_hclust(
      MLGs = list(eMLG=eMLG[, ids, drop = FALSE],map=map),
      map_cl=map_cl[CL_id %in% ids],
      r2_threshold = r2_threshold,
      method = method
    )

  })

  C_eMLG <- do.call(
    cbind,
    lapply(res, `[[`, "C_eMLG")
  )

  map_out <- data.table::rbindlist(
    lapply(res, `[[`, "map"),
    fill = TRUE
  )


  new_ids <- paste0("C_eMLG_", seq_len(ncol(C_eMLG)))

  old_ids <- colnames(C_eMLG)

  colnames(C_eMLG) <- new_ids
  map_out[, old_C_eMLG_id := C_eMLG_id]
  map_out[, C_eMLG_id := new_ids]


  return(list(map_map_eMLG=map_out,C_eMLG=C_eMLG))
}

### example use ###
map_CL <- LD_clustering(ld_decay,map=map_DIEM,ld_th = 0.8,d_th=1e6,l_min=3,cores = 4)

eMLGs <- compute_eMLGs(GTs=GTs_DIEM,map = map_CL[n_loci>=10,],ncores = 8,return_probs = TRUE,input = "hard",cor_th = 0.8)
names(eMLGs)
map_CL[,CL_id := NULL]

table(eMLGs$map_SNP$CL_id %in% eMLGs$map_eMLG$CL_id)

map_CL[eMLGs$map_SNP, `:=`(
  CL_id   = i.CL_id,
  CL_col  = i.CL_col,
  n_loci  = i.n_loci,
  r2_eMLG = i.r2_eMLG
), on = "marker"]


ce <- collapse_eMLGs_by_chr_window(MLGs=eMLGs,map_cl=map_CL,distance_threshold = 5e5,r2_threshold=0.8,method = "complete")


r2 <- suppressWarnings(stats::cor(ce$C_eMLG, use = "pairwise.complete.obs")^2)
r2[!is.finite(r2)] <- 0
diag(r2) <- 1

d <- stats::as.dist(1 - r2)
hc <- stats::hclust(d, method = method)

groups <- split(colnames(ce$C_eMLG), stats::cutree(hc, h = 1 - 0.8))
LD_cls <- groups[lengths(groups)>1]

LD_cl_global_lookup = data.table(LD_cl_global=rep(ce$map$C_eMLG_id,ce$map$n_eMLGs),CL_id=unlist(ce$map$members))

GTs[,map_CL[CL_id %in% LD_cl_global_lookup[LD_cl_global %in% LD_cls[[5]],CL_id],marker]]

image(t(polarize_genotypes(GTs[,map_CL[CL_id %in% LD_cl_global_lookup[LD_cl_global %in% LD_cls[[5]],CL_id],marker]])))




global_C <- collapse_eMLGs_hclust(MLG = ce$C_eMLG, r2_threshold = 0.75,method = "complete")



keep <- which(ce$map[,n_eMLGs>1])
plot(ce$map[keep,min_r2_to_members],apply(ce$C_eML[,keep],2,function(x)cor(round(x),x,use = "pairwise.complete.obs")^2),ylim = c(0.9,1))

which(apply(ce$C_eML,2,function(x)cor(round(x),x,use = "pairwise.complete.obs")^2)<0.5)

ce$map[C_eMLG_id=="C_eMLG_7"]



#eMLGs$eMLG[,ce$map[n_eMLGs>1 & min_r2_to_members<0.75]
image(t(cor(eMLGs$eMLG[,ce$map[n_eMLGs>1 & min_r2_to_members<0.75][,members][[1]]])^2))

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
