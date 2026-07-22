# Helpers -----------------------------------------------------------------
col_vector <- c(
  "#B2DF8A", "#FFD92F", "firebrick", "#33A02C", "#7FC97F", "#CAB2D6",
  "#FB8072", "grey30", "#E6AB02", "#FDC086", "steelblue", "#1F78B4",
  "#FB9A99", "#1B9E77", "#BC80BD", "#E31A1C", "#7570B3", "#A6761D",
  "#A6CEE3", "salmon", "#FFFF33", "forestgreen", "#FDCDAC", "#BF5B17",
  "#A6761D", "#FBB4AE", "#4DAF4A", "#B3E2CD", "#FDDAEC", "#BEBADA",
  "#FFF2AE", "#1F78B4", "#66C2A5", "#F0027F", "#E6AB02", "#E78AC3",
  "#FF7F00", "#8DA0CB", "#6A3D9A", "#B15928", "#E41A1C"
)

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


polarize_genotypes <- function(x, R2_th = 0) {
  x <- as.matrix(x)

  if (ncol(x) <= 1L) return(x)

  ref_idx <- which.max(colSums(!is.na(x)))
  ref <- x[, ref_idx]




  flip <- vapply(seq_len(ncol(x)), function(j) {
    r <- suppressWarnings(
      stats::cor(ref, x[, j], use = "pairwise.complete.obs")
    )
    is.finite(r) && r < R2_th
  }, logical(1))

  if (any(flip)) {
    x[, flip] <- 2 - x[, flip]
  }
  return(x)
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
LD_clustering <- function(
    ld_decay,
    #el_folder = "./EL/",
    map,
    SNPs   = NULL,
    rho_ld = 0.95,
    rho_d  = 0.95,
    ld_th  = NULL,
    d_th   = NULL,
    col_vector = NULL,
    l_min  = 10,
    cores  = 1) {
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

  #ch = names(ld_decay$by_chr)[1]
  data.table::rbindlist(parallel::mclapply(names(ld_decay$by_chr), function(ch) {
    message(ch)
    map <- data.table::copy(map[Chr == ch])
    chr_obj <- ld_decay$by_chr[[ch]]
    if(length(chr_obj$el)==1){
      el <- data.table::fread(chr_obj$el, showProgress = FALSE)
    }else{
      el <- chr_obj$el
    }



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

    if(!is.null(SNPs)){
      ed <- el[
        r2 > ld_th_chr &
          d < d_th_chr &
          (SNP1 %in% SNPs | SNP2 %in% SNPs),
        .(SNP1, SNP2)
      ]
    }else{
      ed <- el[
        r2 > ld_th_chr &
          d < d_th_chr,
        (SNP1 %in% SNPs | SNP2 %in% SNPs),
        .(SNP1, SNP2)
      ]
    }


    map[, `:=`(CL_id = NA_character_, CL_col = "grey")]

    if (nrow(ed) == 0L) {
      return(map)
    }

    g <- igraph::graph_from_data_frame(ed, directed = FALSE)
    comps <- igraph::components(g)
    cls <- split(names(comps$membership), comps$membership)
    cls <- cls[vapply(cls, length, integer(1)) >= l_min]

    if (!length(cls)) {
      return(map)
    }

    cl_cols <- rep(col_vector, length.out = length(cls))

    cl_dt <- data.table::data.table(
      CL_id = rep(paste0(ch, "_", seq_along(cls)), lengths(cls)),
      marker = unlist(cls, use.names = FALSE),
      CL_col = rep(cl_cols, lengths(cls)),
      n_loci = rep(lengths(cls),lengths(cls))
    )

    map[cl_dt, `:=`(CL_id = i.CL_id, CL_col = i.CL_col,n_loci=i.n_loci), on = "marker"]

  },mc.cores=cores))
}

make_eMLGs <- function(GTs,
                       map_cl,
                       l_min = 0,
                       cores = 1) {




  map_cl <- map_cl[!is.na(CL_id) & n_loci >= l_min]

  cls <- split(map_cl$marker, map_cl$CL_id)
  cls <- cls[lengths(cls) >= max(2,l_min)]

  if (!length(cls)) {
    stop("No LD clusters passed l_min.")
  }

  # #input = "rSNP"
  # expected_fun <- switch(
  #   input,
  #   hard = expected_gt_hard,
  #   dosage = expected_gt_dosage,
  #   rSNP = select_representative_snp
  # )
  #select_representative_snp
  process_cluster <- function(markers) {
    select_representative_snp(GTs[, markers, drop = FALSE])
  }


  if(length(cls)>1000){
    batch_size <- 1000

    idx_chunks <- split(seq_along(cls),ceiling(seq_along(cls) / batch_size))

    n_batches <- length(idx_chunks)

    #t1 <- Sys.time()
    # i <- 1
    res <- lapply(seq_along(idx_chunks), function(i) {
      idx <- idx_chunks[[i]]

      message(sprintf(
        "Starting batch %d/%d: entries %d-%d",
        i, n_batches, min(idx), max(idx)
      ))



      #x <- GTs[, cls[idx][[1]], drop = FALSE]
      out <- parallel::mclapply(
        cls[idx],
        process_cluster,
        mc.cores = cores
      )
    })
    res <- unlist(res, recursive = FALSE)
  }else{
    message(sprintf(
      "Starting batch 1: entries %d-%d",1, length(cls)
    ))

    res <- parallel::mclapply(
      cls,
      process_cluster,
      mc.cores = cores
    )
  }

  # t2 <- Sys.time()



  #res$Chr1_1$
  # Representative SNP names
  #res$Chr1_1
  #which(lengths(rep_snps)==0)
  #table(lengths(res))
  #res[1:100]

  cls <- cls[sapply(res,function(x) length(x$snp))>0]
  res <- res[sapply(res,function(x) length(x$snp))>0]
  rep_snps <- vapply(res, `[[`, character(1), "snp")

  # Representative SNP column indices
  rep_indices <- vapply(res, `[[`, integer(1), "index")

  # Correlation of representative SNP with expected genotype
  R2_best <- vapply(res, `[[`, numeric(1), "R2_best")
  R2_E_hard <- vapply(res, `[[`, numeric(1), "R2_E_hard")

  # Expected genotype matrix
  # rows = individuals, columns = clusters
  eMLG <- do.call(cbind, lapply(res, `[[`, "eMLG"))
  colnames(eMLG) <- names(res)

  MLG <- do.call(cbind, lapply(res, `[[`, "MLG"))
  colnames(MLG) <- names(res)

  # Representative SNP genotype matrix (NA replaced by consensus)
  rep_gt_mat <- do.call(cbind, lapply(res, `[[`, "representative"))
  colnames(rep_gt_mat) <- rep_snps

  # Weight matrix: number of SNPs contributing per individual
  # rows = individuals, columns = clusters
  weight_mat <- do.call(cbind, lapply(res, `[[`, "n"))
  colnames(weight_mat) <- names(res)


  map_snp <- data.table::data.table(
    CL_id = rep(names(cls), lengths(cls)),
    marker = unlist(cls, use.names = FALSE)
  )

  map_snp <- add_snp_positions(map_snp, map_cl)

  map_eMLG <- data.table(CL_id = names(cls),
                         Chr = map_snp[match(rep_snps,marker),Chr],
                         n_loci = lengths(cls),
                         rep_snp = rep_snps,
                         R2_best = R2_best,
                         R2_E_hard = R2_E_hard)


  positions <- map_snp[,.(Chr=Chr[1],min_pos=min(Pos),max_pos=max(Pos)),by=CL_id]
  positions[,range_pos:=max_pos-min_pos]

  map_eMLG <- positions[map_eMLG,on="CL_id"]

  list(
    eMLG = eMLG,
    MLG = MLG,
    imp_rep_gt = rep_gt_mat,
    map_SNP = map_snp,
    map_eMLG = map_eMLG,
    weight_mat = weight_mat,
    clusters = cls
  )
}

#x <- GTs[, cls[idx][[1]], drop = FALSE]

expected_gt_hard <- function(x) {
  x <- polarize_genotypes(x)

  n <- rowSums(!is.na(x))
  c1 <- rowSums(x == 1, na.rm = TRUE)
  c2 <- rowSums(x == 2, na.rm = TRUE)

  E <- (c1 + 2 * c2) / n
  E[n == 0] <- NA_real_

  list(
    expected = E,
    n = n
  )
}

select_representative_snp <- function(x) {
  x <- polarize_genotypes(x)

  eg <- expected_gt_hard(x)
  E <- eg$expected
  n <- eg$n

  # Consensus genotype per individual
  consensus <- round(E)

  # Replace missing genotypes with consensus
  na_idx <- is.na(x)
  if (any(na_idx)) {
    row_ids <- row(x)[na_idx]
    x[na_idx] <- consensus[row_ids]
  }

  # Correlation to cluster-wide expected genotype
  cors <- apply(x, 2, function(g) {
    ok <- !is.na(E)
    if (sum(ok) < 3) return(NA_real_)
    cor(g[ok], E[ok])^2
  })

  best <- which.max(cors)

  # Representative SNP genotype with imputed missing values
  rep_gt <- x[, best]

  cluster_gt <- round(E)
  #cluster_gt[n < min_support] <- NA_integer_

  R2_E_hard <- cor(E, cluster_gt, use = "complete.obs")^2
  if(length(unique(na.omit(cluster_gt)))>1){
    list(
      snp = colnames(x)[best],
      index = unname(best),
      MLG = cluster_gt,
      R2_best = cors[best],
      R2_E_hard = R2_E_hard,
      representative = rep_gt,
      eMLG = E,
      n = n
    )
  }else{
    NULL
  }


}

refine_eMLGs_complete <- function(GTs,
                                  MLGs,
                                  input = "dosage",
                                  R2_th = 0.8,
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



  good_ids <- MLGs$map_eMLG[R2_E_hard >= R2_th,CL_id]
  bad_ids  <- MLGs$map_eMLG[R2_E_hard <  R2_th,CL_id]

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

    gr <- stats::cutree(hc, h = 1 - R2_th)

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

  split_bad

  expected_fun <- switch(
    input,
    hard = expected_gt_hard,
    dosage = expected_gt_dosage
  )

  process_cluster <- function(markers) {
    expected_gt_hard(GTs[, markers, drop = FALSE])
  }
  GTs[, markers, drop = FALSE]
  res <- parallel::mclapply(cls_final, process_cluster, mc.cores = cores)
  cor(eMLG[,"Chr11_827"],round(eMLG[,"Chr11_827"]))^2
  cor(res$Chr17_207.1$expected,round(res$Chr17_207.1$expected))^2

  eMLG_split <- do.call(cbind, res)
  colnames(eMLG_split) <- names(split_bad)
  rownames(eMLG_split) <- rownames(GTs)

  r2_final <- apply(eMLG_split, 2, score_eMLG)

  keep <- !is.na(r2_final) & r2_final >= R2_th

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
    R2_th = R2_th
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

  # Standardize position column names
  if ("min_pos" %in% names(map)) data.table::setnames(map, "min_pos", "Pos_min")
  if ("max_pos" %in% names(map)) data.table::setnames(map, "max_pos", "Pos_max")

  data.table::setorderv(map, c("Chr", "Pos_min"))

  map[, gap_bp := Pos_min - data.table::shift(Pos_max), by = Chr]
  map[, new_run := is.na(gap_bp) | gap_bp > distance_threshold, by = Chr]
  map[, run := cumsum(new_run), by = Chr]

  offdiag_stats <- function(r2mat) {
    if (ncol(r2mat) <= 1L) {
      return(list(mean = 1, min = 1))
    }

    vals <- r2mat[upper.tri(r2mat)]
    vals <- vals[is.finite(vals)]

    list(
      mean = mean(vals, na.rm = TRUE),
      min = min(vals, na.rm = TRUE)
    )
  }

  collapse_run <- function(gmap) {
    ids <- intersect(gmap$CL_id, colnames(eMLG))

    if (length(ids) == 0L) return(NULL)

    if (length(ids) == 1L) {
      y <- eMLG[, ids]

      return(list(
        C_eMLG = matrix(y, ncol = 1, dimnames = list(rownames(eMLG), ids)),
        map = data.table::data.table(
          old_C_eMLG_id = ids,
          n_eMLGs = 1L,
          n_loci = gmap[match(ids, CL_id), n_loci],
          mean_r2_to_members = 1,
          min_r2_to_members = 1,
          r2_eMLG = score_eMLG(y),
          members = list(ids)
        )
      ))
    }

    x <- eMLG[, ids, drop = FALSE]

    r2 <- suppressWarnings(stats::cor(x, use = "pairwise.complete.obs")^2)
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
        r2_members <- matrix(1, 1, 1)
      } else {
        w <- map[match(members, CL_id), n_loci]
        w[is.na(w)] <- 1

        E <- weighted_row_mean(
          polarize_genotypes(z),
          w
        )

        y <- round(E)
        y[is.na(E)] <- NA_real_

        r2_members <- suppressWarnings(
          stats::cor(z, use = "pairwise.complete.obs")^2
        )
        r2_members[!is.finite(r2_members)] <- NA_real_
        diag(r2_members) <- 1
      }

      od <- offdiag_stats(r2_members)

      list(
        y = y,
        n_loci = sum(map[match(members, CL_id), n_loci], na.rm = TRUE),
        mean_r2_to_members = od$mean,
        min_r2_to_members = od$min,
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


  chromosomes <- unique(map$Chr)
  n_chr <- length(chromosomes)

  res <- list()

  for (i in seq_along(chromosomes)) {
    chr <- chromosomes[i]

    message(sprintf(
      "Starting chromosome %s (%d/%d)",
      chr, i, n_chr
    ))

    chr_runs <- runs[vapply(runs, function(z) z$Chr[1] == chr, logical(1))]

    res_chr <- lapply(chr_runs, collapse_run)

    res <- c(res, res_chr)
  }

  res <- Filter(Negate(is.null), res)

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

  members <- groups[[1]]
  MLGs$map_eMLG[match(members, C_eMLG_id)]

  score_group <- function(membrs) {
    x <- eMLG[, membrs, drop = FALSE]

    n_loci_group <- sum(
      MLGs$map_eMLG[
        match(membrs, C_eMLG_id),
        n_loci
      ],
      na.rm = TRUE
    )

    if (length(membrs) == 1L) {
      y <- as.numeric(x[, 1])
      r2_members <- 1
    } else {

      w <- MLGs$map_eMLG[
        match(membrs, C_eMLG_id),
        n_loci
      ]

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
      n_loci = n_loci_group,
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
    n_loci = vapply(scores, `[[`, numeric(1), "n_loci"),
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
    clusters = split(lookup$C_eMLG_id, lookup$G_eMLG_id),
    params = list(
      stage = "collapse_eMLGs_global",
      r2_threshold = r2_threshold,
      method = method
    )
  )
}

annotate_C_pair_distance <- function(map_SNP) {
  C_map <- map_SNP[
    ,
    .(
      Chr = Chr[1],
      Pos_min = as.numeric(min(Pos, na.rm = TRUE)),
      Pos_max = as.numeric(max(Pos, na.rm = TRUE)),
      Pos_mid = as.numeric(stats::median(Pos, na.rm = TRUE))
    ),
    by = C_eMLG_id
  ]

  pairs <- data.table::CJ(
    C1 = C_map$C_eMLG_id,
    C2 = C_map$C_eMLG_id
  )[C1 < C2]

  pairs[C_map, `:=`(
    Chr1 = i.Chr,
    Pos_min1 = i.Pos_min,
    Pos_max1 = i.Pos_max
  ), on = .(C1 = C_eMLG_id)]

  pairs[C_map, `:=`(
    Chr2 = i.Chr,
    Pos_min2 = i.Pos_min,
    Pos_max2 = i.Pos_max
  ), on = .(C2 = C_eMLG_id)]

  pairs[
    ,
    same_chr := Chr1 == Chr2
  ]

  pairs[
    ,
    gap_bp := data.table::fifelse(
      same_chr,
      pmax(0, pmax(Pos_min1, Pos_min2) - pmin(Pos_max1, Pos_max2)),
      Inf
    )
  ]

  pairs[
    ,
    overlap := same_chr & gap_bp == 0
  ]

  pairs[]
}
