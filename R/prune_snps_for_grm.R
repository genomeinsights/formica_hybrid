#' Prune markers for GRM construction using LD-network components
#'
#' Constructs a pruned set of representative SNPs for genomic relationship matrix
#' (GRM) estimation by grouping markers into LD-connected components and retaining
#' one representative SNP per component.
#'
#' For each chromosome, an LD threshold is derived from the fitted LD-decay model
#' at a user-specified decay level (`rho_ld`). SNP pairs with LD above this
#' threshold are treated as edges in an undirected graph. Connected components of
#' this graph define LD blocks, and one representative SNP is retained per block.
#' The representative SNP is chosen as the marker with the highest local LD
#' support (`ld_w_col`), which favors central markers within each LD block.
#'
#' This pruning strategy is intended for GRM construction, where the goal is to
#' reduce over-representation of extended haplotypes and long-range LD regions
#' while preserving genome-wide background relatedness structure.
#'
#' @param map A `data.table` containing marker metadata. Must include at least
#'   columns `Chr`, `marker`, and the column specified by `ld_w_col`.
#' @param ld_decay A list-like object containing chromosome-specific LD-decay
#'   results. Must contain:
#'   \describe{
#'     \item{`by_chr`}{A named list with one entry per chromosome. Each entry must
#'     contain an element `el`, giving the path to or a `data.table` of an edge
#'     list with LD values.}
#'     \item{`decay_sum`}{A table containing one row per chromosome and columns
#'     `Chr`, `b`, and `c_pred`, used by [ld_from_rho()] to derive the LD
#'     threshold.}
#'   }
#' @param rho_ld Numeric scalar in `(0, 1)`. Decay level used to derive the LD
#'   threshold for pruning. For example, `rho_ld = 0.5` uses the LD threshold
#'   corresponding to the point where the fitted decay curve reaches 50% of its
#'   scaled decay.
#' @param ld_w_col Character scalar giving the name of the column in `map` used
#'   to rank markers within each LD block. The marker with the highest value is
#'   retained as the representative SNP. Default is `"ld_w_95"`.
#' @param edge_ld_col Character scalar giving the name of the LD column in the
#'   edge list. Default is `"r2"`.
#' @param pos1_col,pos2_col Character scalars giving the names of the position
#'   columns in the edge list. Defaults are `"pos1"` and `"pos2"`.
#' @param marker_sep Character scalar used to construct marker IDs from chromosome
#'   and position. Default is `":"`.
#' @param block_prefix Character scalar used as prefix for generated block IDs.
#'   Default is `"GRM"`.
#' @param show_progress Logical; passed to [data.table::fread()] when edge lists
#'   are stored as file paths. Default is `FALSE`.
#'
#' @return A `data.table` with one row per retained LD block representative and
#'   the following columns:
#'   \describe{
#'     \item{`Chr`}{Chromosome identifier.}
#'     \item{`block_id`}{LD block identifier.}
#'     \item{`core_snp`}{Representative SNP retained for the block.}
#'     \item{`n_snps`}{Number of SNPs in the block.}
#'     \item{`mean_ld_w`}{Mean value of `ld_w_col` across SNPs in the block.}
#'   }
#'
#' @details
#' Isolated SNPs with no LD edges above the threshold are retained as singleton
#' blocks. All markers in `map` are therefore assigned to exactly one block.
#'
#' This function is designed for pruning markers prior to GRM estimation, not for
#' defining outlier regions in the main scan. In practice, it may be useful to
#' compare GRMs generated with and without known inversion regions as a
#' sensitivity analysis.
#'
#' @seealso [ld_from_rho()]
#'
#' @examples
#' \dontrun{
#' pruned_snps <- prune_snps_for_grm(
#'   map = map,
#'   ld_decay = ld_decay,
#'   rho_ld = 0.5
#' )
#' }
#'
#' @export
prune_snps_for_grm <- function(map,
                               ld_decay,
                               rho_ld = 0.5,
                               ld_w_col = "ld_w_95",
                               edge_ld_col = "r2",
                               pos1_col = "pos1",
                               pos2_col = "pos2",
                               marker_sep = ":",
                               block_prefix = NULL,
                               show_progress = FALSE) {
  if (!data.table::is.data.table(map)) {
    map <- data.table::as.data.table(map)
  }

  req_map_cols <- c("Chr", "marker", ld_w_col)
  miss_map <- setdiff(req_map_cols, names(map))
  if (length(miss_map) > 0) {
    stop("`map` is missing required column(s): ",
         paste(miss_map, collapse = ", "))
  }

  if (!is.list(ld_decay) || is.null(ld_decay$by_chr) || is.null(ld_decay$decay_sum)) {
    stop("`ld_decay` must contain components `by_chr` and `decay_sum`.")
  }

  if (!data.table::is.data.table(ld_decay$decay_sum)) {
    decay_sum <- data.table::as.data.table(ld_decay$decay_sum)
  } else {
    decay_sum <- data.table::copy(ld_decay$decay_sum)
  }

  req_decay_cols <- c("Chr", "b", "c_pred")
  miss_decay <- setdiff(req_decay_cols, names(decay_sum))
  if (length(miss_decay) > 0) {
    stop("`ld_decay$decay_sum` is missing required column(s): ",
         paste(miss_decay, collapse = ", "))
  }

  if (!is.numeric(rho_ld) || length(rho_ld) != 1L || is.na(rho_ld) ||
      rho_ld <= 0 || rho_ld >= 1) {
    stop("`rho_ld` must be a single numeric value in (0, 1).")
  }

  chrs <- intersect(names(ld_decay$by_chr), unique(map$Chr))
  if (length(chrs) == 0L) {
    stop("No overlapping chromosomes found between `map$Chr` and `names(ld_decay$by_chr)`.")
  }

  pb <- txtProgressBar(min = 0, max = length(chrs)-1, style = 3)
  setTxtProgressBar(pb, 0)

  out <- lapply(chrs, function(ch) {
    mp <- data.table::copy(map[Chr == ch])

    if (nrow(mp) == 0L) return(NULL)

    pars <- decay_sum[Chr == ch]
    if (nrow(pars) != 1L) {
      stop("Expected exactly one row in `ld_decay$decay_sum` for chromosome ", ch, ".")
    }

    ld_th <- ld_from_rho(
      b = pars$b,
      c = pars$c_pred,
      rho = rho_ld
    )

    chr_obj <- ld_decay$by_chr[[ch]]
    if (is.null(chr_obj$el)) {
      stop("`ld_decay$by_chr[[", ch, "]]` must contain element `el`.")
    }

    el <- chr_obj$el
    if (is.character(el) && length(el) == 1L) {
      el <- data.table::fread(el, showProgress = show_progress)
    } else if (!data.table::is.data.table(el)) {
      el <- data.table::as.data.table(el)
    }

    req_el_cols <- c(edge_ld_col, pos1_col, pos2_col)
    miss_el <- setdiff(req_el_cols, names(el))
    if (length(miss_el) > 0) {
      stop("Edge list for chromosome ", ch, " is missing required column(s): ",
           paste(miss_el, collapse = ", "))
    }

    ed <- el[get(edge_ld_col) > ld_th, .(
      SNP1 = paste(ch, get(pos1_col), sep = marker_sep),
      SNP2 = paste(ch, get(pos2_col), sep = marker_sep)
    )]

    g <- igraph::graph_from_data_frame(
      d = ed,
      directed = FALSE,
      vertices = data.frame(name = mp$marker, stringsAsFactors = FALSE)
    )

    comps <- igraph::components(g)$membership

    mp[, block_id := paste0(
      block_prefix, "_", ch, "_", comps[match(marker, names(comps))]
    )]

    setTxtProgressBar(pb, which(names(ld_decay$by_chr)==chr_obj$decay_sum$Chr))

    mp[, .(
      Chr = ch,
      core_snp = marker[which.max(get(ld_w_col))],
      n_snps = .N,
      mean_ld_w = mean(get(ld_w_col), na.rm = TRUE)
    ), by = block_id]

  })
  close(pb)


  out <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  print(out[, .(
    n_blocks = .N,
    mean_size = mean(n_snps),
    max_size = max(n_snps)
  )])

  return(out)
}
