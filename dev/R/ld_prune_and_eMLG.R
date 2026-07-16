# Combined LD-pruning + eMLG generation, replacing merge_ld_clusters() for
# clusters flagged by high local LD support (ld_w).
#
# Context: the existing pipeline (R/ld_pruning_DIEM.R, formerly baypass.R)
# runs ld_complexity_reduction() once genome-wide, flags clusters containing
# a marker above an ld_w
# threshold, then runs merge_ld_clusters() (complete linkage @ rho=0.5) on
# just the flagged clusters to reunite sliding-window-fragmented blocks --
# and SEPARATELY builds eMLGs from that output and dynamic-cuts them AGAIN
# (average linkage, score_eMLG-gated) for the Ohta-statistic/long-range-LD
# work. That's two independent all-pairs correlation passes over
# overlapping/related data.
#
# This does ONE pass instead: Stage 1 -> flag -> distance-restricted
# dynamic cut directly on the flagged Stage-1 clusters, producing BOTH a
# pruned representative marker set (for LD-pruning/GRM) AND an eMLG matrix
# (for downstream long-range LD analysis) from the same clustering.
#
# Why this is defensible for LD-pruning too, not just eMLG summarization
# (raised and worked through in conversation, not just assumed):
#   - complete linkage's strict "every pairwise r2 must clear the
#     threshold" guarantee exists to avoid pruning away real, independent
#     genetic signal -- appropriate when markers COULD be genuinely
#     unlinked. But in a young hybrid population's low-recombination
#     region, that concern doesn't really apply: with very few actual
#     recombination events having occurred there, the region genuinely
#     behaves as a small number of effectively-independent haplotype
#     blocks, not many independent markers being falsely conflated.
#     Representing it with fewer units isn't discarding signal, it's
#     accurately reflecting the true effective degrees of freedom.
#   - the real risk was merging PHYSICALLY DISTANT clusters whose
#     correlation might reflect intra-chromosomal epistatic selection
#     rather than physical linkage -- exactly the signal later analysis
#     (Ohta statistics) wants to detect, not erase by premature merging.
#     Fixed by restricting merging to physically-contiguous runs
#     (single-linkage BY DISTANCE: consecutive cluster gap <=
#     distance_threshold, default 5e5bp). Note this deliberately does NOT
#     cap a run's total physical span -- a genuine non-recombining block
#     is expected to be one long, closely-spaced, contiguous run that can
#     legitimately span megabases; what should break it into separate runs
#     is a real gap between neighbours, not the block's overall size.
#     (Tried capping total span first; wrong, reverted -- see
#     distance_threshold_semantics in project memory if this gets
#     revisited.)
#
# See dynamic_cut_eMLG.R for the average-vs-single-vs-complete-linkage
# comparison and the score_eMLG/pair_r2 gating rationale/bugs already
# fixed; this file reuses that function unchanged as the within-run
# merging step and adds: flagging, distance-restricted run splitting
# across chromosomes, unflagged cluster passthrough, and representative
# marker selection.
#
# Cost: the dominant cost (all-pairs cor()) scales quadratically with
# cluster count (confirmed directly: ~0.01s at 292 clusters, ~31s at
# 15,000) -- restricting to ld_w-flagged clusters is load-bearing for
# keeping this tractable, not a workaround to relax. At full-chromosome,
# unfiltered scale this would hit the same wall merge_ld_clusters() did
# before min_n_snps filtering.
#
# eMLG computation is itself opt-out/opt-in at two levels, since callers
# doing LD-pruning only don't need eMLGs at all, and long-range-LD analysis
# will apply its own n_loci threshold downstream anyway (too few loci isn't
# informative for detecting extended LD patterns regardless):
#   - compute_unflagged_eMLG: skip eMLG computation for unflagged clusters
#     entirely (usually the large majority of clusters) when only the
#     pruned marker set is needed.
#   - min_n_loci_eMLG: skip storing an eMLG for any group (flagged or
#     unflagged) below this many raw loci -- the pruned representative
#     marker is unaffected either way, this only trims the returned eMLG
#     matrix.

source("./dev/R/dynamic_cut_eMLG.R")

## expected_gt_dosage() is provably the identity for a single-marker
## cluster (polarize_genotypes() short-circuits at ncol<=1, rowMeans() over
## one column returns that column unchanged) -- skip the redundant
## function-call/polarization overhead and use the raw column directly.
## Always correct, no output change, so applied unconditionally rather than
## gated by a parameter.
consensus_dosage <- function(GTs, mk) {
  if (length(mk) == 1) return(GTs[, mk])
  expected_gt_dosage(GTs[, mk, drop = FALSE])
}

## Splits clusters into physically-contiguous runs: single linkage BY
## DISTANCE, i.e. a new run starts whenever the gap to the previous
## (position-sorted) cluster exceeds distance_threshold. A run's TOTAL span
## is deliberately unbounded -- see file header.
split_by_distance <- function(pos_min, pos_max, distance_threshold = 5e5) {
  ord <- order(pos_min)
  gap <- pos_min[ord] - c(NA, pos_max[ord][-length(ord)])
  new_run <- is.na(gap) | gap > distance_threshold
  run <- cumsum(new_run)
  run_id <- integer(length(ord))
  run_id[ord] <- run
  names(run_id) <- names(pos_min)[ord]
  run_id[names(pos_min)]
}

## Representative marker for a final group: among its constituent Stage-1
## clusters, the one whose consensus signal has the highest median r2 to
## the OTHERS in the group (same "highest median r2" principle
## ld_complexity_reduction()/merge_ld_clusters() use), using that cluster's
## own core_snp. Cheap even though computed per-group: groups are small
## (tens of constituent clusters at most, not raw markers).
pick_representative <- function(cl_ids, eMLG, stage1_clusters) {
  if (length(cl_ids) == 1) {
    return(stage1_clusters[CL_id == cl_ids, core_snp])
  }
  sub_r2 <- suppressWarnings(stats::cor(eMLG[, cl_ids, drop = FALSE], use = "pairwise.complete.obs")^2)
  sub_r2[!is.finite(sub_r2)] <- 0
  diag(sub_r2) <- NA
  med_r2 <- apply(sub_r2, 1, stats::median, na.rm = TRUE)
  best_cl <- cl_ids[which.max(med_r2)]
  stage1_clusters[CL_id == best_cl, core_snp]
}

#' Combined LD-pruning + eMLG generation via distance-restricted dynamic cut
#'
#' @param GTs Numeric matrix, individuals x markers.
#' @param stage1 Output of [ld_complexity_reduction()].
#' @param ld_w_col Name of the ld_w column in `stage1$map_snp` used to flag
#'   which clusters need the expensive treatment (a cluster is flagged if
#'   ANY member exceeds `ld_w_threshold`).
#' @param ld_w_threshold Flagging threshold.
#' @param score_threshold Minimum `score_eMLG` a candidate merge's result
#'   must clear (see [dynamic_cut_eMLG()]).
#' @param min_r2 Minimum r2 required directly between the two sides being
#'   merged (see [dynamic_cut_eMLG()]).
#' @param distance_threshold Max consecutive-gap in bp allowed within one
#'   physically-contiguous run (see [split_by_distance()]). Clusters more
#'   than this apart are never merged, regardless of correlation.
#' @param compute_unflagged_eMLG If `FALSE`, skip eMLG computation for
#'   unflagged clusters entirely (usually the large majority) -- for
#'   callers who only need the pruned marker set, not eMLGs. Default `TRUE`.
#' @param min_n_loci_eMLG Skip storing an eMLG for any group (flagged or
#'   unflagged) with fewer than this many raw loci -- the pruned
#'   representative marker is returned regardless, this only trims the
#'   returned `eMLG` matrix. Default `1` (no filtering).
#'
#' @return A list: `eMLG` (matrix, individuals x groups that have one --
#'   see `compute_unflagged_eMLG`/`min_n_loci_eMLG`), `groups` (data.table:
#'   group_id, Chr, representative, n_loci, score, has_eMLG, members --
#'   always lists EVERY group regardless of eMLG filtering), `pruned`
#'   (character vector of representative markers, one per group -- for
#'   LD-pruning use, unaffected by eMLG filtering).
ld_prune_and_eMLG <- function(GTs, stage1, ld_w_col, ld_w_threshold,
                               score_threshold = 0.80, min_r2 = 0.2,
                               distance_threshold = 5e5,
                               compute_unflagged_eMLG = TRUE,
                               min_n_loci_eMLG = 1) {

  map_snp  <- stage1$map_snp
  clusters <- stage1$clusters

  needs_merge_ids <- map_snp[get(ld_w_col) > ld_w_threshold, unique(CL_id)]
  flagged   <- clusters[CL_id %in% needs_merge_ids]
  unflagged <- clusters[!CL_id %in% needs_merge_ids]

  ## unflagged clusters pass straight through unchanged -- their own
  ## core_snp is already the representative regardless of whether an eMLG
  ## is computed for them at all.
  if (nrow(unflagged)) {
    eligible <- compute_unflagged_eMLG & unflagged$n_snps >= min_n_loci_eMLG
    unflagged_for_eMLG <- unflagged[eligible]

    unflagged_eMLG <- if (nrow(unflagged_for_eMLG)) {
      n_u <- nrow(unflagged_for_eMLG)
      message("Computing eMLGs for ", n_u, " unflagged clusters...")
      pb <- utils::txtProgressBar(min = 0, max = n_u, style = 3)
      m_list <- vector("list", n_u)
      for (i in seq_len(n_u)) {
        m_list[[i]] <- consensus_dosage(GTs, unflagged_for_eMLG$members[[i]])
        utils::setTxtProgressBar(pb, i)
      }
      close(pb)
      m <- do.call(cbind, m_list)
      colnames(m) <- paste0("U", unflagged_for_eMLG$CL_id)
      m
    } else {
      matrix(numeric(0), nrow = nrow(GTs), ncol = 0)
    }

    scores <- setNames(rep(NA_real_, nrow(unflagged)), unflagged$CL_id)
    if (nrow(unflagged_for_eMLG)) {
      scores[as.character(unflagged_for_eMLG$CL_id)] <- apply(unflagged_eMLG, 2, score_eMLG)
    }

    unflagged_groups <- data.table::data.table(
      group_id = paste0("U", unflagged$CL_id),
      Chr = unflagged$Chr,
      representative = unflagged$core_snp,
      n_loci = unflagged$n_snps,
      score = scores[as.character(unflagged$CL_id)],
      has_eMLG = eligible,
      members = unflagged$members
    )
  } else {
    unflagged_eMLG <- matrix(numeric(0), nrow = nrow(GTs), ncol = 0)
    unflagged_groups <- data.table::data.table(
      group_id = character(0), Chr = character(0), representative = character(0),
      n_loci = integer(0), score = numeric(0), has_eMLG = logical(0), members = list()
    )
  }

  if (nrow(flagged) == 0) {
    return(list(
      eMLG = unflagged_eMLG,
      groups = unflagged_groups,
      pruned = unflagged_groups$representative
    ))
  }

  ## flagged clusters' own eMLGs are needed as input to the dynamic cut
  ## regardless of size or min_n_loci_eMLG -- even a single-marker flagged
  ## cluster might merge into something larger, so this can't be skipped
  ## the way unflagged eMLGs can be. consensus_dosage() still applies the
  ## n_loci==1 fast path for the ones that stay singletons.
  n_f <- nrow(flagged)
  message("Computing eMLGs for ", n_f, " flagged clusters...")
  pb <- utils::txtProgressBar(min = 0, max = n_f, style = 3)
  flagged_eMLG_list <- vector("list", n_f)
  for (i in seq_len(n_f)) {
    flagged_eMLG_list[[i]] <- consensus_dosage(GTs, flagged$members[[i]])
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  flagged_eMLG <- do.call(cbind, flagged_eMLG_list)
  colnames(flagged_eMLG) <- flagged$CL_id

  n_loci_flagged <- setNames(flagged$n_snps, flagged$CL_id)

  pos_min <- setNames(
    vapply(flagged$members, function(mk) min(map_snp[marker %in% mk, Pos]), numeric(1)),
    flagged$CL_id
  )
  pos_max <- setNames(
    vapply(flagged$members, function(mk) max(map_snp[marker %in% mk, Pos]), numeric(1)),
    flagged$CL_id
  )

  flagged_group_list <- list()

  ## restricted to same chromosome AND same physically-contiguous run --
  ## never merge across chromosomes (no physical basis) or across a real
  ## position gap (possible epistatic-selection signal, not redundancy).
  ## Progress is tracked per chromosome, not per run: run count per
  ## chromosome is a poor proxy for time remaining (cost is dominated by
  ## the largest run's O(n^2) cor(), and a chromosome can have very few but
  ## very large runs -- Chr26's ld_w>0.2 subset was 4 runs from 292
  ## clusters, one of them 263), so a per-chromosome message reports run
  ## count/sizes for context instead of trying to bar-track them.
  chr_levels <- unique(flagged$Chr)
  message("Distance-restricted dynamic cut: ", length(chr_levels), " chromosomes")
  pb <- utils::txtProgressBar(min = 0, max = length(chr_levels), style = 3)

  for (chr_i in seq_along(chr_levels)) {
    ch <- chr_levels[chr_i]
    ids_ch <- as.character(flagged$CL_id[flagged$Chr == ch])
    runs <- split_by_distance(pos_min[ids_ch], pos_max[ids_ch], distance_threshold)
    run_sizes <- table(runs)
    message(
      ch, ": ", length(ids_ch), " flagged clusters -> ", length(run_sizes), " run(s) ",
      "(largest ", max(run_sizes), ")"
    )

    for (r in unique(runs)) {
      run_ids <- names(runs)[runs == r]

      if (length(run_ids) == 1) {
        cid <- run_ids
        flagged_group_list[[length(flagged_group_list) + 1]] <- list(
          Chr = ch, cl_ids = cid, emlg = flagged_eMLG[, cid],
          n_loci = n_loci_flagged[[cid]]
        )
        next
      }

      sub_groups <- dynamic_cut_eMLG(
        flagged_eMLG[, run_ids, drop = FALSE], n_loci_flagged,
        threshold = score_threshold, min_r2 = min_r2
      )
      for (g in sub_groups) {
        flagged_group_list[[length(flagged_group_list) + 1]] <- list(
          Chr = ch, cl_ids = g$members, emlg = g$emlg,
          n_loci = sum(n_loci_flagged[g$members])
        )
      }
    }
    utils::setTxtProgressBar(pb, chr_i)
  }
  close(pb)

  flagged_groups <- data.table::rbindlist(lapply(seq_along(flagged_group_list), function(i) {
    g <- flagged_group_list[[i]]
    data.table::data.table(
      group_id = paste0("F", i),
      Chr = g$Chr,
      representative = pick_representative(g$cl_ids, flagged_eMLG, flagged),
      n_loci = g$n_loci,
      score = score_eMLG(g$emlg),
      has_eMLG = g$n_loci >= min_n_loci_eMLG,
      members = list(unlist(flagged$members[flagged$CL_id %in% g$cl_ids]))
    )
  }))

  ## the merged eMLG signal is already computed as a byproduct of the
  ## dynamic cut regardless of size, so min_n_loci_eMLG only trims which
  ## columns are kept in the returned matrix here, not what gets computed
  keep <- flagged_groups$has_eMLG
  flagged_eMLG_final <- do.call(cbind, lapply(flagged_group_list[keep], `[[`, "emlg"))
  colnames(flagged_eMLG_final) <- flagged_groups$group_id[keep]

  groups <- data.table::rbindlist(list(unflagged_groups, flagged_groups))

  list(
    eMLG = cbind(unflagged_eMLG, flagged_eMLG_final),
    groups = groups,
    pruned = groups$representative
  )
}
