## =========================================================
## Cluster-level (eMLG) parallelism -- consensus assembly
## =========================================================
##
## Companion to parallelism_stats(): running the binomial parallelism test
## directly on SNPs treats every marker as an independent trial, but markers
## within the same Stage 2 LD-cluster (see ld_prune_and_eMLG(), LDscnR) are
## correlated -- a single real (or noisy) sorting pattern in one block can get
## counted as dozens or hundreds of "independent" significant SNPs. The fix is
## to test each CLUSTER once, on its consensus dosage, not each member SNP.
##
## APPROACH (one canonical consensus definition, no re-clustering):
##   A1  run parallelism_stats() per SNP -> the "sorted" markers.
##   A2  build_sorted_eMLG(): assemble a hybrid-side consensus for every
##       Stage-2 cluster those sorted markers touch, REUSING the columns
##       already stored in the frozen clustering (clusters that cleared
##       min_n_loci_eMLG) and BUILDING only the rest fresh, with the identical
##       method (consensus_dosage / build_group_consensus). A cluster below the
##       size floor -- including a singleton -- is just "built" (a 1-marker
##       consensus is that marker's own dosage), so there is no
##       singleton/small/large branching: every unit is "reused" or "built".
##   A3  add a matched parent side (consensus_dosage_matched) and run
##       parallelism_stats() ONCE on the combined matrix -> one BH-FDR over all
##       units, natively coherent.
##
## The canonical eMLG_5loci_0025_cM05.rds stays FROZEN (it is geared to the
## forthcoming Ohta analyses and is the keystone Modules B/C/D join on). A2
## writes a SEPARATE, sorting-derived companion file, so the canonical
## clustering never has to be regenerated to test sorting.
##
## Why the parent side is built separately (not stored in either file): the
## clustering / eMLG is built from HYBRIDS ONLY by design (clustering with
## parents in would let parent-vs-hybrid structure distort every LD estimate).
## But parallelism_stats() needs the parents' allele frequency per unit to
## orient ancestry. consensus_dosage_matched() applies the hybrid-derived
## reference/flip decisions to the parent matrix, so the parent-side consensus
## sits on the SAME orientation as the stored hybrid-side eMLG -- the identical
## flip basis LDscnR used to build it (ref = hybrids).

library(data.table)
library(parallel)


## ---------------------------------------------------------
## consensus_dosage_matched()
##
## The consensus dosage for a set of markers `mk`, computed on GTs_new, but
## using the reference-marker choice and per-marker flip decisions derived
## from GTs_ref instead of recomputing them fresh on GTs_new.
##
## Why this matters: LDscnR's own consensus_dosage()/polarize_genotypes()
## pick a reference marker (most complete column) and flip any other column
## anti-correlated with it -- computed on whatever matrix they're given.
## Calling that fresh on a parent-only subset (30 individuals, strong
## between-species differentiation) can pick a different reference and flip
## different columns than what actually built the existing hybrid-side eMLG
## (164 individuals) -- silently putting the two consensus vectors on
## OPPOSITE orientations for any marker where the two choices disagree, and
## making them impossible to compare or combine.
##
## consensus_dosage_matched(GTs_ref, GTs_ref, mk) reproduces LDscnR's own
## consensus_dosage(GTs_ref, mk) exactly (same reference/flip logic, same
## cor_threshold=0 default) -- this function is a strict generalisation,
## not a different method.
##
## GTs_ref, GTs_new : individuals x markers dosage matrices (can be the
##                     same matrix, or e.g. hybrids vs. parents).
## mk                : character vector of member marker names.
## ---------------------------------------------------------
consensus_dosage_matched <- function(GTs_ref, GTs_new, mk) {
  if (length(mk) == 1) return(GTs_new[, mk])

  x_ref <- as.matrix(GTs_ref[, mk, drop = FALSE])
  ref_idx <- which.max(colSums(!is.na(x_ref)))
  ref <- x_ref[, ref_idx]
  flip <- vapply(seq_len(ncol(x_ref)), function(j) {
    r <- suppressWarnings(stats::cor(ref, x_ref[, j], use = "pairwise.complete.obs"))
    is.finite(r) && r < 0
  }, logical(1))

  x_new <- as.matrix(GTs_new[, mk, drop = FALSE])
  if (any(flip)) x_new[, flip] <- 2 - x_new[, flip]

  E <- rowMeans(x_new, na.rm = TRUE)
  E[is.nan(E)] <- NA_real_
  E
}


## ---------------------------------------------------------
## build_group_consensus()
##
## mclapply wrapper around consensus_dosage_matched() for many groups at
## once. Returns an individuals(GTs_new) x groups matrix.
##
## members_list   : named list, one character vector of marker names per
##                  group (names become the output's column names).
## progress, label : if progress=TRUE, print a per-chunk progress line with
##                  elapsed time + ETA under `label` (dependency-free; ~50
##                  updates max). Off by default.
## ---------------------------------------------------------
build_group_consensus <- function(members_list, GTs_ref, GTs_new, cores = 1,
                                  progress = FALSE, label = "consensus") {
  n <- length(members_list)
  worker <- function(mk) consensus_dosage_matched(GTs_ref, GTs_new, mk)

  if (!progress || n <= 1L) {
    res <- mclapply(members_list, worker, mc.cores = cores)
  } else {
    ## chunk so progress + ETA can be reported between chunks; each chunk is
    ## one mclapply over `cores` workers. Capped at ~50 updates.
    n_chunks <- min(n, 50L)
    idx <- split(seq_len(n), cut(seq_len(n), n_chunks, labels = FALSE))
    res <- vector("list", length(idx))
    t0 <- Sys.time(); done <- 0L
    for (ci in seq_along(idx)) {
      res[[ci]] <- mclapply(members_list[idx[[ci]]], worker, mc.cores = cores)
      done <- done + length(idx[[ci]])
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      message(sprintf("  [%s] %d/%d units (%.0f%%) | %.0fs elapsed | ~%.0fs left",
                      label, done, n, 100 * done / n, el, el * (n - done) / done))
    }
    res <- unlist(res, recursive = FALSE)
  }

  out <- do.call(cbind, res)
  rownames(out) <- rownames(GTs_new)
  colnames(out) <- names(members_list)
  out
}


## ---------------------------------------------------------
## cluster_DI()
##
## Aggregate a marker-level DiagnosticIndex to cluster level, for the min_DI
## gate at the cluster (eMLG) resolution.
##
## groups : ld_prune_and_eMLG() result$groups (needs group_id, representative,
##          members). Joined with on="group_id", so no key is required.
## ids    : group_ids to aggregate; the result is reindexed to this order,
##          with NA for any id not found.
## DI     : marker-named DiagnosticIndex vector.
## di_agg : "max" (default, pipeline-locked) uses the best (max) DI across ALL
##          of a cluster's members, so a cluster is only excluded if NONE of
##          its members were individually diagnostic enough. "representative"
##          uses only the cluster's representative marker's DI -- cheaper, but
##          can gate out a cluster even when the specific SNP that triggered
##          interest in it individually passed DI.
## ---------------------------------------------------------
cluster_DI <- function(groups, ids, DI, di_agg = c("max", "representative")) {
  di_agg <- match.arg(di_agg)
  if (di_agg == "representative") {
    rep_mk <- groups[.(ids), on = "group_id", representative]
    return(setNames(DI[rep_mk], ids))
  }
  ## "max": explode groups to one row per member, look each member's DI up by
  ## name, take the grouped max. all-NA group -> max(numeric(0)) = -Inf -> NA;
  ## every requested id is returned (reindexed by `ids`).
  g <- groups[.(ids), on = "group_id", .(group_id, members)]
  long <- data.table(
    group_id = rep(g$group_id, lengths(g$members)),
    di       = DI[unlist(g$members)]
  )
  agg <- long[, .(di = suppressWarnings(max(di, na.rm = TRUE))), by = group_id]
  agg[!is.finite(di), di := NA_real_]
  setNames(agg$di, agg$group_id)[ids]
}


## ---------------------------------------------------------
## build_sorted_eMLG()
##
## Extend the hybrid-side eMLG consensus to EVERY Stage-2 cluster touched by
## `markers` (the sorted SNPs), reusing the existing consensus columns from a
## frozen clustering where they exist and building the rest fresh with the
## identical method. No re-clustering: the partition (groups) is taken as-is
## from the frozen clustering. There is no singleton/small/large branching --
## a unit is either "reused" (already in `eMLG`) or "built".
##
## groups      : frozen clustering $groups (group_id, representative, n_loci,
##               members). Joined with on="group_id" (no key required).
## eMLG        : frozen clustering $eMLG (hybrid individuals x group_id,
##               >= min_n_loci_eMLG columns only).
## GTs_hybrids : the hybrid dosage matrix the eMLG was built from (rows define
##               the canonical individual order used for the output).
## markers     : sorted SNPs; every cluster containing >= 1 of them is covered.
## cores       : passed to build_group_consensus() for the freshly-built units.
## progress    : if TRUE, report the reuse/build split and per-chunk build
##               progress (passed through to build_group_consensus()).
##
## Returns list(eMLG = <hybrid individuals x unit consensus, reused+built>,
##              groups = <the covered subset of `groups`>,
##              source = <named "reused"/"built" per unit>). Same $eMLG/$groups
## shape as the canonical clustering, minus a parent side (built at test time).
## ---------------------------------------------------------
build_sorted_eMLG <- function(groups, eMLG, GTs_hybrids, markers, cores = 1,
                              progress = FALSE) {
  ind <- rownames(GTs_hybrids)                 # canonical individual order

  marker_to_group <- groups[, .(marker = unlist(members)), by = group_id]
  touched_ids <- unique(
    marker_to_group[.(markers), on = "marker", nomatch = NULL, group_id]
  )

  reuse_ids <- intersect(touched_ids, colnames(eMLG))
  build_ids <- setdiff(touched_ids, reuse_ids)
  if (progress) message(sprintf(
    "  [A2] %d clusters touched by sorted markers: %d reused (frozen) + %d to build",
    length(touched_ids), length(reuse_ids), length(build_ids)))

  parts <- list()
  if (length(reuse_ids)) parts$reuse <- eMLG[ind, reuse_ids, drop = FALSE]
  if (length(build_ids)) {
    bmem <- setNames(groups[.(build_ids), on = "group_id", members], build_ids)
    bcons <- build_group_consensus(bmem, GTs_hybrids, GTs_hybrids, cores = cores,
                                   progress = progress, label = "A2 build")
    parts$build <- bcons[ind, , drop = FALSE]
  }
  eMLG_sorted <- do.call(cbind, parts)

  src <- c(
    setNames(rep("reused", length(reuse_ids)), reuse_ids),
    setNames(rep("built",  length(build_ids)), build_ids)
  )[colnames(eMLG_sorted)]

  list(
    eMLG   = eMLG_sorted,
    groups = groups[.(colnames(eMLG_sorted)), on = "group_id"],
    source = src
  )
}


## ---------------------------------------------------------
## Example usage (not run) -- see dev/R/moduleA_sorting_phenomenon.R for the
## full A1/A2/A3 runner; this is the essence.
## ---------------------------------------------------------
if (FALSE) {

  library(data.table)
  source("dev/R/Ohta.R")
  source("dev/R/parallelism_stats.R")
  source("dev/R/eMLG_parallelism.R")

  e1 <- new.env(); load("./data/hybrids_only_maf005.Rdata", envir = e1)
  GTs_hyb <- e1$GTs_hybrids_005
  e2 <- new.env(); load("./data/hybrids_and_parents_maf005.Rdata", envir = e2)
  GTs_wp <- e2$GTs_with_parents; sample_data <- e2$sample_data_with_parents
  map <- e2$map_hyb_005
  parent_ids  <- sample_data[grepl("_parent$", Population), Sample_ID]
  GTs_parents <- GTs_wp[parent_ids, , drop = FALSE]
  DI_vec <- setNames(map$DiagnosticIndex, map$marker)

  clust  <- readRDS("./data/eMLG_5loci_0025_cM05.rds")
  groups <- clust$groups; eMLG <- clust$eMLG

  aqu_pops <- "aquilonia_parent"; pol_pops <- "polyctena_parent"
  hybrid_pops <- setdiff(unique(sample_data$Population), c(aqu_pops, pol_pops))

  ## A1
  prep <- ohta_fast_prepare(GTs_wp, pops = sample_data$Population)
  snp <- parallelism_stats(prep, hybrid_pops = hybrid_pops, aqu_pops = aqu_pops,
                           pol_pops = pol_pops, DI = DI_vec, min_DI = -15,
                           sort_th = 0.5, fix_th = 0.15, null_prob = 0.5)
  sorted <- snp[differentiated == TRUE & sort_class != "unsorted", marker]

  ## A2
  se <- build_sorted_eMLG(groups, eMLG, GTs_hyb, sorted, cores = 8)

  ## A3
  units <- colnames(se$eMLG)
  umem  <- setNames(groups[.(units), on = "group_id", members], units)
  par_cons  <- build_group_consensus(umem, GTs_hyb, GTs_parents, cores = 8)
  GTs_units <- rbind(se$eMLG, par_cons)
  pops_u    <- sample_data[match(rownames(GTs_units), Sample_ID), Population]
  DI_u      <- cluster_DI(groups, units, DI_vec, di_agg = "max")
  prep_u    <- ohta_fast_prepare(GTs_units, pops = pops_u)
  cl <- parallelism_stats(prep_u, hybrid_pops = hybrid_pops, aqu_pops = aqu_pops,
                          pol_pops = pol_pops, DI = DI_u, min_DI = -15,
                          sort_th = 0.5, fix_th = 0.15, null_prob = 0.5)
  cl[differentiated == TRUE, .N, by = sort_class][order(-N)]

}
