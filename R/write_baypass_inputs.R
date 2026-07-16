# Shared logic for writing BayPass input files (genotype counts, poolsize,
# PC1/PC2 covariates), reusing the already-computed pruned marker set
# rather than re-running LD-pruning/eMLG generation.
#
# Extracted into one function so the poolsize-ordering fix only needs to
# exist in one place: table(pop) sorts populations ALPHABETICALLY, but the
# genotype-count columns and PC1/PC2 covariates (via unique(pop) /
# !duplicated()) are in first-appearance order -- these only coincidentally
# agreed for 2 of 20 populations in this dataset (Aland first, Vuosaari
## alphabetically last too). Confirmed directly: 17 of 20 populations had a
# mismatched poolsize in every BayPass run this pipeline produced before
# this fix (e.g. Tvarminne is 3rd by first appearance but 19th
# alphabetically, shifting everything between them by one). Since BayPass
# uses poolsize to model sampling variance around each population's allele
# frequency estimate, this silently distorted both Omega and the BF/eBPis
# association statistics, not just added noise.

library(data.table)

#' Write BayPass input files (genotype counts, poolsize, PC1/PC2 covariates)
#'
#' @param GTs Genotype matrix, individuals x markers, rows aligned 1:1 with
#'   `sample_data`.
#' @param map data.table with at least a `marker` column, aligned to `GTs`'s
#'   columns.
#' @param sample_data data.table with Population, Sample_ID, PC1, PC2 --
#'   rows aligned 1:1 with `GTs`.
#' @param pruned_markers Character vector of markers to use for the pruned
#'   genotype file (`u_DIEM.geno_pruned`); the full (unpruned) file always
#'   uses every column of `GTs`.
#' @param out_folder Output directory (created if missing).
#' @param exclude_population Population name(s) to drop before writing, or
#'   `NULL` (default) to include everyone.
write_baypass_inputs <- function(GTs, map, sample_data, pruned_markers,
                                  out_folder, exclude_population = NULL) {

  dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)

  stopifnot(
    "pruned_markers has entries not present in map$marker -- stale relative to current data" =
      all(pruned_markers %in% map$marker),
    "GTs and sample_data row counts don't match" = nrow(GTs) == nrow(sample_data)
  )

  keep <- rep(TRUE, nrow(sample_data))
  if (!is.null(exclude_population)) {
    keep <- !sample_data$Population %in% exclude_population
    message(
      "Excluding: ", paste(exclude_population, collapse = ", "),
      " (", sum(!keep), " samples)"
    )
  }
  sample_data_ex <- sample_data[keep]
  GTs_ex <- GTs[keep, ]

  cat("Populations:", uniqueN(sample_data_ex$Population), "\n")

  pop <- sample_data_ex$Population

  message("Writing pruned genotype counts...")
  baypass_pruned <- do.call(cbind, lapply(unique(pop), function(y){
    t(apply(GTs_ex[pop==y, map$marker %in% pruned_markers], 2, function(x){
      c(length(which(x==0))*2+length(which(x==1)),
        length(which(x==2))*2+length(which(x==1))
      )
    }))
  }))
  write.table(baypass_pruned, file=paste0(out_folder,"u_DIEM.geno_pruned"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  ## table(pop)[unique(pop)] forces poolsize into the SAME order as the
  ## genotype columns and PC1/PC2 covariates -- see file header.
  poolsize <- t(as.vector(table(pop)[unique(pop)]))
  write.table(poolsize, file=paste0(out_folder,"u_DIEM.size"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  message("Writing full (unpruned) genotype counts...")
  baypass_full <- do.call(cbind, lapply(unique(pop), function(y){
    t(apply(GTs_ex[pop==y, ], 2, function(x){
      c(length(which(x==0))*2+length(which(x==1)),
        length(which(x==2))*2+length(which(x==1))
      )
    }))
  }))
  write.table(baypass_full, file=paste0(out_folder,"u_DIEM.geno"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  message("Writing PC1/PC2 covariates...")
  pc1_env <- sample_data_ex[!duplicated(Population), PC1]
  write.table(t(as.matrix(pc1_env)), file=paste0(out_folder,"u.PC1"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  pc2_env <- sample_data_ex[!duplicated(Population), PC2]
  write.table(t(as.matrix(pc2_env)), file=paste0(out_folder,"u.PC2"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  stopifnot(
    "genotype and covariate population order mismatch" =
      identical(unique(pop), sample_data_ex[!duplicated(Population), Population]),
    "poolsize count doesn't match number of populations" = ncol(poolsize) == length(unique(pop)),
    "pruned genotype file has wrong number of population columns (expect 2x populations, ref+alt)" =
      ncol(baypass_pruned) == 2 * length(unique(pop)),
    "full genotype file has wrong number of population columns" = ncol(baypass_full) == 2 * length(unique(pop))
  )

  message("Done.")
  cat("Populations (n=", length(unique(pop)), "):\n", sep = "")
  print(unique(pop))
  cat("\nMarkers: ", nrow(baypass_pruned), " (pruned), ", nrow(baypass_full), " (full)\n", sep = "")
  cat("Files written to ", out_folder, ":\n", sep = "")
  print(list.files(out_folder))

  invisible(list(pop_order = unique(pop), n_pruned = nrow(baypass_pruned), n_full = nrow(baypass_full)))
}
