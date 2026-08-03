## =========================================================
## Module B: BayPass input for the eMLG-level climate association
## =========================================================
## Builds a BayPass -countdatafile in which each "marker" is one LD cluster's
## eMLG (its consensus genotype across member markers), rather than an
## individual SNP. This lets us test whether the CLUSTER AS A UNIT tracks the
## climate covariate -- a stronger candidate criterion than ">=5 member SNPs
## individually cross BF(dB)>=15", which the raw candidate list uses (see
## moduleB_climate_candidates.R). Only the primary configuration is built:
## Aland-excluded, and the association is run against the fixed LD-pruned Omega.
##
## Reuses, unchanged, the per-population inputs of the SNP-level primary run,
## because Omega and the PC covariates are properties of the 19 populations, not
## of the marker set:
##   aland_excluded/omega_mat_omega.out  (fixed LD-pruned Omega)
##   aland_excluded/u.PC1, u.PC2         (per-population climate covariates)
##   aland_excluded/u_DIEM.size          (per-population pool sizes)
## The population COLUMN ORDER of u_eMLG.geno must match those files. We rebuild
## the pool sizes from sample_data (first-appearance order, Aland dropped) and
## assert they equal aland_excluded/u_DIEM.size -- if the order disagreed, that
## check fails rather than silently misaligning covariates with genotypes.
##
## eMLG dosages are continuous in [0,2]; per the established convention
## (formica_hybrid_old/R/baypass_eMLG.R) they are round()ed to {0,1,2} per
## individual before the diploid allele-count aggregation.
##
## Writes: aland_excluded_eMLG/u_eMLG.geno  (32,840 rows x 2*19 cols)
##         aland_excluded_eMLG/eMLG_group_order.txt  (group_id per geno row)
##         aland_excluded_eMLG/{omega_mat_omega.out,u.PC1,u.PC2,u_DIEM.size} (copies)
## Run from the formica_hybrid repo root.

suppressMessages(library(data.table))

SRC     <- "aland_excluded"                 # reuse this run's Omega/covariates/poolsize
OUTDIR  <- "aland_excluded_eMLG"
EXCLUDE <- "Aland"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## ---- inputs -------------------------------------------------------------
e <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e)
sd <- as.data.table(e$sample_data)                       # Population, Sample_ID, PC1, PC2
E  <- readRDS("module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds")$eMLG      # inds x has_eMLG clusters, dosage [0,2]

## ---- align eMLG rows to sample_data, drop the excluded population -------
stopifnot("eMLG individuals not all in sample_data" = all(rownames(E) %in% sd$Sample_ID))
sd <- sd[Sample_ID %in% rownames(E)]                     # (all 164; defensive)
sd <- sd[Population != EXCLUDE]                           # -> 154 inds, 19 pops
E  <- E[sd$Sample_ID, , drop = FALSE]                    # reorder eMLG rows to sd order
pop <- sd$Population
pop_order <- unique(pop)                                  # first-appearance order
cat("Populations (n=", length(pop_order), "): ", paste(pop_order, collapse = ", "), "\n", sep = "")
cat("Individuals: ", nrow(E), "   eMLG clusters: ", ncol(E), "\n", sep = "")

## ---- round dosages to {0,1,2}; NA preserved (round(NA)=NA) --------------
Er <- round(E)
stopifnot(all(Er %in% c(0, 1, 2, NA)))

## ---- per-population diploid allele counts (same formula as the SNP run) -
## For each population column-pair: (#0*2 + #1, #2*2 + #1). NAs simply drop
## out of the counts, so a cluster's per-population total = 2 * (# non-NA inds).
geno <- do.call(cbind, lapply(pop_order, function(y) {
  t(apply(Er[pop == y, , drop = FALSE], 2, function(x) {
    c(sum(x == 0, na.rm = TRUE) * 2 + sum(x == 1, na.rm = TRUE),
      sum(x == 2, na.rm = TRUE) * 2 + sum(x == 1, na.rm = TRUE))
  }))
}))
stopifnot(nrow(geno) == ncol(E), ncol(geno) == 2 * length(pop_order))

## ---- pool sizes rebuilt here, then asserted == the reused SRC file ------
poolsize <- as.integer(table(pop)[pop_order])
src_size <- scan(file.path(SRC, "u_DIEM.size"), what = integer(), quiet = TRUE)
stopifnot(
  "rebuilt pool sizes disagree with the reused run's u_DIEM.size -- population ORDER mismatch; reusing its Omega/covariates would misalign" =
    identical(poolsize, src_size)
)
cat("Pool-size order check vs ", SRC, "/u_DIEM.size: OK\n", sep = "")

## ---- write geno + the group_id order (for joining BF back to clusters) --
write.table(geno, file.path(OUTDIR, "u_eMLG.geno"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
writeLines(colnames(E), file.path(OUTDIR, "eMLG_group_order.txt"))

## ---- stage the reused per-population files so the run is self-contained -
for (f in c("omega_mat_omega.out", "u.PC1", "u.PC2", "u_DIEM.size"))
  file.copy(file.path(SRC, f), file.path(OUTDIR, f), overwrite = TRUE)

cat("\nWrote ", OUTDIR, "/u_eMLG.geno (", nrow(geno), " eMLGs x ", ncol(geno),
    " cols) + eMLG_group_order.txt, staged Omega/covariates/poolsize.\n", sep = "")
