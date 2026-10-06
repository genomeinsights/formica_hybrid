# Check per-chromosome AIMs files (input to the SLiM founder set-up) against the
# fixation-level table and LD-cluster table they were built from.
#
# Usage:
#   Rscript check_AIMs.R <AIMs_dir> <fixation_levels.tsv> <group_info.rds> [diem_output.bed]
#
# Part 1 (AIMs files, run BEFORE SLiM) -- must be 100% on every chromosome:
#   1a. every input marker within the chromosome length appears exactly once, at position-1
#   1b. the cluster partition in the AIMs files equals the LD-cluster partition in group_info
#       (independent SNPs = their own singleton cluster)
#   1c. pm10_aq == mean(faqu_fix) and pm10_pol == 1 - mean(fpol_fix) over the cluster's OWN markers
#   1d. bug signature: correlation of pm10_aq with the marker's own value vs with chromosome 1's
#       value at the same row number (the original bug gives ~0 vs strongly non-zero)
#
# Part 2 (optional, AFTER SLiM + diem): given one diem_boot*_output.bed, the simulated parental
#   difference |p_aq - p_pol| per marker should track the INTENDED founder difference (cluster-mean
#   faqu_fix vs 1 - fpol_fix from the input table, independent of what the AIMs files contain)
#   on EVERY chromosome (not just chromosome 1).

options(scipen = 999)
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: Rscript check_AIMs.R <AIMs_dir> <fixation_levels.tsv> <group_info.rds> [diem_output.bed]")
aims_dir <- args[1]; fix_file <- args[2]; grp_file <- args[3]
bed_file <- if (length(args) >= 4) args[4] else NA

TOL <- 1e-9
Ls <- c(16449812, 18558483, 15805178, 16666556, 14957794, 10925039, 13372598, 13028660, 10584038, 13832892, 11830133,
        11260715, 11449430, 7730741, 10965718, 11699404, 11021207, 8644910, 7297899, 9820556, 10314430, 7671390,
        NA, 6269374, 13114028, 9969264, 7047176) - 1

## ---- inputs ----------------------------------------------------------------
fx <- fread(fix_file)
fx[, chr := as.integer(sub("chromosome_", "", chromosome))]
fx[, marker := paste0(chr, ":", as.integer(position))]
stopifnot(!anyDuplicated(fx$marker))
fx <- fx[position - 1 <= Ls[chr]]              # markers SLiM can place

grp <- as.data.table(readRDS(grp_file))
grp[, marker := paste0(as.integer(sub("chromosome_", "", chromosome)), ":", as.integer(position))]

A <- rbindlist(lapply(list.files(aims_dir, pattern = "^AIMs_ch[0-9]+\\.txt$", full.names = TRUE), function(f)
  fread(f, header = FALSE, col.names = c("chr", "pos0", "cluster", "pm10_aq", "pm10_pol"))))
A[, marker := paste0(chr, ":", as.integer(pos0 + 1))]

cat(sprintf("AIMs rows: %d | input markers (within chromosome length): %d | group_info markers: %d\n",
            nrow(A), nrow(fx), nrow(grp)))

## ---- 1a. coverage ------------------------------------------------------------
dup  <- sum(duplicated(A$marker))
miss <- setdiff(fx$marker, A$marker)
extra <- setdiff(A$marker, fx$marker)
cat(sprintf("\n[1a] duplicated AIMs rows: %d | input markers missing from AIMs: %d | AIMs markers not in input: %d\n",
            dup, length(miss), length(extra)))

## ---- 1b. cluster partition ----------------------------------------------------
A <- merge(A, fx[, .(marker, faqu_fix, fpol_fix)], by = "marker", all.x = TRUE)
A[, group_id := grp$group_id[match(marker, grp$marker)]]
A[is.na(group_id), group_id := paste0("indep_", marker)]   # independent SNP = own group
part <- A[, .(n_aims_clusters_per_group = uniqueN(cluster)), by = .(chr, group_id)]
part2 <- A[, .(n_groups_per_aims_cluster = uniqueN(group_id)), by = .(chr, cluster)]
cat(sprintf("[1b] LD groups split over >1 AIMs cluster: %d | AIMs clusters mixing >1 LD group: %d\n",
            sum(part$n_aims_clusters_per_group > 1), sum(part2$n_groups_per_aims_cluster > 1)))

## ---- 1c. values -----------------------------------------------------------------
A[, exp_aq := mean(faqu_fix, na.rm = TRUE), by = .(chr, cluster)]
A[, exp_pol := 1 - mean(fpol_fix, na.rm = TRUE), by = .(chr, cluster)]
A[, ok := abs(pm10_aq - exp_aq) < TOL & abs(pm10_pol - exp_pol) < TOL]

## ---- 1d. bug signature ---------------------------------------------------------
setorder(A, chr, pos0)
A[, row_in_chr := seq_len(.N), by = chr]
c1 <- A[chr == 1]
A[, chr1_same_row := c1$faqu_fix[row_in_chr]]

safe_cor <- function(x, y) { k <- is.finite(x) & is.finite(y); if (sum(k) < 10 || sd(x[k]) == 0 || sd(y[k]) == 0) NA_real_ else cor(x[k], y[k]) }
per_chr <- A[, .(markers = .N,
                 pct_values_correct = round(100 * mean(ok, na.rm = TRUE), 2),
                 cor_own_marker = round(safe_cor(pm10_aq, faqu_fix), 3),
                 cor_chr1_same_row = round(safe_cor(pm10_aq, chr1_same_row), 3)), by = chr]
cat("\n[1c/1d] per chromosome (pct_values_correct must be 100; cor_own_marker should be high and\n",
    "        cor_chr1_same_row near 0 except on chromosome 1):\n", sep = "")
print(per_chr, row.names = FALSE)

pass1 <- dup == 0 && length(miss) == 0 && length(extra) == 0 &&
         all(part$n_aims_clusters_per_group == 1) && all(part2$n_groups_per_aims_cluster == 1) &&
         all(A$ok %in% TRUE)
cat(sprintf("\nPART 1 (AIMs files): %s\n", if (pass1) "PASS" else "FAIL"))

## ---- 2. optional: simulated parents after SLiM + diem --------------------------
if (!is.na(bed_file)) {
  h <- readLines(bed_file, n = 2)
  inds <- strsplit(strsplit(h[2], "\t")[[1]][10], "|", fixed = TRUE)[[1]]
  s <- fread(bed_file, skip = 2, header = FALSE, sep = "\t", select = c(1, 3, 10),
             colClasses = list(character = c(1, 10)), showProgress = FALSE)
  S <- do.call(rbind, strsplit(sub("^S", "", s$V10), "", fixed = TRUE))
  frq <- function(ix) { M <- matrix(suppressWarnings(as.integer(S[, ix])), nrow(S)); rowMeans(M, na.rm = TRUE) / 2 }
  sim <- data.table(marker = paste0(sub("ch", "", s$V1), ":", as.integer(s$V3)),
                    sim_dp = abs(frq(grep("^aq_", inds)) - frq(grep("^pol_", inds))))
  sim <- merge(sim, A[, .(marker, chr, exp_dp = abs(exp_aq - exp_pol))], by = "marker")   # INTENDED values (input table), not the AIMs file
  per_chr2 <- sim[, .(markers = .N, spearman_sim_vs_founder_expectation =
                        round(suppressWarnings(cor(sim_dp, exp_dp, method = "spearman", use = "complete")), 3)), by = chr]
  setorder(per_chr2, chr)
  cat(sprintf("\n[2] %s: simulated parental |p_aq - p_pol| vs intended founder difference (input table)\n", basename(bed_file)))
  cat("    (should be clearly positive on EVERY chromosome; the original bug gives ~0 everywhere except chromosome 1)\n")
  print(per_chr2, row.names = FALSE)
}
