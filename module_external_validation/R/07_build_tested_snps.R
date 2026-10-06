## module_external_validation / 07: choose the SNP to test in the external samples for each
## outlier REGION (Stage-2 group containing >=1 raw-threshold Stage-1 unit) of every analysis.
##
## Hierarchy per region (units ordered by original association statistic, strongest first):
##   tier "best"        : the best marker (the SNP actually tested by BayPass) of the TOP unit, if it is
##                        usable in the external data (MAF>0.1, unique position, <=20% missing)
##   tier "best_other"  : the best marker of another raw-crossing unit of the same region
##   tier "proxy"       : a member of the TOP unit's LD cluster that is usable externally, chosen by
##                        highest |r| to the top unit's best marker in the ORIGINAL hybrids (|r| >= 0.8)
##   NA                 : region not testable
## Orientation caveat: effects are later reported for the VCF ALT allele. Original allele coding is not
## recoverable from the DIEM-derived genotypes, so p_orig_coded (mean dosage/2 of the original coded
## allele) and p_ext_alt are stored to allow orientation to be inferred/checked downstream.
##
## Run from the repo root:  Rscript module_external_validation/R/07_build_tested_snps.R
## Reads : results/outlier_units_vs_vcf.tsv (from 00_), data/ext_geno_maf10.rds, original stage-1/2 + genotypes
## Writes: results/tested_snps.tsv, results/tested_snps_summary.txt

suppressPackageStartupMessages(library(data.table))
OUT <- "module_external_validation"; MAXMISS <- 0.2; MIN_R <- 0.8
units <- fread(file.path(OUT, "results/outlier_units_vs_vcf.tsv"))
ext <- readRDS(file.path(OUT, "data/ext_geno_maf10.rds")); emap <- ext$map
usable <- emap[fmiss21 <= MAXMISS, marker]

## ---- unit -> Stage-2 region (core SNP assignment, as in module_BayPass/R/repeatability_comparison.R)
s1 <- as.data.table(readRDS("module0_ld_pruning/data/pruned_stage1.rds")$clusters)[n_snps >= 5L]
s2 <- readRDS("module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05.rds")
g2 <- as.data.table(s2$groups); m2g <- g2[, .(marker = unlist(members)), by = .(region = group_id)]; setkey(m2g, marker)
u2r <- data.table(group_id = paste0("S1_", s1$CL_id), region = m2g[.(s1$core_snp), region])
stopifnot(!anyNA(u2r$region))
units <- merge(units, u2r, by = "group_id", all.x = TRUE); stopifnot(!anyNA(units$region))

## ---- original genotypes for proxy r and original allele frequency
e <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e); G <- e$GTs_hybrids_005
b <- readRDS("module_manuscript_rho05/data/moduleB_stage1_units_bestsnp.rds")
members <- setNames(b$groups$members, b$groups$group_id)
orig_af <- function(mk) colMeans(G[, mk, drop = FALSE], na.rm = TRUE) / 2

pick <- function(reg) {
  u <- reg[order(-stat)]
  top <- u[1]
  if (top$best_marker %in% usable) return(list(tested = top$best_marker, tier = "best", r = 1, unit = top$group_id))
  for (i in seq_len(nrow(u))[-1]) if (u$best_marker[i] %in% usable) return(list(tested = u$best_marker[i], tier = "best_other", r = 1, unit = u$group_id[i]))
  cand <- intersect(members[[top$group_id]], usable)
  if (length(cand)) {
    r <- suppressWarnings(cor(G[, top$best_marker], G[, cand, drop = FALSE], use = "pairwise.complete.obs"))[1, ]
    r <- r[is.finite(r)]
    if (length(r) && max(abs(r)) >= MIN_R) { k <- names(which.max(abs(r))); return(list(tested = k, tier = "proxy", r = r[[k]], unit = top$group_id)) }
  }
  list(tested = NA_character_, tier = NA_character_, r = NA_real_, unit = top$group_id)
}
res <- units[, { p <- pick(.SD); o <- order(-stat)
                 .(top_unit = group_id[o][1], best_marker_top = best_marker[o][1], top_stat = max(stat), n_units_region = .N,
                   tested_marker = p$tested, tier = p$tier, r_to_best = p$r, tested_unit = p$unit) },
             by = .(analysis, region)]
mi <- match(res$tested_marker, emap$marker)
res[, `:=`(maf21 = emap$maf21[mi], fmiss21 = emap$fmiss21[mi], REF = emap$REF[mi], ALT = emap$ALT[mi])]
## p_ext_alt must be ALT allele frequency (MAF is folded): compute from genotypes
gx <- ext$geno
res[!is.na(tested_marker), p_ext_alt := colMeans(gx[, tested_marker], na.rm = TRUE) / 2]
res[is.na(tested_marker), p_ext_alt := NA_real_]
res[!is.na(tested_marker), p_orig_coded := orig_af(tested_marker)]
res[, `:=`(Chr = sub(":.*", "", tested_marker), Pos = as.integer(sub(".*:", "", tested_marker)))]
fwrite(res, file.path(OUT, "results/tested_snps.tsv"), sep = "\t")

sink(file.path(OUT, "results/tested_snps_summary.txt"), split = TRUE)
cat("Outlier regions per analysis and testability in the external samples\n")
print(res[, .(regions = .N, testable = sum(!is.na(tested_marker)), best = sum(tier == "best", na.rm = TRUE),
              best_other = sum(tier == "best_other", na.rm = TRUE), proxy = sum(tier == "proxy", na.rm = TRUE)), by = analysis])
cat("\nDistinct tested SNPs per analysis:\n"); print(res[!is.na(tested_marker), .(snps = uniqueN(tested_marker)), by = analysis])
cat("\nProxy |r| to original best marker:\n"); print(summary(abs(res[tier == "proxy", r_to_best])))
sink()
