## module_external_validation / 00: audit of the external hybrid VCF against the
## original hybrid data set and the BayPass outlier units.
##
## Questions answered
##  1. Same assembly / coordinates as the original hybrid data?  (site overlap)
##  2. How many of the outlier units' BEST SNPs (the SNP actually tested in the
##     Stage-1 BayPass scan) exist in the external VCF, and how many units have
##     at least one LD-cluster member available as a proxy?
##  3. How many of those are polymorphic enough (MAF in the 21 samples) to test?
##
## Run from the repo root:  Rscript module_external_validation/R/00_audit_external_vcf.R
## Needs bcftools on PATH. Reads (not modified):
##   data/HybridSamples_SNPQ30...FiSeDe.vcf.gz (21 samples; INFO NS/AC/AF are STALE,
##     computed on a 384-sample call -- MAF is therefore recomputed here)
##   module_localscore_crosscheck/data/fullsnp_chr_pos.rds
##   module_manuscript_rho05/data/moduleB_stage1_units_bestsnp.rds
##   module_manuscript_rho05/module_BayPass/R/analysis_registry.R (+ its Stage-1 result files)
## Writes: module_external_validation/data/vcf_sites_maf21.tsv.gz
##         module_external_validation/results/outlier_units_vs_vcf.tsv
##         module_external_validation/results/audit_summary.txt

suppressPackageStartupMessages(library(data.table))
OUT <- "module_external_validation"
VCF <- "data/HybridSamples_SNPQ30.biall.fixedHeader.indDP.minDP8.hwe.376inds.AN10percMiss_FiSeDe.vcf.gz"
sites_f <- file.path(OUT, "data/vcf_sites_maf21.tsv.gz")

if (!file.exists(sites_f)) {
  cmd <- sprintf("bcftools +fill-tags %s -Ou -- -t MAF,F_MISSING 2>/dev/null | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%MAF\\t%%F_MISSING\\n' | sed 's/^chromosome_/Chr/' | gzip > %s", VCF, sites_f)
  stopifnot(system(cmd) == 0)
}
v <- fread(sites_f, col.names = c("Chr", "Pos", "REF", "ALT", "MAF", "FMISS")); v[, snp := paste0(Chr, ":", Pos)]
sink(file.path(OUT, "results/audit_summary.txt"), split = TRUE)
cat("External VCF: ", nrow(v), "biallelic SNPs; MAF in the 21 samples >0:", sum(v$MAF > 0),
    " >0.05:", sum(v$MAF > 0.05), " >0.1:", sum(v$MAF > 0.1), " >0.2:", sum(v$MAF > 0.2), "\n")

## ---- 1. overlap with the original hybrid SNP set --------------------------------
p <- as.data.table(readRDS("module_localscore_crosscheck/data/fullsnp_chr_pos.rds"))
cat("Original hybrid SNPs:", nrow(p), " also in external VCF:", nrow(merge(p, v, by = c("Chr", "Pos"))), "\n")

## ---- 2. outlier units ------------------------------------------------------------
repo <- getwd(); mod <- file.path(repo, "module_manuscript_rho05")
setwd(file.path(repo, "module_manuscript_rho05/module_BayPass")); source("config/paths.R"); source("R/analysis_registry.R"); setwd(repo)
b <- readRDS(file.path(mod, "data/moduleB_stage1_units_bestsnp.rds"))
st <- as.data.table(b$best$stats); gr <- b$groups
mem <- gr[, .(snp = unlist(members)), by = group_id]
ord <- readLines(file.path(mod, "baypass_stage1/aland_excluded_S1units/S1units_group_order.txt"))
n_mem <- mem[, .N, by = group_id]; n_mem <- setNames(n_mem$N, n_mem$group_id)
mem[, in_vcf := snp %in% v$snp]
n_mem_vcf <- mem[, .(N = sum(in_vcf)), by = group_id]; n_mem_vcf <- setNames(n_mem_vcf$N, n_mem_vcf$group_id)
rows <- list()
for (a in c("PC1", "PC2", "mitoC2", "coastal_inland_C2", "heat_tolerance")) {
  A <- ANALYSIS_REGISTRY[[a]]; f <- file.path(repo, A$stage1_file)
  stat <- if (grepl("betai_reg", f)) fread(f, select = "BF(dB)")[[1]] else if (grepl("contrast", f)) fread(f, select = "log10(1/pval)")[[1]] else fread(f)$obs
  stopifnot(length(stat) == length(ord))
  u <- ord[stat >= A$raw_threshold]
  d <- data.table(analysis = a, group_id = u, stat = stat[stat >= A$raw_threshold])
  d$best_marker <- st$best_marker[match(d$group_id, st$group_id)]
  d$best_in_vcf <- d$best_marker %in% v$snp
  d$n_members <- n_mem[d$group_id]
  d$n_members_in_vcf <- n_mem_vcf[d$group_id]
  d$best_MAF21 <- v$MAF[match(d$best_marker, v$snp)]
  rows[[a]] <- d
}
res <- rbindlist(rows)
fwrite(res, file.path(OUT, "results/outlier_units_vs_vcf.tsv"), sep = "\t")
cat("\nRaw-threshold outlier units, by analysis:\n")
print(res[, .(n_units = .N, best_SNP_in_VCF = sum(best_in_vcf), best_in_VCF_and_MAF21_ge_0.1 = sum(best_in_vcf & best_MAF21 >= 0.1, na.rm = TRUE),
              units_with_any_member_in_VCF = sum(n_members_in_vcf > 0)), by = analysis])
sink()
