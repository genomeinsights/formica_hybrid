## module_external_validation / 08: GRM from the LD-pruned (Stage-1 representative) markers and
## EMMAX of the climate covariates (PC1, PC2; WorldClim-consistent projection) on the 21 external
## hybrids.  Design = GWAS with the environmental variable as the "phenotype" and SNP dosage (ALT
## allele) as predictor, relatedness through the GRM (LDscnR::emmax_setup / emmax_fast; REML variance
## components estimated once per phenotype).
##
##  * GRM: Z = (dosage - 2p)/sqrt(2p(1-p)) on mean-imputed dosages of the 155,665 Stage-1 cluster
##    representatives (MAF>0.1, <=20% missing), K = ZZ'/m. Own implementation (imputed entries
##    contribute 0), avoiding the SNPRelate missing.rate default that silently drops markers.
##  * Scan set: all MAF>0.1 SNPs with <=20% missing (genome-wide scan gives lambda_GC and the
##    reference p-value distribution for the tested outlier-region SNPs).
##  * Phenotypes: raw PC and rank-based inverse-normal PC (PC2 has 4 extreme high-leverage samples).
##  * Tested SNPs: results/tested_snps.tsv (07_). Reported: p, GLS effect (per ALT allele), SE,
##    genome-wide empirical percentile of p. Effect SIGN is for the VCF ALT allele; orientation
##    relative to the original BayPass Beta_is is NOT yet resolved (see 07_ header).
##
## Run from the repo root:  Rscript module_external_validation/R/08_grm_emmax.R
## Writes: data/ext_grm_pruned.rds; results/emmax_tested_snps.tsv; results/emmax_scan_<pheno>.rds;
##         results/grm_emmax_summary.txt; results/emmax_qq.png

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
devtools::load_all("~/gitlab/LDscnR/", quiet = TRUE)
OUT <- "module_external_validation"; MAXMISS <- 0.2

d <- readRDS(file.path(OUT, "data/ext_geno_maf10.rds")); G <- d$geno; map <- d$map
keep <- which(map$fmiss21 <= MAXMISS); G <- G[, keep, drop = FALSE]; map <- map[keep]
meanimp <- function(M) { mu <- colMeans(M, na.rm = TRUE); w <- which(is.na(M), arr.ind = TRUE); M[w] <- mu[w[, 2]]; M }
X <- meanimp(G); n <- nrow(X)
stopifnot(all(is.finite(X)))

## ---- GRM ---------------------------------------------------------------------------------
pr <- readLines(file.path(OUT, "data/ext_pruned_markers_maf10.txt")); pr <- intersect(pr, colnames(X))
p <- colMeans(X[, pr]) / 2; Z <- sweep(X[, pr], 2, 2 * p); Z <- sweep(Z, 2, sqrt(2 * p * (1 - p)), "/")
K <- tcrossprod(Z) / length(pr); dimnames(K) <- list(rownames(X), rownames(X))
saveRDS(list(K = K, markers = pr), file.path(OUT, "data/ext_grm_pruned.rds"))

## ---- phenotypes -----------------------------------------------------------------------------
cl <- fread(file.path(OUT, "data/new_sample_climate_pca_wc.tsv"))
cl[, sample := sub("\\.bam$", "", BAM_ID)]; cl <- cl[match(rownames(X), sample)]; stopifnot(!anyNA(cl$PC1))
rint <- function(y) qnorm((rank(y) - 0.5) / length(y))
Y <- list(PC1 = cl$PC1, PC2 = cl$PC2, PC1_rint = rint(cl$PC1), PC2_rint = rint(cl$PC2))

## ---- EMMAX scan (genome-wide) -------------------------------------------------------------------
prep <- emmax_setup(X, K)
beta_se <- function(prep, y, idx) {   # GLS effect for selected columns, mirrors emmax_fast internals
  re <- emma.REMLE(y, prep$Xo, prep$Kn, eig.R = prep$eigR)
  wv <- 1 / sqrt(re$vg * prep$lam + re$ve)
  yt <- as.numeric(crossprod(prep$V, y)) * wv; xo <- prep$xot * wv; Xtw <- prep$Xt[, idx, drop = FALSE] * wv
  a2 <- sum(xo * xo); ra <- yt - xo * (sum(xo * yt) / a2)
  Rb <- Xtw - outer(xo, as.numeric(crossprod(xo, Xtw)) / a2)
  b <- as.numeric(crossprod(ra, Rb)) / colSums(Rb * Rb)
  rss <- sum(ra * ra) - b^2 * colSums(Rb * Rb); se <- sqrt(rss / prep$df2 / colSums(Rb * Rb))
  list(beta = b, se = se, vg = re$vg, ve = re$ve, h2 = re$vg / (re$vg + re$ve))
}
ts <- fread(file.path(OUT, "results/tested_snps.tsv"), na.strings = c("", "NA"))[!is.na(tested_marker)]
scan <- list(); out <- list(); varcomp <- list()
for (ph in names(Y)) {
  pv <- emmax_fast(prep, Y[[ph]]); names(pv) <- colnames(X); scan[[ph]] <- pv
  lam <- median(qchisq(pv, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)
  a <- sub("_rint$", "", ph)
  tt <- ts[analysis == a]
  if (nrow(tt)) {
    idx <- match(tt$tested_marker, colnames(X)); stopifnot(!anyNA(idx))
    bs <- beta_se(prep, Y[[ph]], idx)
    o <- copy(tt)[, `:=`(pheno = ph, beta_alt = bs$beta, se = bs$se, p = pv[idx],
                         pct_genomewide = vapply(pv[idx], function(x) mean(pv <= x, na.rm = TRUE), 0))]
    out[[ph]] <- o
  }
  bs0 <- beta_se(prep, Y[[ph]], 1L)
  varcomp[[ph]] <- data.table(pheno = ph, n_scan = sum(is.finite(pv)), n_nan = sum(!is.finite(pv)), lambda_GC = lam, h2_pop_structure = bs0$h2,
                              frac_p05 = mean(pv < 0.05, na.rm = TRUE), frac_p01 = mean(pv < 0.01, na.rm = TRUE), min_p = min(pv, na.rm = TRUE))
}
saveRDS(scan, file.path(OUT, "results/emmax_scan_all_phenotypes.rds"))
R <- rbindlist(out); fwrite(R, file.path(OUT, "results/emmax_tested_snps.tsv"), sep = "\t")

## ---- equivalence check against LDscnR::emmax() on the tested SNPs ---------------------------------
chk <- local({ tt <- ts[analysis == "PC1"]; idx <- match(tt$tested_marker, colnames(X))
  e1 <- emmax(Y$PC1, X[, idx, drop = FALSE], K)$pval; cor(log10(e1), log10(scan$PC1[idx])) })

## ---- summary ----------------------------------------------------------------------------------------------
sink(file.path(OUT, "results/grm_emmax_summary.txt"), split = TRUE)
cat("n =", n, "individuals; scan SNPs:", ncol(X), " GRM markers:", length(pr), "\n")
cat("GRM: mean diag", round(mean(diag(K)), 3), " mean off-diag", round(mean(K[upper.tri(K)]), 4), " range off-diag", round(range(K[upper.tri(K)]), 3), "\n")
eg <- eigen(K, symmetric = TRUE); cat("GRM eigenvalue share, top 5 PCs:", round(eg$values[1:5] / sum(eg$values), 3), "\n")
pc <- eg$vectors[, 1:3]
cat("cor(GRM PC1..3, PC1_climate):", round(cor(pc, cl$PC1), 2), "| cor(GRM PC1..3, PC2_climate):", round(cor(pc, cl$PC2), 2), "\n")
cat("cor(GRM PC1..3, latitude):", round(cor(pc, cl$latitude), 2), "| cor(GRM PC1..3, longitude):", round(cor(pc, cl$longitude), 2), "\n")
cat("climate PC1 vs PC2 correlation in the 21 samples:", round(cor(cl$PC1, cl$PC2), 2), "\n")
cat("\nEMMAX genome-wide calibration:\n"); print(rbindlist(varcomp), digits = 3)
cat("\nemmax_fast vs emmax() on tested PC1 SNPs: cor(log10 p) =", round(chk, 6), "\n")
cat("\nTested outlier-region SNPs:\n")
print(R[, .(n_tested = .N, n_p_lt_0.05 = sum(p < 0.05), expected_if_null_genomewide = round(.N * unique(varcomp[[pheno[1]]]$frac_p05), 2),
            median_pct_genomewide = round(median(pct_genomewide), 3), min_p = signif(min(p), 3)), by = pheno])
sink()

q <- rbindlist(lapply(names(scan), function(ph) { p <- sort(scan[[ph]]); data.table(pheno = ph, exp = -log10((seq_along(p) - 0.5) / length(p)), obs = -log10(p)) }))
g <- ggplot(q[sample(.N, min(.N, 4e5))], aes(exp, obs)) + geom_abline(linetype = 2) + geom_point(size = .3) + facet_wrap(~pheno) + theme_classic() +
  labs(x = "expected -log10 p", y = "observed -log10 p", title = "EMMAX genome-wide scan, 21 external hybrids")
ggsave(file.path(OUT, "results/emmax_qq.png"), g, width = 8, height = 6, dpi = 200)
