## module_external_validation / 09: infer, for SNPs shared with the original hybrid data, whether the
## VCF ALT allele equals the ORIGINAL CODED allele (DIEM-polarised genotype coding).
##
## Fact established first (07/orientation exploration): in the original data the coded allele (dosage 2)
## is the F. aquilonia-associated allele for 100% of markers with DI > -25 (and 77% at DI > -50), i.e.
## coded dosage ~ "aquilonia-allele count"; p_aqu / p_pol = coded-allele frequency in the 15+15
## allopatric reference parents (data/hybrids_and_parents_maf005.Rdata).
##
## Method (ancestry axis, allele-coding independent):
##  1. SNPs shared by the external VCF set and the original data with |p_aqu - p_pol| >= D_MIN.
##  2. PCA of the 21 external hybrids on these SNPs (PCA is invariant to per-SNP allele coding);
##     PC1 = ancestry axis z_i (sign g unknown).
##  3. Per SNP, cor(ALT dosage, z) and d_j = p_aqu - p_pol give o_j(g) = +1 (ALT = coded) if
##     g * sign(cor) * sign(d) > 0.
##  4. Global sign g from the model fit: for each g, per-individual aquilonia proportion q_i(g) by
##     least squares of coded dosage/2 on p_pol + q d, and SSE(g). The g<->-g symmetry is exact only for
##     SNPs with p_aqu + p_pol = 1, so the SSE (levels of allele frequency) breaks it.
##
## Run from the repo root:  Rscript module_external_validation/R/09_orientation_ancestry.R
## Reads : data/orig_parent_af.rds (made inline below from the parents file), data/ext_geno_maf10.rds
## Writes: data/orientation_ancestry.rds (per-SNP table), results/orientation_ancestry_summary.txt, results/orientation_ancestry.png

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
OUT <- "module_external_validation"; D_MIN <- 0.6; MAXMISS <- 0.2
paf_f <- file.path(OUT, "data/orig_parent_af.rds")
if (!file.exists(paf_f)) {
  e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
  G <- e$GTs_with_parents; sd <- as.data.table(e$sample_data_with_parents); mp <- as.data.table(e$map_hyb_005)
  ia <- which(sd$Population == "aquilonia_parent"); ip <- which(sd$Population == "polyctena_parent")
  mp[, `:=`(p_aqu = colMeans(G[ia, ], na.rm = TRUE) / 2, p_pol = colMeans(G[ip, ], na.rm = TRUE) / 2)]
  saveRDS(mp[, .(marker, Chr, Pos, Polarity, DiagnosticIndex, maf_hyb, p_aqu, p_pol)], paf_f)
}
paf <- readRDS(paf_f)
d <- readRDS(file.path(OUT, "data/ext_geno_maf10.rds")); X <- d$geno; emap <- d$map

sh <- merge(emap[fmiss21 <= MAXMISS, .(marker, REF, ALT)], paf, by = "marker")
sh[, dd := p_aqu - p_pol]
sel <- sh[is.finite(dd) & abs(dd) >= D_MIN]
cat("shared SNPs:", nrow(sh), " strongly differentiated (|d| >=", D_MIN, "):", nrow(sel), "\n")
Xs <- X[, sel$marker]; mu <- colMeans(Xs, na.rm = TRUE); w <- which(is.na(Xs), arr.ind = TRUE); Xs[w] <- mu[w[, 2]]
keepv <- apply(Xs, 2, sd) > 0; sel <- sel[keepv]; Xs <- Xs[, keepv]
pc <- prcomp(Xs, scale. = TRUE); z <- pc$x[, 1]
cat("PC1 variance share:", round(summary(pc)$importance[2, 1], 3), "PC2:", round(summary(pc)$importance[2, 2], 3), "\n")

cr <- as.numeric(cor(Xs, z)); sel[, `:=`(cor_z = cr, abs_cor = abs(cr))]
fit <- function(g) {
  same <- g * sign(sel$cor_z) * sign(sel$dd) > 0                    # ALT == coded allele?
  C <- Xs; C[, !same] <- 2 - C[, !same]                              # coded-allele dosage
  Pp <- sel$p_pol; dv <- sel$dd
  q <- as.numeric((sweep(C / 2, 2, Pp) %*% dv) / sum(dv^2))
  res <- sweep(C / 2, 2, Pp) - outer(q, dv)
  list(same = same, q = q, sse = sum(res^2), mean_q = mean(q), sse_asym = sum(res[, abs(sel$p_aqu + sel$p_pol - 1) > 0.4]^2), n_asym = sum(abs(sel$p_aqu + sel$p_pol - 1) > 0.4))
}
fp <- fit(+1); fm <- fit(-1)
cat(sprintf("\nSSE(g=+1) = %.1f, SSE(g=-1) = %.1f  (ratio %.3f); on asymmetric SNPs (|pA+pP-1|>0.4, n=%d): %.1f vs %.1f\n",
            fp$sse, fm$sse, fp$sse / fm$sse, fp$n_asym, fp$sse_asym, fm$sse_asym))
cat("mean aquilonia proportion q: g=+1:", round(fp$mean_q, 3), " g=-1:", round(fm$mean_q, 3), "\n")
g_best <- if (fp$sse < fm$sse) +1 else -1; fb <- if (g_best == 1) fp else fm
sel[, `:=`(ext_alt_is_coded = fb$same, g = g_best)]
cat("chosen g =", g_best, "; fraction of SNPs with ALT = coded allele:", round(mean(sel$ext_alt_is_coded), 3), "\n")
cat("q_i (chosen g):", round(sort(fb$q), 2), "\n")
## sanity: AF prediction for the strongest SNPs
sel[, p_ext_alt := colMeans(Xs) / 2]
sel[, p_ext_coded := fifelse(ext_alt_is_coded, p_ext_alt, 1 - p_ext_alt)]
sel[, p_pred := p_pol + mean(fb$q) * dd]
cat("cor(ext coded-allele AF, AF predicted from inferred mean ancestry) =", round(cor(sel$p_ext_coded, sel$p_pred), 3),
    " [with g flipped:", round(cor(fifelse(sel$ext_alt_is_coded, 1 - sel$p_ext_coded, sel$p_ext_coded), sel$p_pred), 3), "]\n")
saveRDS(list(table = sel, q = setNames(fb$q, rownames(X)), g = g_best, sse = c(plus = fp$sse, minus = fm$sse), z = setNames(z, rownames(X))),
        file.path(OUT, "data/orientation_ancestry.rds"))
cl <- fread(file.path(OUT, "data/new_sample_climate_pca_wc.tsv"))
cat("\ncor(ancestry q_i, latitude) =", round(cor(fb$q, cl$latitude), 2), " cor(q_i, climate PC1) =", round(cor(fb$q, cl$PC1), 2), " cor(q_i, PC2) =", round(cor(fb$q, cl$PC2), 2), "\n")
png(file.path(OUT, "results/orientation_ancestry.png"), 1500, 700, res = 150); par(mfrow = c(1, 2))
hist(sel$cor_z * sign(sel$dd), breaks = 60, main = "cor(ALT dosage, ancestry axis) x sign(d)", xlab = "")
plot(sort(fb$q), pch = 19, ylab = "aquilonia proportion q_i", xlab = "individual (sorted)", main = "inferred ancestry of the 21 hybrids"); dev.off()
