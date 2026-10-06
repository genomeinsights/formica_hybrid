## module_external_validation / 11: directional concordance of the external EMMAX effects with the
## ORIGINAL BayPass associations, for the SNP tested in each outlier region (PC1, PC2).
##
## 1. Convention of the original Beta_is. The unit genotype fed to BayPass is the best marker's coded dosage
##    (the `flipped` flag only affects how MISSING calls were filled); the population count file lists the
##    dosage-0 allele first. Expectation: Beta_is refers to the NON-coded allele. Checked empirically here:
##    for strongly crossing units (BF >= 15) compare sign(Beta_is) with the sign of the population-level
##    correlation between coded-allele frequency and the covariate (19 populations, Aland excluded).
## 2. Predicted direction of the CODED allele of the tested SNP:
##       e = k * sign(Beta_is[tested_unit]) * sign(r_to_best)
##    k = +1 if Beta_is refers to the coded allele, -1 otherwise (from step 1); r_to_best = signed correlation of
##    the tested SNP's coded dosage with the unit's best marker in the original hybrids (1 for best markers).
## 3. Observed direction: sign(beta_coded) from 08_/10_ (external EMMAX; ALT oriented to coded by Polarity).
##    Covariate signs: the WorldClim-based external PCs were sign-aligned to the stored PC1/PC2.
## Tests per phenotype: one-sided binomial on concordant count; sign-flip permutation of T = sum e_i z_i (z = beta/se).
##
## Run from the repo root:  Rscript module_external_validation/R/11_concordance_test.R
## Writes: results/concordance_tested_snps.tsv, results/concordance_summary.txt, results/concordance.png

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
OUT <- "module_external_validation"; mod <- "module_manuscript_rho05"
S1DIR <- file.path(mod, "baypass_stage1/aland_excluded_S1units")
ord <- readLines(file.path(S1DIR, "S1units_group_order.txt"))
e <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e); G <- e$GTs_hybrids_005; sd <- as.data.table(e$sample_data)
bst <- readRDS(file.path(mod, "data/moduleB_stage1_units_bestsnp.rds"))$best$stats
beta <- list(PC1 = fread(file.path(S1DIR, "PC1_S1units_withOmega_summary_betai_reg.out"), select = c("Beta_is", "BF(dB)")),
             PC2 = fread(file.path(S1DIR, "PC2_S1units_withOmega_summary_betai_reg.out"), select = c("Beta_is", "BF(dB)")))
stopifnot(all(sapply(beta, nrow) == length(ord)))
pops <- sd[Population != "Aland", unique(Population)]
popcov <- sd[Population != "Aland", .(PC1 = PC1[1], PC2 = PC2[1]), by = Population][match(pops, Population)]
rows <- match(sd[Population != "Aland", Sample_ID], rownames(G)); samp_pop <- sd[Population != "Aland", Population]
popfreq <- function(mk) vapply(pops, function(p) mean(G[rows[samp_pop == p], mk], na.rm = TRUE) / 2, 0)

## ---- 1. convention check -----------------------------------------------------------------------------------
sink(file.path(OUT, "results/concordance_summary.txt"), split = TRUE)
conv <- rbindlist(lapply(c("PC1", "PC2"), function(a) {
  idx <- which(beta[[a]]$`BF(dB)` >= 15)
  rbindlist(lapply(idx, function(i) { mk <- bst$best_marker[match(ord[i], bst$group_id)]
    data.table(analysis = a, unit = ord[i], bf = beta[[a]]$`BF(dB)`[i], beta_is = beta[[a]]$Beta_is[i],
               cor_coded_pop = suppressWarnings(cor(popfreq(mk), popcov[[a]]))) }))
}))
conv[, prod_sign := sign(beta_is) * sign(cor_coded_pop)]
cat("1. Beta_is convention. Units with BF >= 15: sign(Beta_is) x sign(cor(coded-allele pop frequency, covariate))\n")
print(conv[, .(n = .N, frac_negative_product = round(mean(prod_sign < 0, na.rm = TRUE), 3), frac_positive_product = round(mean(prod_sign > 0, na.rm = TRUE), 3)), by = analysis])
K <- if (mean(conv$prod_sign < 0, na.rm = TRUE) > 0.5) -1 else +1
cat("=> Beta_is refers to the", if (K == -1) "NON-coded (dosage-0)" else "CODED", "allele; k =", K, "\n\n")

## ---- 2/3. predicted vs observed direction of the coded allele --------------------------------------------------
ts <- fread(file.path(OUT, "results/tested_snps.tsv"), na.strings = c("", "NA"))[!is.na(tested_marker) & analysis %in% c("PC1", "PC2")]
ob <- fread(file.path(OUT, "results/emmax_tested_snps_oriented.tsv"), na.strings = c("", "NA"))
ts[, beta_unit := mapply(function(a, u) beta[[a]]$Beta_is[match(u, ord)], analysis, tested_unit)]
ts[, e_coded := K * sign(beta_unit) * sign(r_to_best)]
D <- merge(ob[, .(pheno, analysis, region, tested_marker, tier, alt_is_coded, beta_coded, se, p)], ts[, .(analysis, region, tested_marker, tested_unit, beta_unit, r_to_best, e_coded)],
           by = c("analysis", "region", "tested_marker"))
D[, `:=`(z = beta_coded / se, concordant = sign(beta_coded) == e_coded)]
fwrite(D, file.path(OUT, "results/concordance_tested_snps.tsv"), sep = "\t")

set.seed(1)
perm_p <- function(e, z, B = 50000) { T0 <- sum(e * z); Tb <- replicate(B, sum(sample(c(-1, 1), length(e), TRUE) * e * z)); (1 + sum(Tb >= T0)) / (B + 1) }
res <- D[, .(n = .N, concordant = sum(concordant), frac = round(mean(concordant), 3),
             p_binom = binom.test(sum(concordant), .N, 0.5, alternative = "greater")$p.value,
             sum_ez = round(sum(e_coded * z), 2), p_perm = perm_p(e_coded, z),
             n_nominal = sum(p < 0.05), concordant_among_nominal = sum(concordant & p < 0.05)), by = pheno]
cat("2. Concordance of external effect direction with the original association (tested SNP per outlier region)\n")
print(res[, `:=`(p_binom = signif(p_binom, 3), p_perm = signif(p_perm, 3))])
cat("\nBy tier (raw PC phenotypes):\n"); print(D[pheno %in% c("PC1", "PC2"), .(n = .N, concordant = sum(concordant)), by = .(pheno, tier)][order(pheno, tier)])
cat("\nPer SNP (raw PC):\n")
print(D[pheno %in% c("PC1", "PC2"), .(pheno, region, tested_marker, tier, e_orig = e_coded, beta_coded = signif(beta_coded, 3), z = round(z, 2), p = signif(p, 3), concordant)][order(pheno, -abs(z))])
sink()
g <- ggplot(D[pheno %in% c("PC1", "PC2")], aes(reorder(tested_marker, e_coded * z), e_coded * z, fill = concordant)) + geom_col() + coord_flip() +
  facet_wrap(~pheno, scales = "free_y") + labs(x = NULL, y = "z (external effect on coded allele) x predicted sign  [>0 = concordant]") +
  scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#D55E00")) + theme_classic(base_size = 9)
ggsave(file.path(OUT, "results/concordance.png"), g, width = 10, height = 6, dpi = 200)
