## module_external_validation / 10: allele-orientation rule and its validation, then application to the
## tested SNPs.
##
## Rule (from 09_): VCF ALT allele == ORIGINAL CODED allele  <=>  DIEM Polarity == 0
## (consistent with DIEM's raw track counting the ALT allele; the project flips Polarity==1 markers,
## coded = 2 - raw). In the original data the coded allele is the F. aquilonia-associated allele for
## 100% of markers with DI > -25.
##
## Validation, independent of how the rule was derived:
##  (a) ancestry-axis inference (09_, uses no Polarity) vs the Polarity rule, by evidence strength |cor|
##  (b) LD-phase check on neighbouring SNP pairs of ANY differentiation level: with orientation o_j
##      (+1 if Polarity 0, else -1), sign(r_ext(j,k)) must equal o_j * o_k * sign(r_orig(j,k)) for pairs in
##      strong LD in both data sets (|r_orig| >= 0.8, |r_ext| >= 0.6). If the rule were wrong for a class of
##      SNPs the agreement would fall toward 50%.
##
## Run from the repo root:  Rscript module_external_validation/R/10_orientation_polarity_rule.R
## Writes: results/orientation_validation.txt, results/emmax_tested_snps_oriented.tsv

suppressPackageStartupMessages(library(data.table))
OUT <- "module_external_validation"
paf <- readRDS(file.path(OUT, "data/orig_parent_af.rds")); setkey(paf, marker)
anc <- readRDS(file.path(OUT, "data/orientation_ancestry.rds"))$table
d <- readRDS(file.path(OUT, "data/ext_geno_maf10.rds")); X <- d$geno; emap <- d$map

sink(file.path(OUT, "results/orientation_validation.txt"), split = TRUE)
## ---- (a) ancestry inference vs Polarity rule, by evidence strength -------------------------------
anc[, rule_alt_is_coded := Polarity == 0]
anc[, bin := cut(abs_cor, c(0, .2, .4, .6, .8, 1), include.lowest = TRUE)]
cat("(a) agreement of ancestry-based orientation (no Polarity used) with the Polarity rule\n")
print(anc[, .(n = .N, agree = round(mean(ext_alt_is_coded == rule_alt_is_coded), 4)), by = bin][order(bin)])
cat("Note: even with a perfect rule the agreement at low |cor| tends to ~50% because the ancestry-based call is itself noise there.\n")
cat("overall at |cor| > 0.6:", round(anc[abs_cor > 0.6, mean(ext_alt_is_coded == rule_alt_is_coded)], 4), " (n =", anc[abs_cor > 0.6, .N], ")\n\n")

## ---- (b) LD-phase check on neighbouring pairs (all shared SNPs, any DI) -----------------------------
e <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e); G <- e$GTs_hybrids_005
sh <- merge(emap[fmiss21 <= 0.2, .(marker, Chr, Pos)], paf[, .(marker, Polarity, DiagnosticIndex)], by = "marker")
sh[, chr_n := as.integer(sub("Chr", "", Chr))]; setorder(sh, chr_n, Pos)
pr <- sh[, .(j = marker[-.N], k = marker[-1], dpos = diff(Pos)), by = chr_n][dpos <= 20000]   # consecutive shared SNPs within 20 kb
set.seed(1); pr <- pr[sample(.N, min(.N, 300000))]
zs <- function(M) { M <- scale(M); M[is.na(M)] <- 0; M }
use <- unique(c(pr$j, pr$k))
Zo <- zs(G[, use]); Ze <- zs(X[, use]); n_o <- nrow(Zo); n_e <- nrow(Ze)
ij <- match(pr$j, use); ik <- match(pr$k, use)
pr[, `:=`(r_o = colSums(Zo[, ij] * Zo[, ik]) / (n_o - 1), r_e = colSums(Ze[, ij] * Ze[, ik]) / (n_e - 1))]
pol <- setNames(sh$Polarity, sh$marker); di <- setNames(sh$DiagnosticIndex, sh$marker)
pr[, `:=`(oj = ifelse(pol[j] == 0, 1, -1), ok = ifelse(pol[k] == 0, 1, -1), di_min = pmin(di[j], di[k]))]
st <- pr[abs(r_o) >= 0.8 & abs(r_e) >= 0.6]
st[, agree := sign(r_e) == oj * ok * sign(r_o)]
st[, agree_norule := sign(r_e) == sign(r_o)]
cat("(b) LD-phase agreement for neighbouring SNP pairs (|r_orig| >= 0.8 & |r_ext| >= 0.6): n pairs =", nrow(st), "\n")
cat("    with Polarity rule applied:", round(mean(st$agree), 4), " | WITHOUT any orientation (raw ALT vs coded):", round(mean(st$agree_norule), 4), "\n")
cat("    pairs where the two SNPs have DIFFERENT polarity (rule matters):", round(mean(st[oj != ok]$agree), 4), "(n =", nrow(st[oj != ok]), ")\n")
cat("    pairs where both SNPs have the same polarity:", round(mean(st[oj == ok]$agree), 4), "(n =", nrow(st[oj == ok]), ")\n")
cat("    non-diagnostic pairs only (both DI <= -25), different polarity:", round(mean(st[oj != ok & di_min <= -25]$agree), 4), "(n =", nrow(st[oj != ok & di_min <= -25]), ")\n\n")
sink()

## ---- apply to the tested SNPs -------------------------------------------------------------------------------
R <- fread(file.path(OUT, "results/emmax_tested_snps.tsv"), na.strings = c("", "NA"))
R[, polarity := paf[tested_marker, Polarity]]
R[, `:=`(alt_is_coded = polarity == 0)]
R[, `:=`(beta_coded = ifelse(alt_is_coded, beta_alt, -beta_alt), p_ext_coded = ifelse(alt_is_coded, p_ext_alt, 1 - p_ext_alt))]
fwrite(R, file.path(OUT, "results/emmax_tested_snps_oriented.tsv"), sep = "\t")
cat("Tested SNPs (PC1/PC2 rows) with ALT == coded allele:", R[pheno %in% c("PC1", "PC2"), sum(alt_is_coded)], "of", R[pheno %in% c("PC1", "PC2"), .N], "\n")

## ---- (c) allele-frequency consistency of the tested SNPs ------------------------------------------------------
## Predicted external coded-allele frequency = p_pol + mean(q) * (p_aqu - p_pol), with mean ancestry q from 09_.
## Only SNPs with |p_aqu - p_pol| >= 0.3 are informative; single-SNP AF at n = 21 is noisy (SD ~0.1), so this is
## supporting evidence only (the LD-phase check above is the strong validation).
qbar <- mean(readRDS(file.path(OUT, "data/orientation_ancestry.rds"))$q)
T <- R[pheno %in% c("PC1", "PC2")]
T[, `:=`(p_aqu = paf[tested_marker, p_aqu], p_pol = paf[tested_marker, p_pol])]
T[, `:=`(dd = p_aqu - p_pol, pred = p_pol + qbar * (p_aqu - p_pol))]
T[, `:=`(err_rule = abs(p_ext_coded - pred), err_flip = abs((1 - p_ext_coded) - pred))]
inf <- T[abs(dd) >= 0.3]
sink(file.path(OUT, "results/orientation_validation.txt"), append = TRUE, split = TRUE)
cat("\n(c) AF consistency of tested SNPs (mean ancestry q =", round(qbar, 3), "): informative SNPs", nrow(inf), "of", nrow(T),
    "; rule closer to predicted AF for", sum(inf$err_rule < inf$err_flip), "; with AF gap > 0.3 between rule and flip:",
    sum(abs(inf$err_flip - inf$err_rule) > 0.3), "(all favour the rule:", all(inf[abs(err_flip - err_rule) > 0.3]$err_rule < inf[abs(err_flip - err_rule) > 0.3]$err_flip), ")\n")
sink()
