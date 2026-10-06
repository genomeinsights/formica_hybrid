## module_external_validation / 12: extend the validation from the single best SNP per outlier region to
## ALL outlier SNPs of the original data.
##
## Outlier SNPs  = every member SNP of every raw-threshold Stage-1 unit of PC1 / PC2 (BF >= 15) in the
##                 original hybrid data (not only the unit's best marker).
## External unit = Stage-1 cluster of the external data (03_, rho = 0.5) containing >= 1 outlier SNP.
## Tested SNP    = the best SNP of that external cluster: the member (restricted to SNPs that also exist in the
##                 original data and have <= 20% missing) with the highest |r| to the cluster consensus
##                 (mean dosage of members sign-aligned to the cluster core), as in LDscnR::eMLG_best_snp.
## Predicted direction of the tested SNP's CODED allele, from the original data:
##    e_j  (outlier SNP j) = k * sign(Beta_is[original unit]) * sign(r_orig(original best marker, j)),  k = -1 (11_)
##    e_best = sign( sum_j |r_orig(best, j)| * e_j * sign(r_orig(best, j)) )  over outlier SNPs j in the cluster with
##             |r_orig(best, j)| >= 0.5; units with no such SNP are dropped.
## Observed: external EMMAX effect of the best SNP on its coded allele (ALT oriented by DIEM Polarity, 10_).
## Dependence: external units inside the same original region (Stage-2 group) are LD-correlated, so inference is
## made at the REGION level (sign-flip permutation of region means of e*z); unit-level fractions are descriptive.
##
## Run from the repo root:  Rscript module_external_validation/R/12_outlier_snps_to_external_units.R
## Writes: results/extunits_tested.tsv, results/extunits_summary.txt, results/extunits.png

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
devtools::load_all("~/gitlab/LDscnR/", quiet = TRUE); source("module_external_validation/R/utils_emmax.R")
OUT <- "module_external_validation"; mod <- "module_manuscript_rho05"; K_ORIG <- -1; R_MIN <- 0.5; MAXMISS <- 0.2
set.seed(1)

## ---- original side ------------------------------------------------------------------------------------------------
S1DIR <- file.path(mod, "baypass_stage1/aland_excluded_S1units"); ord <- readLines(file.path(S1DIR, "S1units_group_order.txt"))
bres <- readRDS(file.path(mod, "data/moduleB_stage1_units_bestsnp.rds")); members <- setNames(bres$groups$members, bres$groups$group_id)
bst <- bres$best$stats
e <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e); G <- e$GTs_hybrids_005
paf <- readRDS(file.path(OUT, "data/orig_parent_af.rds")); setkey(paf, marker)
units <- fread(file.path(OUT, "results/outlier_units_vs_vcf.tsv"))[analysis %in% c("PC1", "PC2")]
s1o <- as.data.table(readRDS("module0_ld_pruning/data/pruned_stage1.rds")$clusters)[n_snps >= 5L]
s2 <- readRDS("module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05.rds"); m2g <- as.data.table(s2$groups)[, .(marker = unlist(members)), by = .(region = group_id)]; setkey(m2g, marker)
u2r <- data.table(group_id = paste0("S1_", s1o$CL_id), region = m2g[.(s1o$core_snp), region]); units <- merge(units, u2r, by = "group_id")
beta <- sapply(c("PC1", "PC2"), function(a) fread(file.path(S1DIR, sprintf("%s_S1units_withOmega_summary_betai_reg.out", a)), select = "Beta_is")[[1]], simplify = FALSE)

osnp <- rbindlist(lapply(seq_len(nrow(units)), function(i) {
  u <- units[i]; mk <- intersect(members[[u$group_id]], colnames(G)); bm <- u$best_marker
  r <- suppressWarnings(cor(G[, bm], G[, mk, drop = FALSE], use = "pairwise.complete.obs"))[1, ]
  data.table(analysis = u$analysis, region = u$region, unit = u$group_id, snp = mk, r_to_unit_best = as.numeric(r),
             beta_unit = beta[[u$analysis]][match(u$group_id, ord)])
}))
osnp[, e_snp := K_ORIG * sign(beta_unit) * sign(r_to_unit_best)]
osnp <- osnp[is.finite(e_snp)]

## ---- external side -------------------------------------------------------------------------------------------------------
d <- readRDS(file.path(OUT, "data/ext_geno_maf10.rds")); Gx <- d$geno; emap <- d$map
s1 <- readRDS(file.path(OUT, "data/ext_stage1_maf10_rho05.rds")); mp <- as.data.table(s1$map_snp); cl <- as.data.table(s1$clusters)
osnp[, ext_cl := mp$CL_id[match(snp, mp$marker)]]
cat_sum <- osnp[, .(outlier_snps = .N, in_external = sum(!is.na(ext_cl)), ext_units = uniqueN(ext_cl, na.rm = TRUE)), by = analysis]

usable <- emap[fmiss21 <= MAXMISS & marker %in% colnames(G), marker]
best_of <- function(cid) {
  mk <- cl$members[[match(cid, cl$CL_id)]]; core <- cl$core_snp[match(cid, cl$CL_id)]
  cand <- intersect(mk, usable); if (!length(cand)) return(NA_character_)
  if (length(mk) == 1L) return(cand[1])
  M <- Gx[, intersect(mk, colnames(Gx)), drop = FALSE]
  ref <- if (core %in% colnames(M)) M[, core] else M[, 1]
  s <- suppressWarnings(sign(cor(M, ref, use = "pairwise.complete.obs"))[, 1]); s[!is.finite(s)] <- 1
  Ma <- M; Ma[, s < 0] <- 2 - M[, s < 0, drop = FALSE]; cons <- rowMeans(Ma, na.rm = TRUE)
  rr <- suppressWarnings(abs(cor(Gx[, cand, drop = FALSE], cons, use = "pairwise.complete.obs"))[, 1]); rr[!is.finite(rr)] <- -1
  cand[which.max(rr)]
}
ucl <- unique(osnp$ext_cl[!is.na(osnp$ext_cl)]); bestsnp <- setNames(vapply(ucl, best_of, ""), ucl)

## ---- predicted direction for each (analysis, external unit) ----------------------------------------------------------
vote <- function(a, cid) {
  b <- bestsnp[[as.character(cid)]]; if (is.na(b)) return(NULL)
  o <- osnp[analysis == a & ext_cl == cid]
  r <- suppressWarnings(cor(G[, b], G[, o$snp, drop = FALSE], use = "pairwise.complete.obs"))[1, ]
  o[, r_best := as.numeric(r)]; o <- o[is.finite(r_best) & abs(r_best) >= R_MIN]; if (!nrow(o)) return(NULL)
  o[, v := abs(r_best) * e_snp * sign(r_best)]
  reg <- o[, .(w = sum(abs(r_best))), by = region][which.max(w), region]
  data.table(analysis = a, ext_unit = cid, tested_marker = b, region = reg, n_outlier_snps = nrow(o), n_orig_regions = uniqueN(o$region),
             e_best = sign(sum(o$v)), vote_agreement = max(mean(sign(o$v) > 0), mean(sign(o$v) < 0)), best_is_outlier = b %in% o$snp)
}
pu <- unique(osnp[!is.na(ext_cl), .(analysis, ext_cl)])
V <- rbindlist(lapply(seq_len(nrow(pu)), function(i) vote(pu$analysis[i], pu$ext_cl[i])))
V <- V[e_best != 0]

## ---- EMMAX on the best SNPs -------------------------------------------------------------------------------------------------
keep <- which(emap$fmiss21 <= MAXMISS); Xs <- meanimp(Gx[, keep, drop = FALSE]); emap2 <- emap[keep]
K <- readRDS(file.path(OUT, "data/ext_grm_pruned.rds"))$K; prep <- emmax_setup(Xs, K)
cl_pc <- fread(file.path(OUT, "data/new_sample_climate_pca_wc.tsv")); cl_pc[, sample := sub("\\.bam$", "", BAM_ID)]; cl_pc <- cl_pc[match(rownames(Xs), sample)]
Y <- list(PC1 = cl_pc$PC1, PC2 = cl_pc$PC2, PC1_rint = rint(cl_pc$PC1), PC2_rint = rint(cl_pc$PC2))
polarity <- paf$Polarity[match(V$tested_marker, paf$marker)]
res <- rbindlist(lapply(names(Y), function(ph) {
  v <- V[analysis == sub("_rint$", "", ph)]; if (!nrow(v)) return(NULL)
  idx <- match(v$tested_marker, colnames(Xs)); stopifnot(!anyNA(idx))
  bs <- beta_se(prep, Y[[ph]], idx); pv <- emmax_fast(prep, Y[[ph]])
  alt_coded <- paf$Polarity[match(v$tested_marker, paf$marker)] == 0
  o <- copy(v)[, `:=`(pheno = ph, beta_coded = ifelse(alt_coded, bs$beta, -bs$beta), se = bs$se, p = pv[idx])]
  o[, `:=`(z = beta_coded / se, concordant = sign(beta_coded) == e_best, maf21 = emap2$maf21[idx])][]
}))
fwrite(res, file.path(OUT, "results/extunits_tested.tsv"), sep = "\t")

## ---- inference at region level ------------------------------------------------------------------------------------------------------
perm_region <- function(dt, B = 50000) {
  R <- dt[, .(T = mean(e_best * z)), by = region]; T0 <- sum(R$T)
  Tb <- replicate(B, sum(sample(c(-1, 1), nrow(R), TRUE) * R$T)); list(T0 = T0, p = (1 + sum(Tb >= T0)) / (B + 1), n_reg = nrow(R), reg_conc = sum(R$T > 0))
}
summ <- res[, { pr <- perm_region(.SD); .(units = .N, regions = pr$n_reg, units_concordant = sum(concordant), frac_units = round(mean(concordant), 3),
                                          regions_concordant = pr$reg_conc, p_regions_binom = binom.test(pr$reg_conc, pr$n_reg, 0.5, alternative = "greater")$p.value,
                                          p_perm_regions = pr$p, units_p05 = sum(p < 0.05), concordant_among_p05 = sum(p < 0.05 & concordant)) }, by = pheno]
sink(file.path(OUT, "results/extunits_summary.txt"), split = TRUE)
cat("Outlier SNPs (members of raw-threshold original units) and their mapping to external Stage-1 units\n"); print(cat_sum)
cat("\nExternal units with a usable best SNP and a predicted direction (|r_orig| >= ", R_MIN, "):\n", sep = ""); print(V[, .(units = .N, regions = uniqueN(region), median_outlier_snps = as.numeric(median(n_outlier_snps)),
    best_is_outlier = sum(best_is_outlier), median_vote_agreement = round(median(vote_agreement), 2), multi_region_units = sum(n_orig_regions > 1)), by = analysis])
cat("\nConcordance of external effect direction with the original association (all outlier SNPs -> external units)\n")
print(summ[, `:=`(p_regions_binom = signif(p_regions_binom, 3), p_perm_regions = signif(p_perm_regions, 3))])
cat("\nBy whether the tested SNP is itself an outlier SNP (raw PC phenotypes):\n"); print(res[pheno %in% c("PC1", "PC2"), .(units = .N, concordant = sum(concordant), frac = round(mean(concordant), 3)), by = .(pheno, best_is_outlier)][order(pheno)])
cat("\nBy vote agreement >= 0.9 (raw PC phenotypes):\n"); print(res[pheno %in% c("PC1", "PC2"), .(units = .N, concordant = sum(concordant), frac = round(mean(concordant), 3)), by = .(pheno, high_agree = vote_agreement >= 0.9)][order(pheno)])
sink()
g <- ggplot(res[pheno %in% c("PC1", "PC2")], aes(e_best * z, fill = concordant)) + geom_histogram(bins = 25) + facet_wrap(~pheno) + geom_vline(xintercept = 0) +
  scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#D55E00")) + labs(x = "external z x predicted sign (>0 concordant)", y = "external units") + theme_classic()
ggsave(file.path(OUT, "results/extunits.png"), g, width = 9, height = 4, dpi = 200)
