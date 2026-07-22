## ==========================================================================
## Member-level robustness check for the sorting x BayPass-outlier overlap.
## Instead of one representative tag SNP per cluster, score EVERY member SNP
## and aggregate to the fraction of (differentiated) members that are
## directionally sorted. Re-test outlier vs background, size-adjusted.
## ==========================================================================

suppressMessages({library(data.table)})
source("dev/R/Ohta.R"); source("dev/R/parallelism_stats.R"); set.seed(1)

e <- new.env(); load("./data/hybrids_and_parents_maf005.Rdata", envir=e)
GTs <- e$GTs_with_parents; sd <- e$sample_data_with_parents; map <- e$map_hyb_005
DI_by <- setNames(map$DiagnosticIndex, map$marker)
aqu <- "aquilonia_parent"; pol <- "polyctena_parent"
hyb <- setdiff(unique(sd$Population), c(aqu, pol))
mo <- new.env(); load("./data/hybrids_only_maf005.Rdata", envir=mo); map_bf <- mo$map_hyb_005

g <- readRDS("data/eMLG_5loci_0025_cM05.rds")$groups[has_eMLG == TRUE]
marker_group <- g[, .(marker = unlist(members), n_loci = n_loci), by = group_id]
cat("has_eMLG clusters:", nrow(g), "| member SNPs:", nrow(marker_group), "\n")

## ---- BayPass outlier clusters (union over 8 configs), as before ----
sig_threshold <- 20; min_n_sig <- 10
configs <- CJ(ps=c("with_aland","aland_excluded"), pc=c("PC1","PC2"), om=c("noOmega","withOmega"))
import_bf <- function(f){ r<-fread(f); stopifnot(identical(r$MRK,seq_len(nrow(r)))); setNames(r$`BF(dB)`, map_bf$marker) }
oid <- lapply(seq_len(nrow(configs)), function(i){ cf<-configs[i]
  bf<-import_bf(sprintf("./%s/%s_DIEM_%s_summary_betai_reg.out",cf$ps,cf$pc,cf$om))
  mg<-copy(marker_group); mg[,BF:=bf[marker]]
  mg[,.(n_sig=sum(BF>=sig_threshold,na.rm=TRUE)),by=group_id][n_sig>=min_n_sig,group_id] })
names(oid) <- configs[, paste(ps,pc,om,sep="_")]
outlier_any <- unique(unlist(oid))
outlier_PC2 <- unique(unlist(oid[grep("PC2",names(oid))]))
cat("outlier clusters (any):", length(outlier_any), "| PC2:", length(outlier_PC2), "\n")

## ---- per-MEMBER sort_class over all member SNPs ----
mk <- unique(marker_group$marker); mk <- intersect(mk, colnames(GTs))
prep <- ohta_fast_prepare(GTs[, mk], pops = sd$Population)
ps <- parallelism_stats(prep, hyb, aqu, pol, DI = DI_by,
                        min_parent_diff = 0.5, sort_th = 0.5, fix_th = 0.1)
ps[, `:=`(m_diff = differentiated,
          m_dir  = sort_class %in% c("aquilonia","polyctena"),
          m_aqu  = sort_class == "aquilonia",
          m_pol  = sort_class == "polyctena")]

mg <- ps[, .(marker, m_diff, m_dir, m_aqu, m_pol, uni_score)][marker_group, on = "marker"]

## ---- aggregate to cluster: fraction of DIFFERENTIATED members directional ----
cl <- mg[, .(
    n_loci      = n_loci[1],
    n_diff      = sum(m_diff, na.rm = TRUE),
    frac_dir    = mean(m_dir[m_diff], na.rm = TRUE),      # among differentiated members
    frac_aqu    = mean(m_aqu[m_diff], na.rm = TRUE),
    frac_pol    = mean(m_pol[m_diff], na.rm = TRUE),
    mean_uni    = mean(uni_score[m_diff], na.rm = TRUE)
  ), by = group_id]
cl[, `:=`(is_outlier   = as.integer(group_id %in% outlier_any),
          is_outlierP2 = as.integer(group_id %in% outlier_PC2))]

## clusters with enough differentiated members to estimate a fraction
d <- cl[n_diff >= 5]
d[, cl_dir := as.integer(frac_dir >= 0.5)]              # majority of members directional
cat(sprintf("\nclusters with >=5 differentiated members: %d | BayPass outliers among them: %d (PC2: %d)\n",
            nrow(d), d[, sum(is_outlier)], d[, sum(is_outlierP2)]))

cat("\n=== member-level directional fraction: outlier vs background ===\n")
print(d[, .(n=.N, med_nloci=as.double(median(n_loci)), med_ndiff=as.double(median(n_diff)),
            mean_frac_dir=round(mean(frac_dir),3),
            median_frac_dir=round(median(frac_dir),3),
            pct_majority_dir=round(100*mean(cl_dir),1)),
        by=.(outlier=ifelse(is_outlier==1,"BayPass outlier","background"))])

## ---- size-adjusted tests ----
cat("\n=== size-adjusted: continuous frac_dir ~ outlier + log(n_loci) ===\n")
m1 <- lm(frac_dir ~ is_outlier + log(n_loci), data = d)
print(round(coef(summary(m1)), 4))
cat("\n=== size-adjusted: majority-directional (logistic) ~ outlier + log(n_loci) ===\n")
m2 <- glm(cl_dir ~ is_outlier + log(n_loci), data = d, family = binomial)
ci <- confint.default(m2, "is_outlier")
cat(sprintf("  OR = %.2f  (95%% CI %.2f-%.2f)  p = %.2g\n",
            exp(coef(m2)["is_outlier"]), exp(ci[1]), exp(ci[2]), summary(m2)$coefficients["is_outlier","Pr(>|z|)"]))
cat("\n=== PC2-only, majority-directional logistic ===\n")
m3 <- glm(cl_dir ~ is_outlierP2 + log(n_loci), data = d, family = binomial)
ci3 <- confint.default(m3, "is_outlierP2")
cat(sprintf("  OR = %.2f  (95%% CI %.2f-%.2f)  p = %.2g\n",
            exp(coef(m3)["is_outlierP2"]), exp(ci3[1]), exp(ci3[2]), summary(m3)$coefficients["is_outlierP2","Pr(>|z|)"]))

## ---- the outlier clusters themselves: per-cluster member consensus ----
cat("\n=== each BayPass-outlier cluster (>=5 diff members): member-level sorting ===\n")
oc <- d[is_outlier == 1][order(-frac_dir)]
oc[, direction := fifelse(frac_aqu > frac_pol, "aquilonia", fifelse(frac_pol > frac_aqu, "polyctena", "mixed"))]
print(oc[, .(group_id, n_loci, n_diff, frac_dir=round(frac_dir,2),
             frac_aqu=round(frac_aqu,2), frac_pol=round(frac_pol,2), direction, PC2=is_outlierP2)])
