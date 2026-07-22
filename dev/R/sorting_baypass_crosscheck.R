## ==========================================================================
## Cross-check: do PREDICTABLY SORTED eMLG clusters co-localise with
## BayPass PC1/PC2 association (climate-adaptation) outlier clusters?
## ==========================================================================
## Sorting (intrinsic dynamics) vs BayPass outliers (extrinsic covariate
## association) -- are they the same regions? Outlier clusters are
## mechanically large (require >=10 significant members), so we control for
## cluster size, exactly as R/diagnostic_index_enrichment.R does.

suppressMessages({library(data.table)})
source("dev/R/Ohta.R"); source("dev/R/parallelism_stats.R"); set.seed(1)

## ---- data ----
e <- new.env(); load("./data/hybrids_and_parents_maf005.Rdata", envir=e)
GTs <- e$GTs_with_parents; sd <- e$sample_data_with_parents; map <- copy(e$map_hyb_005)
DI_by <- setNames(map$DiagnosticIndex, map$marker)
aqu <- "aquilonia_parent"; pol <- "polyctena_parent"
hyb <- setdiff(unique(sd$Population), c(aqu, pol))

## BF alignment uses the map from hybrids_only (identical marker order as the
## BayPass inputs); we key BF by marker name to be robust to ordering.
mo <- new.env(); load("./data/hybrids_only_maf005.Rdata", envir=mo)
map_bf <- mo$map_hyb_005

g_all <- readRDS("data/eMLG_5loci_0025_cM05.rds")$groups
g <- g_all[has_eMLG == TRUE]                    # n_loci >= 5 clusters
marker_group <- g[, .(marker = unlist(members), n_loci = n_loci), by = group_id]

sig_threshold <- 20; min_n_sig <- 10

import_bf <- function(file) {
  res <- fread(file); stopifnot(identical(res$MRK, seq_len(nrow(res))))
  setNames(res$`BF(dB)`, map_bf$marker)
}

## ---- BayPass outlier clusters per config, and union ----
configs <- CJ(ps=c("with_aland","aland_excluded"), pc=c("PC1","PC2"), om=c("noOmega","withOmega"))
outlier_ids <- list()
for (i in seq_len(nrow(configs))) {
  cf <- configs[i]
  bf <- import_bf(sprintf("./%s/%s_DIEM_%s_summary_betai_reg.out", cf$ps, cf$pc, cf$om))
  mg <- copy(marker_group); mg[, BF := bf[marker]]
  ids <- mg[, .(n_sig = sum(BF >= sig_threshold, na.rm=TRUE)), by=group_id][n_sig >= min_n_sig, group_id]
  outlier_ids[[paste(cf$ps,cf$pc,cf$om,sep="_")]] <- ids
}
outlier_any  <- unique(unlist(outlier_ids))
outlier_PC1  <- unique(unlist(outlier_ids[grep("PC1", names(outlier_ids))]))
outlier_PC2  <- unique(unlist(outlier_ids[grep("PC2", names(outlier_ids))]))
cat("outlier clusters: any config =", length(outlier_any),
    "| PC1 =", length(outlier_PC1), "| PC2 =", length(outlier_PC2), "\n")

## ---- sort_class per cluster (from the representative tag SNP) ----
reps <- intersect(g$representative, colnames(GTs))
prep <- ohta_fast_prepare(GTs[, reps], pops = sd$Population)
ps <- parallelism_stats(prep, hyb, aqu, pol, DI = DI_by,
                        min_parent_diff = 0.5, sort_th = 0.5, fix_th = 0.1)

cl <- g[, .(group_id, marker = representative, n_loci)]
cl <- ps[, .(marker, differentiated, sort_class, prop_fixed, uni_score, bi_score)][cl, on = "marker"]
cl[, `:=`(is_outlier   = as.integer(group_id %in% outlier_any),
          is_outlierP1 = as.integer(group_id %in% outlier_PC1),
          is_outlierP2 = as.integer(group_id %in% outlier_PC2),
          sorted       = as.integer(sort_class %in% c("aquilonia","polyctena","bidirectional")),
          directional  = as.integer(sort_class %in% c("aquilonia","polyctena")))]

d <- cl[differentiated == TRUE]                 # clusters testable for sorting
cat(sprintf("\ntestable (differentiated) eMLG clusters: %d | of which BayPass outliers: %d\n",
            nrow(d), d[, sum(is_outlier)]))

## ---- how sorted are outlier vs non-outlier clusters? ----
cat("\n=== raw rates (differentiated clusters) ===\n")
print(d[, .(n=.N, med_nloci=as.double(median(n_loci)),
            pct_sorted=round(100*mean(sorted),1),
            pct_directional=round(100*mean(directional),1),
            pct_bidir=round(100*mean(sort_class=="bidirectional"),2),
            mean_prop_fixed=round(mean(prop_fixed,na.rm=TRUE),3)),
        by=.(outlier = ifelse(is_outlier==1,"BayPass outlier","background"))])

## ---- size-adjusted tests (logistic, control log(n_loci)) ----
adj <- function(y, x) {
  f <- as.formula(paste(y, "~", x, "+ log(n_loci)"))
  m <- glm(f, data = d, family = binomial); cf <- summary(m)$coefficients
  ci <- confint.default(m, x)
  data.table(response=y, predictor=x, OR=exp(cf[x,"Estimate"]),
             OR_lo=exp(ci[1]), OR_hi=exp(ci[2]), p=cf[x,"Pr(>|z|)"])
}
cat("\n=== size-adjusted odds ratios (outlier vs background) ===\n")
res <- rbindlist(list(
  adj("sorted","is_outlier"), adj("directional","is_outlier"),
  adj("directional","is_outlierP1"), adj("directional","is_outlierP2")))
res[, `:=`(OR=round(OR,2), OR_lo=round(OR_lo,2), OR_hi=round(OR_hi,2), p=signif(p,2))]
print(res)

## ---- also: naive (unadjusted) directional overlap, for contrast ----
ft <- fisher.test(table(d$is_outlier, d$directional))
cat(sprintf("\nnaive directional overlap (unadjusted) OR = %.2f, p = %.2g\n", ft$estimate, ft$p.value))

## ---- direction of sorting among outlier clusters (aqu vs pol vs bi) ----
cat("\n=== sort_class of BayPass-outlier clusters ===\n")
print(d[is_outlier==1, .N, by=sort_class][order(-N)])
