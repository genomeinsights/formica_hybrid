## =========================================================================
## DI25 among-region -- pairwise relationships among the four LD estimates.
## =========================================================================
## Scatter-matrix of the per-pair estimates from di25D_emmax.R:
##   Rsq          = EMMAX structure-corrected r^2 (double-LOCO GRM + 10 genome PCs)
##   r2_raw       = plain pairwise r^2 (individuals unrelated, no correction)
##   R_st         = among-population correlation (signed; the D'2st / structure axis)
##   within_pop_r = median signed within-population correlation (structure removed by strata)
## Plotted on a SUBSET so it renders (the full unlinked set is millions of pairs): a random sample PLUS the top tail
## of each estimate PLUS all FDR edges (orange), so no extreme is lost to subsampling.
## Upper-panel correlations are on the FULL data (unbiased), annotated in the title.
## Reads di25D_emmax.rds. Writes Figures/di25D_estimate_pairs.png.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(GGally) })
e <- readRDS("module_di25/among_region_association/data/di25D_emmax.rds")
p <- e$pairs[unlinked == TRUE & is.finite(within_pop_r)]
## ABSOLUTE values: put all four estimates on one |correlation| scale [0,1] for a clean
## magnitude comparison (EMMAX/raw as |r| = sqrt(r^2); R_st / within_pop_r as |r|).
p[, `:=`(aRsq = sqrt(Rsq), araw = sqrt(r2_raw), aRst = abs(R_st), awpr = abs(within_pop_r))]
vars <- c("aRsq", "araw", "aRst", "awpr")
p[, pk := paste(pmin(focal, partner), pmax(focal, partner))]
ek <- paste(pmin(e$edges$cluster1, e$edges$cluster2), pmax(e$edges$cluster1, e$edges$cluster2))

## subset: random 4000 + top 400 |.| of each estimate + all FDR edges
set.seed(1)
rnd   <- p[sample(.N, 4000)]
tails <- rbindlist(lapply(vars, function(v) p[order(-abs(get(v)))][1:400]))
fdr   <- p[pk %in% ek]
sub   <- unique(rbind(rnd, tails, fdr), by = "pk")
sub[, grp := ifelse(pk %in% ek, "FDR edge", "other")]
message(sprintf("[pairs] %d subset points (%d random + tails + %d FDR edges) of %s",
                nrow(sub), nrow(rnd), nrow(fdr), format(nrow(p), big.mark = ",")))

lower_fn <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_point(data = data[data$grp == "other", ], colour = "grey65", size = 0.35, alpha = 0.28) +
    geom_point(data = data[data$grp == "FDR edge", ], colour = "#D55E00", size = 1.9) +
    theme_bw(base_size = 8)
}
diag_fn <- function(data, mapping, ...)
  ggplot(data, mapping) + geom_density(fill = "grey80", linewidth = 0.3) + theme_bw(base_size = 8)

g <- ggpairs(sub, columns = vars, mapping = aes(),
             lower = list(continuous = lower_fn), diag = list(continuous = diag_fn),
             upper = list(continuous = wrap("cor", size = 2.6, use = "complete.obs")),
             columnLabels = c("|r| EMMAX", "|r| raw", "|R_st|", "|within-pop r|")) +
  theme(strip.text = element_text(size = 7))
ttl <- sprintf(paste0("DI25 pairwise LD estimates, |r| (subset: 4k random + top-400 tails + FDR edges [orange]; q<%.2g)\n",
                      "FULL-data:  cor(|r|raw, |R_st|)=%.2f   cor(|r|EMMAX, |R_st|)=%.2f   cor(|r|EMMAX, |within|)=%.2f   cor(|r|EMMAX, |r|raw)=%.2f"),
               e$params$Q_FDR, cor(p$araw, p$aRst), cor(p$aRsq, p$aRst), cor(p$aRsq, p$awpr), cor(p$aRsq, p$araw))

png("module_di25/among_region_association/Figures/di25D_estimate_pairs.png",
    width = 1700, height = 1700, res = 200)
print(g + labs(title = ttl) + theme(plot.title = element_text(size = 8.5)))
dev.off()
cat("[done] Figures/di25D_estimate_pairs.png\n")
