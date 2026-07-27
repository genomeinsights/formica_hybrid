## =========================================================
## MODULE A -- follow-up: DI vs unidirectional (aqu/pol) asymmetry
## =========================================================
## Question raised from the A1 result: the MAF-gated per-SNP sorting is
## POLYCTENA-leaning (17.3% vs 14.3% aquilonia), a flip from both earlier
## threads' ~1.8-2:1 AQUILONIA bias -- but those were measured on a DI-gated
## set. If the DIRECTION of unidirectional sorting depends on DI, the two
## results are just opposite ends of the DI axis, not a contradiction.
##
## This script tests that directly, as a pure post-hoc on data/moduleA_snp.rds
## (per-SNP parallelism_stats output; no genotype matrices, no re-run). It
## reads only that file, so it can run independently of / after the main
## Module A run. "Unidirectional asymmetry" = the signed aquilonia-vs-polyctena
## lean, measured two equivalent ways:
##   * uni_score = (n_aqu - n_pol)/n_obs  per locus (signed; + = aquilonia)
##   * among UNIDIRECTIONAL loci, the fraction classed "aquilonia"
##
## We always work on the differentiated (parent_maf-gated) set, and -- because
## parent_maf is a known confounder of sorting -- also control for it, so the
## DI effect is not just MAF leaking through.

suppressPackageStartupMessages(library(data.table))

snp <- readRDS("data/moduleA_snp.rds")
d   <- snp[differentiated == TRUE & !is.na(sort_class) & !is.na(DI)]
cat("Differentiated, classified loci:", nrow(d), "\n\n")

## ---- 1. headline: does DI predict directional lean? --------------------
## Spearman over ALL classified loci (uni_score is signed lean), and a
## MAF-controlled locus-level view.
cat("Spearman cor(DI, uni_score) [+ = aquilonia lean]:",
    round(cor(d$DI, d$uni_score, method = "spearman", use = "complete.obs"), 3), "\n")
cat("Spearman cor(DI, parent_maf):",
    round(cor(d$DI, d$parent_maf, method = "spearman", use = "complete.obs"), 3), "\n\n")

## ---- 2. DI-decile table ------------------------------------------------
## Per DI decile: composition + the two asymmetry measures + median MAF (so we
## can see whether MAF co-varies with DI across the deciles).
d[, DI_decile := cut(DI, breaks = quantile(DI, 0:10/10, na.rm = TRUE),
                     include.lowest = TRUE, labels = FALSE)]
uni <- c("aquilonia", "polyctena")
di_tab <- d[, .(
  n            = .N,
  DI_lo        = round(min(DI), 1),
  DI_hi        = round(max(DI), 1),
  pct_aqu      = round(100 * mean(sort_class == "aquilonia"), 1),
  pct_pol      = round(100 * mean(sort_class == "polyctena"), 1),
  pct_unsorted = round(100 * mean(sort_class == "unsorted"), 1),
  frac_aqu_of_unidir = round(sum(sort_class == "aquilonia") /
                               max(1L, sum(sort_class %in% uni)), 3),
  mean_uni_score     = round(mean(uni_score, na.rm = TRUE), 3),
  median_parent_maf  = round(median(parent_maf, na.rm = TRUE), 3)
), by = DI_decile][order(DI_decile)]
cat("--- DI decile x directional asymmetry ---\n")
print(di_tab)

## ---- 3. DI-threshold sweep (cumulative DI >= t) ------------------------
## Mirrors "gating at DI >= t": how the aqu/pol balance among unidirectional
## loci moves as the DI floor rises -- the direct bridge to the earlier
## DI-gated threads (they used high-ish DI floors).
thr <- quantile(d$DI, c(.5, .6, .7, .8, .9, .95, .99), na.rm = TRUE)
sweep <- rbindlist(lapply(seq_along(thr), function(i) {
  s <- d[DI >= thr[i]]
  data.table(
    DI_threshold       = round(thr[i], 1),
    quantile           = names(thr)[i],
    n                  = nrow(s),
    n_unidir           = s[sort_class %in% uni, .N],
    frac_aqu_of_unidir = round(s[sort_class %in% uni, mean(sort_class == "aquilonia")], 3),
    mean_uni_score     = round(mean(s$uni_score, na.rm = TRUE), 3)
  )
}))
cat("\n--- DI threshold sweep (loci with DI >= t) ---\n")
print(sweep)

## ---- 4. MAF-controlled test: P(aquilonia | unidirectional) ~ DI + MAF ---
u <- d[sort_class %in% uni]
u[, is_aqu := as.integer(sort_class == "aquilonia")]
fit <- glm(is_aqu ~ DI + parent_maf, data = u, family = binomial())
cat("\n--- logistic P(aquilonia direction) ~ DI + parent_maf (unidirectional loci) ---\n")
print(round(summary(fit)$coefficients, 4))
cat("\n(interpretation: a positive, significant DI coefficient => more-diagnostic\n",
    " loci lean aquilonia even after controlling for parental MAF -- i.e. the\n",
    " polyctena tilt is concentrated at LOW DI, reconciling the thread flip.)\n")

## ========================================================================
## 5. CLUSTER level: direction vs cluster SIZE (and DI), from A3
## ========================================================================
## A3 showed the aqu/pol lean FLIPS with cluster size: small clusters lean
## polyctena, large (>=5-loci) blocks lean aquilonia. Cluster size tracks low
## recombination / high LD, so this is the architecture x direction signal --
## and the likeliest reconciliation of the thread flip (the earlier aquilonia
## bias lived in the big low-recombination blocks). Uses the cluster-level A3
## output (one row per independent unit), which carries n_loci, DI (cluster
## max), parent_maf (cluster), uni_score and sort_class.
cl <- readRDS("data/moduleA_clusters.rds")
c2 <- cl[differentiated == TRUE & !is.na(sort_class)]
cat("\nCluster-level: differentiated, classified units:", nrow(c2), "\n")

## size-class x direction
c2[, size_class := cut(n_loci, breaks = c(0, 1, 4, 9, 49, 199, Inf),
                       labels = c("1", "2-4", "5-9", "10-49", "50-199", "200+"))]
size_tab <- c2[, .(
  n                  = .N,
  pct_unsorted       = round(100 * mean(sort_class == "unsorted"), 1),
  frac_aqu_of_unidir = round(sum(sort_class == "aquilonia") /
                               max(1L, sum(sort_class %in% uni)), 3),
  mean_uni_score     = round(mean(uni_score, na.rm = TRUE), 3),
  median_DI          = round(median(DI, na.rm = TRUE), 1),
  median_parent_maf  = round(median(parent_maf, na.rm = TRUE), 3)
), by = size_class][order(size_class)]
cat("\n--- CLUSTER size-class x directional asymmetry (A3 units) ---\n")
print(size_tab)

## does size predict direction AFTER controlling DI + MAF?
cu <- c2[sort_class %in% uni]
cu[, is_aqu     := as.integer(sort_class == "aquilonia")]
cu[, log2_nloci := log2(n_loci)]
fit_sz <- glm(is_aqu ~ log2_nloci + DI + parent_maf, data = cu, family = binomial())
cat("\n--- logistic P(aquilonia) ~ log2(n_loci) + DI + parent_maf (unidirectional CLUSTERS) ---\n")
print(round(summary(fit_sz)$coefficients, 4))
cat("\n(a positive log2_nloci coefficient => larger clusters lean aquilonia even\n",
    " after DI + MAF, i.e. the size flip is not just DI or MAF in disguise;\n",
    " compare the DI and parent_maf coefficients for their independent pulls.)\n")

saveRDS(list(di_decile = di_tab, di_sweep = sweep, glm_snp = fit,
             size_tab = size_tab, glm_cluster = fit_sz),
        "data/moduleA_di_asymmetry.rds")
cat("\nSaved: data/moduleA_di_asymmetry.rds\n")

## ========================================================================
## 6. Figure 1 (a-c) -- publication panel
## ========================================================================
## Vector PDF (Overleaf/LaTeX-ready) + PNG preview. Categorical colours are the
## Okabe-Ito colourblind-safe set (fixed order); "unsorted" is a deliberate
## neutral grey. One series per panel where the axis already names the entity.
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("patchwork", quietly = TRUE)) {
  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
  dir.create("Figures", showWarnings = FALSE)

  cls_col <- c(aquilonia = "#0072B2", polyctena = "#D55E00",
               unsorted  = "#999999", bidirectional = "#CC79A7")
  th <- theme_classic(base_size = 9) +
    theme(plot.tag = element_text(face = "bold", size = 11),
          axis.title = element_text(size = 8))

  ## (a) genome-wide sort-class proportions among parent-polymorphic loci
  a_dt <- d[!is.na(sort_class), .N, by = sort_class]
  a_dt[, pct := 100 * N / sum(N)]
  a_dt[, sort_class := factor(sort_class,
        levels = c("unsorted", "aquilonia", "polyctena", "bidirectional"))]
  p1a <- ggplot(a_dt, aes(sort_class, pct, fill = sort_class)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.35, size = 2.6) +
    scale_fill_manual(values = cls_col, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(x = NULL, y = "% of parent-polymorphic loci") +
    th + theme(axis.text.x = element_text(angle = 20, hjust = 1))

  ## (b) direction reverses across the diagnostic-index axis
  b_dt <- copy(di_tab)
  b_dt[, DI_mid := (DI_lo + DI_hi) / 2]
  b_dt[, lean := ifelse(frac_aqu_of_unidir >= 0.5, "aquilonia", "polyctena")]
  p1b <- ggplot(b_dt, aes(DI_mid, frac_aqu_of_unidir)) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
    geom_line(colour = "grey45", linewidth = 0.5) +
    geom_point(aes(colour = lean), size = 2.3) +
    scale_colour_manual(values = cls_col, guide = "none") +
    scale_y_continuous(limits = c(0, 1)) +
    annotate("text", x = max(b_dt$DI_mid), y = 0.97, label = "toward aquilonia",
             hjust = 1, size = 2.5, colour = cls_col[["aquilonia"]]) +
    annotate("text", x = max(b_dt$DI_mid), y = 0.03, label = "toward polyctena",
             hjust = 1, size = 2.5, colour = cls_col[["polyctena"]]) +
    labs(x = "Diagnostic index (decile midpoint)",
         y = "Fraction fixing toward aquilonia") + th

  ## (c) sorting collapses in large linkage blocks
  c_dt <- copy(size_tab)
  p1c <- ggplot(c_dt, aes(size_class, pct_unsorted, group = 1)) +
    geom_line(colour = "grey45", linewidth = 0.5) +
    geom_point(size = 2.3, colour = "grey20") +
    geom_text(aes(label = sprintf("%.0f", pct_unsorted)), vjust = -0.8, size = 2.4) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Linkage-block size (n loci)", y = "% of units unsorted") +
    th + theme(axis.text.x = element_text(angle = 25, hjust = 1))

  fig1 <- p1a + p1b + p1c + plot_annotation(tag_levels = "a")
  ggsave("Figures/moduleA_fig1.pdf", fig1, width = 180, height = 70, units = "mm")
  ggsave("Figures/moduleA_fig1.png", fig1, width = 180, height = 70, units = "mm", dpi = 300)
  cat("Saved: Figures/moduleA_fig1.pdf / .png\n")
} else {
  cat("ggplot2/patchwork not available -- skipped Figure 1.\n")
}
