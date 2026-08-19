## =========================================================
## MODULE A (clean _PK) -- genomic architecture of ancestry sorting
## =========================================================
## Reproduces every result and figure in
##   manuscript_notes/moduleB_results_summary.pdf
## (that draft is "Module B"; in the _PK pipeline it is Module A). Built from
## RDS_data/ inputs at the clean _PK conventions: pooled-parental MAF gate 0.15,
## sort_th 0.5, fix_th 0.1 (major-allele floor 0.90), cM05 clustering, DI ungated.
## These settings differ from the current PDF (fix_th 0.15 / cM05), so the printed
## numbers shift slightly -- each cat() prints the regenerated value next to the
## PDF value so the manuscript text can be updated. Consolidates the former
## R/moduleB_architecture.R and the B2/B3 in R_PK/moduleB_architecture.R.
##
## ------------------------------------------------------------------------
## MAP -- PDF part  ->  code section  (grep the tag to jump)
##   [PDF S1] "DI is separable from recombination ..."      tag: == [PDF S1]
##            -> rho(DI,recomb) -0.03, rho(DI,MAF) -0.02,
##               rho(DI,pi) -0.46, rho(DI,dxy) +0.13; Fig 1a
##   [PDF S2] "Sorting is depleted ... in low recombination" tag: == [PDF S2]
##            -> unit 0.04 vs ~0.32, SNP 0.14; magnitude slope +0.05 (z=117); Fig 1b
##   [PDF S3] "Direction is governed by diagnostic index"    tag: == [PDF S3]
##            -> DI +1.57 (z=112), recomb -0.09 (z=-6.5), ~17x weaker; Fig 1c
##   [SUPP]   direction coefficients vs sorting threshold    tag: == [SUPP S3]
##            -> moduleA_supp_th_sensitivity.tex Fig S3
##   assembled 3-panel figure (Fig 1 a/b/c)                  tag: == [FIG 1]
## ------------------------------------------------------------------------
## Reads : data/hybrids_and_parents_maf005.Rdata, module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds,
##         data/Frufa_DTOL_PR.ref_genome.recmap
## Writes: moduleA_sorting/data/moduleA_architecture.rds, moduleA_sorting/data/moduleA_r_{unit,snp}.rds,
##         moduleA_sorting/Figures/moduleA_architecture_fig.{pdf,png}, moduleA_sorting/Figures/moduleA_direction_sweep.png
## Run from the repo root. DESCRIPTIVE only: any recombination effect awaits the
## recombination-matched neutral null (Module E).

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(wesanderson)
})
devtools::load_all("~/gitlab/LDscnR/")            # ohta_fast_prepare() etc. (LDscnR package)
source("moduleA_sorting/R/parallelism_stats.R")   # parallelism_stats()
set.seed(1)

## ---- PARAMETERS (clean _PK conventions) --------------------------------
MIN_PARENT_MAF <- 0.15
FIX_TH         <- 0.15                         # major-allele fixation floor 0.85
SORT_TH        <- 0.6
SORT_RULE      <- "binom"                       # prop_fixed magnitude gate + binomial direction test (see parallelism_stats)
ALPHA          <- 0.05                          # binom direction-test significance (matches the rest of Module A)
SNP_SAMPLE     <- 200000L
SORT_TH_SWEEP  <- c(0.5, 0.6, 0.7, 0.8)        # for the [SUPP S3] direction sweep
CLUSTERING     <- "module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds"
BESTSNP        <- "module0_ld_pruning/data/eMLG_5loci_0025_cM05_bestsnp.rds"
RECMAP         <- "data/Frufa_DTOL_PR.ref_genome.recmap"
RECOMPUTE      <- FALSE

elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))
cache_rds <- function(path, build, label = basename(path), valid = NULL) {
  if (!RECOMPUTE && file.exists(path)) {
    obj <- readRDS(path)
    if (is.null(valid) || isTRUE(valid(obj))) {
      cat(sprintf("      [cache] read  %s\n", path)); return(obj)
    }
    cat(sprintf("      [cache] %s present but stale -> rebuilding\n", path))
  }
  cat(sprintf("      [build] %s ...\n", label))
  obj <- build(); saveRDS(obj, path)
  cat(sprintf("      [cache] wrote %s\n", path)); obj
}

## ---- inputs & per-marker architecture table ----------------------------
e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
GTs <- e$GTs_with_parents; smeta <- e$sample_data_with_parents; map <- copy(e$map_hyb_005)

## per-marker recombination rate (cM/Mb), map-interpolated (exogenous, not size)
rec <- fread(RECMAP); stopifnot(ncol(rec) >= 4)
setnames(rec, 1:4, c("chr", "pos", "cM", "cMMb")); rec[, Chr := sub("chromosome_", "Chr", chr)]
map[, recomb := NA_real_]
for (ch in unique(map$Chr)) {
  r <- rec[Chr == ch]; if (nrow(r) < 2) next
  idx <- map[, which(Chr == ch)]
  map[idx, recomb := approx(r$pos, r$cMMb, xout = map$Pos[idx], rule = 2)$y]
}

## per-marker cluster size (canonical cM05 clustering)
g <- readRDS(CLUSTERING)$groups
## canonical best-SNP: one real marker per cluster (best_marker for has_eMLG,
## else the centrality representative for < min_n_loci_eMLG clusters). The
## LD-reduced-unit architecture is built on this real SNP, not the consensus.
rep_snp_all <- as.data.table(readRDS(BESTSNP)$rep_snp_all)
memb <- data.table(marker = unlist(g$members), cluster_size = rep(g$n_loci, lengths(g$members)))
map <- memb[map, on = "marker"]

## parental frequencies -> pi (within), dxy (absolute), Fst (relative), pooled MAF
aqu <- "aquilonia_parent"; pol <- "polyctena_parent"
aqu_rows <- which(smeta$Population == aqu); pol_rows <- which(smeta$Population == pol)
pa <- (colMeans(GTs[aqu_rows, , drop = FALSE], na.rm = TRUE) / 2)[map$marker]
pp <- (colMeans(GTs[pol_rows, , drop = FALSE], na.rm = TRUE) / 2)[map$marker]
par_pool <- (colMeans(GTs[c(aqu_rows, pol_rows), , drop = FALSE], na.rm = TRUE) / 2)[map$marker]
Hs_a <- 2 * pa * (1 - pa); Hs_p <- 2 * pp * (1 - pp)
map[, pi_within := (Hs_a + Hs_p) / 2]
map[, dxy       := pa * (1 - pp) + pp * (1 - pa)]
pbar <- (pa + pp) / 2; Ht <- 2 * pbar * (1 - pbar)
map[, Fst        := ifelse(Ht > 0, (Ht - (Hs_a + Hs_p) / 2) / Ht, NA_real_)]
map[, parent_maf := pmin(par_pool, 1 - par_pool)]

DI_by     <- setNames(map$DiagnosticIndex, map$marker)
maf_by    <- setNames(map$parent_maf,      map$marker)
recomb_by <- setNames(map$recomb,          map$marker)
cs_by     <- setNames(rep_snp_all$n_loci,  rep_snp_all$rep_snp)   # n_loci keyed by best-SNP marker
hyb       <- setdiff(unique(smeta$Population), c(aqu, pol))

## ========================================================================
## == [PDF S1] "Diagnostic index is separable from recombination and reflects
##             reduced diversity."
##    feeds: rho(DI,recomb), rho(DI,MAF), rho(DI,pi), rho(DI,dxy); Fig 1a
## ========================================================================
d1 <- map[is.finite(recomb) & is.finite(DiagnosticIndex) & !is.na(cluster_size)]
sp <- function(a, b) round(cor(a, b, method = "spearman", use = "complete.obs"), 3)
cat("\n[PDF S1] Spearman correlations (n =", nrow(d1), "markers)  [regenerated | PDF]:\n")
cat("  DI vs recombination :", sp(d1$DiagnosticIndex, d1$recomb),      "| -0.03\n")
cat("  DI vs parental MAF  :", sp(d1$DiagnosticIndex, d1$parent_maf),  "| -0.02\n")
cat("  DI vs pi_within     :", sp(d1$DiagnosticIndex, d1$pi_within),   "| -0.46\n")
cat("  DI vs dxy           :", sp(d1$DiagnosticIndex, d1$dxy),         "| +0.13\n")
d1[, rbin := cut(recomb, quantile(recomb, 0:10 / 10, na.rm = TRUE),
                 include.lowest = TRUE, labels = FALSE)]
arch_tab <- d1[, .(med_recomb = round(median(recomb), 2), n = .N,
                   DI = round(mean(DiagnosticIndex), 1),
                   Fst = round(mean(Fst, na.rm = TRUE), 3),
                   dxy = round(mean(dxy, na.rm = TRUE), 3),
                   pi_within = round(mean(pi_within, na.rm = TRUE), 3),
                   cluster_size = round(mean(cluster_size), 1)), by = rbin][order(rbin)]
cat("\n[PDF S1 / Fig 1a] architecture by recombination decile (rbin 1 = lowest recomb):\n")
print(arch_tab)

## ========================================================================
## == [PDF S2 + S3] classify independent units (LD-pruning reps) and SNPs
##    (the expensive step; cached and reused if the gate/threshold match)
## ========================================================================
run_ps <- function(markers) {
  prep <- ohta_fast_prepare(GTs[, markers], pops = smeta$Population)
  r <- parallelism_stats(prep, hybrid_pops = hyb, aqu_pops = aqu, pol_pops = pol,
                         DI = DI_by, min_DI = NULL, parent_maf = maf_by,
                         min_parent_maf = MIN_PARENT_MAF, sort_th = SORT_TH, fix_th = FIX_TH,
                         sort_rule = SORT_RULE)
  r[, `:=`(recomb = recomb_by[marker], n_loci = cs_by[marker])]; r[]
}
PS_PARAMS <- list(min_parent_maf = MIN_PARENT_MAF, sort_th = SORT_TH, fix_th = FIX_TH, sort_rule = SORT_RULE)
classify_cached <- function(path, mk_fun, label)
  cache_rds(path, label = label, valid = function(o) identical(o$params, PS_PARAMS),
            build = function() list(params = PS_PARAMS, r = run_ps(mk_fun())))$r
reps <- intersect(rep_snp_all$rep_snp, colnames(GTs))            # best-SNP LD-reduced unit set
cat("\n[B2] LD-reduced units:", length(reps), "| SNP sample:", SNP_SAMPLE, "\n")
t0 <- Sys.time()
r_unit <- classify_cached("moduleA_sorting/data/moduleA_r_unit.rds", function() reps,
                          "parallelism_stats (LD-reduced units)")
r_snp  <- classify_cached("moduleA_sorting/data/moduleA_r_snp.rds",
                          function() sample(colnames(GTs), SNP_SAMPLE),
                          "parallelism_stats (SNP sample)")
cat(sprintf("      classification done | %.0fs\n", elapsed(t0)))

## ========================================================================
## == [PDF S2] "Sorting is depleted, not enriched, in low-recombination regions."
##    feeds: unit-level 0.04 vs ~0.32, SNP-level 0.14; magnitude slope; Fig 1b
## ========================================================================
brk <- quantile(map$recomb, 0:10 / 10, na.rm = TRUE); uni <- c("aquilonia", "polyctena")
bin_it <- function(x) x[, rbin := cut(recomb, brk, include.lowest = TRUE, labels = FALSE)][]
r_unit <- bin_it(r_unit); r_snp <- bin_it(r_snp)
summ <- function(x, lab) x[differentiated == TRUE & !is.na(sort_class), .(
  level = lab, med_r = median(recomb, na.rm = TRUE), n = .N,
  frac_sorted = mean(sort_class != "unsorted"),
  frac_aqu_of_unidir = sum(sort_class == "aquilonia") / max(1L, sum(sort_class %in% uni))
), by = rbin][order(rbin)]
su <- summ(r_unit, "LD-reduced unit"); ss <- summ(r_snp, "SNP")
cat("\n[PDF S2 / Fig 1b] fraction sorted by recombination decile  [regenerated | PDF]:\n")
cat("  UNIT decile 1:", round(su$frac_sorted[1], 3), "| 0.06   ;  UNIT deciles 2-10 mean:",
    round(mean(su$frac_sorted[-1]), 3), "| ~0.42\n")
cat("  SNP  decile 1:", round(ss$frac_sorted[1], 3), "| 0.17\n")

## magnitude: prop_fixed increases with recombination (standardised)
du <- r_unit[differentiated == TRUE & is.finite(recomb) & is.finite(DI)]
du[, `:=`(zr  = as.numeric(scale(log10(recomb + 0.1))),
          zDI = as.numeric(scale(DI)),
          zmaf = as.numeric(scale(parent_maf)),
          zcs = as.numeric(scale(log10(n_loci))))]
lm_mag <- lm(prop_fixed ~ zr + zDI, data = du)
cat("\n[PDF S2] magnitude: prop_fixed ~ recombination (standardised), zr coefficient  [PDF: +0.05, z=117]:\n")
print(round(coef(summary(lm_mag))["zr", ], 4))

## ========================================================================
## == [PDF S3] "Direction is governed by diagnostic index; recombination adds
##             only a weak independent pull."
##    feeds: DI +1.57/z=112, recomb -0.09/z=-6.5, ~17x weaker; Fig 1c
## ========================================================================
cu <- du[sort_class %in% uni]; cu[, is_aqu := as.integer(sort_class == "aquilonia")]
glm_dir <- glm(is_aqu ~ zDI + zr + zmaf + zcs, data = cu, family = binomial())
co <- as.data.table(summary(glm_dir)$coefficients, keep.rownames = "term")
setnames(co, c("term", "estimate", "se", "z", "p"))
cat("\n[PDF S3] direction model P(aquilonia) ~ zDI + zr + zmaf + zcs (unidirectional units, n =", nrow(cu), "):\n")
print(co[, .(term, estimate = round(estimate, 3), z = round(z, 1), p = signif(p, 3))])
di_b <- co[term == "zDI", estimate]; r_b <- co[term == "zr", estimate]
cat(sprintf("  DI %+.2f (z=%.0f) [PDF +1.46/125] ; recomb %+.2f (z=%.1f) [PDF -0.09/-7.6] ; DI ~%.0fx recomb [PDF 16x]\n",
            di_b, co[term == "zDI", z], r_b, co[term == "zr", z], abs(di_b / r_b)))

## ========================================================================
## == [SUPP S3] direction coefficients vs sorting threshold
##    feeds: moduleA_supp_th_sensitivity.tex Fig S3 (moduleA_direction_sweep.png).
##    Reclassify unidirectional set at each sort_th (uni_score/bi_score stored),
##    refit with the z-covariates held fixed over the full differentiated set.
## ========================================================================
## Unidirectional set defined by the binom rule (classify_sort), reclassified at
## each sort_th -- consistent with the main direction model above (sort_class %in%
## uni), NOT the superseded |uni_score| >= bi_score component criterion.
dir_sweep <- rbindlist(lapply(SORT_TH_SWEEP, function(th) {
  du2 <- copy(du)
  du2[, sc := classify_sort(n_aqu, n_pol, n_obs, sort_th = th, sort_rule = "binom", alpha = ALPHA)]
  dd  <- du2[sc %in% uni]
  dd[, is_aqu := as.integer(sc == "aquilonia")]
  f  <- glm(is_aqu ~ zDI + zr + zmaf + zcs, data = dd, family = binomial())
  cc <- as.data.table(summary(f)$coefficients, keep.rownames = "term")
  setnames(cc, c("term", "estimate", "se", "z", "p"))
  cc[term != "(Intercept)"][, `:=`(sort_th = th, n_unidir = nrow(dd))][]
}))

## ========================================================================
## == [FIG 1] assembled 3-panel figure (a: architecture; b: sorting vs recomb;
##            c: direction vs recomb) -> moduleA_sorting/Figures/moduleA_architecture_fig.*
## ========================================================================
dir.create("Figures_PK", showWarnings = FALSE)
ok_col  <- c(Fst = "#D55E00", dxy = "#009E73", pi_within = "#0072B2")
lvl_col <- c("LD-reduced unit" = "#315B7D", "SNP" = "#C7CDD2")
tt <- theme_classic(base_size = 9) +
  theme(plot.tag = element_text(face = "bold", size = 11), axis.title = element_text(size = 8),
        legend.position = "bottom", legend.title = element_blank(),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 7))

## (a) Cruickshank-Hahn: standardised Fst / dxy / pi across recombination deciles
zsc <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
cha <- melt(arch_tab[, .(med_recomb, Fst = zsc(Fst), dxy = zsc(dxy), pi_within = zsc(pi_within))],
            id.vars = "med_recomb", variable.name = "metric", value.name = "z")
p_a <- ggplot(cha, aes(med_recomb, z, colour = metric)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) + scale_x_log10() +
  scale_colour_manual(values = ok_col) +
  labs(x = "recombination (cM/Mb, log)", y = "standardised (z)") + tt

## (b) fraction sorted vs recombination: unit vs SNP (the pseudo-replication test)
cmp <- rbind(su, ss); cmp[, level := factor(level, levels = c("SNP", "LD-reduced unit"))]
p_b <- ggplot(cmp, aes(med_r, frac_sorted, colour = level)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) + scale_x_log10() +
  scale_colour_manual(values = lvl_col) +
  labs(x = "recombination (cM/Mb, log)", y = "fraction sorted") + tt

## (c) direction vs recombination at unit level (DI polarity is not a recomb gradient)
p_c <- ggplot(su, aes(med_r, frac_aqu_of_unidir)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 0.6, colour = "grey45") + geom_point(size = 1.8, colour = "grey20") +
  scale_x_log10() + scale_y_continuous(limits = c(0, 1)) +
  labs(x = "recombination (cM/Mb, log)", y = "fraction fixing toward aquilonia") + tt

fig1 <- p_a + p_b + p_c + plot_annotation(tag_levels = "a")
ggsave("moduleA_sorting/Figures/moduleA_architecture_fig.pdf", fig1, width = 180, height = 72, units = "mm")
ggsave("moduleA_sorting/Figures/moduleA_architecture_fig.png", fig1, width = 180, height = 72, units = "mm", dpi = 300)

## [SUPP S3] figure: standardised coefficients vs sorting threshold
pal <- wes_palette("Zissou1", 1, type = "continuous")
dir_sweep[, term := factor(term, levels = c("zDI", "zr", "zmaf", "zcs"),
  labels = c("diagnostic index", "recombination", "parental MAF", "log block size"))]
p_sw <- ggplot(dir_sweep, aes(sort_th, estimate)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_line(colour = "#315B7D") + geom_point(colour = "#315B7D", size = 1.8) +
  geom_errorbar(aes(ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se),
                width = 0.02, colour = "#315B7D") +
  facet_wrap(~ term, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = SORT_TH_SWEEP) +
  labs(x = "sorting threshold", y = "standardised coefficient (95% CI)",
       title = "Module A -- direction-model coefficients vs sorting threshold",
       subtitle = "P(aquilonia) ~ DI + recombination + parental MAF + log block size (unidirectional units)") +
  theme_bw() + theme(strip.background = element_blank())
ggsave("moduleA_sorting/Figures/moduleA_direction_sweep.png", p_sw, width = 12, height = 3.2, dpi = 300)

## ========================================================================
## == [SUPP] within-cluster DI heterogeneity vs LD-cluster size
##    Justifies di_agg = "max": large low-recombination clusters mix diagnostic
##    and non-diagnostic markers, so a cluster's DI is not a single value. Uses
##    per-marker DI and the cluster membership from the same clustering.
## ========================================================================
di_marker <- setNames(map$DiagnosticIndex, map$marker)
clmem <- data.table(group_id = rep(g$group_id, lengths(g$members)),
                    n_loci   = rep(g$n_loci,   lengths(g$members)),
                    DI       = di_marker[unlist(g$members)])
dic <- clmem[, .(n_loci = n_loci[1], di_sd = sd(DI, na.rm = TRUE),
                 di_range = { d <- DI[!is.na(DI)]; if (length(d) >= 2) max(d) - min(d) else NA_real_ },
                 di_mean = mean(DI, na.rm = TRUE), n_di = sum(!is.na(DI))), by = group_id
               ][n_loci >= 2 & n_di >= 2]
dic[, size_class := cut(n_loci, breaks = c(1, 2, 4, 9, 49, 199, Inf),
                        labels = c("2", "3-4", "5-9", "10-49", "50-199", "200+"))]
di_size_tab <- dic[, .(n_clusters = .N, median_di_sd = round(median(di_sd), 2),
                       median_di_range = round(median(di_range), 1),
                       median_di_mean = round(median(di_mean), 1)), by = size_class][order(size_class)]
gw_di_sd <- sd(map$DiagnosticIndex, na.rm = TRUE)
di_sd_size_spearman <- cor(dic$n_loci, dic$di_sd, method = "spearman", use = "complete.obs")
cat("\n[SUPP] within-cluster DI sd by cluster size (di_agg = 'max' justification):\n"); print(di_size_tab)
cat(sprintf("  Spearman cor(n_loci, within-cluster DI sd) = %.3f ; genome-wide marker DI sd = %.1f\n",
            di_sd_size_spearman, gw_di_sd))

p_div <- ggplot(dic, aes(size_class, di_sd)) +
  geom_hline(yintercept = gw_di_sd, linetype = 2, colour = "#D08A45") +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.15, fill = "#8EB4CE", linewidth = 0.3) +
  labs(x = "LD-cluster size (n markers)", y = "within-cluster DI standard deviation",
       title = "Within-cluster diagnostic-index spread vs LD-cluster size",
       subtitle = sprintf("%d multi-marker clusters (%s); dashed = genome-wide marker DI sd (%.1f)",
                          nrow(dic), basename(CLUSTERING), gw_di_sd)) +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
ggsave("moduleA_sorting/Figures/moduleA_di_variance.pdf", p_div, width = 165, height = 100, units = "mm")
ggsave("moduleA_sorting/Figures/moduleA_di_variance.png", p_div, width = 165, height = 100, units = "mm", dpi = 300)

## ---- save the result objects -------------------------------------------
saveRDS(list(correlations = list(DI_recomb = sp(d1$DiagnosticIndex, d1$recomb),
                                 DI_maf = sp(d1$DiagnosticIndex, d1$parent_maf),
                                 DI_pi = sp(d1$DiagnosticIndex, d1$pi_within),
                                 DI_dxy = sp(d1$DiagnosticIndex, d1$dxy)),
             arch_tab = arch_tab, unit_by_recomb = su, snp_by_recomb = ss,
             lm_magnitude = lm_mag, glm_direction = glm_dir,
             direction_coefs = co, direction_sweep = dir_sweep,
             di_size_tab = di_size_tab, di_sd_size_spearman = di_sd_size_spearman,
             params = c(PS_PARAMS, list(clustering = CLUSTERING, snp_sample = SNP_SAMPLE))),
        "moduleA_sorting/data/moduleA_architecture.rds")

cat("\nModule A architecture complete. Outputs:\n",
    "  moduleA_sorting/data/moduleA_architecture.rds\n",
    "  moduleA_sorting/Figures/moduleA_architecture_fig.{pdf,png}   (Fig 1 a/b/c)\n",
    "  moduleA_sorting/Figures/moduleA_direction_sweep.png          (Supp S3)\n",
    "  moduleA_sorting/Figures/moduleA_di_variance.{pdf,png}        (Supp: DI spread vs cluster size)\n")
