## =========================================================
## MODULE B -- genomic architecture of differentiation & sorting
## =========================================================
## Consolidates and updates the former tier1_architecture.R + tier2_sorting_vs_recomb.R
## onto the canonical clustering and the Module A conventions.
##
##   B1  architecture (per marker): DI vs recombination / within-species pi /
##       absolute divergence d_xy / relative F_ST / cluster size. Answers whether
##       DI is a diversity artifact (low pi) or genuine divergence (high d_xy),
##       and whether DI is separable from recombination.
##   B2  sorting vs recombination: fraction sorted by recombination decile at the
##       INDEPENDENT-UNIT level (LD-pruning representatives) vs SNP level. A slope
##       that FLATTENS unit->SNP was the spatial-redundancy artifact; one that
##       PERSISTS is a candidate linked-selection signal.
##   B3  direction x architecture (the centrepiece): does recombination (or size)
##       predict the DIRECTION of sorting BEYOND DI? Module A showed DI governs
##       direction (low DI -> polyctena, high DI -> aquilonia); if low recombination
##       adds an independent aquilonia pull after DI, that is a linked-selection /
##       reduced-local-Ne lead (intrinsic side). If DI absorbs it, direction is a
##       pure diagnostic-index property.
##
## Conventions locked to match Module A: primary gate parent_maf >= 0.15 (pooled
## parental MAF), fix_th 0.15, sort_th 0.5, null_prob 0.5, DI ungated (covariate).
##
## DESCRIPTIVE only: a non-recombining autosomal block does NOT have lower Ne than
## its surroundings under neutrality, so any recombination->sorting/direction effect
## needs the recombination-matched null (Module E) before a causal (linked-selection)
## reading. Units = LD-pruning representatives genome-wide (consensus for all ~474k
## clusters is prohibitive; Module A's consensus results are the sorted-cluster view).
##
## Run from repo root. Reads only GTs/map/recmap/clustering; writes Figures/moduleB_fig2.*
## and data/moduleB_architecture.rds. Sources Ohta.R + parallelism_stats.R.

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("dev/R/Ohta.R"); source("dev/R/parallelism_stats.R")
set.seed(1)

## ---- PARAMETERS (match Module A) ---------------------------------------
MIN_PARENT_MAF <- 0.15
FIX_TH         <- 0.15
SORT_TH        <- 0.5
NULL_PROB      <- 0.5
SNP_SAMPLE     <- 200000L                       # SNP-level comparison sample
CLUSTERING     <- "data/eMLG_5loci_0025_cM05.rds"
RECMAP         <- "data/Frufa_DTOL_PR.ref_genome.recmap"

elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

## ---- inputs -------------------------------------------------------------
e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
GTs <- e$GTs_with_parents; sd <- e$sample_data_with_parents; map <- copy(e$map_hyb_005)

## per-marker recombination rate (cM/Mb), map-interpolated (exogenous, not size)
rec <- fread(RECMAP)
stopifnot(ncol(rec) >= 4)                        # chr, pos, cM, cM/Mb
setnames(rec, 1:4, c("chr", "pos", "cM", "cMMb")); rec[, Chr := sub("chromosome_", "Chr", chr)]
map[, recomb := NA_real_]
for (ch in unique(map$Chr)) {
  r <- rec[Chr == ch]; if (nrow(r) < 2) next
  idx <- map[, which(Chr == ch)]
  map[idx, recomb := approx(r$pos, r$cMMb, xout = map$Pos[idx], rule = 2)$y]
}

## per-marker cluster size (from the canonical clustering)
g <- readRDS(CLUSTERING)$groups
memb <- data.table(marker = unlist(g$members),
                   cluster_size = rep(g$n_loci, lengths(g$members)))
map <- memb[map, on = "marker"]

## parental frequencies -> pi (within), d_xy (absolute), F_ST (relative), pooled MAF
aqu <- "aquilonia_parent"; pol <- "polyctena_parent"
aqu_rows <- which(sd$Population == aqu); pol_rows <- which(sd$Population == pol)
pa <- (colMeans(GTs[aqu_rows, , drop = FALSE], na.rm = TRUE) / 2)[map$marker]
pp <- (colMeans(GTs[pol_rows, , drop = FALSE], na.rm = TRUE) / 2)[map$marker]
par_pool <- (colMeans(GTs[c(aqu_rows, pol_rows), , drop = FALSE], na.rm = TRUE) / 2)[map$marker]
Hs_a <- 2 * pa * (1 - pa); Hs_p <- 2 * pp * (1 - pp)
map[, pi_within := (Hs_a + Hs_p) / 2]
map[, dxy       := pa * (1 - pp) + pp * (1 - pa)]
pbar <- (pa + pp) / 2; Ht <- 2 * pbar * (1 - pbar)
map[, Fst        := ifelse(Ht > 0, (Ht - (Hs_a + Hs_p) / 2) / Ht, NA_real_)]
map[, parent_maf := pmin(par_pool, 1 - par_pool)]

DI_by       <- setNames(map$DiagnosticIndex, map$marker)
maf_by      <- setNames(map$parent_maf,      map$marker)
recomb_by   <- setNames(map$recomb,          map$marker)

## ========================================================================
## B1 -- architecture (per marker)
## ========================================================================
d1 <- map[is.finite(recomb) & is.finite(DiagnosticIndex) & !is.na(cluster_size)]
cat("[B1] markers used:", nrow(d1), "\n")
sp <- function(a, b) round(cor(a, b, method = "spearman", use = "complete.obs"), 3)
cat("\n=== Spearman correlations ===\n")
cat("  DI vs recomb        :", sp(d1$DiagnosticIndex, d1$recomb), "\n")
cat("  DI vs cluster_size  :", sp(d1$DiagnosticIndex, d1$cluster_size), "\n")
cat("  DI vs pi_within     :", sp(d1$DiagnosticIndex, d1$pi_within), "\n")
cat("  DI vs dxy           :", sp(d1$DiagnosticIndex, d1$dxy), "\n")
cat("  DI vs parent_maf    :", sp(d1$DiagnosticIndex, d1$parent_maf), "\n")
cat("  recomb vs cluster   :", sp(d1$recomb, d1$cluster_size), "\n")
cat("  recomb vs pi_within :", sp(d1$recomb, d1$pi_within), "\n")
cat("  recomb vs dxy       :", sp(d1$recomb, d1$dxy), "\n")
cat("  recomb vs Fst       :", sp(d1$recomb, d1$Fst), "\n")

d1[, rbin := cut(recomb, quantile(recomb, 0:10/10, na.rm = TRUE),
                 include.lowest = TRUE, labels = FALSE)]
arch_tab <- d1[, .(med_recomb = round(median(recomb), 2), n = .N,
                   DI = round(mean(DiagnosticIndex), 1),
                   Fst = round(mean(Fst, na.rm = TRUE), 3),
                   dxy = round(mean(dxy, na.rm = TRUE), 3),
                   pi_within = round(mean(pi_within, na.rm = TRUE), 3),
                   cluster_size = round(mean(cluster_size), 1)), by = rbin][order(rbin)]
cat("\n=== means by recombination decile (rbin 1 = lowest recomb) ===\n"); print(arch_tab)

## ========================================================================
## B2 -- sorting vs recombination: unit (LD-pruning rep) vs SNP level
## ========================================================================
reps  <- intersect(g$representative, colnames(GTs))
cs_by <- setNames(g$n_loci, g$representative)
hyb   <- setdiff(unique(sd$Population), c(aqu, pol))
cat("\n[B2] independent units (representatives):", length(reps), "\n")

run_ps <- function(markers) {
  prep <- ohta_fast_prepare(GTs[, markers], pops = sd$Population)
  r <- parallelism_stats(prep, hybrid_pops = hyb, aqu_pops = aqu, pol_pops = pol,
                         DI = DI_by, min_DI = NULL,
                         parent_maf = maf_by, min_parent_maf = MIN_PARENT_MAF,
                         sort_th = SORT_TH, fix_th = FIX_TH, null_prob = NULL_PROB)
  r[, `:=`(recomb = recomb_by[marker], n_loci = cs_by[marker])]
  r
}
t0 <- Sys.time()
r_unit <- run_ps(reps)
r_snp  <- run_ps(sample(colnames(GTs), SNP_SAMPLE))
cat(sprintf("      parallelism_stats (unit + SNP) done | %.0fs\n", elapsed(t0)))

brk <- quantile(map$recomb, 0:10/10, na.rm = TRUE)
uni <- c("aquilonia", "polyctena")
bin_it <- function(x) { x[, rbin := cut(recomb, brk, include.lowest = TRUE, labels = FALSE)][] }
r_unit <- bin_it(r_unit); r_snp <- bin_it(r_snp)
summ <- function(x, lab) x[differentiated == TRUE & !is.na(sort_class), .(
  level = lab, med_r = median(recomb, na.rm = TRUE), n = .N,
  frac_sorted = mean(sort_class != "unsorted"),
  frac_uni    = mean(sort_class %in% uni),
  frac_bi     = mean(sort_class == "bidirectional"),
  frac_aqu_of_unidir = sum(sort_class == "aquilonia") / max(1L, sum(sort_class %in% uni)),
  prop_fixed  = mean(prop_fixed, na.rm = TRUE)), by = rbin][order(rbin)]
su <- summ(r_unit, "unit (eMLG rep)"); ss <- summ(r_snp, "SNP")
cat("\n=== fraction sorted by recomb decile: UNIT level ===\n"); print(su)
cat("\n=== fraction sorted by recomb decile: SNP level ===\n");  print(ss)

## ========================================================================
## B3 -- direction x architecture (does recomb predict direction beyond DI?)
## ========================================================================
du <- r_unit[differentiated == TRUE & is.finite(recomb) & is.finite(DI)]
du[, `:=`(zr = as.numeric(scale(log10(recomb + 0.1))),
          zDI = as.numeric(scale(DI)),
          zmaf = as.numeric(scale(parent_maf)),
          zcs = as.numeric(scale(log10(n_loci))))]

cat("\n=== B3.1 magnitude: prop_fixed ~ recomb + DI (unit level) ===\n")
print(round(coef(summary(lm(prop_fixed ~ zr + zDI, data = du))), 4))

cu <- du[sort_class %in% uni]; cu[, is_aqu := as.integer(sort_class == "aquilonia")]
cat("\n=== B3.2 DIRECTION: P(aquilonia) ~ DI + recomb + parent_maf + log(size) ===\n")
cat("    (unidirectional units; does recomb pull toward aquilonia AFTER DI?)\n")
fitB <- glm(is_aqu ~ zDI + zr + zmaf + zcs, data = cu, family = binomial())
print(round(summary(fitB)$coefficients, 4))
cat("\n(zr = standardised log recombination. A negative, significant zr => LOW\n",
    " recombination independently leans aquilonia after DI/MAF/size -- a linked-\n",
    " selection lead, DESCRIPTIVE pending the Module E null. A ~0 zr => DI absorbs it.)\n")

saveRDS(list(arch_tab = arch_tab, unit_by_recomb = su, snp_by_recomb = ss,
             lm_magnitude = lm(prop_fixed ~ zr + zDI, data = du), glm_direction = fitB),
        "data/moduleB_architecture.rds")

## ========================================================================
## Figure 2 (a-c)
## ========================================================================
dir.create("Figures", showWarnings = FALSE)
ok_col <- c(Fst = "#D55E00", dxy = "#009E73", pi_within = "#0072B2")
lvl_col <- c("unit (eMLG rep)" = "#0072B2", "SNP" = "#999999")
th <- theme_classic(base_size = 9) +
  theme(plot.tag = element_text(face = "bold", size = 11),
        axis.title = element_text(size = 8), legend.position = "bottom",
        legend.title = element_blank(), legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7))

## (a) Cruickshank-Hahn: standardised Fst / dxy / pi across recombination deciles
##     (DI dropped -- it is flat vs recombination; z-scaling only amplified noise)
zsc <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
cha <- melt(arch_tab[, .(med_recomb, Fst = zsc(Fst),
                         dxy = zsc(dxy), pi_within = zsc(pi_within))],
            id.vars = "med_recomb", variable.name = "metric", value.name = "z")
p2a <- ggplot(cha, aes(med_recomb, z, colour = metric)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) + scale_x_log10() +
  scale_colour_manual(values = ok_col) +
  labs(x = "recombination (cM/Mb, log)", y = "standardised (z)") + th

## (b) fraction sorted vs recombination: unit vs SNP (the flattening test)
cmp <- rbindlist(list(su, ss))
p2b <- ggplot(cmp, aes(med_r, frac_sorted, colour = level)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) + scale_x_log10() +
  scale_colour_manual(values = lvl_col) +
  labs(x = "recombination (cM/Mb, log)", y = "fraction sorted") + th

## (c) direction vs recombination at unit level (is the DI flip also a recomb flip?)
p2c <- ggplot(su, aes(med_r, frac_aqu_of_unidir)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 0.6, colour = "grey45") + geom_point(size = 1.8, colour = "grey20") +
  scale_x_log10() + scale_y_continuous(limits = c(0, 1)) +
  labs(x = "recombination (cM/Mb, log)", y = "fraction fixing toward aquilonia") + th

fig2 <- p2a + p2b + p2c + plot_annotation(tag_levels = "a")
ggsave("Figures/moduleB_fig2.pdf", fig2, width = 180, height = 72, units = "mm")
ggsave("Figures/moduleB_fig2.png", fig2, width = 180, height = 72, units = "mm", dpi = 300)
cat("\nSaved: data/moduleB_architecture.rds, Figures/moduleB_fig2.pdf/.png\n")
