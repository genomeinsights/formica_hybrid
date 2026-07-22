## =========================================================
## Tier 2: does sorting concentrate in low-recombination regions,
## at the INDEPENDENT-UNIT (eMLG) level, with DI as a separable covariate?
## =========================================================
## Key design points (from the theory + Tier 1):
##  - Neutral per-UNIT sorting should be ~recombination-independent; the SNP-
##    level recomb-sorting slope is largely the spatial/redundancy artifact
##    (adjacent SNPs redundant in low-recomb blocks). So we compare SNP-level
##    vs unit-level: a slope that FLATTENS at unit level was the artifact; a
##    slope that PERSISTS is a candidate linked-selection / DMI signal.
##  - Use map-based recombination (exogenous), NOT cluster size (collider).
##  - DI is ~orthogonal to recombination (Tier 1), so include it as a
##    separable covariate.
##  - Units = LD-pruning representatives (one tag SNP per cluster), taken from
##    GTs so parents are present for orientation.
##  - NB: no null yet, so this is descriptive; "excess over neutral" awaits
##    the recombination-matched SLiM null.

suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
source("dev/R/Ohta.R"); source("dev/R/parallelism_stats.R")
set.seed(1)

e <- new.env(); load("./data/hybrids_and_parents_maf005.Rdata", envir=e)
GTs <- e$GTs_with_parents; sd <- e$sample_data_with_parents; map <- copy(e$map_hyb_005)

## ---- per-marker recombination ----
rec <- fread("./data/Frufa_DTOL_PR.ref_genome.recmap")
setnames(rec, c("chr","pos","cM","cMMb")); rec[, Chr := sub("chromosome_","Chr",chr)]
map[, recomb := NA_real_]
for (ch in unique(map$Chr)) {
  r <- rec[Chr==ch]; if (nrow(r) < 2) next
  idx <- map[, which(Chr==ch)]
  map[idx, recomb := approx(r$pos, r$cMMb, xout=map$Pos[idx], rule=2)$y]
}
recomb_by <- setNames(map$recomb, map$marker)
DI_by     <- setNames(map$DiagnosticIndex, map$marker)

## ---- independent units = LD-pruning representatives ----
g <- readRDS("data/eMLG_5loci_0025_cM05.rds")$groups
reps <- intersect(g$representative, colnames(GTs))
cs_by <- setNames(g$n_loci, g$representative)      # cluster size per representative
cat("independent units (representatives):", length(reps), "\n")

aqu <- "aquilonia_parent"; pol <- "polyctena_parent"
hyb <- setdiff(unique(sd$Population), c(aqu, pol))

run_ps <- function(markers) {
  prep <- ohta_fast_prepare(GTs[, markers], pops = sd$Population)
  r <- parallelism_stats(prep, hyb, aqu, pol, DI = DI_by,
                         min_parent_diff = 0.5, sort_th = 0.5, fix_th = 0.1)
  r[, recomb := recomb_by[marker]]
  r
}

## unit-level (all representatives) and SNP-level (random 200k SNPs)
r_unit <- run_ps(reps)
r_unit[, n_loci := cs_by[marker]]
r_snp  <- run_ps(sample(colnames(GTs), 200000))

## common recombination deciles (from genome-wide distribution)
brk <- quantile(map$recomb, 0:10/10, na.rm=TRUE)
bin_it <- function(d) { d[, rbin := cut(recomb, brk, include.lowest=TRUE, labels=FALSE)]
                        d[, med_r := median(recomb, na.rm=TRUE), by=rbin][] }
r_unit <- bin_it(r_unit); r_snp <- bin_it(r_snp)

sorted <- function(sc) sc %in% c("aquilonia","polyctena","bidirectional")
summ <- function(d, lab) d[differentiated==TRUE, .(
    level=lab, med_r=median(recomb,na.rm=TRUE), n=.N,
    frac_sorted = mean(sorted(sort_class)),
    frac_uni    = mean(sort_class %in% c("aquilonia","polyctena")),
    frac_bi     = mean(sort_class == "bidirectional"),
    prop_fixed  = mean(prop_fixed, na.rm=TRUE)), by=rbin][order(rbin)]

su <- summ(r_unit, "unit (eMLG)"); ss <- summ(r_snp, "SNP")
cat("\n=== fraction sorted by recomb decile: UNIT level ===\n"); print(su)
cat("\n=== fraction sorted by recomb decile: SNP level ===\n"); print(ss)

## ---- regression at unit level: sorting ~ recomb + DI (separable) ----
du <- r_unit[differentiated==TRUE & is.finite(recomb) & is.finite(DI)]
du[, `:=`(zr = scale(log10(recomb+0.1)), zDI = scale(DI), zcs = scale(log10(n_loci)))]
cat("\n=== unit-level regression: prop_fixed ~ recomb + DI ===\n")
print(round(coef(summary(lm(prop_fixed ~ zr + zDI, data=du))),4))
cat("\n=== add cluster size (collider - for contrast only) ===\n")
print(round(coef(summary(lm(prop_fixed ~ zr + zDI + zcs, data=du))),4))
cat("\n=== directional-only: |uni_score| ~ recomb + DI ===\n")
print(round(coef(summary(lm(abs(uni_score) ~ zr + zDI, data=du))),4))

## ---- plots ----
th <- theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(), legend.position="bottom")
cmp <- rbindlist(list(su, ss))
pA <- ggplot(cmp, aes(med_r, frac_sorted, color=level)) + geom_line(linewidth=1) + geom_point() +
  scale_x_log10() + scale_color_manual(values=c("SNP"="grey55","unit (eMLG)"="firebrick")) +
  labs(x="recombination (cM/Mb, log)", y="fraction sorted",
       title="(A) SNP-level vs independent-unit sorting vs recombination", color=NULL) + th
pB <- ggplot(melt(su[, .(med_r, directional=frac_uni, bidirectional=frac_bi)], id.vars="med_r"),
             aes(med_r, value, color=variable)) + geom_line(linewidth=1) + geom_point() +
  scale_x_log10() + labs(x="recombination (cM/Mb, log)", y="fraction of units",
       title="(B) unit-level: directional vs bidirectional sorting", color=NULL) + th
## DI partial effect: bin DI within recomb tertiles
du[, rtert := cut(recomb, quantile(recomb,0:3/3,na.rm=TRUE), include.lowest=TRUE,
                  labels=c("low r","mid r","high r"))]
du[, DIbin := cut(DI, quantile(DI,0:6/6,na.rm=TRUE), include.lowest=TRUE)]
pd <- du[!is.na(DIbin), .(med_DI=median(DI), prop_fixed=mean(prop_fixed,na.rm=TRUE)), by=.(rtert,DIbin)]
pC <- ggplot(pd, aes(med_DI, prop_fixed, color=rtert)) + geom_line(linewidth=1) + geom_point() +
  labs(x="DI", y="prop_fixed", title="(C) DI effect within recombination tertiles", color=NULL) + th

p <- pA / (pB | pC)
dir.create("Figures", showWarnings=FALSE)
ggsave("Figures/tier2_sorting_vs_recomb.png", p, width=12, height=10, dpi=110)
cat("\nsaved Figures/tier2_sorting_vs_recomb.png\n")
