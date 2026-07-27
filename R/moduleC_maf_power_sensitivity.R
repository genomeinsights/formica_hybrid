## =========================================================
## MODULE C sensitivity: is the DI-enrichment a BayPass-power / allele-frequency
## artifact?
## =========================================================
## CONCERN. BayPass evidence scales with among-population variation in hybrid
## allele frequency: a cluster whose alleles are near-fixed genome-wide in the
## hybrids has little among-population variance and little power (low p_sig); one
## segregating at intermediate frequency has more. If those higher-power clusters
## are ALSO more ancestry-diagnostic, the primary enrichment (frac_hi ~ p_sig)
## could reflect a shared dependence on allele-frequency variation rather than a
## climate <-> diagnostic link.
##
## TEST. Add simple per-cluster power / variation covariates to the PRIMARY
## cluster-level model and watch the p_sig coefficient. If it is stable, the
## power/MAF explanation is closed; if it collapses toward 0, the enrichment was
## largely such an artifact.
##   M0 : frac_hi ~ p_sig + log(n_loci)                       (published primary)
##   M1 :        + maf_hyb        (hybrid MAF, cluster mean)
##   M2 :        + maf_hyb + recomb
##   M3 :        + maf_hyb + recomb + xtx   (BayPass XtX = the direct
##                                           among-population-differentiation /
##                                           power measure, cluster mean)
## The p_sig coefficient is reported on the published scale (per +10 pp of
## significance rate) at each step, with its naive 95% CI, and -- because the
## naive CIs overstate precision here (see doc 05) -- the fully-adjusted (M3)
## coefficient additionally carries a chromosome block-bootstrap 95% CI.
##
## Primary config only (aland_excluded x withOmega), n_loci >= 50, PC1 and PC2 --
## exactly the stratum/config in which the enrichment was reported.
##
## Reads : data/hybrids_only_maf005.Rdata (map only, incl. maf_hyb),
##         data/Frufa_DTOL_PR.ref_genome.recmap, data/eMLG_5loci_0025_cM05.rds,
##         ./aland_excluded/PC{1,2}_DIEM_withOmega_summary_{betai_reg,pi_xtx}.out
## Writes: data/moduleC_maf_power_sensitivity.rds,
##         Figures/moduleC_maf_power_sensitivity.{pdf,png}
## Run from the repo root. DESCRIPTIVE, like the rest of Module C.

suppressMessages({ library(data.table); library(ggplot2) })
set.seed(1)

## ---- parameters (match Module C primary) --------------------------------
SIG_THR    <- 15; DI_TH <- -25; MIN_NLOCI <- 50L
POPSET     <- "aland_excluded"; OMEGA <- "withOmega"
CLUSTERING <- "data/eMLG_5loci_0025_cM05.rds"
RECMAP     <- "data/Frufa_DTOL_PR.ref_genome.recmap"
B          <- 10000L               # block-bootstrap replicates (chromosome blocks)
CI         <- c(0.025, 0.975)
PUBLISHED  <- c(PC1 = 4.457, PC2 = 3.717)   # M0 p_sig*10, for the self-validation gate

pct_ci <- function(x) { x <- x[is.finite(x)]; if (!length(x)) c(NA, NA) else unname(quantile(x, CI)) }

## ---- inputs -------------------------------------------------------------
em <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = em)
map <- as.data.table(em$map_hyb_005)[, .(marker, Chr, Pos, DiagnosticIndex, maf_hyb)]
rm(em); invisible(gc())

rec <- fread(RECMAP); setnames(rec, 1:4, c("chr", "pos", "cM", "cMMb"))
rec[, Chr := sub("chromosome_", "Chr", chr)]
map[, recomb := NA_real_]
for (ch in unique(map$Chr)) {
  r <- rec[Chr == ch]; if (nrow(r) < 2) next
  idx <- map[, which(Chr == ch)]
  map[idx, recomb := approx(r$pos, r$cMMb, xout = map$Pos[idx], rule = 2)$y]
}
maf_by    <- setNames(map$maf_hyb,        map$marker)
recomb_by <- setNames(map$recomb,         map$marker)
DI_by     <- setNames(map$DiagnosticIndex, map$marker)

groups <- readRDS(CLUSTERING)$groups
he  <- groups[has_eMLG == TRUE]
m2g <- he[, .(marker = unlist(members)), by = group_id]
m2g[, `:=`(DI = DI_by[marker], maf_hyb = maf_by[marker], recomb = recomb_by[marker])]

## per-cluster, config-invariant summaries: DI counts + covariate means
base_cl <- m2g[, .(n_hi = sum(DI > DI_TH, na.rm = TRUE), n_di = sum(!is.na(DI)),
                   maf_hyb = mean(maf_hyb, na.rm = TRUE),
                   recomb  = mean(recomb,  na.rm = TRUE)), by = group_id]
## chromosome per cluster (block key) via its representative marker
reppos <- map[data.table(marker = he$representative, group_id = he$group_id),
              on = "marker"][, .(group_id, Chr)]
base_cl <- reppos[base_cl, on = "group_id"]

## ---- per-marker BayPass BF and XtX (positional, guarded) ----------------
import_col <- function(kind, pc, col) {
  f <- sprintf("./%s/%s_DIEM_%s_summary_%s.out", POPSET, pc, OMEGA, kind)
  r <- fread(f)
  stopifnot(nrow(r) == nrow(map), identical(r$MRK, seq_len(nrow(r))), col %in% names(r))
  setNames(r[[col]], map$marker)
}

## ---- assemble the per-cluster modelling frame for one PC ----------------
build_frame <- function(pc) {
  bf  <- import_col("betai_reg", pc, "BF(dB)")
  xtx <- import_col("pi_xtx",    pc, "M_XtX")
  m2g[, `:=`(BF = bf[marker], XT = xtx[marker])]
  a <- m2g[, .(n_loci = .N,
               p_sig = sum(BF >= SIG_THR, na.rm = TRUE) / .N,
               xtx   = mean(XT, na.rm = TRUE)), by = group_id]
  D <- base_cl[a, on = "group_id"]
  D <- D[n_di > 0 & n_loci >= MIN_NLOCI]
  D[, `:=`(frac_hi = n_hi / n_di, logn = log(n_loci))]
  D[stats::complete.cases(D[, .(frac_hi, p_sig, logn, maf_hyb, recomb, xtx)])]
}

## ---- the model ladder; report p_sig coefficient on the published scale ---
LADDER <- list(
  M0 = frac_hi ~ p_sig + logn,
  M1 = frac_hi ~ p_sig + logn + maf_hyb,
  M2 = frac_hi ~ p_sig + logn + maf_hyb + recomb,
  M3 = frac_hi ~ p_sig + logn + maf_hyb + recomb + xtx)

fit_ladder <- function(D, pc) {
  rbindlist(lapply(names(LADDER), function(m) {
    s <- summary(lm(LADDER[[m]], data = D))$coefficients["p_sig", ]
    z <- qnorm(0.975)
    data.table(pc = pc, model = m,
               added = c(M0 = "(baseline)", M1 = "+ maf_hyb", M2 = "+ recomb",
                         M3 = "+ xtx")[m],
               effect = unname(s[1]) * 10,                 # per +10 pp significance rate
               lo = unname(s[1] - z * s[2]) * 10,
               hi = unname(s[1] + z * s[2]) * 10,
               p = unname(s[4]), n = nrow(D))
  }))
}

## ---- chromosome block bootstrap of the M3 p_sig coefficient -------------
boot_M3 <- function(D, pc) {
  X <- model.matrix(LADDER$M3, D); y <- D$frac_hi
  by_block <- split(seq_len(nrow(D)), D$Chr); blocks <- names(by_block)
  full <- .lm.fit(X, y)$coefficients[2] * 10                # col 2 = p_sig
  bs <- numeric(B)
  for (b in seq_len(B)) {
    idx <- unlist(by_block[sample(blocks, length(blocks), replace = TRUE)], use.names = FALSE)
    bs[b] <- tryCatch(.lm.fit(X[idx, , drop = FALSE], y[idx])$coefficients[2] * 10,
                      error = function(e) NA_real_)
  }
  ci <- pct_ci(bs)
  data.table(pc = pc, model = "M3", boot_lo = ci[1], boot_hi = ci[2],
             boot_excl0 = !(ci[1] <= 0 & 0 <= ci[2]), n_fail = sum(!is.finite(bs)))
}

## ========================================================================
## run for PC1 and PC2
## ========================================================================
frames <- lapply(c("PC1", "PC2"), build_frame); names(frames) <- c("PC1", "PC2")
ladder <- rbindlist(Map(fit_ladder, frames, names(frames)))
cat("\n=== p_sig (enrichment) coefficient across the covariate ladder ===\n")
cat("(per +10 pp significance rate; primary config, n_loci >= 50)\n")
print(ladder[, .(pc, model, added, effect = round(effect, 3),
                 naive_lo = round(lo, 3), naive_hi = round(hi, 3), p = signif(p, 3), n)])

## self-validation: M0 must reproduce the published enrichment
m0 <- ladder[model == "M0", setNames(effect, pc)]
cat("\n=== self-validation: M0 vs published ===\n")
for (pc in names(PUBLISHED))
  cat(sprintf("  %s  M0 = %.3f  published = %.3f  %s\n", pc, m0[pc], PUBLISHED[pc],
              ifelse(abs(m0[pc] - PUBLISHED[pc]) <= 0.1, "OK", "MISMATCH")))
if (any(abs(m0[names(PUBLISHED)] - PUBLISHED) > 0.1)) stop("M0 does not reproduce the published enrichment")

message(sprintf("[boot] M3 chromosome block bootstrap, B = %d ...", B))
boot <- rbindlist(Map(boot_M3, frames, names(frames)))
res  <- merge(ladder[model == "M3"], boot, by = c("pc", "model"))
cat("\n=== fully-adjusted (M3) p_sig coefficient: naive vs block-bootstrap CI ===\n")
print(res[, .(pc, effect = round(effect, 3), naive_lo = round(lo, 3), naive_hi = round(hi, 3),
              boot_lo = round(boot_lo, 3), boot_hi = round(boot_hi, 3), boot_excl0, n_fail)])

saveRDS(list(ladder = ladder, m3_boot = boot, frames_n = sapply(frames, nrow),
             params = list(SIG_THR = SIG_THR, DI_TH = DI_TH, MIN_NLOCI = MIN_NLOCI,
                           config = c(POPSET, OMEGA), B = B)),
        "data/moduleC_maf_power_sensitivity.rds")

## ========================================================================
## figure: the p_sig coefficient ladder (naive CI), M3 also with bootstrap CI
## ========================================================================
ladder[, model := factor(model, levels = c("M0", "M1", "M2", "M3"))]
labdt <- unique(ladder[, .(model, added)])
lab <- setNames(sprintf("%s\n%s", labdt$model, labdt$added), as.character(labdt$model))
p <- ggplot(ladder, aes(model, effect)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey65") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, colour = "#0072B2") +
  geom_errorbar(data = res, aes(x = "M3", ymin = boot_lo, ymax = boot_hi),
                width = 0.32, colour = "#D55E00", linewidth = 0.6, inherit.aes = FALSE) +
  geom_point(size = 1.8) +
  facet_wrap(~ pc) +
  scale_x_discrete(labels = lab) +
  labs(x = NULL, y = "enrichment: pp high-DI per +10 pp climate rate",
       title = "Blue = naive 95% CI; orange (M3) = chromosome block-bootstrap 95% CI") +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 7.5), strip.text = element_text(face = "bold"))
dir.create("Figures", showWarnings = FALSE)
ggsave("Figures/moduleC_maf_power_sensitivity.pdf", p, width = 150, height = 80, units = "mm")
ggsave("Figures/moduleC_maf_power_sensitivity.png", p, width = 150, height = 80, units = "mm", dpi = 300)
cat("\nSaved: data/moduleC_maf_power_sensitivity.rds, Figures/moduleC_maf_power_sensitivity.pdf/.png\n")
