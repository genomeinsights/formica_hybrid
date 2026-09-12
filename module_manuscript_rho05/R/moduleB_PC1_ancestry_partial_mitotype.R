## =========================================================
## module_manuscript_rho05 -- does PC1's ancestry link survive controlling
## for mitotype?
## =========================================================
## FULLY SUPERSEDED (2026-09-12 audit): this script has the Omega-null
## variable-substitution bug (null draws replaced PC1, holding ancestry fixed,
## instead of replacing ancestry -- the quantity Omega actually models).
## Nothing downstream reads this script's own output
## (data/moduleB_PC1_ancestry_partial_mitotype.rds, moduleB_PC1_ancestry_
## partial_mitotype.png -- both quarantined to stale_pre_fix_20260912/). Do
## not run this script or cite its numbers; use
## moduleB_ancestry_climate_mitotype_CORRECTED.R instead.
## =========================================================
## Follow-up to moduleB_ancestry_vs_winter_climate.R (PC1 vs ancestry,
## aland_excluded: r=-0.502, block-perm p=0.033) and
## moduleB_mitotype_vs_ancestry_climate.R (mitotype vs ancestry: r=+0.498,
## p~=0.035-0.050; mitotype vs PC1: r=-0.297, not significant). Since
## mitotype correlates with BOTH PC1 and ancestry, PC1's ancestry link could
## partly be riding on the same founder-effect structure rather than an
## independent signal -- tested here via partial correlation, with the SAME
## two non-independence-aware nulls as before (Omega-structured, and
## shared-origin block permutation), applied to the PARTIAL correlation
## statistic itself, not the raw one.
##
## Run from the repo root: Rscript module_manuscript_rho05/R/moduleB_PC1_ancestry_partial_mitotype.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

OMEGA_DIR <- "module_manuscript_rho05/baypass_stage1/aland_excluded"
FIGDIR <- "module_manuscript_rho05/Figures"
NSIM <- 10000; NPERM <- 10000
SEED_SIM <- 2026; SEED_PERM <- 2027

obj <- readRDS("module_manuscript_rho05/data/moduleB_mitotype_vs_ancestry_climate.rds")
dt19 <- copy(obj$dt19)

partial_cor <- function(x, y, z) {
  rxy <- cor(x, y); rxz <- cor(x, z); ryz <- cor(y, z)
  (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
}

raw_r     <- cor(dt19$PC1, dt19$ancestry)
partial_r <- partial_cor(dt19$PC1, dt19$ancestry, dt19$mito01)
cat(sprintf("raw cor(PC1, ancestry)                    = %+.3f\n", raw_r))
cat(sprintf("partial cor(PC1, ancestry | mitotype)     = %+.3f\n", partial_r))
cat(sprintf("(cor(PC1,mito)=%+.3f, cor(ancestry,mito)=%+.3f)\n\n", cor(dt19$PC1, dt19$mito01), cor(dt19$ancestry, dt19$mito01)))

## also via residuals, as a cross-check (should match partial_r) -----------
res_pc1 <- residuals(lm(PC1 ~ mito01, data = dt19))
res_anc <- residuals(lm(ancestry ~ mito01, data = dt19))
cat(sprintf("cross-check via residuals: cor(resid PC1, resid ancestry) = %+.3f\n\n", cor(res_pc1, res_anc)))

## ---- Omega-structured null on the PARTIAL statistic -----------------------
Omega <- as.matrix(fread(file.path(OMEGA_DIR, "omega_mat_omega.out"))); Omega <- (Omega + t(Omega)) / 2
eg <- eigen(Omega, symmetric = TRUE); vals <- pmax(eg$values, 0); P <- nrow(Omega)
stopifnot(P == nrow(dt19))
set.seed(SEED_SIM)
NullMat <- sapply(seq_len(NSIM), function(k) as.numeric(scale(eg$vectors %*% (sqrt(vals) * rnorm(P)))))
## for each null draw (standing in for "a structure-consistent PC1"),
## compute its partial correlation with the OBSERVED ancestry, controlling
## for the OBSERVED mitotype
null_partial <- apply(NullMat, 2, function(v) partial_cor(v, dt19$ancestry, dt19$mito01))
k_om <- sum(abs(null_partial) >= abs(partial_r))
p_om <- (1 + k_om) / (NSIM + 1)
cat(sprintf("Omega-structured null on partial r: null median=%.3f, 95%%=[%.3f,%.3f], k=%d, p=%.4f\n",
            median(null_partial), quantile(null_partial, .025), quantile(null_partial, .975), k_om, p_om))

## ---- shared-origin block permutation on the PARTIAL statistic -------------
unit <- dt19$Population
unit[dt19$Population %in% c("LangholmenW", "LangholmenR")] <- "unit_Lang"
unit[dt19$Population %in% c("Bunkkeri", "Grundsund")]      <- "unit_BunGru"
units <- unique(unit)
set.seed(SEED_PERM)
pc1_by_unit <- setNames(dt19$PC1[match(units, unit)], units)
null_partial_perm <- vapply(seq_len(NPERM), function(i) {
  perm_units <- sample(units)
  perm_pc1 <- setNames(pc1_by_unit[perm_units], units)[unit]
  partial_cor(perm_pc1, dt19$ancestry, dt19$mito01)
}, numeric(1))
k_bp <- sum(abs(null_partial_perm) >= abs(partial_r))
p_bp <- (1 + k_bp) / (NPERM + 1)
cat(sprintf("block-permutation null on partial r (%d units): null median=%.3f, 95%%=[%.3f,%.3f], k=%d, p=%.4f\n",
            length(units), median(null_partial_perm), quantile(null_partial_perm, .025),
            quantile(null_partial_perm, .975), k_bp, p_bp))

## ---- also report for comparison: PC2 (the OTHER climate axis) ------------
partial_r_pc2 <- partial_cor(dt19$PC2, dt19$ancestry, dt19$mito01)
cat(sprintf("\n(for reference) partial cor(PC2, ancestry | mitotype) = %+.3f  (raw was %+.3f)\n",
            partial_r_pc2, cor(dt19$PC2, dt19$ancestry)))

summ <- data.table(
  test = c("raw r", "partial r | mitotype", "Omega-null on partial r", "block-perm on partial r"),
  value = c(round(raw_r, 3), round(partial_r, 3), NA, NA),
  p = c(NA, NA, round(p_om, 4), round(p_bp, 4))
)
cat("\n=== summary ===\n"); print(summ)
saveRDS(list(dt19 = dt19, raw_r = raw_r, partial_r = partial_r, p_omega_null = p_om, p_block_perm = p_bp,
             null_partial = null_partial, null_partial_perm = null_partial_perm),
        "module_manuscript_rho05/data/moduleB_PC1_ancestry_partial_mitotype.rds")

## ---- figure: raw vs partial, and null distributions -----------------------
p1 <- ggplot(dt19, aes(PC1, ancestry, colour = Mitotype)) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey60", linewidth = 0.4) +
  geom_point(size = 2) +
  labs(title = sprintf("Raw: r=%.3f", raw_r), x = "PC1", y = "ancestry") +
  theme_classic(base_size = 9) + theme(legend.position = "bottom")
dt19b <- copy(dt19); dt19b[, `:=`(res_pc1 = res_pc1, res_anc = res_anc)]
p2 <- ggplot(dt19b, aes(res_pc1, res_anc, colour = Mitotype)) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey60", linewidth = 0.4) +
  geom_point(size = 2) +
  labs(title = sprintf("Partial (| mitotype): r=%.3f, p_Omega=%.3f, p_block=%.3f", partial_r, p_om, p_bp),
       x = "PC1 residual (after mitotype)", y = "ancestry residual (after mitotype)") +
  theme_classic(base_size = 9) + theme(legend.position = "bottom")
p_combo <- p1 + p2 + plot_layout(guides = "collect") + plot_annotation(
  title = "Does PC1's ancestry link survive controlling for mitotype?")
ggsave(file.path(FIGDIR, "moduleB_PC1_ancestry_partial_mitotype.png"), p_combo, width = 9, height = 4.5, dpi = 200)
cat("\nSaved figure + moduleB_PC1_ancestry_partial_mitotype.rds\n")
