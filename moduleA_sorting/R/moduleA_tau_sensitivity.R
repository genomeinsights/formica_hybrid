## =========================================================
## MODULE A -- sorting-threshold (tau) sensitivity series
## =========================================================
## Predefined series tau in {0.5 relaxed, 0.6 primary, 0.8 stringent}; everything
## else fixed (fix_th 0.15, min_parent_maf 0.15, cM05 clustering, sort_rule "binom"
## at alpha 0.05, DI ungated). Because the per-unit fixation counts are
## tau-INDEPENDENT, this is a cheap post-hoc reclassification of stored counts -- no
## re-run of the expensive prep. It reports, across the three tau, exactly the four
## questions the manuscript asks:
##   * does the DIRECTION / broad pattern persist?          (direction ratio, zDI sign)
##   * do EFFECT SIZES stay comparable?                     (direction glm zDI/zr + CI)
##   * does only PREVALENCE change, as expected?            (sort_class counts)
##   * does anything appear ONLY at the relaxed threshold?  (unresolved prevalence)
## tau = 0.6 is the reported analysis; 0.5 / 0.8 show what is gained / lost. We never
## pick the strongest tau. Results that appear only at tau = 0.5 (e.g. the
## unresolved peaks) are threshold-sensitive / exploratory; results surviving
## tau = 0.8 are the strongest evidence (on a smaller set, so wider CIs).
##
## Reads : moduleA_sorting/data/moduleA_snp.rds (per-SNP counts),
##         moduleA_sorting/data/moduleA_r_unit.rds (LD-reduced unit counts +
##         recomb/n_loci/parent_maf, from moduleA_architecture.R),
##         moduleA_sorting/data/moduleA_cluster_sorting_counts.rds (has_eMLG counts,
##         if moduleA_cluster_sorting.R has been re-run).
## Writes: moduleA_sorting/data/moduleA_tau_sensitivity.rds,
##         moduleA_sorting/Figures/moduleA_tau_sensitivity.png
## Run from the repo root.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("moduleA_sorting/R/parallelism_stats.R")   # classify_sort(), MODULEA_TAU_SERIES, tau_stamp()
RULE <- "binom"; ALPHA <- 0.05
CLASSLV <- c("aquilonia", "polyctena", "unresolved", "ambiguous", "unsorted")

## ---- inputs --------------------------------------------------------------
snp  <- as.data.table(readRDS("moduleA_sorting/data/moduleA_snp.rds"))
unit <- as.data.table(readRDS("moduleA_sorting/data/moduleA_r_unit.rds")$r)  # cache = list(params, r)
emlg_f <- "moduleA_sorting/data/moduleA_cluster_sorting_counts.rds"
emlg <- if (file.exists(emlg_f)) as.data.table(readRDS(emlg_f)) else NULL

## z-covariates for the direction / magnitude models (match moduleA_architecture.R)
du <- unit[differentiated == TRUE & is.finite(recomb) & is.finite(DI)]
du[, `:=`(zr   = as.numeric(scale(log10(recomb + 0.1))),
          zDI  = as.numeric(scale(DI)),
          zmaf = as.numeric(scale(parent_maf)),
          zcs  = as.numeric(scale(log10(n_loci))))]

## ---- prevalence across tau (SNP / LD-reduced unit / has_eMLG) -------------
prev_one <- function(dt, tau, level) {
  cl <- classify_sort(dt$n_aqu, dt$n_pol, dt$n_obs, tau, RULE, ALPHA)
  as.data.table(as.list(table(factor(cl, levels = CLASSLV))))[
    , `:=`(level = level, tau = tau_stamp(tau))][]
}
res_prev <- rbindlist(c(
  lapply(MODULEA_TAU_SERIES, function(t) prev_one(snp[differentiated == TRUE & n_obs > 0], t, "SNP")),
  lapply(MODULEA_TAU_SERIES, function(t) prev_one(du, t, "unit")),
  if (!is.null(emlg)) lapply(MODULEA_TAU_SERIES, function(t) prev_one(emlg[differentiated == TRUE & n_obs > 0], t, "eMLG"))
), fill = TRUE)
setcolorder(res_prev, c("level", "tau"))
res_prev[, `:=`(unidir = aquilonia + polyctena, aqu_pol_ratio = round(aquilonia / polyctena, 2))]

## ---- tau-INDEPENDENT magnitude slope (reported once) ---------------------
mag <- summary(lm(prop_fixed ~ zr + zDI, data = du))$coefficients

## ---- effect sizes across tau: direction glm + low-recomb depletion -------
du[, rdec := cut(recomb, quantile(recomb, 0:10 / 10, na.rm = TRUE), include.lowest = TRUE, labels = FALSE)]
res_eff <- rbindlist(lapply(MODULEA_TAU_SERIES, function(tau) {
  du[, sc := classify_sort(n_aqu, n_pol, n_obs, tau, RULE, ALPHA)]
  cu <- du[sc %in% c("aquilonia", "polyctena")][, is_aqu := as.integer(sc == "aquilonia")]
  g  <- glm(is_aqu ~ zDI + zr + zmaf + zcs, data = cu, family = binomial())
  co <- as.data.table(summary(g)$coefficients, keep.rownames = "term")
  setnames(co, c("term", "est", "se", "z", "p"))
  ci <- confint.default(g)
  du[, sorted := !is.na(sc) & sc != "unsorted"]
  data.table(tau = tau_stamp(tau), n_unidir = nrow(cu),
             zDI = co[term == "zDI", est], zDI_z = co[term == "zDI", z],
             zDI_lo = ci["zDI", 1], zDI_hi = ci["zDI", 2],
             zr = co[term == "zr", est], zr_z = co[term == "zr", z],
             DI_vs_recomb = abs(co[term == "zDI", est] / co[term == "zr", est]),
             sorted_dec1 = du[rdec == 1, mean(sorted)],
             sorted_rest = du[rdec > 1, mean(sorted)])
}))
res_eff[, `:=`(depletion = round(sorted_rest - sorted_dec1, 3))]

cat("=== PREVALENCE across tau (only prevalence should change) ===\n")
print(res_prev[, .(level, tau, aquilonia, polyctena, unresolved, ambiguous, unsorted, aqu_pol_ratio)])
cat(sprintf("\ntau-independent magnitude slope prop_fixed ~ recomb (zr): %+.3f (t = %.0f)\n",
            mag["zr", "Estimate"], mag["zr", "t value"]))
cat("\n=== EFFECT SIZES across tau (direction glm on unidirectional units; low-recomb depletion) ===\n")
print(res_eff[, .(tau, n_unidir, zDI = round(zDI, 2), zDI_CI = sprintf("[%.2f,%.2f]", zDI_lo, zDI_hi),
                  zDI_z = round(zDI_z), zr = round(zr, 2), DI_vs_recomb = round(DI_vs_recomb),
                  sorted_dec1 = round(sorted_dec1, 3), sorted_rest = round(sorted_rest, 3))])
cat("\nRELAXED-ONLY (threshold-sensitive / exploratory): unresolved prevalence",
    paste(sprintf("%s=%d", res_prev[level == "SNP", tau], res_prev[level == "SNP", unresolved]), collapse = " -> "), "\n")

saveRDS(list(prevalence = res_prev, effects = res_eff, mag_slope = mag,
             params = list(series = MODULEA_TAU_SERIES, primary = MODULEA_TAU_PRIMARY,
                           rule = RULE, alpha = ALPHA, fix_th = 0.15, min_parent_maf = 0.15)),
        "moduleA_sorting/data/moduleA_tau_sensitivity.rds")

## ---- figure: prevalence (SNP) + direction effect size, across tau --------
pv <- melt(res_prev[level == "SNP"], id.vars = "tau",
           measure.vars = c("aquilonia", "polyctena", "unresolved", "unsorted"),
           variable.name = "class", value.name = "n")
pv[, class := factor(class, levels = c("aquilonia", "polyctena", "unresolved", "unsorted"))]
cols <- c(aquilonia = "#2c7fb8", polyctena = "#d95f0e", unresolved = "#31a354", unsorted = "#bdbdbd")
pA <- ggplot(pv, aes(tau, n, fill = class)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = cols) + scale_y_continuous(labels = scales::comma) +
  labs(title = "A. Per-SNP prevalence by threshold", x = NULL, y = "differentiated SNPs", fill = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "top")
pB <- ggplot(res_eff, aes(tau, zDI)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_errorbar(aes(ymin = zDI_lo, ymax = zDI_hi), width = .15) +
  geom_point(size = 3, colour = "#2166AC") +
  geom_text(aes(label = sprintf("z=%.0f, n=%s", zDI_z, formatC(n_unidir, big.mark = ",", format = "d"))),
            vjust = -1.1, size = 3) +
  labs(title = "B. Direction effect size: DI -> P(aquilonia)", subtitle = "logistic coefficient (95% CI) on unidirectional units",
       x = expression(tau), y = "zDI coefficient") +
  theme_minimal(base_size = 10)
ggsave("moduleA_sorting/Figures/moduleA_tau_sensitivity.png", pA / pB, width = 8, height = 8, dpi = 150)
cat("\nSaved: moduleA_sorting/data/moduleA_tau_sensitivity.rds, moduleA_sorting/Figures/moduleA_tau_sensitivity.png\n")
