## =========================================================
## MODULE C reformulation: size-normalised climate association
## (removes the size-gated threshold artifact)
## =========================================================
## PROBLEM this fixes. The outlier definition ">= MIN_N_SIG members with
## BF(dB) >= SIG_THR" is an ABSOLUTE COUNT, so eligibility scales with cluster
## size. Measured on this data:
##     9 loci  (background median) -> needs 56% of members significant
##    13 loci  (75th pct)          -> 38%
##    21 loci  (90th pct)          -> 24%
##    64-90    (outlier median)    -> 6-8%
## and 88% of background clusters have <20 loci (needing >25%). So the "same"
## criterion is in fact far stricter for small clusters. Adjusting with
## log(n_loci) assumes a smooth log-linear size effect and cannot represent that
## hard, highly nonlinear ELIGIBILITY threshold -- the model extrapolates the
## outlier coefficient into a size range where outliers structurally cannot occur.
## (At the original 20/10 it was worse: >=10 significant members is literally
## impossible for the 56.5% of clusters with <10 loci.)
##
## FIX (primary): never binarise on a size-gated count. Use the per-cluster
## significance RATE p_sig = n_sig / n_loci -- size-normalised by construction --
## as the measure of climate-association strength, keeping log(n_loci) as a
## covariate for the independent size->DI / size->sorting effects. Also report the
## binomial EXCESS z (observed n_sig vs its size expectation) as an alternative
## parameterisation, and the old binary is_outlier side by side, so we can see
## exactly how much of the result depended on the gated count.
##
## ROBUSTNESS: repeat restricted to clusters where outlier status is realistically
## attainable (n_loci >= 20, >= 50), i.e. enforcing overlap rather than
## extrapolating.
##
## Reads: canonical clustering, the 8 BayPass BF files, map (marker order for the
## positional BF join + DI), and the CACHED consensus sort_class
## (data/moduleC_C3_cl.rds). No consensus rebuild, no hybrids_only load.

suppressMessages({ library(data.table) })

SIG_THR <- 15; MIN_N_SIG <- 5      # must match moduleC_sorting_climate.R / C2
DI_TH   <- -25                      # C2's high-DI definition
CLUSTERING <- "data/eMLG_5loci_0025_cM05.rds"
CKPT <- "data/moduleC_C3_cl.rds"
uni <- c("aquilonia", "polyctena")

## ---- inputs -------------------------------------------------------------
e2 <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
map <- e2$map_hyb_005                       # marker order (positional BF join) + DI
rm(e2); invisible(gc())
groups <- readRDS(CLUSTERING)$groups
cl     <- readRDS(CKPT)                     # cached consensus sort_class per cluster
message("[rate] inputs loaded")

he  <- groups[has_eMLG == TRUE]
m2g <- he[, .(marker = unlist(members)), by = group_id]
m2g[, DI := setNames(map$DiagnosticIndex, map$marker)[marker]]

import_bf <- function(f) {
  r <- fread(f)
  stopifnot(nrow(r) == nrow(map), identical(r$MRK, seq_len(nrow(r))))
  setNames(r$`BF(dB)`, map$marker)
}
cfg <- CJ(popset = c("with_aland", "aland_excluded"), pc = c("PC1", "PC2"),
          omega = c("noOmega", "withOmega"))

## per-cluster DI counts (constant across configs)
di_cl <- m2g[, .(n_hi = sum(DI > DI_TH, na.rm = TRUE), n_di = sum(!is.na(DI))), by = group_id]

## ---- per-config cluster table: n_sig, rate, excess z, binary flag --------
cluster_table <- function(bf) {
  m2g[, BF := bf[marker]]
  a <- m2g[, .(n_loci = .N, n_sig = sum(BF >= SIG_THR, na.rm = TRUE)), by = group_id]
  p0 <- mean(m2g$BF >= SIG_THR, na.rm = TRUE)          # background significance rate
  a[, `:=`(p_sig = n_sig / n_loci,
           z_excess = (n_sig - n_loci * p0) / sqrt(n_loci * p0 * (1 - p0)),
           is_outlier = as.integer(n_sig >= MIN_N_SIG))]
  a <- di_cl[a, on = "group_id"]
  cl[, .(group_id, differentiated, directional, sort_class)][a, on = "group_id"]
}

## ---- the two tests, each with the three predictor parameterisations ------
## DI: CLUSTER-LEVEL -- each cluster is ONE observation. A marker-level binomial
## treats LD-correlated members as independent trials (dispersion ~64 on these
## data), understating the SE ~8-fold; that is precisely the pseudo-replication
## this pipeline collapses everywhere else (Module A). PRIMARY = unweighted
## cluster-level linear model on the high-DI fraction. CONSERVATIVE CHECK =
## size-weighted quasibinomial (overdispersion-corrected).
fit_di <- function(D, pred) {
  d <- copy(D[n_di > 0])[, frac_hi := n_hi / n_di]
  m1 <- lm(as.formula(sprintf("frac_hi ~ %s + log(n_loci)", pred)), data = d)
  c1 <- summary(m1)$coefficients
  m2 <- glm(as.formula(sprintf("cbind(n_hi, n_di - n_hi) ~ %s + log(n_loci)", pred)),
            data = d, family = quasibinomial)
  c2 <- summary(m2)$coefficients
  rbind(
    data.table(test = "DI cluster-level", predictor = pred,
               est = c1[pred, "Estimate"], p = c1[pred, "Pr(>|t|)"], n = nrow(d)),
    data.table(test = "DI quasibinomial", predictor = pred,
               est = c2[pred, "Estimate"], p = c2[pred, "Pr(>|t|)"], n = nrow(d))
  )
}
## sorting: cluster-level logistic on directional (differentiated clusters only)
fit_sort <- function(D, pred) {
  d <- D[differentiated == TRUE & !is.na(directional)]
  f <- as.formula(sprintf("directional ~ %s + log(n_loci)", pred))
  m <- glm(f, data = d, family = binomial); cf <- summary(m)$coefficients
  data.table(test = "sorting-overlap", predictor = pred,
             est = cf[pred, "Estimate"], p = cf[pred, "Pr(>|z|)"], n = nrow(d))
}
## report on a common scale: OR per +10 percentage points of significance rate
## (p_sig), per +1 SD of excess z, and per the binary flag as-is.
## Common reporting scale. The cluster-level DI model is linear in the high-DI
## FRACTION, so its effect is reported in PERCENTAGE POINTS; the logistic /
## quasibinomial models are reported as ODDS RATIOS. Both are expressed per +10
## percentage points of significance rate (p_sig), or per unit for the others.
scale_est <- function(dt) {
  lin <- dt$test == "DI cluster-level"
  dt[, effect := fifelse(lin,
                         fifelse(predictor == "p_sig", est * 10, est * 100),
                         fifelse(predictor == "p_sig", exp(est * 0.1), exp(est)))]
  dt[, scale := fifelse(lin, "pp high-DI", "odds ratio")]
  dt[, per := fifelse(predictor == "p_sig", "per +10pp sig rate",
              fifelse(predictor == "z_excess", "per +1 excess z", "outlier vs not"))]
  dt[, .(test, predictor, per, scale, effect = round(effect, 3), p = signif(p, 3), n)]
}

run_config <- function(i, min_nloci = 0) {
  c <- cfg[i]
  bf <- import_bf(sprintf("./%s/%s_DIEM_%s_summary_betai_reg.out", c$popset, c$pc, c$omega))
  D <- cluster_table(bf)[n_loci >= min_nloci]
  res <- rbindlist(lapply(c("p_sig", "z_excess", "is_outlier"),
                          function(p) rbind(fit_di(D, p), fit_sort(D, p))))
  cbind(c[, .(popset, pc, omega)], scale_est(res))
}

## ========================================================================
## 1 -- full set (all has_eMLG clusters)
## ========================================================================
message("[rate] fitting all 8 configs, full set ...")
full <- rbindlist(lapply(seq_len(nrow(cfg)), run_config))
cat("\n=== all has_eMLG clusters: size-normalised association (p_sig) ===\n")
print(full[predictor == "p_sig", .(popset, pc, omega, test, scale, effect, p)])
cat("\n--- PRIMARY config (aland_excluded x withOmega), with p-values ---\n")
print(full[popset == "aland_excluded" & omega == "withOmega"])

## ========================================================================
## 2 -- ROBUSTNESS: restrict to clusters where outlier status is attainable
## ========================================================================
for (mn in c(20, 50)) {
  message(sprintf("[rate] robustness: n_loci >= %d ...", mn))
  rb <- rbindlist(lapply(seq_len(nrow(cfg)), run_config, min_nloci = mn))
  cat(sprintf("\n=== ROBUSTNESS: clusters with n_loci >= %d only ===\n", mn))
  print(rb[popset == "aland_excluded" & omega == "withOmega"])
  assign(sprintf("rob%d", mn), rb)
}

## ========================================================================
## 3 -- dose-response figure: ancestry-diagnostic content vs climate-association
##      strength, in size-comparable strata (primary config)
## ========================================================================
suppressMessages(library(ggplot2))
message("[rate] building dose-response figure ...")
prim <- cfg[popset == "aland_excluded" & omega == "withOmega"]
dr <- rbindlist(lapply(seq_len(nrow(prim)), function(i) {
  c <- prim[i]
  bf <- import_bf(sprintf("./%s/%s_DIEM_%s_summary_betai_reg.out", c$popset, c$pc, c$omega))
  cluster_table(bf)[, pc := c$pc][]
}))
gw <- 100 * mean(map$DiagnosticIndex > DI_TH, na.rm = TRUE)   # genome-wide baseline

## p_sig is heavily zero-inflated (most clusters have no member at BF >= SIG_THR),
## so quantile breaks over all clusters collapse to a single bin. Bin the zeros
## separately ("no climate signal"), then split the positive rates by quantile.
bin_dr <- function(d, lab, nb = 5) {
  d <- d[n_di > 0 & is.finite(p_sig)]
  d[, bin := NA_character_]
  d[p_sig == 0, bin := "0"]
  pos <- d[p_sig > 0]
  if (nrow(pos) > nb) {
    br <- unique(quantile(pos$p_sig, 0:nb / nb, na.rm = TRUE))
    if (length(br) > 1) d[p_sig > 0, bin := as.character(cut(p_sig, br, include.lowest = TRUE))]
  }
  d[!is.na(bin), .(x = 100 * median(p_sig), y = 100 * sum(n_hi) / sum(n_di), n = .N),
    by = .(pc, bin)][order(x)][, stratum := lab][]
}
dr_plot <- rbind(bin_dr(dr[n_loci >= 20], "clusters >= 20 loci"),
                 bin_dr(dr[n_loci >= 50], "clusters >= 50 loci"))

dir.create("Figures", showWarnings = FALSE)
p_dr <- ggplot(dr_plot, aes(x, y, colour = stratum)) +
  geom_hline(yintercept = gw, linetype = 2, colour = "grey60") +
  annotate("text", x = Inf, y = gw, label = "genome-wide", hjust = 1.05, vjust = -0.6,
           size = 2.3, colour = "grey45") +
  geom_line(linewidth = 0.5) + geom_point(aes(size = n)) +
  scale_colour_manual(values = c("clusters >= 20 loci" = "#0072B2",
                                 "clusters >= 50 loci" = "#D55E00")) +
  scale_size_continuous(range = c(1, 3.2), guide = "none") +
  facet_wrap(~ pc) +
  labs(x = expression("climate association: % of member loci with BF(dB) " >= " 15"),
       y = "% of member loci with DI > -25", colour = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7), strip.text = element_text(size = 8),
        axis.title = element_text(size = 8))
ggsave(sprintf("Figures/moduleC_dose_response_%d_%d.pdf", MIN_N_SIG, SIG_THR),
       p_dr, width = 150, height = 78, units = "mm")
ggsave(sprintf("Figures/moduleC_dose_response_%d_%d.png", MIN_N_SIG, SIG_THR),
       p_dr, width = 150, height = 78, units = "mm", dpi = 300)

saveRDS(list(full = full, rob20 = rob20, rob50 = rob50, dose_response = dr_plot),
        sprintf("data/moduleC_rate_based_%d_%d.rds", MIN_N_SIG, SIG_THR))
cat(sprintf("\nSaved: data/moduleC_rate_based_%d_%d.rds, Figures/moduleC_dose_response_%d_%d.pdf/.png\n",
            MIN_N_SIG, SIG_THR, MIN_N_SIG, SIG_THR))
