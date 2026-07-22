library(data.table)
library(ggplot2)

# ------------------------------------------------------------
# DiagnosticIndex enrichment in BayPass outlier association regions
# ------------------------------------------------------------
# Checks whether the "outlier" LD clusters identified from BayPass PC1/PC2
# association results (see R/analyse_baypass.R: clusters with >=10
# individually significant, BF(dB)>=20, member loci) are enriched for
# high-DiagnosticIndex loci relative to the rest of the (eMLG-filtered)
# genome.
#
# DiagnosticIndex is DIEM's per-marker ancestry-informativeness statistic
# (attached to map_hyb_005 during parsing, see R/LD_decay_from_DIEM.R).
# DiagnosticIndex > -25 is the same threshold already used elsewhere in
# this pipeline (R/LD_decay_from_DIEM.R:104).
#
# IMPORTANT confound, checked and corrected for here: the outlier
# definition REQUIRES n_loci >= 10 by construction (can't have 10
# significant members in a smaller cluster), so outlier clusters are
# mechanically size-biased relative to the full background. Separately,
# cluster size itself is associated with DiagnosticIndex genome-wide
# (independent of any BayPass result -- checked directly: Spearman r=0.23
# between cluster size and % high-DiagnosticIndex members; marker-level
# logistic regression OR=1.135 per doubling of cluster size, p~0),
# consistent with larger non-recombining regions being more likely to
# harbor multiple ancestry-diagnostic loci (and plausibly Dobzhansky-Muller
# incompatibilities) in the first place. So part of the naive outlier-vs-
# background enrichment could be a pure size artifact rather than evidence
# that BayPass is finding real signal.
#
# A simple size-floor match (background restricted to n_loci>=10 too)
# turns out NOT to be enough: outlier clusters are still far larger even
# within that restricted pool (checked directly on one config: outlier
# median 147 loci / mean 743, vs. size-matched background median 15 /
# mean 28 -- a ~25-30x gap in typical size that a >=10 floor alone doesn't
# remove). The primary analysis here instead uses a marker-level logistic
# regression with log(n_loci) as a continuous covariate alongside outlier
# status, which properly accounts for this rather than just truncating the
# small end of the size distribution.

message("=== Loading data ===")
load("./data/hybrids_only_maf005.Rdata")
DI <- map_hyb_005[TRUE, .(marker, DiagnosticIndex)]

eMLG_result <- readRDS("./data/eMLG_5loci_0025_cM05.rds")
has_eMLG_groups <- eMLG_result$groups[has_eMLG == TRUE]
marker_group <- has_eMLG_groups[TRUE, .(marker = unlist(members), n_loci = n_loci), by = group_id]
marker_group <- DI[marker_group, on = "marker"]

import_bf <- function(file, n_expected) {
  res <- fread(file)
  stopifnot(
    "summary_betai_reg.out row count doesn't match map_hyb_005 -- MRK-based positional join would silently misalign" =
      nrow(res) == n_expected,
    "MRK isn't 1:N in order -- positional join assumption doesn't hold for this file" =
      identical(res$MRK, seq_len(nrow(res)))
  )
  res$`BF(dB)`
}

sig_threshold <- 20
min_n_sig_loci <- 10
di_threshold <- -25

# ------------------------------------------------------------
# Step 0: is cluster size associated with DiagnosticIndex at all, genome-
# wide, independent of any BayPass result? (confirms the confound is real
# before worrying about correcting for it)
# ------------------------------------------------------------
marker_group[, hi := DiagnosticIndex > di_threshold]
m0 <- glm(hi ~ log(n_loci), data = marker_group, family = binomial)
size_effect_or <- unname(exp(coef(m0)["log(n_loci)"] * log(2)))
message(sprintf(
  "Cluster-size effect on DiagnosticIndex (genome-wide, independent of BayPass): OR = %.3f per doubling of cluster size, p = %.2e",
  size_effect_or, summary(m0)$coefficients["log(n_loci)", "Pr(>|z|)"]
))

# ------------------------------------------------------------
# Per-config enrichment: naive (unadjusted) vs. size-adjusted
# ------------------------------------------------------------
enrichment_for_config <- function(population_set, pc, omega) {
  file <- paste0("./", population_set, "/", pc, "_DIEM_", omega, "_summary_betai_reg.out")
  bf <- import_bf(file, nrow(map_hyb_005))

  dt <- copy(marker_group)
  dt[, BF := bf[match(marker, DI$marker)]]

  sig_ids <- dt[, .(n_sig = sum(BF >= sig_threshold, na.rm = TRUE)), by = group_id][n_sig >= min_n_sig_loci, group_id]
  dt[, is_outlier := as.integer(group_id %in% sig_ids)]
  outlier <- dt[is_outlier == 1L]
  background <- dt[is_outlier == 0L]

  ## naive: raw 2x2 Fisher's exact test, outlier vs. ALL background
  n_out <- nrow(outlier); n_out_hi <- sum(outlier$DiagnosticIndex > di_threshold, na.rm = TRUE)
  n_bg <- nrow(background); n_bg_hi <- sum(background$DiagnosticIndex > di_threshold, na.rm = TRUE)
  ft <- fisher.test(matrix(c(n_out_hi, n_out - n_out_hi, n_bg_hi, n_bg - n_bg_hi), nrow = 2))

  ## size-adjusted: marker-level logistic regression, outlier status +
  ## log(cluster size) -- the outlier coefficient is the enrichment that
  ## survives controlling for cluster size continuously
  m <- glm(hi ~ is_outlier + log(n_loci), data = dt, family = binomial)
  cf <- summary(m)$coefficients
  ## Wald CI (confint.default), not profile-likelihood confint() -- the
  ## latter is far too slow on ~500k-row models and the two agree closely
  ## at this sample size anyway.
  ci <- confint.default(m, "is_outlier")

  data.table(
    population_set = population_set, pc = pc, omega = omega,
    n_clusters = uniqueN(sig_ids), n_outlier_markers = n_out,
    pct_outlier = 100 * n_out_hi / n_out, pct_background = 100 * n_bg_hi / n_bg,
    naive_OR = unname(ft$estimate), naive_OR_lo = ft$conf.int[1], naive_OR_hi = ft$conf.int[2], naive_p = ft$p.value,
    adj_OR = unname(exp(cf["is_outlier", "Estimate"])),
    adj_OR_lo = exp(ci[1]), adj_OR_hi = exp(ci[2]),
    adj_p = cf["is_outlier", "Pr(>|z|)"]
  )
}

message("=== Computing naive and size-adjusted enrichment across all available result sets ===")
configs <- CJ(
  population_set = c("with_aland", "aland_excluded"),
  pc = c("PC1", "PC2"),
  omega = c("noOmega", "withOmega"),
  sorted = FALSE
)
enrichment <- configs[, enrichment_for_config(population_set, pc, omega), by = .I]
enrichment[, I := NULL]
enrichment[, naive_sig := fifelse(naive_p < 0.001, "***", fifelse(naive_p < 0.01, "**", fifelse(naive_p < 0.05, "*", "n.s.")))]
enrichment[, adj_sig := fifelse(adj_p < 0.001, "***", fifelse(adj_p < 0.01, "**", fifelse(adj_p < 0.05, "*", "n.s.")))]

genome_wide_pct <- 100 * mean(DI$DiagnosticIndex > di_threshold, na.rm = TRUE)

dir.create("./data/", showWarnings = FALSE)
fwrite(enrichment, "./data/diagnostic_index_enrichment.csv")
message("Saved: ./data/diagnostic_index_enrichment.csv")
print(enrichment[, .(population_set, pc, omega, n_clusters, n_outlier_markers,
                      pct_outlier, pct_background, naive_OR, naive_sig, adj_OR, adj_sig)])

# ------------------------------------------------------------
# Figure: naive vs. size-adjusted odds ratio, per config -- the point of
# this figure IS the correction, so both estimates are shown side by side
# rather than only the "final" one.
# ------------------------------------------------------------
enrichment[, config_label := paste(pc, population_set, omega)]
enrichment[, config_label := factor(config_label, levels = rev(config_label))]

plot_dt <- rbindlist(list(
  enrichment[, .(config_label, method = "Naive (unadjusted)", or = naive_OR, or_lo = naive_OR_lo, or_hi = naive_OR_hi, p = naive_p)],
  enrichment[, .(config_label, method = "Size-adjusted", or = adj_OR, or_lo = adj_OR_lo, or_hi = adj_OR_hi, p = adj_p)]
))
plot_dt[, method := factor(method, levels = c("Naive (unadjusted)", "Size-adjusted"))]

p_forest <- ggplot(plot_dt, aes(or, config_label, color = method)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey40") +
  geom_errorbar(aes(xmin = pmax(or_lo, 0.1), xmax = pmin(or_hi, 100)),
                position = position_dodge(width = 0.5), width = 0.25, orientation = "y") +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  scale_x_log10() +
  scale_color_manual(values = c("Naive (unadjusted)" = "grey55", "Size-adjusted" = "#a83a2b"), name = NULL) +
  labs(
    x = "Odds ratio (outlier cluster vs. background), log scale", y = NULL,
    title = "DiagnosticIndex enrichment: naive vs. cluster-size-adjusted",
    subtitle = "Dashed line: no enrichment (OR=1).\nSize-adjusted = marker-level logistic regression controlling for log(cluster size).",
    caption = "Naive test ignores that outlier clusters are mechanically larger (require >=10 significant members). Adjusted test controls for this."
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(size = 13), plot.caption = element_text(hjust = 0, size = 9))

ggsave("./Figures/diagnostic_index_enrichment_forest.png", p_forest, width = 9.5, height = 5.5, dpi = 150)
message("Saved: ./Figures/diagnostic_index_enrichment_forest.png")

# ------------------------------------------------------------
# Figure: proportions -- outlier cluster vs. background range, per config
# (unchanged from before; still useful as the plain-language "how common
# is this" view, alongside the OR-based confound correction above)
# ------------------------------------------------------------
background_range <- range(enrichment$pct_background)

p_proportions <- ggplot(enrichment, aes(config_label, pct_outlier)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = background_range[1], ymax = background_range[2],
           fill = "grey70", alpha = 0.4) +
  geom_hline(yintercept = genome_wide_pct, linetype = 3, color = "grey30", linewidth = 0.5) +
  geom_col(fill = "#a83a2b", width = 0.6) +
  geom_text(aes(label = naive_sig, y = pct_outlier + 3), size = 4) +
  coord_flip() +
  labs(
    x = NULL, y = "% of outlier-cluster markers with DiagnosticIndex > -25",
    title = "Raw proportions: outlier clusters vs. background",
    subtitle = "Grey band: background range across configs (7.4-8.0%).\nDotted line: genome-wide baseline (4.5%).",
    caption = paste0(
      "Outlier cluster = LD cluster with >=10 individually significant (BF(dB)>=20) member loci.\n",
      "Significance stars are from the NAIVE (unadjusted) test -- see the size-adjusted forest plot for the confound-corrected result."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 13), plot.caption = element_text(hjust = 0, size = 9))

ggsave("./Figures/diagnostic_index_enrichment_proportions.png", p_proportions, width = 9, height = 6, dpi = 150)
message("Saved: ./Figures/diagnostic_index_enrichment_proportions.png")
