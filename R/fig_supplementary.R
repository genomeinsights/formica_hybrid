## =========================================================
## SUPPLEMENTARY FIGURES (restyled to the manuscript palette)
## =========================================================
## Regenerates the key supplementary figures from saved analysis outputs, in the
## same palette/theme as manuscript/figures/main_ancestry_sorting. Non-key figures
## (LD decay, ROC, tracks, stage1-vs-combined, fidelity hist, consensus-vs-rep,
## ancestry confound, manhattans) are copied as-is elsewhere.
##
## Reads : data/moduleA_snp.rds, data/moduleA_di_asymmetry.rds,
##         data/moduleB_architecture.rds, data/moduleC_rate_based_5_15.rds,
##         data/moduleC_maf_power_sensitivity.rds, data/sensitivity_block_bootstrap.rds
## Writes: manuscript/figures/figS07_sorting.{pdf,png},
##         figS08_architecture.*, figS11_dose_response.*, figS12_maf_power.*,
##         figS14_block_bootstrap.*
## Run from the repo root. Re-plots existing estimates only.

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

Accent <- "#315B7D"; MidBlue <- "#8EB4CE"; Muted <- "#5B6570"
Aqu <- "#3A7D77"; Pol <- "#D08A45"; Grey <- "#C7CDD2"
th <- theme_classic(base_size = 9) +
  theme(plot.tag = element_text(face = "bold", size = 11, colour = Accent),
        plot.title = element_text(size = 10, face = "bold", colour = Accent),
        axis.title = element_text(size = 8), legend.position = "none",
        legend.title = element_blank(), legend.text = element_text(size = 6.5),
        legend.key.size = unit(3, "mm"),
        legend.background = element_rect(fill = scales::alpha("white", 0.7), colour = NA),
        plot.margin = margin(4, 6, 4, 4))
sv <- function(f, p, w, h) { ggsave(sprintf("manuscript/figures/%s.pdf", f), p, width = w, height = h, units = "mm")
  ggsave(sprintf("manuscript/figures/%s.png", f), p, width = w, height = h, units = "mm", dpi = 300) }
dir.create("manuscript/figures", showWarnings = FALSE, recursive = TRUE)

## ========================================================================
## figS07 -- Module A: proportions, DI-decile direction, size washout
## ========================================================================
# NOTE: the Module A sorting figure (proportions / DI-decile / size washout) was
# removed from the manuscript set: its content is in main figure panels A-B and in
# supplementary Table S3 (size aggregation + dilution). moduleA_fig1 remains in /doc.

## ========================================================================
## figS08 -- Module B architecture: differentiation deciles + direction
## (the SNP-vs-unit panel was dropped: it duplicates main figure panel C)
## ========================================================================
b <- readRDS("data/moduleB_architecture.rds")
arch <- as.data.table(b$arch_tab)
zsc <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
cha <- melt(arch[, .(med_recomb, Fst = zsc(Fst), dxy = zsc(dxy), pi_within = zsc(pi_within))],
            id.vars = "med_recomb", variable.name = "metric", value.name = "z")
s8a <- ggplot(cha, aes(med_recomb, z, colour = metric)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) + scale_x_log10() +
  scale_colour_manual(values = c(Fst = Pol, dxy = Aqu, pi_within = Accent),
                      labels = c(Fst = expression(F[ST]), dxy = expression(d[xy]), pi_within = expression(pi))) +
  labs(title = "A", x = "recombination (cM/Mb, log)", y = "standardised (z)") +
  th + theme(legend.position = c(0.72, 0.20), legend.direction = "horizontal")
s8c <- ggplot(as.data.table(b$unit_by_recomb), aes(med_r, frac_aqu_of_unidir)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = Muted) +
  geom_line(linewidth = 0.6, colour = Muted) + geom_point(size = 1.7, colour = Accent) +
  scale_x_log10() + ylim(0, 1) +
  labs(title = "B", x = "recombination (cM/Mb, log)", y = "fraction toward aquilonia (unit)") + th
sv("figS08_architecture", (s8a | s8c), 120, 62)

## ========================================================================
## figS11 -- climate dose-response: DI content vs climate rate, by size stratum
## ========================================================================
dr <- as.data.table(readRDS("data/moduleC_rate_based_5_15.rds")$dose_response)
s11 <- ggplot(dr, aes(x, y, colour = stratum)) +
  geom_line(linewidth = 0.5) + geom_point(aes(size = n)) +
  scale_colour_manual(values = c("clusters >= 20 loci" = MidBlue, "clusters >= 50 loci" = Accent)) +
  scale_size_continuous(range = c(1, 3.2), guide = "none") + facet_wrap(~ pc) +
  labs(x = "climate association: % of members with BF(dB) >= 15",
       y = "% of members with DI > -25", colour = NULL) +
  th + theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
sv("figS11_dose_response", s11, 150, 82)

## ========================================================================
## figS12 -- climate MAF/recomb/XtX sensitivity ladder
## ========================================================================
mp <- readRDS("data/moduleC_maf_power_sensitivity.rds")
ld <- as.data.table(mp$ladder); bo <- as.data.table(mp$m3_boot)
ld[, model := factor(model, levels = c("M0", "M1", "M2", "M3"))]
labdt <- unique(ld[, .(model, added)]); labs <- setNames(sprintf("%s\n%s", labdt$model, labdt$added), as.character(labdt$model))
s12 <- ggplot(ld, aes(model, effect)) +
  geom_hline(yintercept = 0, linetype = 2, colour = Muted) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, colour = MidBlue) +
  geom_errorbar(data = bo, aes(x = "M3", ymin = boot_lo, ymax = boot_hi), width = 0.32,
                colour = Accent, linewidth = 0.6, inherit.aes = FALSE) +
  geom_point(size = 1.8) + facet_wrap(~ pc) + scale_x_discrete(labels = labs) +
  labs(x = NULL, y = "enrichment: pp high-DI per +10 pp climate rate",
       title = "Blue = naive 95% CI; dark (M3) = chromosome block-bootstrap 95% CI") +
  th + theme(plot.title = element_text(size = 7), strip.text = element_text(face = "bold"))
sv("figS12_maf_power", s12, 150, 82)

## ========================================================================
## figS14 -- block-bootstrap CIs for the five central coefficients
## ========================================================================
res <- as.data.table(readRDS("data/sensitivity_block_bootstrap.rds")$results)
lvl <- c("model-based", "chromosome", "10 cM")
fig <- melt(res, id.vars = "target", measure.vars = patterns("_lo$", "_hi$"),
            variable.name = "method", value.name = c("lo", "hi"))
fig[, method := factor(lvl[method], levels = lvl)][, est := res$estimate[match(target, res$target)]]
s14 <- ggplot(fig, aes(est, method)) +
  geom_vline(xintercept = 0, linetype = 2, colour = Muted) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = method), height = 0.25, linewidth = 0.6) +
  geom_point(size = 1.5) + facet_wrap(~ target, scales = "free_x", ncol = 1) +
  scale_colour_manual(values = setNames(c(Grey, Accent, Pol), lvl), guide = "none") +
  labs(x = "coefficient (published scale)", y = NULL) +
  th + theme(strip.text = element_text(size = 8, face = "bold"))
sv("figS14_block_bootstrap", s14, 150, 165)

## ========================================================================
## figS06 -- consensus fidelity histograms, dropping the ld_w = 0.2 threshold
## (regenerated from the clustering runs so the four panels fit the page; the
## full five-threshold version remains in /doc)
## ========================================================================
res <- readRDS("data/results_min_loci5.rds")
fl <- rbindlist(lapply(res, function(r) { g <- r$groups[startsWith(r$groups$group_id, "F")]
  g[, .(score, th, n_loci)] }))
fl <- fl[th != 0.2]                                   # drop the 0.2 panel
fl[, th_label := sprintf("ld_w threshold = %s   (n = %s, mean = %.3f)",
     th, format(.N, big.mark = ","), mean(score)), by = th]
fl[, th_label := factor(th_label, levels = unique(th_label[order(th)]))]
fl[, cluster_type := fifelse(n_loci == 1, "singleton (no merge possible)", "merged (n_loci > 1)")]
hb <- fl[, { br <- seq(0.80, 1, by = 0.01); h <- hist(score, breaks = br, plot = FALSE)
  .(bin_mid = h$mids, count = h$counts) }, by = .(th_label, cluster_type)]
hb <- fl[, .(N = .N), by = th_label][hb, on = "th_label"][, pct := 100 * count / N]
s6 <- ggplot(hb, aes(bin_mid, pct, fill = cluster_type)) +
  geom_col(colour = "white", linewidth = 0.1, width = 0.01) +
  geom_vline(xintercept = 0.80, linetype = 2, colour = Pol, linewidth = 0.5) +
  scale_fill_manual(values = c("singleton (no merge possible)" = Grey, "merged (n_loci > 1)" = Accent),
                    name = NULL) +
  facet_wrap(~ th_label, ncol = 1) + coord_cartesian(xlim = c(0.80, 1)) +
  labs(x = expression(score[eMLG] == cor(round(x), x)^2),
       y = "% of that run's flagged clusters") +
  theme_classic(base_size = 9) +
  theme(legend.position = "top", legend.text = element_text(size = 7),
        strip.text = element_text(size = 7.5, face = "bold"), axis.title = element_text(size = 8))
sv("figS06_fidelity_hist", s6, 120, 150)

cat("Saved figS06, figS08, figS11, figS12, figS14 to manuscript/figures/\n")
