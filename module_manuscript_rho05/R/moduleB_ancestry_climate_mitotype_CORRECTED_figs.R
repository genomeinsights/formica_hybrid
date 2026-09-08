suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
FIGDIR <- "module_manuscript_rho05/Figures"
obj <- readRDS("module_manuscript_rho05/data/moduleB_ancestry_climate_mitotype_CORRECTED.rds")
dt17 <- obj$dt17; res_climate <- obj$res_climate; res_mito <- obj$res_mito

## Fig A: ancestry vs climate, 17 units
m <- melt(dt17, id.vars = c("unit_id", "ancestry"), measure.vars = c("PC1", "PC2", "bio6", "bio11"),
          variable.name = "variable", value.name = "value")
m[res_climate, on = "variable", `:=`(r = i.r_17unit, p_om = i.p_omega_null, p_bp = i.p_block_perm, p_bh = i.p_block_perm_BH)]
m[, lab := sprintf("%s: r=%+.2f, p_Om=%.3f, p_blk=%.3f, p_BH=%.3f", variable, r, p_om, p_bp, p_bh)]
pA <- ggplot(m, aes(ancestry, value)) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey70", linewidth = 0.5) +
  geom_point(size = 1.8, colour = "#315B7D") +
  geom_text(aes(label = unit_id), size = 1.9, vjust = -0.7, colour = "grey30") +
  facet_wrap(~ lab, scales = "free_y") +
  labs(x = "genome-wide aquilonia ancestry (17 independent lineage units)", y = "value",
       title = "CORRECTED: ancestry vs climate, 17 lineage units") +
  theme_classic(base_size = 8.5)
ggsave(file.path(FIGDIR, "moduleB_ancestry_vs_climate_17units_CORRECTED.png"), pA, width = 9, height = 5, dpi = 200)

## Fig B: mitotype vs ancestry/climate, 17 units
dt17b <- copy(dt17); dt17b[, Mitotype := ifelse(mito01 == 1, "Faquilonia", "Fpolyctena")]
m2 <- melt(dt17b, id.vars = c("unit_id", "Mitotype"), measure.vars = c("ancestry", "PC1", "PC2", "bio6", "bio11"),
           variable.name = "variable", value.name = "value")
m2[res_mito, on = "variable", `:=`(r = i.r_17unit, p_om = i.p_omega_null, p_bp = i.p_block_perm, p_bh = i.p_block_perm_BH)]
m2[, lab := sprintf("%s: r=%+.2f, p_blk=%.3f, p_BH=%.3f%s", variable, r, p_bp, p_bh,
                    ifelse(is.na(p_om), "", sprintf(", p_Om=%.3f", p_om)))]
pB <- ggplot(m2, aes(Mitotype, value)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, colour = "grey60") +
  geom_jitter(width = 0.08, size = 1.6, colour = "#315B7D") +
  facet_wrap(~ lab, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL, title = "CORRECTED: mitotype vs ancestry/climate, 17 lineage units") +
  theme_classic(base_size = 8)
ggsave(file.path(FIGDIR, "moduleB_mitotype_vs_ancestry_climate_17units_CORRECTED.png"), pB, width = 13, height = 3.6, dpi = 200)

## Fig C: PC1-ancestry, raw vs partial | mitotype, 17 units
res_pc1 <- residuals(lm(PC1 ~ mito01, data = dt17))
res_anc <- residuals(lm(ancestry ~ mito01, data = dt17))
dt17c <- copy(dt17b); dt17c[, `:=`(res_pc1 = res_pc1, res_anc = res_anc)]
p1 <- ggplot(dt17c, aes(PC1, ancestry, colour = Mitotype)) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey60", linewidth = 0.4) +
  geom_point(size = 2) +
  labs(title = sprintf("Raw (17 units): r=%.3f", obj$raw_r17), x = "PC1", y = "ancestry") +
  theme_classic(base_size = 9) + theme(legend.position = "bottom")
p2 <- ggplot(dt17c, aes(res_pc1, res_anc, colour = Mitotype)) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey60", linewidth = 0.4) +
  geom_point(size = 2) +
  labs(title = sprintf("Partial (| mitotype): r=%.3f, p_Om=%.3f, p_blk=%.3f", obj$partial_r17, obj$p_om_partial, obj$p_perm_partial),
       x = "PC1 residual (after mitotype)", y = "ancestry residual (after mitotype)") +
  theme_classic(base_size = 9) + theme(legend.position = "bottom")
p_combo <- p1 + p2 + plot_layout(guides = "collect") +
  plot_annotation(title = "CORRECTED: does PC1's ancestry link survive controlling for mitotype? (17 units)")
ggsave(file.path(FIGDIR, "moduleB_PC1_ancestry_partial_mitotype_17units_CORRECTED.png"), p_combo, width = 9, height = 4.5, dpi = 200)

cat("wrote 3 corrected figures\n")
