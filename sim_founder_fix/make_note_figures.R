## =========================================================================
## sim_founder_fix -- figures for NOTE_parent_ld (current vs empirical vs mosaic LD)
## Needs the check_parent_ld.R outputs:
##   sim_founder_fix/out/check_current.rds     (bed mode, existing diem_outs_demo replicates)
##   sim_founder_fix/out/check_mosaic_l05.rds  (vcf mode, make_mosaic_founders.R LAMBDA 0.5)
##   sim_founder_fix/out/check_mosaic_l1.rds   (vcf mode, LAMBDA 1)
## and, for the hybrid-level preview, module_allele_specific_sorting/Figures/explore_sim_vs_emp_decay_cM.png
## Run from the formica_hybrid repo root.
## =========================================================================
source("sim_founder_fix/parent_ld_lib.R")
FIG <- "sim_founder_fix/note_parent_ld/figures"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

src <- c(current = "current simulations", mosaic_l05 = "mosaic founders, 0.5 switches/cM",
         mosaic_l1 = "mosaic founders, 1 switch/cM")
r <- lapply(names(src), function(k) readRDS(sprintf("sim_founder_fix/out/check_%s.rds", k))); names(r) <- names(src)
d <- rbind(r$current$compare[, .(species, cls, source = "empirical parents", r2 = r2_emp, lo = NA_real_, hi = NA_real_)],
           rbindlist(lapply(names(src), function(k) r[[k]]$compare[, .(species, cls, source = src[[k]], r2 = r2_test, lo, hi)])))
d[, source := factor(source, levels = c("empirical parents", src))]
pal <- c("black", "#d95f02", "#9ecae1", "#3182bd")

## Figure 1: decay with genetic distance
p1 <- ggplot(d[!grepl("sim cluster", cls)], aes(cls, r2, colour = source, group = source)) +
  geom_line() + geom_point(size = 1.4) +
  facet_wrap(~ species, labeller = as_labeller(c(aquilonia = "F. aquilonia parents", polyctena = "F. polyctena parents"))) +
  scale_y_log10() + scale_colour_manual(values = pal, name = NULL) +
  labs(x = "genetic distance between SNPs", y = expression("mean within-species "*r^2*" (log)")) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
ggsave(file.path(FIG, "fig1_decay.pdf"), p1, width = 9, height = 4.6)

## Figure 2: close pairs inside vs between the simulation's LD clusters
d2 <- d[grepl("sim cluster", cls)]
d2[, cls := factor(cls, levels = c("<5kb same sim cluster", "<5kb different sim cluster"),
                   labels = c("same simulation cluster", "different simulation clusters"))]
p2 <- ggplot(d2, aes(cls, r2, fill = source)) +
  geom_col(position = position_dodge(0.85), width = 0.8) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.85), width = 0.25) +
  facet_wrap(~ species, labeller = as_labeller(c(aquilonia = "F. aquilonia parents", polyctena = "F. polyctena parents"))) +
  scale_fill_manual(values = pal, name = NULL) +
  labs(x = "SNP pairs less than 5 kb apart", y = expression("mean within-species "*r^2)) +
  theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2))
ggsave(file.path(FIG, "fig2_close_pairs.pdf"), p2, width = 9, height = 4.4)

## Figure 3: hybrid-level preview from the current simulations (copied)
f3 <- "module_allele_specific_sorting/Figures/explore_sim_vs_emp_decay_cM.png"
if (file.exists(f3)) file.copy(f3, file.path(FIG, "fig3_hybrid_preview.png"), overwrite = TRUE)

## numbers for the note
tab <- dcast(d[, .(species, cls, source, r2)], species + cls ~ source, value.var = "r2")
info <- rbind(r$current$emp_info[, source := "empirical parents"],
              rbindlist(lapply(names(src), function(k) r[[k]]$test_info[, source := src[[k]]])), fill = TRUE)
saveRDS(list(table = tab, info = info, checks = lapply(r, `[[`, "checks")), file.path(FIG, "note_numbers.rds"))
print(tab, digits = 3); print(info, digits = 3)
