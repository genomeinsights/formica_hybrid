## =============================================================================
## module_manuscript_rho05 -- high-DI Fst stratified by sorting class and
## prop_fixed, complementing the per-locus Fst violin (di25_fst_vs_di_violin_rho05.R)
## =============================================================================
## USER HYPOTHESIS (2026-09-13): high DI converts genome-wide variation in
## population ancestry into broadly elevated Fst; local sorting affects how
## EXTREME and DIRECTIONAL that differentiation becomes, but it is not its
## sole source. Demonstrated directly here by stratifying per-locus Fst
## (already computed, MAF-folded, in di25_fst_vs_di_violin_rho05.rds) by
## Module A's sort_class and continuous prop_fixed -- computed FRESH on the
## exact same 395,996 MAF-gated units (not reused from moduleA_stage1_cluster_
## sorting.rds, which only covers the 18,361-unit Stage-1-direct subset) via
## the identical parallelism_stats() recipe moduleA_stage1_cluster_sorting.R
## uses, just scaled to the full LD-reduced unit set. All 20 hybrid
## populations used (Aland included) -- this is a pure Fst/sorting question,
## not tied to the BayPass/Omega population universe (Issue 5 in
## AUDIT_FIXES.md), matching di25_fst_vs_di_rho05.R's own population handling.
##
## High-DI cutoff: DI > -25 (the established DI25 convention).
##
## Run from the repo root:
##   Rscript module_manuscript_rho05/R/di25_fst_vs_di_sorting_stratified_rho05.R
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2) })
devtools::load_all("~/gitlab/LDscnR/")
source("moduleA_sorting/R/parallelism_stats.R")   # parallelism_stats(), classify_sort()

VIOLIN_RDS <- "module_manuscript_rho05/data/di25_fst_vs_di_violin_rho05.rds"
OUTRDS     <- "module_manuscript_rho05/data/di25_fst_vs_di_sorting_stratified_rho05.rds"
FIGDIR     <- "module_manuscript_rho05/Figures"
OUTPNG     <- file.path(FIGDIR, "di25_fst_vs_di_sorting_stratified_rho05.png")
OUTPDF     <- sub("\\.png$", ".pdf", OUTPNG)
DI_HIGH_CUTOFF <- -25   # established DI25 convention
MIN_PARENT_MAF <- 0.15; SORT_TH <- 0.6; FIX_TH <- 0.15; SORT_RULE <- "binom"; ALPHA <- 0.05
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

vio <- readRDS(VIOLIN_RDS)
units <- copy(vio$units)   # 395,996 MAF-gated units: group_id, best, DI, bin, pmaf, fst_locus
message(sprintf("[fst-sort] %d MAF-gated units from the violin analysis; %d have DI>%d (high-DI subset)",
                nrow(units), sum(units$DI > DI_HIGH_CUTOFF), DI_HIGH_CUTOFF))

## ---- genotypes at each unit's marker, ALL 20 hybrid pops + parents --------
e2 <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
GTs_wp <- e2$GTs_with_parents; sample_data <- e2$sample_data_with_parents
aqu_pops <- "aquilonia_parent"; pol_pops <- "polyctena_parent"
hybrid_pops <- setdiff(unique(sample_data$Population), c(aqu_pops, pol_pops))   # all 20 -- no Aland exclusion
parent_ids  <- sample_data[grepl("_parent$", Population), Sample_ID]
GTs_parents <- GTs_wp[parent_ids, , drop = FALSE]
hybrid_ids  <- setdiff(rownames(GTs_wp), parent_ids)
GTs_hybrid  <- GTs_wp[hybrid_ids, , drop = FALSE]

stopifnot("some unit markers missing from GTs_with_parents" = all(units$best %in% colnames(GTs_wp)))
hyb_units <- GTs_hybrid[, units$best, drop = FALSE]
par_units <- GTs_parents[, units$best, drop = FALSE]
GTs_units  <- rbind(hyb_units, par_units)
pops_units <- sample_data[match(rownames(GTs_units), Sample_ID), Population]
DI_units   <- setNames(units$DI, units$best)
maf_units  <- setNames(units$pmaf, units$best)   # already correctly FOLDED (see Fst-script fix)

message(sprintf("[fst-sort] classifying %d units (all 20 hybrid pops + parents, %d individuals) ...", nrow(units), nrow(GTs_units)))
t0 <- Sys.time()
prep_units <- ohta_fast_prepare(GTs_units, pops = pops_units)
ps <- parallelism_stats(prep_units, hybrid_pops = hybrid_pops, aqu_pops = aqu_pops,
                        pol_pops = pol_pops, DI = DI_units, min_DI = NULL,
                        parent_maf = maf_units, min_parent_maf = MIN_PARENT_MAF,
                        sort_th = SORT_TH, fix_th = FIX_TH, sort_rule = SORT_RULE, alpha = ALPHA)
message(sprintf("      done | %.0fs", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

setnames(ps, "marker", "best")
uni_cls <- c("aquilonia", "polyctena")
ps[, sort_class := NA_character_]
ok <- ps$differentiated & ps$n_obs > 0 & !is.na(ps$uni_score)
ps[ok, sort_class := classify_sort(n_aqu, n_pol, n_obs, sort_th = SORT_TH, sort_rule = SORT_RULE, alpha = ALPHA)]
## BUGFIX (caught from the first render): classify_sort() also returns
## "unsorted" (differentiated, observed, but below the fixation threshold --
## the LARGEST category, not an edge case) and "ambiguous" (binom rule,
## sorted but too few fixed individuals for power) -- both were falling
## through to raw pass-through before, and "unsorted" was then silently
## dropped to NA because it wasn't in lab_order below.
ps[, sort_label := fifelse(!differentiated, "not differentiated",
                    fifelse(is.na(sort_class), "differentiated, unclassified",
                    fifelse(sort_class %in% uni_cls, paste0("sorted: ", sort_class),
                    fifelse(sort_class == "unresolved", "differentiated, unresolved",
                    fifelse(sort_class == "unsorted", "differentiated, not fixed",
                    fifelse(sort_class == "ambiguous", "differentiated, low power", sort_class))))))]

units <- ps[, .(best, differentiated, prop_fixed, uni_score, sort_class, sort_label)][units, on = "best"]
stopifnot("row count changed after sorting-classification join" = nrow(units) == nrow(vio$units))
message("[fst-sort] sort_label counts (all MAF-gated units):")
print(units[, .N, by = sort_label][order(-N)])

saveRDS(list(units = units, di_high_cutoff = DI_HIGH_CUTOFF,
             sort_params = list(min_parent_maf = MIN_PARENT_MAF, sort_th = SORT_TH,
                                fix_th = FIX_TH, sort_rule = SORT_RULE, alpha = ALPHA)), OUTRDS)

## ---- high-DI subset for the figure -----------------------------------------
hi <- units[DI > DI_HIGH_CUTOFF & is.finite(fst_locus)]
message(sprintf("[fst-sort] %d units with DI>%d used for the figure", nrow(hi), DI_HIGH_CUTOFF))
lab_order <- c("not differentiated", "differentiated, unclassified", "differentiated, not fixed",
               "differentiated, low power", "differentiated, unresolved",
               "sorted: aquilonia", "sorted: polyctena")
stopifnot("a sort_label value is missing from lab_order" = all(unique(hi$sort_label) %in% lab_order))
hi[, sort_label := factor(sort_label, levels = intersect(lab_order, unique(sort_label)))]

n_lab <- hi[, .N, by = sort_label]
lab_with_n <- setNames(sprintf("%s\n(n=%s)", n_lab$sort_label, format(n_lab$N, big.mark = ",")), n_lab$sort_label)

## ---- panel A: Fst by sort_class (violin+box); panel B: Fst vs prop_fixed ---
pA <- ggplot(hi, aes(sort_label, fst_locus, fill = sort_label)) +
  geom_violin(alpha = 0.3, scale = "width", trim = TRUE, linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, outlier.alpha = 0.3, fill = "white", linewidth = 0.4) +
  scale_x_discrete(labels = lab_with_n) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = NULL, y = expression("per-locus " * F[ST]),
       title = "a  Fst by sorting class", subtitle = sprintf("DI>%d, MAF>=%.2f", DI_HIGH_CUTOFF, MIN_PARENT_MAF)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1))

pB <- ggplot(hi, aes(prop_fixed, fst_locus)) +
  geom_point(size = 0.4, alpha = 0.15, colour = "#377EB8") +
  geom_smooth(method = "loess", colour = "#D95F02", linewidth = 0.9, se = TRUE) +
  labs(x = "prop_fixed (sorting magnitude)", y = expression("per-locus " * F[ST]),
       title = "b  Fst vs sorting magnitude", subtitle = sprintf("DI>%d, MAF>=%.2f, differentiated units only", DI_HIGH_CUTOFF, MIN_PARENT_MAF)) +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())

library(patchwork)
full <- pA + pB + plot_layout(widths = c(1.3, 1))
ggsave(OUTPDF, full, width = 12, height = 5.2)
ggsave(OUTPNG, full, width = 12, height = 5.2, dpi = 200)
cat("saved:", OUTPNG, "\n")
