## =========================================================
## Module C prerequisite: PC1/PC2 provenance & the PC2-ancestry confound
## =========================================================
## Manuscript-relevant. Establishes that the BayPass covariates PC1/PC2 are
## per-population CLIMATE axes (not a genetic-structure artifact of the same
## genotypes), quantifies how far each is confounded with genome-wide ancestry,
## and documents how that confound is controlled -- so the sorting x BayPass
## cross-check (Module C) is a genuine extrinsic (climate) test, not a circular
## genotype-vs-structure one.
##
## FINDINGS / DECISIONS (see numbers printed below):
##  1. PC1/PC2 are CONSTANT within every population (within-pop SD = 0) => they
##     are per-population covariates, NOT a per-individual genetic PCA of these
##     genotypes. This rules out the worst-case circularity. (Confirmed by PK:
##     they are climate variables; which climate variables load on PC1 vs PC2 is
##     TBD and needed to name the axes.)
##  2. Ancestry link (cor with per-population aquilonia ancestry, ~49k diagnostic
##     markers): PC1 weak overall (-0.23) though Aland-masked (-0.50 without
##     Aland); PC2 moderate (+0.57). A real but partial confound -- sorting IS
##     ancestry fixation, so a PC2 climate outlier could co-localise with sorting
##     via shared ancestry rather than climate. (Both PCs therefore need the
##     structure control, not just PC2.)
##  3. LEVERAGE: the PC2-ancestry link is driven by ALAND (+0.57 -> +0.38 when
##     dropped; extreme low PC2 AND lowest ancestry). SIELVA does NOT drive it --
##     dropping Sielva makes it WORSE (+0.57 -> +0.64; its extreme PC2 sits
##     against the ancestry trend). A residual structural correlation (~0.35)
##     persists without either extreme.
##  4. CONTROL: BayPass's Omega (the among-population covariance) conditions out
##     exactly this shared-ancestry structure, and it is present in BOTH runs.
##     with_aland/run_baypass.sh step 1 estimates Omega on the LD-PRUNED set and
##     the withOmega association reuses it (cleaner, LD-deflated); noOmega lets
##     BayPass estimate its own internally. => PRIMARY = withOmega; noOmega is a
##     robustness contrast. ALAND leverage is handled by the aland_excluded runs.
##     (Caveat: the pre-estimated Omega came from an EARLIER complexity-reduction
##     pruning [pruned_markers.rds], to be regenerated with the final _cM05
##     pruning later; adequate for now -- PK.)
##  5. NO new BayPass run without Sielva: the evidence (point 3) shows it would
##     not help this confound. sielva_excluded/ inputs exist but were left unrun
##     for this reason. Sielva's extreme PC2 does give it leverage on the PC2
##     ASSOCIATION itself, but the >=10-significant-member outlier definition and
##     Omega already guard against single-population artifacts; revisit only if a
##     specific PC2 outlier region looks Sielva-driven.
##
## So Module C uses runs we already have: withOmega as primary (structure-
## controlled), aland_excluded as the leverage robustness, noOmega/with_aland as
## uncontrolled comparators. The PC2-sorting overlap counts as climate-specific
## only if it survives withOmega AND aland_excluded.

suppressMessages({ library(data.table); library(ggplot2) })

e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
sd  <- as.data.table(e$sample_data_with_parents)
GT  <- e$GTs_with_parents
map <- copy(e$map_hyb_005)

aqu <- "aquilonia_parent"; pol <- "polyctena_parent"
h   <- sd[!grepl("_parent$", Population)]
hyb_pops <- unique(h$Population)

## ---- 1. per-population covariate check (constant within population?) ----
const <- h[, .(sd_PC1 = sd(PC1), sd_PC2 = sd(PC2)), by = Population]
cat(sprintf("[1] max within-population SD:  PC1 = %.4f  PC2 = %.4f\n",
            max(const$sd_PC1), max(const$sd_PC2)),
    "    (0 => per-population covariate, NOT a per-individual genetic PCA)\n")

## ---- 2. per-population aquilonia ANCESTRY, oriented by parents ----
## orient each marker so dosage counts the aquilonia allele, then average the
## per-population aquilonia-allele frequency over DIAGNOSTIC markers (parents
## clearly differentiated, |dp| >= 0.5, and DI above the enrichment threshold).
pa_v <- colMeans(GT[sd$Population == aqu, ], na.rm = TRUE) / 2
pp_v <- colMeans(GT[sd$Population == pol, ], na.rm = TRUE) / 2
map[, `:=`(pa = pa_v[marker], pp = pp_v[marker])]
map[, sign_aqu := sign(pa - pp)]
diag_mk <- map[DiagnosticIndex > -25 & abs(pa - pp) >= 0.5 & sign_aqu != 0, marker]
diag_mk <- intersect(diag_mk, colnames(GT))
s_aqu   <- setNames(map$sign_aqu, map$marker)[diag_mk]
cat("[2] diagnostic markers used for the ancestry estimate:", length(diag_mk), "\n")

anc <- vapply(hyb_pops, function(P) {
  f <- colMeans(GT[h[Population == P, Sample_ID], diag_mk, drop = FALSE], na.rm = TRUE) / 2
  mean(ifelse(s_aqu > 0, f, 1 - f), na.rm = TRUE)          # aquilonia-allele freq
}, numeric(1))

dt <- h[, .(PC1 = PC1[1], PC2 = PC2[1]), by = Population]
dt[, ancestry := anc[Population]]

## ---- 3. correlations with ancestry ----
cat(sprintf("\n[3] cor(PC1, ancestry) = %+.3f   cor(PC2, ancestry) = %+.3f\n",
            cor(dt$PC1, dt$ancestry), cor(dt$PC2, dt$ancestry)))

## ---- 4. leverage: is the PC2-ancestry link driven by Aland / Sielva? ----
cat("\n[4] leverage (drop high-leverage populations):\n")
lev <- function(drop, lab) {
  d <- dt[!Population %in% drop]
  cat(sprintf("    %-20s cor(PC1)=%+.3f  cor(PC2)=%+.3f  (n=%d)\n",
              lab, cor(d$PC1, d$ancestry), cor(d$PC2, d$ancestry), nrow(d)))
}
lev(character(0), "all populations")
lev("Aland", "drop Aland")
lev("Sielva", "drop Sielva")
lev(c("Aland", "Sielva"), "drop Aland+Sielva")

## ---- 5. figure: PC vs ancestry, extremes highlighted ----
dir.create("Figures", showWarnings = FALSE)
m <- melt(dt, id.vars = c("Population", "ancestry"),
          measure.vars = c("PC1", "PC2"), variable.name = "PC", value.name = "score")
m[, extreme := Population %in% c("Aland", "Sielva")]
p <- ggplot(m, aes(ancestry, score)) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey70", linewidth = 0.5) +
  geom_point(aes(colour = extreme), size = 1.8) +
  geom_text(data = m[extreme == TRUE], aes(label = Population),
            vjust = -0.7, size = 2.4) +
  scale_colour_manual(values = c("FALSE" = "grey40", "TRUE" = "#D55E00"), guide = "none") +
  facet_wrap(~ PC, scales = "free_y") +
  labs(x = "genome-wide aquilonia ancestry (per population)", y = "climate PC score") +
  theme_classic(base_size = 9)
ggsave("Figures/moduleC_ancestry_confound.png", p, width = 150, height = 70, units = "mm", dpi = 300)

saveRDS(dt, "data/moduleC_ancestry_confound.rds")
cat("\nSaved: data/moduleC_ancestry_confound.rds, Figures/moduleC_ancestry_confound.png\n")
