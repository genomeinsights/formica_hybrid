## =========================================================
## MODULE C (C3) -- sorting x climate: do sorted clusters co-localise with
## BayPass PC1/PC2 climate-association outliers?
## =========================================================
## The extrinsic-side test of the headline question. Companion scripts:
##   C1  R/analyse_baypass.R              -- outlier definition + Manhattans
##   C2  R/diagnostic_index_enrichment.R  -- DI-enrichment of outlier clusters
##   this file (C3)                       -- sorting x outlier overlap
##   C4  (later)                          -- PC1/PC2 axis identity + reconciliation
##                                           (needs the climate-variable loadings)
##
## What is new here vs the old sorting_baypass_crosscheck{,_memberlevel}.R:
##  - sort_class is the eMLG CONSENSUS call for EVERY has_eMLG cluster (stored
##    hybrid-side eMLG + a matched parent-side consensus), not the representative
##    tag SNP -- per the eMLG-vs-representative validation (representatives agree
##    on direction but mis-call magnitude for >=5-loci blocks; the overlap is a
##    magnitude question, so use the consensus).
##  - the sorting gate is the pipeline-locked parent_maf >= 0.15 (was the old
##    min_parent_diff=0.5 / fix_th=0.1); null_prob 0.5, sort_th 0.5, di_agg "max".
##  - the ancestry confound (PC2 ~ ancestry +0.57; see moduleC_ancestry_confound.R)
##    is handled by BayPass's Omega, present in BOTH runs. PRIMARY = withOmega
##    (Omega pre-estimated on the LD-pruned set) x aland_excluded (removes the
##    Aland leverage point). The overlap counts as climate-specific only if it
##    survives that; noOmega / with_aland are uncontrolled comparators. All 8
##    configs are reported as a robustness grid.
##
## DESCRIPTIVE co-localisation (a climate lead), not a causal claim.
## Run from repo root. Reads BayPass .out (all present), canonical clustering,
## GTs, moduleA_snp.rds (member-level), C2's enrichment CSV (Fig 3a).

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("dev/R/Ohta.R"); source("dev/R/parallelism_stats.R"); source("dev/R/eMLG_parallelism.R")
set.seed(1)

## ---- PARAMETERS (locked, match Module A) --------------------------------
MIN_PARENT_MAF <- 0.15; FIX_TH <- 0.15; SORT_TH <- 0.5; NULL_PROB <- 0.5
DI_AGG <- "max"; CORES <- 8
## BayPass outlier definition: >= MIN_N_SIG members with BF(dB) >= SIG_THR.
## RELAXED from the original 20/10 to 15/5 to raise the number of PUTATIVE outlier
## clusters (aland_excluded PC2 had only 2-4, leaving the sorting overlap
## untestable). Trade-off: false positives under a relaxed rule are SIZE-BIASED --
## the chance of >=5 members crossing grows with cluster size and is amplified by
## within-cluster LD -- and cluster size independently predicts DI, so they can ADD
## spurious enrichment rather than merely dilute it. Read the SIZE-ADJUSTED
## estimates, not the naive ones. MUST match R/diagnostic_index_enrichment.R (C2),
## or Fig 3a and 3b describe different outlier sets.
SIG_THR <- 15; MIN_N_SIG <- 5
TAG <- sprintf("%d_%d", MIN_N_SIG, SIG_THR)   # output suffix so both settings coexist
CLUSTERING <- "data/eMLG_5loci_0025_cM05.rds"
uni <- c("aquilonia", "polyctena")
elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

## ---- inputs -------------------------------------------------------------
e2 <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
GTs_wp <- e2$GTs_with_parents; sample_data <- e2$sample_data_with_parents; map <- e2$map_hyb_005
e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1); GTs_hyb <- e1$GTs_hybrids_005

clust  <- readRDS(CLUSTERING); groups <- clust$groups; eMLG <- clust$eMLG
DI_vec <- setNames(map$DiagnosticIndex, map$marker)
aqu <- "aquilonia_parent"; pol <- "polyctena_parent"
hyb <- setdiff(unique(sample_data$Population), c(aqu, pol))
parent_ids  <- sample_data[grepl("_parent$", Population), Sample_ID]
GTs_parents <- GTs_wp[parent_ids, , drop = FALSE]
message("[C3] inputs loaded (", nrow(GTs_wp), " indiv x ", ncol(GTs_wp), " markers)")

## ========================================================================
## C3.1 -- BayPass outlier clusters (has_eMLG cluster with >=10 BF>=20 members)
## ========================================================================
he  <- groups[has_eMLG == TRUE]
m2g <- he[, .(marker = unlist(members)), by = group_id]     # member -> cluster
import_bf <- function(f) {
  r <- fread(f)
  stopifnot(nrow(r) == nrow(map), identical(r$MRK, seq_len(nrow(r))))
  setNames(r$`BF(dB)`, map$marker)
}
outlier_of <- function(bf) {
  m2g[, BF := bf[marker]]                       # ONE vectorised name lookup (NOT per group)
  m2g[, .(n_sig = sum(BF >= SIG_THR, na.rm = TRUE)), by = group_id][n_sig >= MIN_N_SIG, group_id]
}
cfg <- CJ(popset = c("with_aland", "aland_excluded"), pc = c("PC1", "PC2"),
          omega = c("noOmega", "withOmega"))
oid <- lapply(seq_len(nrow(cfg)), function(i) {
  c <- cfg[i]; outlier_of(import_bf(sprintf("./%s/%s_DIEM_%s_summary_betai_reg.out", c$popset, c$pc, c$omega)))
})
names(oid) <- cfg[, paste(popset, pc, omega, sep = "_")]
message("[C3.1] outlier clusters computed")
cat("[C3.1] outlier clusters per config:\n"); print(sapply(oid, length))
## primary structure-controlled sets, and per-PC unions (for robustness)
prim_PC1 <- oid[["aland_excluded_PC1_withOmega"]]
prim_PC2 <- oid[["aland_excluded_PC2_withOmega"]]

## ========================================================================
## C3.2 -- consensus sort_class for EVERY has_eMLG cluster (CHECKPOINTED:
##         the parent build is the only slow step, so cache its result and skip
##         it on re-runs -- a downstream bug must never discard it again).
## ========================================================================
## NB deliberately NOT tagged by TAG: the consensus sort_class does not depend on
## the BayPass outlier definition, so it is reused across threshold settings --
## which is exactly what makes re-running at new thresholds cheap.
CKPT <- "data/moduleC_C3_cl.rds"
if (file.exists(CKPT)) {
  message("[C3.2] loading cached consensus sort_class (skipping the parent build): ", CKPT)
  cl <- readRDS(CKPT)
} else {
  has_ids <- colnames(eMLG)                                  # >=5-loci clusters
  umem <- setNames(groups[.(has_ids), on = "group_id", members], has_ids)
  message("[C3.2] building matched parent-side consensus for ", length(has_ids), " has_eMLG clusters ...")
  t0 <- Sys.time()
  par_cons <- build_group_consensus(umem, GTs_hyb, GTs_parents, cores = CORES, progress = TRUE, label = "C3 parent")
  message(sprintf("      parent consensus done | %.0fs", elapsed(t0)))
  GTs_units  <- rbind(eMLG[, has_ids, drop = FALSE], par_cons)
  pops_units <- sample_data[match(rownames(GTs_units), Sample_ID), Population]
  DI_units   <- cluster_DI(groups, has_ids, DI_vec, di_agg = DI_AGG)
  maf_units  <- { pf <- colMeans(par_cons, na.rm = TRUE) / 2; pmin(pf, 1 - pf) }
  prep_units <- ohta_fast_prepare(GTs_units, pops = pops_units)
  ps <- parallelism_stats(prep_units, hybrid_pops = hyb, aqu_pops = aqu, pol_pops = pol,
                          DI = DI_units, min_DI = NULL, parent_maf = maf_units, min_parent_maf = MIN_PARENT_MAF,
                          sort_th = SORT_TH, fix_th = FIX_TH, null_prob = NULL_PROB)
  setnames(ps, "marker", "group_id")
  cl <- groups[.(has_ids), on = "group_id", .(group_id, n_loci)][
    ps[, .(group_id, differentiated, sort_class, DI, prop_fixed, uni_score)], on = "group_id"]
  cl[, `:=`(directional = as.integer(sort_class %in% uni),
            sorted      = as.integer(sort_class %in% c(uni, "bidirectional")))]
  saveRDS(cl, CKPT)
  message("[C3.2] cached consensus sort_class -> ", CKPT)
}

## ========================================================================
## C3.3 -- overlap: are outlier clusters more likely to be sorted? (size-adj.)
## ========================================================================
overlap_or <- function(ids, resp = "directional") {
  d <- cl[differentiated == TRUE]; d[, is_out := as.integer(group_id %in% ids)]
  if (sum(d$is_out) < 3) return(data.table(OR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, n_out = sum(d$is_out)))
  m <- glm(as.formula(paste(resp, "~ is_out + log(n_loci)")), data = d, family = binomial)
  ci <- confint.default(m, "is_out"); cf <- summary(m)$coefficients
  data.table(OR = exp(cf["is_out","Estimate"]), lo = exp(ci[1]), hi = exp(ci[2]),
             p = cf["is_out","Pr(>|z|)"], n_out = sum(d$is_out))
}
cat("\n[C3.3] size-adjusted directional-overlap OR across all 8 configs (robustness grid):\n")
grid <- cbind(cfg, rbindlist(lapply(names(oid), function(nm) overlap_or(oid[[nm]]))))
grid[, `:=`(OR = round(OR, 2), lo = round(lo, 2), hi = round(hi, 2), p = signif(p, 2))]
print(grid[order(pc, popset, omega)])
cat("\nPRIMARY (aland_excluded x withOmega), directional overlap:\n")
cat("  PC1:"); print(overlap_or(prim_PC1)); cat("  PC2:"); print(overlap_or(prim_PC2))

## ---- member-level robustness (fraction of differentiated members directional) ----
snp <- readRDS("data/moduleA_snp.rds")
mg  <- snp[, .(marker, m_diff = differentiated, m_dir = sort_class %in% uni)][m2g, on = "marker"]
memb <- mg[, .(n_loci = .N, n_diff = sum(m_diff, na.rm = TRUE),
               frac_dir = mean(m_dir[m_diff], na.rm = TRUE)), by = group_id][n_diff >= 5]
memb_or <- function(ids) {
  d <- copy(memb); d[, is_out := as.integer(group_id %in% ids)]
  m <- lm(frac_dir ~ is_out + log(n_loci), data = d); cf <- summary(m)$coefficients
  data.table(beta = round(cf["is_out","Estimate"], 3), p = signif(cf["is_out","Pr(>|t|)"], 2), n_out = sum(d$is_out))
}
cat("\n[C3.3b] member-level frac_directional ~ outlier + log(n_loci) (primary configs):\n")
cat("  PC1:"); print(memb_or(prim_PC1)); cat("  PC2:"); print(memb_or(prim_PC2))

## ========================================================================
## C3.4 -- direction (aqu/pol) and DI of sorted outlier clusters
## ========================================================================
cat("\n[C3.4] sort_class of PRIMARY outlier clusters (differentiated):\n")
cat("  PC1:\n"); print(cl[differentiated == TRUE & group_id %in% prim_PC1, .N, by = sort_class][order(-N)])
cat("  PC2:\n"); print(cl[differentiated == TRUE & group_id %in% prim_PC2, .N, by = sort_class][order(-N)])
cat("\nmedian DI of sorted vs unsorted PC2-outlier clusters:\n")
print(cl[differentiated == TRUE & group_id %in% prim_PC2,
         .(n = .N, median_DI = round(median(DI, na.rm = TRUE), 1)),
         by = .(sorted = ifelse(directional == 1, "directional", "not"))])

## outlier-vs-background cluster size: did relaxing the definition narrow the gap
## the size adjustment has to extrapolate across? (was outlier ~147 vs bg ~15)
cat("\n[C3.3c] outlier vs background cluster size (median n_loci):\n")
szg <- function(ids, lab) {
  d <- cl[differentiated == TRUE]
  data.table(set = lab, n_out = sum(d$group_id %in% ids),
             med_nloci_outlier = as.double(median(d[group_id %in% ids, n_loci])),
             med_nloci_background = as.double(median(d[!group_id %in% ids, n_loci])))
}
print(rbind(szg(prim_PC1, "PC1 primary"), szg(prim_PC2, "PC2 primary")))

saveRDS(list(cl = cl, grid = grid, memb = memb,
             outliers = list(PC1 = prim_PC1, PC2 = prim_PC2)),
        sprintf("data/moduleC_sorting_climate_%s.rds", TAG))

## ========================================================================
## Figure 3 (a-c)
## ========================================================================
dir.create("Figures", showWarnings = FALSE)
th <- theme_classic(base_size = 8) +
  theme(plot.tag = element_text(face = "bold", size = 10),
        plot.title = element_text(size = 8.5, hjust = 0),
        axis.title = element_text(size = 8), axis.text = element_text(size = 6.5),
        strip.text = element_text(size = 7.5),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 7), legend.key.size = unit(3, "mm"),
        plot.margin = margin(4, 9, 2, 4))
pc_col <- c(PC1 = "#0072B2", PC2 = "#D55E00")
om_lab <- function(o) ifelse(o == "withOmega", "wOm", "noOm")

## (a) DI-enrichment (from C2): size-adjusted OR, PC1 vs PC2 (aland_excluded)
## match C2's parameterised name; fall back to the legacy unsuffixed file (= 20/10)
enr_f <- sprintf("data/diagnostic_index_enrichment%s.csv", TAG)
if (!file.exists(enr_f)) enr_f <- "data/diagnostic_index_enrichment.csv"
p3a <- if (file.exists(enr_f)) {
  enr <- fread(enr_f)[population_set == "aland_excluded"][, ylab := paste(pc, om_lab(omega))]
  ggplot(enr, aes(adj_OR, ylab)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
    geom_errorbarh(aes(xmin = adj_OR_lo, xmax = adj_OR_hi), height = 0.2, colour = "grey30") +
    geom_point(size = 1.8, colour = "grey20") +
    labs(x = "DI-enrichment OR", y = NULL, title = "DI-enrichment") + th
} else patchwork::plot_spacer()

## (b) sorting-overlap OR, PC1 vs PC2 (all configs; primary = larger point)
ov <- copy(grid)[, `:=`(ylab = paste0(ifelse(popset == "with_aland", "Al-in ", "Al-out "), om_lab(omega)),
                        prim = popset == "aland_excluded" & omega == "withOmega")]
p3b <- ggplot(ov, aes(OR, ylab)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, na.rm = TRUE, colour = "grey30") +
  geom_point(aes(size = prim), na.rm = TRUE, colour = "grey20") +
  scale_size_manual(values = c(`FALSE` = 1.4, `TRUE` = 2.7), guide = "none") +
  scale_x_log10(breaks = c(0.3, 3, 30)) + facet_wrap(~ pc) +
  labs(x = "sorting-overlap OR (log)", y = NULL, title = "Sorting overlap") + th

## (c) direction of sorted primary-outlier clusters
dir_dt <- rbind(
  cl[differentiated == TRUE & group_id %in% prim_PC1, .(pc = "PC1", sort_class)],
  cl[differentiated == TRUE & group_id %in% prim_PC2, .(pc = "PC2", sort_class)]
)[!is.na(sort_class)][, .N, by = .(pc, sort_class)]
dir_dt[, sort_class := factor(sort_class, levels = c("unsorted", "aquilonia", "polyctena", "bidirectional"))]
scol <- c(aquilonia = "#0072B2", polyctena = "#D55E00", unsorted = "#999999", bidirectional = "#CC79A7")
p3c <- ggplot(dir_dt, aes(pc, N, fill = sort_class)) +
  geom_col(position = "fill", width = 0.7) + scale_fill_manual(values = scol) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "proportion", title = "Direction") + th

fig3 <- p3a + p3b + p3c + plot_layout(widths = c(1, 1.5, 0.85), guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
ggsave(sprintf("Figures/moduleC_fig3_%s.pdf", TAG), fig3, width = 210, height = 82, units = "mm")
ggsave(sprintf("Figures/moduleC_fig3_%s.png", TAG), fig3, width = 210, height = 82, units = "mm", dpi = 300)
cat(sprintf("\nSaved: data/moduleC_sorting_climate_%s.rds, Figures/moduleC_fig3_%s.pdf/.png\n", TAG, TAG))
