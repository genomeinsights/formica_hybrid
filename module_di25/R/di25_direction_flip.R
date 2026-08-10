## =========================================================
## module_di25 (high-DI analyses) -- the SNP -> eMLG ancestry-sorting DIRECTION FLIP
## =========================================================
## Quantifies the flip: per-SNP the most strongly sorted high-DI loci lean
## F. polyctena (~34% aquilonia-directed at tau = 0.8); collapsed to independent
## LD units the bias inverts to a stable ~80% F. aquilonia. Additive: reads the
## outputs of di25_sorting.R (which it does NOT modify) and computes one new
## series -- an "LD-thinned" control (one representative SNP per cluster, NO
## consensus) that separates two explanations for the flip:
##   * thinned tracks eMLG  -> the flip is INDEPENDENCE (pseudoreplication removed);
##   * thinned tracks perSNP -> the flip is the eMLG CONSENSUS construction.
##
## Run from the repo root:  Rscript module_di25/R/di25_direction_flip.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
devtools::load_all("~/gitlab/LDscnR/")             # ohta_fast_prepare()
source("moduleA_sorting/R/parallelism_stats.R")    # parallelism_stats(), classify_sort()

## ---- parameters (identical to di25_sorting.R) ---------------------------
TAU_GRID    <- c(0.5, 0.6, 0.7, 0.8)
FIX_TH      <- 0.15      # phi = 0.85
MIN_PMAF    <- 0.15
SORT_RULE   <- "binom"
ALPHA       <- 0.05
CM_STAMP    <- "cM5"
TAU_PANEL_C <- 0.6       # tau for the mechanism panel
OUTDIR <- "module_di25/data"; FIGDIR <- "module_di25/Figures"

COL_AQU <- "#21918C"; COL_POL <- "#D3C93B"; COL_HET <- "#440154"; COL_AMB <- "#BDBDBD"
SER_COL <- c("per-SNP" = "grey55", "SNP, LD-thinned" = "#E69F00",
             "per-eMLG" = "#315B7D", "per-eMLG (consensus only)" = "#C2549D")
SER_LTY <- c("per-SNP" = 1, "SNP, LD-thinned" = 1, "per-eMLG" = 1, "per-eMLG (consensus only)" = 2)
SER_SHP <- c("per-SNP" = 16, "SNP, LD-thinned" = 17, "per-eMLG" = 16, "per-eMLG (consensus only)" = 15)

## ---- inputs -------------------------------------------------------------
sweep   <- as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_sweep.rds")))
ps_snp  <- as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_snp.rds")))
ps_emlg <- as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_emlg.rds")))
res <- readRDS(file.path(OUTDIR, sprintf("di25_clustering_%s.rds", CM_STAMP))); g <- res$groups
inp <- readRDS(file.path(OUTDIR, "di25_inputs.rds"))
e2  <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd  <- e2$sample_data_with_parents
DI_vec <- setNames(e2$map_hyb_005$DiagnosticIndex, e2$map_hyb_005$marker)

## ---- rebuild the 194-individual matrix EXACTLY as di25_sorting.R --------
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)
GTs_all <- GTs_all[rownames(GTs_all) %in% sd$Sample_ID, ]
pops     <- sd$Population[match(rownames(GTs_all), sd$Sample_ID)]
aqu_pops <- "aquilonia_parent"; pol_pops <- "polyctena_parent"
hybrid_pops <- setdiff(unique(pops), c(aqu_pops, pol_pops))
parent_rows <- grepl("_parent$", pops)
stopifnot(nrow(GTs_all) == 194L, sum(!parent_rows) == 164L,
          length(hybrid_pops) == 20L, sum(parent_rows) == 30L)
par_freq_snp <- colMeans(GTs_all[parent_rows, , drop = FALSE], na.rm = TRUE) / 2
pmaf_snp     <- pmin(par_freq_snp, 1 - par_freq_snp)

## ---- local tally helper (must reproduce di25_sorting.R::tally_level) ----
tally_level <- function(ps, level) {
  base <- ps[differentiated == TRUE & n_obs > 0]
  rbindlist(lapply(TAU_GRID, function(tau) {
    cls <- classify_sort(base$n_aqu, base$n_pol, base$n_obs,
                         sort_th = tau, sort_rule = SORT_RULE, alpha = ALPHA)
    n_aqu <- sum(cls == "aquilonia"); n_pol <- sum(cls == "polyctena")
    n_unres <- sum(cls == "unresolved"); n_amb <- sum(cls == "ambiguous")
    n_sorted <- n_aqu + n_pol + n_unres + n_amb
    data.table(level = level, tau = tau, n_differentiated = nrow(base),
               n_sorted = n_sorted, pct_sorted = 100 * n_sorted / nrow(base),
               toward_aqu = n_aqu, toward_pol = n_pol,
               dir_unresolved = n_unres, ambiguous = n_amb,
               pct_aqu_of_resolved = 100 * n_aqu / (n_aqu + n_pol))
  }))
}
## correctness check: reproduce the stored sweep for SNP and eMLG
chk <- rbind(tally_level(ps_snp, "SNP"), tally_level(ps_emlg, "eMLG"))
setkey(chk, level, tau); setkey(sweep, level, tau)
repro_ok <- isTRUE(all.equal(chk, sweep[, names(chk), with = FALSE], check.attributes = FALSE))
cat("[check] local tally reproduces di25_sorting_sweep.rds (SNP+eMLG): ", repro_ok, "\n")

## =========================================================================
## LD-thinned control: one representative SNP per cluster, NO consensus
## =========================================================================
reps <- g$representative
stopifnot(!anyDuplicated(reps), all(reps %in% colnames(GTs_all)))
GTs_thin <- GTs_all[, reps, drop = FALSE]                 # 194 x ~11,052 real SNPs
cat(sprintf("[thinned] %d individuals x %d representative SNPs (one per cluster)\n",
            nrow(GTs_thin), ncol(GTs_thin)))
prep_thin <- ohta_fast_prepare(GTs_thin, pops = pops)
ps_thin <- parallelism_stats(prep_thin, hybrid_pops = hybrid_pops,
                             aqu_pops = aqu_pops, pol_pops = pol_pops,
                             fix_th = FIX_TH, DI = DI_vec[reps], min_DI = NULL,   # per-SNP DI, not max-over-members
                             parent_maf = pmaf_snp[reps], min_parent_maf = MIN_PMAF,
                             sort_rule = SORT_RULE, alpha = ALPHA)
thin_sweep <- tally_level(ps_thin, "SNP, LD-thinned")

## supplementary: eMLG restricted to the 4,010 true consensus units
emlg_cons_sweep <- tally_level(ps_emlg[is_emlg == TRUE], "eMLG (consensus only)")

flip <- rbind(copy(sweep)[, level := fifelse(level == "SNP", "per-SNP", "per-eMLG")],
              thin_sweep,
              copy(emlg_cons_sweep)[, level := "per-eMLG (consensus only)"])
LVLS <- c("per-SNP", "SNP, LD-thinned", "per-eMLG", "per-eMLG (consensus only)")
flip[, level := factor(level, levels = LVLS)]
saveRDS(list(flip = flip, repro_ok = repro_ok),
        file.path(OUTDIR, "di25_direction_flip.rds"))

## Wilson CI for a proportion (shown on every series except the pseudoreplicated per-SNP)
wilson <- function(k, n) { z <- 1.959964; p <- k/n; d <- 1 + z^2/n
  ctr <- (p + z^2/(2*n))/d; hw <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/d
  list(lo = 100*pmax(0, ctr-hw), hi = 100*pmin(1, ctr+hw)) }
flip[, nres := toward_aqu + toward_pol]
flip[level != "per-SNP", `:=`(ci_lo = wilson(toward_aqu, nres)$lo, ci_hi = wilson(toward_aqu, nres)$hi)]

## =========================================================================
## Panel a -- pct aquilonia of resolved vs tau (three series)
## =========================================================================
XLIM <- c(0.46, 0.84)
pA <- ggplot(flip, aes(tau, pct_aqu_of_resolved, colour = level, linetype = level, shape = level)) +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey70") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.012, na.rm = TRUE, linetype = 1) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
  ## call out the coincidence directly, with an arrow to the overlapping lines
  annotate("segment", x = 0.585, y = 88, xend = 0.60, yend = 82.5,
           arrow = arrow(length = unit(0.16, "cm")), colour = "grey35", linewidth = 0.4) +
  annotate("label", x = 0.5, y = 90.5, hjust = 0, size = 3.0, fontface = "bold",
           colour = "grey15", label.size = 0, fill = "white",
           label = "thinned tracks eMLG (lines coincide):\nthe flip is independence, not consensus") +
  annotate("text", x = 0.46, y = 27, hjust = 0, size = 2.7, colour = "grey45",
           label = "no CI on per-SNP: markers are not independent (per-SNP n = 2,268-9,017)") +
  scale_colour_manual(values = SER_COL, name = NULL) +
  scale_linetype_manual(values = SER_LTY, name = NULL) +
  scale_shape_manual(values = SER_SHP, name = NULL) +
  scale_x_continuous(breaks = TAU_GRID, limits = XLIM) +
  labs(x = NULL, y = "% of resolved loci toward F. aquilonia",
       title = "a  The direction flip",
       subtitle = "above 50% = aquilonia-biased, below = polyctena-biased") +
  coord_cartesian(ylim = c(25, 94)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

## n as a labelled table row beneath panel a (the small-denominator unit series)
ntab <- flip[level %in% c("per-eMLG", "per-eMLG (consensus only)"), .(level, tau, nres)]
ntab[, row := fifelse(level == "per-eMLG", 2L, 1L)]
pN <- ggplot(ntab, aes(tau, row, label = nres, colour = level)) +
  geom_text(size = 3.1, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = SER_COL) +
  scale_x_continuous(breaks = TAU_GRID, limits = XLIM,
                     labels = function(b) sprintf("%.1f", b)) +
  scale_y_continuous(breaks = c(1, 2), limits = c(0.5, 2.5),
                     labels = c("consensus units  n", "eMLG units  n")) +
  labs(x = expression(sorting~threshold~tau), y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 8, hjust = 1))

## composition strip: sorted-set make-up per level x tau
comp <- melt(flip[, .(level, tau, toward_aqu, toward_pol, dir_unresolved, ambiguous)],
             id.vars = c("level","tau"), variable.name = "cls", value.name = "n")
comp[, cls := factor(cls, levels = c("toward_aqu","dir_unresolved","ambiguous","toward_pol"),
                     labels = c("aquilonia","unresolved","ambiguous","polyctena"))]
pStrip <- ggplot(comp, aes(factor(tau), n, fill = cls)) +
  geom_col(position = "fill", width = 0.8) +
  facet_wrap(~ level, nrow = 1) +
  scale_fill_manual(values = c(aquilonia = COL_AQU, unresolved = COL_HET,
                               ambiguous = COL_AMB, polyctena = COL_POL), name = NULL) +
  labs(x = expression(tau), y = "sorted make-up") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom",
        strip.text = element_text(size = 8), panel.grid = element_blank())

## =========================================================================
## Panel b (mechanism) -- polyctena signal is a few big blocks
## =========================================================================
mc <- data.table(group_id = rep(g$group_id, lengths(g$members)), marker = unlist(g$members))
base_s <- ps_snp[differentiated == TRUE & n_obs > 0]
base_s[, cls := classify_sort(n_aqu, n_pol, n_obs, sort_th = TAU_PANEL_C,
                              sort_rule = SORT_RULE, alpha = ALPHA)]
dir_snp <- mc[base_s[cls %in% c("aquilonia","polyctena"), .(marker, cls)], on = "marker"]
percluster <- dir_snp[, .(n_snp = .N), by = .(cls, group_id)][order(cls, -n_snp)]
cumdt <- percluster[, .(group_id, rank = seq_len(.N),
                        cum_share = cumsum(n_snp) / sum(n_snp)), by = cls]
top10 <- cumdt[rank == 10, .(cls, top10_share = round(cum_share, 3))]
cat("[panel b] share of each direction's sorted SNPs from the top-10 clusters:\n"); print(top10)

## identity of the top polyctena / aquilonia clusters
memb_pos <- data.table(group_id = rep(g$group_id, lengths(g$members)),
                       marker = unlist(g$members))
memb_pos[, `:=`(chr = sub(":.*","",marker), pos = as.integer(sub(".*:","",marker)))]
top_clusters <- percluster[, head(.SD, 5), by = cls][
  g[, .(group_id, n_loci, Chr)], on = "group_id", nomatch = 0][order(cls, -n_snp)]
top_clusters <- memb_pos[, .(pos_lo = min(pos), pos_hi = max(pos)), by = group_id][
  top_clusters, on = "group_id"]
cat("[panel b] top clusters per direction (expect Chr5/25/26 for polyctena):\n")
print(top_clusters[, .(cls, group_id, Chr, n_loci, n_snp, pos_lo, pos_hi)])

pB <- ggplot(cumdt, aes(rank, cum_share, colour = cls)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c(aquilonia = COL_AQU, polyctena = COL_POL), name = NULL) +
  scale_x_log10() + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  annotate("text", x = 1.15, y = 0.93, hjust = 0, vjust = 1, size = 3.0, colour = "grey25",
    label = sprintf("top-10 clusters supply\n%.0f%% of polyctena- vs %.0f%% of\naquilonia-directed sorted SNPs\n(polyctena blocks: Chr 5, 25, 26)",
                    100*top10[cls=="polyctena", top10_share], 100*top10[cls=="aquilonia", top10_share])) +
  labs(x = "cluster rank (log; clusters ordered by sorted-SNP count)",
       y = "cumulative share of direction's sorted SNPs",
       title = sprintf("b  Why: polyctena signal is a few blocks (tau = %.1f)", TAU_PANEL_C)) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

## =========================================================================
## assemble + save
## =========================================================================
## Main talk figure: panel a (with an n table row) beside the mechanism panel b.
## The composition strip is deliberately NOT here -- it is backup-slide material.
left <- pA / pN + plot_layout(heights = c(6.5, 1))
fig  <- (left | pB) + plot_layout(widths = c(1.25, 1)) +
  plot_annotation(title = "Ancestry sorting flips direction from SNPs to independent LD units",
                  subtitle = "SNP and unit sorted-sets are not nested subsets -- a shift in composition, not a reclassification of the same loci",
                  theme = theme(plot.title = element_text(face = "bold", size = 12)))
ggsave(file.path(FIGDIR, "di25_direction_flip.png"), fig, width = 13, height = 6, dpi = 200)

## Backup slide: the sorted-set composition strip on its own.
ggsave(file.path(FIGDIR, "di25_direction_flip_composition.png"),
       pStrip + labs(title = "Backup: sorted-set composition per level x tau",
                     subtitle = "pct toward aquilonia (panel a) excludes the unresolved/ambiguous slices shown here"),
       width = 11, height = 3.4, dpi = 200)
cat("[direction-flip] wrote main figure + composition-strip backup\n")

## =========================================================================
## acceptance checks + interpretation
## =========================================================================
cat("\n===== ACCEPTANCE CHECKS =====\n")
cat("1. tally reproduces stored sweep (SNP+eMLG): ", repro_ok, "\n")
cat(sprintf("2. sorted%% at tau=0.6: per-SNP %.1f (exp ~12.6) | per-eMLG %.1f (exp ~7.2)\n",
            sweep[level=="SNP" & tau==0.6, pct_sorted], sweep[level=="eMLG" & tau==0.6, pct_sorted]))
cat(sprintf("3. pct_aqu at tau=0.8: per-SNP %.1f (exp ~34) | per-eMLG %.1f (exp ~80)\n",
            sweep[level=="SNP" & tau==0.8, pct_aqu_of_resolved], sweep[level=="eMLG" & tau==0.8, pct_aqu_of_resolved]))
cat("4. individuals 194 / 164 hybrids / 20 pops / 30 parents: OK (asserted)\n")
cat(sprintf("5. thinned matrix %d cols; all colnames are per-SNP markers: %s\n",
            ncol(GTs_thin), all(colnames(GTs_thin) %in% ps_snp$marker)))

cat("\n----- THINNED CONTROL (the interpretation decider) -----\n")
print(dcast(flip, tau ~ level, value.var = "pct_aqu_of_resolved")[order(tau)])
cat("\neMLG consensus-only (4,010 units) pct_aqu_of_resolved by tau:\n")
print(emlg_cons_sweep[, .(tau, nres = toward_aqu + toward_pol, pct_aqu = round(pct_aqu_of_resolved,1))])
d08 <- flip[tau == 0.8]
cat(sprintf("\nAt tau=0.8: per-SNP %.0f%%, LD-thinned %.0f%%, per-eMLG %.0f%% toward aquilonia.\n",
            d08[level=="per-SNP", pct_aqu_of_resolved], d08[level=="SNP, LD-thinned", pct_aqu_of_resolved],
            d08[level=="per-eMLG", pct_aqu_of_resolved]))
cat("[direction-flip] done\n")
