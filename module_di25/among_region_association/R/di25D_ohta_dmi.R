## =========================================================================
## DI25 among-region two-locus Ohta test -- HIGH-DI-ONLY re-run of Module D's Ohta arm.
## =========================================================================
## Same question as moduleD_ohta_dmi.R, restricted to the high-DI diagnostic markers
## and their FROM-SCRATCH DI25 LD partition (di25_clustering_cM5): do UNLINKED DI25
## units show elevated AMONG-POPULATION two-locus LD -- the signature of epistasis /
## Dobzhansky-Muller incompatibilities -- beyond the per-locus ancestry sorting?
##
## LOGIC. 20 independent replicate hybrid populations sort the same admixture. For an
## UNLINKED pair, independent neutral drift -> ~0 SYSTEMATIC among-population covariance
## of ancestry frequencies; a DMI makes some ancestry combinations unfit -> correlated
## sorting in the SAME direction across replicates -> non-zero systematic covariance.
##
## PRIMARY STATISTIC = R_st (cor of the 20 per-population mean dosages). Its exact Ohta
## counterpart is D'2st (`Dp2st`), the SYSTEMATIC among-population LD (LD of the
## population-averaged gamete frequencies), NOT D2st. This was measured directly on this
## data: over the full range of random unlinked pairs cor(R_st^2, Dp2st) = 0.74 while
## cor(R_st^2, D2st) = 0.02. D2st is the among-population *co-differentiation* variance
## (product of each locus's marginal freqs vs total); with every DI25 unit strongly
## structured it is large and ~constant (~0.16) regardless of association, so it is NOT
## the epistasis axis. The epistasis-relevant axis is the SYSTEMATIC one (Dp2st ~ R_st^2),
## and here it is ~0 (median 0.008) while the co-differentiation D2st is large -- i.e.
## drift/structure, no systematic parallel LD, pre-null.
##
## OHTA DRIFT-vs-SELECTION HEURISTIC (compare WITHIN each family; the primed terms carry
## a /4 so raw D2st-vs-Dp2st magnitude is not comparable). Epistatic selection: D2is>D2st
## AND Dp2st>Dp2is (LD within demes, systematic/directional). Drift/structure: D2st>D2is
## AND Dp2is>Dp2st. Stored as Ohta_D = D2is-D2st and Ohta_Dprime = Dp2st-Dp2is (both > 0
## under selection). This data: D2st>D2is and Dp2is>Dp2st -> both families read drift.
## NB even an elevated Dp2st does not by itself prove selection -- neutral admixture with
## shared founding also makes systematic among-pop ancestry LD -- so Dp2st/R_st is the
## SCREENING axis and the call is deferred to the recombination-matched neutral null.
##
## DIFFERENCES FROM Module D's canonical run (all baked into the read-out):
##   * UNITS = di25D_units.rds (165 hybrids x 11,052 high-DI units, aqu-oriented). ALL
##     units enter -- eMLG CONSENSUS for >2-marker clusters AND the REPRESENTATIVE SNP
##     for 1-2-marker clusters (the full DI25 LD partition), not only the has_eMLG set.
##   * The `differentiated` gate is near-vacuous here (11,042/11,052 TRUE) because every
##     DI25 marker is already ancestry-diagnostic; it drops only the 10 low-parent-MAF
##     units. Scope is therefore ~all units, so the pre-filter runs over ~all pairs.
##   * Ohta functions come from LDscnR (devtools::load_all), as in Module A -- no local
##     Ohta.R copy.
##
## INFERENCE STATUS (unchanged). DESCRIPTIVE distribution-generator, NOT a candidate
## caller: it emits the observed R_st / D'2st distributions + the full per-pair Ohta
## decomposition, and -- for THIS pipeline -- supplies the extreme +/-R_st pairs that
## seed the EMMAX arm's focal set (di25D_emmax.R), mirroring how moduleD_emmax.R read
## moduleD_ohta.rds. Candidate DMIs are defined downstream vs a recomb-matched null.
##
## Run from repo root. Reads di25D_units.rds; writes data/di25D_ohta.rds (per-pair
## decomposition + R_st distributions + top-R_st ranked view) and Figures/di25D_fig4.{pdf,png}.

suppressMessages({ library(data.table); library(ggplot2); library(patchwork); library(parallel) })
devtools::load_all("~/gitlab/LDscnR/")                       # ohta_fast_prepare(), dstat_unphased_scan()
source("module_di25/among_region_association/R/di25D_paralogy.R")   # flag_paralogy(), moduleD_cluster_het()
set.seed(1)

## ---- PARAMETERS ---------------------------------------------------------
UNITS      <- "module_di25/among_region_association/data/di25D_units.rds"
OUT        <- "module_di25/among_region_association/data/di25D_ohta.rds"
FIG        <- "module_di25/among_region_association/Figures/di25D_fig4"
CORES      <- 8
LINK_CM    <- 10                # unlinked = different chr OR same chr > LINK_CM cM (matches Module D)
DROP_SIELVA <- FALSE            # near-F1 pop; TRUE = robustness rerun (default keeps all 20 pops)
RST_THR    <- 0.7               # |R_st| that sends an unlinked pair to the exact Ohta scan (reuse in E)
PREFILTER_BLOCK <- 500          # cluster columns per crossprod block
TOP_N      <- 2000              # size of the ranked top-R_st view saved for inspection
PARALOGY_R <- 0.9               # |within-pop r| above which a survivor pair is a duplicate
elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

## =========================================================================
## Reusable core (identical algorithm to Module D; kept local for a self-contained
## module and any future DI25 Module-E null). ------------------------------
## =========================================================================

## Per-population MEAN dosage matrix (pops x units), continuous consensus (no rounding).
## Columns monomorphic ACROSS populations (zero variance over pop means) are dropped and
## reported in attr(,'dropped').
di25D_pop_freq_matrix <- function(GT, pops, cols) {
  pops <- as.character(pops); lv <- sort(unique(pops))
  F <- suppressWarnings(vapply(lv, function(p) colMeans(GT[pops == p, cols, drop = FALSE], na.rm = TRUE),
              numeric(length(cols))))
  F <- t(F); colnames(F) <- cols
  sdv <- apply(F, 2, sd, na.rm = TRUE); keep <- is.finite(sdv) & sdv > 0
  out <- F[, keep, drop = FALSE]; attr(out, "dropped") <- cols[!keep]; out
}

## Blocked among-population correlation R_st = cor(F) over pops. Linked = same chr AND
## |dcM| <= link_cm; otherwise unlinked. Returns $edges (unlinked |R_st|>=rst_thr, i<j),
## $hist (linked vs unlinked R_st distribution), $decay (same-chr mean R_st^2 by cM) and
## $baseline_rst2 (inter-chromosomal mean). No K x K materialisation.
di25D_prefilter <- function(Fmat, chr, cm, rst_thr = 0.7, link_cm = 10, block = 500,
                            hist_breaks = seq(-1, 1, by = 0.05),
                            cm_breaks = c(seq(0, 40, by = 2), Inf)) {
  K  <- ncol(Fmat)
  Fs <- scale(Fmat, center = TRUE, scale = FALSE)
  nrm <- sqrt(colSums(Fs^2)); Fs <- sweep(Fs, 2, nrm, "/")
  chr <- as.character(chr); cm <- as.numeric(cm)
  nb  <- length(hist_breaks) - 1L; nbc <- length(cm_breaks) - 1L
  h_link <- numeric(nb); h_unlink <- numeric(nb)
  s_rst2 <- numeric(nbc); n_cm <- numeric(nbc); base_s <- 0; base_n <- 0
  binsum <- function(x, b) { o <- numeric(nbc); t <- tapply(x, b, sum); o[as.integer(names(t))] <- t; o }
  edges <- vector("list", ceiling(K / block)); bi <- 0L
  for (start in seq(1L, K, by = block)) {
    bi  <- bi + 1L; idx <- start:min(start + block - 1L, K)
    Rb  <- crossprod(Fs[, idx, drop = FALSE], Fs)
    gi <- rep(idx, times = K); gj <- rep(seq_len(K), each = length(idx)); rv <- as.vector(Rb)
    up <- gj > gi & is.finite(rv)
    if (any(up)) {
      gi <- gi[up]; gj <- gj[up]; rv <- rv[up]
      same <- chr[gi] == chr[gj]; dcm <- abs(cm[gi] - cm[gj])
      link <- same & is.finite(dcm) & dcm <= link_cm
      hb <- cut(rv, hist_breaks, include.lowest = TRUE, labels = FALSE)
      h_link   <- h_link   + tabulate(hb[link],  nbins = nb)
      h_unlink <- h_unlink + tabulate(hb[!link], nbins = nb)
      sel <- !link & abs(rv) >= rst_thr
      if (any(sel)) edges[[bi]] <- data.table(i = gi[sel], j = gj[sel], R_st = rv[sel])
      sc <- same & is.finite(dcm)
      if (any(sc)) { cb <- cut(dcm[sc], cm_breaks, include.lowest = TRUE, right = FALSE, labels = FALSE)
        r2 <- rv[sc]^2; s_rst2 <- s_rst2 + binsum(r2, cb); n_cm <- n_cm + tabulate(cb, nbins = nbc) }
      ic <- !same; if (any(ic)) { base_s <- base_s + sum(rv[ic]^2); base_n <- base_n + sum(ic) }
    }
  }
  edges <- rbindlist(edges)
  list(edges = edges,
       hist  = data.table(mid = (head(hist_breaks, -1) + tail(hist_breaks, -1)) / 2,
                          linked = h_link, unlinked = h_unlink),
       decay = data.table(cm_lo = head(cm_breaks, -1), cm_hi = tail(cm_breaks, -1),
                          mean_rst2 = s_rst2 / pmax(n_cm, 1), n = n_cm),
       baseline_rst2 = if (base_n > 0) base_s / base_n else NA_real_)
}

## Exact Ohta decomposition for column-index pairs, on ROUNDED consensus dosages.
di25D_scan <- function(prep, pairs, cores = 1)
  dstat_unphased_scan(pairs = pairs, prep = prep, tot_maf = 0, pop_maf = 0, cores = cores)

## =========================================================================
## D1 -- load the assembled DI25 units, populations, per-unit metadata
## =========================================================================
u <- readRDS(UNITS)
eMLG   <- u$dosage                     # 165 hybrids x 11,052 high-DI units, aqu-oriented 0..2
groups <- u$groups
cl_gate <- u$gate[, .(group_id, differentiated, DI, sort_class)]
pops_all <- as.character(u$pops)
chr_of <- u$chr_of; cm_of <- u$cm_of; marker_Ho <- u$marker_Ho
stopifnot(length(pops_all) == nrow(eMLG), !anyNA(pops_all))

if (DROP_SIELVA) {
  keep_ind <- pops_all != "Sielva"
  eMLG <- eMLG[keep_ind, , drop = FALSE]; pops_all <- pops_all[keep_ind]
}
message(sprintf("[D1] %d hybrid individuals x %d DI25 units | %d populations%s",
                nrow(eMLG), ncol(eMLG), length(unique(pops_all)),
                if (DROP_SIELVA) " (Sielva dropped)" else ""))
message(sprintf("[D1] genetic positions attached (%d units with finite cM)", sum(is.finite(cm_of))))

## =========================================================================
## D2 -- scope: differentiated units (near-vacuous here; drops only low-parent-MAF units)
## =========================================================================
diff_ids <- cl_gate[differentiated == TRUE, group_id]
scope    <- intersect(colnames(eMLG), diff_ids)          # order = dosage column order
message(sprintf("[D2] differentiated units in scope: %d (of %d units)", length(scope), ncol(eMLG)))

## =========================================================================
## D3 -- pre-filter (cheap R_st over ALL unlinked scope pairs) -> exact Ohta on tail
## =========================================================================
t0 <- Sys.time()
Fmat <- di25D_pop_freq_matrix(eMLG, pops_all, scope)
dropped <- attr(Fmat, "dropped"); scope_kept <- colnames(Fmat)
chr_scope <- chr_of[scope_kept]; cm_scope <- cm_of[scope_kept]
message(sprintf("[D3] pop-frequency matrix %d x %d built (%d monomorphic-across-pop dropped) | %.0fs",
                nrow(Fmat), ncol(Fmat), length(dropped), elapsed(t0)))

t0 <- Sys.time()
pf <- di25D_prefilter(Fmat, chr = chr_scope, cm = cm_scope, rst_thr = RST_THR,
                      link_cm = LINK_CM, block = PREFILTER_BLOCK)
n_pairs_total <- pf$hist[, sum(linked + unlinked)]
message(sprintf("[D3] unlinked = different chr or same chr > %g cM; %s linked (<=%g cM) reference pairs",
                LINK_CM, format(pf$hist[, sum(linked)], big.mark = ","), LINK_CM))
message(sprintf("[D3] pre-filter done: %s scope pairs scanned, %d unlinked survivors (|R_st|>=%.2f) | %.0fs",
                format(n_pairs_total, big.mark = ","), nrow(pf$edges), RST_THR, elapsed(t0)))

## exact Ohta decomposition on the surviving unlinked pairs (ROUNDED consensus classes)
prep <- ohta_fast_prepare(round(eMLG[, scope_kept, drop = FALSE]), pops = pops_all)
t0 <- Sys.time()
ohta <- di25D_scan(prep, pairs = as.matrix(pf$edges[, .(i, j)]), cores = CORES)
message(sprintf("[D3] exact Ohta scan on %d unlinked survivors | %.0fs", nrow(ohta), elapsed(t0)))

## assemble per-pair table
pairs_dt <- cbind(
  pf$edges,
  data.table(cluster1 = scope_kept[pf$edges$i], cluster2 = scope_kept[pf$edges$j],
             Chr1 = chr_scope[pf$edges$i], Chr2 = chr_scope[pf$edges$j]),
  ohta)
pairs_dt[, Ohta_D := D2is - D2st]         # unprimed contrast (> 0 under epistatic selection)
pairs_dt[, Ohta_Dprime := Dp2st - Dp2is]  # PRIMED/systematic contrast (> 0 under selection) -- the R_st axis
di_v <- setNames(cl_gate$DI, cl_gate$group_id); sc_v <- setNames(cl_gate$sort_class, cl_gate$group_id)
pairs_dt[, `:=`(DI1 = di_v[cluster1], DI2 = di_v[cluster2],
                sort1 = sc_v[cluster1], sort2 = sc_v[cluster2])]

## cross-chromosome paralogy / duplication filter (shared with the EMMAX arm)
t0 <- Sys.time()
het_of <- moduleD_cluster_het(groups, scope_kept, marker_Ho)
pairs_dt <- flag_paralogy(pairs_dt, "cluster1", "cluster2", eMLG, pops_all,
                          het_of = het_of, thr = PARALOGY_R, cores = CORES)
message(sprintf("[D3] paralogy filter: %d/%d survivor pairs flagged as duplicates (|within-pop r| > %.2f) | %.0fs",
                sum(pairs_dt$paralog), nrow(pairs_dt), PARALOGY_R, elapsed(t0)))

## =========================================================================
## D4 -- descriptive summary (NO candidate network: deferred to the null)
## =========================================================================
tailf <- function(t) pf$hist[mid >= t, .(unlinked = sum(unlinked), linked = sum(linked))]
tail_contrast <- rbindlist(lapply(c(0.5, 0.7, 0.9), function(t) {
  h <- tailf(t)
  data.table(R_st_cut = t, frac_unlinked = h$unlinked / pf$hist[, sum(unlinked)],
             frac_linked = h$linked / pf$hist[, sum(linked)])
}))[, ratio_link_unlink := frac_linked / frac_unlinked]

decomp <- pairs_dt[, .(D2is = median(D2is, na.rm = TRUE), D2st = median(D2st, na.rm = TRUE),
                       Dp2st = median(Dp2st, na.rm = TRUE), Dp2is = median(Dp2is, na.rm = TRUE),
                       Ohta_D = median(Ohta_D, na.rm = TRUE), Ohta_Dprime = median(Ohta_Dprime, na.rm = TRUE),
                       frac_D2st_gt_D2is = mean(D2st > D2is, na.rm = TRUE),      # >0.5 -> drift (unprimed)
                       frac_Dp2st_gt_Dp2is = mean(Dp2st > Dp2is, na.rm = TRUE))] # >0.5 -> selection (primed)

ranked <- head(pairs_dt[order(-abs(R_st))], TOP_N)[, .(
  cluster1, cluster2, Chr1, Chr2, R_st, Dp2st, Dp2is, D2st, D2is, Ohta_Dprime, Ohta_D, Fst_mean,
  DI1, DI2, sort1, sort2, within_pop_r, concordance, paralog)]

cat("\n[D4] unlinked-vs-linked R_st tail contrast:\n"); print(tail_contrast)
cat("\n[D4] Ohta decomposition on unlinked survivors (medians):\n"); print(decomp)
cat(sprintf("\n[D4] Ohta drift-vs-selection: unprimed %s (D2st>D2is frac %.2f); primed %s (Dp2st>Dp2is frac %.2f)\n",
            ifelse(decomp$frac_D2st_gt_D2is   > 0.5, "DRIFT", "selection"), decomp$frac_D2st_gt_D2is,
            ifelse(decomp$frac_Dp2st_gt_Dp2is > 0.5, "selection", "DRIFT"), decomp$frac_Dp2st_gt_Dp2is))
cat("     epistasis axis = SYSTEMATIC Dp2st (~ R_st^2); ~0 here -> the co-differentiation\n")
cat("     D2st is large under structure, so this is drift/structure, NOT parallel epistatic LD.\n")
cat(sprintf("\n[D4] among-pop LD (mean R_st^2) vs genetic distance | inter-chromosomal baseline = %.4f:\n",
            pf$baseline_rst2))
print(pf$decay[is.finite(cm_hi)][, .(cM = sprintf("%g-%g", cm_lo, cm_hi),
                                     mean_rst2 = round(mean_rst2, 4), n)])

## =========================================================================
## D5 -- outputs + figure
## =========================================================================
saveRDS(list(pairs = pairs_dt, rst_hist = pf$hist, ranked = ranked,
             tail_contrast = tail_contrast, decomp = decomp,
             rst_decay = pf$decay, baseline_rst2 = pf$baseline_rst2,
             scope = scope_kept, dropped = dropped, n_pairs_total = n_pairs_total,
             params = list(RST_THR = RST_THR, LINK_CM = LINK_CM, DROP_SIELVA = DROP_SIELVA,
                           TOP_N = TOP_N, units = UNITS)),
        OUT)
message("[D5] saved ", OUT)

dir.create(dirname(FIG), showWarnings = FALSE)
th <- theme_classic(base_size = 8) +
  theme(plot.tag = element_text(face = "bold", size = 10),
        plot.title = element_text(size = 8.5, hjust = 0),
        axis.title = element_text(size = 8), axis.text = element_text(size = 6.5),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 7), legend.key.size = unit(3, "mm"),
        plot.margin = margin(4, 9, 2, 4))

hh <- melt(pf$hist, id.vars = "mid", variable.name = "pairtype", value.name = "n")
hh[, frac := n / sum(n), by = pairtype]
p4a <- ggplot(hh, aes(mid, frac, colour = pairtype)) +
  geom_vline(xintercept = c(-RST_THR, RST_THR), linetype = 2, colour = "grey70") +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = c(linked = "#D55E00", unlinked = "#0072B2"),
                      labels = c(linked = sprintf("linked (<=%g cM)", LINK_CM), unlinked = "unlinked")) +
  scale_y_sqrt() +
  labs(x = expression(R[st]~"(among-population freq. correlation)"),
       y = "fraction of pairs (sqrt)", title = "a  Unlinked vs linked among-pop LD (DI25)") + th

## Epistasis axis = the SYSTEMATIC (primed) family: Dp2st (systematic among-pop LD, the
## R_st counterpart) vs Dp2is (non-systematic scatter). Points ABOVE the diagonal =
## systematic LD dominates = selection-like; here they hug the x-axis (Dp2st≈0) = drift.
p4b <- if (nrow(pairs_dt) > 0) {
  ggplot(pairs_dt, aes(Dp2is, Dp2st)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70") +
    geom_point(alpha = 0.25, size = 0.4, colour = "#0072B2") +
    labs(x = "non-systematic  D'2is", y = "systematic  D'2st  (~ R_st^2)",
         title = "b  Systematic vs non-systematic LD (D'2st≈0)") + th
} else patchwork::plot_spacer()

p4c <- if (pf$decay[is.finite(cm_hi) & n > 0, .N] > 0) {
  dd <- pf$decay[is.finite(cm_hi) & n > 0][, cm_mid := (cm_lo + cm_hi) / 2]
  ggplot(dd, aes(cm_mid, mean_rst2)) +
    geom_hline(yintercept = pf$baseline_rst2, linetype = 3, colour = "grey50") +
    geom_vline(xintercept = LINK_CM, linetype = 2, colour = "grey70") +
    geom_line(linewidth = 0.5, colour = "#D55E00") + geom_point(size = 0.7, colour = "#D55E00") +
    annotate("text", x = LINK_CM + 1, y = max(dd$mean_rst2), hjust = 0, size = 2.2,
             label = "inter-chr baseline", colour = "grey50", vjust = 3) +
    labs(x = "genetic distance (cM)", y = expression("among-pop LD  ("*R[st]^2*")"),
         title = "c  LD decay vs cM (same chr)") + th
} else patchwork::plot_spacer()

fig4 <- p4a + p4b + p4c + plot_layout(widths = c(1.1, 1, 1))
ggsave(paste0(FIG, ".pdf"), fig4, width = 210, height = 78, units = "mm")
ggsave(paste0(FIG, ".png"), fig4, width = 210, height = 78, units = "mm", dpi = 300)
cat(sprintf("\n[done] %s, %s.{pdf,png} | %d unlinked survivor pairs (%d clean)\n",
            OUT, FIG, nrow(pairs_dt), if (nrow(pairs_dt) > 0) sum(!pairs_dt$paralog) else 0L))
