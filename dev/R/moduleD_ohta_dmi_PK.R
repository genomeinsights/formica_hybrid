## =========================================================
## MODULE D (D1-D5) -- among-region two-locus Ohta test: do UNLINKED eMLG
## clusters show elevated AMONG-POPULATION two-locus LD -- the direct signature
## of epistasis / Dobzhansky-Muller incompatibilities (DMIs), as opposed to the
## generic per-locus ancestry sorting Modules A/B measured?
## =========================================================
## The INTRINSIC-side test of the headline question (C is the extrinsic side).
##
## LOGIC (why the AMONG-population component is the DMI signature HERE).
##   We have 20 INDEPENDENT, non-migrating replicate hybrid populations, each
##   independently sorting parental ancestry out of the same F.polyctena x
##   F.aquilonia admixture. For a pair of UNLINKED clusters (different chromosome):
##     - under independent neutral drift the two loci sort INDEPENDENTLY across
##       replicates, so the expected covariance of their per-population ancestry
##       frequencies is ~0;
##     - a DMI makes certain ancestry COMBINATIONS unfit, so the two loci sort in
##       CORRELATED directions across replicates -> non-zero among-population
##       covariance of allele frequencies.
##   PRIMARY STATISTIC = that among-population allele-frequency covariance, measured
##   by R_st = cor of the 20 per-population mean dosages, and confirmed by its exact
##   Ohta counterpart D2st (the among-population variance component). NB the run
##   settled a subtlety: Ohta's LITERAL D'2st (`Dp2st`, the systematic LD of the
##   averaged population) is NOT the same quantity and turns out ~0 here (large D2st
##   + Dp2st ~ 0 is Ohta's own signature of correlated allele-frequency divergence
##   among populations, i.e. structure/drift, rather than systematic epistatic LD).
##   So R_st / D2st are the primary axis; Dp2st ~ 0 and Ohta_D = D2is - D2st are
##   reported as decomposition findings, not ranking axes. The full decomposition
##   (D2it, D2is, D2st, Dp2st, Dp2is) is always stored.
##
## SCOPE (settled with PK). All-pairs over 32,840 has_eMLG clusters is ~5e8
##   comparisons -> infeasible for the exact Ohta dstat. Chosen design:
##     (i)  restrict to DIFFERENTIATED clusters (parent_maf >= 0.15; = 27,223).
##          This is an INFORMATIVENESS gate (a precondition for reading ancestry),
##          NOT the sorting OUTCOME, so it does not condition on the thing tested.
##     (ii) CHEAP PRE-FILTER over all ~3.5e8 unlinked differentiated pairs: R_st
##          from the 20 x K per-population mean-dosage matrix (continuous consensus,
##          no rounding) -- a fast linear proxy of the very statistic of interest.
##          Keep the top tail (|R_st| >= RST_THR).
##     (iii) EXACT Ohta decomposition (dstat_unphased_scan) only on the surviving
##          unlinked tail -> confirm the among-population component (D2st) and
##          characterise the within-population part (D2is) and the averaged LD (Dp2st).
##   LINKED (same-chromosome) pairs are retained as an internal POSITIVE CONTROL:
##   physical linkage SHOULD raise LD, so the interesting anomaly is UNLINKED pairs
##   carrying among-population LD.
##
## INFERENCE STATUS. This is a DESCRIPTIVE DISTRIBUTION-GENERATOR, not a candidate
##   caller. The run showed the raw among-population LD tail is essentially IDENTICAL
##   for unlinked and linked pairs (structure/shared-ancestry-axis dominated), so
##   NOTHING is pair-specific before the null and a candidate-DMI network here would
##   be premature (it was a 2,564-cluster hairball). Module D therefore emits the
##   observed R_st / D2st DISTRIBUTIONS (and the full per-pair decomposition), and
##   CANDIDATE DMIs are defined downstream as pairs whose among-population LD EXCEEDS
##   Module E's recombination-matched neutral null. The scan is factored into
##   `moduleD_pop_freq_matrix()`, `moduleD_prefilter()`, `moduleD_scan()` (top of
##   file) so Module E sources THIS file and runs the IDENTICAL pre-filter + scan on
##   simulated genotypes -- the null plugs in with no re-implementation. IMPORTANT
##   for E: apply the SAME RST_THR.
##
## Run from repo root. Reads canonical clustering, hybrids_only sample_data,
## moduleC_C3_cl.rds (the differentiated / DI gate). Writes data/moduleD_ohta.rds
## (per-pair decomposition + R_st distributions + top-R_st ranked view) and
## Figures/moduleD_fig4.{pdf,png}.

suppressMessages({ library(data.table); library(ggplot2); library(patchwork); library(parallel) })
source("dev/R/Ohta.R")   # ohta_fast_prepare(), dstat_unphased_scan() (uses parallel::mclapply)
set.seed(1)

## ---- PARAMETERS ---------------------------------------------------------
CLUSTERING <- "data/eMLG_5loci_0025_cM05.rds"
CL_GATE    <- "data/moduleC_C3_cl.rds"   # per-cluster differentiated / DI / sort_class
CORES      <- 8
DROP_SIELVA <- FALSE   # near-F1 pop (ancestry SD ~0.001): TRUE = robustness rerun.
                       # It mostly adds a ~0.5-frequency leverage point to the
                       # among-pop covariance; default keeps all 20 pops (as A/B/C).
## Pre-filter: how large an among-population frequency correlation sends a pair to
## the exact Ohta scan. |R_st| >= 0.7 on 20 populations is ~0.2% of pairs under
## independence -> ~7e5 candidate unlinked pairs (dstat feasible in minutes on
## CORES). Lower it to widen the net (more compute); Module E MUST reuse this value.
RST_THR    <- 0.7
PREFILTER_BLOCK <- 500          # cluster columns per crossprod block (~2 GB transient
                                # peak at K~27k; lower it if memory-constrained)
TOP_N      <- 2000              # size of the ranked top-R_st view saved for inspection
                                # (a convenience list, NOT a candidate call -- see header)
elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

## =========================================================================
## Reusable core (Module E sources these) -----------------------------------
## =========================================================================

## Per-population MEAN dosage matrix (pops x clusters), continuous consensus (no
## rounding -- frequencies want the full information). Columns that are monomorphic
## ACROSS populations (zero variance over the 20 pop means -> undefined correlation)
## are dropped and reported in `attr(,'dropped')`.
moduleD_pop_freq_matrix <- function(GT, pops, cols) {
  pops <- as.character(pops)
  lv   <- sort(unique(pops))
  ## suppressWarnings: a pop with all-NA members for a cluster yields NaN (that pair's
  ## R_st is then NaN and dropped downstream by is.finite) -- not an error.
  F <- suppressWarnings(vapply(lv, function(p) colMeans(GT[pops == p, cols, drop = FALSE], na.rm = TRUE),
              numeric(length(cols))))                     # length(cols) x nPops
  F <- t(F)                                               # nPops x nClusters
  colnames(F) <- cols
  sdv <- apply(F, 2, sd, na.rm = TRUE)
  keep <- is.finite(sdv) & sdv > 0
  out <- F[, keep, drop = FALSE]
  attr(out, "dropped") <- cols[!keep]
  out
}

## Blocked among-population correlation R_st = cor(F) over pops, returning
##   $edges : data.table(i, j, R_st) for pairs with |R_st| >= rst_thr, chr[i]!=chr[j],
##            i<j (indices into colnames(Fmat)); linked pairs are for the histogram
##            only and are NOT returned as edges.
##   $hist  : binned distribution of R_st for linked vs unlinked pairs (all pairs),
##            for the figure -- accumulated without materialising the full K x K matrix.
## Correlation across the (few) populations is computed as a crossproduct of
## column-standardised frequency vectors: scale each cluster's pop-vector to
## mean 0 / unit L2 norm, then crossprod == Pearson r.
moduleD_prefilter <- function(Fmat, chr, rst_thr = 0.7, block = 500,
                              hist_breaks = seq(-1, 1, by = 0.05)) {
  K  <- ncol(Fmat)
  Fs <- scale(Fmat, center = TRUE, scale = FALSE)         # centre over pops
  nrm <- sqrt(colSums(Fs^2))
  Fs <- sweep(Fs, 2, nrm, "/")                            # unit norm -> crossprod = r
  chr <- as.character(chr)
  nb  <- length(hist_breaks) - 1L
  h_link <- numeric(nb); h_unlink <- numeric(nb)
  edges <- vector("list", ceiling(K / block))
  bi <- 0L
  for (start in seq(1L, K, by = block)) {
    bi  <- bi + 1L
    idx <- start:min(start + block - 1L, K)
    Rb  <- crossprod(Fs[, idx, drop = FALSE], Fs)         # |idx| x K  (r of block vs all)
    ## histogram over off-diagonal pairs with global-j > global-i (each pair once)
    gi <- rep(idx, times = K)
    gj <- rep(seq_len(K), each = length(idx))
    rv <- as.vector(Rb)
    up <- gj > gi & is.finite(rv)
    if (any(up)) {
      link <- chr[gi[up]] == chr[gj[up]]
      hb <- cut(rv[up], hist_breaks, include.lowest = TRUE, labels = FALSE)
      h_link   <- h_link   + tabulate(hb[link],  nbins = nb)
      h_unlink <- h_unlink + tabulate(hb[!link], nbins = nb)
      ## surviving UNLINKED edges
      sel <- up & !link & abs(rv) >= rst_thr
      if (any(sel)) edges[[bi]] <- data.table(i = gi[sel], j = gj[sel], R_st = rv[sel])
    }
  }
  edges <- rbindlist(edges)
  list(edges = edges,
       hist  = data.table(mid = (head(hist_breaks, -1) + tail(hist_breaks, -1)) / 2,
                          linked = h_link, unlinked = h_unlink))
}

## Exact Ohta decomposition for a set of column-index pairs, on ROUNDED consensus
## dosages (the gametic tabulation needs integer 0/1/2 classes). `prep` is
## ohta_fast_prepare(round(units), pops); `pairs` is a 2-col matrix of column indices
## into prep$G. Gating is done upstream, so tot_maf = pop_maf = 0.
moduleD_scan <- function(prep, pairs, cores = 1) {
  dstat_unphased_scan(pairs = pairs, prep = prep, tot_maf = 0, pop_maf = 0, cores = cores)
}

## =========================================================================
## D1 -- assemble hybrid consensus units, populations, per-cluster metadata
## =========================================================================
clust  <- readRDS(CLUSTERING)
eMLG   <- clust$eMLG                       # 164 hybrids x 32,840 has_eMLG clusters, 0..2
groups <- clust$groups
cl     <- readRDS(CL_GATE)                 # per-cluster differentiated / DI / sort_class

e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1)
sample_data <- as.data.table(e1$sample_data)
pops_all <- sample_data[match(rownames(eMLG), Sample_ID), Population]
stopifnot(!anyNA(pops_all), length(pops_all) == nrow(eMLG))

if (DROP_SIELVA) {
  keep_ind <- pops_all != "Sielva"
  eMLG <- eMLG[keep_ind, , drop = FALSE]; pops_all <- pops_all[keep_ind]
}
message(sprintf("[D1] %d hybrid individuals x %d clusters | %d populations%s",
                nrow(eMLG), ncol(eMLG), length(unique(pops_all)),
                if (DROP_SIELVA) " (Sielva dropped)" else ""))

## chromosome + DI + differentiated flag, aligned to eMLG columns
chr_of  <- setNames(as.character(groups$Chr), groups$group_id)
cl_gate <- cl[, .(group_id, differentiated, DI, sort_class)]

## =========================================================================
## D2 -- scope: differentiated clusters (the informativeness gate)
## =========================================================================
diff_ids <- cl_gate[differentiated == TRUE, group_id]
scope    <- intersect(colnames(eMLG), diff_ids)          # order = eMLG column order
message(sprintf("[D2] differentiated clusters in scope: %d (of %d has_eMLG)",
                length(scope), ncol(eMLG)))

## =========================================================================
## D3 -- pre-filter (cheap R_st over ALL unlinked scope pairs) -> exact Ohta on tail
## =========================================================================
t0 <- Sys.time()
Fmat <- moduleD_pop_freq_matrix(eMLG, pops_all, scope)
dropped <- attr(Fmat, "dropped")
scope_kept <- colnames(Fmat)                              # monomorphic-across-pop dropped
chr_scope  <- chr_of[scope_kept]
message(sprintf("[D3] pop-frequency matrix %d x %d built (%d monomorphic-across-pop dropped) | %.0fs",
                nrow(Fmat), ncol(Fmat), length(dropped), elapsed(t0)))

t0 <- Sys.time()
pf <- moduleD_prefilter(Fmat, chr = chr_scope, rst_thr = RST_THR, block = PREFILTER_BLOCK)
n_pairs_total <- pf$hist[, sum(linked + unlinked)]
message(sprintf("[D3] pre-filter done: %s scope pairs scanned, %d unlinked survivors (|R_st|>=%.2f) | %.0fs",
                format(n_pairs_total, big.mark = ","), nrow(pf$edges), RST_THR, elapsed(t0)))

## exact Ohta decomposition on the surviving unlinked pairs
## prep on ROUNDED consensus (integer genotype classes), columns = scope_kept order
prep <- ohta_fast_prepare(round(eMLG[, scope_kept, drop = FALSE]), pops = pops_all)
t0 <- Sys.time()
ohta <- moduleD_scan(prep, pairs = as.matrix(pf$edges[, .(i, j)]), cores = CORES)
message(sprintf("[D3] exact Ohta scan on %d unlinked survivors | %.0fs", nrow(ohta), elapsed(t0)))

## assemble per-pair table
pairs_dt <- cbind(
  pf$edges,
  data.table(cluster1 = scope_kept[pf$edges$i], cluster2 = scope_kept[pf$edges$j],
             Chr1 = chr_scope[pf$edges$i], Chr2 = chr_scope[pf$edges$j]),
  ohta
)
pairs_dt[, Ohta_D := D2is - D2st]                         # prototype's secondary statistic
## dstat_unphased_scan returns Dp2st (= D'2st, among-pop) and Dp2is (= D'2is, within-pop)
## annotate with the gate metadata of each partner (DI / sort_class)
di_v <- setNames(cl_gate$DI, cl_gate$group_id); sc_v <- setNames(cl_gate$sort_class, cl_gate$group_id)
pairs_dt[, `:=`(DI1 = di_v[cluster1], DI2 = di_v[cluster2],
                sort1 = sc_v[cluster1], sort2 = sc_v[cluster2])]

## =========================================================================
## D4 -- descriptive summary (NO candidate network: deferred to Module E)
## =========================================================================
## The interesting quantity for E is the OBSERVED among-population LD distribution.
## Here we (i) quantify the unlinked-vs-linked tail contrast (is the among-pop LD
## linkage-specific, or a genome-wide ancestry-axis property?), (ii) summarise the
## Ohta decomposition on survivors (structure vs systematic-epistasis reading), and
## (iii) save a ranked top-R_st view for inspection. CANDIDATE DMIs = pairs whose
## among-population LD exceeds Module E's null; that call is NOT made here.

## (i) tail contrast: fraction of pairs above R_st cutoffs, unlinked vs linked
tailf <- function(t) pf$hist[mid >= t, .(unlinked = sum(unlinked), linked = sum(linked))]
tail_contrast <- rbindlist(lapply(c(0.5, 0.7, 0.9), function(t) {
  h <- tailf(t)
  data.table(R_st_cut = t,
             frac_unlinked = h$unlinked / pf$hist[, sum(unlinked)],
             frac_linked   = h$linked   / pf$hist[, sum(linked)])
}))[, ratio_link_unlink := frac_linked / frac_unlinked]

## (ii) decomposition medians on unlinked survivors
decomp <- pairs_dt[, .(D2is = median(D2is, na.rm = TRUE), D2st = median(D2st, na.rm = TRUE),
                       Dp2st = median(Dp2st, na.rm = TRUE), Dp2is = median(Dp2is, na.rm = TRUE),
                       Ohta_D = median(Ohta_D, na.rm = TRUE),
                       frac_D2st_gt_D2is = mean(D2st > D2is, na.rm = TRUE),
                       frac_Dp2st_gt_Dp2is = mean(Dp2st > Dp2is, na.rm = TRUE))]

## (iii) ranked top-R_st view (a convenience list for eyeballing, NOT a candidate call)
ranked <- head(pairs_dt[order(-abs(R_st))], TOP_N)[, .(
  cluster1, cluster2, Chr1, Chr2, R_st, D2st, D2is, Dp2st, Ohta_D, Fst_mean,
  DI1, DI2, sort1, sort2)]

cat("\n[D4] unlinked-vs-linked R_st tail contrast:\n"); print(tail_contrast)
cat("\n[D4] Ohta decomposition on unlinked survivors (medians):\n"); print(decomp)
cat(sprintf("\n[D4] reading: D2st >> D2is (frac %.2f) with Dp2st ~ 0 (frac Dp2st>Dp2is %.2f) => the\n",
            decomp$frac_D2st_gt_D2is, decomp$frac_Dp2st_gt_Dp2is))
cat("     among-population LD is correlated allele-frequency DIVERGENCE (structure/\n")
cat("     shared ancestry axis), not systematic epistasis. Excess-over-null = E.\n")

## =========================================================================
## D5 -- outputs + figure
## =========================================================================
saveRDS(list(pairs = pairs_dt, rst_hist = pf$hist, ranked = ranked,
             tail_contrast = tail_contrast, decomp = decomp,
             scope = scope_kept, dropped = dropped, n_pairs_total = n_pairs_total,
             params = list(RST_THR = RST_THR, DROP_SIELVA = DROP_SIELVA,
                           TOP_N = TOP_N, clustering = CLUSTERING)),
        "data/moduleD_ohta.rds")
message("[D5] saved data/moduleD_ohta.rds")

dir.create("Figures", showWarnings = FALSE)
th <- theme_classic(base_size = 8) +
  theme(plot.tag = element_text(face = "bold", size = 10),
        plot.title = element_text(size = 8.5, hjust = 0),
        axis.title = element_text(size = 8), axis.text = element_text(size = 6.5),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 7), legend.key.size = unit(3, "mm"),
        plot.margin = margin(4, 9, 2, 4))

## (a) THE KEY CONTRAST: among-population frequency correlation (R_st), unlinked vs
## linked, over ALL scope pairs. Linked (positive control) carries the physical-LD
## tail; an UNLINKED tail beyond it = candidate correlated among-population sorting.
hh <- melt(pf$hist, id.vars = "mid", variable.name = "pairtype", value.name = "n")
hh[, frac := n / sum(n), by = pairtype]
p4a <- ggplot(hh, aes(mid, frac, colour = pairtype)) +
  geom_vline(xintercept = c(-RST_THR, RST_THR), linetype = 2, colour = "grey70") +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = c(linked = "#D55E00", unlinked = "#0072B2")) +
  scale_y_sqrt() +
  labs(x = expression(R[st]~"(among-population freq. correlation)"),
       y = "fraction of pairs (sqrt)", title = "a  Unlinked vs linked among-pop LD") + th

## (b) Ohta decomposition on survivors: among-population (D2st) vs within-population
## (D2is). Points ABOVE the diagonal = among-pop LD dominates = structure/drift
## signature (correlated allele-freq divergence), not systematic within-deme epistasis.
p4b <- if (nrow(pairs_dt) > 0) {
  ggplot(pairs_dt, aes(D2is, D2st)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70") +
    geom_point(alpha = 0.25, size = 0.4, colour = "#0072B2") +
    labs(x = "within-pop  D2is", y = "among-pop  D2st",
         title = "b  Among- vs within-pop (D'2st≈0)") + th
} else patchwork::plot_spacer()

## (c) pre-filter validation: the cheap R_st proxy tracks the exact among-pop
## component D2st (so the genome-wide R_st distribution E is compared against is a
## faithful stand-in for the exact statistic).
p4c <- if (nrow(pairs_dt) > 0) {
  ggplot(pairs_dt, aes(abs(R_st), D2st)) +
    geom_point(alpha = 0.2, size = 0.4, colour = "#0072B2") +
    labs(x = expression("|"*R[st]*"|  (pre-filter proxy)"), y = "among-pop  D2st",
         title = "c  Proxy vs exact among-pop LD") + th
} else patchwork::plot_spacer()

fig4 <- p4a + p4b + p4c + plot_layout(widths = c(1.1, 1, 1)) +
  plot_annotation(tag_levels = NULL)
ggsave("Figures/moduleD_fig4.pdf", fig4, width = 210, height = 78, units = "mm")
ggsave("Figures/moduleD_fig4.png", fig4, width = 210, height = 78, units = "mm", dpi = 300)
cat("\n[done] data/moduleD_ohta.rds, Figures/moduleD_fig4.{pdf,png}\n")
