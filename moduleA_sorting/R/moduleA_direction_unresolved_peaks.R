## =========================================================
## MODULE A -- diagnosing the per-SNP "direction-unresolved" manhattan peaks
## =========================================================
## The binom sorting manhattan (moduleA_fig_sorting_manhattan.R) shows local
## peaks in the per-SNP BIDIRECTIONAL fraction at tau = 0.5. Several coincide with
## high-ld_w regions, so they could be single large LD clusters counted many times
## rather than genomic regions of genuinely heterogeneous ancestry outcomes. This
## script tests that, and answers four questions:
##   (1) Do the peaks persist when counted per LD-REDUCED UNIT (one representative
##       per cluster) instead of per SNP?
##   (2) What is the n_fixed distribution inside the peaks? (n_fixed = 10-11 means
##       prop_fixed 0.50-0.55, which drops out at the primary tau = 0.6 gate.)
##   (3) List the peak intervals.
##   (4) Do the same populations sort in opposite directions across loci within a
##       peak, and do different peaks share the same population split or not?
##
## Bidirectional is defined here at tau = 0.5 (prop_fixed >= 0.5 & p_binom >= ALPHA
## & n_fixed >= N_POW) -- the threshold at which the peaks are visible -- NOT the
## stored tau = 0.6 sort_class.
##
## Reads : moduleA_sorting/data/moduleA_snp.rds  (per-SNP binom stats: prop_fixed,
##         p_binom, n_fixed, uni_score, ...), module0_ld_pruning/data/
##         eMLG_5loci_0025_cM05.rds (LD clusters + representatives),
##         data/hybrids_and_parents_maf005.Rdata (genotypes, populations)
## Writes: moduleA_sorting/data/moduleA_direction_unresolved_peaks.csv (peak summary),
##         moduleA_sorting/Figures/moduleA_direction_unresolved_peaks.png (pop x peak split)
## Run from the repo root.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

TAU   <- 0.5                                # magnitude threshold at which peaks appear
ALPHA <- 0.05                              # binom direction-test significance
N_POW <- ceiling(log(ALPHA/2)/log(0.5))    # smallest testable n_fixed (=6 at 0.05)
FIXTH <- 0.15                              # per-population near-fixation tolerance (Module A)
WIN   <- 100000L                           # 100-kb windows
QUAL_FRAC <- 0.08; QUAL_NBI <- 3L; PEAK_MINBI <- 10L   # peak-calling thresholds
N_TOP_PEAKS <- 10L                         # peaks carried into the population analysis

## ---- per-SNP, direction-unresolved at tau = 0.5 ---------------------------------
snp <- as.data.table(readRDS("moduleA_sorting/data/moduleA_snp.rds"))
snp <- snp[differentiated == TRUE & is.finite(uni_score) & !is.na(p_binom)]
snp[, c("Chr", "Pos") := tstrsplit(marker, ":", fixed = TRUE)][, Pos := as.integer(Pos)]
snp[, `:=`(win = Pos %/% WIN,
           is_bi = prop_fixed >= TAU & p_binom >= ALPHA & n_fixed >= N_POW)]

## marker -> LD cluster + is-representative flag
groups <- readRDS("module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds")$groups
snp[groups[, .(marker = unlist(members)), by = group_id], on = "marker", group_id := i.group_id]
snp[, is_rep := marker %in% unique(groups$representative)]
e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
snp[as.data.table(e$map_hyb_005)[, .(marker, ld_w_095)], on = "marker", ld_w := i.ld_w_095]

## ---- (3) peaks: contiguous 100-kb windows enriched for direction-unresolved -----
w <- snp[, .(n = .N, n_bi = sum(is_bi)), by = .(Chr, win)][, bi_frac := n_bi / n]
w[, qual := n_bi >= QUAL_NBI & bi_frac >= QUAL_FRAC]
setorder(w, Chr, win)
w[, run := cumsum(!qual | c(TRUE, diff(win) > 2)), by = Chr]
peaks <- w[qual == TRUE, .(win0 = min(win), win1 = max(win), n_bi = sum(n_bi),
                           n_diff = sum(n)), by = .(Chr, run)][n_bi >= PEAK_MINBI]
peaks[, `:=`(mb0 = round(win0 * WIN / 1e6, 2), mb1 = round((win1 + 1) * WIN / 1e6, 2))]
peaks[, peak := sprintf("%s:%.1f-%.1f", Chr, mb0, mb1)]
setorder(peaks, -n_bi)

snp[, peak := NA_character_]
for (i in seq_len(nrow(peaks)))
  snp[Chr == peaks$Chr[i] & win >= peaks$win0[i] & win <= peaks$win1[i], peak := peaks$peak[i]]
bi <- snp[is_bi == TRUE & !is.na(peak)]

## ---- (1) per-SNP vs LD-reduced recount + LD context ----------------------
red <- bi[, .(n_bi_SNP = .N, n_clusters = uniqueN(group_id), n_bi_rep = sum(is_rep),
              LDinflation = round(.N / uniqueN(group_id), 1),
              mean_ldw = round(mean(ld_w, na.rm = TRUE), 2)), by = peak]
tab <- peaks[, .(peak, bi_frac = round(n_bi / n_diff, 3))][red, on = "peak"][order(-n_bi_SNP)]
fwrite(tab, "moduleA_sorting/data/moduleA_direction_unresolved_peaks.csv")
cat("=== (1)+(3) PEAK INTERVALS: per-SNP vs LD-reduced (representative) ===\n")
print(tab, row.names = FALSE)
cat(sprintf("\nTOTAL across peaks: %d per-SNP direction-unresolved -> %d distinct LD clusters -> %d direction-unresolved representatives\n",
            nrow(bi), uniqueN(bi$group_id), sum(bi$is_rep)))

## ---- (2) n_fixed distribution within peaks -------------------------------
cat("\n=== (2) n_fixed of direction-unresolved SNPs within peaks (n_obs = 20) ===\n")
print(bi[, .N, keyby = n_fixed])
cat(sprintf("  n_fixed in {10,11}: %d of %d (%.0f%%) -- these have prop_fixed 0.50-0.55 and drop out at tau = 0.6\n",
            bi[n_fixed %in% 10:11, .N], nrow(bi), 100 * bi[n_fixed %in% 10:11, .N] / nrow(bi)))

## ---- (4) population groupings (top peaks) --------------------------------
GT <- e$GTs_with_parents; sd <- as.data.table(e$sample_data_with_parents); sd <- sd[match(rownames(GT), Sample_ID)]
is_aqu_par <- sd$Population == "aquilonia_parent"; is_pol_par <- sd$Population == "polyctena_parent"
hybpop <- setdiff(unique(sd$Population), c("aquilonia_parent", "polyctena_parent"))
sel <- head(peaks$peak, N_TOP_PEAKS)
bip <- bi[peak %in% sel]; mk <- bip$marker

sub <- GT[, mk, drop = FALSE]                       # orient so dosage 2 = aquilonia allele
sgn <- sign(colMeans(sub[is_aqu_par, , drop = FALSE], na.rm = TRUE) -
            colMeans(sub[is_pol_par, , drop = FALSE], na.rm = TRUE))
for (j in seq_along(mk)) if (isTRUE(sgn[j] < 0)) sub[, j] <- 2 - sub[, j]
popfreq  <- sapply(hybpop, function(p) colMeans(sub[sd$Population == p, , drop = FALSE], na.rm = TRUE) / 2)
dir_call <- ifelse(popfreq >= 1 - FIXTH, 1L, ifelse(popfreq <= FIXTH, -1L, 0L))  # +1 aqu / -1 pol / 0
dc <- as.data.table(dir_call); dc[, marker := mk]; dc <- dc[bip[, .(marker, peak)], on = "marker"]
long <- melt(dc, id.vars = c("marker", "peak"), variable.name = "Population", value.name = "call")
grid <- long[, .(net = mean(call), n_aqu = sum(call == 1), n_pol = sum(call == -1)), by = .(peak, Population)]

## within-peak coherence: does every population sort the same way across the peak's loci?
coh <- long[call != 0, .(agree = max(mean(call == 1), mean(call == -1))), by = .(peak, Population)][
  , .(coherence = round(mean(agree), 2)), by = peak][order(-coherence)]
cat("\n=== (4) within-peak coherence (1 = each population always sorts the same way in the peak) ===\n")
print(coh, row.names = FALSE)

## cross-peak: do different peaks share the same population split?
m <- as.matrix(dcast(grid, Population ~ peak, value.var = "net")[, -1])
cat("\ncross-peak correlation of population net-direction vectors (top peaks):\n")
print(round(cor(m, use = "pairwise.complete.obs"), 2))

## heatmap: population x peak, coloured by net direction
ord <- grid[, .(m = mean(net, na.rm = TRUE)), by = Population][order(m), Population]
grid[, `:=`(Population = factor(Population, levels = ord), peak = factor(peak, levels = sel))]
p <- ggplot(grid, aes(peak, Population, fill = net)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%d/%d", n_aqu, n_pol)), size = 2.3) +
  scale_fill_gradient2(low = "#d95f0e", mid = "#f7f7f7", high = "#2c7fb8", midpoint = 0,
                       limits = c(-1, 1), name = "net\n(+aqu / -pol)") +
  labs(title = "Population direction within per-SNP direction-unresolved peaks (tau = 0.5)",
       subtitle = "tile = mean per-locus fixation direction; label = #loci aquilonia-fixed / polyctena-fixed",
       x = "peak (direction-unresolved region)", y = "hybrid population (ordered by mean net)") +
  theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("moduleA_sorting/Figures/moduleA_direction_unresolved_peaks.png", p, width = 9, height = 6, dpi = 150)
cat("\nSaved: moduleA_sorting/Figures/moduleA_direction_unresolved_peaks.png,",
    "moduleA_sorting/data/moduleA_direction_unresolved_peaks.csv\n")
