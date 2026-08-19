## =============================================================================
## module_di25 -- is the empirical ld_w EXCESS over the neutral sim associated
##                with sorting and/or BDMI regions?
## =============================================================================
## residual = empirical windowed-median ld_w_0.95  -  across-1000-rep simulated
## median. Because both sides share the marker panel and recombination landscape,
## a POSITIVE residual = empirical local LD in excess of what neutral drift on the
## same map produces (the "points above the 1:1 line" in the scatter).
##
## For each 500kb window we attach:
##   * sorted_frac -- fraction of differentiated diagnostic SNPs that are ancestry
##                    SORTED (any direction) at tau=0.6 (module-canonical classify_sort)
##   * bdmi_frac   -- fraction of diagnostic SNPs inside a BDMI candidate interval
##                    (colleague's liftoff BEDs), at a chosen X^2 cutoff
## and ask whether residual rises with either.
##
## Caveat printed with the result: windows are spatially autocorrelated and both
## excess LD and sorting concentrate in the same big low-recombination blocks, so
## the Spearman is descriptive; a per-chromosome circular-rotation null (as in
## bdmi_sorting_overlap.R) would be the formal test.
##
## Run from the formica_hybrid repo root:
##   Rscript module_di25/R/di25_ldw_residual_vs_sorting_bdmi.R [bdmi_cutoff_tag]
##     bdmi_cutoff_tag default "13_003"  (moderate; 215 intervals)
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()

args    <- commandArgs(trailingOnly = TRUE)
CUT     <- if (length(args) >= 1) args[1] else "13_003"
ENV     <- "module_di25/data/di25_ldw_envelope.rds"
SORTRDS <- "module_di25/data/di25_sorting_snp.rds"
BED     <- sprintf("data/liftoff_Frufa_DTOL_PR/bdmi_candidates.cutoff_%s.liftoff.Frufa_DTOL_PR.bed", CUT)
OUTPNG  <- "module_di25/Figures/di25_ldw_residual_vs_sorting_bdmi.png"
OUTPDF  <- sub("\\.png$", ".pdf", OUTPNG)
WIN_BP  <- 500000L
TAU     <- 0.6

## ---- windows with empirical + simulated-median ld_w ------------------------
pl <- as.data.table(readRDS(ENV)$pl)[is.finite(emp) & is.finite(sim_med)]
pl[, residual := emp - sim_med]

## ---- per-window sorted fraction (tau=0.6) ----------------------------------
ps  <- as.data.table(readRDS(SORTRDS))
ps[, `:=`(chr_id = as.integer(sub("Chr", "", sub(":.*", "", marker))),
          pos    = as.integer(sub(".*:", "", marker)))]
ok  <- ps$differentiated & ps$n_obs > 0 & is.finite(ps$uni_score)
ps[, sorted := FALSE]
ps[ok, sorted := classify_sort(n_aqu, n_pol, n_obs, sort_th = TAU,
                               sort_rule = "binom", alpha = 0.05) != "unsorted"]

## ---- per-SNP BDMI membership -----------------------------------------------
bdmi <- fread(BED, header = FALSE, col.names = c("chr", "start", "end"))
bdmi[, chr_id := as.integer(sub("chromosome_", "", chr))]
setkey(bdmi, chr_id, start, end)
ps[, `:=`(s = pos, e = pos)]
ov  <- foverlaps(ps[, .(chr_id, s, e, marker)], bdmi[, .(chr_id, start, end)],
                 by.x = c("chr_id", "s", "e"), type = "within", nomatch = 0L)
ps[, in_bdmi := marker %in% ov$marker]

## ---- aggregate to windows and merge ----------------------------------------
ps[, win := pos %/% WIN_BP]
wagg <- ps[, .(n_snp = .N, n_diff = sum(differentiated),
               sorted_frac = sum(sorted) / max(1L, sum(differentiated)),
               bdmi_frac   = mean(in_bdmi)),
           by = .(chr_id, win)]
d <- merge(pl, wagg, by = c("chr_id", "win"), all.x = TRUE)
d[is.na(sorted_frac), sorted_frac := 0]; d[is.na(bdmi_frac), bdmi_frac := 0]

## ---- association stats ------------------------------------------------------
sp <- function(x, y) cor(x, y, method = "spearman", use = "complete.obs")
above <- d[residual > 0]; below <- d[residual <= 0]
d[, excess := emp > sim_hi]                          # above the 95% envelope (sim_hi carried from pl)

cat(sprintf("\n=== ld_w residual (emp - sim_median) vs sorting / BDMI [cutoff %s] ===\n", CUT))
cat(sprintf("windows: %d  (above 1:1 line: %d, below: %d; above 95%% envelope: %d)\n",
            nrow(d), nrow(above), nrow(below), d[, sum(excess)]))
cat(sprintf("Spearman(residual, sorted_frac) = %+.3f\n", sp(d$residual, d$sorted_frac)))
cat(sprintf("Spearman(residual, bdmi_frac)   = %+.3f\n", sp(d$residual, d$bdmi_frac)))
cat("-- restricted to non-saturated windows (sim_med < 0.5), removing the big LD blocks --\n")
dl <- d[sim_med < 0.5]
cat(sprintf("   n=%d  Spearman(residual, sorted_frac) = %+.3f | (residual, bdmi_frac) = %+.3f\n",
            nrow(dl), sp(dl$residual, dl$sorted_frac), sp(dl$residual, dl$bdmi_frac)))
grp <- d[, .(n = .N, mean_sorted = mean(sorted_frac), mean_bdmi = mean(bdmi_frac)),
         by = .(above_line = residual > 0)]
cat("-- group means --\n"); print(grp)
cat("-- above vs within/below the 95% envelope --\n")
print(d[, .(n = .N, mean_sorted = mean(sorted_frac), mean_bdmi = mean(bdmi_frac)), by = excess])
cat("\nNOTE: windows are spatially autocorrelated; this is descriptive. A per-chromosome\n",
    "circular-rotation null (cf. bdmi_sorting_overlap.R) is the formal enrichment test.\n")

## ---- figure: the scatter, coloured by sorting then by BDMI -----------------
base <- function(fill, name, lab)
  ggplot(d, aes(sim_med, emp)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
    geom_point(aes(colour = .data[[fill]]), size = 1.1, alpha = 0.85, stroke = 0) +
    scale_colour_viridis_c(name = name, option = "C") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = expression("simulated median  " * ld[w][0.95]),
         y = expression("empirical  " * ld[w][0.95]), title = lab) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
                                     legend.position = "right")
p <- base("sorted_frac", "sorted\nfraction", "coloured by ancestry-sorted fraction") +
     base("bdmi_frac", "BDMI\nfraction", sprintf("coloured by BDMI overlap (cutoff %s)", CUT)) +
     plot_annotation(title = "Empirical ld_w excess over neutral sim vs sorting / BDMI",
                     subtitle = "one point per 500kb window; points above the dashed 1:1 line = empirical LD above the neutral expectation")
ggsave(OUTPDF, p, width = 11, height = 5.3)
ggsave(OUTPNG, p, width = 11, height = 5.3, dpi = 150)
cat("saved:", OUTPNG, "\n")
