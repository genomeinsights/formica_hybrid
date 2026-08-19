## =============================================================================
## module_di25 -- is the empirical ld_w EXCESS over neutral associated with
##                genotype-CLIMATE association (local adaptation)?
## =============================================================================
## residual = empirical windowed-median ld_w_0.95 - across-1000-rep simulated
## median  (empirical local LD beyond neutral drift on the same map).
##
## Per-window climate-association strength from the module B BayPass GEA
## (moduleB_snp_bf.rds): per-SNP Bayes Factors on the two raw environmental PCs
## (BF1 = climate PC1, BF2 = climate PC2; structure-controlled via Omega). Each
## 500kb window gets clim = max over its SNPs of max(BF1, BF2) [dB], the strongest
## genotype-climate association in the window (n15 = count of SNPs with BF >= 15).
##
## Test: does residual (excess LD) rise with clim? Spearman + per-chromosome
## circular-rotation null (phase-randomise the climate track vs the fixed residual
## track, within chromosome, two-sided). For context, the same correlation is run
## for raw empirical ld_w and for the neutral simulated ld_w -- if climate tracks
## the RESIDUAL specifically (and not the neutral landscape), that points at
## adaptation rather than architecture.
##
## Run from the formica_hybrid repo root:
##   Rscript module_di25/R/di25_ldw_residual_vs_climate.R
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
set.seed(1)
ENV     <- "module_di25/data/di25_ldw_envelope.rds"
BF      <- "moduleB_climate_GEA/data/moduleB_snp_bf.rds"
OUTPNG  <- "module_di25/Figures/di25_ldw_residual_vs_climate.png"
OUTPDF  <- sub("\\.png$", ".pdf", OUTPNG)
WIN_BP  <- 500000L
N_PERM  <- 2000L
BF_STRONG <- 15

## ---- residual per window ---------------------------------------------------
pl <- as.data.table(readRDS(ENV)$pl)[is.finite(emp) & is.finite(sim_med)]
pl[, `:=`(residual = emp - sim_med, excess = emp > sim_hi)]

## ---- per-window climate-association strength -------------------------------
bf <- as.data.table(readRDS(BF))
bf[, `:=`(chr_id = as.integer(sub("Chr", "", Chr)), win = Pos %/% WIN_BP,
          bfmax = pmax(BF1, BF2))]
clim <- bf[, .(clim = max(bfmax), n15 = sum(bfmax >= BF_STRONG), n_bf = .N),
           by = .(chr_id, win)]
d <- merge(pl, clim, by = c("chr_id", "win"))
setorder(d, chr_id, win)

## ---- correlations ----------------------------------------------------------
sp <- function(x, y) suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
r_res  <- sp(d$residual, d$clim)
r_emp  <- sp(d$emp,      d$clim)
r_sim  <- sp(d$sim_med,  d$clim)
r_res15 <- sp(d$residual, d$n15)

## per-chromosome circular-rotation null for residual vs clim (two-sided) ------
chr_idx <- split(seq_len(nrow(d)), d$chr_id)
rot_clim <- function() {
  cr <- d$clim
  for (ix in chr_idx) { k <- length(ix)
    if (k > 1) cr[ix] <- d$clim[ix][ (seq_len(k) - 1 + sample.int(k, 1)) %% k + 1 ] }
  sp(d$residual, cr)
}
null_r <- replicate(N_PERM, rot_clim())
p_res  <- (1 + sum(abs(null_r - mean(null_r)) >= abs(r_res - mean(null_r)))) / (N_PERM + 1)

cat("\n=== ld_w residual (excess over neutral) vs genotype-climate association ===\n")
cat(sprintf("windows: %d  (above 95%% envelope: %d)\n", nrow(d), d[, sum(excess)]))
cat(sprintf("Spearman(residual, clim BF_max) = %+.3f   (rotation p = %.4f)\n", r_res, p_res))
cat(sprintf("Spearman(residual, n[BF>=%d])    = %+.3f\n", BF_STRONG, r_res15))
cat("-- context (same climate track vs raw LD levels) --\n")
cat(sprintf("Spearman(empirical ld_w, clim)  = %+.3f\n", r_emp))
cat(sprintf("Spearman(neutral   ld_w, clim)  = %+.3f\n", r_sim))
cat("-- climate strength by window class --\n")
print(d[, .(n = .N, med_clim = round(median(clim), 2), mean_n15 = round(mean(n15), 2)),
         by = .(class = fifelse(excess, "above 95% envelope",
                        fifelse(residual > 0, "above 1:1 line", "below 1:1 line")))])
cat("\nNOTE: windows are spatially autocorrelated; rotation p respects that. clim is\n",
    "empirical-only (no neutral analogue), so the emp/sim rows show whether climate\n",
    "tracks LD level generally vs the excess-over-neutral residual specifically.\n")

## ---- figure ----------------------------------------------------------------
pA <- ggplot(d, aes(clim, residual)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_vline(xintercept = BF_STRONG, linetype = 3, colour = "grey65") +
  geom_point(aes(colour = excess), size = 1.1, alpha = 0.8, stroke = 0) +
  geom_smooth(method = "loess", se = FALSE, span = 1, colour = "#d95f02", linewidth = 0.7) +
  scale_colour_manual(values = c("FALSE" = "grey45", "TRUE" = "#d73027"),
                      labels = c("within/below envelope", "above 95% envelope"), name = NULL) +
  labs(x = "per-window climate association  max BF (dB)  [dotted = BF 15]",
       y = expression("ld_w residual  (empirical - neutral median)"),
       title = sprintf("excess LD vs climate association  (Spearman %+.2f, rot p=%.3f)", r_res, p_res)) +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "top")

pB <- ggplot(d, aes(fifelse(excess, "above\n95% env", fifelse(residual > 0, "above\n1:1", "below\n1:1")), clim)) +
  geom_boxplot(outlier.size = 0.5, fill = "grey90", width = 0.6) +
  geom_hline(yintercept = BF_STRONG, linetype = 3, colour = "grey55") +
  labs(x = NULL, y = "per-window climate max BF (dB)",
       title = "climate association by LD-excess class") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())

p <- pA + pB + plot_layout(widths = c(2, 1)) +
  plot_annotation(title = "Empirical ld_w excess over neutral vs genotype-climate association (local adaptation)",
                  subtitle = "residual = empirical minus 1000-rep neutral median; climate = BayPass BF on raw environmental PC1/PC2")
ggsave(OUTPDF, p, width = 12, height = 5.2)
ggsave(OUTPNG, p, width = 12, height = 5.2, dpi = 150)
saveRDS(list(d = d, r_res = r_res, p_res = p_res, r_emp = r_emp, r_sim = r_sim,
             r_res15 = r_res15), "module_di25/data/di25_ldw_residual_vs_climate.rds")
cat("saved:", OUTPNG, "\n")
