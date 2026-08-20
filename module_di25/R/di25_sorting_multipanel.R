## =============================================================================
## module_di25 -- combined sorting multipanel (self-contained, titles stripped)
## =============================================================================
## Rebuilds all four panels so the descriptive in-plot titles/subtitles can be
## dropped (they belong in the caption); layout puts (a) next to (c).
##   a  ancestry sorting: observed vs 1000-rep simulated null   (from di25_sorting_emp_vs_sim.rds)
##   b  sorting vs recombination across tau                     (from di25_recomb_tau_sweep.rds)
##   c  rotation-null overlap of sorted SNPs with BDMI regions  (recomputed; cutoff 5, mirrors bdmi_sorting_null_hist.R)
##   d  BDMI regions vs sorting vs low recombination (circos)   (recomputed; cutoff 5, mirrors bdmi_sorting_circos.R)
## a/b are native ggplot; c/d are base-R rendered title-less to cached PNGs
## (module_di25/data/mp_*.png) and placed via cowplot::draw_image. The heavy c/d
## renders are cached -- set RECOMPUTE=1 to force-rebuild them.
##
## Run from the formica_hybrid repo root:
##   Rscript module_di25/R/di25_sorting_multipanel.R
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(cowplot); library(magick) })
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()
source("module_di25/R/diem_circos_core.R")          # render_circos_raster()

FIG <- "module_di25/Figures"; DATA <- "module_di25/data"
OUTPNG <- file.path(FIG, "di25_sorting_multipanel.png"); OUTPDF <- sub("png$", "pdf", OUTPNG)
C_PNG  <- file.path(DATA, "mp_bdmi_hist_notitle.png")
D_PNG  <- file.path(DATA, "mp_bdmi_circos_notitle.png")
RECOMP <- nzchar(Sys.getenv("RECOMPUTE"))
BEDDIR <- "data/liftoff_Frufa_DTOL_PR"

merge_iv <- function(s, e) {
  o <- order(s); s <- s[o]; e <- e[o]; cs <- s[1L]; ce <- e[1L]; oS <- numeric(0); oE <- numeric(0)
  for (i in seq_along(s)[-1L]) if (s[i] <= ce) ce <- max(ce, e[i]) else { oS <- c(oS, cs); oE <- c(oE, ce); cs <- s[i]; ce <- e[i] }
  list(s = c(oS, cs), e = c(oE, ce))
}
in_intervals <- function(q, iv) { if (!length(iv$s)) return(logical(length(q)))
  (findInterval(q, as.vector(rbind(iv$s, iv$e))) %% 2L) == 1L }

## ============================ panel a (ggplot) ==============================
o   <- readRDS(file.path(DATA, "di25_sorting_emp_vs_sim.rds"))
sim <- as.data.table(o$sim); emp <- as.data.table(o$emp)
pA <- ggplot(sim, aes(factor(tau), pct_sorted)) +
  geom_jitter(width = 0.18, height = 0, colour = "#1b9e77", alpha = 0.25, size = 0.5) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = NA, colour = "#0b3d2e") +
  geom_point(data = emp, aes(factor(tau), pct_sorted), colour = "#d95f02", size = 3.2, shape = 18) +
  labs(x = expression("sorting threshold  " * tau), y = "% SNPs sorted") +
  theme_bw(base_size = 13) +
  ## large left margin pushes the y-axis right so the top-left panel label clears it
  theme(panel.grid.minor = element_blank(), plot.margin = margin(t = 10, r = 5, b = 3, l = 34))

## ============================ panel b (ggplot) ==============================
rs   <- readRDS(file.path(DATA, "di25_recomb_tau_sweep.rds"))
dec  <- as.data.table(rs$deciles); sweep <- as.data.table(rs$sweep)
ne   <- as.data.table(rs$n_units_per_decile)[order(rbin)]
de   <- dec[level == "eMLG"]; de[, tau := factor(tau)]
med_lab <- de[, .(m = med_recomb[1]), by = rbin][order(rbin)]$m
ymax <- max(de$frac_sorted); bs <- ymax * 0.98 / max(ne$N); N_MIN <- 100L
lo_bins <- ne[N < N_MIN, rbin]
pB1 <- ggplot(de, aes(rbin)) +
  { if (length(lo_bins)) annotate("rect", xmin = min(lo_bins) - 0.5, xmax = max(lo_bins) + 0.5, ymin = -Inf, ymax = Inf, fill = "grey95") } +
  geom_col(data = ne, aes(rbin, N * bs), fill = "grey86", width = 0.85, inherit.aes = FALSE) +
  geom_line(aes(y = frac_sorted, colour = tau), linewidth = 0.9) +
  geom_point(aes(y = frac_sorted, colour = tau), size = 1.8) +
  scale_colour_viridis_d(end = 0.9, name = expression(tau)) +
  scale_x_continuous(breaks = 1:10, labels = med_lab) +
  scale_y_continuous(name = "fraction sorted (eMLG units)",
                     sec.axis = sec_axis(~ . / bs, name = "n independent units")) +
  labs(x = "recombination decile (median cM/Mb)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.margin = margin(t = 14, r = 5, b = 3, l = 8))
pB2 <- ggplot(sweep, aes(factor(tau), coef_recomb)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_errorbar(aes(ymin = boot_lo, ymax = boot_hi), width = 0.18, colour = "#315B7D") +
  geom_point(size = 2.6, colour = "#315B7D") +
  labs(x = expression(tau), y = "recomb coef (logit / log10 cM/Mb)") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
pB <- plot_grid(pB1, pB2, ncol = 2, rel_widths = c(2, 1))

## ================= panel c: rotation-null histograms (base R) ================
if (RECOMP || !file.exists(C_PNG)) {
  set.seed(1); N_PERM <- 5000L; TAU <- 0.6; CUTOFF_K <- 5L
  ps  <- readRDS(file.path(DATA, "di25_sorting_snp.rds"))
  chr <- sub(":.*", "", ps$marker); pos <- as.integer(sub(".*:", "", ps$marker))
  ok  <- ps$differentiated & ps$n_obs > 0 & is.finite(ps$uni_score)
  cls <- rep("unsorted", nrow(ps))
  cls[ok] <- classify_sort(ps$n_aqu[ok], ps$n_pol[ok], ps$n_obs[ok], sort_th = TAU, sort_rule = "binom", alpha = 0.05)
  snp <- data.table(chr = chr, pos = pos, sorted = cls != "unsorted", pol = cls == "polyctena", aqu = cls == "aquilonia")
  setkey(snp, chr, pos); chr_len <- snp[, .(len = max(pos)), by = chr]
  bedf <- list.files(BEDDIR, sprintf("^bdmi_candidates\\.cutoff_%d_.*\\.bed$", CUTOFF_K))
  x2 <- as.numeric(sub("^(0)(\\d+)$", "0.\\2", sub(".*cutoff_\\d+_(\\d+)\\..*", "\\1", bedf)))
  bed <- fread(file.path(BEDDIR, bedf), header = FALSE, col.names = c("chr", "start", "end"))
  bed[, chr := sub("chromosome_", "Chr", chr)]; bed <- bed[chr %in% chr_len$chr]
  ivs <- lapply(split(bed, bed$chr), function(b) merge_iv(b$start, b$end))
  snp[, inb := FALSE]; for (cc in names(ivs)) { idx <- snp$chr == cc; snp$inb[idx] <- in_intervals(snp$pos[idx], ivs[[cc]]) }
  obs <- c(any = sum(snp$inb & snp$sorted), pol = sum(snp$inb & snp$pol), aqu = sum(snp$inb & snp$aqu)); n_in <- sum(snp$inb)
  null <- matrix(0L, N_PERM, 3, dimnames = list(NULL, c("any", "pol", "aqu")))
  for (p in seq_len(N_PERM)) { s <- 0L; po <- 0L; a <- 0L
    for (cc in names(ivs)) { idx <- which(snp$chr == cc); if (!length(idx)) next
      L <- chr_len[chr == cc, len]; off <- sample.int(L, 1L)
      m <- in_intervals(((snp$pos[idx] - 1L + off) %% L) + 1L, ivs[[cc]])
      s <- s + sum(m & snp$sorted[idx]); po <- po + sum(m & snp$pol[idx]); a <- a + sum(m & snp$aqu[idx]) }
    null[p, ] <- c(s, po, a) }
  emp_p <- sapply(1:3, function(j) (1 + sum(null[, j] >= obs[j])) / (1 + N_PERM)); fold <- obs / colMeans(null)
  labs <- c(any = "any sorted", pol = "toward F. polyctena", aqu = "toward F. aquilonia")
  cols <- c(any = "#315B7D", pol = "#D3C93B", aqu = "#21918C")
  png(C_PNG, width = 2400, height = 820, res = 300, type = "cairo")
  op <- par(mfrow = c(1, 3), mar = c(4.2, 4.0, 2.4, 1.0), mgp = c(2.3, 0.7, 0))
  for (j in 1:3) { v <- null[, j]; ob <- obs[j]
    xr <- range(v, ob); br <- seq(xr[1] - 0.5, xr[2] + 0.5, length.out = 40); h <- hist(v, breaks = br, plot = FALSE)
    plot(h, col = "grey85", border = "grey70", main = labs[j], cex.main = 1.1,
         xlab = "sorted SNPs inside BDMI regions", ylab = "null permutations", xlim = c(min(br), max(xr[2], ob) * 1.02))
    abline(v = ob, col = cols[j], lwd = 2.5); ymx <- par("usr")[4]; gap <- 0.02 * diff(par("usr")[1:2])
    p_str <- if (emp_p[j] <= 1 / (N_PERM + 1)) sprintf("<%.1g", 1 / N_PERM) else sprintf("%.4f", emp_p[j])
    text(ob - gap, ymx * 0.96, sprintf("observed = %d", ob), col = cols[j], adj = c(1, 1), cex = 0.95, font = 2)
    text(ob - gap, ymx * 0.78, sprintf("null mean = %.0f\nfold = %.1fx\np = %s", mean(v), fold[j], p_str), adj = c(1, 1), cex = 0.85, col = "grey20") }
  par(op); dev.off()
  message(sprintf("[c] cutoff 5 (X2=%.3f): obs any/pol/aqu = %d/%d/%d", x2, obs[1], obs[2], obs[3]))
}

## ================= panel d: BDMI x sorting x low-recomb circos (base R) ======
if (RECOMP || !file.exists(D_PNG)) {
  TAU_GRID <- c(0.5, 0.6, 0.7, 0.8); CUTOFF_K <- 5L; LD_A <- 1e-5
  BDMI_COL <- "#E08214"; BDMI_OVL <- "#4575B4"; LOWREC_COL <- "#D01C8B"
  SORT_PAL <- c("#F4F4F4", "#21918C", "#D3C93B", "#440154"); CLS_CODE <- c(unsorted = 0L, aquilonia = 1L, polyctena = 2L, unresolved = 3L, ambiguous = 3L)
  ps <- readRDS(file.path(DATA, "di25_sorting_snp.rds"))
  chr_num <- as.integer(sub("Chr", "", sub(":.*", "", ps$marker))); pos <- as.integer(sub(".*:", "", ps$marker))
  ord <- order(chr_num, pos); ps <- ps[ord]; chr_num <- chr_num[ord]; pos <- pos[ord]
  ok <- ps$differentiated & ps$n_obs > 0 & is.finite(ps$uni_score)
  code <- vapply(TAU_GRID, function(tau) { v <- integer(nrow(ps))
    cl <- classify_sort(ps$n_aqu[ok], ps$n_pol[ok], ps$n_obs[ok], sort_th = tau, sort_rule = "binom", alpha = 0.05)
    v[which(ok)] <- CLS_CODE[cl]; v }, integer(nrow(ps)))
  bedf <- list.files(BEDDIR, sprintf("^bdmi_candidates\\.cutoff_%d_.*\\.bed$", CUTOFF_K))
  bed <- fread(file.path(BEDDIR, bedf), header = FALSE, col.names = c("chr", "start", "end")); bed[, chrn := as.integer(sub("chromosome_", "", chr))]
  inb <- logical(nrow(ps))
  for (cc in unique(bed$chrn)) { iv <- merge_iv(bed[chrn == cc]$start, bed[chrn == cc]$end); idx <- chr_num == cc; inb[idx] <- in_intervals(pos[idx], iv) }
  sorted06 <- code[, 2] != 0L; bdmi_state <- integer(nrow(ps)); bdmi_state[inb & !sorted06] <- 1L; bdmi_state[inb & sorted06] <- 2L
  ld <- readRDS(file.path(DATA, "di25_ld_decay.rds")); lowrec <- logical(nrow(ps))
  for (cc in unique(chr_num)) { dec2 <- ld$by_chr[[paste0("Chr", cc)]]$decay; if (is.null(dec2) || !nrow(dec2)) next
    low <- dec2[is.finite(a) & a < LD_A]; if (!nrow(low)) next; iv <- merge_iv(low$start, low$end); idx <- chr_num == cc; lowrec[idx] <- in_intervals(pos[idx], iv) }
  png(D_PNG, width = 2600, height = 2600, res = 300, type = "cairo")
  op <- par(mar = c(0, 0, 0, 0), xpd = NA); plot.new(); plot.window(xlim = c(-1, 1), ylim = c(-1, 1), asp = 1)
  render_circos_raster(code, chr_num, palette = SORT_PAL, title = "", r_in = 0.30, r_out = 0.82,
    ring_labels = sprintf("tau=%.1f", TAU_GRID), draw_chr_labels = FALSE, new_device = FALSE, add = TRUE, bg_col = "#FFFFFF")
  render_circos_raster(matrix(bdmi_state, ncol = 1L), chr_num, palette = c(NA, BDMI_COL, BDMI_OVL), r_in = 0.835, r_out = 0.87, ring_sep = FALSE, draw_chr_labels = FALSE, add = TRUE, bg_col = NA)
  render_circos_raster(matrix(as.integer(lowrec), ncol = 1L), chr_num, palette = c(NA, LOWREC_COL), r_in = 0.885, r_out = 0.92, ring_sep = FALSE, draw_chr_labels = TRUE, chr_label_cex = 0.5, add = TRUE, bg_col = NA)
  legend("bottomleft", bty = "n", cex = 0.8, border = NA, inset = c(0.0, 0.0),
    legend = c("toward F. aquilonia", "toward F. polyctena", "direction unresolved", "unsorted",
               "BDMI region (no overlap)", "BDMI region overlapping sorted", "low recombination (LD-decay ~ 0)"),
    fill = c(SORT_PAL[2], SORT_PAL[3], SORT_PAL[4], SORT_PAL[1], BDMI_COL, BDMI_OVL, LOWREC_COL))
  par(op); dev.off(); message("[d] wrote circos")
}

## ================================ compose ===================================
## drop images slightly from the top so the bold panel labels sit in a clear strip
pC <- ggdraw() + draw_image(C_PNG, y = 0, height = 0.92)
pD <- ggdraw() + draw_image(D_PNG, y = 0, height = 0.96)
## circos = panel a (left, large); sorting = b, histograms = c (top-right row,
## kept short so b is ~as tall as the c histograms); recomb = d fills below.
row_bc <- plot_grid(pA, pC, ncol = 2, rel_widths = c(1, 2), labels = c("b", "c"), label_size = 24, label_fontface = "bold")
right  <- plot_grid(row_bc, pB, ncol = 1, rel_heights = c(1, 1.8), labels = c("", "d"), label_size = 24, label_fontface = "bold")
full   <- plot_grid(pD, right, ncol = 2, rel_widths = c(1, 1.3), labels = c("a", ""), label_size = 24, label_fontface = "bold")
full   <- full + theme(plot.margin = margin(12, 6, 4, 6))   # outer headroom so top-edge labels aren't clipped
ggsave(OUTPNG, full, width = 17, height = 6.8, dpi = 200, bg = "white")
ggsave(OUTPDF, full, width = 17, height = 6.8, device = cairo_pdf, bg = "white")
cat("saved:", OUTPNG, "\n       ", OUTPDF, "\n")
