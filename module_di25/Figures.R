## =========================================================
## module_di25 (high-DI analyses) -- MANUSCRIPT figures
## =========================================================
## Home for publication-ready figures. Unlike the exploratory scripts in R/,
## these carry NO descriptive titles/subtitles -- only axes, legend, panel tags
## and minimal condition labels; the wording lives in the manuscript caption.
## Each figure is a self-contained fig_*() section that writes a PDF (+PNG) to
## Figures/. Reads the cached di25_* outputs; runs nothing heavy.
##
## Run from the repo root:  Rscript module_di25/Figures.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()

## ---- shared palette / paths --------------------------------------------
COL_AQU <- "#21918C"; COL_POL <- "#D3C93B"           # F. aquilonia / F. polyctena
OUTDIR <- "module_di25/data"; FIGDIR <- "module_di25/Figures"
theme_ms <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.tag = element_text(face = "bold"),
        legend.key.width = unit(1.4, "cm"))          # wide keys so line types read

save_fig <- function(fig, name, width, height) {
  ggsave(file.path(FIGDIR, paste0(name, ".pdf")), fig, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(FIGDIR, paste0(name, ".png")), fig, width = width, height = height, dpi = 300,
         device = "png", type = "cairo")
  cat("[Figures.R] wrote", name, ".pdf/.png\n")
}

## =========================================================================
## FIG: directional ancestry sorting vs recombination, tau = 0.6 | 0.8
## =========================================================================
fig_recomb_directional <- function() {
  rec <- fread("data/Frufa_DTOL_PR.ref_genome.recmap"); setnames(rec, 1:4, c("chr","p","cM","cMMb"))
  rec[, Chr := sub("chromosome_", "Chr", chr)]
  rc <- function(ch, pos) { o <- rep(NA_real_, length(pos))
    for (cc in unique(ch)) { r <- rec[Chr == paste0("Chr", cc)]; if (nrow(r) < 2) next
      i <- which(ch == cc); o[i] <- approx(r$p, r$cMMb, xout = pos[i], rule = 2)$y }; o }
  wilson <- function(k, n) { z <- 1.959964; p <- k/n; d <- 1 + z^2/n
    ctr <- (p + z^2/(2*n))/d; hw <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/d
    list(lo = pmax(0, ctr-hw), hi = pmin(1, ctr+hw)) }
  prep <- function(dt, is_snp) {
    if (is_snp) dt[, `:=`(chr = as.integer(sub("Chr","",sub(":.*","",marker))), pos = as.integer(sub(".*:","",marker)))]
    else { g <- readRDS(file.path(OUTDIR, "di25_clustering_cM5.rds"))$groups
      rmk <- g$representative[match(dt$group_id, g$group_id)]
      dt[, `:=`(chr = as.integer(sub("Chr","",sub(":.*","",rmk))), pos = as.integer(sub(".*:","",rmk)))] }
    dt[, recomb := rc(chr, pos)]
    dt <- dt[differentiated == TRUE & n_obs > 0 & is.finite(recomb)]
    dt[, level := if (is_snp) "SNP" else "eMLG"][]
  }
  s <- prep(as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_snp.rds"))),  TRUE)
  e <- prep(as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_emlg.rds"))), FALSE)
  rbrk <- quantile(s$recomb, 0:10/10, na.rm = TRUE)
  s[, rdec := cut(recomb, rbrk, include.lowest = TRUE, labels = FALSE)]
  e[, rdec := cut(recomb, rbrk, include.lowest = TRUE, labels = FALSE)]

  dir_long <- function(tau) {
    cl <- function(dt) classify_sort(dt$n_aqu, dt$n_pol, dt$n_obs, sort_th = tau, sort_rule = "binom", alpha = 0.05)
    a <- rbind(cbind(s[, .(level, rdec)], dir = cl(s)), cbind(e[, .(level, rdec)], dir = cl(e)))
    ag <- a[!is.na(rdec), .(n = .N, n_aqu = sum(dir == "aquilonia"), n_pol = sum(dir == "polyctena")),
            by = .(level, rdec)]
    lg <- melt(ag, id.vars = c("level","rdec","n"), measure.vars = c("n_aqu","n_pol"),
               variable.name = "species", value.name = "k")
    lg[, `:=`(frac = k/n, species = fifelse(species == "n_aqu", "aquilonia", "polyctena"))]
    lg[level == "eMLG", `:=`(lo = wilson(k, n)$lo, hi = wilson(k, n)$hi)][]
  }
  d6 <- dir_long(0.6); d8 <- dir_long(0.8)

  ne      <- unique(e[!is.na(rdec), .(n = .N), by = rdec][order(rdec)])   # tau-independent unit counts
  lo_bins <- ne[n < 100L, rdec]
  xlab    <- s[!is.na(rdec), round(median(recomb), 1), by = rdec][order(rdec)]$V1
  ymax <- max(c(d6$hi, d6$frac, d8$hi, d8$frac), na.rm = TRUE)
  bs   <- ymax * 0.98 / max(ne$n)

  panel <- function(long, tau, secondary) {
    p <- ggplot(long, aes(rdec, frac, colour = species, linetype = level)) +
      { if (length(lo_bins)) annotate("rect", xmin = min(lo_bins)-0.5, xmax = max(lo_bins)+0.5,
                                      ymin = -Inf, ymax = Inf, fill = "grey95") } +
      geom_col(data = ne, aes(rdec, n * bs), fill = "grey86", width = 0.85, inherit.aes = FALSE) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, na.rm = TRUE, linetype = 1, show.legend = FALSE) +
      geom_line(linewidth = 0.9) + geom_point(size = 1.7) +
      annotate("text", x = 5.5, y = ymax * 0.99, size = 3.5, parse = TRUE,
               label = sprintf("phantom() >= %d*'%% of populations near-fixed'", round(tau * 100))) +
      scale_colour_manual(values = c(aquilonia = COL_AQU, polyctena = COL_POL),
                          breaks = c("aquilonia", "polyctena"),
                          labels = c(expression(italic("F. aquilonia")), expression(italic("F. polyctena"))),
                          name = "sorted toward") +
      scale_linetype_manual(values = c(eMLG = 1, SNP = 2), breaks = c("eMLG", "SNP"), name = "level") +
      guides(linetype = guide_legend(override.aes = list(linewidth = 0.7)),
             colour   = guide_legend(override.aes = list(linetype = 0, size = 3))) +
      scale_x_continuous(breaks = 1:10, labels = xlab) +
      coord_cartesian(ylim = c(0, ymax)) +
      labs(x = "recombination decile (median cM/Mb)") + theme_ms
    if (secondary)
      p + scale_y_continuous(name = NULL, sec.axis = sec_axis(~ . / bs, name = "n independent eMLG units")) +
          theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank())
    else
      p + scale_y_continuous(name = "fraction of differentiated loci\nsorted toward the species")
  }

  ## panel c: directional recombination coefficient per tau, with a genomic
  ## block bootstrap over chromosomes (units are only ~independent). One logistic
  ## per (tau, species): P(sorted toward species) ~ log10(recomb). Negative =
  ## that direction's sorting is concentrated at low recombination.
  lr <- log10(e$recomb + 0.1); chr_idx <- split(seq_len(nrow(e)), e$chr); chrs <- names(chr_idx)
  fit_slope <- function(y, idx) {
    if (sum(y[idx]) < 3L) return(NA_real_)
    suppressWarnings(tryCatch(glm.fit(cbind(1, lr[idx]), y[idx], family = binomial())$coefficients[2],
                              error = function(e) NA_real_))
  }
  set.seed(1)
  coefd <- rbindlist(lapply(c(0.5, 0.6, 0.7, 0.8), function(tau) {
    cl <- classify_sort(e$n_aqu, e$n_pol, e$n_obs, sort_th = tau, sort_rule = "binom", alpha = 0.05)
    rbindlist(lapply(c("aquilonia", "polyctena"), function(sp) {
      y <- as.integer(cl == sp)
      bc <- vapply(1:400, function(b)
        fit_slope(y, unlist(chr_idx[sample(chrs, length(chrs), TRUE)], use.names = FALSE)), numeric(1))
      data.table(tau = tau, species = sp, est = fit_slope(y, seq_along(y)),
                 lo = quantile(bc, .025, na.rm = TRUE), hi = quantile(bc, .975, na.rm = TRUE))
    }))
  }))
  cat("[Figures.R] directional recomb coefficients:\n"); print(coefd[, .(tau, species, est = round(est,2), lo = round(lo,2), hi = round(hi,2))])

  pd <- position_dodge(0.45)
  pc <- ggplot(coefd, aes(factor(tau), est, colour = species)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, position = pd, show.legend = FALSE) +
    geom_point(size = 2.4, position = pd) +
    scale_colour_manual(values = c(aquilonia = COL_AQU, polyctena = COL_POL),
                        breaks = c("aquilonia", "polyctena"),
                        labels = c(expression(italic("F. aquilonia")), expression(italic("F. polyctena"))),
                        name = "sorted toward") +
    guides(colour = "none") +                       # same species colours as a/b -> reuse that legend
    labs(x = "sorting threshold (fraction of populations near-fixed)",
         y = "recombination coefficient\n(logit per log10 cM/Mb)") + theme_ms

  ## tag the whole top row as "a" (two sub-panels), the coefficient panel as "b"
  fig <- ((panel(d6, 0.6, FALSE) + labs(tag = "a")) | panel(d8, 0.8, TRUE)) / (pc + labs(tag = "b")) +
    plot_layout(guides = "collect", heights = c(1.35, 1)) & theme(legend.position = "bottom")
  save_fig(fig, "fig_recomb_directional", width = 11, height = 8)
}

fig_recomb_directional()
