## =========================================================================
## module_di25 -- TALK figures (Figures_Jonna)
## =========================================================================
## Standalone single-panel figures for slides, each rendered on its own so it
## can be placed/arranged independently. tau = 0.8, phi = 0.85 throughout.
##   * fig2a_recomb_tau08      -- fraction sorted toward each species by
##       recombination decile (eMLG solid + SNP dashed; grey bars = n eMLG units).
##   * fig3a_clustersize_tau08 -- fraction of units sorted by LD-cluster size
##       class (grey bars = n eMLG units).
## Both mirror panel a of the corresponding main-text figure in Figures.R
## (fig_recomb_directional / fig_size_directional), lifted out standalone with a
## title and both y-axes shown. Saved to module_di25/Figures_Jonna/.
##
## Run from the repo root:  Rscript module_di25/Figures_Jonna.R
## =========================================================================

suppressMessages({ library(data.table); library(ggplot2);library(patchwork) })
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()
COL_AQU <- "#21918C"; COL_POL <- "#D3C93B"
OUTDIR <- "module_di25/data"; FIGDIR <- "module_di25/Figures_Jonna"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
theme_talk <- theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))
wilson <- function(k, n) { z <- 1.959964; p <- k/n; d <- 1 + z^2/n
  ctr <- (p + z^2/(2*n))/d; hw <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/d; list(lo = pmax(0, ctr-hw), hi = pmin(1, ctr+hw)) }
save2 <- function(p, name, w, h) {
  ggsave(file.path(FIGDIR, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIGDIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, device = "png", type = "cairo")
  cat("[Figures_Jonna] wrote", name, "\n") }

## =========================================================================
## Fig 2a (tau = 0.8): ancestry sorting vs recombination decile
## =========================================================================
rec <- fread("data/Frufa_DTOL_PR.ref_genome.recmap"); setnames(rec, 1:4, c("chr","p","cM","cMMb"))
rec[, Chr := sub("chromosome_", "Chr", chr)]
rc <- function(ch, pos) { o <- rep(NA_real_, length(pos))
  for (cc in unique(ch)) { r <- rec[Chr == paste0("Chr", cc)]; if (nrow(r) < 2) next
    i <- which(ch == cc); o[i] <- approx(r$p, r$cMMb, xout = pos[i], rule = 2)$y }; o }
prep <- function(dt, is_snp) {
  if (is_snp) dt[, `:=`(chr = as.integer(sub("Chr","",sub(":.*","",marker))), pos = as.integer(sub(".*:","",marker)))]
  else { g <- readRDS(file.path(OUTDIR, "di25_clustering_cM5.rds"))$groups
    rmk <- g$representative[match(dt$group_id, g$group_id)]
    dt[, `:=`(chr = as.integer(sub("Chr","",sub(":.*","",rmk))), pos = as.integer(sub(".*:","",rmk)))] }
  dt[, recomb := rc(chr, pos)]; dt <- dt[differentiated == TRUE & n_obs > 0 & is.finite(recomb)]
  dt[, level := if (is_snp) "SNP" else "eMLG"][] }
s <- prep(as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_snp.rds"))),  TRUE)
e <- prep(as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_emlg.rds"))), FALSE)
rbrk <- quantile(s$recomb, 0:10/10, na.rm = TRUE)
s[, rdec := cut(recomb, rbrk, include.lowest = TRUE, labels = FALSE)]
e[, rdec := cut(recomb, rbrk, include.lowest = TRUE, labels = FALSE)]
clr <- function(dt) classify_sort(dt$n_aqu, dt$n_pol, dt$n_obs, sort_th = 0.8, sort_rule = "binom", alpha = 0.05)
a  <- rbind(cbind(s[, .(level, rdec)], dir = clr(s)), cbind(e[, .(level, rdec)], dir = clr(e)))
ag <- a[!is.na(rdec), .(n = .N, n_aqu = sum(dir == "aquilonia"), n_pol = sum(dir == "polyctena")), by = .(level, rdec)]
d8 <- melt(ag, id.vars = c("level","rdec","n"), measure.vars = c("n_aqu","n_pol"), variable.name = "species", value.name = "k")
d8[, `:=`(frac = k/n, species = fifelse(species == "n_aqu", "aquilonia", "polyctena"))]
d8[level == "eMLG", `:=`(lo = wilson(k, n)$lo, hi = wilson(k, n)$hi)]
ne <- unique(e[!is.na(rdec), .(n = .N), by = rdec][order(rdec)]); lo_bins <- ne[n < 100L, rdec]
xlab <- s[!is.na(rdec), round(median(recomb), 1), by = rdec][order(rdec)]$V1
ymax <- max(d8$hi, d8$frac, na.rm = TRUE); bs <- ymax * 0.98 / max(ne$n)
p2a <- ggplot(d8, aes(rdec, frac, colour = species, linetype = level)) +
  { if (length(lo_bins)) annotate("rect", xmin = min(lo_bins)-0.5, xmax = max(lo_bins)+0.5, ymin = -Inf, ymax = Inf, fill = "grey95") } +
  geom_col(data = ne, aes(rdec, n * bs), fill = "grey86", width = 0.85, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, na.rm = TRUE, linetype = 1, show.legend = FALSE) +
  geom_line(linewidth = 1) + geom_point(size = 2.2) +
  scale_colour_manual(values = c(aquilonia = COL_AQU, polyctena = COL_POL), breaks = c("aquilonia", "polyctena"),
                      labels = c(expression(italic("F. aquilonia")), expression(italic("F. polyctena"))), name = "sorted toward") +
  scale_linetype_manual(values = c(eMLG = 1, SNP = 2), breaks = c("eMLG", "SNP"), name = "level") +
  guides(linetype = guide_legend(override.aes = list(linewidth = 0.7)), colour = guide_legend(override.aes = list(linetype = 0, size = 3))) +
  scale_x_continuous(breaks = 1:10, labels = xlab) +
  scale_y_continuous(name = "fraction of units sorted\ntoward the species", labels = scales::percent,
                     sec.axis = sec_axis(~ . / bs, name = "n eMLG units")) +
  coord_cartesian(ylim = c(0, ymax)) +
  labs(x = "recombination decile (median cM/Mb)",
       title = expression(paste("Ancestry sorting vs recombination  (", tau, " = 0.8)"))) + theme_talk
save2(p2a, "fig2a_recomb_tau08", 7.5, 6)


## =========================================================================
## Fig 3a (tau = 0.8): ancestry sorting vs LD-cluster size class
## =========================================================================
u <- as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_emlg.rds")))
u <- u[differentiated == TRUE & n_obs > 0 & !is.na(n_loci)]
SLV <- c("1", "2", "3-4", "5-9", "10-49", "50-199", "200+")
u[, sclass := cut(n_loci, c(0, 1, 2, 4, 9, 49, 199, Inf), labels = SLV)]
cl8 <- classify_sort(u$n_aqu, u$n_pol, u$n_obs, sort_th = 0.8, sort_rule = "binom", alpha = 0.05)
d <- data.table(sclass = as.character(u$sclass), dir = fifelse(cl8 == "aquilonia", "aquilonia", fifelse(cl8 == "polyctena", "polyctena", "uns")))
tot <- d[, .(n = .N), by = sclass]; agg <- d[dir != "uns", .N, by = .(sclass, dir)]
full <- CJ(sclass = SLV, dir = c("aquilonia", "polyctena")); agg <- agg[full, on = .(sclass, dir)][is.na(N), N := 0L]
fa <- tot[agg, on = "sclass"][, `:=`(frac = N/n, lo = wilson(N, n)$lo, hi = wilson(N, n)$hi)]
fa[, `:=`(sclass = factor(sclass, levels = SLV), dir = factor(dir, levels = c("aquilonia", "polyctena")))]
bar <- tot[, .(sclass = factor(sclass, levels = SLV), n)]
ymax3 <- max(fa$hi, fa$frac, na.rm = TRUE); bs3 <- ymax3 * 0.98 / max(bar$n)
p3a <- ggplot(fa, aes(sclass, frac, colour = dir, group = dir)) +
  geom_col(data = bar, aes(sclass, n * bs3), fill = "grey86", width = 0.85, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, show.legend = FALSE) +
  geom_line(linewidth = 1) + geom_point(size = 2.2) +
  scale_colour_manual(values = c(aquilonia = COL_AQU, polyctena = COL_POL), breaks = c("aquilonia", "polyctena"),
                      labels = c(expression(italic("F. aquilonia")), expression(italic("F. polyctena"))), name = "sorted toward") +
  scale_y_continuous(name = "fraction of units sorted\ntoward the species", labels = scales::percent,
                     sec.axis = sec_axis(~ . / bs3, name = "n eMLG units")) +
  coord_cartesian(ylim = c(0, ymax3)) +
  labs(x = "LD-cluster size (n markers)",
       title = expression(paste("Ancestry sorting vs LD-cluster size  (", tau, " = 0.8)"))) +
  theme_talk + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save2(p3a, "fig3a_clustersize_tau08", 7.5, 6)
save2(p2a | p3a, "fig23a_tau08", 15, 6)

## =========================================================================
## keep exact copies of the separately-generated talk circos figures here too,
## so all talk figures live in one folder. They are produced by their own
## scripts (rerun those first to refresh): diem_circos_compare.R,
## di25_population_fixation.R, diem_circos_sorting_sweep_slide.R.
## =========================================================================
for (f in c("diem_circos_compare_snp_vs_ldreduced.png",
            "di25_popfix_snp.png",
            "diem_circos_sorting_sweep_slide.png")) {
  src <- file.path("module_di25/Figures", f)
  if (file.exists(src)) { file.copy(src, file.path(FIGDIR, f), overwrite = TRUE)
    cat("[Figures_Jonna] copied", f, "\n")
  } else cat("[Figures_Jonna] (source missing, skipped)", f, "\n")
}

