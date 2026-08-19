## =========================================================
## module_di25 (high-DI analyses) -- sorting vs recombination, split by DIRECTION
## =========================================================
## The recombination-effect figure (di25_architecture.R style) with the sorted
## fraction split by direction: line type = toward F. aquilonia (solid) vs toward
## F. polyctena (dashed); colour = level (eMLG vs SNP). Grey bars = n independent
## eMLG units per decile; Wilson CIs on the eMLG lines; low-n deciles shaded.
## phi = 0.85, tau = 0.6. Additive: reads di25_sorting_{snp,emlg}.rds only.
##
## Run from the repo root:  Rscript module_di25/R/di25_recomb_directional.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2) })
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()

TAU <- as.numeric(Sys.getenv("TAU", "0.6")); N_MIN <- 100L    # env-overridable threshold
STAMP <- sprintf("tau%02d", round(TAU * 10))
OUTDIR <- "module_di25/data"; FIGDIR <- "module_di25/Figures"
COL_AQU <- "#21918C"; COL_POL <- "#D3C93B"          # species colours

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
  dt[, dir := classify_sort(n_aqu, n_pol, n_obs, sort_th = TAU, sort_rule = "binom", alpha = 0.05)]
  dt[, level := if (is_snp) "SNP" else "eMLG"][]
}
s <- prep(as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_snp.rds"))),  TRUE)
e <- prep(as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_emlg.rds"))), FALSE)
rbrk <- quantile(s$recomb, 0:10/10, na.rm = TRUE)
s[, rdec := cut(recomb, rbrk, include.lowest = TRUE, labels = FALSE)]
e[, rdec := cut(recomb, rbrk, include.lowest = TRUE, labels = FALSE)]

agg <- rbind(s[, .(level, rdec, recomb, dir)], e[, .(level, rdec, recomb, dir)])[!is.na(rdec),
  .(med = median(recomb), n = .N, n_aqu = sum(dir == "aquilonia"), n_pol = sum(dir == "polyctena")),
  by = .(level, rdec)][order(level, rdec)]
long <- melt(agg, id.vars = c("level","rdec","med","n"), measure.vars = c("n_aqu","n_pol"),
             variable.name = "species", value.name = "k")
long[, `:=`(frac = k/n, species = fifelse(species == "n_aqu", "aquilonia", "polyctena"))]
long[level == "eMLG", `:=`(lo = wilson(k, n)$lo, hi = wilson(k, n)$hi)]
ne <- agg[level == "eMLG"]; lo_bins <- ne[n < N_MIN, rdec]
saveRDS(list(agg = agg, long = long), file.path(OUTDIR, sprintf("di25_recomb_directional_%s.rds", STAMP)))

ymax <- max(long$hi, long$frac, na.rm = TRUE); bs <- ymax*0.98 / max(ne$n)
xlab <- s[!is.na(rdec), round(median(recomb),1), by = rdec][order(rdec)]$V1
p <- ggplot(long, aes(rdec, frac, colour = species, linetype = level)) +
  { if (length(lo_bins)) annotate("rect", xmin = min(lo_bins)-0.5, xmax = max(lo_bins)+0.5, ymin = -Inf, ymax = Inf, fill = "grey95") } +
  geom_col(data = ne, aes(rdec, n*bs), fill = "grey86", width = 0.85, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, na.rm = TRUE, linetype = 1, show.legend = FALSE) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.7) +
  scale_colour_manual(values = c(aquilonia = COL_AQU, polyctena = COL_POL),
                      labels = c("F. aquilonia", "F. polyctena"), name = "sorted toward") +
  scale_linetype_manual(values = c(eMLG = 1, SNP = 2), name = "level") +
  scale_x_continuous(breaks = 1:10, labels = xlab) +
  scale_y_continuous(name = "fraction of differentiated loci sorted toward the species",
                     sec.axis = sec_axis(~ ./bs, name = "n independent eMLG units (bars)")) +
  labs(x = "recombination decile  (median cM/Mb)",
       title = sprintf("Ancestry sorting vs recombination, by direction (high-DI loci, tau = %.1f)", TAU),
       subtitle = "colour = species sorted toward; line type = level (solid eMLG / dashed SNP). Bars = n units; low-n deciles shaded; CIs on eMLG.") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 8.5))
ggsave(file.path(FIGDIR, sprintf("di25_recomb_directional_%s.png", STAMP)), p, width = 9.5, height = 5.6, dpi = 200)
cat("[recomb-directional] wrote", file.path(FIGDIR, sprintf("di25_recomb_directional_%s.png", STAMP)), "\n")
