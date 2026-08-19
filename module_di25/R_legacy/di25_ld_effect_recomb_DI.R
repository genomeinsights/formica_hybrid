## =========================================================
## module_di25 (high-DI analyses) -- the effect of LD-clustering on sorting results
## =========================================================
## Amount of sorting toward each species vs recombination decile, split by DI bin,
## for SNP-based vs cluster-based (eMLG) results. Shows the LD-pseudoreplication
## effect directly: at the SNP level low-recombination regions inflate the sorted
## counts (direction set by DI); collapsing to independent units removes it.
## phi = 0.85, tau = 0.6. Additive: reads di25_sorting_{snp,emlg}.rds only.
##
## Run from the repo root:  Rscript module_di25/R/di25_ld_effect_recomb_DI.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2) })
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()

TAU <- as.numeric(Sys.getenv("TAU", "0.6"))          # env-overridable sorting threshold
STAMP <- sprintf("tau%02d", round(TAU * 10))
COL_AQU <- "#21918C"; COL_POL <- "#D3C93B"
OUTDIR <- "module_di25/data"; FIGDIR <- "module_di25/Figures"

## ---- recombination interpolation ---------------------------------------
rec <- fread("data/Frufa_DTOL_PR.ref_genome.recmap"); setnames(rec, 1:4, c("chr","p","cM","cMMb"))
rec[, Chr := sub("chromosome_", "Chr", chr)]
rc <- function(ch, pos) { o <- rep(NA_real_, length(pos))
  for (cc in unique(ch)) { r <- rec[Chr == paste0("Chr", cc)]; if (nrow(r) < 2) next
    i <- which(ch == cc); o[i] <- approx(r$p, r$cMMb, xout = pos[i], rule = 2)$y }; o }

## ---- load + annotate (recomb, DI, sort direction at tau) ----------------
prep <- function(dt, is_snp) {
  if (is_snp) { dt[, `:=`(chr = as.integer(sub("Chr","",sub(":.*","",marker))), pos = as.integer(sub(".*:","",marker)))] }
  else {
    g <- readRDS(file.path(OUTDIR, "di25_clustering_cM5.rds"))$groups
    rm <- g$representative[match(dt$group_id, g$group_id)]
    dt[, `:=`(chr = as.integer(sub("Chr","",sub(":.*","",rm))), pos = as.integer(sub(".*:","",rm)))]
  }
  dt[, recomb := rc(chr, pos)]
  dt <- dt[differentiated == TRUE & n_obs > 0 & is.finite(recomb) & is.finite(DI)]
  cl <- classify_sort(dt$n_aqu, dt$n_pol, dt$n_obs, sort_th = TAU, sort_rule = "binom", alpha = 0.05)
  dt[, dir := cl]
  dt[, level := if (is_snp) "SNP-based" else "cluster-based (eMLG)"]
  dt[]
}
s <- prep(as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_snp.rds"))),  TRUE)
e <- prep(as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_emlg.rds"))), FALSE)

## ---- shared bins: recomb deciles (from SNPs) + DI terciles --------------
rbrk  <- quantile(s$recomb, 0:10/10, na.rm = TRUE)
dibrk <- quantile(s$DI, c(0, 1/3, 2/3, 1), na.rm = TRUE); dibrk[c(1, length(dibrk))] <- c(-Inf, Inf)
di_lab <- c(sprintf("low DI  (%.0f to %.0f)\nmore polyctena-like", dibrk[1] <- min(s$DI), quantile(s$DI,1/3)),
            sprintf("mid DI  (%.0f to %.0f)", quantile(s$DI,1/3), quantile(s$DI,2/3)),
            sprintf("high DI  (%.0f to %.0f)\nmore aquilonia-like", quantile(s$DI,2/3), max(s$DI)))
for (d in list(s, e)) {
  d[, rdec := cut(recomb, rbrk, include.lowest = TRUE, labels = FALSE)]
  d[, dibin := cut(DI, dibrk, labels = di_lab)]
}

## ---- fraction sorted toward each species per (level, DI bin, decile) ----
cols <- c("level", "dibin", "rdec", "recomb", "dir")
agg <- rbind(s[, ..cols], e[, ..cols])[!is.na(rdec) & !is.na(dibin),
  .(med_recomb = median(recomb), n = .N,
    aquilonia = mean(dir == "aquilonia"), polyctena = mean(dir == "polyctena")),
  by = .(level, dibin, rdec)]
long <- melt(agg, id.vars = c("level","dibin","rdec","med_recomb","n"),
             measure.vars = c("aquilonia","polyctena"), variable.name = "species", value.name = "frac")
long[, level := factor(level, levels = c("SNP-based", "cluster-based (eMLG)"))]
long[, dibin := factor(dibin, levels = di_lab)]
saveRDS(agg, file.path(OUTDIR, sprintf("di25_ld_effect_recomb_DI_%s.rds", STAMP)))

## ---- figure ------------------------------------------------------------
xlab_med <- s[!is.na(rdec), round(median(recomb),1), by = rdec][order(rdec)]$V1
MIN_N <- 15L                                          # lines only through adequately-sampled cells
ymax <- max(long$frac, na.rm = TRUE)
bs   <- ymax / max(agg$n)                             # scale n-bars to the fraction axis (shared across panels)
p <- ggplot(long, aes(rdec)) +
  geom_col(data = agg, aes(rdec, n * bs), fill = "grey87", width = 0.85, inherit.aes = FALSE) +
  geom_line(data = long[n >= MIN_N], aes(y = frac, colour = species), linewidth = 0.8) +
  geom_point(aes(y = frac, colour = species), size = 1.2) +
  facet_grid(dibin ~ level) +
  scale_colour_manual(values = c(aquilonia = COL_AQU, polyctena = COL_POL),
                      labels = c("toward F. aquilonia", "toward F. polyctena"), name = NULL) +
  scale_x_continuous(breaks = 1:10, labels = xlab_med) +
  scale_y_continuous(name = "fraction of differentiated loci sorted toward the species",
                     sec.axis = sec_axis(~ . / bs, name = "grey bars: n loci per cell  (SNPs left / clusters right)")) +
  labs(x = "recombination decile  (median cM/Mb)",
       title = sprintf("Effect of LD-clustering on ancestry sorting (phi = 0.85, tau = %.1f)", TAU),
       subtitle = "SNP sorting inflates at low recombination (left); grey bars show the SNP counts collapse to almost no clusters there (right).\nRows = DI bins (direction gradient); lines drawn only where a cell has >=15 loci.") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text.y = element_text(size = 8),
        plot.subtitle = element_text(size = 8.5), axis.text.x = element_text(size = 7))
ggsave(file.path(FIGDIR, sprintf("di25_ld_effect_recomb_DI_%s.png", STAMP)), p, width = 10, height = 7.5, dpi = 200)
cat("[ld-effect] wrote", file.path(FIGDIR, sprintf("di25_ld_effect_recomb_DI_%s.png", STAMP)), "\n")
print(agg[, .(n_loci = sum(n)), by = .(level, dibin)])
