## =========================================================
## module_di25 -- ternary view of the recombination / sorting / BDMI interplay
## =========================================================
## A ternary plot needs a 3-part composition summing to 1. Recombination, sorting
## and BDMI are not naturally compositional, so we build one: bin diagnostic SNPs
## by recombination, and within each bin take the loci that are ACTIVE (sorted or
## BDMI-candidate) and split them into three exhaustive classes ---
##   sorted-only : both : BDMI-only  --- which sum to 1.
## Each point is one recombination bin, placed in the triangle by that mix and
## COLOURED by its recombination rate; the points are joined in recombination
## order so the path from low- to high-recombination bins is visible.
##
##   corner (bottom-left)  = sorted only   (sorted, not a BDMI candidate)
##   corner (top)          = both          (sorted AND BDMI candidate)
##   corner (bottom-right) = BDMI only     (BDMI candidate, not sorted)
##
## Reads sorting toward LOW recombination and BDMIs toward HIGH recombination as a
## drift from the sorted-only corner to the BDMI-only corner along the gradient.
##
## Run from repo root:  Rscript module_di25/R/bdmi_sorting_ternary.R
## =========================================================

suppressMessages(library(data.table))
source("moduleA_sorting/R/parallelism_stats.R")

TAU    <- 0.6
A_CUT  <- as.integer(Sys.getenv("CUTOFF", "13"))   # cutoff-k index (5 -> X2=0.05, 13 -> X2=0.03)
NBIN   <- 15L
BEDDIR <- "data/liftoff_Frufa_DTOL_PR"
OUTPNG <- sprintf("module_di25/Figures/bdmi_sorting_ternary_cutoff%d.png", A_CUT)

rec <- fread("data/Frufa_DTOL_PR.ref_genome.recmap"); setnames(rec, 1:4, c("chr","p","cM","cMMb"))
rec[, Chr := sub("chromosome_", "Chr", chr)]
rc <- function(ch, pos) { o <- rep(NA_real_, length(pos))
  for (cc in unique(ch)) { r <- rec[Chr == paste0("Chr", cc)]; if (nrow(r) < 2) next
    i <- which(ch == cc); o[i] <- approx(r$p, r$cMMb, xout = pos[i], rule = 2)$y }; o }
merge_iv <- function(s, e) { o <- order(s); s <- s[o]; e <- e[o]
  cs <- s[1L]; ce <- e[1L]; oS <- numeric(0); oE <- numeric(0)
  for (i in seq_along(s)[-1L]) { if (s[i] <= ce) ce <- max(ce, e[i]) else { oS <- c(oS,cs); oE <- c(oE,ce); cs <- s[i]; ce <- e[i] } }
  list(s = c(oS, cs), e = c(oE, ce)) }
in_iv <- function(q, iv) { if (!length(iv$s)) return(logical(length(q))); (findInterval(q, as.vector(rbind(iv$s, iv$e))) %% 2L) == 1L }

ps <- as.data.table(readRDS("module_di25/data/di25_sorting_snp.rds"))
ps[, `:=`(chr = as.integer(sub("Chr","",sub(":.*","",marker))), pos = as.integer(sub(".*:","",marker)))]
ps[, recomb := rc(chr, pos)]
ps <- ps[differentiated == TRUE & n_obs > 0 & is.finite(recomb) & is.finite(uni_score)]
ps[, sorted := classify_sort(n_aqu, n_pol, n_obs, sort_th = TAU, sort_rule = "binom", alpha = 0.05) != "unsorted"]

beds <- sort(list.files(BEDDIR, pattern = "^bdmi_candidates\\.cutoff_.*\\.bed$"))
cut_k <- as.integer(sub(".*cutoff_(\\d+)_.*", "\\1", beds))
x2v   <- as.numeric(sub("^(0)(\\d+)$", "0.\\2", sub(".*cutoff_\\d+_(\\d+)\\..*", "\\1", beds)))
bi <- which(cut_k == A_CUT); bed <- fread(file.path(BEDDIR, beds[bi]), header = FALSE, col.names = c("chr","start","end"))
bed[, cn := as.integer(sub("chromosome_","",chr))]
ps[, inb := FALSE]; for (cc in unique(bed$cn)) { iv <- merge_iv(bed[cn==cc]$start, bed[cn==cc]$end); idx <- ps$chr == cc; ps$inb[idx] <- in_iv(ps$pos[idx], iv) }

## recombination bins -> composition among active loci
qb <- quantile(ps$recomb, 0:NBIN/NBIN, na.rm = TRUE)
ps[, rbin := cut(recomb, qb, include.lowest = TRUE, labels = FALSE)]
comp <- ps[!is.na(rbin), .(rmed = median(recomb),
                           s_only = sum(sorted & !inb), both = sum(sorted & inb), b_only = sum(!sorted & inb)),
           by = rbin][order(rbin)]
comp[, active := s_only + both + b_only]
comp <- comp[active > 0]
comp[, `:=`(a = s_only/active, m = both/active, c = b_only/active)]   # sorted-only / both / BDMI-only

## barycentric -> cartesian (L=sorted-only (0,0), R=BDMI-only (1,0), T=both (0.5,h))
h <- sqrt(3)/2
comp[, `:=`(x = c*1 + m*0.5, y = m*h)]

## colour by recombination rate (log scale) -> viridis
pal <- hcl.colors(100, "viridis")
rr  <- log10(comp$rmed + 0.1); col_i <- pmax(1, pmin(100, round(1 + 99 * (rr - min(rr)) / diff(range(rr)))))
comp[, col := pal[col_i]]

png(OUTPNG, width = 1900, height = 1900, res = 300, type = "cairo")
par(mar = c(1, 1, 3, 1), xpd = NA)
plot(NA, xlim = c(-0.12, 1.12), ylim = c(-0.12, h + 0.14), asp = 1, axes = FALSE, xlab = "", ylab = "",
     main = "")
title(main = sprintf("Recombination gradient of the sorting/BDMI mix  (X^2=%.2f)", x2v[bi]),
      cex.main = 0.92, line = 0.5)
## triangle + light gridlines
polygon(c(0, 1, 0.5), c(0, 0, h), border = "grey55", lwd = 1.4)
for (f in c(.25, .5, .75)) {
  segments(f*0.5, f*h, 1 - f*0.5, f*h, col = "grey88")              # horizontal (const 'both')
  segments(f, 0, 0.5 + f*0.5, (1-f)*h, col = "grey88")              # const BDMI-only
  segments(1-f, 0, 0.5 - f*0.5, (1-f)*h, col = "grey88")            # const sorted-only
}
## corner labels
text(0, -0.04, "sorted only", adj = c(0.5, 1), font = 2, cex = 0.85, col = "#8A3A0E")
text(1, -0.04, "BDMI only",   adj = c(0.5, 1), font = 2, cex = 0.85, col = "#B4661E")
text(0.5, h + 0.05, "both (sorted & BDMI)", adj = c(0.5, 0), font = 2, cex = 0.85, col = "#21918C")
## path + points (joined in recombination order)
lines(comp$x, comp$y, col = "grey60", lwd = 1)
points(comp$x, comp$y, pch = 21, bg = comp$col, col = "grey25", cex = 1.4 + 2.2 * sqrt(comp$active / max(comp$active)))
## endpoints annotated
text(comp$x[1], comp$y[1], sprintf("%.1f cM/Mb", comp$rmed[1]), pos = 2, cex = 0.7, col = "grey30")
n <- nrow(comp); text(comp$x[n], comp$y[n], sprintf("%.0f cM/Mb", comp$rmed[n]), pos = 4, cex = 0.7, col = "grey30")
## colour legend (recombination)
xl <- 0.86; yl <- seq(h*0.55, h*0.95, length.out = 100)
for (i in 1:99) rect(xl, yl[i], xl + 0.03, yl[i+1], col = pal[i], border = NA)
text(xl + 0.05, yl[100], sprintf("%.0f", 10^max(rr)), adj = c(0,1), cex = 0.65, col = "grey30")
text(xl + 0.05, yl[1],   sprintf("%.1f", 10^min(rr)), adj = c(0,0), cex = 0.65, col = "grey30")
text(xl + 0.015, yl[100] + 0.03, "cM/Mb", adj = c(0.5, 0), cex = 0.65, col = "grey30")
dev.off()
cat("[out] wrote", OUTPNG, "\n")
cat("[ternary] bins low->high recomb (sorted-only / both / BDMI-only):\n")
print(comp[, .(rmed = round(rmed,1), active, sorted_only = round(a,2), both = round(m,2), bdmi_only = round(c,2))])
