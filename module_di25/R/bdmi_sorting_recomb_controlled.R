## =========================================================
## module_di25 -- BDMI<->sorting association, controlling for recombination
## =========================================================
## Figs 4-6 show BDMI candidate regions overlap ancestry-sorted loci; Fig 7 shows
## BDMIs and sorting have OPPOSITE recombination tendencies (sorting low-recomb,
## BDMIs high-recomb). So: does the BDMI<->sorting association survive controlling
## for recombination rate? Because BDMIs sit in high-recomb / low-sorting
## territory, adjusting for recombination should if anything STRENGTHEN it.
##
## Panel a (intuitive): P(sorted) across recombination deciles, split by whether
##   the SNP is inside a BDMI region -- at one representative cutoff. If the
##   inside line sits above the outside line within every decile, the association
##   is not a recombination artefact.
## Panel b (test): the BDMI term of a logistic regression of P(sorted) on BDMI
##   membership, UNADJUSTED  (sorted ~ inBDMI)  vs  recombination-ADJUSTED
##   (sorted ~ inBDMI + log10 recomb), per X^2 cutoff, chromosome block bootstrap.
##
## Run from repo root:  Rscript module_di25/R/bdmi_sorting_recomb_controlled.R
## =========================================================

suppressMessages(library(data.table))
source("moduleA_sorting/R/parallelism_stats.R")   # classify_sort()

set.seed(1)
NB      <- 1000L
TAU     <- 0.6
A_CUT   <- 13L                                     # panel-a cutoff (X2=0.03, good coverage)
BEDDIR  <- "data/liftoff_Frufa_DTOL_PR"
OUTPNG  <- "module_di25/Figures/bdmi_sorting_recomb_controlled.png"
OUTRDS  <- "module_di25/data/bdmi_sorting_recomb_controlled.rds"

## ---- diagnostic SNPs: recombination + sorted status ------------------------
rec <- fread("data/Frufa_DTOL_PR.ref_genome.recmap"); setnames(rec, 1:4, c("chr","p","cM","cMMb"))
rec[, Chr := sub("chromosome_", "Chr", chr)]
rc <- function(ch, pos) { o <- rep(NA_real_, length(pos))
  for (cc in unique(ch)) { r <- rec[Chr == paste0("Chr", cc)]; if (nrow(r) < 2) next
    i <- which(ch == cc); o[i] <- approx(r$p, r$cMMb, xout = pos[i], rule = 2)$y }; o }
wilson <- function(k, n) { z <- 1.959964; if (n == 0) return(c(NA, NA)); p <- k/n; d <- 1 + z^2/n
  ctr <- (p + z^2/(2*n))/d; hw <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/d; c(max(0, ctr-hw), min(1, ctr+hw)) }
merge_iv <- function(s, e) { o <- order(s); s <- s[o]; e <- e[o]
  cs <- s[1L]; ce <- e[1L]; oS <- numeric(0); oE <- numeric(0)
  for (i in seq_along(s)[-1L]) { if (s[i] <= ce) ce <- max(ce, e[i]) else { oS <- c(oS,cs); oE <- c(oE,ce); cs <- s[i]; ce <- e[i] } }
  list(s = c(oS, cs), e = c(oE, ce)) }
in_iv <- function(q, iv) { if (!length(iv$s)) return(logical(length(q)))
  (findInterval(q, as.vector(rbind(iv$s, iv$e))) %% 2L) == 1L }

ps <- as.data.table(readRDS("module_di25/data/di25_sorting_snp.rds"))
ps[, `:=`(chr = as.integer(sub("Chr","",sub(":.*","",marker))), pos = as.integer(sub(".*:","",marker)))]
ps[, recomb := rc(chr, pos)]
ps <- ps[differentiated == TRUE & n_obs > 0 & is.finite(recomb) & is.finite(uni_score)]
ps[, `:=`(lr = log10(recomb + 0.1),
          sorted = as.integer(classify_sort(n_aqu, n_pol, n_obs, sort_th = TAU, sort_rule = "binom", alpha = 0.05) != "unsorted"))]

beds <- sort(list.files(BEDDIR, pattern = "^bdmi_candidates\\.cutoff_.*\\.bed$"))
cut_k <- as.integer(sub(".*cutoff_(\\d+)_.*", "\\1", beds))
x2v   <- as.numeric(sub("^(0)(\\d+)$", "0.\\2", sub(".*cutoff_\\d+_(\\d+)\\..*", "\\1", beds)))
o <- order(cut_k); beds <- beds[o]; cut_k <- cut_k[o]; x2v <- x2v[o]
member <- function(i) { bed <- fread(file.path(BEDDIR, beds[i]), header = FALSE, col.names = c("chr","start","end"))
  bed[, cn := as.integer(sub("chromosome_","",chr))]; inb <- logical(nrow(ps))
  for (cc in unique(bed$cn)) { iv <- merge_iv(bed[cn==cc]$start, bed[cn==cc]$end); idx <- ps$chr == cc
    inb[idx] <- in_iv(ps$pos[idx], iv) }; inb }

## ---- panel-a data: P(sorted) by recomb decile, inside vs outside BDMI -------
rbrk <- quantile(ps$recomb, 0:10/10, na.rm = TRUE)
ps[, rdec := cut(recomb, rbrk, include.lowest = TRUE, labels = FALSE)]
xlab <- ps[!is.na(rdec), round(median(recomb), 1), by = rdec][order(rdec)]$V1
inbA <- member(which(cut_k == A_CUT))
adat <- rbindlist(lapply(c(TRUE, FALSE), function(side) {
  d <- ps[inbA == side]
  d[, .(grp = side, n = .N, k = sum(sorted)), by = rdec][order(rdec)]
}))
adat[, `:=`(frac = k/n)]; adat[, c("lo","hi") := as.list(wilson(k, n)), by = .(grp, rdec)]

## ---- panel-b: BDMI log-odds on sorting, unadjusted vs recomb-adjusted -------
chr_idx <- split(seq_len(nrow(ps)), ps$chr); chrs <- names(chr_idx)
## return c(unadj_bdmi, adj_bdmi) from ONE index resample
fit_two <- function(y, inb, idx) {
  yy <- y[idx]; b1 <- inb[idx]
  if (sum(b1) < 3L || sum(yy[b1]) < 2L || sum(yy[!b1]) < 2L) return(c(NA, NA))
  u <- suppressWarnings(tryCatch(glm.fit(cbind(1, b1), yy, family = binomial())$coefficients[2], error = function(e) NA_real_))
  a <- suppressWarnings(tryCatch(glm.fit(cbind(1, b1, ps$lr[idx]), yy, family = binomial())$coefficients[2], error = function(e) NA_real_))
  c(u, a)
}
## the block bootstrap is the only slow part; reuse it if already computed
## (delete module_di25/data/bdmi_sorting_recomb_controlled.rds to force a re-run)
if (file.exists(OUTRDS)) {
  res <- readRDS(OUTRDS); cat("[cache] reusing", OUTRDS, "-- only re-plotting\n")
} else {
  res <- vector("list", length(cut_k))
  for (i in seq_along(cut_k)) {
    inb <- member(i); y <- ps$sorted
    est <- fit_two(y, inb, seq_len(nrow(ps)))
    bs  <- vapply(seq_len(NB), function(b) fit_two(y, inb, unlist(chr_idx[sample(chrs, length(chrs), TRUE)], use.names = FALSE)), numeric(2))
    res[[i]] <- data.table(x2 = x2v[i], n_in = sum(inb),
      unadj = est[1], unadj_lo = quantile(bs[1,], .025, na.rm = TRUE), unadj_hi = quantile(bs[1,], .975, na.rm = TRUE),
      adj   = est[2], adj_lo   = quantile(bs[2,], .025, na.rm = TRUE), adj_hi   = quantile(bs[2,], .975, na.rm = TRUE))
    cat(sprintf("[cutoff X2=%.4f] BDMI log-odds on sorting: unadjusted %.2f [%.2f,%.2f]  recomb-adjusted %.2f [%.2f,%.2f]\n",
                x2v[i], res[[i]]$unadj, res[[i]]$unadj_lo, res[[i]]$unadj_hi, res[[i]]$adj, res[[i]]$adj_lo, res[[i]]$adj_hi))
  }
  res <- rbindlist(res); saveRDS(res, OUTRDS)
}

## ---- figure (no panel titles; wording lives in the manuscript caption) -----
png(OUTPNG, width = 2100, height = 2000, res = 300, type = "cairo")
op <- par(mfrow = c(2, 1), mar = c(4.4, 4.6, 1.6, 1.2), mgp = c(2.6, 0.7, 0))

## panel a
ci <- c("in" = "#8A3A0E", "out" = "grey55")
with(adat[grp == TRUE], { plot(rdec, frac, type = "n", xaxt = "n", ylim = c(0, max(adat$hi, na.rm = TRUE)),
  xlab = "recombination decile (median cM/Mb)", ylab = "P(ancestry-sorted)") })
mtext("a", side = 3, line = 0.2, adj = 0, font = 2, cex = 1.1)
axis(1, at = 1:10, labels = xlab)
for (side in c(TRUE, FALSE)) { d <- adat[grp == side]; col <- if (side) ci["in"] else ci["out"]
  segments(d$rdec, d$lo, d$rdec, d$hi, col = col); lines(d$rdec, d$frac, col = col, lwd = 2)
  points(d$rdec, d$frac, col = col, pch = 19) }
legend("topright", bty = "n", lwd = 2, pch = 19, col = ci[c("in","out")],
       legend = c("inside a BDMI region", "outside"))

## panel b
xr <- range(x2v)
with(res, { plot(x2, adj, type = "n", log = "x", ylim = range(0, unadj_lo, unadj_hi, adj_lo, adj_hi, na.rm = TRUE),
  xlab = expression("X"^2*" cutoff  ("%<-%" more permissive; more stringent "%->%")"),
  ylab = "BDMI effect on P(sorted)\n(logit; log-odds)") })
mtext("b", side = 3, line = 0.2, adj = 0, font = 2, cex = 1.1)
abline(h = 0, lty = 2, col = "grey70")
dx <- 1.02
with(res, segments(x2/dx, unadj_lo, x2/dx, unadj_hi, col = "grey55"))
with(res, { points(x2/dx, unadj, pch = 17, col = "grey45"); lines(x2/dx, unadj, col = "grey55", lty = 3) })
with(res, segments(x2*dx, adj_lo, x2*dx, adj_hi, col = "#8A3A0E"))
with(res, { points(x2*dx, adj, pch = 19, col = "#8A3A0E"); lines(x2*dx, adj, col = "#8A3A0E") })
legend("topleft", bty = "n", pch = c(17, 19), lty = c(3, 1), col = c("grey45", "#8A3A0E"),
       legend = c("unadjusted", "recombination-adjusted"))
par(op); dev.off()
cat("[out] wrote", OUTPNG, "and", OUTRDS, "\n")
