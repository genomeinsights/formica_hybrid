## =========================================================
## module_di25 -- rotation-null histogram for the BDMI/sorting overlap test
## =========================================================
## Visualises the significance test behind Figure 5: for ONE X^2 cutoff, the
## per-chromosome circular-rotation null distribution of the number of SORTED
## diagnostic SNPs falling inside the BDMI candidate regions, with the OBSERVED
## count marked. Three panels: any-direction sorted, toward polyctena, toward
## aquilonia. More permutations than the sweep (one cutoff only) for a smooth null.
##
##   CUTOFF=k       which bed (default 5 -> X2=0.05, a stringent, strongly-enriched set)
##   EXCLUDE_CHR=   comma-separated chromosomes to drop (e.g. 5,25,26)
##   N_PERM=        permutations (default 5000)
## Run from repo root:  CUTOFF=5 Rscript module_di25/R/bdmi_sorting_null_hist.R
## =========================================================

suppressMessages(library(data.table))
source("moduleA_sorting/R/parallelism_stats.R")   # classify_sort()

set.seed(1)
CUTOFF_K <- as.integer(Sys.getenv("CUTOFF", "5"))
N_PERM   <- as.integer(Sys.getenv("N_PERM", "5000"))
TAU      <- 0.6
BEDDIR   <- "data/liftoff_Frufa_DTOL_PR"
EXCL <- Sys.getenv("EXCLUDE_CHR", "")
excl_num <- if (nzchar(EXCL)) as.integer(strsplit(EXCL, ",")[[1]]) else integer(0)
sfx  <- if (length(excl_num)) paste0("_no", paste(excl_num, collapse = "_")) else ""

## ---- diagnostic SNPs: position + sort class at TAU -------------------------
ps  <- readRDS("module_di25/data/di25_sorting_snp.rds")
chr <- sub(":.*", "", ps$marker); pos <- as.integer(sub(".*:", "", ps$marker))
ok  <- ps$differentiated & ps$n_obs > 0 & is.finite(ps$uni_score)
cls <- rep("unsorted", nrow(ps))
cls[ok] <- classify_sort(ps$n_aqu[ok], ps$n_pol[ok], ps$n_obs[ok],
                         sort_th = TAU, sort_rule = "binom", alpha = 0.05)
snp <- data.table(chr = chr, pos = pos, sorted = cls != "unsorted",
                  pol = cls == "polyctena", aqu = cls == "aquilonia")
if (length(excl_num)) snp <- snp[!chr %in% paste0("Chr", excl_num)]
setkey(snp, chr, pos)
chr_len <- snp[, .(len = max(pos)), by = chr]

merge_iv <- function(s, e) {
  o <- order(s); s <- s[o]; e <- e[o]
  cs <- s[1L]; ce <- e[1L]; outS <- numeric(0); outE <- numeric(0)
  for (i in seq_along(s)[-1L]) {
    if (s[i] <= ce) ce <- max(ce, e[i])
    else { outS <- c(outS, cs); outE <- c(outE, ce); cs <- s[i]; ce <- e[i] }
  }
  list(s = c(outS, cs), e = c(outE, ce))
}
in_intervals <- function(qpos, iv) {
  if (!length(iv$s)) return(logical(length(qpos)))
  brk <- as.vector(rbind(iv$s, iv$e)); (findInterval(qpos, brk) %% 2L) == 1L
}

## ---- load the chosen cutoff's BDMI regions ---------------------------------
bedf <- list.files(BEDDIR, pattern = sprintf("^bdmi_candidates\\.cutoff_%d_.*\\.bed$", CUTOFF_K))
stopifnot(length(bedf) == 1L)
x2  <- as.numeric(sub("^(0)(\\d+)$", "0.\\2", sub(".*cutoff_\\d+_(\\d+)\\..*", "\\1", bedf)))
bed <- fread(file.path(BEDDIR, bedf), header = FALSE, col.names = c("chr", "start", "end"))
bed[, chr := sub("chromosome_", "Chr", chr)]; bed <- bed[chr %in% chr_len$chr]
ivs <- lapply(split(bed, bed$chr), function(b) merge_iv(b$start, b$end))

## observed
snp[, inb := FALSE]
for (cc in names(ivs)) { idx <- snp$chr == cc; snp$inb[idx] <- in_intervals(snp$pos[idx], ivs[[cc]]) }
obs <- c(any = sum(snp$inb & snp$sorted), pol = sum(snp$inb & snp$pol), aqu = sum(snp$inb & snp$aqu))
n_in <- sum(snp$inb)

## rotation null
null <- matrix(0L, N_PERM, 3, dimnames = list(NULL, c("any", "pol", "aqu")))
for (p in seq_len(N_PERM)) {
  s <- 0L; po <- 0L; a <- 0L
  for (cc in names(ivs)) {
    idx <- which(snp$chr == cc); if (!length(idx)) next
    L <- chr_len[chr == cc, len]; off <- sample.int(L, 1L)
    m <- in_intervals(((snp$pos[idx] - 1L + off) %% L) + 1L, ivs[[cc]])
    s <- s + sum(m & snp$sorted[idx]); po <- po + sum(m & snp$pol[idx]); a <- a + sum(m & snp$aqu[idx])
  }
  null[p, ] <- c(s, po, a)
}
emp_p <- sapply(1:3, function(j) (1 + sum(null[, j] >= obs[j])) / (1 + N_PERM))
fold  <- obs / colMeans(null)

## ---- figure: one histogram panel per category ------------------------------
OUTPNG <- sprintf("module_di25/Figures/bdmi_sorting_null_hist_cutoff%d%s.png", CUTOFF_K, sfx)
labs <- c(any = "any sorted", pol = "toward F. polyctena", aqu = "toward F. aquilonia")
cols <- c(any = "#315B7D", pol = "#D3C93B", aqu = "#21918C")
png(OUTPNG, width = 2400, height = 900, res = 300, type = "cairo")
op <- par(mfrow = c(1, 3), mar = c(4.2, 4.0, 3.0, 1.0), mgp = c(2.3, 0.7, 0), oma = c(0, 0, 2.0, 0))
for (j in 1:3) {
  v <- null[, j]; ob <- obs[j]
  xr <- range(v, ob); br <- seq(xr[1] - 0.5, xr[2] + 0.5, length.out = 40)
  h <- hist(v, breaks = br, plot = FALSE)
  plot(h, col = "grey85", border = "grey70", main = labs[j], cex.main = 1.0,
       xlab = "sorted SNPs inside BDMI regions", ylab = "null permutations",
       xlim = c(min(br), max(xr[2], ob) * 1.02))
  abline(v = ob, col = cols[j], lwd = 2.5)
  ## annotations right-justified just LEFT of the observed line so they never overlap it
  ymax <- par("usr")[4]; gap <- 0.02 * diff(par("usr")[1:2])
  p_str <- if (emp_p[j] <= 1 / (N_PERM + 1)) sprintf("<%.1g", 1 / N_PERM) else sprintf("%.4f", emp_p[j])
  text(ob - gap, ymax * 0.96, sprintf("observed = %d", ob), col = cols[j], adj = c(1, 1), cex = 0.9, font = 2)
  text(ob - gap, ymax * 0.78,
       sprintf("null mean = %.0f\nfold = %.1fx\np = %s", mean(v), fold[j], p_str),
       adj = c(1, 1), cex = 0.8, col = "grey20")
}
mtext(sprintf("Rotation-null overlap of sorted diagnostic SNPs with BDMI regions  (X^2 cutoff %.3f, %d regions, %d SNPs inside%s)",
              x2, nrow(bed), n_in, if (length(excl_num)) sprintf(", excl Chr%s", paste(excl_num, collapse = ",")) else ""),
      outer = TRUE, cex = 0.72, font = 2)
par(op); dev.off()
cat(sprintf("[hist] cutoff %d (X2=%.3f)%s: obs any/pol/aqu = %d/%d/%d ; fold %.1f/%.1f/%.1f ; p %.4f/%.4f/%.4f\n",
            CUTOFF_K, x2, sfx, obs[1], obs[2], obs[3], fold[1], fold[2], fold[3], emp_p[1], emp_p[2], emp_p[3]))
cat("[out] wrote", OUTPNG, "\n")
