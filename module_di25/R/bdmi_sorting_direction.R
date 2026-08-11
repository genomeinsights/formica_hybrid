## =========================================================
## module_di25 -- BDMI regions vs SORTING DIRECTION (concordance)
## =========================================================
## The overlap test (bdmi_sorting_overlap.R) shows diagnostic SNPs inside BDMI
## candidate regions are enriched for SORTED loci, more strongly toward polyctena.
## This asks the sharper directional question: CONDITIONAL on a locus being
## directionally sorted, are the ones inside BDMI regions disproportionately
## POLYCTENA-directional relative to the genome-wide aquilonia:polyctena ratio?
## If BDMI incompatibilities drive polyctena fixation, in-region sorted loci
## should lean polyctena beyond chance positioning.
##
## Statistic per X^2 cutoff: pol fraction among in-region directional-sorted loci
##   f_in = pol_in / (pol_in + aqu_in),  vs genome-wide f0 = pol / (pol + aqu).
## Odds ratio OR = (pol_in/aqu_in) / (pol_out/aqu_out).
## Significance: per-chromosome circular-rotation null (same as the overlap test)
## -- randomise only the BDMI phase, keep each chromosome's sorted-direction
## pattern fixed; one-sided p that f_in exceeds the null.
##
##   EXCLUDE_CHR=5,25,26 drops the three big polyctena-block chromosomes.
## Run from repo root:  EXCLUDE_CHR=5,25,26 Rscript module_di25/R/bdmi_sorting_direction.R
## =========================================================

suppressMessages(library(data.table))
source("moduleA_sorting/R/parallelism_stats.R")   # classify_sort()

set.seed(1)
N_PERM  <- 2000L
TAU     <- 0.6
BEDDIR  <- "data/liftoff_Frufa_DTOL_PR"
EXCL <- Sys.getenv("EXCLUDE_CHR", "")
excl_num <- if (nzchar(EXCL)) as.integer(strsplit(EXCL, ",")[[1]]) else integer(0)
sfx     <- if (length(excl_num)) paste0("_no", paste(excl_num, collapse = "_")) else ""
OUTRDS  <- sprintf("module_di25/data/bdmi_sorting_direction%s.rds", sfx)
OUTPNG  <- sprintf("module_di25/Figures/bdmi_sorting_direction%s.png", sfx)

## ---- diagnostic SNPs: position + direction at TAU --------------------------
ps  <- readRDS("module_di25/data/di25_sorting_snp.rds")
chr <- sub(":.*", "", ps$marker); pos <- as.integer(sub(".*:", "", ps$marker))
ok  <- ps$differentiated & ps$n_obs > 0 & is.finite(ps$uni_score)
cls <- rep("unsorted", nrow(ps))
cls[ok] <- classify_sort(ps$n_aqu[ok], ps$n_pol[ok], ps$n_obs[ok],
                         sort_th = TAU, sort_rule = "binom", alpha = 0.05)
snp <- data.table(chr = chr, pos = pos,
                  pol = cls == "polyctena", aqu = cls == "aquilonia")
if (length(excl_num)) snp <- snp[!chr %in% paste0("Chr", excl_num)]
setkey(snp, chr, pos)
chr_len <- snp[, .(len = max(pos)), by = chr]

f0 <- sum(snp$pol) / (sum(snp$pol) + sum(snp$aqu))   # genome-wide polyctena fraction among directional-sorted
cat(sprintf("[snp] %d SNPs; directional-sorted pol %d / aqu %d -> genome-wide pol fraction f0 = %.3f%s\n",
            as.integer(nrow(snp)), as.integer(sum(snp$pol)), as.integer(sum(snp$aqu)), f0,
            if (length(excl_num)) sprintf(" (excl Chr%s)", paste(excl_num, collapse = ",")) else ""))

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

## ---- per-cutoff direction concordance + rotation null ----------------------
beds <- sort(list.files(BEDDIR, pattern = "^bdmi_candidates\\.cutoff_.*\\.bed$"))
cut_k  <- as.integer(sub(".*cutoff_(\\d+)_.*", "\\1", beds))
x2_val <- as.numeric(sub("^(0)(\\d+)$", "0.\\2", sub(".*cutoff_\\d+_(\\d+)\\..*", "\\1", beds)))
ord <- order(cut_k); beds <- beds[ord]; cut_k <- cut_k[ord]; x2_val <- x2_val[ord]

res <- vector("list", length(beds))
for (i in seq_along(beds)) {
  bed <- fread(file.path(BEDDIR, beds[i]), header = FALSE, col.names = c("chr", "start", "end"))
  bed[, chr := sub("chromosome_", "Chr", chr)]
  bed <- bed[chr %in% chr_len$chr]
  ivs <- lapply(split(bed, bed$chr), function(b) merge_iv(b$start, b$end))

  snp[, inb := FALSE]
  for (cc in names(ivs)) { idx <- snp$chr == cc; snp$inb[idx] <- in_intervals(snp$pos[idx], ivs[[cc]]) }
  pol_in <- sum(snp$inb & snp$pol);  aqu_in <- sum(snp$inb & snp$aqu)
  pol_ou <- sum(!snp$inb & snp$pol); aqu_ou <- sum(!snp$inb & snp$aqu)
  f_in   <- pol_in / max(pol_in + aqu_in, 1L)
  OR     <- (pol_in / max(aqu_in, 1L)) / (pol_ou / max(aqu_ou, 1L))

  null_fin <- numeric(N_PERM)
  for (p in seq_len(N_PERM)) {
    pin <- 0L; ain <- 0L
    for (cc in names(ivs)) {
      idx <- which(snp$chr == cc); if (!length(idx)) next
      L <- chr_len[chr == cc, len]; off <- sample.int(L, 1L)
      qp <- ((snp$pos[idx] - 1L + off) %% L) + 1L
      m  <- in_intervals(qp, ivs[[cc]])
      pin <- pin + sum(m & snp$pol[idx]); ain <- ain + sum(m & snp$aqu[idx])
    }
    null_fin[p] <- pin / max(pin + ain, 1L)
  }
  p_emp <- (1 + sum(null_fin >= f_in)) / (1 + N_PERM)

  res[[i]] <- data.table(cutoff_k = cut_k[i], x2 = x2_val[i], n_region = nrow(bed),
                         pol_in = pol_in, aqu_in = aqu_in, f_in = f_in, f0 = f0,
                         OR = OR, null_fin_mean = mean(null_fin), p = p_emp)
  cat(sprintf("[cutoff %2d  X2=%.5f] sorted-in pol/aqu = %3d/%3d  pol-frac %.3f vs f0 %.3f  OR %.2f  p=%.4f\n",
              cut_k[i], x2_val[i], pol_in, aqu_in, f_in, f0, OR, p_emp))
}
res <- rbindlist(res); saveRDS(res, OUTRDS); cat("[out] wrote", OUTRDS, "\n")

## ---- figure ----------------------------------------------------------------
xlab_e <- expression(""%<-%"more permissive    "*"X"^2*" cutoff"*"    more stringent"%->%"")
png(OUTPNG, width = 2100, height = 1600, res = 300, type = "cairo")
op <- par(mfrow = c(2, 1), mar = c(4.4, 5.4, 2.0, 1.2), mgp = c(2.6, 0.7, 0))
plot(res$x2, res$f_in, type = "b", pch = 19, log = "x", col = "#B4661E",
     xlab = xlab_e, ylab = "polyctena fraction among\nsorted SNPs in BDMI regions",
     ylim = range(0, res$f_in, f0))
abline(h = f0, lty = 2, col = "grey45")
text(max(res$x2), f0, sprintf("genome-wide = %.2f", f0), adj = c(1, -0.4), cex = 0.7, col = "grey35")
plot(res$x2, -log10(res$p), type = "b", pch = 19, log = "x", col = "#B4661E",
     xlab = xlab_e, ylab = expression(-log[10]~"(rotation-null p)"),
     ylim = c(0, max(3.2, -log10(res$p))))
abline(h = -log10(0.05), lty = 2, col = "grey55")
text(min(res$x2), -log10(0.05), "p = 0.05", pos = 3, cex = 0.7, col = "grey45")
par(op); dev.off(); cat("[out] wrote", OUTPNG, "\n")
