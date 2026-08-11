## =========================================================
## module_di25 -- do BDMI candidate regions overlap sorted regions?
## =========================================================
## A colleague provided per-cutoff BDMI (Bateson-Dobzhansky-Muller
## incompatibility) candidate regions from male hybrids of one F. aquilonia x
## F. polyctena population, lifted onto the Frufa_DTOL_PR reference:
##   data/liftoff_Frufa_DTOL_PR/bdmi_candidates.cutoff_<k>_<X2>.liftoff.*.bed
## Each .bed is the set of merged genomic INTERVALS (nodes) whose loci take part
## in >=1 significant DMI edge at that X^2 cutoff. Smaller X^2 cutoff (higher k)
## = more permissive = more/larger regions.
##
## Question: are the diagnostic SNPs that fall inside these BDMI regions more
## often ANCESTRY-SORTED (toward F. aquilonia / F. polyctena) than expected?
##
## Test: per-chromosome CIRCULAR-ROTATION null. Each chromosome's diagnostic-SNP
## sort labels and its BDMI interval coverage are held fixed; only the PHASE of
## the BDMI regions relative to the SNPs is randomised (shift the SNP coordinates
## by a random per-chromosome offset, wrap within the chromosome, recount).
## A positive result therefore means the *real* regions genuinely sit on sorted
## loci -- not merely that both like the same chromosomes or the big blocks.
##
## Run from repo root:  Rscript module_di25/R/bdmi_sorting_overlap.R
## =========================================================

suppressMessages(library(data.table))
source("moduleA_sorting/R/parallelism_stats.R")   # classify_sort()

set.seed(1)
N_PERM  <- 2000L
TAU     <- 0.6                                     # primary sorting operating point
PHI     <- 0.85                                    # (fixed; matches the tau-sweep circos)
BEDDIR  <- "data/liftoff_Frufa_DTOL_PR"
## EXCLUDE_CHR: comma-separated chromosome numbers dropped from BOTH the SNP pool
## and the BDMI regions (robustness check). e.g. EXCLUDE_CHR=5,25,26 removes the
## three big low-recombination polyctena blocks' chromosomes.
EXCL <- Sys.getenv("EXCLUDE_CHR", "")
excl_num <- if (nzchar(EXCL)) as.integer(strsplit(EXCL, ",")[[1]]) else integer(0)
sfx     <- if (length(excl_num)) paste0("_no", paste(excl_num, collapse = "_")) else ""
OUTRDS  <- sprintf("module_di25/data/bdmi_sorting_overlap%s.rds", sfx)
OUTPNG  <- sprintf("module_di25/Figures/bdmi_sorting_overlap%s.png", sfx)

## ---- 1. diagnostic SNPs: position + sorted status at TAU --------------------
ps  <- readRDS("module_di25/data/di25_sorting_snp.rds")
chr <- sub(":.*", "", ps$marker)
pos <- as.integer(sub(".*:", "", ps$marker))
ok  <- ps$differentiated & ps$n_obs > 0 & is.finite(ps$uni_score)

cls <- rep("unsorted", nrow(ps))
cls[ok] <- classify_sort(ps$n_aqu[ok], ps$n_pol[ok], ps$n_obs[ok],
                         sort_th = TAU, sort_rule = "binom", alpha = 0.05)
sorted <- cls != "unsorted"                        # any direction
sort_pol <- cls == "polyctena"
sort_aqu <- cls == "aquilonia"

snp <- data.table(chr = chr, pos = pos, sorted = sorted,
                  pol = sort_pol, aqu = sort_aqu)
if (length(excl_num)) {                            # drop excluded chromosomes
  drop <- snp$chr %in% paste0("Chr", excl_num)
  cat(sprintf("[excl] dropping %s: %d SNPs removed\n",
              paste0("Chr", excl_num, collapse = ","), as.integer(sum(drop))))
  snp <- snp[!drop]
  cls <- cls[!drop]; sort_pol <- snp$pol; sort_aqu <- snp$aqu
}
setkey(snp, chr, pos)
chr_levels <- sort(unique(snp$chr))
## per-chromosome span used for the circular rotation (cover SNPs and regions)
chr_len <- snp[, .(len = max(pos)), by = chr]

overall_sorted <- mean(snp$sorted)
cat(sprintf("[snp] %d diagnostic SNPs; sorted at tau=%.1f: %d (%.1f%%)  [pol %d / aqu %d / unresolved %d]\n",
            as.integer(nrow(snp)), TAU, as.integer(sum(snp$sorted)), 100 * overall_sorted,
            as.integer(sum(sort_pol)), as.integer(sum(sort_aqu)), as.integer(sum(cls == "unresolved"))))

## collapse a set of [start,end] intervals into sorted, non-overlapping ones
## (the beds are mostly merged, but a few cutoffs still contain overlaps).
merge_iv <- function(s, e) {
  o <- order(s); s <- s[o]; e <- e[o]
  cs <- s[1L]; ce <- e[1L]; outS <- numeric(0); outE <- numeric(0)
  for (i in seq_along(s)[-1L]) {
    if (s[i] <= ce) ce <- max(ce, e[i])
    else { outS <- c(outS, cs); outE <- c(outE, ce); cs <- s[i]; ce <- e[i] }
  }
  list(s = c(outS, cs), e = c(outE, ce))
}
## membership test against sorted, non-overlapping [start,end] intervals:
## findInterval on the interleaved boundary vector; inside == odd index.
in_intervals <- function(qpos, iv) {
  if (!length(iv$s)) return(logical(length(qpos)))
  brk <- as.vector(rbind(iv$s, iv$e))              # s1,e1,s2,e2,... (non-decreasing)
  (findInterval(qpos, brk) %% 2L) == 1L
}

## ---- 2. per-cutoff overlap + rotation null ---------------------------------
beds <- sort(list.files(BEDDIR, pattern = "^bdmi_candidates\\.cutoff_.*\\.bed$"))
cut_k  <- as.integer(sub(".*cutoff_(\\d+)_.*", "\\1", beds))
x2_str <- sub(".*cutoff_\\d+_(\\d+)\\..*", "\\1", beds)
x2_val <- as.numeric(sub("^(0)(\\d+)$", "0.\\2", x2_str))   # 00275 -> 0.00275
ord <- order(cut_k); beds <- beds[ord]; cut_k <- cut_k[ord]; x2_val <- x2_val[ord]

res <- vector("list", length(beds))
for (i in seq_along(beds)) {
  bed <- fread(file.path(BEDDIR, beds[i]), header = FALSE,
               col.names = c("chr", "start", "end"))
  bed[, chr := sub("chromosome_", "Chr", chr)]
  bed <- bed[chr %in% chr_levels]

  ## merged intervals per chromosome (reused by observed pass and every permutation)
  ivs <- lapply(split(bed, bed$chr), function(b) merge_iv(b$start, b$end))

  ## observed: mark each diagnostic SNP as inside a BDMI region
  snp[, inb := FALSE]
  for (cc in names(ivs)) {
    idx <- snp$chr == cc
    snp$inb[idx] <- in_intervals(snp$pos[idx], ivs[[cc]])
  }
  n_in       <- sum(snp$inb)
  obs_sorted <- sum(snp$inb & snp$sorted)
  obs_pol    <- sum(snp$inb & snp$pol)
  obs_aqu    <- sum(snp$inb & snp$aqu)

  ## null: rotate SNP coordinates within each chromosome, recount overlaps
  null_sorted <- integer(N_PERM); null_pol <- integer(N_PERM); null_aqu <- integer(N_PERM)
  for (p in seq_len(N_PERM)) {
    ns <- 0L; np <- 0L; na <- 0L
    for (cc in names(ivs)) {
      idx <- which(snp$chr == cc)
      if (!length(idx)) next
      L <- chr_len[chr == cc, len]
      off <- sample.int(L, 1L)
      qp  <- ((snp$pos[idx] - 1L + off) %% L) + 1L   # circular shift of SNP coords
      inb <- in_intervals(qp, ivs[[cc]])
      ns <- ns + sum(inb & snp$sorted[idx])
      np <- np + sum(inb & snp$pol[idx])
      na <- na + sum(inb & snp$aqu[idx])
    }
    null_sorted[p] <- ns; null_pol[p] <- np; null_aqu[p] <- na
  }

  emp_p <- function(obsv, nullv) (1 + sum(nullv >= obsv)) / (1 + length(nullv))
  fold  <- function(obsv, nullv) obsv / max(mean(nullv), .Machine$double.eps)

  res[[i]] <- data.table(
    cutoff_k = cut_k[i], x2 = x2_val[i],
    n_region = nrow(bed), n_snp_in = n_in,
    frac_in_sorted = obs_sorted / max(n_in, 1L),
    overall_sorted = overall_sorted,
    obs_sorted = obs_sorted, null_sorted_mean = mean(null_sorted),
    fold_sorted = fold(obs_sorted, null_sorted), p_sorted = emp_p(obs_sorted, null_sorted),
    obs_pol = obs_pol, fold_pol = fold(obs_pol, null_pol), p_pol = emp_p(obs_pol, null_pol),
    obs_aqu = obs_aqu, fold_aqu = fold(obs_aqu, null_aqu), p_aqu = emp_p(obs_aqu, null_aqu))
  cat(sprintf("[cutoff %2d  X2=%.5f] %4d regions, %5d SNPs inside | sorted in-frac %.3f vs %.3f  fold %.2f  p=%.4f  (pol fold %.2f p=%.4f / aqu fold %.2f p=%.4f)\n",
              cut_k[i], x2_val[i], nrow(bed), n_in,
              res[[i]]$frac_in_sorted, overall_sorted, res[[i]]$fold_sorted, res[[i]]$p_sorted,
              res[[i]]$fold_pol, res[[i]]$p_pol, res[[i]]$fold_aqu, res[[i]]$p_aqu))
}
res <- rbindlist(res)
saveRDS(res, OUTRDS)
cat("[out] wrote", OUTRDS, "\n")

## ---- 3. summary figure -----------------------------------------------------
plot_overlap <- function(res, outpng) {
  xlab_e <- expression(""%<-%"more permissive    "*"X"^2*" cutoff"*"    more stringent"%->%"")
  png(outpng, width = 2100, height = 1600, res = 300, type = "cairo")
  op <- par(mfrow = c(2, 1), mar = c(4.4, 5.4, 1.8, 1.2), mgp = c(2.6, 0.7, 0))

  ## panel a: fold-enrichment of sorted SNPs inside BDMI regions vs cutoff
  plot(res$x2, res$fold_sorted, type = "b", pch = 19, log = "x", col = "#315B7D",
       xlab = xlab_e, ylab = "fold-enrichment of\nsorted SNPs in BDMI regions",
       ylim = range(1, res$fold_sorted, res$fold_pol, res$fold_aqu))
  abline(h = 1, lty = 2, col = "grey55")
  lines(res$x2, res$fold_pol, type = "b", pch = 17, col = "#D3C93B")
  lines(res$x2, res$fold_aqu, type = "b", pch = 15, col = "#21918C")
  legend("topleft", bty = "n", cex = 0.8,
         legend = c("any sorted", "toward polyctena", "toward aquilonia"),
         col = c("#315B7D", "#D3C93B", "#21918C"), pch = c(19, 17, 15), lty = 1)

  ## panel b: rotation-null p-value vs cutoff
  plot(res$x2, -log10(res$p_sorted), type = "b", pch = 19, log = "x", col = "#315B7D",
       xlab = xlab_e, ylab = expression(-log[10]~"(rotation-null p)"),
       ylim = c(0, max(3.2, -log10(res$p_sorted))))
  abline(h = -log10(0.05), lty = 2, col = "grey55")
  lines(res$x2, -log10(res$p_pol), type = "b", pch = 17, col = "#D3C93B")
  lines(res$x2, -log10(res$p_aqu), type = "b", pch = 15, col = "#21918C")
  text(min(res$x2), -log10(0.05), "p = 0.05", pos = 3, cex = 0.7, col = "grey45")
  par(op); dev.off()
}
plot_overlap(res, OUTPNG)
cat("[out] wrote", OUTPNG, "\n")
