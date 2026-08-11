## =========================================================
## module_di25 -- BDMI candidate regions overlaid on the sorting circos
## =========================================================
## The tau-sweep sorting circos (diem_circos_sorting_sweep.R) with an extra OUTER
## ring marking the diagnostic-SNP units that fall inside a colleague's BDMI
## candidate regions at one X^2 cutoff. Lets you SEE whether the BDMI regions
## (orange ticks) sit on the ancestry-sorted units (coloured rings) rather than
## the unsorted background (near-white).
##
## Geometry is identical to the sweep, so the BDMI ring aligns unit-for-unit:
## units are ordered by (chromosome, position) and each chromosome's arc is
## proportional to its diagnostic-SNP COUNT (not bp), so a BDMI interval is
## mapped to the diagnostic SNPs whose positions fall inside it.
##
##   CUTOFF env var selects the bed by its k index (default 14 -> X2=0.00275,
##   the cutoff the colleague plotted).  SWEEP_LEVEL is SNP (only SNP supported).
##
## Run from repo root:  CUTOFF=14 Rscript module_di25/R/bdmi_sorting_circos.R
## =========================================================

suppressMessages(library(data.table))
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()
source("module_di25/R/diem_circos_core.R")          # render_circos_raster()

CUTOFF_K <- as.integer(Sys.getenv("CUTOFF", "14"))
TAU_GRID <- c(0.5, 0.6, 0.7, 0.8)
BEDDIR   <- "data/liftoff_Frufa_DTOL_PR"
BDMI_COL <- "#B4661E"                                # orange, matches colleague's arcs

SORT_PAL <- c("#F4F4F4", "#21918C", "#D3C93B", "#440154")
CLS_CODE <- c(unsorted = 0L, aquilonia = 1L, polyctena = 2L, unresolved = 3L, ambiguous = 3L)

## ---- sorting units, in the sweep's own order -------------------------------
ps  <- readRDS("module_di25/data/di25_sorting_snp.rds")
chr_num <- as.integer(sub("Chr", "", sub(":.*", "", ps$marker)))
pos     <- as.integer(sub(".*:", "", ps$marker))
ord <- order(chr_num, pos); ps <- ps[ord]; chr_num <- chr_num[ord]; pos <- pos[ord]

ok <- ps$differentiated & ps$n_obs > 0 & is.finite(ps$uni_score)
code <- vapply(TAU_GRID, function(tau) {
  v <- integer(nrow(ps))
  cl <- classify_sort(ps$n_aqu[ok], ps$n_pol[ok], ps$n_obs[ok],
                      sort_th = tau, sort_rule = "binom", alpha = 0.05)
  v[which(ok)] <- CLS_CODE[cl]; v
}, integer(nrow(ps)))

## ---- BDMI membership per unit at the chosen cutoff -------------------------
bedf <- list.files(BEDDIR, pattern = sprintf("^bdmi_candidates\\.cutoff_%d_.*\\.bed$", CUTOFF_K))
stopifnot(length(bedf) == 1L)
x2 <- as.numeric(sub("^0", "0.", sub(".*cutoff_\\d+_(\\d+)\\..*", "\\1", bedf)))
bed <- fread(file.path(BEDDIR, bedf), header = FALSE, col.names = c("chr", "start", "end"))
bed[, chrn := as.integer(sub("chromosome_", "", chr))]

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
inb <- logical(nrow(ps))
for (cc in unique(bed$chrn)) {
  iv  <- merge_iv(bed[chrn == cc]$start, bed[chrn == cc]$end)
  idx <- chr_num == cc
  inb[idx] <- in_intervals(pos[idx], iv)
}
bdmi_code <- matrix(as.integer(inb), ncol = 1L)      # units x 1 ring
n_in <- sum(inb); n_in_sorted <- sum(inb & code[, 2] != 0L)   # tau=0.6 column
message(sprintf("[circos] cutoff %d (X2=%.5f): %d regions, %d/%d in-region units sorted (%.0f%%)",
                CUTOFF_K, x2, nrow(bed), n_in_sorted, n_in, 100 * n_in_sorted / n_in))

## ---- render: tau-sweep base (inner) + BDMI ring (outer) --------------------
OUTPNG <- sprintf("module_di25/Figures/bdmi_sorting_circos_cutoff%d.png", CUTOFF_K)
ttl <- sprintf("BDMI candidate regions vs ancestry sorting  (X^2 cutoff = %.5f, %d regions)",
               x2, nrow(bed))
png(OUTPNG, width = 2600, height = 2600, res = 300, type = "cairo")
op <- par(mar = c(0, 0, 2, 0), xpd = NA); on.exit({ par(op); dev.off() }, add = TRUE)
plot.new(); plot.window(xlim = c(-1, 1), ylim = c(-1, 1), asp = 1)

render_circos_raster(
  code, chr_num, palette = SORT_PAL, title = "", r_in = 0.30, r_out = 0.85,
  ring_labels = sprintf("tau=%.1f", TAU_GRID), draw_chr_labels = FALSE,
  new_device = FALSE, add = TRUE, bg_col = "#FFFFFF")

render_circos_raster(
  bdmi_code, chr_num, palette = c(NA, BDMI_COL), r_in = 0.865, r_out = 0.905,
  ring_sep = FALSE, draw_chr_labels = TRUE, chr_label_cex = 0.5,
  add = TRUE, bg_col = NA)

title(main = ttl, cex.main = 0.9, line = 0.2)
legend("bottomleft", bty = "n", cex = 0.72, border = NA, inset = c(0.02, 0.02),
       legend = c("toward F. aquilonia", "toward F. polyctena", "direction unresolved",
                  "unsorted", "BDMI candidate region"),
       fill = c(SORT_PAL[2], SORT_PAL[3], SORT_PAL[4], SORT_PAL[1], BDMI_COL))
message("[circos] wrote ", OUTPNG)
