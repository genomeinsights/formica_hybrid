## =========================================================
## module_di25 (high-DI analyses) -- empirical vs SIMULATED DIEM circos (side by side)
## =========================================================
## One figure, two panels on the SAME 51,612 DI25 diagnostic markers (26 chr, no
## chr23) and the SAME 195-individual design (165 hybrids + 15 aq + 15 pol parents):
##   (a) empirical   -- module_di25 per-SNP DIEM  (data/species_diagnostic_markers_DI25_20pops.tsv.gz)
##   (b) simulated   -- ONE bootstrap replicate    (formica_hybrid/data/diem_outs_demo/diem_boot<N>_output.bed)
##
## The simulated set is a *separate realisation* with its own individual labels
## (hyb_aland_.., aq_sielva_.., pol_karsi_..), so panels cannot share an individual
## order.  Each panel is therefore ordered by its OWN per-SNP hybrid index (inner
## ring = most F. aquilonia, outer = most F. polyctena).  The comparison is
## structural: does the neutral simulation reproduce the sorting bands / genomic
## architecture seen empirically?
##
## Simulated markers carry DIEM's arbitrary per-marker polarity, so each is
## re-oriented here using its own aq/pol parents (aqu hom -> 1, het -> 2, pol hom -> 3),
## matching the canonical coding of the empirical TSV.
##
## Run from the formica_hybrid repo root:
##   Rscript module_di25/R/diem_circos_compare_emp_vs_sim.R [rep]
##   (rep defaults to 1 -> diem_boot1_output.bed)
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")
source("module_di25/R/diem_circos_core.R")

args   <- commandArgs(trailingOnly = TRUE)
REP    <- if (length(args) >= 1) as.integer(args[1]) else 1L
TSV    <- "data/species_diagnostic_markers_DI25_20pops.tsv.gz"
BED    <- sprintf("../formica_hybrid/data/diem_outs_demo/diem_boot%d_output.bed", REP)
if (!file.exists(BED)) BED <- sprintf("data/diem_outs_demo/diem_boot%d_output.bed", REP)
OUTPNG <- sprintf("module_di25/Figures/diem_circos_compare_emp_vs_sim_rep%d.png", REP)

## ---- helper: per-SNP hybrid index -> individual order (inner = most aqu) -----
hybrid_order <- function(gt) {                 # gt: markers x individuals, 1/2/3, 0 miss
  hi <- apply(gt, 2, function(x) { x <- x[x > 0]; if (!length(x)) NA_real_ else mean((x - 1) / 2) })
  colnames(gt)[order(hi)]
}

## ---- (a) empirical per-SNP --------------------------------------------------
message("[emp-vs-sim] empirical per-SNP panel")
d       <- fread(TSV)
emp_chr <- as.integer(sub("chromosome_", "", d$chromosome))
emp_gt  <- as.matrix(d[, !c("chromosome", "position"), with = FALSE])   # 1/2/3, 0=miss
storage.mode(emp_gt) <- "integer"
ord_e   <- order(emp_chr, d$position)
emp_gt  <- emp_gt[ord_e, ]; emp_chr <- emp_chr[ord_e]
emp_gt  <- emp_gt[, hybrid_order(emp_gt)]

## ---- (b) simulated per-SNP (one DIEM bootstrap replicate) -------------------
message("[emp-vs-sim] simulated per-SNP panel  <- ", BED)

## individual names live in the '#Chrom ...' header line, field 10, pipe-separated
hdr      <- readLines(BED, n = 2)
sim_inds <- strsplit(strsplit(hdr[2], "\t")[[1]][10], "\\|")[[1]]
n_ind    <- length(sim_inds)

## data: V1 = chrom (ch#), V3 = end/position, V10 = 'S' + one state char per individual
sim      <- fread(BED, skip = 2, header = FALSE, sep = "\t",
                  select = c(1, 3, 10), colClasses = list(character = c(1, 10)))
sim_chr  <- as.integer(sub("ch", "", sim$V1))
sim_pos  <- as.integer(sim$V3)
geno     <- sub("^S", "", sim$V10)                       # drop the leading 'S'

## explode the state strings into a markers x individuals character matrix
S <- matrix(unlist(strsplit(geno, "", fixed = TRUE), use.names = FALSE),
            nrow = length(geno), byrow = TRUE)
stopifnot(ncol(S) == n_ind)

## DIEM states: '0'/'2' = the two homozygotes, '1' = het, '_' or '.' = missing.
## dosage in {0,1,2}; NA = missing.
dos <- matrix(NA_integer_, nrow(S), ncol(S))
dos[S == "0"] <- 0L; dos[S == "1"] <- 1L; dos[S == "2"] <- 2L
colnames(dos) <- sim_inds

## orient each marker by its parents so low dosage = F. aquilonia
aq  <- grep("^aq_",  sim_inds); pol <- grep("^pol_", sim_inds)
m_aq  <- rowMeans(dos[, aq,  drop = FALSE], na.rm = TRUE)
m_pol <- rowMeans(dos[, pol, drop = FALSE], na.rm = TRUE)
flip  <- which(m_aq > m_pol)                             # aq is the high-dosage homozygote -> flip
dos[flip, ] <- 2L - dos[flip, ]

sim_gt <- dos + 1L                                       # 0/1/2 dosage -> 1/2/3 code
sim_gt[is.na(sim_gt)] <- 0L                              # missing -> 0
storage.mode(sim_gt) <- "integer"

ord_s  <- order(sim_chr, sim_pos)
sim_gt <- sim_gt[ord_s, ]; sim_chr <- sim_chr[ord_s]
sim_gt <- sim_gt[, hybrid_order(sim_gt)]

## ---- draw both panels into one wide canvas ---------------------------------
message("[emp-vs-sim] rendering -> ", OUTPNG)
png(OUTPNG, width = 5400, height = 3000, res = 300)
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
render_diem_circos(emp_gt, emp_chr, new_device = FALSE,
                   title = sprintf("a   Empirical: %s diagnostic SNPs x %d individuals",
                                   format(nrow(emp_gt), big.mark = ","), ncol(emp_gt)),
                   cex_main = 1.0, chr_label_cex = 0.85)
render_diem_circos(sim_gt, sim_chr, new_device = FALSE,
                   title = sprintf("b   Simulated (rep %d): %s SNPs x %d individuals",
                                   REP, format(nrow(sim_gt), big.mark = ","), ncol(sim_gt)),
                   cex_main = 1.0, chr_label_cex = 0.85)
mtext(paste("DIEM ancestry (teal = F. aquilonia, dark = het, yellow = F. polyctena)  |",
            "rings inner to outer: most aquilonia to most polyctena  |  each panel ordered by its own hybrid index"),
      outer = TRUE, cex = 0.8, line = 0.3)
dev.off()
message("[emp-vs-sim] done")
