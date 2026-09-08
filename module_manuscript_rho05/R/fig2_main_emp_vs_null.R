## =========================================================
## module_manuscript_rho05 -- Fig 2 (main text): observed vs neutral-null
## ancestry sorting, assembled composite (a) empirical circos, (b) simulated
## (SLiM) circos, (c) % sorted vs 1000-replicate null
## =========================================================
## This is the missing assembly script identified when mapping manuscript
## figures to producing code: panels (a)/(b) and (c) existed as two separate
## images with no combining script. Both are SNP-level and therefore
## UNAFFECTED by the min_r2 choice (module0/di25 clustering never enters
## here) -- this figure is the same under min_r2=0.2 and min_r2_rho=0.5.
##
## Panel sources (logic reused, not duplicated):
##   (a)/(b) module_di25/R_legacy/diem_circos_compare_emp_vs_sim.R
##            -- empirical per-SNP DIEM circos vs ONE SLiM bootstrap replicate
##               (data/diem_outs_demo/diem_boot<REP>_output.bed)
##   (c)      module_di25/R/di25_sorting_emp_vs_sim.R's panel A
##            -- % diagnostic SNPs sorted, observed vs the 1000-replicate
##               SLiM-null distribution, reusing its saved
##               module_di25/data/di25_sorting_emp_vs_sim.rds (no recompute)
##
## Run from the repo root:  Rscript module_manuscript_rho05/R/fig2_main_emp_vs_null.R [rep]
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(magick) })
devtools::load_all("~/gitlab/LDscnR/")
source("module_di25/R/diem_circos_core.R")

args   <- commandArgs(trailingOnly = TRUE)
REP    <- if (length(args) >= 1) as.integer(args[1]) else 1L
TSV    <- "data/species_diagnostic_markers_DI25_20pops.tsv.gz"
BED    <- sprintf("data/diem_outs_demo/diem_boot%d_output.bed", REP)
SIM_RDS <- "module_di25/data/di25_sorting_emp_vs_sim.rds"
FIGDIR <- "module_manuscript_rho05/Figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
TMP_AB <- file.path(FIGDIR, sprintf("_tmp_fig2_ab_rep%d.png", REP))
TMP_C  <- file.path(FIGDIR, "_tmp_fig2_c.png")
OUTPNG <- file.path(FIGDIR, "fig2_main_emp_vs_null.png")

## ---- helper: per-SNP hybrid index -> individual order (inner = most aqu) -----
hybrid_order <- function(gt) {
  hi <- apply(gt, 2, function(x) { x <- x[x > 0]; if (!length(x)) NA_real_ else mean((x - 1) / 2) })
  colnames(gt)[order(hi)]
}

## ============================================================
## panels (a) empirical, (b) simulated -- circos, into one temp PNG
## ============================================================
message("[fig2-main] (a) empirical per-SNP panel")
d       <- fread(TSV)
emp_chr <- as.integer(sub("chromosome_", "", d$chromosome))
emp_gt  <- as.matrix(d[, !c("chromosome", "position"), with = FALSE])
storage.mode(emp_gt) <- "integer"
ord_e   <- order(emp_chr, d$position)
emp_gt  <- emp_gt[ord_e, ]; emp_chr <- emp_chr[ord_e]
emp_gt  <- emp_gt[, hybrid_order(emp_gt)]

message("[fig2-main] (b) simulated per-SNP panel  <- ", BED)
hdr      <- readLines(BED, n = 2)
sim_inds <- strsplit(strsplit(hdr[2], "\t")[[1]][10], "\\|")[[1]]
n_ind    <- length(sim_inds)
sim      <- fread(BED, skip = 2, header = FALSE, sep = "\t",
                  select = c(1, 3, 10), colClasses = list(character = c(1, 10)))
sim_chr  <- as.integer(sub("ch", "", sim$V1))
sim_pos  <- as.integer(sim$V3)
geno     <- sub("^S", "", sim$V10)
S <- matrix(unlist(strsplit(geno, "", fixed = TRUE), use.names = FALSE),
            nrow = length(geno), byrow = TRUE)
stopifnot(ncol(S) == n_ind)
dos <- matrix(NA_integer_, nrow(S), ncol(S))
dos[S == "0"] <- 0L; dos[S == "1"] <- 1L; dos[S == "2"] <- 2L
colnames(dos) <- sim_inds
aq  <- grep("^aq_",  sim_inds); pol <- grep("^pol_", sim_inds)
m_aq  <- rowMeans(dos[, aq,  drop = FALSE], na.rm = TRUE)
m_pol <- rowMeans(dos[, pol, drop = FALSE], na.rm = TRUE)
flip  <- which(m_aq > m_pol)
dos[flip, ] <- 2L - dos[flip, ]
sim_gt <- dos + 1L
sim_gt[is.na(sim_gt)] <- 0L
storage.mode(sim_gt) <- "integer"
ord_s  <- order(sim_chr, sim_pos)
sim_gt <- sim_gt[ord_s, ]; sim_chr <- sim_chr[ord_s]
sim_gt <- sim_gt[, hybrid_order(sim_gt)]

png(TMP_AB, width = 5400, height = 3000, res = 300)
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
render_diem_circos(emp_gt, emp_chr, new_device = FALSE,
                   title = sprintf("a   Empirical: %s diagnostic SNPs x %d individuals",
                                   format(nrow(emp_gt), big.mark = ","), ncol(emp_gt)),
                   cex_main = 1.0, chr_label_cex = 0.85)
render_diem_circos(sim_gt, sim_chr, new_device = FALSE,
                   title = sprintf("b   Simulated (SLiM, rep %d): %s SNPs x %d individuals",
                                   REP, format(nrow(sim_gt), big.mark = ","), ncol(sim_gt)),
                   cex_main = 1.0, chr_label_cex = 0.85)
mtext("DIEM ancestry (teal = F. aquilonia, dark = het, yellow = F. polyctena)  |  rings inner to outer: most aquilonia to most polyctena",
      outer = TRUE, cex = 0.8, line = 0.3)
dev.off()

## ============================================================
## panel (c) -- % sorted vs null, reusing the saved comparison (no recompute)
## ============================================================
message("[fig2-main] (c) % sorted vs null panel  <- ", SIM_RDS)
obj <- readRDS(SIM_RDS)
sim_tally <- obj$sim; emp_tally <- obj$emp
TAU_GRID <- c(0.5, 0.6, 0.7, 0.8)
pC <- ggplot(sim_tally, aes(factor(tau), pct_sorted)) +
  geom_jitter(width = 0.18, height = 0, colour = "#1b9e77", alpha = 0.25, size = 0.8) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = NA, colour = "#0b3d2e", linewidth = 0.8) +
  geom_point(data = emp_tally, aes(factor(tau), pct_sorted), colour = "#d95f02", size = 5, shape = 18) +
  labs(x = expression("sorting threshold  " * tau), y = "% diagnostic SNPs sorted",
       title = sprintf("c   Observed vs %d-replicate SLiM null", obj$n_rep)) +
  theme_bw(base_size = 15) + theme(panel.grid.minor = element_blank())
ggsave(TMP_C, pC, width = 5.5, height = 5.5, dpi = 300)

## ============================================================
## stitch a/b (wide) above c (narrower, centred) into one composite
## ============================================================
message("[fig2-main] stitching -> ", OUTPNG)
img_ab <- image_read(TMP_AB)
img_c  <- image_read(TMP_C)
w_ab <- image_info(img_ab)$width
img_c_rs <- image_resize(img_c, paste0(round(w_ab * 0.42), "x"))
composite <- image_append(c(img_ab, image_border(img_c_rs, "white", paste0(round((w_ab - image_info(img_c_rs)$width)/2), "x20"))),
                          stack = TRUE)
image_write(composite, OUTPNG)
file.remove(TMP_AB, TMP_C)
message("[fig2-main] done")
