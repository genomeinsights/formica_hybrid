## =========================================================
## module_di25 (high-DI analyses) -- ancestry-sorting tau-sweep circos (SLIDE)
## =========================================================
## Slide-optimised twin of diem_circos_sorting_sweep.R: the same four-ring tau
## sweep (tau = 0.5, 0.6, 0.7, 0.8 inner->outer; phi = 0.85 fixed), but larger
## fonts, the cairo device (so the greek-tau / em-dash glyphs render), short
## "tau = x" ring labels, and higher resolution.
##
## Each ring classifies all 51,612 diagnostic SNPs at that tau: a SNP is
## coloured by sort direction (aquilonia teal / polyctena yellow / direction
## unresolved purple) or left grey (unsorted). Angle = genomic position;
## radius = tau ring. As tau rises outward, isolated sorted SNPs drop out and
## only the large robustly-sorted low-recombination blocks persist
## (chr26/5/6 -> polyctena, chr16/15 -> aquilonia). NB the within-ring dashed
## look of an isolated SNP is undersampling (each SNP is sub-pixel among 51k);
## the white gaps between rings are the deliberate ring separators.
##
## Run from the repo root:  Rscript module_di25/R/diem_circos_sorting_sweep_slide.R
## =========================================================

suppressMessages(library(data.table))
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()
source("module_di25/R/diem_circos_core.R")          # render_circos_raster()

TAU_GRID <- c(0.5, 0.6, 0.7, 0.8)
OUTPNG   <- "module_di25/Figures/diem_circos_sorting_sweep_slide.png"
TAU <- intToUtf8(0x3c4); EMD <- intToUtf8(0x2014)   # tau, em-dash (keeps the source ASCII)

## sort-class palette (code 0 unsorted .. 3 unresolved)
SORT_PAL <- c("#F4F4F4",   # 0 unsorted / not differentiated
              "#21918C",   # 1 toward F. aquilonia
              "#D3C93B",   # 2 toward F. polyctena
              "#440154")   # 3 direction unresolved
CLS_CODE <- c(unsorted = 0L, aquilonia = 1L, polyctena = 2L, unresolved = 3L, ambiguous = 3L)

## ---- per-SNP sorting counts + genomic positions -------------------------
ps <- as.data.table(readRDS("module_di25/data/di25_sorting_snp.rds"))
chr_num <- as.integer(sub("Chr", "", sub(":.*", "", ps$marker)))
pos     <- as.integer(sub(".*:", "", ps$marker))
ord <- order(chr_num, pos); ps <- ps[ord]; chr_num <- chr_num[ord]

## ---- classify every SNP at each tau -> code matrix (SNPs x rings) -------
ok <- ps$differentiated & ps$n_obs > 0 & is.finite(ps$uni_score)
code <- vapply(TAU_GRID, function(tau) {
  v <- integer(nrow(ps))
  cl <- classify_sort(ps$n_aqu[ok], ps$n_pol[ok], ps$n_obs[ok],
                      sort_th = tau, sort_rule = "binom", alpha = 0.05)
  v[which(ok)] <- CLS_CODE[cl]
  v
}, integer(nrow(ps)))

## ---- render (slide fonts) ----------------------------------------------
message("[sorting-sweep-slide] rendering -> ", OUTPNG)
render_circos_raster(
  code, chr_num, palette = SORT_PAL, outpng = OUTPNG, npx = 3800, res = 200,
  title = sprintf("Ancestry-sorting threshold sweep %s %s diagnostic SNPs",
                  EMD, format(nrow(ps), big.mark = ",")),
  cex_main = 1.8, main_line = 1.0,
  ring_labels = sprintf("%s = %.1f", TAU, TAU_GRID), ring_label_cex = 1.35,
  chr_label_cex = 1.5,
  legend_labels = c("toward F. aquilonia", "toward F. polyctena",
                    "direction unresolved", "unsorted"),
  legend_cols = SORT_PAL[c(2, 3, 4, 1)], legend_cex = 1.4
)
message("[sorting-sweep-slide] done (inner->outer rings: tau ",
        paste(TAU_GRID, collapse = ", "), ")")
