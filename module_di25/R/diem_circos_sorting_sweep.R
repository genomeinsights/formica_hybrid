## =========================================================
## module_di25 (high-DI analyses) -- ancestry-sorting tau-sweep circos
## =========================================================
## The tau sweep drawn AROUND the genome, in the same circular chromosome layout
## as the DIEM circos -- but the concentric rings are the four sorting thresholds
## tau = 0.5, 0.6, 0.7, 0.8 (inner -> outer), and each unit is coloured by its
## sort class at that tau. phi = 0.85 fixed. No DIEM genotype data here.
##
##   colour: toward F. aquilonia (purple) / toward F. polyctena (yellow) /
##           direction unresolved (teal) / unsorted (near-white)
##
## As tau increases outward, colour thins out (fewer units sorted); direction
## colour shows where each ancestry prevails and how the balance shifts with tau.
##
## LEVEL = "SNP" (51,612 markers, mirrors the flagship DIEM circos) or
##         "eMLG" (11,052 LD-reduced units).
##
## Run from the repo root:  Rscript module_di25/R/diem_circos_sorting_sweep.R
## =========================================================

suppressMessages(library(data.table))
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()
source("module_di25/R/diem_circos_core.R")  # render_circos_raster()

LEVEL    <- Sys.getenv("SWEEP_LEVEL", "SNP")         # "SNP" or "eMLG" (env-overridable)
TAU_GRID <- c(0.5, 0.6, 0.7, 0.8)
OUTPNG   <- sprintf("module_di25/Figures/diem_circos_sorting_sweep_%s.png", tolower(LEVEL))

## sort-class palette (code 0 unsorted .. 3 unresolved)
SORT_PAL <- c("#F4F4F4",   # 0 unsorted / not differentiated
              "#21918C",   # 1 toward F. aquilonia (purple)
              "#D3C93B",   # 2 toward F. polyctena (yellow)
              "#440154")   # 3 direction unresolved (teal)
CLS_CODE <- c(unsorted = 0L, aquilonia = 1L, polyctena = 2L, unresolved = 3L, ambiguous = 3L)

## ---- load the per-unit sorting counts + genomic positions ---------------
if (LEVEL == "SNP") {
  ps <- readRDS("module_di25/data/di25_sorting_snp.rds")
  chr_num <- as.integer(sub("Chr", "", sub(":.*", "", ps$marker)))
  pos     <- as.integer(sub(".*:", "", ps$marker))
} else {
  ps <- readRDS("module_di25/data/di25_sorting_emlg.rds")
  g  <- readRDS("module_di25/data/di25_clustering_cM5.rds")$groups
  rep_mk  <- g$representative[match(ps$group_id, g$group_id)]
  chr_num <- as.integer(sub("Chr", "", sub(":.*", "", rep_mk)))
  pos     <- as.integer(sub(".*:", "", rep_mk))
}

## order units by chromosome then position
ord <- order(chr_num, pos); ps <- ps[ord]; chr_num <- chr_num[ord]

## ---- classify every unit at each tau -> code matrix (units x rings) ------
ok <- ps$differentiated & ps$n_obs > 0 & is.finite(ps$uni_score)
code <- vapply(TAU_GRID, function(tau) {
  v <- integer(nrow(ps))                                   # 0 = unsorted / not differentiated
  cl <- classify_sort(ps$n_aqu[ok], ps$n_pol[ok], ps$n_obs[ok],
                      sort_th = tau, sort_rule = "binom", alpha = 0.05)
  v[which(ok)] <- CLS_CODE[cl]
  v
}, integer(nrow(ps)))                                      # units x length(TAU_GRID)

## ---- render -------------------------------------------------------------
n_lab <- if (LEVEL == "SNP") "diagnostic SNPs" else "LD-reduced units"
ttl <- sprintf("Ancestry-sorting tau sweep (phi = 0.85) -- %s %s",
               format(nrow(ps), big.mark = ","), n_lab)
message("[sorting-sweep] rendering -> ", OUTPNG)
render_circos_raster(
  code, chr_num, palette = SORT_PAL, outpng = OUTPNG, title = ttl,
  ring_labels = sprintf("tau=%.1f", TAU_GRID),
  legend_labels = c("toward F. aquilonia", "toward F. polyctena",
                    "direction unresolved", "unsorted"),
  legend_cols = SORT_PAL[c(2, 3, 4, 1)]
)
message("[sorting-sweep] done (inner->outer rings: tau ",
        paste(TAU_GRID, collapse = ", "), ")")
