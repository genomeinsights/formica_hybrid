## =========================================================
## module_di25 (high-DI analyses) -- LD-reduced DIEM circos + sort-class outer track
## =========================================================
## The 5 cM LD-reduced DIEM circos (per-eMLG genotype rings), with an OUTER
## annotation ring giving each unit's ancestry-sorting class at tau = 0.6
## (phi = 0.85). Genotype rings and the sort track share the SAME unit order
## (representative chr/pos), so the outer class aligns radially with the unit
## below it.
##
## Run from the repo root:  Rscript module_di25/R/diem_circos_di25_sortclass.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()
source("module_di25/R/diem_circos_core.R")

TAU_TRACK <- 0.6
OUTPNG    <- "module_di25/Figures/diem_circos_di25_eMLG_sortclass_cM5.png"

## sort-class palette (matches the tau-sweep figure)
SORT_PAL <- c("#F4F4F4", "#21918C", "#D3C93B", "#440154")   # unsorted/aqu/pol/unresolved
CLS_CODE <- c(unsorted = 0L, aquilonia = 1L, polyctena = 2L, unresolved = 3L, ambiguous = 3L)

## ---- inputs (clustering, prepared genotypes, sorting counts) -------------
g   <- readRDS("module_di25/data/di25_clustering_cM5.rds")$groups
inp <- readRDS("module_di25/data/di25_inputs.rds")
ps  <- readRDS("module_di25/data/di25_sorting_emlg.rds")
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)
faqu <- grep("^Faqu", rownames(GTs_all)); fpol <- grep("^Fpol", rownames(GTs_all))

## shared unit order: representative chr/pos
rep_chr <- as.integer(sub("Chr", "", sub(":.*", "", g$representative)))
rep_pos <- as.integer(sub(".*:", "", g$representative))
ord_u   <- order(rep_chr, rep_pos)
chr_num <- rep_chr[ord_u]

## ---- genotype rings (per-eMLG, oriented to aquilonia) -------------------
is_emlg <- g$n_loci > 2
D <- vapply(seq_len(nrow(g)), function(i) {
  if (is_emlg[i]) consensus_dosage(GTs_all, g$members[[i]]) else GTs_all[, g$representative[i]]
}, numeric(nrow(GTs_all)))
D <- t(D)                                                    # units x individuals
flip <- which(rowMeans(D[, faqu, drop = FALSE], na.rm = TRUE) <
              rowMeans(D[, fpol, drop = FALSE], na.rm = TRUE))
D[flip, ] <- 2 - D[flip, ]
code_geno <- 3L - round(D); code_geno[is.na(code_geno)] <- 0L
storage.mode(code_geno) <- "integer"
colnames(code_geno) <- rownames(GTs_all)
code_geno <- code_geno[ord_u, ]                              # unit order
hi <- apply(code_geno, 2, function(x) { x <- x[x > 0]; if (!length(x)) NA_real_ else mean((x - 1) / 2) })
code_geno <- code_geno[, order(hi)]                          # inner ring = most aquilonia

## ---- sort-class outer track (same unit order) --------------------------
psm <- ps[match(g$group_id, group_id)]                       # align ps to g order
ok  <- psm$differentiated & psm$n_obs > 0 & is.finite(psm$uni_score)
track <- integer(nrow(g))
cl <- classify_sort(psm$n_aqu[ok], psm$n_pol[ok], psm$n_obs[ok],
                    sort_th = TAU_TRACK, sort_rule = "binom", alpha = 0.05)
track[which(ok)] <- CLS_CODE[cl]
track <- matrix(track[ord_u], ncol = 1)                      # units x 1, unit order
storage.mode(track) <- "integer"

## ---- draw: genotype inside, sort-class ring outside --------------------
message("[sortclass] rendering -> ", OUTPNG)
png(OUTPNG, width = 2600, height = 2600, res = 300)
render_diem_circos(code_geno, chr_num, new_device = FALSE, draw_chr_labels = FALSE,
                   r_in = 0.30, r_out = 0.80, cex_main = 0.85,
                   title = sprintf("LD-reduced DIEM (5 cM) + sort class @ tau = %.1f (phi = 0.85)  |  %s units",
                                   TAU_TRACK, format(nrow(g), big.mark = ",")))
render_circos_raster(track, chr_num, palette = SORT_PAL, add = TRUE, bg_col = NA,
                     r_in = 0.825, r_out = 0.90, ring_sep = FALSE, draw_chr_labels = TRUE,
                     ring_labels = sprintf("sort tau=%.1f", TAU_TRACK),
                     legend_labels = c("toward F. aquilonia", "toward F. polyctena",
                                       "direction unresolved", "unsorted"),
                     legend_cols = SORT_PAL[c(2, 3, 4, 1)])
dev.off()
message("[sortclass] done")
