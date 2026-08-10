## =========================================================
## module_di25 (high-DI analyses) -- DIEM circos on the from-scratch DI25 clustering
## =========================================================
## LD-reduced DIEM circos built on the high-DI-ONLY clustering produced by
## di25_ld_clustering.R (5 cM cap), plotting ALL units:
##   * clusters with > 2 markers  -> eMLG consensus dosage
##   * clusters with 1-2 markers  -> representative SNP genotype
##
## Unlike diem_circos_eMLG.R (which used the canonical all-marker eMLGs), every
## unit here summarises high-DI variation only.
##
## Consensus is computed over the COMBINED 165 hybrids + 30 parents in one pass
## (one polarization per unit), so hybrid and parent sides share a sign
## convention; each unit is then oriented to F. aquilonia via the parent rows.
##
## Run from the repo root:  Rscript module_di25/R/diem_circos_di25_eMLG.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")                 # consensus_dosage(), expected_gt_dosage()
source("module_di25/R/diem_circos_core.R")

CM_STAMP <- "cM5"
CLUST    <- sprintf("module_di25/data/di25_clustering_%s.rds", CM_STAMP)
INPUTS   <- "module_di25/data/di25_inputs.rds"
OUTPNG   <- "module_di25/Figures/diem_circos_di25_eMLG_cM5.png"
DI_LABEL <- -25

## ---- load clustering + prepared 012 genotypes ---------------------------
message("[di25-circos] loading ", CLUST)
res <- readRDS(CLUST); g <- res$groups
inp <- readRDS(INPUTS)
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)             # 195 individuals x markers (012, marker colnames)
faqu <- grep("^Faqu", rownames(GTs_all)); fpol <- grep("^Fpol", rownames(GTs_all))
message(sprintf("  %d units | %d individuals (%d hybrids + %d parents)",
                nrow(g), nrow(GTs_all), nrow(inp$GTs_hyb), nrow(inp$GTs_par)))

## ---- one dosage vector per unit (eMLG consensus or representative SNP) ---
is_emlg <- g$n_loci > 2
message(sprintf("  %d eMLG units (>2 markers) + %d representative-SNP units",
                sum(is_emlg), sum(!is_emlg)))
D <- vapply(seq_len(nrow(g)), function(i) {
  if (is_emlg[i]) consensus_dosage(GTs_all, g$members[[i]])   # 0..2, polarized
  else            GTs_all[, g$representative[i]]               # 0/1/2 representative SNP
}, numeric(nrow(GTs_all)))                                     # individuals x units
D <- t(D)                                                      # units x individuals

## ---- orient each unit to F. aquilonia via the parent side ---------------
maqu <- rowMeans(D[, faqu, drop = FALSE], na.rm = TRUE)
mpol <- rowMeans(D[, fpol, drop = FALSE], na.rm = TRUE)
flip <- which(maqu < mpol)
D[flip, ] <- 2 - D[flip, ]                                     # 2 = aquilonia everywhere

## ---- recode to shared DIEM scheme (1 aqu, 2 het, 3 pol; 0 missing) -------
code <- 3L - round(D)                                          # 2->1 aqu, 1->2 het, 0->3 pol
code[is.na(code)] <- 0L
storage.mode(code) <- "integer"

## ---- position by representative marker; order units & individuals -------
rep_chr <- as.integer(sub("Chr", "", sub(":.*", "", g$representative)))
rep_pos <- as.integer(sub(".*:", "", g$representative))
ord_u   <- order(rep_chr, rep_pos)
code    <- code[ord_u, ]; chr_num <- rep_chr[ord_u]

hi <- apply(code, 2, function(x) { x <- x[x > 0]; if (!length(x)) NA_real_ else mean((x - 1) / 2) })
code <- code[, order(hi)]                                      # inner = most aquilonia

## ---- render -------------------------------------------------------------
n_miss <- sum(code == 0L)
ttl <- sprintf("LD-reduced (DI = %d, %s): %s units [%s eMLG + %s rep-SNP], %s missing (%.1f%%)",
               DI_LABEL, sub("cM", "", CM_STAMP) |> (\(x) paste0(x, " cM"))(),
               format(nrow(code), big.mark = ","),
               format(sum(is_emlg), big.mark = ","), format(sum(!is_emlg), big.mark = ","),
               format(n_miss, big.mark = ","), 100 * n_miss / length(code))

message("[di25-circos] rendering -> ", OUTPNG)
render_diem_circos(code, chr_num, OUTPNG, title = ttl, cex_main = 0.8)
message("[di25-circos] done")
