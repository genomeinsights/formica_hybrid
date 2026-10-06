## module_external_validation / 06: does the lower, flatter-trend LD-decay rate seen in the external
## 21-sample fit (a_pred ~0.855 x original) come from n = 21 + the MAF>0.1 marker set, or is it
## a property of the external samples?
##
## Design: draw 21 hybrids from the ORIGINAL data (one per population, 20 populations, + 1 extra
## random individual), keep SNPs with MAF > 0.1 IN THE SUBSAMPLE (as for the external VCF), refit
## compute_LD_decay with the external settings (n_win_decay = 100, slide = 100, corr, seed set).
## Two marker variants per replicate:
##   "all"     : every SNP with subsample-MAF > 0.1 (denser than the external set)
##   "thinned" : random subset of the same size as the external set (230,392 SNPs), so slide_bp /
##               marker density match (LDscnR docs: a depends on density because slide is in SNPs)
## Compared with: full original (a_pred median 0.00174), external fit (a 0.00154, a_pred 0.00148,
## b 0.214, c 0.896).
##
## Run from the repo root:  Rscript module_external_validation/R/06_subsample_refit_original.R [n_rep]
## Writes: results/subsample_refit_original.rds, results/subsample_refit_original.tsv, .txt

suppressPackageStartupMessages({ library(data.table); library(SNPRelate); library(parallel); library(ggplot2) })
devtools::load_all("~/gitlab/LDscnR/", quiet = TRUE)
OUT <- "module_external_validation"
NREP <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) as.integer(a[1]) else 5L }
N_EXT <- 230392L; SLIDE <- 100L; CORES <- 8L; MAFMIN <- 0.1

e <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e)
G <- e$GTs_hybrids_005; map <- as.data.table(e$map_hyb_005)[, .(marker, Chr, Pos)]; sd <- as.data.table(e$sample_data)
stopifnot(identical(colnames(G), map$marker), all(rownames(G) %in% sd$Sample_ID))
full <- as.data.table(readRDS("module0_ld_pruning/data/ld_decay_DIEM_100w.rds")$decay_sum)
ext  <- as.data.table(readRDS(file.path(OUT, "data/ld_decay_ext_maf10_slide100.rds"))$decay_sum)

summ <- function(x, rep, variant, nm, ids) {
  d <- as.data.table(x$decay_sum)
  cmp <- merge(d[, .(Chr, a)], full[, .(Chr, a_full = a)], by = "Chr")
  data.table(rep = rep, variant = variant, n_snps = nm, background_b = unique(round(d$b, 4)), c_med = median(d$c),
             a_med = median(d$a), a_pred_med = median(d$a_pred), slide_bp_med = median(d$slide_bp),
             a_pred_ratio_to_full = median(d$a_pred) / median(full$a_pred),
             cor_a_vs_full = cor(cmp$a, cmp$a_full), spearman_a_vs_full = cor(cmp$a, cmp$a_full, method = "spearman"))
}
res <- list()
for (rep in seq_len(NREP)) {
  set.seed(1000 + rep)
  one <- sd[, .SD[sample(.N, 1)], by = Population]$Sample_ID
  ids <- c(one, sample(setdiff(sd$Sample_ID, one), 21L - length(one)))
  g <- G[ids, , drop = FALSE]
  maf <- colMeans(g, na.rm = TRUE) / 2; maf <- pmin(maf, 1 - maf)
  keep <- which(maf > MAFMIN)
  for (variant in c("all", "thinned")) {
    idx <- if (variant == "all") keep else sort(sample(keep, min(N_EXT, length(keep))))
    gf <- file.path(tempdir(), sprintf("sub_%s_%d.gds", variant, rep))
    gds <- create_gds_from_geno(g[, idx], map[idx], gf)
    t0 <- Sys.time()
    x <- compute_LD_decay(gds, n_win_decay = 100, slide = SLIDE, ld_method = "corr", keep_el = FALSE,
                          cores = CORES, min_maf_decay = MAFMIN, seed = rep)
    snpgdsClose(gds); file.remove(gf)
    r <- summ(x, rep, variant, length(idx), ids); res[[length(res) + 1]] <- r
    message(sprintf("rep %d %s: n_snps=%d  a_pred_ratio=%.3f  b=%.3f  (%.1f min)", rep, variant, length(idx), r$a_pred_ratio_to_full, r$background_b,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
R <- rbindlist(res); saveRDS(R, file.path(OUT, "results/subsample_refit_original.rds")); fwrite(R, file.path(OUT, "results/subsample_refit_original.tsv"), sep = "\t")
ext_row <- data.table(rep = NA, variant = "EXTERNAL", n_snps = N_EXT, background_b = round(unique(ext$b), 4), c_med = median(ext$c), a_med = median(ext$a),
                      a_pred_med = median(ext$a_pred), slide_bp_med = median(ext$slide_bp), a_pred_ratio_to_full = median(ext$a_pred) / median(full$a_pred))
sink(file.path(OUT, "results/subsample_refit_original.txt"), split = TRUE)
cat("Full original: a median", signif(median(full$a), 3), " a_pred median", signif(median(full$a_pred), 3), " b", round(unique(full$b), 3), " c", round(median(full$c), 3), " slide_bp(100) n/a\n\n")
print(rbind(R, ext_row, fill = TRUE), digits = 4)
cat("\nMean over replicates by variant:\n"); print(R[, lapply(.SD, mean), by = variant, .SDcols = c("n_snps", "background_b", "c_med", "a_med", "a_pred_med", "slide_bp_med", "a_pred_ratio_to_full", "cor_a_vs_full")], digits = 4)
sink()
