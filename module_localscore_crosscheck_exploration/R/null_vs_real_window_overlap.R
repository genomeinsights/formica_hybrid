## =========================================================
## module_localscore_crosscheck -- do the null draws' significant windows
## land on the SAME genomic regions as the real covariates' windows?
## =========================================================
## Prompted by a visual impression that the full-SNP null demonstration
## figures (local_score_null_fullsnp_manhattan.png) "look similar" to the
## real observed-data figure (local_score_manhattan.png). Checked directly
## by exact physical overlap between window sets, plus a position-
## randomization test (shuffle each null window to a random position on
## its own chromosome, preserving window size and chromosome; 1,000 sims)
## to ask whether any observed overlap exceeds what similar genome
## coverage would produce by chance alone.
##
## Reads : module_localscore_crosscheck/data/localscore_{PC1,PC2,bio_winter,mitoC2}.rds
##         module_localscore_crosscheck/data/localscore_null_fullsnp_continuous_draw74.rds
##         module_localscore_crosscheck/data/localscore_null_fullsnp_mitoC2_draw212.rds
## Writes: module_localscore_crosscheck/data/null_vs_real_window_overlap.tsv
##
## Run from the repo root:
##   Rscript module_localscore_crosscheck/R/null_vs_real_window_overlap.R
## =========================================================

suppressMessages(library(data.table))
DATA <- "module_localscore_crosscheck/data"
load("data/hybrids_only_maf005.Rdata")
chr_lens <- map_hyb_005[, .(len = max(Pos)), by = Chr]
genome_len <- sum(chr_lens$len)

## ---- per-SNP LD context: ld_w_095 (LD weight) and Stage-1/Stage-2 cluster
## size, used below to test whether large/overlapping windows specifically
## sit in high-LD, low-recombination genomic regions -----------------------
map_snp <- readRDS("module0_ld_pruning/data/pruned_stage1.rds")$map_snp
s2ref <- readRDS("module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05.rds")
g2ref <- as.data.table(s2ref$groups)
marker2s2ref <- g2ref[, .(marker = unlist(members)), by = .(s2_group = group_id, s2_nloci = n_loci)]
setkey(marker2s2ref, marker); setkey(map_snp, marker)
map_snp[marker2s2ref, on = "marker", `:=`(s2_group = i.s2_group, s2_nloci = i.s2_nloci)]

annotate_windows <- function(win_dt) {
  win_dt <- copy(win_dt)
  win_dt[, `:=`(mean_ldw095 = NA_real_, max_s2_nloci = NA_real_, n_snp_in_win = NA_integer_)]
  for (i in seq_len(nrow(win_dt))) {
    sub <- map_snp[Chr == win_dt$chr[i] & Pos >= win_dt$beg[i] & Pos <= win_dt$end[i]]
    if (nrow(sub) == 0) next
    win_dt[i, `:=`(mean_ldw095 = mean(sub$ld_w_095, na.rm = TRUE),
                   max_s2_nloci = max(sub$s2_nloci, na.rm = TRUE), n_snp_in_win = nrow(sub))]
  }
  win_dt
}

overlap_count <- function(real, null) {
  if (nrow(real) == 0 || nrow(null) == 0) return(0L)
  hits <- 0L
  for (i in seq_len(nrow(null))) {
    ov <- real[chr == null$chr[i] & beg <= null$end[i] & end >= null$beg[i]]
    if (nrow(ov) > 0) hits <- hits + 1L
  }
  hits
}

randomization_test <- function(real, null, n_sim = 1000L, seed = 1L) {
  set.seed(seed)
  counts <- integer(n_sim)
  for (s in seq_len(n_sim)) {
    shuf <- copy(null)
    for (i in seq_len(nrow(shuf))) {
      L <- shuf$end[i] - shuf$beg[i]
      maxlen <- chr_lens[Chr == shuf$chr[i], len]
      if (maxlen > L) {
        newbeg <- sample.int(maxlen - L, 1)
        shuf[i, `:=`(beg = newbeg, end = newbeg + L)]
      }
    }
    counts[s] <- overlap_count(real, shuf)
  }
  counts
}

run_one <- function(real_tag, real_rds, null_tag, null_rds, n_sim = 1000L) {
  real <- as.data.table(readRDS(real_rds)$significant.windows)
  null <- as.data.table(readRDS(null_rds)$significant.windows)
  obs <- overlap_count(real, null)
  real_bp <- sum(real$end - real$beg); null_bp <- sum(null$end - null$beg)
  message(sprintf("[%s vs %s] %d real windows (%.2f%% genome), %d null windows (%.2f%% genome), observed overlap=%d",
                  real_tag, null_tag, nrow(real), 100 * real_bp / genome_len,
                  nrow(null), 100 * null_bp / genome_len, obs))
  if (nrow(null) < 3) {
    message("  (too few null windows for a meaningful randomization test)")
    return(data.table(real = real_tag, null = null_tag, n_real_windows = nrow(real), n_null_windows = nrow(null),
                      real_pct_genome = 100 * real_bp / genome_len, null_pct_genome = 100 * null_bp / genome_len,
                      observed_overlap = obs, null_of_null_mean = NA_real_, null_of_null_lo = NA_real_,
                      null_of_null_hi = NA_real_, emp_p = NA_real_))
  }
  sim_counts <- randomization_test(real, null, n_sim = n_sim)
  emp_p <- mean(sim_counts >= obs)
  message(sprintf("  position-randomized null-of-null: mean=%.2f 95%%=[%.0f,%.0f] max=%.0f | empirical p(>=obs)=%.4f",
                  mean(sim_counts), quantile(sim_counts, 0.025), quantile(sim_counts, 0.975), max(sim_counts), emp_p))
  data.table(real = real_tag, null = null_tag, n_real_windows = nrow(real), n_null_windows = nrow(null),
            real_pct_genome = 100 * real_bp / genome_len, null_pct_genome = 100 * null_bp / genome_len,
            observed_overlap = obs, null_of_null_mean = mean(sim_counts),
            null_of_null_lo = quantile(sim_counts, 0.025, names = FALSE),
            null_of_null_hi = quantile(sim_counts, 0.975, names = FALSE), emp_p = emp_p)
}

res <- rbindlist(list(
  run_one("PC1", file.path(DATA, "localscore_PC1.rds"), "null_draw74", file.path(DATA, "localscore_null_fullsnp_continuous_draw74.rds")),
  run_one("PC2", file.path(DATA, "localscore_PC2.rds"), "null_draw74", file.path(DATA, "localscore_null_fullsnp_continuous_draw74.rds")),
  run_one("bio_winter", file.path(DATA, "localscore_bio_winter.rds"), "null_draw74", file.path(DATA, "localscore_null_fullsnp_continuous_draw74.rds")),
  run_one("mitoC2", file.path(DATA, "localscore_mitoC2.rds"), "null_draw212", file.path(DATA, "localscore_null_fullsnp_mitoC2_draw212.rds"))
))
fwrite(res, file.path(DATA, "null_vs_real_window_overlap.tsv"), sep = "\t")
cat("\n"); print(res)
cat("\n[overlap] wrote module_localscore_crosscheck/data/null_vs_real_window_overlap.tsv\n")

## ========================================================================
## Why do null and real windows overlap more than a size-matched random
## shuffle predicts? Test whether window SIZE itself (and hence overlap
## probability) is driven by local LD/recombination context (ld_w_095,
## Stage-2 cluster size) rather than by covariate-specific signal -- true
## for BOTH real and null windows alike, this would explain excess overlap
## as a shared LD-architecture artifact of the method, not of either
## specific covariate.
## ========================================================================
real_mito <- annotate_windows(as.data.table(readRDS(file.path(DATA, "localscore_mitoC2.rds"))$significant.windows))
null_mito <- annotate_windows(as.data.table(readRDS(file.path(DATA, "localscore_null_fullsnp_mitoC2_draw212.rds"))$significant.windows))
null_mito[, is_overlap := FALSE]
for (i in seq_len(nrow(null_mito))) {
  ov <- real_mito[chr == null_mito$chr[i] & beg <= null_mito$end[i] & end >= null_mito$beg[i]]
  if (nrow(ov) > 0) null_mito[i, is_overlap := TRUE]
}

allw <- rbind(real_mito, null_mito, fill = TRUE)
allw[, size_bp := end - beg]
cor_size_ldw <- cor(log10(allw$size_bp), allw$mean_ldw095, use = "complete.obs")
cor_size_s2  <- cor(log10(allw$size_bp), allw$max_s2_nloci, use = "complete.obs")
cat(sprintf("\n=== window size vs. LD context (mitoC2 real+null windows pooled, n=%d) ===\n", nrow(allw)))
cat(sprintf("cor(log10(window size bp), mean ld_w_095) = %.3f\n", cor_size_ldw))
cat(sprintf("cor(log10(window size bp), max Stage-2 cluster n_loci) = %.3f\n", cor_size_s2))

null_mito[, size_kb := (end - beg) / 1000]
by_overlap <- null_mito[, .(n = .N, mean_size_kb = mean(size_kb), mean_ldw095 = mean(mean_ldw095, na.rm = TRUE),
                            mean_max_s2_nloci = mean(max_s2_nloci, na.rm = TRUE)), by = is_overlap]
cat("\n=== null mitoC2 windows: overlapping a real window vs not ===\n")
print(by_overlap)

fwrite(allw, file.path(DATA, "null_vs_real_window_ld_context.tsv"), sep = "\t")
fwrite(by_overlap, file.path(DATA, "null_overlap_vs_ld_context_summary.tsv"), sep = "\t")
cat("\n[overlap] wrote null_vs_real_window_ld_context.tsv and null_overlap_vs_ld_context_summary.tsv\n")
