## =========================================================
## module_localscore_crosscheck -- independent local-score computation on
## the 40 full-SNP null sweep outputs (10 replicates x continuous/mitoC2 x
## structured/unstructured).
## =========================================================
## Fresh implementation (not sourced from the exploratory module). Uses
## BayPass's own compute.local.scores() (Fariello et al. 2017 / Bonhomme
## et al. 2019), same call convention as validated by the exploratory
## module's usage (confirmed by direct inspection: xi=1,
## pval.local.score.thres=0.01, min.maf=0.2, min.nsnp=1e4, MAF filtering is
## internal to compute.local.scores() via snp.pi, no separate pre-filter,
## per-chromosome looping is internal to the function -- one call per run
## covers the whole genome).
##
## Reads (per run, all under RAW_DIR):
##   <run>_summary_betai_reg.out   (continuous: MRK, BF(dB))
##   <run>_summary_contrast.out    (mitoC2: MRK, log10(1/pval))
##   <run>_summary_pi_xtx.out      (both modes: MRK, M_P -- required by
##                                  compute.local.scores() for MAF filtering)
##   fullsnp_chr_pos.rds           (Chr, Pos per SNP, row-order-aligned to
##                                  u_DIEM.geno / all BayPass MRK indices)
## Writes (per run): <run>_localscores.rds (list: res.local.scores,
##   significant.windows), plus a combined manifest/summary table.
##
## Designed to run standalone on mini2 against RAW_DIR=/Volumes/T9/... (no
## need to transfer the ~97MB pi_xtx files back to the laptop). Also runs
## locally if RAW_DIR/BAYPASS_UTILS are edited to local paths.
## =========================================================

suppressMessages(library(data.table))

RAW_DIR       <- Sys.getenv("LOCALSCORE_RAW_DIR", "/Volumes/T9/module_localscore_crosscheck_full_snp_null10/raw")
POS_FILE      <- Sys.getenv("LOCALSCORE_POS_FILE", "fullsnp_chr_pos.rds")
BAYPASS_UTILS <- Sys.getenv("LOCALSCORE_BAYPASS_UTILS", "~/baypass_public/utils/baypass_utils.R")
OUT_DIR       <- Sys.getenv("LOCALSCORE_OUT_DIR", file.path(RAW_DIR, "localscores"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(BAYPASS_UTILS)
pos <- readRDS(POS_FILE)
stopifnot(ncol(pos) == 2, nrow(pos) == 1114423)

NULLTYPES <- c("structured", "unstructured")
DRAWS     <- 1:10

run_continuous <- function(nulltype, draw) {
  tag <- sprintf("%s_continuous_draw%02d", nulltype, draw)
  bf <- fread(file.path(RAW_DIR, sprintf("%s_summary_betai_reg.out", tag)), select = c("MRK", "BF(dB)"))
  pi <- fread(file.path(RAW_DIR, sprintf("%s_summary_pi_xtx.out", tag)), select = c("MRK", "M_P"))
  stopifnot(all(bf$MRK == seq_len(nrow(bf))), all(pi$MRK == seq_len(nrow(pi))),
            nrow(bf) == nrow(pos), nrow(pi) == nrow(pos))
  res <- compute.local.scores(snp.position = pos, snp.pi = pi$M_P, snp.bf = bf$`BF(dB)`,
                               xi = 1, pval.local.score.thres = 0.01, min.maf = 0.2, min.nsnp = 1e4)
  saveRDS(res, file.path(OUT_DIR, sprintf("%s_localscores.rds", tag)))
  list(tag = tag, mode = "continuous", nulltype = nulltype, draw = draw,
       n_windows = if (is.null(res$significant.windows)) 0L else nrow(res$significant.windows))
}

run_mitoC2 <- function(nulltype, draw) {
  tag <- sprintf("%s_mitoC2_draw%02d", nulltype, draw)
  pv <- fread(file.path(RAW_DIR, sprintf("%s_summary_contrast.out", tag)), select = c("MRK", "log10(1/pval)"))
  pi <- fread(file.path(RAW_DIR, sprintf("%s_summary_pi_xtx.out", tag)), select = c("MRK", "M_P"))
  stopifnot(all(pv$MRK == seq_len(nrow(pv))), all(pi$MRK == seq_len(nrow(pi))),
            nrow(pv) == nrow(pos), nrow(pi) == nrow(pos))
  res <- compute.local.scores(snp.position = pos, snp.pi = pi$M_P, snp.pvalue = pv$`log10(1/pval)`,
                               xi = 1, pval.local.score.thres = 0.01, min.maf = 0.2, min.nsnp = 1e4)
  saveRDS(res, file.path(OUT_DIR, sprintf("%s_localscores.rds", tag)))
  list(tag = tag, mode = "mitoC2", nulltype = nulltype, draw = draw,
       n_windows = if (is.null(res$significant.windows)) 0L else nrow(res$significant.windows))
}

summary_rows <- list()
for (nt in NULLTYPES) for (d in DRAWS) {
  message("[localscore] ", nt, " continuous draw ", d)
  summary_rows[[length(summary_rows) + 1]] <- run_continuous(nt, d)
  message("[localscore] ", nt, " mitoC2 draw ", d)
  summary_rows[[length(summary_rows) + 1]] <- run_mitoC2(nt, d)
}

summary_dt <- rbindlist(summary_rows)
setorder(summary_dt, mode, nulltype, draw)
fwrite(summary_dt, file.path(OUT_DIR, "localscore_window_counts.tsv"), sep = "\t")
saveRDS(summary_dt, file.path(OUT_DIR, "localscore_window_counts.rds"))

cat("\n=== window count summary ===\n")
print(summary_dt[, .(mean_n_windows = mean(n_windows), sd_n_windows = sd(n_windows),
                      min_n_windows = min(n_windows), max_n_windows = max(n_windows)),
                  by = .(mode, nulltype)])
message("[localscore] done: 40 runs, output in ", OUT_DIR)
