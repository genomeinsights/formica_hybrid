## =========================================================================
## module_population_partitioning -- candidate-locus follow-up 20: item 5,
## matched-Stage-1-unit-set null comparison.
##
## For each target's RAW candidate set (the most complete, statistically
## meaningful set per target -- floor/top-N sets are too small for a
## matched-null comparison to be informative and are skipped here), draw
## B=1000 matched SETS of real Stage-1-direct units: one match per
## candidate, drawn from an admissible pool (all 18,361 units minus the
## union of ALL FOUR targets' raw candidates, so a locus can never stand in
## as "background" for a scan that itself flagged it), matched on DI
## (+-2), folded parental MAF (+-0.05), recombination rate (+-50%), cluster
## size (+-3 markers), same-chromosome preferred with genome-wide fallback
## if the same-chromosome pool is empty. Matching is WITHOUT REPLACEMENT
## within one draw (no Stage-1 unit stands in for two different candidates
## in the same draw) but WITH REPLACEMENT ACROSS draws (each of the 1000
## draws is independent). Every matched-null unit's profile is a REAL,
## OBSERVED profile of a real physically-linked genomic unit -- nothing is
## synthetically permuted, satisfying "physical linkage is not broken by
## independently permuting population values".
##
## Two statistics are compared, observed-candidate-set vs B=1000-draw null,
## on BOTH raw and structure-adjusted (PC1-residual) profiles:
##   (i)  the same near/far/cross-chr mean |r| contrast from script 18/19
##   (ii) "structure_dominance" = leading-eigenvalue fraction of the locus x
##        locus correlation matrix -- directly operationalizes the main
##        inferential question (near 1 = one dominant genome-wide division;
##        well below 1, comparable PC1/PC2/PC3 eigenvalues = several
##        reproducible, distinguishable groups).
##
## Run from the formica_hybrid repo root, after 19_structure_adjustment.R:
##   Rscript module_population_partitioning/R/20_matched_null.R
## Reads : module_population_partitioning/data/followup/{15,18,19}_*.rds
## Writes: module_population_partitioning/data/followup/20_matched_null.rds
##         module_population_partitioning/Figures/followup/20_matched_null_<target>.png
## =========================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
t_start <- Sys.time()
OUTDIR <- "module_population_partitioning/data/followup"
FIGDIR <- "module_population_partitioning/Figures/followup"
r15 <- readRDS(file.path(OUTDIR, "15_candidate_data.rds"))
r18 <- readRDS(file.path(OUTDIR, "18_similarity_raw.rds"))
r19 <- readRDS(file.path(OUTDIR, "19_structure_adjustment.rds"))
u <- r15$u; Fmat_all <- r15$Fmat_all; Resid1 <- r19$Resid1
stopifnot("candidate unit table must have exactly 18,361 rows -- rerun 15_candidate_import.R" = nrow(u) == 18361L)

DI_TOL <- 2; MAF_TOL <- 0.05; RECOMB_TOL_PCT <- 0.5; NLOCI_TOL <- 3
B_MATCH <- 1000; NEAR_BP <- 1e5
set.seed(20260914)
theme_ms <- theme_bw(base_size = 11) + theme(strip.background = element_blank(), panel.grid.minor = element_blank())

all_raw_candidates <- unique(unlist(r15$cand_ids))
pool <- u[!group_id %in% all_raw_candidates]
cat(sprintf("[matched-null] admissible pool: %d/%d Stage-1-direct units (excludes %d raw candidates across all 4 targets)\n",
            nrow(pool), nrow(u), length(all_raw_candidates)))

## ---------------------------------------------------------------------
## structure_dominance: leading-eigenvalue fraction of the locus x locus
## correlation matrix
## ---------------------------------------------------------------------
structure_dominance <- function(Fsub) {
  keep <- colSums(!is.na(Fsub)) > 0
  Fsub <- Fsub[, keep, drop = FALSE]
  if (ncol(Fsub) < 3) return(NA_real_)
  Rm <- suppressWarnings(cor(Fsub, use = "pairwise.complete.obs"))
  if (anyNA(Rm)) return(NA_real_)
  ev <- eigen(Rm, symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(ev, 0)
  ev[1] / sum(ev)
}
dcat_stat <- function(Fsub, sub_u) {
  n <- ncol(Fsub)
  Rm <- suppressWarnings(cor(Fsub, use = "pairwise.complete.obs"))
  ii <- rep(seq_len(n), times = n); jj <- rep(seq_len(n), each = n); keep <- ii < jj; ii <- ii[keep]; jj <- jj[keep]
  same_chr <- sub_u$Chr[ii] == sub_u$Chr[jj]
  dist_bp <- ifelse(same_chr, abs(sub_u$Pos[ii] - sub_u$Pos[jj]), NA_real_)
  dcat <- ifelse(!same_chr, "cross-chr", ifelse(dist_bp <= NEAR_BP, "near", "far"))
  absr <- abs(Rm[cbind(ii, jj)])
  c(near = mean(absr[dcat == "near"], na.rm = TRUE), far = mean(absr[dcat == "far"], na.rm = TRUE),
    cross = mean(absr[dcat == "cross-chr"], na.rm = TRUE))
}

## ---------------------------------------------------------------------
## precompute, per candidate, its admissible pool (same-chr + genome-wide
## fallback) -- fixed across draws, only the within-draw "already used" set
## changes
## ---------------------------------------------------------------------
precompute_pools <- function(cand_ids) {
  cand <- u[match(cand_ids, group_id)]
  lapply(seq_len(nrow(cand)), function(i) {
    c1 <- cand[i]
    elig <- pool[abs(DI - c1$DI) <= DI_TOL & abs(pmaf - c1$pmaf) <= MAF_TOL &
                abs(recomb - c1$recomb) <= RECOMB_TOL_PCT * c1$recomb & abs(n_loci - c1$n_loci) <= NLOCI_TOL]
    list(same_chr = elig[Chr == c1$Chr, group_id], genome_wide = elig$group_id)
  })
}
match_one_set <- function(elig_list) {
  used <- character(0); matched <- character(length(elig_list))
  for (i in seq_along(elig_list)) {
    cand_pool <- setdiff(elig_list[[i]]$same_chr, used)
    if (!length(cand_pool)) cand_pool <- setdiff(elig_list[[i]]$genome_wide, used)
    if (!length(cand_pool)) { matched[i] <- NA_character_; next }
    pick <- cand_pool[sample.int(length(cand_pool), 1)]
    matched[i] <- pick; used <- c(used, pick)
  }
  matched
}

## ---------------------------------------------------------------------
## per-target matched-null run (raw candidate set only)
## ---------------------------------------------------------------------
run_matched_null <- function(tgt) {
  cand_ids <- r15$cand_ids[[tgt]]
  elig_list <- precompute_pools(cand_ids)
  n_no_match <- sum(vapply(elig_list, function(e) length(e$genome_wide) == 0, logical(1)))
  cat(sprintf("[matched-null] %-10s: %d candidates, %d with zero admissible pool members (excluded from all draws)\n",
              tgt, length(cand_ids), n_no_match))

  null_stat <- matrix(NA_real_, B_MATCH, 4, dimnames = list(NULL, c("near", "far", "cross", "structure_dominance")))
  null_stat_adj <- matrix(NA_real_, B_MATCH, 4, dimnames = list(NULL, c("near", "far", "cross", "structure_dominance")))
  n_matched_vec <- integer(B_MATCH)
  for (b in seq_len(B_MATCH)) {
    matched <- match_one_set(elig_list)
    matched <- matched[!is.na(matched)]
    n_matched_vec[b] <- length(matched)
    if (length(matched) < 3) next
    sub_u <- u[match(matched, group_id)]
    Fsub <- Fmat_all[, matched, drop = FALSE]
    ds <- dcat_stat(Fsub, sub_u)
    null_stat[b, ] <- c(ds, structure_dominance(Fsub))
    Fsub_adj <- Resid1[, matched, drop = FALSE]
    ds_adj <- dcat_stat(Fsub_adj, sub_u)
    null_stat_adj[b, ] <- c(ds_adj, structure_dominance(Fsub_adj))
  }
  cat(sprintf("    -> %d/%d draws had >=3 matched loci; median matched = %d/%d candidates\n",
              sum(n_matched_vec >= 3), B_MATCH, median(n_matched_vec), length(cand_ids)))

  ## observed candidate-set values, raw (from script 18) and adjusted (script 19)
  dcat_to_key <- c("same-chr near (<=100kb)" = "near", "same-chr far" = "far", "cross-chr" = "cross")
  obs_raw <- r18$sim_results[[paste0(tgt, "_raw")]]
  obs_adj <- r19$adj_results_pc1[[paste0(tgt, "_raw")]]
  obs_raw_dcat <- if (!is.null(obs_raw)) setNames(obs_raw$obs$mean_absr, dcat_to_key[obs_raw$obs$dcat]) else c(near = NA, far = NA, cross = NA)
  obs_adj_dcat <- if (!is.null(obs_adj)) setNames(obs_adj$obs$mean_absr, dcat_to_key[obs_adj$obs$dcat]) else c(near = NA, far = NA, cross = NA)
  obs_raw_dcat <- obs_raw_dcat[c("near", "far", "cross")]; obs_adj_dcat <- obs_adj_dcat[c("near", "far", "cross")]
  obs_sd_raw <- structure_dominance(Fmat_all[, cand_ids, drop = FALSE])
  obs_sd_adj <- structure_dominance(Resid1[, cand_ids, drop = FALSE])

  p_value <- function(obs_v, null_v) mean(null_v >= obs_v, na.rm = TRUE)
  summary_tab <- data.table(
    target = tgt,
    statistic = c("cross-chr mean|r|", "same-chr far mean|r|", "same-chr near mean|r|", "structure_dominance"),
    obs_raw = c(obs_raw_dcat["cross"], obs_raw_dcat["far"], obs_raw_dcat["near"], obs_sd_raw),
    null_mean_raw = c(mean(null_stat[, "cross"], na.rm = TRUE), mean(null_stat[, "far"], na.rm = TRUE),
                      mean(null_stat[, "near"], na.rm = TRUE), mean(null_stat[, "structure_dominance"], na.rm = TRUE)),
    null_lo_raw = c(quantile(null_stat[, "cross"], 0.025, na.rm = TRUE), quantile(null_stat[, "far"], 0.025, na.rm = TRUE),
                    quantile(null_stat[, "near"], 0.025, na.rm = TRUE), quantile(null_stat[, "structure_dominance"], 0.025, na.rm = TRUE)),
    null_hi_raw = c(quantile(null_stat[, "cross"], 0.975, na.rm = TRUE), quantile(null_stat[, "far"], 0.975, na.rm = TRUE),
                    quantile(null_stat[, "near"], 0.975, na.rm = TRUE), quantile(null_stat[, "structure_dominance"], 0.975, na.rm = TRUE)),
    p_raw = c(p_value(obs_raw_dcat["cross"], null_stat[, "cross"]), p_value(obs_raw_dcat["far"], null_stat[, "far"]),
             p_value(obs_raw_dcat["near"], null_stat[, "near"]), p_value(obs_sd_raw, null_stat[, "structure_dominance"])),
    obs_adj = c(obs_adj_dcat["cross"], obs_adj_dcat["far"], obs_adj_dcat["near"], obs_sd_adj),
    null_mean_adj = c(mean(null_stat_adj[, "cross"], na.rm = TRUE), mean(null_stat_adj[, "far"], na.rm = TRUE),
                      mean(null_stat_adj[, "near"], na.rm = TRUE), mean(null_stat_adj[, "structure_dominance"], na.rm = TRUE)),
    p_adj = c(p_value(obs_adj_dcat["cross"], null_stat_adj[, "cross"]), p_value(obs_adj_dcat["far"], null_stat_adj[, "far"]),
             p_value(obs_adj_dcat["near"], null_stat_adj[, "near"]), p_value(obs_sd_adj, null_stat_adj[, "structure_dominance"])))
  print(summary_tab)

  ## figure: observed vs null distribution for cross-chr |r| and structure_dominance (raw)
  df1 <- data.table(stat = null_stat[, "cross"]); df1 <- df1[!is.na(stat)]
  p1 <- ggplot(df1, aes(stat)) + geom_histogram(bins = 40, fill = "grey75") +
    geom_vline(xintercept = obs_raw_dcat["cross"], colour = "firebrick", linewidth = 1) +
    labs(title = "cross-chr mean |r|: candidate (red) vs matched-null", x = "mean |r|", y = NULL) + theme_ms
  df2 <- data.table(stat = null_stat[, "structure_dominance"]); df2 <- df2[!is.na(stat)]
  p2 <- ggplot(df2, aes(stat)) + geom_histogram(bins = 40, fill = "grey75") +
    geom_vline(xintercept = obs_sd_raw, colour = "firebrick", linewidth = 1) +
    labs(title = "structure_dominance: candidate (red) vs matched-null", x = "leading-eigenvalue fraction", y = NULL) + theme_ms
  fig <- p1 + p2
  ggsave(file.path(FIGDIR, sprintf("20_matched_null_%s.png", tgt)), fig, width = 10, height = 4.5, dpi = 200)

  list(target = tgt, n_candidates = length(cand_ids), n_no_match = n_no_match, n_matched_vec = n_matched_vec,
      null_stat = null_stat, null_stat_adj = null_stat_adj, summary_tab = summary_tab)
}

matched_null_results <- lapply(names(r15$TARGETS), run_matched_null)
names(matched_null_results) <- names(r15$TARGETS)

result <- list(matched_null_results = matched_null_results, DI_TOL = DI_TOL, MAF_TOL = MAF_TOL,
              RECOMB_TOL_PCT = RECOMB_TOL_PCT, NLOCI_TOL = NLOCI_TOL, B_MATCH = B_MATCH,
              session_info = sessionInfo(), run_time = Sys.time(),
              elapsed_secs = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
saveRDS(result, file.path(OUTDIR, "20_matched_null.rds"))
cat(sprintf("\n[matched-null] saved -> %s (elapsed %.1fs)\n", file.path(OUTDIR, "20_matched_null.rds"), result$elapsed_secs))
