## =========================================================================
## module_population_partitioning -- follow-up 12: ancestry-run proxy for
## recent introgression. Does NOT modify any existing pp_*.R script or
## output, or R/10_geographic_prediction.R / R/11_individual_influence.R.
##
## Question: do individuals/populations carry unusually long runs of
## ancestry-homozygous genotype along the genetic map -- a signature
## consistent with recent backcrossing -- and do population-level run
## metrics predict the other follow-up results (individual influence,
## residual-profile magnitude, geographic residual axes, contribution to
## high-DI differentiation)?
##
## *** Authoritative local-ancestry tract calls: NONE FOUND ***
## Searched the repository (grep across all *.R, case-insensitive, for
## "local.ancestry", "ancestry.hmm"/"ancestry_hmm", "loter", "elai",
## "tract.call", "ancestry.tract", "hapmix", "rfmix", "lamp") before writing
## any code, per the brief. Every hit is an incidental mention of
## "ancestry-tract" as a THEORETICAL quantity in comments about demographic
## simulation (e.g. dev/R/moduleD_ohta_dmi.R's Haldane-mapping-function
## comment, dev/R/moduleE_analyze_sweep.R's "tract clock" comment) -- none is
## an implementation of, or a saved output from, an actual local-ancestry
## caller (no RFMix/ELAI/Loter/Ancestry_HMM usage anywhere). This module
## therefore implements ONLY a conservative, explicitly-labelled ANCESTRY-RUN
## PROXY on unphased genotypes, never called an "ancestry tract" in any
## output. See the closing note for what this can and cannot distinguish.
##
## *** DATASET NOTE (same as R/10, R/11) *** uses the current rho05 primary
## dataset (20,807 DI25 units), not the legacy 11,052-unit lineage.
##
## Inspected module_di25/R_legacy/di25_pruning_test.R for implementation
## ideas (not reused as authoritative, per the brief) -- notably its
## explicit finding that raw tract/run length ALONE cannot separate "F1 +
## scattered genotyping noise" from "early backcross" (an F1 with ~20%
## homozygous units by chance already produces short median runs by pure
## scatter); it addresses this with a spatial clustering (Wald-Wolfowitz
## runs) test, which is NOT implemented here (out of scope for this pass)
## but is flagged as a natural follow-up rather than silently ignored.
##
## Run detection: uses the SAME representative-SNP/unit-level oriented
## genotypes as every other analysis in this module (one marker per
## LD-reduced unit, NOT all 51,612 raw DI25 SNPs) -- using raw SNPs would
## just re-introduce the LD pseudoreplication the whole module's unit
## construction exists to remove. Units ordered by GENETIC (cM) position
## (module_population_partitioning/data/pp_recombination.rds, itself
## interpolated from data/Frufa_DTOL_PR.ref_genome.recmap), not physical
## position. A run is a maximal set of CONSECUTIVE (in cM order) units with
## IDENTICAL oriented genotype state (homozygous-aquilonia / homozygous-
## polyctena / heterozygous); heterozygous runs are their own category, not
## merged with either homozygous state. Missing genotypes BREAK a run in
## the primary analysis (no bridging) -- implemented by recoding NA to a
## per-position-unique sentinel before rle(), which forces a break at every
## missing call. A one-marker-gap-bridging rule is examined separately as a
## named sensitivity analysis, never the primary result.
##
## Run from the formica_hybrid repo root, after pp_recombination.R,
## R/10_geographic_prediction.R and R/11_individual_influence.R:
##   Rscript module_population_partitioning/R/12_ancestry_run_proxy.R
## Writes: module_population_partitioning/data/followup/12_ancestry_run_proxy.rds
##         module_population_partitioning/Figures/followup/12_*.png
## =========================================================================
suppressMessages(library(data.table))
t_start <- Sys.time()
OUTDIR <- "module_population_partitioning/data/followup"
FIGDIR <- "module_population_partitioning/Figures/followup"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------------
## 1. reconstruct individual-level oriented genotypes (identical
##    reconstruction to R/11_individual_influence.R -- G was not saved to
##    disk there either; re-derived and re-verified here independently)
## ---------------------------------------------------------------------
devtools::load_all("~/gitlab/LDscnR/", quiet = TRUE)
obj <- readRDS("module_population_partitioning/data/pp_units_Fmat.rds")
rc <- readRDS("module_population_partitioning/data/pp_recombination.rds")   # for cM_pos
u <- copy(rc$u); Fmat <- obj$Fmat; setDT(u); setorder(u, ChrNum, Pos)
stopifnot(all(u$group_id %in% colnames(Fmat)))
Fmat <- Fmat[, u$group_id, drop = FALSE]   # align column order to u (cM-ordered per chromosome)

inp <- readRDS("module_di25/data/di25_inputs.rds")
e2 <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd <- e2$sample_data_with_parents
GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)
GTs_all <- GTs_all[rownames(GTs_all) %in% sd$Sample_ID, ]
pops_all <- sd$Population[match(rownames(GTs_all), sd$Sample_ID)]
aqu_pops <- "aquilonia_parent"; pol_pops <- "polyctena_parent"
parent_rows <- grepl("_parent$", pops_all)

E <- GTs_all[, u$unit_marker, drop = FALSE]; colnames(E) <- u$group_id
prep <- ohta_fast_prepare(E, pops = pops_all)
P <- prep$pop_means / 2
sign_aqu <- sign(P[aqu_pops, ] - P[pol_pops, ])
flip <- which(sign_aqu < 0); undef <- which(is.na(sign_aqu) | sign_aqu == 0)
G <- E[!parent_rows, , drop = FALSE] / 2
if (length(flip)) G[, flip] <- 1 - G[, flip]
if (length(undef)) G[, undef] <- NA_real_
pops_hyb <- pops_all[!parent_rows]
stopifnot(identical(colnames(G), u$group_id))
chk <- colMeans(G[pops_hyb == "Karsikas", ], na.rm = TRUE)
stopifnot(max(abs(chk - Fmat["Karsikas", ]), na.rm = TRUE) < 1e-8)
cat(sprintf("[run] %d hybrid individuals x %d units, oriented genotypes verified against Fmat\n", nrow(G), ncol(G)))

## ---------------------------------------------------------------------
## 2. run detection: per individual x chromosome, units in cM order
## ---------------------------------------------------------------------
STATE <- ifelse(G > 0.9, "aqu", ifelse(G < 0.1, "pol", ifelse(!is.na(G), "het", NA_character_)))
dimnames(STATE) <- dimnames(G)
chrs <- unique(u$Chr)

detect_runs <- function(state_row, chr_vec, cm_vec, bridge_one_gap = FALSE) {
  out <- vector("list", length(unique(chr_vec)))
  for (k in seq_along(unique(chr_vec))) {
    ch <- unique(chr_vec)[k]
    idx <- which(chr_vec == ch)               # already cM-ordered (u was setorder'd upstream)
    s <- state_row[idx]; cm <- cm_vec[idx]
    if (bridge_one_gap) {
      ## sensitivity only: a single NA flanked by two markers of the SAME state
      ## is recoded to that state (bridges exactly one missing/discordant call)
      for (i in seq_along(s)) {
        if (is.na(s[i]) && i > 1 && i < length(s) && !is.na(s[i - 1]) && !is.na(s[i + 1]) && s[i - 1] == s[i + 1]) s[i] <- s[i - 1]
      }
    }
    key <- ifelse(is.na(s), paste0("NA_", seq_along(s)), s)   # forces a break at every remaining NA
    r <- rle(key)
    ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
    keep_run <- !startsWith(r$values, "NA_")
    if (any(keep_run)) {
      out[[k]] <- data.table(Chr = ch, state = r$values[keep_run], n_markers = r$lengths[keep_run],
                             cM_start = cm[starts[keep_run]], cM_end = cm[ends[keep_run]])
    }
  }
  rbindlist(out)
}

t0 <- Sys.time()
all_runs <- rbindlist(lapply(rownames(G), function(id) {
  r <- detect_runs(STATE[id, ], u$Chr, u$cM_pos)
  if (nrow(r)) r[, Sample_ID := id]
  r
}))
all_runs[, run_cM := cM_end - cM_start]
cat(sprintf("[run] primary (no gap-bridging) run detection: %d runs across %d individuals, %.1fs\n",
            nrow(all_runs), length(rownames(G)), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat("[run] run-state counts:\n"); print(all_runs[, .N, by = state])

## sensitivity: one-marker-gap bridging rule
all_runs_bridged <- rbindlist(lapply(rownames(G), function(id) {
  r <- detect_runs(STATE[id, ], u$Chr, u$cM_pos, bridge_one_gap = TRUE)
  if (nrow(r)) r[, Sample_ID := id]
  r
}))
all_runs_bridged[, run_cM := cM_end - cM_start]
cat(sprintf("[run] sensitivity (one-marker-gap bridging): %d runs (vs %d primary)\n", nrow(all_runs_bridged), nrow(all_runs)))

## sensitivity: stricter marker-informativeness threshold (current_map_DI > -15
## instead of the full DI25 panel's DI > -25 ascertainment)
strict_units <- u[!is.na(current_map_DI) & current_map_DI > -15, group_id]
cat(sprintf("[run] marker-informativeness sensitivity: %d/%d units retained at current_map_DI > -15\n", length(strict_units), nrow(u)))
u_strict <- u[group_id %in% strict_units]; setorder(u_strict, ChrNum, Pos)
G_strict <- G[, u_strict$group_id, drop = FALSE]
STATE_strict <- STATE[, u_strict$group_id, drop = FALSE]
all_runs_strict <- rbindlist(lapply(rownames(G_strict), function(id) {
  r <- detect_runs(STATE_strict[id, ], u_strict$Chr, u_strict$cM_pos)
  if (nrow(r)) r[, Sample_ID := id]
  r
}))
all_runs_strict[, run_cM := cM_end - cM_start]
cat(sprintf("[run] informativeness-sensitivity run detection: %d runs\n", nrow(all_runs_strict)))

## ---------------------------------------------------------------------
## 3. per-individual summaries (primary, un-bridged, full DI25 panel)
## ---------------------------------------------------------------------
callable_cM <- u[, .(span = max(cM_pos) - min(cM_pos)), by = Chr][, sum(span)]
cat(sprintf("\n[run] callable genetic map (fixed, panel-wide): %.1f cM\n", callable_cM))
LONG_THRESH <- c(0.5, 1, 2, 5)   # several prespecified thresholds, per the brief -- no single arbitrary cut

summarize_indiv <- function(runs_dt) {
  homo <- runs_dt[state %in% c("aqu", "pol")]
  out <- rbindlist(lapply(unique(runs_dt$Sample_ID), function(id) {
    ra <- homo[Sample_ID == id & state == "aqu", run_cM]
    rp <- homo[Sample_ID == id & state == "pol", run_cM]
    frac_long <- sapply(LONG_THRESH, function(th) sum(homo[Sample_ID == id & run_cM >= th, run_cM]) / callable_cM)
    names(frac_long) <- paste0("frac_long_", LONG_THRESH, "cM")
    as.data.table(c(list(Sample_ID = id, n_runs_aqu = length(ra), n_runs_pol = length(rp),
                        median_aqu_cM = if (length(ra)) median(ra) else NA_real_,
                        p90_aqu_cM = if (length(ra)) quantile(ra, 0.9) else NA_real_,
                        max_aqu_cM = if (length(ra)) max(ra) else NA_real_,
                        median_pol_cM = if (length(rp)) median(rp) else NA_real_,
                        p90_pol_cM = if (length(rp)) quantile(rp, 0.9) else NA_real_,
                        max_pol_cM = if (length(rp)) max(rp) else NA_real_),
                   as.list(frac_long)))
  }))
  out
}
indiv_runs <- summarize_indiv(all_runs)
indiv_runs[, Population := pops_hyb[match(Sample_ID, rownames(G))]]
setorder(indiv_runs, -frac_long_2cM)
cat("\n[run] per-individual run summary, top 10 by fraction of genome in runs >=2cM:\n")
print(indiv_runs[1:10, .(Sample_ID, Population, n_runs_aqu, n_runs_pol, median_aqu_cM, median_pol_cM, frac_long_2cM)])

## ---------------------------------------------------------------------
## 4. population-level run metrics (populations are the biological
##    replicate here, per the brief -- n=20)
## ---------------------------------------------------------------------
pop_runs <- indiv_runs[, .(n_indiv = .N,
                           mean_median_aqu_cM = mean(median_aqu_cM, na.rm = TRUE),
                           mean_median_pol_cM = mean(median_pol_cM, na.rm = TRUE),
                           mean_frac_long_1cM = mean(frac_long_1cM, na.rm = TRUE),
                           mean_frac_long_2cM = mean(frac_long_2cM, na.rm = TRUE),
                           mean_frac_long_5cM = mean(frac_long_5cM, na.rm = TRUE)),
                       by = Population]
setorder(pop_runs, -mean_frac_long_2cM)
cat("\n[run] population-level run metrics (n=20 populations, ordered by mean_frac_long_2cM):\n")
print(pop_runs)

## ---------------------------------------------------------------------
## 5. relate population-level run metrics to the 4 required targets
## ---------------------------------------------------------------------
infl <- readRDS("module_population_partitioning/data/followup/11_individual_influence.rds")
pop_influence <- infl$indiv_tab[, .(mean_influence = mean(influence_rms)), by = Population]

rez <- readRDS("module_population_partitioning/data/pp_residual_ancestry.rds")
pop_resid_mag <- data.table(Population = rownames(rez$Resid), resid_rms = sqrt(rowMeans(rez$Resid^2, na.rm = TRUE)))

geo <- readRDS("module_population_partitioning/data/followup/10_geographic_prediction.rds")
pop_pca <- geo$pca_scores[, .(pop, PC1, PC2)]; setnames(pop_pca, "pop", "Population")

fbar <- colMeans(Fmat, na.rm = TRUE)
pop_contrib <- data.table(Population = rownames(Fmat), contrib_diff = rowMeans(sweep(Fmat, 2, fbar, "-")^2, na.rm = TRUE))

targets <- Reduce(function(a, b) merge(a, b, by = "Population"), list(pop_runs, pop_influence, pop_resid_mag, pop_pca, pop_contrib))
stopifnot(nrow(targets) == 20)

metric_col <- "mean_frac_long_2cM"   # primary run-metric used against all 4 targets
test_targets <- c(mean_influence = "individual influence (Analysis 2)", resid_rms = "residual-profile magnitude",
                  PC1 = "geographic residual PC1 (Analysis 1)", PC2 = "geographic residual PC2 (Analysis 1)",
                  contrib_diff = "contribution to high-DI differentiation")
assoc <- rbindlist(lapply(names(test_targets), function(tc) {
  rho <- cor(targets[[metric_col]], targets[[tc]], method = "spearman")
  loo <- sapply(seq_len(20), function(i) cor(targets[[metric_col]][-i], targets[[tc]][-i], method = "spearman"))
  data.table(target = test_targets[tc], rho = rho, loo_min = min(loo), loo_max = max(loo))
}))
cat(sprintf("\n[run] population-level run metric (%s) vs the 4 required targets (n=20, Spearman + leave-one-pop-out range):\n", metric_col))
print(assoc)

## ---------------------------------------------------------------------
## 6. figures
## ---------------------------------------------------------------------
suppressMessages(library(ggplot2))
theme_ms <- theme_bw(base_size = 12) + theme(strip.background = element_blank(), panel.grid.minor = element_blank())

homo_runs <- all_runs[state %in% c("aqu", "pol")]
fig_dist <- ggplot(homo_runs, aes(run_cM, fill = state)) +
  geom_histogram(position = "identity", alpha = 0.55, bins = 60) +
  scale_fill_manual(values = c(aqu = "#21918C", pol = "#D3C93B"), labels = c(aqu = "F. aquilonia", pol = "F. polyctena")) +
  labs(x = "run length (cM)", y = "count", fill = NULL,
      title = "Continuous distribution of ancestry-homozygous run lengths (all individuals pooled)") +
  theme_ms
ggsave(file.path(FIGDIR, "12_run_length_distribution.png"), fig_dist, width = 8, height = 5.5, dpi = 200)

indiv_runs[, Population := factor(Population, levels = pop_runs[order(mean_frac_long_2cM), Population])]
fig_by_pop <- ggplot(indiv_runs, aes(Population, frac_long_2cM)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") + geom_jitter(width = 0.15, size = 1.6, colour = "grey30") +
  labs(x = NULL, y = "fraction of genetic map in runs >=2cM", title = "Long-run fraction by population") +
  theme_ms + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(FIGDIR, "12_long_run_fraction_by_population.png"), fig_by_pop, width = 10, height = 5.5, dpi = 200)

assoc[, target := factor(target, levels = target)]
fig_assoc <- ggplot(assoc, aes(target, rho)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 3, colour = "#7570b3") + geom_errorbar(aes(ymin = loo_min, ymax = loo_max), width = 0.15, colour = "#7570b3") +
  labs(x = NULL, y = "Spearman rho (population-level, n=20)\npoint = full sample, bars = leave-one-population-out range",
      title = "Population-level long-run fraction vs the 4 required targets") +
  theme_ms + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(FIGDIR, "12_run_metric_associations.png"), fig_assoc, width = 8, height = 6, dpi = 200)
cat("\n[run] figures saved: 12_run_length_distribution.png, 12_long_run_fraction_by_population.png, 12_run_metric_associations.png\n")

## ---------------------------------------------------------------------
## 7. save
## ---------------------------------------------------------------------
cat("\n*** CAUTION (required statement, not a finding) ***\n")
cat("Unphased ancestry runs on ~20,800 LD-reduced units cannot distinguish all recent\n")
cat("backcross histories from other explanations (e.g. an F1-like individual can show\n")
cat("short homozygous runs by pure chance without any backcrossing; conversely phase\n")
cat("errors/genotyping noise can fragment a genuine long run). This is a conservative\n")
cat("PROXY, not a validated ancestry-tract call. If an authoritative tract caller\n")
cat("becomes available, retain this analysis as a sensitivity check only.\n")

result <- list(all_runs = all_runs, all_runs_bridged_sensitivity = all_runs_bridged,
               all_runs_informativeness_sensitivity = all_runs_strict,
               indiv_runs = indiv_runs, pop_runs = pop_runs, targets = targets, assoc = assoc,
               callable_cM = callable_cM, long_thresholds = LONG_THRESH,
               caveat = "conservative unphased-genotype ancestry-run PROXY, not a validated local-ancestry tract call; no authoritative caller found in this repository",
               dataset_note = "uses the current rho05 (20,807-unit) dataset",
               session_info = sessionInfo(), run_time = Sys.time(), elapsed_secs = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
saveRDS(result, file.path(OUTDIR, "12_ancestry_run_proxy.rds"))
cat(sprintf("\n[run] saved -> %s (elapsed %.1fs)\n", file.path(OUTDIR, "12_ancestry_run_proxy.rds"), result$elapsed_secs))
