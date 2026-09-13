## =========================================================================
## module_population_partitioning -- follow-up 12: UNPHASED GENOTYPE-STATE
## RUN PROXY (NOT an ancestry-tract / local-ancestry analysis -- see caveat
## below and throughout). Does NOT modify any existing pp_*.R script or
## output, or R/10_geographic_prediction.R / R/11_individual_influence.R.
##
## REVISED (audit response item 5) after a first pass that had several
## fixable problems, now addressed:
##   - runs no longer bridge arbitrarily large genetic-map gaps between
##     same-state markers (a maximum-gap threshold now breaks a run even
##     when state is unchanged -- see part 2 below);
##   - run length no longer silently includes unobserved genetic interval:
##     it is only ever the span between two markers actually accepted into
##     the same segment;
##   - callable genetic length is computed SEPARATELY PER INDIVIDUAL
##     (depends on that individual's own missing calls), not one fixed
##     panel-wide number;
##   - single-marker segments (0 cM by construction) are explicitly flagged
##     and distinguished from multi-marker (>=2 markers) runs throughout;
##   - marker count and genetic length are reported jointly, never length
##     alone;
##   - the gap-length distribution is inspected BEFORE choosing the primary
##     gap threshold (part 1), with >=2 additional thresholds as named
##     sensitivity analyses;
##   - population-level associations are now explicitly flagged as
##     non-independent (item 10) -- run metrics, residual magnitude, and
##     differentiation contribution are all computed from the same
##     underlying genotypes/allele frequencies, so a correlation between
##     them is not 3 separate lines of evidence;
##   - nothing here is interpreted as evidence FOR or AGAINST recent
##     backcrossing -- this proxy cannot resolve introgression timing (item
##     9); it is a purely descriptive genotype-state summary.
##
## Question this DESCRIBES (not resolves): how much of each individual's/
## population's genome is covered by short vs. longer runs of consecutive,
## same-parental-state, genetic-map-ordered LD-reduced-unit genotypes?
##
## *** Authoritative local-ancestry tract calls: NONE FOUND ***
## Searched the repository (grep across all *.R, case-insensitive, for
## "local.ancestry", "ancestry.hmm"/"ancestry_hmm", "loter", "elai",
## "tract.call", "ancestry.tract", "hapmix", "rfmix", "lamp") before writing
## any code. Every hit is an incidental mention of "ancestry-tract" as a
## THEORETICAL quantity in demographic-simulation comments (e.g.
## dev/R/moduleD_ohta_dmi.R's Haldane-mapping-function comment) -- no
## implementation of, or saved output from, an actual local-ancestry caller
## exists anywhere in this repository.
##
## *** WHAT THIS IS AND IS NOT (required statement) ***
## This is an UNPHASED GENOTYPE-STATE RUN PROXY on LD-reduced representative
## SNPs. It is NOT a probabilistic local-ancestry call, NOT phased, and the
## representative-SNP marker set (one marker per LD cluster, chosen for the
## population-partitioning analysis, not for ancestry-tract inference) is
## not an ideal marker set for this purpose either. A "run" here means
## "consecutive same-state observed genotypes within a chosen genetic-map
## gap tolerance" -- nothing more. It cannot: distinguish recent
## backcrossing from old admixture followed by drift, correct for phase
## ambiguity, or be interpreted as excluding (or supporting) any specific
## introgression timing. Predominantly short runs are NOT evidence against
## recent backcrossing (an F1 individual, or an old, fully-recombined
## admixed individual, can both show short runs for different reasons that
## this proxy cannot tell apart).
##
## *** DATASET NOTE (same as R/10, R/11) *** uses the current rho05 primary
## dataset (20,807 DI25 units), not the legacy 11,052-unit lineage.
##
## Inspected module_di25/R_legacy/di25_pruning_test.R for implementation
## ideas (not reused as authoritative): its explicit finding that raw
## run/tract length ALONE cannot separate "F1 + scattered genotyping noise"
## from "early backcross" motivates the caveat above; its spatial
## (Wald-Wolfowitz) clustering test is NOT implemented here (out of scope).
##
## Run detection: uses the SAME representative-SNP/unit-level oriented
## genotypes as every other analysis in this module (one marker per
## LD-reduced unit, NOT all 51,612 raw DI25 SNPs). Units ordered by GENETIC
## (cM) position. A run is a maximal set of consecutive OBSERVED
## (non-missing) units with IDENTICAL genotype state (homozygous-aquilonia /
## homozygous-polyctena / heterozygous; heterozygous kept separate from
## either homozygous state), where the genetic distance between consecutive
## observed units does not exceed a maximum-gap threshold -- exceeding it
## breaks the run even when the state either side is identical.
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
## 0. reconstruct individual-level oriented genotypes (identical
##    reconstruction to R/11_individual_influence.R -- G was not saved to
##    disk there either; re-derived and re-verified here independently)
## ---------------------------------------------------------------------
devtools::load_all("~/gitlab/LDscnR/", quiet = TRUE)
obj <- readRDS("module_population_partitioning/data/pp_units_Fmat.rds")
rc <- readRDS("module_population_partitioning/data/pp_recombination.rds")   # for cM_pos
u <- copy(rc$u); Fmat <- obj$Fmat; setDT(u); setorder(u, ChrNum, Pos)
stopifnot("unit table must have exactly the 20,807 DI25 rho05 units (min_r2_rho=0.5) -- check for a legacy/stale input" = nrow(u) == 20807L,
         "all unit_marker/group_id values must be present in Fmat's columns" = all(u$group_id %in% colnames(Fmat)))
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
## 1. inspect the adjacent-OBSERVED-marker cM-gap distribution BEFORE
##    choosing the primary maximum-gap threshold (panel-wide, ignoring any
##    individual's own missingness -- missingness is negligible, ~0.1%, so
##    this is a good guide to the threshold; per-individual gaps used for
##    the actual run detection below account for each individual's own
##    missing calls exactly)
## ---------------------------------------------------------------------
panel_gaps <- u[, .(gap = diff(cM_pos)), by = Chr]$gap
cat("\n[run] adjacent-marker cM-gap distribution (panel-wide, for threshold selection):\n")
print(summary(panel_gaps))
gap_pctile <- quantile(panel_gaps, c(0.5, 0.75, 0.9, 0.95, 0.99, 0.999))
cat("percentiles:\n"); print(round(gap_pctile, 3))
for (th in c(0.1, 0.5, 1, 2)) cat(sprintf("  fraction of adjacent gaps > %.1fcM: %.2f%%\n", th, 100 * mean(panel_gaps > th)))
cat("\n[run] PRIMARY max-gap threshold chosen: 1.0 cM (sits at the ~99th percentile of the panel-wide\n")
cat("    adjacent-marker gap distribution above -- only exceptionally sparse regions are excluded from\n")
cat("    bridging; sensitivity thresholds 0.5cM (stricter) and 2.0cM (more permissive) below).\n")

## ---------------------------------------------------------------------
## 2. run detection: per individual x chromosome, OBSERVED units only, in
##    cM order. Breaks a run when state changes OR the gap to the next
##    observed unit exceeds max_gap (even if state is unchanged).
## ---------------------------------------------------------------------
STATE <- ifelse(G > 0.9, "aqu", ifelse(G < 0.1, "pol", ifelse(!is.na(G), "het", NA_character_)))
dimnames(STATE) <- dimnames(G)
chrs <- unique(u$Chr)

detect_runs <- function(state_row, chr_vec, cm_vec, max_gap) {
  out <- vector("list", length(unique(chr_vec)))
  for (k in seq_along(unique(chr_vec))) {
    ch <- unique(chr_vec)[k]
    idx <- which(chr_vec == ch)                    # already cM-ordered
    s <- state_row[idx]; cm <- cm_vec[idx]
    obs <- which(!is.na(s))                         # OBSERVED positions only -- missingness never bridged
    if (length(obs) == 0) next
    s_obs <- s[obs]; cm_obs <- cm[obs]
    gap <- c(Inf, diff(cm_obs))                      # gap[i] = genetic distance from obs i-1 to obs i
    state_change <- c(TRUE, s_obs[-1] != s_obs[-length(s_obs)])
    break_here <- state_change | (gap > max_gap)     # AUDIT FIX: gap alone also breaks a run now
    run_id <- cumsum(break_here)
    dt <- data.table(run_id = run_id, state = s_obs, cm = cm_obs)
    out[[k]] <- dt[, .(Chr = ch, state = state[1], n_markers = .N, cM_start = min(cm), cM_end = max(cm)), by = run_id][, run_id := NULL]
  }
  rbindlist(out)
}

run_all_individuals <- function(max_gap, label) {
  t0 <- Sys.time()
  runs <- rbindlist(lapply(rownames(G), function(id) {
    r <- detect_runs(STATE[id, ], u$Chr, u$cM_pos, max_gap)
    if (nrow(r)) r[, Sample_ID := id]
    r
  }))
  runs[, run_cM := cM_end - cM_start]
  runs[, is_single_marker := n_markers == 1L]
  cat(sprintf("[run] %s (max_gap=%.1fcM): %d segments across %d individuals, %.1fs\n",
              label, max_gap, nrow(runs), length(rownames(G)), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  runs
}

all_runs <- run_all_individuals(1.0, "PRIMARY")
cat("[run] segment-state counts (primary):\n"); print(all_runs[, .N, by = state])
cat(sprintf("[run] single-marker segments: %d/%d (%.1f%%) -- these have run_cM=0 by construction, reported\n",
            sum(all_runs$is_single_marker), nrow(all_runs), 100 * mean(all_runs$is_single_marker)))
cat("    separately from multi-marker (>=2 observed markers) runs throughout, not conflated with them.\n")

all_runs_gap05 <- run_all_individuals(0.5, "sensitivity (stricter gap)")
all_runs_gap2 <- run_all_individuals(2.0, "sensitivity (more permissive gap)")

## sensitivity: stricter marker-informativeness threshold (current_map_DI > -15
## instead of the full DI25 panel's DI > -25 ascertainment)
strict_units <- u[!is.na(current_map_DI) & current_map_DI > -15, group_id]
cat(sprintf("\n[run] marker-informativeness sensitivity: %d/%d units retained at current_map_DI > -15\n", length(strict_units), nrow(u)))
u_strict <- u[group_id %in% strict_units]; setorder(u_strict, ChrNum, Pos)
G_strict <- G[, u_strict$group_id, drop = FALSE]; STATE_strict <- STATE[, u_strict$group_id, drop = FALSE]
all_runs_strict <- rbindlist(lapply(rownames(G_strict), function(id) {
  r <- detect_runs(STATE_strict[id, ], u_strict$Chr, u_strict$cM_pos, 1.0)
  if (nrow(r)) r[, Sample_ID := id]
  r
}))
all_runs_strict[, run_cM := cM_end - cM_start]
cat(sprintf("[run] informativeness-sensitivity run detection (max_gap=1.0cM): %d segments\n", nrow(all_runs_strict)))

## ---------------------------------------------------------------------
## 3. per-individual summaries (primary: max_gap=1.0cM, full DI25 panel).
##    callable_cM is now PER INDIVIDUAL (sum of accepted, i.e. <=max_gap,
##    genetic intervals between that individual's own observed calls) --
##    not one fixed panel-wide number.
## ---------------------------------------------------------------------
callable_per_indiv <- function(max_gap) {
  sapply(rownames(G), function(id) {
    s <- STATE[id, ]
    total <- 0
    for (ch in chrs) {
      idx <- which(u$Chr == ch); cm <- u$cM_pos[idx]; si <- s[idx]
      obs <- which(!is.na(si)); if (length(obs) < 2) next
      g <- diff(cm[obs])
      total <- total + sum(g[g <= max_gap])
    }
    total
  })
}
callable_cM_indiv <- callable_per_indiv(1.0)
cat(sprintf("\n[run] per-individual callable genetic length (max_gap=1.0cM): median %.1f cM, range [%.1f, %.1f] (panel-wide upper bound for reference: %.1f cM)\n",
            median(callable_cM_indiv), min(callable_cM_indiv), max(callable_cM_indiv),
            u[, .(span = max(cM_pos) - min(cM_pos)), by = Chr][, sum(span)]))

LONG_THRESH <- c(0.5, 1, 2, 5)   # several prespecified thresholds -- no single arbitrary cut

summarize_indiv <- function(runs_dt, callable_vec) {
  homo <- runs_dt[state %in% c("aqu", "pol")]
  multi <- homo[is_single_marker == FALSE]   # runs supported by >=2 markers -- the ones with a genuine cM length
  out <- rbindlist(lapply(names(callable_vec), function(id) {
    ra <- multi[Sample_ID == id & state == "aqu", run_cM]
    rp <- multi[Sample_ID == id & state == "pol", run_cM]
    n_single_aqu <- homo[Sample_ID == id & state == "aqu" & is_single_marker, .N]
    n_single_pol <- homo[Sample_ID == id & state == "pol" & is_single_marker, .N]
    cal <- callable_vec[id]
    frac_long <- sapply(LONG_THRESH, function(th) sum(multi[Sample_ID == id & run_cM >= th, run_cM]) / cal)
    names(frac_long) <- paste0("frac_long_", LONG_THRESH, "cM")
    as.data.table(c(list(Sample_ID = id, callable_cM = cal,
                        n_multimarker_runs_aqu = length(ra), n_single_marker_aqu = n_single_aqu,
                        n_multimarker_runs_pol = length(rp), n_single_marker_pol = n_single_pol,
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
indiv_runs <- summarize_indiv(all_runs, callable_cM_indiv)
indiv_runs[, Population := pops_hyb[match(Sample_ID, rownames(G))]]
setorder(indiv_runs, -frac_long_2cM)
cat("\n[run] per-individual summary (marker count AND genetic length reported jointly), top 10 by fraction in runs >=2cM:\n")
print(indiv_runs[1:10, .(Sample_ID, Population, callable_cM, n_multimarker_runs_aqu, n_multimarker_runs_pol,
                         median_aqu_cM, median_pol_cM, frac_long_2cM)])

## ---------------------------------------------------------------------
## 4. population-level run metrics (populations are the biological
##    replicate here -- n=20)
## ---------------------------------------------------------------------
pop_runs <- indiv_runs[, .(n_indiv = .N, mean_callable_cM = mean(callable_cM),
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
## 5. relate population-level run metrics to the 4 required targets.
##    AUDIT FIX (item 10): these are explicitly NOT independent evidence --
##    the run metric, residual-profile magnitude, and differentiation
##    contribution are all computed from the SAME underlying oriented
##    genotypes/Fmat, so a correlation between them reflects shared
##    measurement, not 3 separately-corroborating observations. Stated here
##    and repeated in every downstream consumer of this object.
## ---------------------------------------------------------------------
infl <- readRDS("module_population_partitioning/data/followup/11_individual_influence.rds")
pop_influence <- infl$indiv_tab[, .(mean_influence = mean(influence_rms)), by = Population]

rez <- readRDS("module_population_partitioning/data/pp_residual_ancestry.rds")
stopifnot("residual matrix must have exactly the 20,807 DI25 rho05 units" = ncol(rez$Resid) == 20807L)
pop_resid_mag <- data.table(Population = rownames(rez$Resid), resid_rms = sqrt(rowMeans(rez$Resid^2, na.rm = TRUE)))

geo <- readRDS("module_population_partitioning/data/followup/10_geographic_prediction.rds")
pop_pca <- geo$pca_scores[, .(pop, PC1, PC2)]; setnames(pop_pca, "pop", "Population")

fbar <- colMeans(Fmat, na.rm = TRUE)
pop_contrib <- data.table(Population = rownames(Fmat), contrib_diff = rowMeans(sweep(Fmat, 2, fbar, "-")^2, na.rm = TRUE))

targets <- Reduce(function(a, b) merge(a, b, by = "Population"), list(pop_runs, pop_influence, pop_resid_mag, pop_pca, pop_contrib))
stopifnot(nrow(targets) == 20)

metric_col <- "mean_frac_long_2cM"
test_targets <- c(mean_influence = "individual influence (Analysis 2)", resid_rms = "residual-profile magnitude (NOT independent -- same genotypes)",
                  PC1 = "geographic residual PC1 (Analysis 1)", PC2 = "geographic residual PC2 (Analysis 1)",
                  contrib_diff = "contribution to high-DI differentiation (NOT independent -- same genotypes)")
assoc <- rbindlist(lapply(names(test_targets), function(tc) {
  rho <- cor(targets[[metric_col]], targets[[tc]], method = "spearman")
  loo <- sapply(seq_len(20), function(i) cor(targets[[metric_col]][-i], targets[[tc]][-i], method = "spearman"))
  data.table(target = test_targets[tc], rho = rho, loo_min = min(loo), loo_max = max(loo))
}))
cat(sprintf("\n[run] population-level run metric (%s) vs the 4 required targets (n=20, Spearman + leave-one-pop-out range):\n", metric_col))
cat("    CAUTION: resid_rms and contrib_diff are computed from the SAME Fmat/Resid genotypes as the run\n")
cat("    metric itself -- these two associations are NOT independent corroboration, just shared measurement.\n")
print(assoc)

## ---------------------------------------------------------------------
## 6. figures
## ---------------------------------------------------------------------
suppressMessages(library(ggplot2))
theme_ms <- theme_bw(base_size = 12) + theme(strip.background = element_blank(), panel.grid.minor = element_blank())

multi_runs <- all_runs[state %in% c("aqu", "pol") & !is_single_marker]
fig_dist <- ggplot(multi_runs, aes(run_cM, fill = state)) +
  geom_histogram(position = "identity", alpha = 0.55, bins = 60) +
  scale_fill_manual(values = c(aqu = "#21918C", pol = "#D3C93B"), labels = c(aqu = "F. aquilonia", pol = "F. polyctena")) +
  labs(x = "run length (cM), multi-marker runs only", y = "count", fill = NULL,
      title = "Genotype-state run lengths (multi-marker runs; single-marker\nsegments excluded here -- see text for their separate count)") +
  theme_ms
ggsave(file.path(FIGDIR, "12_run_length_distribution.png"), fig_dist, width = 8, height = 5.5, dpi = 200)

indiv_runs[, Population := factor(Population, levels = pop_runs[order(mean_frac_long_2cM), Population])]
fig_by_pop <- ggplot(indiv_runs, aes(Population, frac_long_2cM)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") + geom_jitter(width = 0.15, size = 1.6, colour = "grey30") +
  labs(x = NULL, y = "fraction of (individual) callable genetic map in runs >=2cM",
      title = "Long-run fraction by population (exploratory genotype-state\nrun proxy, NOT an ancestry-tract inference)") +
  theme_ms + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(FIGDIR, "12_long_run_fraction_by_population.png"), fig_by_pop, width = 10, height = 5.5, dpi = 200)

assoc[, target := factor(target, levels = target)]
fig_assoc <- ggplot(assoc, aes(target, rho)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 3, colour = "#7570b3") + geom_errorbar(aes(ymin = loo_min, ymax = loo_max), width = 0.15, colour = "#7570b3") +
  labs(x = NULL, y = "Spearman rho (population-level, n=20)\npoint = full sample, bars = leave-one-population-out range",
      title = "Population-level long-run fraction vs the 4 required targets\n(NOT independent tests -- several share the same underlying genotypes)") +
  theme_ms + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(FIGDIR, "12_run_metric_associations.png"), fig_assoc, width = 8, height = 6.5, dpi = 200)
cat("\n[run] figures saved: 12_run_length_distribution.png, 12_long_run_fraction_by_population.png, 12_run_metric_associations.png\n")

## ---------------------------------------------------------------------
## 7. save
## ---------------------------------------------------------------------
cat("\n*** REQUIRED CAVEAT (not a finding) ***\n")
cat("This is a conservative, unphased GENOTYPE-STATE RUN PROXY, not a validated local-ancestry\n")
cat("tract call. It cannot distinguish recent backcrossing from old admixture followed by drift,\n")
cat("and predominantly short runs are NOT interpreted here as evidence against recent backcrossing\n")
cat("-- this proxy simply cannot resolve introgression timing either way. If an authoritative tract\n")
cat("caller becomes available, retain this analysis as an exploratory sensitivity check only.\n")

result <- list(all_runs = all_runs, all_runs_gap05_sensitivity = all_runs_gap05, all_runs_gap2_sensitivity = all_runs_gap2,
               all_runs_informativeness_sensitivity = all_runs_strict,
               panel_gap_distribution = panel_gaps, primary_max_gap_cM = 1.0,
               callable_cM_indiv = callable_cM_indiv,
               indiv_runs = indiv_runs, pop_runs = pop_runs, targets = targets, assoc = assoc,
               long_thresholds = LONG_THRESH,
               caveat = "conservative unphased GENOTYPE-STATE RUN PROXY (not 'ancestry tract'); cannot resolve introgression timing/recency; no authoritative local-ancestry caller found in this repository; population-level associations with resid_rms/contrib_diff are NOT independent (shared genotypes)",
               dataset_note = "uses the current rho05 (20,807-unit) dataset",
               session_info = sessionInfo(), run_time = Sys.time(), elapsed_secs = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
saveRDS(result, file.path(OUTDIR, "12_ancestry_run_proxy.rds"))
cat(sprintf("\n[run] saved -> %s (elapsed %.1fs)\n", file.path(OUTDIR, "12_ancestry_run_proxy.rds"), result$elapsed_secs))
