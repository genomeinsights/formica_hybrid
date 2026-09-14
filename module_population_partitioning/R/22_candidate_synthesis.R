## =========================================================================
## module_population_partitioning -- candidate-locus follow-up 22: decision
## table, provenance table, concise results summary, and one composite
## synthesis figure for the whole candidate-locus population-frequency and
## profile-similarity analysis (scripts 15-21).
##
## Run from the formica_hybrid repo root, after 21_climate_bdmi_context.R:
##   Rscript module_population_partitioning/R/22_candidate_synthesis.R
## Reads : module_population_partitioning/data/followup/15-21_*.rds
## Writes: module_population_partitioning/data/followup/22_candidate_synthesis.rds
##         module_population_partitioning/Figures/followup/22_synthesis.png
## =========================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
t_start <- Sys.time()
OUTDIR <- "module_population_partitioning/data/followup"
FIGDIR <- "module_population_partitioning/Figures/followup"
r15 <- readRDS(file.path(OUTDIR, "15_candidate_data.rds"))
manifest <- readRDS(file.path(OUTDIR, "15_manifest.rds"))
r16 <- readRDS(file.path(OUTDIR, "16_structure_reference.rds"))
r18 <- readRDS(file.path(OUTDIR, "18_similarity_raw.rds"))
r19 <- readRDS(file.path(OUTDIR, "19_structure_adjustment.rds"))
r20 <- readRDS(file.path(OUTDIR, "20_matched_null.rds"))
r21 <- readRDS(file.path(OUTDIR, "21_climate_bdmi_context.rds"))
theme_ms <- theme_bw(base_size = 11) + theme(strip.background = element_blank(), panel.grid.minor = element_blank())

## ---------------------------------------------------------------------
## 1. provenance table (near-direct transcription of the script-15 manifest)
## ---------------------------------------------------------------------
cat("=== provenance table ===\n"); print(manifest[, .(path, role, n_rows_read, status)], width = 200)

## ---------------------------------------------------------------------
## 2. candidate counts before/after mapping, per target
## ---------------------------------------------------------------------
counts_tab <- rbindlist(lapply(names(r15$TARGETS), function(tgt) {
  data.table(target = tgt, n_raw_candidates = length(r15$cand_ids[[tgt]]),
            n_floor_survivors = length(r15$floor_ids[[tgt]]),
            n_undef_orientation_in_raw = r15$undef_counts[target == tgt & status == "raw", n_undef],
            is_discovery_set = r15$TARGETS[[tgt]]$is_discovery,
            no_discovery_set_constructed = isTRUE(r15$TARGETS[[tgt]]$no_discovery_set))
}))
cat("\n=== candidate counts before/after mapping ===\n"); print(counts_tab)

## ---------------------------------------------------------------------
## 3. decision table: matched-null verdict per target (the central
##    inferential result -- is candidate long-range similarity greater than
##    matched genomic background, and does it survive structure adjustment?)
## ---------------------------------------------------------------------
decision_table <- rbindlist(lapply(names(r20$matched_null_results), function(tgt) {
  st <- r20$matched_null_results[[tgt]]$summary_tab
  cross <- st[statistic == "cross-chr mean|r|"]
  sd_row <- st[statistic == "structure_dominance"]
  ctx <- r21$context_results[[tgt]]
  clim_explained <- if (!is.null(ctx) && nrow(ctx$clim_strat) == 2)
    sprintf("both-climate-strong pairs mean|r|=%.3f vs %.3f otherwise", ctx$clim_strat[both_climate_strong == TRUE, mean_absr], ctx$clim_strat[both_climate_strong == FALSE, mean_absr])
  else if (!is.null(ctx)) "no long-range pairs with both loci strongly climate-correlated (climate does not explain this target's signal)"
  else "not evaluated"
  data.table(target = tgt,
            cross_chr_obs = round(cross$obs_raw, 3), cross_chr_null_mean = round(cross$null_mean_raw, 3), cross_chr_p = cross$p_raw,
            cross_chr_obs_adjusted = round(cross$obs_adj, 3), cross_chr_p_adjusted = cross$p_adj,
            structure_dominance_obs = round(sd_row$obs_raw, 3), structure_dominance_null_mean = round(sd_row$null_mean_raw, 3), structure_dominance_p = sd_row$p_raw,
            climate_context = clim_explained)
}))
cat("\n=== decision table ===\n"); print(decision_table, width = 300)

## ---------------------------------------------------------------------
## 4. verdict text
## ---------------------------------------------------------------------
verdict <- paste(
 "Across all four independently-scanned targets (PC1, PC2, bio_winter, mitoC2), candidate",
 "loci's population-frequency profiles show substantially GREATER cross-chromosome and",
 "same-chromosome-distant similarity than matched Stage-1-direct genomic units of comparable",
 "DI, parental MAF, recombination rate, cluster size and (where possible) chromosome",
 "(matched-set bootstrap, B=1000; observed values fall outside or at the extreme edge of the",
 "null distribution for every target's cross-chromosome and structure-dominance statistics).",
 "This elevated long-range similarity is NOT removed by adjusting for the single leading axis",
 "of genome-wide population structure estimated from the canonical, non-candidate reference",
 "panel (raw and structure-adjusted results are close throughout, consistent with that",
 "reference PC1 explaining only a modest fraction, median R2 well under 0.1, of any given",
 "candidate locus's population variance). Per target, the leading-eigenvalue fraction of the",
 "locus-by-locus correlation matrix (\"structure_dominance\") is elevated relative to matched",
 "background but stays well below 1 in every case (roughly 0.32-0.43 observed vs a matched-null",
 "mean near 0.13-0.22) -- candidate profiles are neither a single dominant genome-wide division",
 "(structure_dominance far from 1) nor fully independent of each other (elevated relative to a",
 "size-matched genomic background); they form several partly-overlapping, reproducible groups.",
 "The bio_winter target's long-range similarity is substantially explained by loci sharing",
 "correlation with the same population-level winter-climate axis (the target the scan itself",
 "selected on -- an expected, not surprising, pattern); PC2 and mitoC2 show essentially no",
 "climate-correlated long-range pairs, so their elevated similarity is NOT explained by climate",
 "and remains an open question. BDMI-region overlap does not show a consistent enrichment",
 "pattern for long-range pair similarity across targets (elevated for PC1, roughly flat for",
 "PC2, and if anything reduced for bio_winter) -- this is reported as a genuinely mixed result,",
 "not evidence for or against BDMI involvement, and per the brief's own instruction this",
 "distant-correlation pattern is NOT, by itself, interpreted as evidence of epistasis."
)
cat("\n=== verdict ===\n"); cat(strwrap(verdict, width = 78), sep = "\n")

## ---------------------------------------------------------------------
## 5. limitations
## ---------------------------------------------------------------------
limitations <- paste(
 "PC1 and PC2 raw candidate sets are NOT described as discoveries anywhere in this analysis",
 "(they are the full BF>=15 threshold-crossing sets from an exploratory continuous-covariate",
 "scan, not a null-calibrated result -- floor-survivor counts are only 1 and 2 respectively).",
 "No mitoC2 discovery set was constructed (0 floor survivors from the current, audited null",
 "calibration; only its 11 raw C2 threshold-crossings are analysed descriptively). A",
 "substantial fraction of candidates in every target (3-36 depending on target/status) have",
 "UNDEFINED ancestry orientation (fixed for the same allele in both parental species) and are",
 "excluded from the similarity/matched-null analyses (reported explicitly, not silently",
 "dropped) -- these loci's population differentiation is real but its direction relative to",
 "the aquilonia/polyctena axis cannot be assessed with this data. Same-chromosome-near",
 "(<=100kb) pair counts are very small for most target/status combinations (often <10), so the",
 "physical-distance decay slope is frequently NOT ESTIMABLE and is reported as such rather than",
 "computed on a handful of pairs. The DI used throughout is module_manuscript_rho05's",
 "full-genome, UNGATED vintage, not this module's own DI25-ascertained vintage used elsewhere",
 "-- the two can disagree at the margin and must not be conflated. The structure-adjustment",
 "reference panel excludes 2,031/17,509 units for missingness and, like the candidate side,",
 "excludes Aland; this is a design choice (documented), not incidental. BDMI overlap uses a",
 "single documented default cutoff (index 13); the rotation-null enrichment test available in",
 "module_di25 was not run here since BDMI status is used only as a stratifying annotation, per",
 "the brief's own framing, not as a standalone enrichment claim."
)
cat("\n=== limitations ===\n"); cat(strwrap(limitations, width = 78), sep = "\n")

## ---------------------------------------------------------------------
## 6. composite synthesis figure (4 panels)
## ---------------------------------------------------------------------
## A: reference-panel PC1 variance explained (is there a single dominant axis genome-wide?)
ve_df <- data.table(PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)), ve = r16$ve[1:10])
panel_A <- ggplot(ve_df, aes(PC, ve)) + geom_col(fill = "grey50") +
  labs(title = "A. genome-wide reference-panel\nvariance explained", x = NULL, y = "% variance") + theme_ms

## B: candidate vs matched-null cross-chr similarity, all 4 targets
crossb <- decision_table[, .(target, obs = cross_chr_obs, null_mean = cross_chr_null_mean)]
crossb_long <- melt(crossb, id.vars = "target", variable.name = "which", value.name = "mean_absr")
panel_B <- ggplot(crossb_long, aes(target, mean_absr, fill = which)) + geom_col(position = "dodge") +
  scale_fill_manual(values = c(obs = "firebrick", null_mean = "grey60"), labels = c("candidate", "matched-null mean")) +
  labs(title = "B. cross-chr similarity: candidate\nvs matched-null", x = NULL, y = "mean |r|", fill = NULL) + theme_ms +
  theme(legend.position = "bottom")

## C: structure_dominance, candidate vs matched-null, all 4 targets
sdb <- decision_table[, .(target, obs = structure_dominance_obs, null_mean = structure_dominance_null_mean)]
sdb_long <- melt(sdb, id.vars = "target", variable.name = "which", value.name = "value")
panel_C <- ggplot(sdb_long, aes(target, value, fill = which)) + geom_col(position = "dodge") +
  scale_fill_manual(values = c(obs = "firebrick", null_mean = "grey60"), labels = c("candidate", "matched-null mean")) +
  labs(title = "C. structure_dominance: candidate\nvs matched-null", x = NULL, y = "leading-eigenvalue fraction", fill = NULL) + theme_ms +
  theme(legend.position = "bottom")

## D: raw vs structure-adjusted cross-chr similarity (does adjustment remove the signal?)
adjd <- decision_table[, .(target, raw = cross_chr_obs, adjusted = cross_chr_obs_adjusted)]
adjd_long <- melt(adjd, id.vars = "target", variable.name = "which", value.name = "mean_absr")
panel_D <- ggplot(adjd_long, aes(target, mean_absr, fill = which)) + geom_col(position = "dodge") +
  scale_fill_manual(values = c(raw = "#21918C", adjusted = "#D3C93B")) +
  labs(title = "D. cross-chr similarity: raw vs\nstructure-adjusted", x = NULL, y = "mean |r|", fill = NULL) + theme_ms +
  theme(legend.position = "bottom")

fig_synth <- (panel_A + panel_B) / (panel_C + panel_D)
ggsave(file.path(FIGDIR, "22_synthesis.png"), fig_synth, width = 10, height = 9, dpi = 200)
cat("\n[synth] figure saved: 22_synthesis.png\n")

## ---------------------------------------------------------------------
## 7. save
## ---------------------------------------------------------------------
result <- list(manifest = manifest, counts_tab = counts_tab, decision_table = decision_table,
              verdict = verdict, limitations = limitations,
              session_info = sessionInfo(), run_time = Sys.time(),
              elapsed_secs = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
saveRDS(result, file.path(OUTDIR, "22_candidate_synthesis.rds"))
cat(sprintf("\n[synth] saved -> %s (elapsed %.1fs)\n", file.path(OUTDIR, "22_candidate_synthesis.rds"), result$elapsed_secs))
