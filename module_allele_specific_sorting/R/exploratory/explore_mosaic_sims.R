## =========================================================================
## module_allele_specific_sorting -- EXPLORATORY: empirical data vs the NEW neutral
## simulations with mosaic founders (sim_founder_fix/slim/; run on mini1/mini2).
##
## Replicate k = run k of every founding group (paired populations lanR/lanW and
## bun/grund share founders, as in the original design). Each population's hybrid
## sample (new queens at cycle 125; Sielva at cycle 11) is downsampled to the
## empirical sample size. Simulated alleles use the empirical allele coding (the
## founders copy empirical genotypes), so the empirical parental orientation applies.
## Statistics, computed with IDENTICAL code for the empirical data and each replicate:
##   F_ST (W&C) at DI25 units and at near-neutral units (DI <= -90; calibration anchor)
##   % directionally sorted DI25 units (tau 0.6, phi 0.85, binom alpha 0.05)
##   Figure-1 profile by cM bin + unlinked: ceiling, realised fraction, share beyond
##   1/(n_pop - 1), ancestry-adjusted within-population LD -- for DI25 and near-neutral units
## Output: data/explore_mosaic_sims.rds, Figures/explore_mosaic_sims_*.png
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/exploratory/explore_mosaic_sims.R <RESULTS_DIR> [N_CORES]
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
source("moduleA_sorting/R/parallelism_stats.R")          # classify_sort()
args <- commandArgs(trailingOnly = TRUE)
RES_DIR <- args[1]; N_CORES <- if (length(args) >= 2) as.integer(args[2]) else 6L
SEED <- 1L

## empirical sample sizes per simulated population code (Beatriz's bootstrap TARGETS)
TARGETS <- c(aland = 10, katis = 5, lanR = 9, lanW = 20, svan1 = 3, svan2 = 7, tvar = 10, bun = 10, grund = 10,
             pik = 10, nyr1 = 5, nyr2 = 5, heina = 10, pari = 6, hiiv = 6, vuos = 10, kumm = 10, karsi = 10,
             jarven = 4, sielva = 5)

## ---- marker sets ----------------------------------------------------------------------------
U  <- load_units(); u25 <- readRDS(file.path(OUT_DATA, "01_units.rds"))
GG <- load_oriented_genotypes(U)                                  # empirical DI25, oriented
inp <- readRDS("module_di25/data/di25_inputs.rds")
sgn <- sign(colMeans(inp$GTs_par[grepl("^Faqu", rownames(inp$GTs_par)), u25$unit_marker], na.rm = TRUE) -
            colMeans(inp$GTs_par[grepl("^Fpol", rownames(inp$GTs_par)), u25$unit_marker], na.rm = TRUE))
neu <- readRDS(file.path(OUT_DATA, "06_lowDI.rds"))$units[, .(marker = rep_snp, Chr, Pos, cM_pos)]
e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
sdp <- e$sample_data_with_parents
Gneu_emp <- e$GTs_with_parents[rownames(GG$G), neu$marker, drop = FALSE]; rm(e); invisible(gc())
sets <- list(DI25 = u25[, .(marker = unit_marker, Chr, Pos, cM_pos)], neutral = neu)

source("module_allele_specific_sorting/R/sim_stats_lib.R")   # wc_fst, profile, sorting_pct, all_stats, read_sim

## ---- empirical -------------------------------------------------------------------------------------
G25_emp_raw <- inp$GTs_hyb[rownames(GG$G), u25$unit_marker]
emp <- all_stats(G25_emp_raw, Gneu_emp, GG$pop, "empirical")
cat("[mosaic] empirical summary:\n"); print(emp$summary, digits = 3)

## ---- simulations -----------------------------------------------------------------------------------
files <- list.files(RES_DIR, pattern = "^females_ckl[0-9]+_.*\\.vcf\\.gz$", full.names = TRUE)
fi <- data.table(file = files, tag = sub("^females_ckl[0-9]+_(.*)\\.vcf\\.gz$", "\\1", basename(files)))
fi[, `:=`(pop = sub("_[0-9]+$", "", tag), run = as.integer(sub(".*_", "", tag)))]
runs <- fi[, .(n_pops = uniqueN(pop)), by = run][n_pops == length(TARGETS), run]
cat(sprintf("[mosaic] %d result files; complete replicates (all 20 populations): %d\n", nrow(fi), length(runs)))
all_markers <- c(u25$unit_marker, neu$marker)
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
sim <- parallel::mclapply(runs, mc.cores = N_CORES, mc.set.seed = TRUE, FUN = function(k) {
  parts <- lapply(names(TARGETS), function(p) read_sim(fi[pop == p & run == k, file], TARGETS[[p]]))
  pop <- rep(names(TARGETS), sapply(parts, nrow)); G <- do.call(rbind, parts)
  all_stats(G[, u25$unit_marker], G[, neu$marker], pop, sprintf("sim run %02d", k))
})
sim_summary <- rbindlist(lapply(sim, `[[`, "summary")); sim_prof <- rbindlist(lapply(sim, `[[`, "profile"))
cat("\n[mosaic] simulated replicates:\n"); print(sim_summary, digits = 3)
saveRDS(list(empirical = emp, sim_summary = sim_summary, sim_profile = sim_prof, runs = runs),
        file.path(OUT_DATA, "explore_mosaic_sims.rds"))

## ---- figures ----------------------------------------------------------------------------------------------
lab <- function(b) factor(ifelse(b == 99L, "unlinked", CM_LABELS[b]), levels = c(CM_LABELS, "unlinked"))
sp <- melt(sim_prof, id.vars = c("set", "bin", "data"), measure.vars = c("ceiling", "share", "within_LD"), variable.name = "qty")[
  , .(mean = mean(value, na.rm = TRUE), lo = quantile(value, 0.025, na.rm = TRUE), hi = quantile(value, 0.975, na.rm = TRUE)), by = .(set, bin, qty)][
  , data := sprintf("neutral simulations, mosaic founders (%d replicates)", length(runs))]
ep <- melt(emp$profile, id.vars = c("set", "bin", "data"), measure.vars = c("ceiling", "share", "within_LD"), variable.name = "qty")[
  , .(set, bin, qty, mean = value, lo = NA_real_, hi = NA_real_, data = "empirical")]
pd <- rbind(ep, sp)[, bin := lab(bin)]
pd[, qty := factor(qty, levels = c("ceiling", "share", "within_LD"),
                   labels = c("ceiling (max among-pop LD, log)", "share of ceiling beyond chance", "within-population LD (log)"))]
pd[, set := factor(set, levels = c("DI25", "neutral"), labels = c("ancestry-informative (DI > -25)", "near-neutral (DI <= -90)"))]
mk <- function(q, logy) {
  p <- ggplot(pd[qty == q], aes(bin, mean, colour = data, group = data)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = data), colour = NA, alpha = 0.2, show.legend = FALSE) +
    geom_line() + geom_point(size = 1.1) + facet_wrap(~ set, nrow = 1) +
    scale_colour_manual(values = c("black", "#3182bd"), aesthetics = c("colour", "fill"), name = NULL) +
    labs(x = NULL, y = NULL, title = q) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom", plot.title = element_text(size = 10))
  if (logy) p + scale_y_log10() else p
}
p <- patchwork::wrap_plots(mk(levels(pd$qty)[1], TRUE), mk(levels(pd$qty)[2], FALSE), mk(levels(pd$qty)[3], TRUE), ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG, "explore_mosaic_sims_profile.png"), p, width = 10, height = 11, dpi = 200)

ss <- melt(sim_summary, id.vars = "data", variable.name = "stat")
es <- melt(emp$summary, id.vars = "data", variable.name = "stat")
p2 <- ggplot(ss, aes(stat, value)) + geom_boxplot(outlier.size = 0.6, colour = "#3182bd") +
  geom_point(data = es, colour = "black", shape = 18, size = 3.5) + facet_wrap(~ stat, scales = "free", nrow = 1) +
  labs(x = NULL, y = NULL, caption = "boxes: simulated replicates; black diamond: empirical") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(OUT_FIG, "explore_mosaic_sims_summary.png"), p2, width = 13, height = 3.6, dpi = 200)
cat("[mosaic] done\n")
