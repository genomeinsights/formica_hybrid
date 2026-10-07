## =========================================================================
## module_allele_specific_sorting -- EXPLORATORY: calibration grid of the neutral
## mosaic-founder model (phased founders; sim_founder_fix/slim/run_calib_job.sh).
##
## Each setting (K, initialN) = several independent runs (one hybrid population each),
## simulated on chromosomes 1-6 and sampled at several cycles. For every setting x cycle
## the runs are pooled as populations (N_PER females each) and the statistics of
## explore_mosaic_sims.R are computed (sim_stats_lib.R). The empirical data are
## restricted to the same chromosomes. Targets: near-neutral F_ST (drift) and the
## distance decay of within-population LD at DI25 units (time since admixture).
## Output: data/explore_calib_sims.rds, Figures/explore_calib_sims_*.png
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/explore_calib_sims.R <RESULTS_DIR> [N_CORES] [N_PER=10]
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
source("moduleA_sorting/R/parallelism_stats.R")          # classify_sort()
args <- commandArgs(trailingOnly = TRUE)
RES_DIR <- args[1]; N_CORES <- if (length(args) >= 2) as.integer(args[2]) else 6L
N_PER <- if (length(args) >= 3) as.integer(args[3]) else 10L
SEED <- 1L; CHRS <- paste0("Chr", 1:6)

## ---- marker sets, restricted to the simulated chromosomes ---------------------------------------
U  <- load_units(); u25 <- readRDS(file.path(OUT_DATA, "01_units.rds"))
GG <- load_oriented_genotypes(U)
inp <- readRDS("module_di25/data/di25_inputs.rds")
sgn_all <- sign(colMeans(inp$GTs_par[grepl("^Faqu", rownames(inp$GTs_par)), u25$unit_marker], na.rm = TRUE) -
                colMeans(inp$GTs_par[grepl("^Fpol", rownames(inp$GTs_par)), u25$unit_marker], na.rm = TRUE))
keep25 <- u25$Chr %in% CHRS
u25 <- u25[keep25]; sgn <- sgn_all[keep25]
neu <- readRDS(file.path(OUT_DATA, "06_lowDI.rds"))$units[Chr %in% CHRS, .(marker = rep_snp, Chr, Pos, cM_pos)]
e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
Gneu_emp <- e$GTs_with_parents[rownames(GG$G), neu$marker, drop = FALSE]; rm(e); invisible(gc())
sets <- list(DI25 = u25[, .(marker = unit_marker, Chr, Pos, cM_pos)], neutral = neu)
all_markers <- c(u25$unit_marker, neu$marker)
source("module_allele_specific_sorting/R/sim_stats_lib.R")
cat(sprintf("[calib] chromosomes %s: %d DI25 units, %d near-neutral units\n",
            paste(range(as.integer(sub("Chr", "", CHRS))), collapse = "-"), nrow(u25), nrow(neu)))

## ---- empirical (same chromosomes) ---------------------------------------------------------------
emp <- all_stats(inp$GTs_hyb[rownames(GG$G), u25$unit_marker], Gneu_emp, GG$pop, "empirical")
cat("[calib] empirical:\n"); print(emp$summary, digits = 3)

## ---- calibration runs ---------------------------------------------------------------------------
files <- list.files(RES_DIR, pattern = "^females_ckl[0-9]+_.*\\.vcf\\.gz$", full.names = TRUE)
fi <- data.table(file = files)
fi[, `:=`(cycle = as.integer(sub("^females_ckl([0-9]+)_.*", "\\1", basename(file))),
          tag = sub("^females_ckl[0-9]+_(.*)\\.vcf\\.gz$", "\\1", basename(file)))]
fi[, `:=`(setting = sub("_[0-9]+$", "", tag), run = as.integer(sub(".*_", "", tag)))]
cells <- fi[, .(n_runs = .N), by = .(setting, cycle)][n_runs >= 4][order(setting, cycle)]
cat("[calib] setting x cycle cells with >= 4 runs:\n"); print(cells)
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
res <- parallel::mclapply(seq_len(nrow(cells)), mc.cores = N_CORES, mc.set.seed = TRUE, FUN = function(i) {
  cl <- cells[i]; ff <- fi[setting == cl$setting & cycle == cl$cycle][order(run)]
  parts <- lapply(ff$file, read_sim, n = N_PER)
  pop <- rep(sprintf("r%02d", ff$run), sapply(parts, nrow)); G <- do.call(rbind, parts)
  s <- all_stats(G[, u25$unit_marker], G[, neu$marker], pop, cl$setting)
  s$summary[, `:=`(cycle = cl$cycle, n_runs = cl$n_runs)]; s$profile[, `:=`(cycle = cl$cycle)]
  s
})
cal_summary <- rbindlist(lapply(res, `[[`, "summary")); cal_prof <- rbindlist(lapply(res, `[[`, "profile"))
cat("\n[calib] simulated:\n"); print(cal_summary[order(data, cycle)], digits = 3)
saveRDS(list(empirical = emp, summary = cal_summary, profile = cal_prof, cells = cells, chroms = CHRS, n_per = N_PER),
        file.path(OUT_DATA, "explore_calib_sims.rds"))

## ---- figures ----------------------------------------------------------------------------------------
lab <- function(b) factor(ifelse(b == 99L, "unlinked", CM_LABELS[b]), levels = c(CM_LABELS, "unlinked"))
fs <- melt(cal_summary, id.vars = c("data", "cycle"), measure.vars = c("fst_DI25_pooled", "fst_neutral_pooled", "pct_sorted"), variable.name = "stat")
fe <- melt(emp$summary, id.vars = "data", measure.vars = c("fst_DI25_pooled", "fst_neutral_pooled", "pct_sorted"), variable.name = "stat")
p1 <- ggplot(fs, aes(cycle, value, colour = data)) + geom_line() + geom_point() +
  geom_hline(data = fe, aes(yintercept = value), linetype = 2) + facet_wrap(~ stat, scales = "free_y", nrow = 1) +
  scale_x_log10() + labs(x = "cycle (log)", y = NULL, colour = "setting", caption = "dashed: empirical (chromosomes 1-6)")
ggsave(file.path(OUT_FIG, "explore_calib_sims_summary.png"), p1, width = 12, height = 3.8, dpi = 200)

pd <- rbind(melt(cal_prof, id.vars = c("set", "bin", "data", "cycle"), measure.vars = c("share", "within_LD"), variable.name = "qty"),
            melt(emp$profile, id.vars = c("set", "bin", "data"), measure.vars = c("share", "within_LD"), variable.name = "qty")[, cycle := NA_integer_])
pd[, bin := lab(bin)][, curve := ifelse(data == "empirical", "empirical", paste("cycle", cycle))]
pd[, curve := factor(curve, levels = c(paste("cycle", sort(unique(cal_prof$cycle))), "empirical"))]
for (q in c("within_LD", "share")) {
  sims <- pd[qty == q & data != "empirical"]; em <- pd[qty == q & data == "empirical"]
  em <- em[, .(setting = unique(sims$data)), by = .(set, bin, value, curve)]
  p <- ggplot(sims[, setting := data], aes(bin, value, colour = curve, group = curve)) + geom_line() + geom_point(size = 0.8) +
    geom_line(data = em, colour = "black", linewidth = 0.9) + facet_grid(set ~ setting, scales = "free_y") +
    scale_colour_viridis_d(end = 0.9, name = NULL) + labs(x = NULL, y = q, caption = "black: empirical (chromosomes 1-6)") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  if (q == "within_LD") p <- p + scale_y_log10()
  ggsave(file.path(OUT_FIG, sprintf("explore_calib_sims_%s.png", q)), p, width = 13, height = 6, dpi = 200)
}
cat("[calib] done\n")
