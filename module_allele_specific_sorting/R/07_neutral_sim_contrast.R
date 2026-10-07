## =========================================================================
## module_allele_specific_sorting -- 07: neutral simulations, near-neutral vs
## ancestry-informative loci.
##
## Question: can a neutral hybrid model that reproduces the differentiation of
## near-neutral loci also reproduce that of ancestry-informative (DI25) loci?
##
## Simulations (sim_founder_fix/, SLiM; see sim_founder_fix/slim/run_calib_job.sh):
## haplodiploid model of Portinha et al., neutral (no QTN targets, no DMI effect),
## founders = genotype mosaics of the Beagle-phased empirical parents (DI25 panel +
## near-neutral SNPs), chromosomes 1-6, unscaled recombination map. Grid: carrying
## capacity K (6,250 / 12,500) x founding number (100 / 1,000), each run one
## independent hybrid population sampled at cycles 60, 125, 250, 500 and 1000. A
## grid cell = one setting x cycle, its runs pooled as populations.
##
## Statistics, identical code for empirical (chromosomes 1-6) and simulated data,
## separately for DI25 units and near-neutral units (DI <= -90; 06_lowDI.rds):
##   F_ST     : Weir & Cockerham, pooled over loci (sim_stats_lib.R::wc_fst)
##   % fixed  : share of unit x population combinations monomorphic in a sample of
##              N_FIX = 9 individuals (frequency <= 0.05 or >= 0.95, i.e. 0 or 1 at
##              n = 9); every population is downsampled to 9 (smaller ones skipped)
##              so the statistic does not depend on sample size
##   excess within-population LD at 0.05-0.2 cM (ancestry-adjusted, both loci
##              segregating; minus the between-chromosome level), reported to
##              document the model's LD limitation, not used as evidence
## Best-fitting cell = the one whose near-neutral F_ST is closest (log ratio) to
## the empirical value.
## Output: data/07_neutral_sim.rds, Figures/07_neutral_sim.{pdf,png}
## Run from the formica_hybrid repo root (after 01 and 06):
##   Rscript module_allele_specific_sorting/R/07_neutral_sim_contrast.R [RESULTS_DIR=sim_founder_fix/out/calib] [N_CORES=6]
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
args <- commandArgs(trailingOnly = TRUE)
RES_DIR <- if (length(args) >= 1) args[1] else "sim_founder_fix/out/calib"
N_CORES <- if (length(args) >= 2) as.integer(args[2]) else 6L
SETTINGS <- c("K6k_N100", "K6k_N1000", "K12k_N100", "K12k_N1000")   # unscaled map
SET_LAB <- c(K6k_N100 = "K 6,250, 100 founders", K6k_N1000 = "K 6,250, 1,000 founders",
             K12k_N100 = "K 12,500, 100 founders", K12k_N1000 = "K 12,500, 1,000 founders")
CHRS <- paste0("Chr", 1:6); N_PER <- 10L; N_FIX <- 9L; MIN_RUNS <- 4L; SEED <- 1L

## ---- marker sets (chromosomes 1-6) and empirical data --------------------------------------
U <- load_units(); u25 <- readRDS(file.path(OUT_DATA, "01_units.rds"))
GG <- load_oriented_genotypes(U)
inp <- readRDS("module_di25/data/di25_inputs.rds")
sgn_all <- sign(colMeans(inp$GTs_par[grepl("^Faqu", rownames(inp$GTs_par)), u25$unit_marker], na.rm = TRUE) -
                colMeans(inp$GTs_par[grepl("^Fpol", rownames(inp$GTs_par)), u25$unit_marker], na.rm = TRUE))
k25 <- u25$Chr %in% CHRS; u25 <- u25[k25]; sgn <- sgn_all[k25]
neu <- readRDS(file.path(OUT_DATA, "06_lowDI.rds"))$units[Chr %in% CHRS, .(marker = rep_snp, Chr, Pos, cM_pos)]
e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
Gneu_emp <- e$GTs_with_parents[rownames(GG$G), neu$marker, drop = FALSE]; rm(e); invisible(gc())
sets <- list(DI25 = u25[, .(marker = unit_marker, Chr, Pos, cM_pos)], neutral = neu)
all_markers <- c(u25$unit_marker, neu$marker)
source("module_allele_specific_sorting/R/sim_stats_lib.R")   # wc_fst, profile, read_sim (globals above)
cat(sprintf("[07] chromosomes 1-6: %d DI25 units, %d near-neutral units\n", nrow(u25), nrow(neu)))

pct_fixed <- function(G, pop, markers) {
  keep <- names(which(table(pop) >= N_FIX))
  Fm <- t(sapply(keep, function(p) { r <- which(pop == p); r <- r[sample.int(length(r), N_FIX)]
    colMeans(G[r, markers, drop = FALSE], na.rm = TRUE) / 2 }))
  100 * mean(Fm <= 0.05 | Fm >= 0.95, na.rm = TRUE)
}
stats <- function(G25raw, Gneu, pop) {
  Gor <- G25raw; Gor[, which(sgn < 0)] <- 2L - Gor[, which(sgn < 0)]; Gor[, which(is.na(sgn) | sgn == 0)] <- NA
  ld <- function(prof) prof[bin == 4L, within_LD] - prof[bin == 99L, within_LD]          # 0.05-0.2 cM minus unlinked
  G <- cbind(G25raw, Gneu)
  data.table(partition = c("DI25", "neutral"),
             fst = c(wc_fst(G25raw, pop)$pooled, wc_fst(Gneu, pop)$pooled),
             pct_fixed = c(pct_fixed(G, pop, u25$unit_marker), pct_fixed(G, pop, neu$marker)),
             excess_ld = c(ld(profile(Gor, pop, sets$DI25, Gor)), ld(profile(Gneu, pop, sets$neutral, Gor))))
}
set.seed(SEED)
emp <- stats(inp$GTs_hyb[rownames(GG$G), u25$unit_marker], Gneu_emp, GG$pop)[, `:=`(setting = "empirical", cycle = NA_integer_, n_runs = NA_integer_)]
cat("[07] empirical:\n"); print(emp, digits = 3)

## ---- simulations ----------------------------------------------------------------------------
fi <- data.table(file = list.files(RES_DIR, "^females_ckl[0-9]+_.*\\.vcf\\.gz$", full.names = TRUE))
fi[, `:=`(cycle = as.integer(sub(".*ckl([0-9]+)_.*", "\\1", basename(file))),
          setting = sub("^females_ckl[0-9]+_(.*)_[0-9]+\\.vcf\\.gz$", "\\1", basename(file)),
          run = as.integer(sub(".*_([0-9]+)\\.vcf\\.gz$", "\\1", basename(file))))]
fi <- fi[setting %in% SETTINGS]
cells <- fi[, .(n_runs = .N), by = .(setting, cycle)][n_runs >= MIN_RUNS][order(setting, cycle)]
cat(sprintf("[07] %d grid cells (setting x cycle) with >= %d runs\n", nrow(cells), MIN_RUNS))
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
sim <- rbindlist(parallel::mclapply(seq_len(nrow(cells)), mc.cores = N_CORES, mc.set.seed = TRUE, FUN = function(i) {
  cl <- cells[i]; ff <- fi[setting == cl$setting & cycle == cl$cycle][order(run)]
  parts <- lapply(ff$file, read_sim, n = N_PER)
  pop <- rep(sprintf("r%02d", ff$run), sapply(parts, nrow)); G <- do.call(rbind, parts)
  stats(G[, u25$unit_marker], G[, neu$marker], pop)[, `:=`(setting = cl$setting, cycle = cl$cycle, n_runs = cl$n_runs)]
}))
res <- rbind(emp, sim)

## ---- best-fitting cell and grid range ---------------------------------------------------------
fit <- sim[partition == "neutral", .(setting, cycle, d = abs(log(fst / emp[partition == "neutral", fst])))][which.min(d)]
best <- sim[setting == fit$setting & cycle == fit$cycle]
rng <- sim[, .(fst_min = min(fst), fst_max = max(fst), fix_min = min(pct_fixed), fix_max = max(pct_fixed),
               ld_min = min(excess_ld), ld_max = max(excess_ld)), by = partition]
summ <- merge(merge(emp[, .(partition, emp_fst = fst, emp_fixed = pct_fixed, emp_ld = excess_ld)],
                    best[, .(partition, best_fst = fst, best_fixed = pct_fixed, best_ld = excess_ld)], by = "partition"), rng, by = "partition")
cat(sprintf("\n[07] best-fitting cell (near-neutral F_ST): %s, cycle %d (%d runs)\n", fit$setting, fit$cycle, best$n_runs[1]))
print(summ, digits = 3)
saveRDS(list(result = res, summary = summ, best = fit[, .(setting, cycle)], cells = cells, chroms = CHRS,
             n_per = N_PER, n_fix = N_FIX, settings = SET_LAB), file.path(OUT_DATA, "07_neutral_sim.rds"))

## ---- figure: F_ST and % fixed by cycle, per partition ----------------------------------------
pd <- melt(sim, id.vars = c("partition", "setting", "cycle"), measure.vars = c("fst", "pct_fixed"), variable.name = "stat")
pe <- melt(emp, id.vars = "partition", measure.vars = c("fst", "pct_fixed"), variable.name = "stat")
lab_part <- c(DI25 = "ancestry-informative (DI > -25)", neutral = "near-neutral (DI <= -90)")
lab_stat <- c(fst = "F[ST]", pct_fixed = "'% unit x population fixed'")
for (d in list(pd, pe)) d[, `:=`(partition = factor(lab_part[partition], levels = lab_part),
                                 stat = factor(lab_stat[as.character(stat)], levels = lab_stat))]
pd[, setting := factor(SET_LAB[setting], levels = SET_LAB)]
p <- ggplot(pd, aes(cycle, value, colour = setting)) + geom_line() + geom_point(size = 1.2) +
  geom_hline(data = pe, aes(yintercept = value), linetype = 2) +
  facet_grid(stat ~ partition, scales = "free_y", labeller = labeller(stat = label_parsed, partition = label_value)) +
  scale_x_log10(breaks = c(60, 125, 250, 500, 1000)) + scale_colour_viridis_d(end = 0.85, name = NULL) +
  labs(x = "simulation cycles since hybridisation (log scale)", y = NULL) + theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG, "07_neutral_sim.pdf"), p, width = 8, height = 5.5)
ggsave(file.path(OUT_FIG, "07_neutral_sim.png"), p, width = 8, height = 5.5, dpi = 200)
cat("[07] done\n")
