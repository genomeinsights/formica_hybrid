## =============================================================================
## module_di25 -- ancestry sorting: empirical vs the 1000-replicate simulated NULL
## =============================================================================
## Estimates per-SNP ancestry sorting EXACTLY as di25_sorting.R does (Module A
## conventions: ohta_fast_prepare -> parallelism_stats, fix_th=0.15 / phi=0.85,
## sort_rule="binom", alpha=0.05, orient by parents, no differentiation re-gate),
## but on each of the 1000 DIEM bootstrap replicates
## (data/diem_outs_demo/diem_boot<N>_output.bed), and compares the OBSERVED
## empirical sorting to the simulated null distribution.
##
## The demo sim is a matched design: 20 hybrid populations (near-identical sizes)
## + 15 aq + 15 pol parents, so sorting is estimated per population identically.
## DIEM 0/1/2 states are used as dosage; orientation is by parents, so DIEM
## polarity is irrelevant. Metrics per tau in {0.5,0.6,0.7,0.8}:
##   pct_sorted, toward_aqu, toward_pol, pct_aqu_of_resolved.
##
## Per-SNP only here (the exact observed-data estimator). Per-eMLG would require
## re-clustering each replicate (heavy) -- offered as a follow-up.
##
## Resumable + parallel (per-rep cache module_di25/data/sorting_sim_cache/rep<N>.rds).
##   Rscript module_di25/R/di25_sorting_emp_vs_sim.R [nreps] [workers]
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(parallel); library(patchwork) })
devtools::load_all("~/gitlab/LDscnR/")
source("moduleA_sorting/R/parallelism_stats.R")

args    <- commandArgs(trailingOnly = TRUE)
NREPS   <- if (length(args) >= 1) as.integer(args[1]) else 1000L
WORKERS <- if (length(args) >= 2) as.integer(args[2]) else 9L
BEDFMT    <- "data/diem_outs_demo/diem_boot%d_output.bed"
CACHE_DIR <- "module_di25/data/sorting_sim_cache"
SWEEP_EMP <- "module_di25/data/di25_sorting_sweep.rds"
OUTRDS    <- "module_di25/data/di25_sorting_emp_vs_sim.rds"
OUTPNG    <- "module_di25/Figures/di25_sorting_emp_vs_sim.png"
OUTPDF    <- sub("\\.png$", ".pdf", OUTPNG)
TAU_GRID  <- c(0.5, 0.6, 0.7, 0.8)
FIX_TH    <- 0.15; SORT_RULE <- "binom"; ALPHA <- 0.05
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- tally sorting for one prepared replicate ------------------------------
tally_reps <- function(ps, REP) {
  base <- ps[differentiated == TRUE & n_obs > 0]
  n_diff <- nrow(base)
  rbindlist(lapply(TAU_GRID, function(tau) {
    cls <- classify_sort(base$n_aqu, base$n_pol, base$n_obs,
                         sort_th = tau, sort_rule = SORT_RULE, alpha = ALPHA)
    n_aqu <- sum(cls == "aquilonia"); n_pol <- sum(cls == "polyctena")
    n_sorted <- sum(cls != "unsorted")
    data.table(rep = REP, tau = tau, n_diff = n_diff, n_sorted = n_sorted,
               pct_sorted = 100 * n_sorted / n_diff,
               toward_aqu = n_aqu, toward_pol = n_pol,
               pct_aqu_of_resolved = 100 * n_aqu / (n_aqu + n_pol))
  }))
}

## ---- per-replicate sorting (cached) ----------------------------------------
## CAVEAT (this null is PRELIMINARY):
##  (1) Ascertainment NOT matched -- every simulated marker enters the null; the
##      empirical DI > -25 diagnostic ascertainment is not recreated on the
##      simulations, so the null panel is not ascertainment-matched to the observed.
##  (2) Cache keyed by REPLICATE INDEX ONLY -- existence-based, not parameter-
##      stamped. If any sorting parameter or the empirical panel changes, delete
##      CACHE_DIR to force a recompute; the cache cannot detect the change itself.
process_rep <- function(REP) {
  cache <- file.path(CACHE_DIR, sprintf("rep%d.rds", REP))
  if (file.exists(cache)) return(invisible(TRUE))
  bed <- sprintf(BEDFMT, REP); if (!file.exists(bed)) return(invisible(FALSE))
  out <- tryCatch({
    hdr  <- readLines(bed, n = 2)
    inds <- strsplit(strsplit(hdr[2], "\t")[[1]][10], "\\|")[[1]]
    sim  <- fread(bed, skip = 2, header = FALSE, sep = "\t",
                  select = c(1, 3, 10), colClasses = list(character = c(1, 10)),
                  showProgress = FALSE)
    markers <- paste0("Chr", sub("ch", "", sim$V1), ":", sim$V3)
    geno <- sub("^S", "", sim$V10)
    S <- matrix(unlist(strsplit(geno, "", fixed = TRUE), use.names = FALSE),
                nrow = length(geno), byrow = TRUE)
    dos <- matrix(NA_integer_, nrow(S), ncol(S))
    dos[S == "0"] <- 0L; dos[S == "1"] <- 1L; dos[S == "2"] <- 2L
    ## individuals x markers dosage; marker names on cols, population labels from names
    G <- t(dos); colnames(G) <- markers               # inds x markers
    pops <- ifelse(grepl("^aq_",  inds), "aquilonia_parent",
             ifelse(grepl("^pol_", inds), "polyctena_parent",
                    sub("^hyb_(.*)_[0-9]+$", "\\1", inds)))
    hybrid_pops <- setdiff(unique(pops), c("aquilonia_parent", "polyctena_parent"))
    prep <- ohta_fast_prepare(G, pops = pops)
    ps <- parallelism_stats(prep, hybrid_pops = hybrid_pops,
                            aqu_pops = "aquilonia_parent", pol_pops = "polyctena_parent",
                            fix_th = FIX_TH, DI = NULL, min_DI = NULL,
                            parent_maf = NULL, min_parent_maf = NULL,
                            sort_rule = SORT_RULE, alpha = ALPHA)
    saveRDS(tally_reps(ps, REP), cache); TRUE
  }, error = function(e) { message(sprintf("  rep%d FAILED: %s", REP, conditionMessage(e))); FALSE })
  invisible(out)
}

reps <- seq_len(NREPS)
todo <- reps[!file.exists(file.path(CACHE_DIR, sprintf("rep%d.rds", reps)))]
message(sprintf("[sorting] %d reps requested, %d cached, %d to compute (workers=%d)",
                NREPS, NREPS - length(todo), length(todo), WORKERS))
if (length(todo)) invisible(mclapply(todo, process_rep, mc.cores = WORKERS, mc.preschedule = FALSE))

## ---- aggregate the null + compare to observed ------------------------------
done <- reps[file.exists(file.path(CACHE_DIR, sprintf("rep%d.rds", reps)))]
sim  <- rbindlist(lapply(done, function(r) readRDS(file.path(CACHE_DIR, sprintf("rep%d.rds", r)))))
emp  <- as.data.table(readRDS(SWEEP_EMP))[level == "SNP",
          .(tau, pct_sorted, toward_aqu, toward_pol, pct_aqu_of_resolved)]

metrics <- c("pct_sorted", "pct_aqu_of_resolved")
cmp <- rbindlist(lapply(TAU_GRID, function(t) rbindlist(lapply(metrics, function(m) {
  nullv <- sim[tau == t][[m]]; nullv <- nullv[is.finite(nullv)]; obs <- emp[tau == t][[m]]
  data.table(tau = t, metric = m, observed = obs,
             null_med = median(nullv), null_lo = quantile(nullv, .025), null_hi = quantile(nullv, .975),
             null_sd = sd(nullv),
             z = (obs - mean(nullv)) / sd(nullv),
             p_ge = (1 + sum(nullv >= obs)) / (length(nullv) + 1),
             p_le = (1 + sum(nullv <= obs)) / (length(nullv) + 1),
             n_rep = length(nullv))
}))))
saveRDS(list(sim = sim, emp = emp, cmp = cmp, n_rep = length(done)), OUTRDS)

cat(sprintf("\n=== per-SNP ancestry sorting: observed vs %d-rep simulated null ===\n", length(done)))
cat("--- pct_sorted (fraction of diagnostic SNPs sorted) ---\n")
print(cmp[metric == "pct_sorted", .(tau, observed = round(observed, 2),
      null_med = round(null_med, 2), null_95 = sprintf("[%.2f,%.2f]", null_lo, null_hi),
      z = round(z, 1), p_ge = signif(p_ge, 3))])
cat("--- pct_aqu_of_resolved (directional bias toward aquilonia) ---\n")
print(cmp[metric == "pct_aqu_of_resolved", .(tau, observed = round(observed, 1),
      null_med = round(null_med, 1), null_95 = sprintf("[%.1f,%.1f]", null_lo, null_hi),
      z = round(z, 1), p_ge = signif(p_ge, 3), p_le = signif(p_le, 3))])

## ---- figure ----------------------------------------------------------------
## Panel A: % sorted, observed vs the (near-zero) null jitter + 95% band.
## Panel B: raw number of sorted SNPs per replicate in the null vs the observed
## count (log1p), to make the degenerate null legible.
emp_ns <- as.data.table(readRDS(SWEEP_EMP))[level == "SNP", .(tau, n_sorted)]
pA <- ggplot(sim, aes(factor(tau), pct_sorted)) +
  geom_jitter(width = 0.18, height = 0, colour = "#1b9e77", alpha = 0.25, size = 0.8) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = NA, colour = "#0b3d2e", linewidth = 0.8) +
  geom_point(data = emp, aes(factor(tau), pct_sorted), colour = "#d95f02", size = 5, shape = 18) +  # bigger observed diamonds for slides
  labs(x = expression("sorting threshold  " * tau), y = "% diagnostic SNPs sorted") +
  theme_bw(base_size = 17) + theme(panel.grid.minor = element_blank())
pB <- ggplot(sim, aes(factor(tau), n_sorted + 1)) +
  geom_jitter(width = 0.18, height = 0, colour = "#1b9e77", alpha = 0.25, size = 0.8) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = NA, colour = "#0b3d2e", linewidth = 0.8) +
  geom_point(data = emp_ns, aes(factor(tau), n_sorted + 1), colour = "#d95f02", size = 5, shape = 18) +  # bigger observed diamonds for slides
  scale_y_log10() +
  labs(x = expression("sorting threshold  " * tau), y = "sorted SNP count + 1 (log scale)") +
  theme_bw(base_size = 17) + theme(panel.grid.minor = element_blank())
p <- pA + pB + patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 19))
ggsave(OUTPDF, p, width = 11, height = 5.2); ggsave(OUTPNG, p, width = 11, height = 5.2, dpi = 200)
cat("saved:", OUTPNG, "\n")
