## =============================================================================
## module_di25 -- ld_w_095 windowed-median Manhattan: empirical vs 1000-REP ENVELOPE
## =============================================================================
## Extends di25_ldw_manhattan_emp_vs_sim.R from ONE replicate to the full set of
## 1000 DIEM bootstrap replicates (data/diem_outs_demo/diem_boot<N>_output.bed).
## For every replicate we compute the SAME DI25 ld_w_095 statistic (hybrids-only
## LD-decay slide=100 corr -> compute_ld_w(0.95)) and its per-chromosome windowed
## median. Across the 1000 replicates each genomic window then gets a simulated
## distribution; the figure overlays the empirical line on the simulated median
## and 2.5-97.5% envelope.
##
## Resumable + parallel: each replicate's windowed medians are cached to
## module_di25/data/ldw_sim_cache/rep<N>.rds and skipped if present, so the job
## can be interrupted/re-run. Replicates run in parallel (WORKERS forks), each
## with per-replicate LD compute cores = 1.
##
## Run from the formica_hybrid repo root:
##   Rscript module_di25/R/di25_ldw_manhattan_envelope.R [nreps] [workers]
##     nreps   : number of replicates to use (default 1000)
##     workers : parallel forks (default 9)
## =============================================================================
suppressMessages({
  library(data.table); library(ggplot2); library(SNPRelate); library(parallel)
})
devtools::load_all("~/gitlab/LDscnR/")

args    <- commandArgs(trailingOnly = TRUE)
NREPS   <- if (length(args) >= 1) as.integer(args[1]) else 1000L
WORKERS <- if (length(args) >= 2) as.integer(args[2]) else 9L

STAGE1    <- "module_di25/data/di25_stage1.rds"
BEDFMT    <- "data/diem_outs_demo/diem_boot%d_output.bed"
CACHE_DIR <- "module_di25/data/ldw_sim_cache"
GDS_TMP   <- "/private/tmp/claude-539526166/-Users-petrikem-gitlab-LDscnR/63fce5ef-28e0-4215-afc5-5b7620c2c579/scratchpad"
OUTRDS    <- "module_di25/data/di25_ldw_envelope.rds"
OUTPNG    <- "module_di25/Figures/di25_ldw_manhattan_envelope.png"
OUTPDF    <- sub("\\.png$", ".pdf", OUTPNG)

RHO_LDW <- 0.95
WIN_BP  <- 500000L
MIN_N   <- 5L
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- windowed median helper (shared with the single-rep script) ------------
win_median <- function(dt) {
  d <- dt[is.finite(ld_w)]
  d[, win := (Pos %/% WIN_BP)]
  d[, .(mid = (win[1] * WIN_BP + WIN_BP / 2), n = .N, ld_w_med = median(ld_w)),
    by = .(chr_id, Chr, win)][n >= MIN_N]
}

## ---- per-replicate ld_w_095 -> windowed medians (cached) -------------------
process_rep <- function(REP) {
  cache <- file.path(CACHE_DIR, sprintf("rep%d.rds", REP))
  if (file.exists(cache)) return(invisible(TRUE))
  bed <- sprintf(BEDFMT, REP)
  if (!file.exists(bed)) return(invisible(FALSE))
  out <- tryCatch({
    hdr  <- readLines(bed, n = 2)
    inds <- strsplit(strsplit(hdr[2], "\t")[[1]][10], "\\|")[[1]]
    sim  <- fread(bed, skip = 2, header = FALSE, sep = "\t",
                  select = c(1, 3, 10), colClasses = list(character = c(1, 10)),
                  showProgress = FALSE)
    chr_id <- as.integer(sub("ch", "", sim$V1)); pos <- as.integer(sim$V3)
    geno   <- sub("^S", "", sim$V10)
    S <- matrix(unlist(strsplit(geno, "", fixed = TRUE), use.names = FALSE),
                nrow = length(geno), byrow = TRUE)
    dos <- matrix(NA_integer_, nrow(S), ncol(S))
    dos[S == "0"] <- 0L; dos[S == "1"] <- 1L; dos[S == "2"] <- 2L
    colnames(dos) <- inds
    hyb <- grep("^hyb_", inds)
    map <- data.table(Chr = paste0("Chr", chr_id), Pos = pos,
                      marker = paste0("Chr", chr_id, ":", pos), chr_id = chr_id)
    ord <- order(map$chr_id, map$Pos); map <- map[ord]; dos <- dos[ord, , drop = FALSE]
    GTs_hyb <- t(dos[, hyb, drop = FALSE])
    gpath <- file.path(GDS_TMP, sprintf("di25_sim_rep%d.gds", REP))
    gds <- create_gds_from_geno(geno = GTs_hyb, map = map, gds_path = gpath)
    ld_decay <- compute_LD_decay(gds, keep_el = TRUE, slide = 100,
                                 ld_method = "corr", cores = 1)
    snpgdsClose(gds); unlink(gpath)
    lw  <- compute_ld_w(ld_decay, rho = RHO_LDW, cores = 1)
    ids <- unlist(lapply(ld_decay$by_chr, function(o) o$snp_ids), use.names = FALSE)
    map[, ld_w := lw[match(marker, ids)]]
    wm <- win_median(map[, .(Chr, Pos, ld_w, chr_id)])
    wm[, rep := REP]
    saveRDS(wm, cache)
    TRUE
  }, error = function(e) { message(sprintf("  rep%d FAILED: %s", REP, conditionMessage(e))); FALSE })
  invisible(out)
}

## ---- run (parallel, resumable) ---------------------------------------------
reps  <- seq_len(NREPS)
todo  <- reps[!file.exists(file.path(CACHE_DIR, sprintf("rep%d.rds", reps)))]
message(sprintf("[envelope] %d reps requested, %d already cached, %d to compute (workers=%d)",
                NREPS, NREPS - length(todo), length(todo), WORKERS))
if (length(todo)) {
  invisible(mclapply(todo, process_rep, mc.cores = WORKERS, mc.preschedule = FALSE))
}

## ---- aggregate simulated windows across replicates -------------------------
done  <- reps[file.exists(file.path(CACHE_DIR, sprintf("rep%d.rds", reps)))]
message(sprintf("[envelope] aggregating %d cached replicates", length(done)))
sim   <- rbindlist(lapply(done, function(r) readRDS(file.path(CACHE_DIR, sprintf("rep%d.rds", r)))))
env   <- sim[, .(sim_med = median(ld_w_med),
                 sim_lo  = quantile(ld_w_med, 0.025),
                 sim_hi  = quantile(ld_w_med, 0.975),
                 n_rep   = .N),
             by = .(chr_id, Chr, win, mid)]

## ---- empirical windowed median (same statistic, cached ld_w_095) -----------
emp <- as.data.table(readRDS(STAGE1)$map)[, .(Chr, Pos, ld_w = ld_w_095)]
emp[, chr_id := as.integer(sub("Chr", "", Chr))]
emp_wm <- win_median(emp)[, .(chr_id, Chr, win, mid, emp = ld_w_med)]

pl <- merge(env, emp_wm, by = c("chr_id", "Chr", "win", "mid"), all = TRUE)
setorder(pl, chr_id, win)
pl[, Chr := factor(Chr, levels = paste0("Chr", sort(unique(pl$chr_id))))]
saveRDS(list(pl = pl, sim = sim, n_rep = length(done),
             win_bp = WIN_BP, rho = RHO_LDW), OUTRDS)

## ---- plot: empirical line over simulated median + 95% envelope -------------
p <- ggplot(pl, aes(mid / 1e6)) +
  geom_ribbon(aes(ymin = sim_lo, ymax = sim_hi, fill = "simulated 95% envelope"), alpha = 0.4) +
  geom_line(aes(y = sim_med, colour = "simulated median"), linewidth = 0.6) +
  geom_line(aes(y = emp, colour = "empirical"), linewidth = 0.8) +
  facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = c("empirical" = "#d95f02", "simulated median" = "#1b9e77"), name = NULL) +
  scale_fill_manual(values = c("simulated 95% envelope" = "#66c2a5"), name = NULL) +
  labs(x = "position (Mbp)", y = expression("windowed median local LD  " * ld[w][0.95])) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(), panel.spacing = unit(1, "pt"),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10, angle = 90, face = "bold"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "top", legend.text = element_text(size = 14),
        legend.key.width = unit(1.4, "cm"))
ggsave(OUTPDF, p, width = 15, height = 5)
ggsave(OUTPNG, p, width = 15, height = 5, dpi = 150)

## ---- summary: coverage of the empirical line by the simulated envelope ------
cov <- pl[is.finite(emp) & is.finite(sim_lo)]
cat(sprintf("\n=== ld_w_0.95 envelope (%d replicates) ===\n", length(done)))
cat(sprintf("windows compared: %d\n", nrow(cov)))
cat(sprintf("empirical within simulated 95%% envelope: %.1f%% of windows\n",
            100 * cov[, mean(emp >= sim_lo & emp <= sim_hi)]))
cat(sprintf("empirical ABOVE envelope: %.1f%%  |  BELOW: %.1f%%\n",
            100 * cov[, mean(emp > sim_hi)], 100 * cov[, mean(emp < sim_lo)]))
cat(sprintf("mean windowed median  emp=%.3f  sim=%.3f\n",
            cov[, mean(emp)], cov[, mean(sim_med)]))
cat("saved: ", OUTRDS, "\n       ", OUTPNG, "\n       ", OUTPDF, "\n")
