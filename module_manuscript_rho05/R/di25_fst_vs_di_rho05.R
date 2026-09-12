## =============================================================================
## module_manuscript_rho05 -- Fig [FST_simulations], rho05 version
## =============================================================================
## AUDIT FIX (2026-09-12, Issue 2): the previous version used only
## `b$stats$best_marker` (17,509 rows -- the large clusters that have an
## eMLG/best-SNP representation), omitting the 643,877 singleton and
## small-cluster representatives. The COMPLETE LD-reduced representation is
## `b$rep_snp_all` (661,386 rows: `rep_snp` = best_marker where has_eMLG==TRUE,
## else the centrality `representative` for the small/singleton clusters) --
## verified: 661,386 unique group_id, 661,386 unique rep_snp, no NAs, no dups.
## This version uses one marker for EVERY LD-reduced unit, labels the series
## "all LD-reduced units" (not "all SNPs"), and reports the large-cluster
## best-SNP vs small/singleton representative-SNP counts separately.
##
## DECISION (pipeline author, 2026-09-12): apply MIN_PARENT_MAF=0.15 as the
## PRIMARY gate on the empirical Fst-by-DI-bin curve (matching Module A's
## locked convention pipeline-wide); the ungated (all parental MAF) curve is
## now reported explicitly as a sensitivity check, alongside the existing
## MAF-stratified breakdown. The sim/neutral comparison band is NOT
## MAF-gated (it answers a different question -- whether empirical
## differentiation at each DI bin exceeds neutral expectations -- and every
## sim-overlap marker already carries its own real pmaf value if a MAF view
## of the sim comparison is wanted later).
##
## Adapted from module_di25/R/di25_fst_vs_di.R / the previous rho05 version.
## Because the marker universe is now completely different (661,386 vs
## 17,509 units -> a different best/representative-SNP set per cluster), the
## per-replicate simulation cache CANNOT be reused -- this reruns the
## 1000-replicate Fst computation against the corrected unit set, in a FRESH
## cache directory (fst_sim_cache_full_rho05_v2/) with an explicit fingerprint
## per cache entry (replicate number, marker-list hash, DI-bin hash, parameter
## hash, source sim-file identity, script version) -- aggregation ABORTS if
## fewer than the requested replicates complete or any fingerprint mismatches.
##
## Everything else (fixed DI bins, background LD, Weir & Cockerham estimator)
## is unchanged from the previous version.
##
## Run from the repo root:
##   Rscript module_manuscript_rho05/R/di25_fst_vs_di_rho05.R [nreps] [workers]   (default 1000, 9)
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(parallel); library(digest) })

args    <- commandArgs(trailingOnly = TRUE)
NREPS   <- if (length(args) >= 1) as.integer(args[1]) else 1000L
WORKERS <- if (length(args) >= 2) as.integer(args[2]) else 9L
BEDFMT    <- "data/diem_outs_demo/diem_boot%d_output.bed"
CACHE_DIR <- "module_manuscript_rho05/data/fst_sim_cache_full_rho05_v2"   # FRESH -- do not reuse the old cache
OUTRDS    <- "module_manuscript_rho05/data/di25_fst_vs_di_rho05.rds"
FIGDIR    <- "module_manuscript_rho05/Figures"
OUTPNG    <- file.path(FIGDIR, "di25_fst_vs_di_rho05.png")
OUTPDF    <- sub("\\.png$", ".pdf", OUTPNG)
DI_BREAKS <- c(-Inf, -90, -75, -60, -50, -40, -30, -25, -20, -15, Inf)
BIN_LAB   <- levels(cut(0, DI_BREAKS)); N_BIN <- length(BIN_LAB)
BG_NSUB   <- 2000L; BG_Q <- 0.95; MIN_SIM_UNITS <- 30L
MIN_PARENT_MAF <- 0.15   # AUDIT FIX Issue 2 decision: primary gate, matching Module A's locked convention
SCRIPT_VERSION <- "di25_fst_vs_di_rho05_AUDITFIX_2026-09-12"
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

## ---- Weir & Cockerham 1984 per-locus a and (a+b+c) --------------------------
wc_ac <- function(G, pop) {
  levs <- unique(pop); M <- ncol(G)
  N <- P <- H <- matrix(0, length(levs), M)
  for (k in seq_along(levs)) {
    g <- G[pop == levs[k], , drop = FALSE]
    n <- colSums(!is.na(g)); s <- colSums(g, na.rm = TRUE); het <- colSums(g == 1, na.rm = TRUE)
    N[k, ] <- n; P[k, ] <- ifelse(n > 0, s / (2 * n), 0); H[k, ] <- ifelse(n > 0, het / n, 0)
  }
  C <- colSums(N); sumN2 <- colSums(N^2); r <- colSums(N > 0)
  nbar <- C / r; nc <- (C - sumN2 / C) / (r - 1)
  pbar <- colSums(N * P) / C; hbar <- colSums(N * H) / C
  s2  <- colSums(N * sweep(P, 2, pbar)^2) / ((r - 1) * nbar)
  msp <- pbar * (1 - pbar)
  a  <- (nbar / nc) * (s2 - (1 / (nbar - 1)) * (msp - ((r - 1) / r) * s2 - 0.25 * hbar))
  b  <- (nbar / (nbar - 1)) * (msp - ((r - 1) / r) * s2 - ((2 * nbar - 1) / (4 * nbar)) * hbar)
  cc <- 0.5 * hbar
  list(a = a, abc = a + b + cc)
}
fst_ratio  <- function(ac, ok = TRUE) { ok <- ok & is.finite(ac$a) & is.finite(ac$abc)
  s <- sum(ac$abc[ok]); if (s <= 0) NA_real_ else sum(ac$a[ok]) / s }
fst_by_bin <- function(ac, bin) vapply(seq_len(N_BIN), function(d) fst_ratio(ac, bin == d), numeric(1))
bg_ld <- function(G, chr) {
  chrs <- unique(chr); Mtot <- ncol(G)
  pool <- unlist(lapply(chrs, function(ch) { ix <- which(chr == ch)
    sample(ix, min(length(ix), ceiling(BG_NSUB * length(ix) / Mtot))) }))
  R <- suppressWarnings(cor(G[, pool], use = "pairwise.complete.obs"))^2
  inter <- outer(chr[pool], chr[pool], "!=") & upper.tri(R)
  as.numeric(quantile(R[inter], BG_Q, na.rm = TRUE))
}

## ---- FULL-data LD-reduced units: ONE MARKER PER UNIT (all 661,386) ----------
## AUDIT FIX (Issue 2): b$rep_snp_all, not b$stats$best_marker. rep_snp is the
## best-SNP for the 17,509 large (has_eMLG) clusters, the centrality
## representative for the 643,877 small/singleton clusters -- verified
## unique group_id, unique rep_snp, no NA, no dup (checked interactively
## before writing this fix).
b   <- readRDS("module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05_bestsnp.rds")
e   <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e)
G0  <- e$GTs_hybrids_005
mDI <- as.data.table(e$map_hyb_005)
sd  <- as.data.table(e$sample_data)

r <- b$rep_snp_all
stopifnot("AUDIT FIX (Issue 2): expected exactly 661,386 LD-reduced units" = nrow(r) == 661386,
          "duplicate group_id in rep_snp_all"  = !anyDuplicated(r$group_id),
          "duplicate rep_snp marker in rep_snp_all" = !anyDuplicated(r$rep_snp),
          "NA rep_snp in rep_snp_all" = !anyNA(r$rep_snp))
units <- data.table(group_id = r$group_id, best = r$rep_snp,
                    rep_type = ifelse(r$has_eMLG, "best_snp_large_cluster", "representative_small_or_singleton"))
units[, DI := mDI$DiagnosticIndex[match(best, mDI$marker)]]
units <- units[is.finite(DI)]
units[, bin := cut(DI, DI_BREAKS, labels = FALSE)]
stopifnot("units lost uniqueness after DI join" = !anyDuplicated(units$group_id))
n_best <- units[rep_type == "best_snp_large_cluster", .N]
n_repr <- units[rep_type == "representative_small_or_singleton", .N]
message(sprintf("[fst-rho05] %d full-data LD-reduced units (rho05): %d best-SNP (large cluster, has eMLG) + %d representative-SNP (small/singleton); %d fixed DI bins (%s .. %s)",
                nrow(units), n_best, n_repr, N_BIN, BIN_LAB[1], BIN_LAB[N_BIN]))

## ---- empirical Fst per fixed DI bin + background LD (ALL units, ungated) ----
pop_emp <- sd$Population[match(rownames(G0), sd$Sample_ID)]
Gb      <- G0[, units$best]
ac_emp  <- wc_ac(Gb, pop_emp)
emp_fst_ungated <- fst_by_bin(ac_emp, units$bin)
message(sprintf("[fst-rho05] empirical overall Fst (ALL LD-reduced units, ungated) = %.3f (neutral bg DI<-90 %.3f -> diagnostic %.3f)",
                fst_ratio(ac_emp), emp_fst_ungated[1], emp_fst_ungated[N_BIN]))
emp_chr <- as.integer(sub("Chr", "", sub(":.*", "", colnames(G0))))
set.seed(1); emp_bg <- bg_ld(G0, emp_chr)
message(sprintf("[fst-rho05] empirical background LD = %.4f", emp_bg))

## ---- parental MAF (joined explicitly by marker name) ------------------------
ep <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = ep)
par_rows <- grepl("_parent$", ep$sample_data_with_parents$Population)
pf <- colMeans(ep$GTs_with_parents[par_rows, units$best, drop = FALSE], na.rm = TRUE) / 2
units[, pmaf := pf[match(best, names(pf))]]   # explicit marker-name join, not positional
stopifnot("pmaf join lost rows" = length(pf) == nrow(units))
if (anyNA(units$pmaf)) message(sprintf("[fst-rho05] note: %d units have NA parental MAF (no non-missing parent genotype calls at that marker) -- excluded from the MAF>=%.2f gate and from MAF-stratified bins",
                                       sum(is.na(units$pmaf)), MIN_PARENT_MAF))
rm(ep); gc()

## ---- AUDIT FIX (Issue 2, decision): primary Fst curve gated on parental MAF ----
units_primary <- units[pmaf >= MIN_PARENT_MAF]
n_best_p <- units_primary[rep_type == "best_snp_large_cluster", .N]
n_repr_p <- units_primary[rep_type == "representative_small_or_singleton", .N]
message(sprintf("[fst-rho05] PRIMARY gate parental MAF>=%.2f: %d of %d units retained (%d best-SNP + %d representative-SNP)",
                MIN_PARENT_MAF, nrow(units_primary), nrow(units), n_best_p, n_repr_p))
Gb_p     <- G0[, units_primary$best]
ac_emp_p <- wc_ac(Gb_p, pop_emp)
emp_fst_primary <- fst_by_bin(ac_emp_p, units_primary$bin)
message(sprintf("[fst-rho05] empirical overall Fst (PRIMARY, parental MAF>=%.2f) = %.3f",
                MIN_PARENT_MAF, fst_ratio(ac_emp_p)))

MAF_BREAKS <- c(0, 0.1, 0.2, 0.35, 0.5); MAF_LAB <- levels(cut(0, MAF_BREAKS, include.lowest = TRUE))
units[, mstr := cut(pmaf, MAF_BREAKS, include.lowest = TRUE)]
MIN_STRAT <- 30L
strat <- rbindlist(lapply(MAF_LAB, function(ms) data.table(mstr = ms, bin = seq_len(N_BIN),
  n   = vapply(seq_len(N_BIN), function(d) sum(!is.na(units$mstr) & units$mstr == ms & units$bin == d), integer(1)),
  fst = vapply(seq_len(N_BIN), function(d) fst_ratio(ac_emp, !is.na(units$mstr) & units$mstr == ms & units$bin == d), numeric(1)))))
strat <- strat[n >= MIN_STRAT]; strat[, mstr := factor(mstr, levels = MAF_LAB)]
message(sprintf("[fst-rho05] parental-MAF strata (sensitivity, all units): %s", paste(MAF_LAB, collapse = ", ")))

rm(G0, Gb, Gb_p, ac_emp_p, e); gc()

## ---- sim: LD-reduced units (full, ungated universe) that fall in the DI25 sim panel ----
## Ungated deliberately -- the neutral-sim comparison asks whether empirical
## differentiation at each DI bin exceeds neutral expectations, a different
## question from the parental-MAF gate on the primary empirical curve; every
## sim-overlap marker already carries its own real pmaf (joined above) if a
## MAF-conditioned view of this comparison is wanted later.
sim_mk  <- { hdr <- readLines(sprintf(BEDFMT, 1), n = 2); s1 <- fread(sprintf(BEDFMT, 1), skip = 2,
              header = FALSE, sep = "\t", select = c(1, 3), colClasses = list(character = 1), showProgress = FALSE)
            paste0("Chr", sub("ch", "", s1$V1), ":", s1$V3) }
ov <- units[best %in% sim_mk]
message(sprintf("[fst-rho05] %d LD-reduced units overlap the DI25 sim panel (bins %s) -- AUDIT FIX expects ~16,616 (was 1,511 pre-fix)",
                nrow(ov), paste(range(ov$bin), collapse = "-")))

## ---- fingerprinting for the fresh, per-replicate cache ----------------------
marker_hash  <- digest(ov$best, algo = "md5")     # ov's row order is fixed for this run
dibin_hash   <- digest(ov$bin,  algo = "md5")
params_hash  <- digest(list(DI_BREAKS = DI_BREAKS, MAF_BREAKS = MAF_BREAKS,
                            MIN_PARENT_MAF = MIN_PARENT_MAF, BG_NSUB = BG_NSUB, BG_Q = BG_Q,
                            BEDFMT = BEDFMT), algo = "md5")
run_fingerprint <- list(marker_hash = marker_hash, dibin_hash = dibin_hash,
                        params_hash = params_hash, script_version = SCRIPT_VERSION)
message(sprintf("[fst-rho05] cache fingerprint: marker=%s dibin=%s params=%s",
                substr(marker_hash, 1, 8), substr(dibin_hash, 1, 8), substr(params_hash, 1, 8)))

process_rep <- function(REP) {
  cache <- file.path(CACHE_DIR, sprintf("rep%d.rds", REP))
  if (file.exists(cache)) {
    ck <- tryCatch(readRDS(cache), error = function(e) NULL)
    if (!is.null(ck) && !is.null(ck$fingerprint) &&
        identical(ck$fingerprint[c("marker_hash", "dibin_hash", "params_hash", "script_version")],
                  run_fingerprint))
      return(invisible(TRUE))
    message(sprintf("  rep%d: cache fingerprint mismatch or malformed -- recomputing", REP))
  }
  bed <- sprintf(BEDFMT, REP); if (!file.exists(bed)) return(invisible(FALSE))
  bed_md5 <- unname(tools::md5sum(bed))
  out <- tryCatch({
    hdr  <- readLines(bed, n = 2)
    inds <- strsplit(strsplit(hdr[2], "\t")[[1]][10], "|", fixed = TRUE)[[1]]
    sim  <- fread(bed, skip = 2, header = FALSE, sep = "\t", select = c(1, 3, 10),
                  colClasses = list(character = c(1, 10)), showProgress = FALSE)
    markers <- paste0("Chr", sub("ch", "", sim$V1), ":", sim$V3)
    chr_id  <- as.integer(sub("ch", "", sim$V1))
    S <- matrix(unlist(strsplit(sub("^S", "", sim$V10), "", fixed = TRUE), use.names = FALSE), nrow = nrow(sim), byrow = TRUE)
    dos <- matrix(NA_integer_, nrow(S), ncol(S)); dos[S == "0"] <- 0L; dos[S == "1"] <- 1L; dos[S == "2"] <- 2L
    hyb <- grep("^hyb_", inds); pop <- sub("^hyb_(.*)_[0-9]+$", "\\1", inds[hyb])
    G <- t(dos[, hyb, drop = FALSE]); colnames(G) <- markers
    mi  <- match(ov$best, colnames(G)); keep <- !is.na(mi)
    ac  <- wc_ac(G[, mi[keep], drop = FALSE], pop)
    fst <- fst_by_bin(ac, ov$bin[keep]); overall <- fst_ratio(ac)
    set.seed(100 + REP); bg <- bg_ld(G, chr_id)
    saveRDS(list(fst = fst, overall = overall, bg = bg, rep = REP,
                 fingerprint = c(run_fingerprint, list(sim_file = normalizePath(bed), sim_file_md5 = bed_md5))),
            paste0(cache, ".tmp"))
    file.rename(paste0(cache, ".tmp"), cache)   # atomic write
    TRUE
  }, error = function(e) { message(sprintf("  rep%d FAILED: %s", REP, conditionMessage(e))); FALSE })
  invisible(out)
}
reps <- seq_len(NREPS)
already_ok <- vapply(reps, function(REP) {
  cache <- file.path(CACHE_DIR, sprintf("rep%d.rds", REP))
  if (!file.exists(cache)) return(FALSE)
  ck <- tryCatch(readRDS(cache), error = function(e) NULL)
  !is.null(ck) && !is.null(ck$fingerprint) &&
    identical(ck$fingerprint[c("marker_hash", "dibin_hash", "params_hash", "script_version")], run_fingerprint)
}, logical(1))
todo <- reps[!already_ok]
message(sprintf("[fst-rho05] %d reps requested, %d already cached (fingerprint OK), %d to compute (workers=%d)",
                NREPS, NREPS - length(todo), length(todo), WORKERS))
if (length(todo)) invisible(mclapply(todo, process_rep, mc.cores = WORKERS, mc.preschedule = FALSE))

## ---- aggregate: ABORT if incomplete or any fingerprint mismatch -------------
done <- reps[file.exists(file.path(CACHE_DIR, sprintf("rep%d.rds", reps)))]
sm   <- lapply(done, function(r) readRDS(file.path(CACHE_DIR, sprintf("rep%d.rds", r))))
fp_ok <- vapply(sm, function(x) !is.null(x$fingerprint) &&
                  identical(x$fingerprint[c("marker_hash", "dibin_hash", "params_hash", "script_version")], run_fingerprint),
                logical(1))
if (length(done) < NREPS)
  stop(sprintf("ABORT: only %d/%d replicates completed (aggregation requires exactly %d). Rerun this script to fill in the remainder -- already-cached, fingerprint-matching replicates are reused automatically.",
               length(done), NREPS, NREPS))
if (!all(fp_ok))
  stop(sprintf("ABORT: %d cached replicate(s) have a fingerprint mismatching the current run (marker set / DI bins / params / script version changed since they were cached). Delete %s and rerun, or investigate why a stale replicate survived.",
               sum(!fp_ok), CACHE_DIR))
message(sprintf("[fst-rho05] all %d replicates present and fingerprint-verified -- aggregating", NREPS))

simF <- do.call(rbind, lapply(sm, function(x) x$fst))
simO <- vapply(sm, function(x) x$overall, numeric(1))
simB <- vapply(sm, function(x) x$bg, numeric(1))
n_emp <- units[, .N, by = bin][order(bin)]
n_emp_p <- units_primary[, .N, by = bin][order(bin)]
n_ov  <- ov[, .N, by = bin][order(bin)]
env  <- data.table(bin = seq_len(N_BIN), DI_bin = BIN_LAB,
                   emp = emp_fst_primary, emp_ungated = emp_fst_ungated,
                   n_emp_units = n_emp_p$N[match(seq_len(N_BIN), n_emp_p$bin)],
                   n_emp_units_ungated = n_emp$N[match(seq_len(N_BIN), n_emp$bin)],
                   n_sim_units = n_ov$N[match(seq_len(N_BIN), n_ov$bin)],
                   sim_med = apply(simF, 2, median, na.rm = TRUE),
                   sim_lo  = apply(simF, 2, quantile, 0.025, na.rm = TRUE),
                   sim_hi  = apply(simF, 2, quantile, 0.975, na.rm = TRUE))
env[is.na(n_emp_units) | n_emp_units < MIN_SIM_UNITS, `:=`(emp = NA)]
env[is.na(n_sim_units) | n_sim_units < MIN_SIM_UNITS, `:=`(sim_med = NA, sim_lo = NA, sim_hi = NA)]
neutral <- c(med = median(simO), lo = quantile(simO, 0.025), hi = quantile(simO, 0.975))
saveRDS(list(env = env, strat = strat, neutral = neutral, emp_bg = emp_bg, sim_bg = simB,
             n_rep = length(done), di_breaks = DI_BREAKS, maf_breaks = MAF_BREAKS,
             min_parent_maf_primary = MIN_PARENT_MAF,
             n_units_total = nrow(units), n_units_best_snp = n_best, n_units_representative = n_repr,
             n_units_primary = nrow(units_primary), n_units_primary_best_snp = n_best_p,
             n_units_primary_representative = n_repr_p, n_sim_overlap = nrow(ov),
             run_fingerprint = run_fingerprint), OUTRDS)

cat(sprintf("\n=== [rho05, AUDIT FIX] among-population Fst by DI bin (empirical, PRIMARY parental MAF>=%.2f, vs %d-rep high-DI neutral sim) ===\n",
            MIN_PARENT_MAF, length(done)))
print(env[, .(DI_bin, n_emp_units, emp = round(emp, 3), emp_ungated = round(emp_ungated, 3), n_sim_units, sim_med = round(sim_med, 3))])
cat(sprintf("neutral sim Fst (over its high-DI units): %.3f [%.3f, %.3f]\n", neutral[1], neutral[2], neutral[3]))
cat(sprintf("background LD (inter-chr r^2, q%.2f): empirical %.4f | sim %.4f [%.4f, %.4f]\n",
            BG_Q, emp_bg, median(simB), quantile(simB, 0.025), quantile(simB, 0.975)))
cat(sprintf("units: %d total (%d best-SNP + %d representative) -> %d retained at MAF>=%.2f (%d best-SNP + %d representative) -> %d overlap the sim panel\n",
            nrow(units), n_best, n_repr, nrow(units_primary), MIN_PARENT_MAF, n_best_p, n_repr_p, nrow(ov)))

## ---- figure -----------------------------------------------------------------
p <- ggplot(env, aes(bin)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = neutral[2], ymax = neutral[3], fill = "#66c2a5", alpha = 0.30) +
  geom_hline(yintercept = neutral[1], colour = "#1b9e77", linetype = 2, linewidth = 0.7) +
  geom_line(data = strat, aes(bin, fst, linetype = mstr), colour = "#d95f02", linewidth = 0.55, na.rm = TRUE) +
  geom_line(aes(y = emp_ungated, colour = "empirical, ungated (sensitivity)"), linewidth = 0.6, linetype = "dashed", na.rm = TRUE) +
  geom_line(aes(y = emp, colour = "empirical, LD-reduced units (parental MAF>=0.15)"), linewidth = 1.0, na.rm = TRUE) +
  geom_point(aes(y = emp, colour = "empirical, LD-reduced units (parental MAF>=0.15)"), size = 2.4, na.rm = TRUE) +
  geom_point(aes(y = sim_med, colour = "high-DI neutral sim"), size = 2.2, na.rm = TRUE) +
  geom_errorbar(aes(ymin = sim_lo, ymax = sim_hi, colour = "high-DI neutral sim"), width = 0.15, na.rm = TRUE) +
  annotate("text", x = 1, y = neutral[1], label = "neutral sim (high-DI)", colour = "#1b7f63",
           hjust = 0, vjust = -0.6, size = 3.4) +
  scale_x_continuous(breaks = seq_len(N_BIN), labels = BIN_LAB) +
  scale_colour_manual(values = c("empirical, LD-reduced units (parental MAF>=0.15)" = "#d95f02",
                                "empirical, ungated (sensitivity)" = "#7570b3",
                                "high-DI neutral sim" = "#1b9e77"), name = NULL) +
  scale_linetype_manual(values = c("[0,0.1]" = "dotted", "(0.1,0.2]" = "dotdash",
                                   "(0.2,0.35]" = "longdash", "(0.35,0.5]" = "22"), name = "parental MAF stratum (all units, sensitivity)") +
  labs(x = "DiagnosticIndex bin  (left = near-neutral background, right = ancestry-informative)",
       y = expression("among-population " * F[ST] * "  (Weir & Cockerham), rho05 LD-reduced units"),
       caption = sprintf("All %s LD-reduced units (best-SNP for large clusters, representative-SNP for small/singleton); primary curve gated at parental MAF>=%.2f.",
                          format(nrow(units), big.mark = ","), MIN_PARENT_MAF)) +
  guides(colour = guide_legend(order = 1, nrow = 2), linetype = guide_legend(order = 2, nrow = 1)) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), legend.position = "top", legend.box = "vertical",
        legend.margin = margin(1, 1, 1, 1), legend.spacing.y = unit(1, "pt"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(8, 12, 4, 6), plot.caption = element_text(size = 8, colour = "grey40"))
ggsave(OUTPDF, p, width = 8.6, height = 6.0); ggsave(OUTPNG, p, width = 8.6, height = 6.0, dpi = 200)
cat("saved:", OUTPNG, "\n")
