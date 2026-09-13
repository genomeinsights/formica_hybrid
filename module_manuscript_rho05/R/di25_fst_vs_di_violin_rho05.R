## =============================================================================
## module_manuscript_rho05 -- per-locus Fst distribution by DI bin (violin +
## boxplot), complementing the pooled Fst-vs-DI curve (di25_fst_vs_di_rho05.R)
## =============================================================================
## The existing Fig [FST_simulations] plots ONE pooled Fst per DI bin (Weir &
## Cockerham's a/(a+b+c) summed across all loci in the bin) -- this hides the
## per-locus VARIATION within a bin entirely. Not every locus in a
## high-DI bin is expected to be "sorted" (directionally fixed between
## aquilonia/polyctena), so per-locus Fst should vary a lot within a bin, and
## conceivably show two modes: a neutral-like cluster (near the simulated
## background) and a differentiated cluster (near/above the empirical
## pooled value) -- this script visualises that directly.
##
## Per-locus Fst is NOT a new computation -- wc_ac()'s `a` and `abc` are
## already per-locus vectors (di25_fst_vs_di_rho05.R pools them across a bin
## via sum(a)/sum(abc) to get ONE value; here they're just left unpooled:
## fst_locus = a / abc, one value per LD-reduced unit). Reuses the identical
## unit construction (all 661,386 LD-reduced units, rep_snp_all; parental
## MAF>=0.15 primary gate; same fixed DI bins) as the audited
## di25_fst_vs_di_rho05.R, and the already-computed neutral-simulation
## reference band from its saved output (no simulation rerun needed).
##
## Run from the repo root:
##   Rscript module_manuscript_rho05/R/di25_fst_vs_di_violin_rho05.R
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(parallel); library(digest) })

FST_RDS   <- "module_manuscript_rho05/data/di25_fst_vs_di_rho05.rds"
BEDFMT    <- "data/diem_outs_demo/diem_boot%d_output.bed"
SIM_LOCUS_CACHE <- "module_manuscript_rho05/data/fst_sim_locus_cache_rho05"
OUTRDS    <- "module_manuscript_rho05/data/di25_fst_vs_di_violin_rho05.rds"
FIGDIR    <- "module_manuscript_rho05/Figures"
OUTPNG    <- file.path(FIGDIR, "di25_fst_vs_di_violin_rho05.png")
OUTPDF    <- sub("\\.png$", ".pdf", OUTPNG)
DI_BREAKS <- c(-Inf, -90, -75, -60, -50, -40, -30, -25, -20, -15, Inf)
BIN_LAB   <- levels(cut(0, DI_BREAKS)); N_BIN <- length(BIN_LAB)
MIN_PARENT_MAF <- 0.15
WORKERS   <- 9L
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SIM_LOCUS_CACHE, showWarnings = FALSE, recursive = TRUE)

fst_prev <- readRDS(FST_RDS)
stopifnot("MIN_PARENT_MAF must match the audited di25_fst_vs_di_rho05.R" =
            fst_prev$min_parent_maf_primary == MIN_PARENT_MAF,
          identical(fst_prev$di_breaks, DI_BREAKS))
neutral <- fst_prev$neutral   # med/lo/hi of the POOLED per-replicate Fst (kept for reference only)

## ---- USER FIX (2026-09-13): the simulated reference must be a PER-LOCUS
## distribution to be comparable with the empirical per-locus violins -- the
## previous version used `overall` (one multilocus-POOLED Fst per replicate),
## whose narrowness came from averaging across ~16,616 loci per replicate,
## not from neutral loci genuinely having low locus-to-locus variance. Fixed
## by re-parsing all 1000 simulation replicates and keeping wc_ac()'s
## per-locus a/abc UNPOOLED (one value per locus per replicate, exactly like
## the empirical fst_locus below), cached per replicate (fingerprinted,
## atomic writes) since this reruns the same genotype parsing as the pooled
## simulation (~8 min for all 1000 reps at 9 workers, observed empirically --
## cheap enough to redo properly rather than reuse the pooled cache).

## ---- Weir & Cockerham 1984 per-locus a and (a+b+c) -- IDENTICAL to di25_fst_vs_di_rho05.R --
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

## ---- FULL-data LD-reduced units (all 661,386; best-SNP or representative) --
b   <- readRDS("module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05_bestsnp.rds")
e   <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e)
G0  <- e$GTs_hybrids_005
mDI <- as.data.table(e$map_hyb_005)
sd  <- as.data.table(e$sample_data)

r <- b$rep_snp_all
stopifnot(nrow(r) == 661386, !anyDuplicated(r$group_id), !anyDuplicated(r$rep_snp), !anyNA(r$rep_snp))
units <- data.table(group_id = r$group_id, best = r$rep_snp,
                    rep_type = ifelse(r$has_eMLG, "best_snp_large_cluster", "representative_small_or_singleton"))
units[, DI := mDI$DiagnosticIndex[match(best, mDI$marker)]]
units <- units[is.finite(DI)]
units[, bin := cut(DI, DI_BREAKS, labels = FALSE)]

## ---- per-locus simulated (neutral) Fst: re-parse all sim replicates, keep
## wc_ac()'s per-locus a/abc UNPOOLED (see USER FIX note above) ---------------
sim_mk <- { hdr <- readLines(sprintf(BEDFMT, 1), n = 2); s1 <- fread(sprintf(BEDFMT, 1), skip = 2,
             header = FALSE, sep = "\t", select = c(1, 3), colClasses = list(character = 1), showProgress = FALSE)
           paste0("Chr", sub("ch", "", s1$V1), ":", s1$V3) }
ov <- units[best %in% sim_mk]
message(sprintf("[fst-violin] %d units overlap the DI25 sim panel (per-locus sim extraction)", nrow(ov)))

marker_hash <- digest(ov$best, algo = "md5")
params_hash <- digest(list(BEDFMT = BEDFMT), algo = "md5")
run_fp <- list(marker_hash = marker_hash, params_hash = params_hash,
               script_version = "di25_fst_vs_di_violin_rho05_locus_v1")

process_rep_locus <- function(REP) {
  cache <- file.path(SIM_LOCUS_CACHE, sprintf("rep%d.rds", REP))
  if (file.exists(cache)) {
    ck <- tryCatch(readRDS(cache), error = function(e) NULL)
    if (!is.null(ck) && !is.null(ck$fingerprint) &&
        identical(ck$fingerprint[c("marker_hash", "params_hash", "script_version")], run_fp))
      return(invisible(TRUE))
  }
  bed <- sprintf(BEDFMT, REP); if (!file.exists(bed)) return(invisible(FALSE))
  out <- tryCatch({
    hdr  <- readLines(bed, n = 2)
    inds <- strsplit(strsplit(hdr[2], "\t")[[1]][10], "|", fixed = TRUE)[[1]]
    sim  <- fread(bed, skip = 2, header = FALSE, sep = "\t", select = c(1, 3, 10),
                  colClasses = list(character = c(1, 10)), showProgress = FALSE)
    markers <- paste0("Chr", sub("ch", "", sim$V1), ":", sim$V3)
    S <- matrix(unlist(strsplit(sub("^S", "", sim$V10), "", fixed = TRUE), use.names = FALSE), nrow = nrow(sim), byrow = TRUE)
    dos <- matrix(NA_integer_, nrow(S), ncol(S)); dos[S == "0"] <- 0L; dos[S == "1"] <- 1L; dos[S == "2"] <- 2L
    hyb <- grep("^hyb_", inds); pop <- sub("^hyb_(.*)_[0-9]+$", "\\1", inds[hyb])
    G <- t(dos[, hyb, drop = FALSE]); colnames(G) <- markers
    mi  <- match(ov$best, colnames(G)); keep <- !is.na(mi)
    ac  <- wc_ac(G[, mi[keep], drop = FALSE], pop)
    fst_locus_sim <- ac$a / ac$abc   # UNPOOLED -- one value per overlap locus, this replicate
    saveRDS(list(fst_locus = fst_locus_sim, rep = REP, fingerprint = run_fp), paste0(cache, ".tmp"))
    file.rename(paste0(cache, ".tmp"), cache)
    TRUE
  }, error = function(e) { message(sprintf("  sim-locus rep%d FAILED: %s", REP, conditionMessage(e))); FALSE })
  invisible(out)
}
reps <- seq_len(fst_prev$n_rep)
already_ok <- vapply(reps, function(REP) {
  cache <- file.path(SIM_LOCUS_CACHE, sprintf("rep%d.rds", REP))
  if (!file.exists(cache)) return(FALSE)
  ck <- tryCatch(readRDS(cache), error = function(e) NULL)
  !is.null(ck) && !is.null(ck$fingerprint) && identical(ck$fingerprint[c("marker_hash", "params_hash", "script_version")], run_fp)
}, logical(1))
todo <- reps[!already_ok]
message(sprintf("[fst-violin] per-locus sim: %d/%d reps already cached, %d to compute (workers=%d)",
                sum(already_ok), length(reps), length(todo), WORKERS))
if (length(todo)) invisible(mclapply(todo, process_rep_locus, mc.cores = WORKERS, mc.preschedule = FALSE))

done_files <- list.files(SIM_LOCUS_CACHE, pattern = "^rep[0-9]+\\.rds$", full.names = TRUE)
stopifnot("not all per-locus sim replicates completed" = length(done_files) == length(reps))
simO_locus <- unlist(lapply(done_files, function(f) readRDS(f)$fst_locus), use.names = FALSE)
n_nonfinite_sim <- sum(!is.finite(simO_locus))
simO_locus <- simO_locus[is.finite(simO_locus)]
message(sprintf("[fst-violin] per-locus simulated Fst: %d values from %d replicates x ~%d loci (%d non-finite excluded; median %.4f, cf. old pooled-per-replicate median %.4f)",
                length(simO_locus), length(reps), nrow(ov), n_nonfinite_sim, median(simO_locus), neutral["med"]))

## ---- parental MAF, primary gate (identical convention to di25_fst_vs_di_rho05.R) --
ep <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = ep)
par_rows <- grepl("_parent$", ep$sample_data_with_parents$Population)
## USER FIX (2026-09-13): fold to minor-allele frequency -- pf alone is a raw
## allele frequency (range [0,1]) and was letting e.g. pf=0.95 (true MAF 0.05)
## pass the "MAF>=0.15" gate. Same bug/fix as di25_fst_vs_di_rho05.R.
pf <- colMeans(ep$GTs_with_parents[par_rows, units$best, drop = FALSE], na.rm = TRUE) / 2
pf_folded <- pmin(pf, 1 - pf)
units[, pmaf := pf_folded[match(best, names(pf_folded))]]
rm(ep); gc()
units_primary <- units[pmaf >= MIN_PARENT_MAF]
message(sprintf("[fst-violin] %d units total -> %d retained at parental MAF>=%.2f (primary gate)",
                nrow(units), nrow(units_primary), MIN_PARENT_MAF))

## ---- per-locus Fst (UNPOOLED -- this is the only difference from the pooled curve) --
pop_emp <- sd$Population[match(rownames(G0), sd$Sample_ID)]
Gb      <- G0[, units_primary$best]
ac_emp  <- wc_ac(Gb, pop_emp)
units_primary[, fst_locus := ac_emp$a / ac_emp$abc]
n_nonfinite <- sum(!is.finite(units_primary$fst_locus))
if (n_nonfinite > 0) message(sprintf("[fst-violin] %d of %d units have non-finite per-locus Fst (zero-variance loci) -- excluded",
                                     n_nonfinite, nrow(units_primary)))
units_primary <- units_primary[is.finite(fst_locus)]
rm(G0, Gb, ac_emp, e); gc()

units_primary[, DI_bin := factor(BIN_LAB[bin], levels = BIN_LAB)]
n_per_bin <- units_primary[, .N, by = DI_bin]
message("[fst-violin] units per DI bin (primary, MAF-gated):")
print(n_per_bin[order(DI_bin)])

saveRDS(list(units = units_primary, neutral = neutral, simO_locus = simO_locus, n_sim_reps = length(reps),
             n_sim_overlap_units = nrow(ov), di_breaks = DI_BREAKS,
             min_parent_maf_primary = MIN_PARENT_MAF, n_per_bin = n_per_bin), OUTRDS)

## ---- figure: simulated (neutral) PER-LOCUS violin/box to the LEFT of the empirical DI-bin ones --
SIM_LAB <- "Simulated\n(neutral, per-locus)"
plot_dt <- rbind(
  data.table(DI_bin = factor(SIM_LAB, levels = c(SIM_LAB, BIN_LAB)), fst_locus = simO_locus, kind = "Simulated"),
  units_primary[, .(DI_bin = factor(as.character(DI_bin), levels = c(SIM_LAB, BIN_LAB)), fst_locus, kind = "Empirical")]
)
p <- ggplot(plot_dt, aes(DI_bin, fst_locus, fill = kind, colour = kind)) +
  geom_violin(alpha = 0.25, scale = "width", linewidth = 0.4, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, outlier.alpha = 0.3, fill = "white", linewidth = 0.4) +
  geom_vline(xintercept = 1.5, linetype = 3, colour = "grey50") +
  scale_fill_manual(values = c(Simulated = "#1b9e77", Empirical = "#d95f02"), guide = "none") +
  scale_colour_manual(values = c(Simulated = "#1b9e77", Empirical = "#d95f02"), guide = "none") +
  labs(x = "DiagnosticIndex bin  (left = near-neutral background, right = ancestry-informative)",
       y = expression("per-locus " * F[ST] * "  (Weir & Cockerham, unpooled)"),
       title = "Per-locus Fst distribution by DI bin (rho05, all LD-reduced units)",
       subtitle = sprintf("%s empirical units (MAF>=%.2f); leftmost = %s per-locus sim values (%d reps)",
                          format(nrow(units_primary), big.mark = ","), MIN_PARENT_MAF,
                          format(length(simO_locus), big.mark = ","), length(reps))) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(8, 12, 4, 6))
ggsave(OUTPDF, p, width = 9.5, height = 5.8)
ggsave(OUTPNG, p, width = 9.5, height = 5.8, dpi = 200)
cat("saved:", OUTPNG, "\n")
