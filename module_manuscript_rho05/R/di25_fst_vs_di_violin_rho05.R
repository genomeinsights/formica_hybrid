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
suppressMessages({ library(data.table); library(ggplot2) })

FST_RDS   <- "module_manuscript_rho05/data/di25_fst_vs_di_rho05.rds"
OUTRDS    <- "module_manuscript_rho05/data/di25_fst_vs_di_violin_rho05.rds"
FIGDIR    <- "module_manuscript_rho05/Figures"
OUTPNG    <- file.path(FIGDIR, "di25_fst_vs_di_violin_rho05.png")
OUTPDF    <- sub("\\.png$", ".pdf", OUTPNG)
DI_BREAKS <- c(-Inf, -90, -75, -60, -50, -40, -30, -25, -20, -15, Inf)
BIN_LAB   <- levels(cut(0, DI_BREAKS)); N_BIN <- length(BIN_LAB)
MIN_PARENT_MAF <- 0.15
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

fst_prev <- readRDS(FST_RDS)
stopifnot("MIN_PARENT_MAF must match the audited di25_fst_vs_di_rho05.R" =
            fst_prev$min_parent_maf_primary == MIN_PARENT_MAF,
          identical(fst_prev$di_breaks, DI_BREAKS))
neutral <- fst_prev$neutral   # med/lo/hi, from the 1000-rep high-DI neutral sim (scalar, not per-bin)

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

## ---- parental MAF, primary gate (identical convention to di25_fst_vs_di_rho05.R) --
ep <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = ep)
par_rows <- grepl("_parent$", ep$sample_data_with_parents$Population)
pf <- colMeans(ep$GTs_with_parents[par_rows, units$best, drop = FALSE], na.rm = TRUE) / 2
units[, pmaf := pf[match(best, names(pf))]]
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

saveRDS(list(units = units_primary, neutral = neutral, di_breaks = DI_BREAKS,
             min_parent_maf_primary = MIN_PARENT_MAF, n_per_bin = n_per_bin), OUTRDS)

## ---- figure: violin + boxplot per DI bin, neutral-sim band for reference --
p <- ggplot(units_primary, aes(DI_bin, fst_locus)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = neutral["lo"], ymax = neutral["hi"],
           fill = "#66c2a5", alpha = 0.30) +
  geom_hline(yintercept = neutral["med"], colour = "#1b9e77", linetype = 2, linewidth = 0.6) +
  geom_violin(fill = "#d95f02", alpha = 0.25, colour = "#d95f02", scale = "width", linewidth = 0.4, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, outlier.alpha = 0.3, fill = "white", linewidth = 0.4) +
  annotate("text", x = 1, y = neutral["med"], label = "neutral sim (high-DI)", colour = "#1b7f63",
           hjust = 0, vjust = -0.6, size = 3.2) +
  labs(x = "DiagnosticIndex bin  (left = near-neutral background, right = ancestry-informative)",
       y = expression("per-locus " * F[ST] * "  (Weir & Cockerham, unpooled)"),
       title = "Per-locus Fst distribution by DI bin (rho05, all LD-reduced units)",
       subtitle = sprintf("%s units, parental MAF>=%.2f; shaded = neutral high-DI simulation 95%% interval (pooled, %d reps)",
                          format(nrow(units_primary), big.mark = ","), MIN_PARENT_MAF, fst_prev$n_rep)) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(8, 12, 4, 6))
ggsave(OUTPDF, p, width = 9.5, height = 5.8)
ggsave(OUTPNG, p, width = 9.5, height = 5.8, dpi = 200)
cat("saved:", OUTPNG, "\n")
