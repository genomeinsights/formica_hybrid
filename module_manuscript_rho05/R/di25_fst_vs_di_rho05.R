## =============================================================================
## module_manuscript_rho05 -- Fig [FST_simulations], rho05 version
## =============================================================================
## Adapted from module_di25/R/di25_fst_vs_di.R. Swaps the FULL-data LD-reduced
## unit set from module0's canonical (min_r2=0.2) best-SNP object to the
## rho05 (min_r2_rho=0.5) one -- 661,386 vs the old 1,114,423-marker-derived
## clustering means a DIFFERENT set of best-SNP markers, so the per-replicate
## simulation cache CANNOT be reused (ov$best differs) -- this reruns the
## 1000-replicate Fst computation against the new unit set, cached separately
## under fst_sim_cache_full_rho05/.
##
## Everything else (fixed DI bins, MAF stratification, background LD,
## Weir & Cockerham estimator) is unchanged from the canonical script.
##
## Run from the repo root:
##   Rscript module_manuscript_rho05/R/di25_fst_vs_di_rho05.R [nreps] [workers]   (default 1000, 9)
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(parallel) })

args    <- commandArgs(trailingOnly = TRUE)
NREPS   <- if (length(args) >= 1) as.integer(args[1]) else 1000L
WORKERS <- if (length(args) >= 2) as.integer(args[2]) else 9L
BEDFMT    <- "data/diem_outs_demo/diem_boot%d_output.bed"
CACHE_DIR <- "module_manuscript_rho05/data/fst_sim_cache_full_rho05"
OUTRDS    <- "module_manuscript_rho05/data/di25_fst_vs_di_rho05.rds"
FIGDIR    <- "module_manuscript_rho05/Figures"
OUTPNG    <- file.path(FIGDIR, "di25_fst_vs_di_rho05.png")
OUTPDF    <- sub("\\.png$", ".pdf", OUTPNG)
DI_BREAKS <- c(-Inf, -90, -75, -60, -50, -40, -30, -25, -20, -15, Inf)
BIN_LAB   <- levels(cut(0, DI_BREAKS)); N_BIN <- length(BIN_LAB)
BG_NSUB   <- 2000L; BG_Q <- 0.95; MIN_SIM_UNITS <- 30L
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

## ---- FULL-data LD-reduced units (best SNPs), rho05 + FIXED DI bins ----------
b   <- readRDS("module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05_bestsnp.rds")
e   <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e)
G0  <- e$GTs_hybrids_005
mDI <- as.data.table(e$map_hyb_005)
sd  <- as.data.table(e$sample_data)
units <- data.table(best = b$stats$best_marker,
                    DI = mDI$DiagnosticIndex[match(b$stats$best_marker, mDI$marker)])
units <- units[is.finite(DI)]
units[, bin := cut(DI, DI_BREAKS, labels = FALSE)]
message(sprintf("[fst-rho05] %d full-data best-SNP units (rho05); %d fixed DI bins (%s .. %s)",
                nrow(units), N_BIN, BIN_LAB[1], BIN_LAB[N_BIN]))

## ---- empirical Fst per fixed DI bin + background LD --------------------------
pop_emp <- sd$Population[match(rownames(G0), sd$Sample_ID)]
Gb      <- G0[, units$best]
ac_emp  <- wc_ac(Gb, pop_emp)
emp_fst <- fst_by_bin(ac_emp, units$bin)
message(sprintf("[fst-rho05] empirical overall Fst = %.3f (neutral bg DI<-90 %.3f -> diagnostic %.3f)",
                fst_ratio(ac_emp), emp_fst[1], emp_fst[N_BIN]))
emp_chr <- as.integer(sub("Chr", "", sub(":.*", "", colnames(G0))))
set.seed(1); emp_bg <- bg_ld(G0, emp_chr)
message(sprintf("[fst-rho05] empirical background LD = %.4f", emp_bg))

## ---- empirical Fst stratified by POOLED PARENTAL MAF -------------------------
ep <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = ep)
par_rows <- grepl("_parent$", ep$sample_data_with_parents$Population)
pf <- colMeans(ep$GTs_with_parents[par_rows, units$best, drop = FALSE], na.rm = TRUE) / 2
units[, pmaf := pmin(pf, 1 - pf)]
rm(ep); gc()
MAF_BREAKS <- c(0, 0.1, 0.2, 0.35, 0.5); MAF_LAB <- levels(cut(0, MAF_BREAKS, include.lowest = TRUE))
units[, mstr := cut(pmaf, MAF_BREAKS, include.lowest = TRUE)]
MIN_STRAT <- 30L
strat <- rbindlist(lapply(MAF_LAB, function(ms) data.table(mstr = ms, bin = seq_len(N_BIN),
  n   = vapply(seq_len(N_BIN), function(d) sum(!is.na(units$mstr) & units$mstr == ms & units$bin == d), integer(1)),
  fst = vapply(seq_len(N_BIN), function(d) fst_ratio(ac_emp, !is.na(units$mstr) & units$mstr == ms & units$bin == d), numeric(1)))))
strat <- strat[n >= MIN_STRAT]; strat[, mstr := factor(mstr, levels = MAF_LAB)]
message(sprintf("[fst-rho05] parental-MAF strata: %s", paste(MAF_LAB, collapse = ", ")))

rm(G0, Gb, ac_emp, e); gc()

## ---- sim: best-SNP units that fall in the DI25 sim panel ---------------------
sim_mk  <- { hdr <- readLines(sprintf(BEDFMT, 1), n = 2); s1 <- fread(sprintf(BEDFMT, 1), skip = 2,
              header = FALSE, sep = "\t", select = c(1, 3), colClasses = list(character = 1), showProgress = FALSE)
            paste0("Chr", sub("ch", "", s1$V1), ":", s1$V3) }
ov <- units[best %in% sim_mk]
message(sprintf("[fst-rho05] %d best-SNP units overlap the DI25 sim panel (bins %s)",
                nrow(ov), paste(range(ov$bin), collapse = "-")))

process_rep <- function(REP) {
  cache <- file.path(CACHE_DIR, sprintf("rep%d.rds", REP)); if (file.exists(cache)) return(invisible(TRUE))
  bed <- sprintf(BEDFMT, REP); if (!file.exists(bed)) return(invisible(FALSE))
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
    saveRDS(list(fst = fst, overall = overall, bg = bg, rep = REP), cache); TRUE
  }, error = function(e) { message(sprintf("  rep%d FAILED: %s", REP, conditionMessage(e))); FALSE })
  invisible(out)
}
reps <- seq_len(NREPS)
todo <- reps[!file.exists(file.path(CACHE_DIR, sprintf("rep%d.rds", reps)))]
message(sprintf("[fst-rho05] %d reps requested, %d cached, %d to compute (workers=%d)", NREPS, NREPS - length(todo), length(todo), WORKERS))
if (length(todo)) invisible(mclapply(todo, process_rep, mc.cores = WORKERS, mc.preschedule = FALSE))

## ---- aggregate --------------------------------------------------------------
done <- reps[file.exists(file.path(CACHE_DIR, sprintf("rep%d.rds", reps)))]
sm   <- lapply(done, function(r) readRDS(file.path(CACHE_DIR, sprintf("rep%d.rds", r))))
simF <- do.call(rbind, lapply(sm, function(x) x$fst))
simO <- vapply(sm, function(x) x$overall, numeric(1))
simB <- vapply(sm, function(x) x$bg, numeric(1))
n_emp <- units[, .N, by = bin][order(bin)]
n_ov  <- ov[, .N, by = bin][order(bin)]
env  <- data.table(bin = seq_len(N_BIN), DI_bin = BIN_LAB,
                   emp = emp_fst, n_emp_units = n_emp$N[match(seq_len(N_BIN), n_emp$bin)],
                   n_sim_units = n_ov$N[match(seq_len(N_BIN), n_ov$bin)],
                   sim_med = apply(simF, 2, median, na.rm = TRUE),
                   sim_lo  = apply(simF, 2, quantile, 0.025, na.rm = TRUE),
                   sim_hi  = apply(simF, 2, quantile, 0.975, na.rm = TRUE))
env[is.na(n_sim_units) | n_sim_units < MIN_SIM_UNITS, `:=`(sim_med = NA, sim_lo = NA, sim_hi = NA)]
neutral <- c(med = median(simO), lo = quantile(simO, 0.025), hi = quantile(simO, 0.975))
saveRDS(list(env = env, strat = strat, neutral = neutral, emp_bg = emp_bg, sim_bg = simB,
             n_rep = length(done), di_breaks = DI_BREAKS, maf_breaks = MAF_BREAKS), OUTRDS)

cat(sprintf("\n=== [rho05] among-population Fst by DI bin (empirical vs %d-rep high-DI neutral sim) ===\n", length(done)))
print(env[, .(DI_bin, n_emp_units, emp = round(emp, 3), n_sim_units, sim_med = round(sim_med, 3))])
cat(sprintf("neutral sim Fst (over its high-DI units): %.3f [%.3f, %.3f]\n", neutral[1], neutral[2], neutral[3]))
cat(sprintf("background LD (inter-chr r^2, q%.2f): empirical %.4f | sim %.4f [%.4f, %.4f]\n",
            BG_Q, emp_bg, median(simB), quantile(simB, 0.025), quantile(simB, 0.975)))

## ---- figure -----------------------------------------------------------------
p <- ggplot(env, aes(bin)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = neutral[2], ymax = neutral[3], fill = "#66c2a5", alpha = 0.30) +
  geom_hline(yintercept = neutral[1], colour = "#1b9e77", linetype = 2, linewidth = 0.7) +
  geom_line(data = strat, aes(bin, fst, linetype = mstr), colour = "#d95f02", linewidth = 0.55, na.rm = TRUE) +
  geom_line(aes(y = emp, colour = "empirical (all SNPs)"), linewidth = 1.0) +
  geom_point(aes(y = emp, colour = "empirical (all SNPs)"), size = 2.4) +
  geom_point(aes(y = sim_med, colour = "high-DI neutral sim"), size = 2.2, na.rm = TRUE) +
  geom_errorbar(aes(ymin = sim_lo, ymax = sim_hi, colour = "high-DI neutral sim"), width = 0.15, na.rm = TRUE) +
  annotate("text", x = 1, y = neutral[1], label = "neutral sim (high-DI)", colour = "#1b7f63",
           hjust = 0, vjust = -0.6, size = 3.4) +
  scale_x_continuous(breaks = seq_len(N_BIN), labels = BIN_LAB) +
  scale_colour_manual(values = c("empirical (all SNPs)" = "#d95f02", "high-DI neutral sim" = "#1b9e77"), name = NULL) +
  scale_linetype_manual(values = c("[0,0.1]" = "dotted", "(0.1,0.2]" = "dotdash",
                                   "(0.2,0.35]" = "longdash", "(0.35,0.5]" = "22"), name = "parental MAF") +
  labs(x = "DiagnosticIndex bin  (left = near-neutral background, right = ancestry-informative)",
       y = expression("among-population " * F[ST] * "  (Weir & Cockerham), rho05 units")) +
  guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2, nrow = 1)) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), legend.position = "top", legend.box = "vertical",
        legend.margin = margin(1, 1, 1, 1), legend.spacing.y = unit(1, "pt"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(8, 12, 4, 6))
ggsave(OUTPDF, p, width = 8.6, height = 5.6); ggsave(OUTPNG, p, width = 8.6, height = 5.6, dpi = 200)
cat("saved:", OUTPNG, "\n")
