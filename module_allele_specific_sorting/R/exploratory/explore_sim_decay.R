## =========================================================================
## module_allele_specific_sorting -- EXPLORATORY (not in the document, not committed):
## Figure-1 statistics on the EXISTING neutral simulations (data/diem_outs_demo).
##
## CAVEAT: these replicates were generated with the founder-frequency indexing bug
## (sim_founder_fix/NOTE_for_Beatriz.md): founders on every chromosome except Chr1
## received Chr1's fixation levels. This is a preview of the shape, not a valid null.
##
## Per replicate, exactly as for the empirical data (01-03):
##   units = the empirical 20,807 representative SNPs (present in the simulations);
##   orientation from the replicate's OWN simulated parents (aq_/pol_ samples);
##   Fmat, G (Nei), W&C F_ST, LOCO-ancestry-residualised profiles (pp_residualize_ancestry.R
##   logic), ancestry-adjusted within-population LD (00_utils.R::within_pop_Z(adj = TRUE));
##   within-chromosome pairs binned by cM (zero-map-distance pairs excluded) and bp;
##   unlinked class via the exact algebraic shortcut (02);
##   realised-fraction baseline from per-unit population-label permutations.
## Output: data/explore_sim_decay.rds, Figures/explore_sim_vs_emp_decay_cM.png
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/exploratory/explore_sim_decay.R [N_REPS] [N_CORES] [N_PERM]
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
args <- commandArgs(trailingOnly = TRUE)
N_REPS  <- if (length(args) >= 1) as.integer(args[1]) else 20L
N_CORES <- if (length(args) >= 2) as.integer(args[2]) else 8L
N_PERM  <- if (length(args) >= 3) as.integer(args[3]) else 5L
BEDFMT  <- "data/diem_outs_demo/diem_boot%d_output.bed"

u0 <- readRDS(file.path(OUT_DATA, "01_units.rds"))[, .(group_id, Chr, Pos, cM_pos, unit_marker)]

read_rep <- function(r) {
  f <- sprintf(BEDFMT, r)
  h <- readLines(f, n = 2)
  inds <- strsplit(strsplit(h[2], "\t")[[1]][10], "|", fixed = TRUE)[[1]]
  s <- fread(f, skip = 2, header = FALSE, sep = "\t", select = c(1, 3, 10),
             colClasses = list(character = c(1, 10)), showProgress = FALSE)
  mk <- paste0("Chr", sub("ch", "", s$V1), ":", s$V3)
  keep <- match(u0$unit_marker, mk)
  S <- do.call(rbind, strsplit(sub("^S", "", s$V10[keep[!is.na(keep)]]), "", fixed = TRUE))
  G <- matrix(NA_integer_, length(inds), nrow(u0))
  G[, !is.na(keep)] <- t(matrix(suppressWarnings(as.integer(S)), nrow(S)))
  list(G = G, inds = inds)
}

wc_fst <- function(G, pop) {
  levs <- unique(pop); M <- ncol(G)
  N <- P2 <- H <- matrix(0, length(levs), M)
  for (k in seq_along(levs)) {
    gk <- G[pop == levs[k], , drop = FALSE]
    n <- colSums(!is.na(gk)); N[k, ] <- n
    P2[k, ] <- ifelse(n > 0, colSums(gk, na.rm = TRUE) / (2 * n), 0); H[k, ] <- ifelse(n > 0, colSums(gk == 1, na.rm = TRUE) / n, 0)
  }
  C <- colSums(N); r <- colSums(N > 0); nbar <- C / r; nc <- (C - colSums(N^2) / C) / (r - 1)
  pbar <- colSums(N * P2) / C; hbar <- colSums(N * H) / C
  s2 <- colSums(N * sweep(P2, 2, pbar)^2) / ((r - 1) * nbar); msp <- pbar * (1 - pbar)
  a <- (nbar / nc) * (s2 - (1 / (nbar - 1)) * (msp - ((r - 1) / r) * s2 - 0.25 * hbar))
  b <- (nbar / (nbar - 1)) * (msp - ((r - 1) / r) * s2 - ((2 * nbar - 1) / (4 * nbar)) * hbar)
  ifelse(a + b + 0.5 * hbar > 0, a / (a + b + 0.5 * hbar), NA_real_)
}

chrs <- unique(u0$Chr); cols <- split(seq_len(nrow(u0)), factor(u0$Chr, levels = chrs))
pairs <- lapply(cols, function(ix) {
  m <- length(ix); ut <- which(upper.tri(diag(m)))
  i <- ix[row(diag(m))[ut]]; j <- ix[col(diag(m))[ut]]
  dcm <- abs(u0$cM_pos[j] - u0$cM_pos[i])
  cm <- cut(dcm, CM_BREAKS, labels = FALSE, right = FALSE); cm[which(dcm == 0)] <- NA
  list(ut = ut, i = i, j = j, cm = cm)
})

one_rep <- function(r) {
  x <- read_rep(r); inds <- x$inds
  hyb <- grepl("^hyb_", inds); pop <- sub("^hyb_(.*)_[0-9]+$", "\\1", inds[hyb])
  sgn <- sign(colMeans(x$G[grepl("^aq_", inds), , drop = FALSE], na.rm = TRUE) -
              colMeans(x$G[grepl("^pol_", inds), , drop = FALSE], na.rm = TRUE))
  G <- x$G[hyb, , drop = FALSE]; G[, which(sgn < 0)] <- 2L - G[, which(sgn < 0)]; G[, which(is.na(sgn) | sgn == 0)] <- NA
  pops <- sort(unique(pop))
  Fm <- t(sapply(pops, function(p) colMeans(G[pop == p, , drop = FALSE], na.rm = TRUE) / 2))
  fst <- wc_fst(G, pop)
  Fbar <- colMeans(Fm); g <- apply(Fm, 2, function(v) mean((v - mean(v))^2)) / (Fbar * (1 - Fbar))
  g[!is.finite(g)] <- NA
  ## LOCO-ancestry residualised profiles (lm with intercept on H_loco, per unit)
  Res <- matrix(NA_real_, nrow(Fm), ncol(Fm))
  S_p <- rowSums(Fm, na.rm = TRUE); N_p <- rowSums(!is.na(Fm))
  for (k in seq_along(cols)) {
    ix <- cols[[k]]; Hc <- (S_p - rowSums(Fm[, ix, drop = FALSE], na.rm = TRUE)) / (N_p - rowSums(!is.na(Fm[, ix, drop = FALSE])))
    X <- cbind(1, Hc); Pm <- diag(nrow(X)) - X %*% solve(crossprod(X), t(X))
    Fx <- Fm[, ix, drop = FALSE]; ok <- colSums(is.na(Fx)) == 0
    Res[, ix[ok]] <- Pm %*% Fx[, ok, drop = FALSE]
  }
  Za <- among_pop_Z(Fm); Zr <- among_pop_Z(Res)
  Zw <- within_pop_Z(list(G = G, pop = pop), u0, adj = TRUE)

  within <- function(A, R, W) rbindlist(lapply(seq_along(cols), function(k) {
    p <- pairs[[k]]; ix <- cols[[k]]
    ca <- crossprod(A[, ix, drop = FALSE])[p$ut]; cr <- crossprod(R[, ix, drop = FALSE])[p$ut]
    w2 <- if (is.null(W)) NA_real_ else crossprod(W[, ix, drop = FALSE])[p$ut]^2
    gg <- g[p$i] * g[p$j]; num <- ca^2 * gg; okn <- !is.na(num)
    data.table(cm = p$cm, num = ifelse(okn, num, 0), den = ifelse(okn, gg, 0), ceil = ifelse(is.na(gg), 0, gg), nceil = !is.na(gg),
               sc = ifelse(is.na(ca), 0, ca), nc = !is.na(ca), sr = ifelse(is.na(cr), 0, cr), nr = !is.na(cr),
               sw = ifelse(is.na(w2), 0, w2), nw = !is.na(w2))[!is.na(cm), lapply(.SD, sum), by = cm]
  }))[, lapply(.SD, sum), by = cm]
  unlinked <- function(A, R, W) {
    half <- function(tot, parts) (tot - parts) / 2
    A0 <- A; A0[is.na(A0)] <- 0; R0 <- R; R0[is.na(R0)] <- 0
    okA <- colSums(is.na(A)) == 0; okR <- colSums(is.na(R)) == 0
    gw <- ifelse(is.na(g) | !okA, 0, g)
    M  <- lapply(cols, function(ix) A0[, ix, drop = FALSE] %*% (gw[ix] * t(A0[, ix, drop = FALSE])))
    sa <- lapply(cols, function(ix) rowSums(A0[, ix, drop = FALSE])); sr <- lapply(cols, function(ix) rowSums(R0[, ix, drop = FALSE]))
    gs <- sapply(cols, function(ix) sum(gw[ix])); na <- sapply(cols, function(ix) sum(okA[ix])); nr <- sapply(cols, function(ix) sum(okR[ix]))
    out <- data.table(cm = 99L,
      num = half(sum(Reduce(`+`, M)^2), sum(sapply(M, function(m) sum(m^2)))),
      den = half(sum(gs)^2, sum(gs^2)), ceil = half(sum(gs)^2, sum(gs^2)), nceil = half(sum(na)^2, sum(na^2)),
      sc = half(sum(Reduce(`+`, sa)^2), sum(sapply(sa, function(v) sum(v^2)))), nc = half(sum(na)^2, sum(na^2)),
      sr = half(sum(Reduce(`+`, sr)^2), sum(sapply(sr, function(v) sum(v^2)))), nr = half(sum(nr)^2, sum(nr^2)),
      sw = NA_real_, nw = NA_real_)
    if (!is.null(W)) {
      W0 <- W; W0[is.na(W0)] <- 0; okW <- colSums(is.na(W)) == 0
      Mw <- lapply(cols, function(ix) tcrossprod(W0[, ix, drop = FALSE])); nwc <- sapply(cols, function(ix) sum(okW[ix]))
      sw_u <- half(sum(Reduce(`+`, Mw)^2), sum(sapply(Mw, function(m) sum(m^2)))); nw_u <- half(sum(nwc)^2, sum(nwc^2))
      out[, `:=`(sw = sw_u, nw = nw_u)]   # locals renamed: inside j, `nw` would mean the column
    }
    out
  }
  agg <- rbind(within(Za, Zr, Zw), unlinked(Za, Zr, Zw))
  ## baseline for the realised fraction: per-unit population-label permutations
  npop <- nrow(Za)
  base <- rowMeans(sapply(seq_len(N_PERM), function(b) {
    P <- replicate(ncol(Za), sample.int(npop))
    pz <- function(Z) matrix(Z[cbind(as.vector(P), rep(seq_len(ncol(Z)), each = npop))], npop)
    a <- rbind(within(pz(Za), pz(Zr), NULL), unlinked(pz(Za), pz(Zr), NULL)); setorder(a, cm)
    a$num / a$den
  }))
  setorder(agg, cm)
  agg[, .(rep = r, cm, ceiling = ceil / nceil, among_LD = num / nceil, realised = num / den, baseline = base,
          realised_share = (num / den - base) / (1 - base), conc = sc / nc, conc_resid = sr / nr,
          within_LD_adj = sw / nw, median_fst = median(fst, na.rm = TRUE), n_sorted_like = NA_integer_)]
}

RNGkind("L'Ecuyer-CMRG"); set.seed(1)
t0 <- Sys.time()
res <- rbindlist(parallel::mclapply(seq_len(N_REPS), one_rep, mc.cores = N_CORES, mc.set.seed = TRUE))
cat(sprintf("[sim] %d replicates in %.1f min\n", N_REPS, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
res[, bin := factor(ifelse(cm == 99L, "unlinked", CM_LABELS[pmin(cm, length(CM_LABELS))]), levels = c(CM_LABELS, "unlinked"))]
long <- melt(res, id.vars = c("rep", "bin"), measure.vars = c("ceiling", "among_LD", "realised", "realised_share",
                                                               "conc", "conc_resid", "within_LD_adj"), variable.name = "qty")
sim <- long[, .(sim_mean = mean(value, na.rm = TRUE), sim_lo = quantile(value, 0.025, na.rm = TRUE),
                sim_hi = quantile(value, 0.975, na.rm = TRUE)), by = .(bin, qty)]
emp <- readRDS(file.path(OUT_DATA, "03_decay.rds"))$cm[, .(bin = factor(as.character(bin), levels = levels(sim$bin)), qty, emp = mean)]
cmp <- merge(sim, emp, by = c("bin", "qty"), all.x = TRUE); setorder(cmp, qty, bin)
cat(sprintf("[sim] median unit F_ST: simulations %.3f (range over replicates %.3f-%.3f); empirical 0.251\n",
            mean(res$median_fst), min(res$median_fst), max(res$median_fst)))
print(dcast(cmp[qty %in% c("ceiling", "realised", "realised_share", "conc_resid", "within_LD_adj")],
            bin ~ qty, value.var = c("emp", "sim_mean")), digits = 3)
saveRDS(list(per_rep = res, compare = cmp, n_reps = N_REPS, n_perm = N_PERM), file.path(OUT_DATA, "explore_sim_decay.rds"))

## ---- figure: empirical vs simulations, same layout as Figure 1 ------------------------
lab <- c(ceiling = "ceiling (max among-pop LD)", among_LD = "among-population LD", within_LD_adj = "within-population LD",
         realised = "realised fraction", conc_resid = "ancestry-profile similarity (residualised)")
d <- melt(cmp[qty %in% names(lab)], id.vars = c("bin", "qty", "sim_lo", "sim_hi"), measure.vars = c("emp", "sim_mean"),
          variable.name = "data")
d[, data := factor(data, levels = c("emp", "sim_mean"), labels = c("empirical", "neutral simulations (buggy founders)"))]
d[data != "neutral simulations (buggy founders)", `:=`(sim_lo = NA, sim_hi = NA)]
mkp <- function(q, logy) {
  p <- ggplot(d[qty %in% q], aes(bin, value, colour = data, linetype = qty, group = interaction(data, qty))) +
    geom_ribbon(aes(ymin = sim_lo, ymax = sim_hi, fill = data), colour = NA, alpha = 0.2, show.legend = FALSE) +
    geom_line() + geom_point(size = 1.2) +
    scale_colour_manual(values = c("black", "#d95f02"), aesthetics = c("colour", "fill"), name = NULL) +
    scale_linetype_discrete(labels = lab, name = NULL) +
    labs(x = "genetic distance between units", y = NULL) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom", legend.box = "vertical")
  if (logy) p + scale_y_log10() else p
}
p <- patchwork::wrap_plots(mkp(c("ceiling", "among_LD", "within_LD_adj"), TRUE) + ggtitle("a  differentiation vs LD (log)"),
                           mkp(c("realised", "conc_resid"), FALSE) + ggtitle("b  shared population partitions"), nrow = 1) +
  patchwork::plot_annotation(caption = sprintf("Simulations: mean and 95%% range over %d replicates; CAVEAT founder-frequency bug in these replicates", N_REPS))
ggsave(file.path(OUT_FIG, "explore_sim_vs_emp_decay_cM.png"), p, width = 12, height = 5.5, dpi = 200)
cat("[sim] done\n")
