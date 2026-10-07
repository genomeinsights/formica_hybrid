## =========================================================================
## module_allele_specific_sorting -- shared statistics for empirical vs simulated data
## (explore_mosaic_sims.R, explore_calib_sims.R). Uses globals defined by the caller:
## u25 (DI25 units), sets (list(DI25, neutral) marker tables), sgn (DI25 orientation),
## all_markers (marker order for read_sim).
## =========================================================================
## ---- statistics for one data set ----------------------------------------------------------------
wc_fst <- function(G, pop) {
  levs <- unique(pop); M <- ncol(G); N <- P2 <- H <- matrix(0, length(levs), M)
  for (k in seq_along(levs)) { gk <- G[pop == levs[k], , drop = FALSE]; n <- colSums(!is.na(gk)); N[k, ] <- n
    P2[k, ] <- ifelse(n > 0, colSums(gk, na.rm = TRUE) / (2 * n), 0); H[k, ] <- ifelse(n > 0, colSums(gk == 1, na.rm = TRUE) / n, 0) }
  C <- colSums(N); r <- colSums(N > 0); nbar <- C / r; nc <- (C - colSums(N^2) / C) / (r - 1)
  pbar <- colSums(N * P2) / C; hbar <- colSums(N * H) / C
  s2 <- colSums(N * sweep(P2, 2, pbar)^2) / ((r - 1) * nbar); msp <- pbar * (1 - pbar)
  a <- (nbar / nc) * (s2 - (1 / (nbar - 1)) * (msp - ((r - 1) / r) * s2 - 0.25 * hbar))
  b <- (nbar / (nbar - 1)) * (msp - ((r - 1) / r) * s2 - ((2 * nbar - 1) / (4 * nbar)) * hbar)
  list(per_locus = ifelse(a + b + 0.5 * hbar > 0, a / (a + b + 0.5 * hbar), NA_real_), pooled = sum(a, na.rm = TRUE) / sum(a + b + 0.5 * hbar, na.rm = TRUE))
}
half <- function(tot, parts) (tot - parts) / 2
## G: individuals x markers of one set (dosage), H25: individuals x DI25 units oriented dosage (hybrid index)
profile <- function(G, pop, set, H25) {
  pops <- sort(unique(pop))
  Fm <- t(sapply(pops, function(p) colMeans(G[pop == p, , drop = FALSE], na.rm = TRUE) / 2))
  Fbar <- colMeans(Fm, na.rm = TRUE)
  g <- apply(Fm, 2, function(x) mean((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE)) / (Fbar * (1 - Fbar)); g[!is.finite(g)] <- NA
  Za <- among_pop_Z(Fm)
  ## ancestry-adjusted within-population LD (hybrid index from DI25 units on other chromosomes)
  X <- G / 2
  for (p in pops) { r <- pop == p; X[r, ] <- sweep(X[r, , drop = FALSE], 2, colMeans(X[r, , drop = FALSE], na.rm = TRUE)) }
  obs <- !is.na(H25); S_all <- rowSums(H25 / 2, na.rm = TRUE); N_all <- rowSums(obs)
  for (ch in unique(set$Chr)) {
    u <- which(u25$Chr == ch)
    h <- (S_all - rowSums(H25[, u, drop = FALSE] / 2, na.rm = TRUE)) / (N_all - rowSums(obs[, u, drop = FALSE]))
    for (p in pops) { r <- pop == p; h[r] <- h[r] - mean(h[r]) }
    cc <- which(set$Chr == ch); Xc <- X[, cc, drop = FALSE]; ok <- !is.na(Xc); Xc0 <- Xc; Xc0[!ok] <- 0
    X[, cc] <- Xc - outer(h, colSums(Xc0 * h) / colSums(ok * h^2))
  }
  X[is.na(X)] <- 0; nw <- sqrt(colSums(X^2)); nw[nw == 0] <- NA; Zw <- sweep(X, 2, nw, "/")
  cols <- split(seq_len(nrow(set)), factor(set$Chr, levels = unique(set$Chr)))
  w <- rbindlist(lapply(cols, function(ix) {
    m <- length(ix); ut <- which(upper.tri(diag(m))); i <- ix[row(diag(m))[ut]]; j <- ix[col(diag(m))[ut]]
    dcm <- abs(set$cM_pos[j] - set$cM_pos[i]); cm <- cut(dcm, CM_BREAKS, labels = FALSE, right = FALSE); cm[which(dcm == 0)] <- NA
    ca <- crossprod(Za[, ix, drop = FALSE])[ut]; w2 <- crossprod(Zw[, ix, drop = FALSE])[ut]^2; gg <- g[i] * g[j]
    num <- ca^2 * gg; ok <- !is.na(num)
    data.table(bin = cm, num = ifelse(ok, num, 0), den = ifelse(ok, gg, 0), ceil = ifelse(is.na(gg), 0, gg), nceil = !is.na(gg),
               sw = ifelse(is.na(w2), 0, w2), nw = !is.na(w2))[!is.na(bin), lapply(.SD, sum), by = bin]
  }))[, lapply(.SD, sum), by = bin]
  A0 <- Za; A0[is.na(A0)] <- 0; okA <- colSums(is.na(Za)) == 0; gw <- ifelse(is.na(g) | !okA, 0, g)
  M <- lapply(cols, function(ix) A0[, ix, drop = FALSE] %*% (gw[ix] * t(A0[, ix, drop = FALSE])))
  gs <- sapply(cols, function(ix) sum(gw[ix])); ng <- sapply(cols, function(ix) sum(!is.na(g[ix])))
  W0 <- Zw; W0[is.na(W0)] <- 0; okW <- colSums(is.na(Zw)) == 0
  Mw <- lapply(cols, function(ix) tcrossprod(W0[, ix, drop = FALSE])); nwc <- sapply(cols, function(ix) sum(okW[ix]))
  unl <- data.table(bin = 99L, num = half(sum(Reduce(`+`, M)^2), sum(sapply(M, function(m) sum(m^2)))),
                    den = half(sum(gs)^2, sum(gs^2)), ceil = half(sum(gs)^2, sum(gs^2)), nceil = half(sum(ng)^2, sum(ng^2)),
                    sw = half(sum(Reduce(`+`, Mw)^2), sum(sapply(Mw, function(m) sum(m^2)))), nw = half(sum(nwc)^2, sum(nwc^2)))
  b <- 1 / (length(pops) - 1)
  rbind(w, unl)[, .(bin, ceiling = ceil / nceil, realised = num / den, share = (num / den - b) / (1 - b), within_LD = sw / nw)]
}
sorting_pct <- function(Gor, pop) {
  pops <- sort(unique(pop))
  Fm <- t(sapply(pops, function(p) colMeans(Gor[pop == p, , drop = FALSE], na.rm = TRUE) / 2))
  n_aqu <- colSums(Fm >= 0.85, na.rm = TRUE); n_pol <- colSums(Fm <= 0.15, na.rm = TRUE); n_obs <- colSums(!is.na(Fm))
  cls <- classify_sort(n_aqu, n_pol, n_obs, sort_th = 0.6, sort_rule = "binom", alpha = 0.05)
  c(pct_sorted = 100 * mean(cls %in% c("aquilonia", "polyctena")), pct_aqu_of_sorted = 100 * mean(cls[cls %in% c("aquilonia", "polyctena")] == "aquilonia"))
}
all_stats <- function(G25raw, Gneu, pop, label) {
  Gor <- G25raw; Gor[, which(sgn < 0)] <- 2L - Gor[, which(sgn < 0)]; Gor[, which(is.na(sgn) | sgn == 0)] <- NA
  f25 <- wc_fst(G25raw, pop); fne <- wc_fst(Gneu, pop)
  s <- sorting_pct(Gor, pop)
  summ <- data.table(data = label, fst_DI25_median = median(f25$per_locus, na.rm = TRUE), fst_DI25_pooled = f25$pooled,
                     fst_neutral_median = median(fne$per_locus, na.rm = TRUE), fst_neutral_pooled = fne$pooled,
                     pct_sorted = s[["pct_sorted"]], pct_aqu_of_sorted = s[["pct_aqu_of_sorted"]])
  prof <- rbind(profile(Gor, pop, sets$DI25, Gor)[, set := "DI25"], profile(Gneu, pop, sets$neutral, Gor)[, set := "neutral"])
  list(summary = summ, profile = prof[, data := label])
}

read_sim <- function(f, n) {
  v <- fread(cmd = sprintf("gzip -dc %s | grep -v '^##'", f), sep = "\t", colClasses = "character")
  mk <- paste0("Chr", sub("^ch", "", v[["#CHROM"]]), ":", v$POS)
  keep <- which(mk %in% all_markers & grepl("MT=10", v$INFO, fixed = TRUE))
  smp <- names(v)[-(1:9)]; smp <- sample(smp, min(n, length(smp)))
  D <- sapply(smp, function(s) { x <- v[[s]][keep]; as.integer(substr(x, 1, 1)) + as.integer(substr(x, 3, 3)) })
  G <- matrix(0L, length(smp), length(all_markers), dimnames = list(smp, all_markers))   # absent site = allele 0
  G[, mk[keep]] <- t(D)
  G
}
