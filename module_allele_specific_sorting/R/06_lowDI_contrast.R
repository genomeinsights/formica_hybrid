## =========================================================================
## module_allele_specific_sorting -- 06: near-neutral contrast (DI <= -90)
##
## Drift and demography act on every locus regardless of DI; ancestry sorting can
## only act where the parents differ. Near-neutral loci (parental diagnostic index
## DI <= -90) are therefore an internal control: do they show the same distance
## profile of shared population partitions as the ancestry-informative DI25 units
## (03), or only a lower level of differentiation?
##
## Units: full-genome rho05 LD units (module0_ld_pruning_rho05, best/representative
## SNP) with DI <= -90 and pooled parental minor-allele frequency >= 0.15 (folded),
## the gate of the full-genome F_ST analysis. Orientation to F. aquilonia is
## meaningless when the parents barely differ, so only sign-free statistics are
## used: ceiling G_i G_j, realised fraction R = sum c^2 G_i G_j / sum G_i G_j, its
## share beyond the population-label permutation baseline S = (R - b)/(1 - b), and
## within-population LD r_w^2 adjusted for each individual's leave-one-chromosome-out
## hybrid index (estimated from the DI25 units, as in 01). Same cM bins (pairs with
## identical map positions excluded), unlinked shortcut, 500 permutations and
## chromosome-block bootstrap as 01-03.
##
## Output: data/06_lowDI.rds, Figures/06_lowDI_contrast.{png,pdf}
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/06_lowDI_contrast.R [N_PERM] [N_CORES]
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
args <- commandArgs(trailingOnly = TRUE)
N_PERM  <- if (length(args) >= 1) as.integer(args[1]) else 500L
N_CORES <- if (length(args) >= 2) as.integer(args[2]) else 8L
DI_MAX <- -90; PMAF_MIN <- 0.15; B <- 2000L; SEED <- 1L

## ---- near-neutral units -------------------------------------------------------------
e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
map <- as.data.table(e$map_hyb_005); sdp <- as.data.table(e$sample_data_with_parents)
reps <- readRDS("module0_ld_pruning_rho05/data/eMLG_5loci_0025_cM05_rho05_bestsnp.rds")$rep_snp_all
reps <- reps[, .(group_id, rep_snp, n_loci)]
reps[, DI := map$DiagnosticIndex[match(rep_snp, map$marker)]]
cand <- reps[is.finite(DI) & DI <= DI_MAX]
pops_all <- sdp$Population[match(rownames(e$GTs_with_parents), sdp$Sample_ID)]
par <- grepl("_parent$", pops_all)
Gc <- e$GTs_with_parents[, cand$rep_snp, drop = FALSE]
rm(e); invisible(gc())
pf <- colMeans(Gc[par, , drop = FALSE], na.rm = TRUE) / 2
cand[, pmaf := pmin(pf, 1 - pf)]
lo <- cand[pmaf >= PMAF_MIN]
lo[, Chr := sub(":.*", "", rep_snp)][, Pos := as.integer(sub(".*:", "", rep_snp))]
lo[, ChrNum := as.integer(sub("Chr", "", Chr))]; setorder(lo, ChrNum, Pos)
Gh <- Gc[!par, lo$rep_snp, drop = FALSE]; pop <- pops_all[!par]
cat(sprintf("[06] near-neutral units: %d (DI <= %d: %d; of which pooled parental MAF >= %.2f: %d); %d hybrids in %d populations\n",
            nrow(lo), DI_MAX, nrow(cand), PMAF_MIN, nrow(lo), nrow(Gh), length(unique(pop))))

## genetic positions (00_utils convention: off-map -> NA)
gm <- fread("data/Frufa_DTOL_PR.ref_genome.recmap"); gm[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
lo[, cM_pos := NA_real_]
for (ch in unique(lo$Chr)) {
  g <- gm[Chr == ch & !is.na(cM)]; setorder(g, pos); idx <- which(lo$Chr == ch)
  if (nrow(g) < 2) next
  v <- approxfun(g$pos, g$cM, rule = 2)(lo$Pos[idx]); v[lo$Pos[idx] < min(g$pos) | lo$Pos[idx] > max(g$pos)] <- NA
  lo$cM_pos[idx] <- v
}

## ---- among-population statistics ----------------------------------------------------
hpops <- sort(unique(pop))
Fm <- t(sapply(hpops, function(p) colMeans(Gh[pop == p, , drop = FALSE], na.rm = TRUE) / 2))
Fbar <- colMeans(Fm, na.rm = TRUE)
g <- apply(Fm, 2, function(x) mean((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE)) / (Fbar * (1 - Fbar))
g[!is.finite(g)] <- NA
Za <- among_pop_Z(Fm)

## ---- within-population LD, hybrid index from the DI25 units --------------------------
U25 <- load_units(); GG25 <- load_oriented_genotypes(U25)
ind <- intersect(rownames(GG25$G), rownames(Gh))
stopifnot(length(ind) == nrow(Gh))
Gh <- Gh[ind, , drop = FALSE]; pop <- pop[match(ind, rownames(Gc)[!par])]
H25 <- GG25$G[ind, , drop = FALSE] / 2; obs25 <- !is.na(H25)
S_all <- rowSums(H25, na.rm = TRUE); N_all <- rowSums(obs25)
X <- Gh / 2
for (p in hpops) { r <- pop == p; X[r, ] <- sweep(X[r, , drop = FALSE], 2, colMeans(X[r, , drop = FALSE], na.rm = TRUE)) }
for (ch in unique(lo$Chr)) {
  u25 <- which(U25$u$Chr == ch)
  h <- (S_all - rowSums(H25[, u25, drop = FALSE], na.rm = TRUE)) / (N_all - rowSums(obs25[, u25, drop = FALSE]))
  for (p in hpops) { r <- pop == p; h[r] <- h[r] - mean(h[r]) }
  cc <- which(lo$Chr == ch); Xc <- X[, cc, drop = FALSE]; ok <- !is.na(Xc); Xc0 <- Xc; Xc0[!ok] <- 0
  X[, cc] <- Xc - outer(h, colSums(Xc0 * h) / colSums(ok * h^2))
}
X[is.na(X)] <- 0
nw <- sqrt(colSums(X^2)); nw[nw == 0] <- NA
Zw <- sweep(X, 2, nw, "/")
rm(GG25, H25, X); invisible(gc())

## ---- pair sums per chromosome x cM bin, and unlinked ----------------------------------
chrs <- unique(lo$Chr); cols <- split(seq_len(nrow(lo)), factor(lo$Chr, levels = chrs))
pairs <- lapply(cols, function(ix) {
  m <- length(ix); ut <- which(upper.tri(diag(m)))
  i <- ix[row(diag(m))[ut]]; j <- ix[col(diag(m))[ut]]
  dcm <- abs(lo$cM_pos[j] - lo$cM_pos[i])
  cm <- cut(dcm, CM_BREAKS, labels = FALSE, right = FALSE); cm[which(dcm == 0)] <- NA
  list(ut = ut, cm = cm, gg = g[i] * g[j])
})
half <- function(tot, parts) (tot - parts) / 2
sums <- function(A, W = NULL) {
  w <- rbindlist(lapply(seq_along(cols), function(k) {
    p <- pairs[[k]]; ca <- crossprod(A[, cols[[k]], drop = FALSE])[p$ut]
    w2 <- if (is.null(W)) NA_real_ else crossprod(W[, cols[[k]], drop = FALSE])[p$ut]^2
    num <- ca^2 * p$gg; okn <- !is.na(num)
    data.table(Chr = chrs[k], bin = p$cm, num = ifelse(okn, num, 0), den = ifelse(okn, p$gg, 0),
               ceil = ifelse(is.na(p$gg), 0, p$gg), nceil = !is.na(p$gg),
               sw = ifelse(is.na(w2), 0, w2), nw = !is.na(w2))[!is.na(bin), lapply(.SD, sum), by = .(Chr, bin)]
  }))
  A0 <- A; A0[is.na(A0)] <- 0; okA <- colSums(is.na(A)) == 0
  gw <- ifelse(is.na(g) | !okA, 0, g); gn <- !is.na(g)
  M <- lapply(cols, function(ix) A0[, ix, drop = FALSE] %*% (gw[ix] * t(A0[, ix, drop = FALSE])))
  gs <- sapply(cols, function(ix) sum(gw[ix])); ng <- sapply(cols, function(ix) sum(gn[ix]))
  ## per chromosome-pair unlinked sums are needed for the bootstrap: keep the matrices
  unl <- list(M = M, gs = gs, ng = ng)
  if (!is.null(W)) {
    W0 <- W; W0[is.na(W0)] <- 0; okW <- colSums(is.na(W)) == 0
    unl$Mw <- lapply(cols, function(ix) tcrossprod(W0[, ix, drop = FALSE])); unl$nwc <- sapply(cols, function(ix) sum(okW[ix]))
  }
  list(within = w, unl = unl)
}
obs <- sums(Za, Zw)

## ---- permutation baseline for the realised fraction ------------------------------------
npop <- nrow(Za)
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
base <- parallel::mclapply(seq_len(N_PERM), mc.cores = N_CORES, mc.set.seed = TRUE, FUN = function(b) {
  P <- replicate(ncol(Za), sample.int(npop))
  Zp <- matrix(Za[cbind(as.vector(P), rep(seq_len(ncol(Za)), each = npop))], npop)
  s <- sums(Zp)
  wb <- s$within[, .(r = sum(num) / sum(den)), by = bin]
  M <- s$unl$M; gs <- s$unl$gs
  rbind(wb, data.table(bin = 99L, r = half(sum(Reduce(`+`, M)^2), sum(sapply(M, function(m) sum(m^2)))) / half(sum(gs)^2, sum(gs^2))))
})
base <- rbindlist(base)[, .(b = mean(r), b_lo = quantile(r, 0.025), b_hi = quantile(r, 0.975)), by = bin]

## ---- estimates + chromosome-block bootstrap -------------------------------------------
draws <- chrom_draws(chrs, B, SEED)
est_within <- function(w, num, den) {
  s <- tapply(w[[num]], w$Chr, sum)[chrs]; n <- tapply(w[[den]], w$Chr, sum)[chrs]; s[is.na(s)] <- 0; n[is.na(n)] <- 0
  bt <- vapply(draws, function(x) sum(x * s) / sum(x * n), numeric(1))
  c(mean = sum(s) / sum(n), lo = unname(quantile(bt, 0.025)), hi = unname(quantile(bt, 0.975)))
}
## unlinked: pair (a,b) contributes <M_a, M_b>; weight w_a w_b, a != b
pairmat <- function(L) { k <- length(L); m <- matrix(0, k, k); for (a in seq_len(k)) for (b2 in seq_len(k)) if (a != b2) m[a, b2] <- sum(L[[a]] * L[[b2]]); m }
est_unl <- function(Snum, Sden) {
  bt <- vapply(draws, function(x) { W <- outer(x, x); sum(W * Snum) / sum(W * Sden) }, numeric(1))
  c(mean = sum(Snum) / sum(Sden), lo = unname(quantile(bt, 0.025)), hi = unname(quantile(bt, 0.975)))
}
U <- obs$unl
S_num <- pairmat(U$M); S_den <- outer(U$gs, U$gs); diag(S_den) <- 0
S_ceil_n <- outer(U$ng, U$ng); diag(S_ceil_n) <- 0
S_w <- pairmat(U$Mw); S_wn <- outer(U$nwc, U$nwc); diag(S_wn) <- 0

res <- rbindlist(lapply(c(sort(unique(obs$within$bin)), 99L), function(bn) {
  if (bn == 99L) {
    rr <- est_unl(S_num, S_den); ce <- est_unl(S_den, S_ceil_n); ww <- est_unl(S_w, S_wn)
  } else {
    w <- obs$within[bin == bn]
    rr <- est_within(w, "num", "den"); ce <- est_within(w, "ceil", "nceil"); ww <- est_within(w, "sw", "nw")
  }
  data.table(bin = bn, qty = c("realised", "ceiling", "within_LD_adj"),
             mean = c(rr["mean"], ce["mean"], ww["mean"]), lo = c(rr["lo"], ce["lo"], ww["lo"]), hi = c(rr["hi"], ce["hi"], ww["hi"]))
}))
res <- merge(res, base, by = "bin", all.x = TRUE)
sh <- res[qty == "realised", .(bin, qty = "realised_share", mean = (mean - b) / (1 - b), lo = (lo - b) / (1 - b), hi = (hi - b) / (1 - b),
                               b = NA_real_, b_lo = NA_real_, b_hi = NA_real_)]
res <- rbind(res, sh)
res[, bin := factor(ifelse(bin == 99L, "unlinked", CM_LABELS[bin]), levels = c(CM_LABELS, "unlinked"))]
res[, set := sprintf("near-neutral (DI <= %d)", DI_MAX)]

d25 <- readRDS(file.path(OUT_DATA, "03_decay.rds"))$cm[qty %in% c("realised", "realised_share", "ceiling", "within_LD_adj"),
                                                        .(bin, qty, mean, lo, hi)][, set := "ancestry-informative (DI > -25)"]
cmp <- rbind(d25, res[, .(bin, qty, mean, lo, hi, set)])
cat(sprintf("\n[06] mean Nei-type G: near-neutral %.3f vs DI25 %.3f\n", mean(g, na.rm = TRUE), mean(U25$u$G_nei, na.rm = TRUE)))
cat("\n[06] DI25 vs near-neutral, by genetic distance:\n")
print(dcast(cmp[, .(bin, qty, set, v = sprintf("%.4f [%.4f,%.4f]", mean, lo, hi))], qty + bin ~ set, value.var = "v"))
saveRDS(list(units = lo, result = res, compare = cmp, baseline = base, n_perm = N_PERM, DI_MAX = DI_MAX, PMAF_MIN = PMAF_MIN),
        file.path(OUT_DATA, "06_lowDI.rds"))

## ---- figure ----------------------------------------------------------------------------
lab <- c(ceiling = "a  ceiling: max among-population LD (log)", realised_share = "b  share of the ceiling realised beyond chance",
         within_LD_adj = "c  within-population LD (log)")
pd <- cmp[qty %in% names(lab)][, qty := factor(qty, levels = names(lab), labels = lab)]
## log scale only for panels a and c
p_log <- function(q) ggplot(pd[qty == lab[[q]]], aes(bin, mean, colour = set, group = set)) +
  geom_line() + geom_point(size = 1.3) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
  scale_y_log10() + scale_colour_manual(values = c("#d95f02", "#1b9e77"), name = NULL) +
  labs(x = "genetic distance between units", y = NULL, title = lab[[q]]) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom", plot.title = element_text(size = 10))
p_lin <- ggplot(pd[qty == lab[["realised_share"]]], aes(bin, mean, colour = set, group = set)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_line() + geom_point(size = 1.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
  scale_colour_manual(values = c("#d95f02", "#1b9e77"), name = NULL) +
  labs(x = "genetic distance between units", y = NULL, title = lab[["realised_share"]]) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom", plot.title = element_text(size = 10))
pp <- patchwork::wrap_plots(p_log("ceiling"), p_lin, p_log("within_LD_adj"), nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG, "06_lowDI_contrast.png"), pp, width = 13, height = 4.8, dpi = 200)
ggsave(file.path(OUT_FIG, "06_lowDI_contrast.pdf"), pp, width = 13, height = 4.8)
cat("[06] done\n")
