## =========================================================================
## module_allele_specific_sorting -- 04: long-range (unlinked) associations (A3)
##
## Epistatic selection (e.g. BDMIs) does not weaken with physical distance, so
## it predicts associations between SPECIFIC unlinked loci, not a general
## elevation of short-range LD. Such pairs are a sparse subset of ~2e8
## cross-chromosome pairs, so the test is on the TAILS, not the mean (the mean
## residual long-range concordance is ~0; module_population_partitioning).
##
## Statistics, over every cross-chromosome unit pair:
##   conc_resid  among-population concordance of ancestry-residualised profiles
##               (BDMI partners resolved into compatible combinations should
##               distinguish the SAME populations -> positive tail)
##   conc_raw    same on raw profiles -- contrast only: includes the shared
##               genome-wide ancestry gradient (which a polygenic, aquilonia-
##               biased BDMI network would ALSO produce -- not separable here)
##   r_w_adj     within-population LD after residualising on individual LOCO
##               hybrid index (conspecific-ancestry association within
##               populations while loci still segregate -> positive tail)
## Orientation: + = same parental ancestry at both loci.
##
## Null (empirical, no simulations): chromosome-wise permutation. For each
## chromosome independently, permute population labels (among-population
## statistics) or individuals within populations (within-population
## statistic), applied jointly to all units on that chromosome. This keeps
## every within-chromosome structure and each unit's marginal distribution,
## and destroys only associations between chromosomes.
##
## Reported: observed vs null tail counts (positive and negative) across a
## threshold grid; positive-minus-negative asymmetry; per-unit hub counts vs
## the null maximum; enrichment of tail-pair endpoints in the colleague's BDMI
## regions (nodes only -- the beds carry no partner information) and in sorted
## units; top pairs with leave-one-population-out robustness.
##
## Inputs : module_population_partitioning/data (via 00_utils.R), data/01_units.rds
## Outputs: data/04_longrange_tail.rds, Figures/04_longrange_tail.{png,pdf}
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/04_longrange_tail.R [B_NULL] [N_CORES]
## Permutation scans run in parallel (parallel::mclapply, L'Ecuyer-CMRG streams,
## reproducible for a given SEED and N_CORES).
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
args <- commandArgs(trailingOnly = TRUE)
B_NULL <- if (length(args) >= 1) as.integer(args[1]) else 50L
N_CORES <- if (length(args) >= 2) as.integer(args[2]) else 1L
SEED <- 1L
BREAKS <- seq(-1, 1, by = 0.005)
T_STAR <- c(conc_resid = 0.80, conc_raw = 0.80, r_w_adj = 0.35)     # primary tail thresholds
T_GRID <- list(conc_resid = seq(0.5, 0.95, 0.05), conc_raw = seq(0.5, 0.95, 0.05),
               r_w_adj = seq(0.15, 0.6, 0.05))
BDMI_CUTOFF <- 13L        # documented default (module_population_partitioning Analysis 6)

U  <- load_units()
u  <- readRDS(file.path(OUT_DATA, "01_units.rds"))
stopifnot(identical(u$group_id, U$u$group_id))
GG <- load_oriented_genotypes(U)
Z <- list(conc_resid = among_pop_Z(U$Resid),
          conc_raw   = among_pop_Z(U$Fmat),
          r_w_adj    = within_pop_Z(GG, u, adj = TRUE))
pop <- GG$pop; rm(GG); invisible(gc())

annot <- list(bdmi = bdmi_membership(u, BDMI_CUTOFF), sorted = u$sorted)
cat(sprintf("[04] units in BDMI regions (cutoff %d): %d / %d; sorted: %d\n",
            BDMI_CUTOFF, sum(annot$bdmi), nrow(u), sum(annot$sorted)))

chr_order <- unique(u$Chr)                         # u is ordered by ChrNum, Pos
chr_cols  <- split(u$idx, factor(u$Chr, levels = chr_order))

## ---- one full cross-chromosome scan ------------------------------------------
scan_cross <- function(Zm, t_star, keep_pairs = FALSE) {
  nb <- length(BREAKS) - 1L
  hist <- numeric(nb); hub <- numeric(ncol(Zm))
  ann_end <- setNames(numeric(length(annot)), names(annot)); n_tail <- 0
  kept <- list()
  for (a in seq_along(chr_order)[-length(chr_order)]) {
    rows <- chr_cols[[a]]; cols <- unlist(chr_cols[(a + 1L):length(chr_order)], use.names = FALSE)
    C <- crossprod(Zm[, rows, drop = FALSE], Zm[, cols, drop = FALSE])
    C[!is.finite(C)] <- NA
    v <- C[!is.na(C)]
    hist <- hist + tabulate(findInterval(pmin(pmax(v, -1), 1), BREAKS, rightmost.closed = TRUE), nb)
    M <- !is.na(C) & C >= t_star
    rs <- rowSums(M); cs <- colSums(M)
    hub[rows] <- hub[rows] + rs; hub[cols] <- hub[cols] + cs
    n_tail <- n_tail + sum(rs)
    for (k in names(annot)) ann_end[k] <- ann_end[k] + sum(rs * annot[[k]][rows]) + sum(cs * annot[[k]][cols])
    if (keep_pairs) {
      w <- which(M | (!is.na(C) & C <= -t_star), arr.ind = TRUE)
      if (nrow(w)) kept[[length(kept) + 1L]] <- data.table(i = rows[w[, 1]], j = cols[w[, 2]], r = C[w])
    }
  }
  list(hist = hist, hub = hub, n_tail = n_tail,
       ann_frac = ann_end / (2 * n_tail), pairs = if (keep_pairs) rbindlist(kept) else NULL)
}

## tail counts from a histogram: positive (r >= T) and negative (r <= -T)
tail_counts <- function(hist, tgrid) {
  lo <- BREAKS[-length(BREAKS)]; hi <- BREAKS[-1]
  rbindlist(lapply(tgrid, function(t) data.table(thr = t, pos = sum(hist[lo >= t - 1e-9]),
                                                 neg = sum(hist[hi <= -t + 1e-9]))))
}

## ---- chromosome-wise permutations --------------------------------------------
perm_among <- function(Zm) {
  for (cc in chr_cols) Zm[, cc] <- Zm[sample.int(nrow(Zm)), cc]
  Zm
}
perm_within <- function(Zm) {
  grp <- split(seq_along(pop), pop)
  for (cc in chr_cols) {
    p <- seq_along(pop)
    for (g in grp) p[g] <- g[sample.int(length(g))]
    Zm[, cc] <- Zm[p, cc]
  }
  Zm
}

## ---- observed + null ------------------------------------------------------------
res <- list()
for (st in names(Z)) {
  cat(sprintf("\n[04] %s: observed scan\n", st))
  obs <- scan_cross(Z[[st]], T_STAR[[st]], keep_pairs = TRUE)
  permf <- if (st == "r_w_adj") perm_within else perm_among
  RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
  nulls <- parallel::mclapply(seq_len(B_NULL), mc.cores = N_CORES, mc.set.seed = TRUE, FUN = function(b) {
    s <- scan_cross(permf(Z[[st]]), T_STAR[[st]])
    if (b %% 10 == 0) cat(sprintf("[04] %s: null %d/%d\n", st, b, B_NULL))
    list(tails = tail_counts(s$hist, T_GRID[[st]]), max_hub = max(s$hub),
         ann_frac = s$ann_frac, n_tail = s$n_tail)
  })
  tob <- tail_counts(obs$hist, T_GRID[[st]])
  tnull <- rbindlist(lapply(seq_along(nulls), function(b) nulls[[b]]$tails[, rep := b]))
  tsum <- tnull[, .(pos_null_mean = mean(pos), pos_null_lo = quantile(pos, 0.025), pos_null_hi = quantile(pos, 0.975),
                    neg_null_mean = mean(neg), neg_null_lo = quantile(neg, 0.025), neg_null_hi = quantile(neg, 0.975),
                    asym_null_lo = quantile(pos - neg, 0.025), asym_null_hi = quantile(pos - neg, 0.975)), by = thr]
  tsum <- merge(tob, tsum, by = "thr")
  tsum[, `:=`(p_pos = (1 + sapply(thr, function(t) sum(tnull[thr == t]$pos >= tob[thr == t]$pos))) / (B_NULL + 1),
              p_neg = (1 + sapply(thr, function(t) sum(tnull[thr == t]$neg >= tob[thr == t]$neg))) / (B_NULL + 1))]

  ann_null <- rbindlist(lapply(nulls, function(n) as.list(n$ann_frac)))
  ann <- rbindlist(lapply(names(annot), function(k) data.table(
    annotation = k, frac_all_units = mean(annot[[k]]), frac_tail_endpoints = obs$ann_frac[[k]],
    null_mean = mean(ann_null[[k]], na.rm = TRUE), null_lo = quantile(ann_null[[k]], 0.025, na.rm = TRUE),
    null_hi = quantile(ann_null[[k]], 0.975, na.rm = TRUE))))
  max_hub_null <- vapply(nulls, `[[`, numeric(1), "max_hub")

  ## leave-one-population-out robustness of the kept (|r| >= T*) observed pairs
  kp <- obs$pairs
  if (nrow(kp)) {
    if (st == "r_w_adj") {
      loo <- sapply(unique(pop), function(p) {
        Zs <- Z[[st]][pop != p, , drop = FALSE]; Zs <- sweep(Zs, 2, sqrt(colSums(Zs^2)), "/")
        colSums(Zs[, kp$i, drop = FALSE] * Zs[, kp$j, drop = FALSE])
      })
    } else {
      M0 <- if (st == "conc_resid") U$Resid else U$Fmat
      loo <- sapply(seq_len(nrow(M0)), function(p) {
        Zs <- among_pop_Z(M0[-p, , drop = FALSE])
        colSums(Zs[, kp$i, drop = FALSE] * Zs[, kp$j, drop = FALSE])
      })
    }
    loo <- matrix(loo, nrow = nrow(kp))
    kp[, r_loo_weakest := ifelse(r > 0, apply(loo, 1, min), apply(loo, 1, max))]
    kp[, robust := abs(r_loo_weakest) >= T_STAR[[st]]]
    kp[, `:=`(unit_i = u$group_id[i], chr_i = u$Chr[i], pos_i = u$Pos[i], fst_i = u$FST[i], sort_i = u$sort_class[i], bdmi_i = annot$bdmi[i],
              unit_j = u$group_id[j], chr_j = u$Chr[j], pos_j = u$Pos[j], fst_j = u$FST[j], sort_j = u$sort_class[j], bdmi_j = annot$bdmi[j])]
  }
  cat(sprintf("[04] %s @ T*=%.2f: observed +tail %s, -tail %s; robust (LOO) +%d / -%d; max hub %d vs null max-hub 95%% %.0f\n",
              st, T_STAR[[st]], format(tob[abs(thr - T_STAR[[st]]) < 1e-9]$pos, big.mark = ","),
              format(tob[abs(thr - T_STAR[[st]]) < 1e-9]$neg, big.mark = ","),
              if (nrow(kp)) kp[r > 0 & robust == TRUE, .N] else 0L, if (nrow(kp)) kp[r < 0 & robust == TRUE, .N] else 0L,
              as.integer(max(obs$hub)), quantile(max_hub_null, 0.95)))
  print(tsum[, .(thr, pos, pos_null_mean, p_pos, neg, neg_null_mean, p_neg)], digits = 3)
  print(ann, digits = 3)

  res[[st]] <- list(tails = tsum, tails_null = tnull, annotation = ann, hub = obs$hub,
                    max_hub_null = max_hub_null, pairs = kp, hist = obs$hist, t_star = T_STAR[[st]])
}
res$meta <- list(B_NULL = B_NULL, breaks = BREAKS, bdmi_cutoff = BDMI_CUTOFF, seed = SEED)
saveRDS(res, file.path(OUT_DATA, "04_longrange_tail.rds"))

## ---- figure ------------------------------------------------------------------------
lab <- c(conc_resid = "among-pop concordance, ancestry-residualised",
         conc_raw = "among-pop concordance, raw", r_w_adj = "within-pop LD, hybrid-index adj.")
tl <- rbindlist(lapply(names(Z), function(st) {
  d <- res[[st]]$tails
  rbind(d[, .(stat = st, thr, tail = "positive (same ancestry)", obs = pos, mean = pos_null_mean, lo = pos_null_lo, hi = pos_null_hi)],
        d[, .(stat = st, thr, tail = "negative (opposite ancestry)", obs = neg, mean = neg_null_mean, lo = neg_null_lo, hi = neg_null_hi)])
}))
p <- ggplot(tl, aes(thr)) +
  geom_ribbon(aes(ymin = pmax(lo, 0.5), ymax = pmax(hi, 0.5), fill = tail), alpha = 0.25) +
  geom_line(aes(y = pmax(mean, 0.5), colour = tail), linetype = 2) +
  geom_point(aes(y = pmax(obs, 0.5), colour = tail), size = 1.4) +
  facet_wrap(~ stat, scales = "free", labeller = as_labeller(lab)) +
  scale_y_log10() +
  scale_colour_manual(values = c("#d95f02", "#1b9e77"), aesthetics = c("colour", "fill")) +
  labs(x = "threshold |r|", y = "cross-chromosome pairs beyond threshold (log)", colour = NULL, fill = NULL,
       caption = "points: observed; dashed + band: chromosome-wise permutation null (mean, 95%)") +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG, "04_longrange_tail.png"), p, width = 12, height = 4.5, dpi = 200)
ggsave(file.path(OUT_FIG, "04_longrange_tail.pdf"), p, width = 12, height = 4.5)
cat("[04] done\n")
