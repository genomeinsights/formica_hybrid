## =========================================================================
## module_allele_specific_sorting -- 03: is sorting allele-specific or
## haplotype-wide? (A2)
##
## Anchors  : directionally sorted units (sort_class aquilonia / polyctena,
##            tau = 0.6, phi = 0.85, binom, alpha = 0.05 -- locked convention).
## Controls : unsorted units matched to each anchor (K nearest, with
##            replacement) on F_ST, log local recombination rate, parental
##            differentiation (parent_diff) and local unit density (units within
##            +-100 kb, which sets how many neighbours a unit can have).
##            Matching on F_ST matters: concordance rises with F_ST genome-wide
##            (module_population_partitioning, rho = 0.14), so an unmatched
##            comparison would confound "sorted" with "differentiated".
## Statistic: mean over each focal unit's neighbours, per distance bin (<= 2 Mb),
##            of conc, conc_resid, r_w, r_w_adj, r_w_poly (signed; oriented to
##            aquilonia, so + = same-ancestry association).
##            r_w_poly = pooled within-population LD over ONLY the populations in
##            which both loci segregate. The pooled r_w is diluted whenever one
##            locus is monomorphic in a population (it adds nothing to the
##            covariance while the other locus's variance still enters the
##            denominator); sorted units are near-fixed in most populations, so
##            r_w is biased downward for them specifically. r_w_poly is the fair
##            anchor-vs-control comparison; r_w / r_w_adj are kept for reference.
## Prediction: allele-specific sorting -> anchor decay ~ control decay.
##             Sweep / haplotype-wide sorting -> anchors show elevated
##             concordance / LD out to larger distances.
## Uncertainty: anchor - control difference, chromosome-block bootstrap.
##
## Inputs : data/01_units.rds, data/01_pairs_le2Mb.rds
## Outputs: data/03_anchor_decay.rds, Figures/03_sorted_anchor_decay.{png,pdf}
## Run from the formica_hybrid repo root.
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
K_CTRL <- 5L; B <- 2000L; SEED <- 1L
NEAR_BREAKS <- BP_BREAKS[BP_BREAKS <= 2e6]; NEAR_LABELS <- BP_LABELS[seq_len(length(NEAR_BREAKS) - 1)]
STATS <- c("conc", "conc_resid", "r_w", "r_w_adj", "r_w_poly")
STATS_PAIR <- setdiff(STATS, "r_w_poly")

u  <- readRDS(file.path(OUT_DATA, "01_units.rds"))
pr <- readRDS(file.path(OUT_DATA, "01_pairs_le2Mb.rds"))
stopifnot(nrow(u) == N_UNITS)

## both directions: each pair contributes to the neighbourhood of i and of j
nb <- rbind(pr[, c(list(focal = i, nbr = j), .SD), .SDcols = c("Chr", "dist_bp", STATS_PAIR)],
            pr[, c(list(focal = j, nbr = i), .SD), .SDcols = c("Chr", "dist_bp", STATS_PAIR)])
nb[, bin := cut(dist_bp, NEAR_BREAKS, labels = NEAR_LABELS, right = FALSE)]
nb <- nb[!is.na(bin)]
rm(pr); invisible(gc())

## local density covariate
dens <- nb[dist_bp <= 1e5, .N, by = focal]
u[, density100k := 0L][dens, on = .(idx = focal), density100k := i.N]

## ---- matching ----------------------------------------------------------------
cov_cols <- c("FST", "log_recomb", "parent_diff", "density100k")
u[, log_recomb := log10(recomb_rate)]
ok <- complete.cases(u[, ..cov_cols]) & is.finite(u$log_recomb)
anc <- u[ok & sorted]; ctl_pool <- u[ok & sort_class == "unsorted"]
sc  <- scale(rbind(anc[, ..cov_cols], ctl_pool[, ..cov_cols]))
Xa  <- sc[seq_len(nrow(anc)), , drop = FALSE]; Xc <- sc[-seq_len(nrow(anc)), , drop = FALSE]
set.seed(SEED)
match_idx <- lapply(seq_len(nrow(Xa)), function(k) {
  d <- colSums((t(Xc) - Xa[k, ])^2)
  ctl_pool$idx[order(d)[seq_len(K_CTRL)]]
})
ctl <- data.table(anchor = rep(anc$idx, each = K_CTRL), focal = unlist(match_idx))

balance <- rbind(anc[, c(list(set = "anchors (sorted)"), lapply(.SD, mean)), .SDcols = cov_cols],
                 u[ctl$focal, c(list(set = "matched controls"), lapply(.SD, mean)), .SDcols = cov_cols],
                 ctl_pool[, c(list(set = "all unsorted"), lapply(.SD, mean)), .SDcols = cov_cols])
cat("[03] covariate balance (means):\n"); print(balance, digits = 3)
cat(sprintf("[03] %d anchors (%d aquilonia, %d polyctena); %d control draws from %d unique unsorted units\n",
            nrow(anc), sum(anc$sort_class == "aquilonia"), sum(anc$sort_class == "polyctena"),
            nrow(ctl), uniqueN(ctl$focal)))

## ---- r_w_poly for the focal units' neighbourhoods (populations where both segregate) ----
nb <- nb[focal %in% c(anc$idx, ctl$focal)]
U  <- load_units(); stopifnot(identical(U$u$group_id, u$group_id))
GG <- load_oriented_genotypes(U)
Xp <- lapply(split(seq_along(GG$pop), GG$pop), function(r) {
  X <- GG$G[r, , drop = FALSE] / 2
  X <- sweep(X, 2, colMeans(X, na.rm = TRUE)); poly <- colSums(X^2, na.rm = TRUE) > 0
  X[is.na(X)] <- 0
  list(X = X, poly = poly)
})
rm(GG, U); invisible(gc())
nb[, r_w_poly := NA_real_]
chunks <- split(seq_len(nrow(nb)), ceiling(seq_len(nrow(nb)) / 2e5))
for (ch in chunks) {
  fi <- nb$focal[ch]; ni <- nb$nbr[ch]
  num <- dx <- dy <- numeric(length(ch))
  for (xp in Xp) {
    both <- xp$poly[fi] & xp$poly[ni]
    a <- xp$X[, fi, drop = FALSE]; b <- xp$X[, ni, drop = FALSE]
    num <- num + both * colSums(a * b); dx <- dx + both * colSums(a^2); dy <- dy + both * colSums(b^2)
  }
  set(nb, ch, "r_w_poly", ifelse(dx > 0 & dy > 0, num / sqrt(dx * dy), NA_real_))
}
rm(Xp); invisible(gc())

## ---- per-focal neighbourhood means ------------------------------------------
focal_means <- nb[, c(lapply(.SD, mean, na.rm = TRUE), list(n_nbr = .N)),
                  by = .(focal, Chr, bin), .SDcols = STATS]

## group tables: anchors once each; controls weighted by how often they were drawn
grp <- rbindlist(list(data.table(focal = anc$idx, set = "anchor", w = 1, dir = anc$sort_class),
                      ctl[, .(w = .N), by = focal][, `:=`(set = "control", dir = NA_character_)]),
                 use.names = TRUE)
grp_dir <- rbindlist(list(grp[set == "anchor"],
                          ctl[, .(focal, dir = u$sort_class[anchor])][, .(w = .N), by = .(focal, dir)][, set := "control"]),
                     use.names = TRUE)
fm <- merge(focal_means, grp, by = "focal", allow.cartesian = TRUE)

## ---- chromosome-block bootstrap of anchor - control ------------------------
chrs <- sort(unique(u$Chr))
set.seed(SEED)
draws <- replicate(B, tabulate(match(sample(chrs, length(chrs), replace = TRUE), chrs), length(chrs)),
                   simplify = FALSE)
summarise_diff <- function(d) {
  rbindlist(lapply(STATS, function(st) {
    d[!is.na(get(st)), {
      sa <- tapply(w * get(st) * (set == "anchor"),  Chr, sum)[chrs]; na <- tapply(w * (set == "anchor"),  Chr, sum)[chrs]
      scn <- tapply(w * get(st) * (set == "control"), Chr, sum)[chrs]; nc <- tapply(w * (set == "control"), Chr, sum)[chrs]
      sa[is.na(sa)] <- 0; na[is.na(na)] <- 0; scn[is.na(scn)] <- 0; nc[is.na(nc)] <- 0
      bt <- vapply(draws, function(x) sum(x * sa) / sum(x * na) - sum(x * scn) / sum(x * nc), numeric(1))
      list(stat = st, anchor = sum(sa) / sum(na), control = sum(scn) / sum(nc),
           diff = sum(sa) / sum(na) - sum(scn) / sum(nc),
           lo = quantile(bt, 0.025, na.rm = TRUE), hi = quantile(bt, 0.975, na.rm = TRUE))
    }, by = bin]
  }))
}
res_all <- summarise_diff(fm)[, dir := "all sorted"]
fm_dir  <- merge(focal_means, grp_dir[!is.na(dir)], by = "focal", allow.cartesian = TRUE)
res_dir <- rbindlist(lapply(c("aquilonia", "polyctena"), function(dd) summarise_diff(fm_dir[dir == dd])[, dir := dd]))
res <- rbind(res_all, res_dir)
res[, bin := factor(bin, levels = NEAR_LABELS)]

cat("\n[03] anchor - control difference by distance (all sorted anchors):\n")
print(res[dir == "all sorted", .(stat, bin, anchor = round(anchor, 4), control = round(control, 4),
                                 diff = sprintf("%.4f [%.4f,%.4f]", diff, lo, hi))][order(stat, bin)])

saveRDS(list(result = res, balance = balance, matches = ctl, K_CTRL = K_CTRL, B = B),
        file.path(OUT_DATA, "03_anchor_decay.rds"))

## ---- figure --------------------------------------------------------------------
lab <- c(conc = "concordance", conc_resid = "concordance, ancestry-residualised",
         r_w = "within-pop LD (r_w, pooled)", r_w_adj = "within-pop LD, hybrid-index adj.",
         r_w_poly = "within-pop LD, both segregating")
long <- melt(res, id.vars = c("stat", "bin", "dir", "lo", "hi", "diff"), measure.vars = c("anchor", "control"),
             variable.name = "set", value.name = "mean")
p1 <- ggplot(long[dir == "all sorted"], aes(bin, mean, colour = set, group = set)) +
  geom_line() + geom_point(size = 1.3) +
  facet_wrap(~ stat, nrow = 1, scales = "free_y", labeller = as_labeller(lab)) +
  scale_colour_manual(values = c(anchor = "#d95f02", control = "grey40"),
                      labels = c(anchor = "sorted anchors", control = "F_ST-matched unsorted")) +
  labs(x = NULL, y = "mean over neighbours", colour = NULL) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
p2 <- ggplot(res, aes(bin, diff, colour = dir, group = dir)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(position = position_dodge(0.4)) + geom_point(size = 1.3, position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, position = position_dodge(0.4)) +
  facet_wrap(~ stat, nrow = 1, scales = "free_y", labeller = as_labeller(lab)) +
  scale_colour_manual(values = c("all sorted" = "black", aquilonia = "#1b9e77", polyctena = "#7570b3")) +
  labs(x = "distance to neighbour", y = "anchor - control", colour = NULL) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
if (requireNamespace("patchwork", quietly = TRUE)) {
  p <- patchwork::wrap_plots(p1, p2, ncol = 1)
  ggsave(file.path(OUT_FIG, "03_sorted_anchor_decay.png"), p, width = 12, height = 7.5, dpi = 200)
  ggsave(file.path(OUT_FIG, "03_sorted_anchor_decay.pdf"), p, width = 12, height = 7.5)
}
cat("[03] done\n")
