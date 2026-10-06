## =========================================================================
## module_allele_specific_sorting -- 04: do sorted loci share their sorting
## with nearby markers more than equally differentiated unsorted loci?
##
## Anchors : directionally sorted units (aquilonia / polyctena; tau = 0.6,
##           phi = 0.85, binom, alpha = 0.05 -- locked convention).
## Controls: up to K = 5 unsorted units per anchor, nearest on standardised
##           F_ST, log10 local recombination rate, parental differentiation and
##           local unit density (units within +-100 kb), within a caliper of
##           1 SD on EVERY covariate. Anchors without any control inside the
##           caliper are dropped and summarised separately. Caliper choice
##           (checked 0.25 / 0.5 / 1 / none): tighter calipers dropped 37% / 10%
##           of anchors -- disproportionately the large sorted LD blocks -- and
##           gave WORSE balance (|SMD| up to 0.19 / 0.11); 1 SD keeps 98% with
##           all |SMD| <= 0.04. Matching on F_ST matters: shared partitions
##           increase with F_ST genome-wide.
## Statistics, per focal locus and distance bin (<= 2 Mb), from the focal
## locus's representative SNP to each neighbouring marker:
##   similarity  correlation of the two oriented 20-population ancestry profiles
##   r_w_poly    within-population LD over the populations where BOTH markers
##               segregate (pooled r_w is diluted for sorted loci, which are
##               monomorphic in most populations), after regressing out each
##               individual's leave-one-chromosome-out hybrid index (as r_w_adj in 01)
## Two marker sets (audit item 1):
##   LD-reduced  neighbours = representative SNPs of OTHER units only
##   unpruned    neighbours = ALL 51,612 DI25 SNPs within 2 Mb, including SNPs in
##               the focal locus's own LD cluster -- where an extended haplotype
##               carrying the sorted allele would be seen
## Plus the within-unit haplotype extent (SNPs per unit, physical span) as a
## matched outcome.
## Inference (audit item 2): per-anchor matched-set difference (anchor minus
## the mean of its own controls), averaged over anchors; chromosome-block
## bootstrap resampling ANCHOR chromosomes with each matched set kept intact.
##
## Inputs : data/01_units.rds, module_di25/data/di25_inputs.rds,
##          data/hybrids_and_parents_maf005.Rdata, DI25 rho05 clustering (spans)
## Outputs: data/04_neighbourhood.rds, Figures/04_neighbourhood.{png,pdf},
##          Figures/04_balance.{png,pdf}
## Run from the formica_hybrid repo root.
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
K_CTRL <- 5L; CALIPER <- 1; MAXD <- 2e6; B <- 2000L; SEED <- 1L
NEAR_BREAKS <- c(0, 1e3, 5e3, 2e4, 1e5, 5e5, 2e6)
NEAR_LABELS <- c("0-1kb", "1-5kb", "5-20kb", "20-100kb", "100-500kb", "0.5-2Mb")
COVS <- c("FST", "log_recomb", "parent_diff", "density100k")

u <- readRDS(file.path(OUT_DATA, "01_units.rds"))
stopifnot(nrow(u) == N_UNITS)

## ---- SNP-level oriented genotypes and population profiles (all 51,612) ---------
inp <- readRDS("module_di25/data/di25_inputs.rds")
e2  <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd  <- e2$sample_data_with_parents
G_all <- rbind(inp$GTs_hyb, inp$GTs_par)
G_all <- G_all[rownames(G_all) %in% sd$Sample_ID, , drop = FALSE]
pops  <- sd$Population[match(rownames(G_all), sd$Sample_ID)]
sgn <- sign(colMeans(G_all[pops == "aquilonia_parent", , drop = FALSE], na.rm = TRUE) -
            colMeans(G_all[pops == "polyctena_parent", , drop = FALSE], na.rm = TRUE))
hyb <- !grepl("_parent$", pops)
G <- G_all[hyb, , drop = FALSE]; pop <- pops[hyb]; rm(G_all); invisible(gc())
G[, sgn < 0] <- 2L - G[, sgn < 0]
G[, is.na(sgn) | sgn == 0] <- NA
hpops <- sort(unique(pop))
P <- t(sapply(hpops, function(p) colMeans(G[pop == p, , drop = FALSE], na.rm = TRUE) / 2))
Z <- among_pop_Z(P)
map <- as.data.table(inp$map)[, `:=`(col = .I, Pos = as.integer(Pos))]
stopifnot(identical(map$marker, colnames(G)))

## guard: unit representative SNP profiles reproduce the unit Fmat used elsewhere
Fm <- readRDS(file.path("module_population_partitioning/data", "pp_units_Fmat.rds"))$Fmat
uc <- match(u$unit_marker, map$marker); stopifnot(!anyNA(uc))
dd <- abs(P[rownames(Fm), uc] - Fm); dd <- dd[is.finite(dd)]
stopifnot("SNP-level orientation does not reproduce the unit Fmat" = max(dd) < 1e-8)

## within-population genotypes for r_w_poly, aligned with 00_utils.R::within_pop_Z(adj = TRUE):
## (1) x = oriented dosage / 2, centred within each population;
## (2) for SNPs on chromosome c, regress x (pooled over populations, no intercept) on the
##     individual's leave-one-chromosome-out hybrid index h_c -- mean oriented dosage over
##     the UNITS' representative SNPs on all other chromosomes (same construction as 01),
##     centred within populations -- and keep the residual;
## (3) missing -> 0. Segregation (the "poly" restriction) is judged on the UNADJUSTED
##     centred genotypes, so it still means "polymorphic in that population".
X <- G / 2
for (p in hpops) { r <- pop == p; X[r, ] <- sweep(X[r, , drop = FALSE], 2, colMeans(X[r, , drop = FALSE], na.rm = TRUE)) }
seg_by_pop <- lapply(hpops, function(p) colSums(X[pop == p, , drop = FALSE]^2, na.rm = TRUE) > 0)
Gu <- G[, uc, drop = FALSE] / 2; obs_u <- !is.na(Gu)
S_all <- rowSums(Gu, na.rm = TRUE); N_all <- rowSums(obs_u)
for (ch in unique(map$Chr)) {
  uchr <- which(u$Chr == ch)
  h <- (S_all - rowSums(Gu[, uchr, drop = FALSE], na.rm = TRUE)) / (N_all - rowSums(obs_u[, uchr, drop = FALSE]))
  for (p in hpops) { r <- pop == p; h[r] <- h[r] - mean(h[r]) }
  cols <- map[Chr == ch, col]
  Xc <- X[, cols, drop = FALSE]; ok <- !is.na(Xc); Xc0 <- Xc; Xc0[!ok] <- 0
  bb <- colSums(Xc0 * h) / colSums(ok * h^2)
  X[, cols] <- Xc - outer(h, bb)
}
X[is.na(X)] <- 0
Xp <- lapply(seq_along(hpops), function(k) list(X = X[pop == hpops[k], , drop = FALSE], seg = seg_by_pop[[k]]))
rm(G, Gu, X); invisible(gc())

## ---- covariates and matching ------------------------------------------------------
u[, col := uc]
u[, density100k := vapply(seq_len(.N), function(k) sum(abs(Pos[Chr == Chr[k]] - Pos[k]) <= 1e5) - 1L, integer(1))]
u[, log_recomb := log10(recomb_rate)]
g <- readRDS("module_di25_rho05/data/di25_clustering_cM5_rho05.rds")$groups
span <- g[, .(group_id, span_kb = vapply(members, function(m) { p <- as.numeric(sub(".*:", "", m)); (max(p) - min(p)) / 1e3 }, numeric(1)))]
u[span, on = "group_id", span_kb := i.span_kb]
u[, log2_size := log2(n_loci_g)]

ok <- complete.cases(u[, ..COVS]) & is.finite(u$log_recomb)
anc_all <- u[ok & sorted]; pool <- u[ok & sort_class == "unsorted"]
mu <- colMeans(rbind(anc_all[, ..COVS], pool[, ..COVS])); sdv <- apply(rbind(anc_all[, ..COVS], pool[, ..COVS]), 2, sd)
Sa <- sweep(sweep(as.matrix(anc_all[, ..COVS]), 2, mu), 2, sdv, "/")
Sp <- sweep(sweep(as.matrix(pool[, ..COVS]), 2, mu), 2, sdv, "/")
matches <- rbindlist(lapply(seq_len(nrow(Sa)), function(k) {
  dev <- abs(sweep(Sp, 2, Sa[k, ]))
  inside <- which(rowSums(dev <= CALIPER) == length(COVS))
  if (!length(inside)) return(NULL)
  d2 <- rowSums(dev[inside, , drop = FALSE]^2)
  data.table(anchor = anc_all$idx[k], control = pool$idx[inside[order(d2)[seq_len(min(K_CTRL, length(inside)))]]])
}))
n_ctrl <- matches[, .N, by = anchor]
cat(sprintf("[04] %d sorted units with covariates; %d matched (>=1 control within %.2f SD caliper), %d dropped; controls per anchor: %s\n",
            nrow(anc_all), nrow(n_ctrl), CALIPER, nrow(anc_all) - nrow(n_ctrl),
            paste(names(table(n_ctrl$N)), table(n_ctrl$N), sep = ":", collapse = " ")))

## balance: standardised mean differences (anchors vs all unsorted, vs matched controls)
smd <- function(a, b, w = NULL) {
  if (is.null(w)) w <- rep(1, length(b))
  mb <- sum(w * b) / sum(w)
  (mean(a) - mb) / sqrt((var(a) + var(rep(b, times = w))) / 2)
}
anc <- u[idx %in% n_ctrl$anchor]
unmatched <- anc_all[!idx %in% n_ctrl$anchor,
                     .(group_id, Chr, Pos, sort_class, FST, recomb_rate, n_loci_g, span_kb)]
cat(sprintf("[04] unmatched anchors (no unsorted unit within the caliper): %d; mean %.0f SNPs, mean span %.0f kb (matched anchors: %.1f SNPs, %.1f kb)\n",
            nrow(unmatched), mean(unmatched$n_loci_g), mean(unmatched$span_kb),
            mean(anc$n_loci_g), mean(anc$span_kb)))
print(unmatched[order(-span_kb)][1:min(10, .N)], digits = 3)
cw  <- matches[, .(w = .N), by = control]
bal_vars <- c(COVS, "log2_size", "span_kb")
balance <- rbindlist(lapply(bal_vars, function(v) data.table(
  variable = v, role = if (v %in% COVS) "matched covariate" else "haplotype extent (outcome)",
  mean_anchor = mean(anc[[v]]), mean_all_unsorted = mean(pool[[v]]),
  mean_controls = weighted.mean(u[[v]][cw$control], cw$w),
  smd_before = smd(anc[[v]], pool[[v]]),
  smd_after = smd(anc[[v]], u[[v]][cw$control], cw$w))))
cat("[04] balance (standardised mean differences):\n"); print(balance, digits = 3)

## ---- neighbourhood statistics for every focal locus ------------------------------------
focal <- u[idx %in% c(anc$idx, matches$control)]
snp_by_chr <- split(map[, .(Chr, Pos, col)], map$Chr)
unit_cols  <- split(u[, .(Chr, Pos, col)], u$Chr)
neighbours <- function(set) rbindlist(lapply(seq_len(nrow(focal)), function(k) {
  f <- focal[k]; cand <- set[[f$Chr]]
  d <- abs(cand$Pos - f$Pos); keep <- d <= MAXD & cand$col != f$col
  if (!any(keep)) return(NULL)
  data.table(focal = f$idx, fcol = f$col, ncol = cand$col[keep], dist = d[keep])
}))
pair_stats <- function(nb) {
  nb[, bin := cut(dist, NEAR_BREAKS, labels = NEAR_LABELS, right = TRUE, include.lowest = TRUE)]
  nb[, similarity := colSums(Z[, fcol, drop = FALSE] * Z[, ncol, drop = FALSE])]
  nb[, r_w_poly := NA_real_]
  for (ch in split(seq_len(nrow(nb)), ceiling(seq_len(nrow(nb)) / 2e5))) {
    fi <- nb$fcol[ch]; ni <- nb$ncol[ch]; num <- dx <- dy <- numeric(length(ch))
    for (xp in Xp) {
      both <- xp$seg[fi] & xp$seg[ni]
      a <- xp$X[, fi, drop = FALSE]; b <- xp$X[, ni, drop = FALSE]
      num <- num + both * colSums(a * b); dx <- dx + both * colSums(a^2); dy <- dy + both * colSums(b^2)
    }
    set(nb, ch, "r_w_poly", ifelse(dx > 0 & dy > 0, num / sqrt(dx * dy), NA_real_))
  }
  nb[, .(similarity = mean(similarity, na.rm = TRUE), r_w_poly = mean(r_w_poly, na.rm = TRUE), n_nbr = .N),
     by = .(focal, bin)]
}
fm <- rbind(pair_stats(neighbours(unit_cols))[, marker_set := "LD-reduced units"],
            pair_stats(neighbours(snp_by_chr))[, marker_set := "unpruned SNPs"])
## within-unit haplotype extent as an additional "bin"-free outcome
ext <- melt(focal[, .(focal = idx, log2_size, span_kb)], id.vars = "focal", variable.name = "stat")
cat(sprintf("[04] focal loci: %d anchors + %d distinct controls\n", nrow(anc), uniqueN(matches$control)))

## ---- matched-set differences + anchor-chromosome bootstrap -------------------------------
chrs <- sort(unique(u$Chr))
draws <- chrom_draws(chrs, B, SEED)
anchor_chr <- setNames(match(u$Chr, chrs), u$idx)
matched_diff <- function(vals, by_cols) {
  ## vals: focal, <by_cols>, value ; returns per-anchor differences
  a <- merge(matches, vals, by.x = "control", by.y = "focal", allow.cartesian = TRUE)
  cmean <- a[!is.na(value), .(control_mean = mean(value)), by = c("anchor", by_cols)]
  av <- vals[focal %in% anc$idx]; setnames(av, "focal", "anchor")
  m <- merge(av[!is.na(value)], cmean, by = c("anchor", by_cols))
  m[, diff := value - control_mean][, chr_k := anchor_chr[as.character(anchor)]]
  m
}
summarise_diff <- function(m, by_cols) {
  m[, {
    w0 <- tabulate(chr_k, length(chrs))
    sd_ <- tapply(diff, factor(chr_k, levels = seq_along(chrs)), sum); sd_[is.na(sd_)] <- 0
    bt <- vapply(draws, function(w) sum(w * sd_) / sum(w * w0), numeric(1))
    list(anchor = mean(value), control = mean(control_mean), diff = mean(diff),
         lo = quantile(bt, 0.025, na.rm = TRUE), hi = quantile(bt, 0.975, na.rm = TRUE), n_anchors = .N)
  }, by = by_cols]
}
long <- melt(fm, id.vars = c("focal", "bin", "marker_set", "n_nbr"), measure.vars = c("similarity", "r_w_poly"),
             variable.name = "stat", value.name = "value")
res <- summarise_diff(matched_diff(long[, .(focal, marker_set, bin, stat, value)], c("marker_set", "bin", "stat")),
                      c("marker_set", "bin", "stat"))
res[, bin := factor(bin, levels = NEAR_LABELS)]; setorder(res, marker_set, stat, bin)
res_ext <- summarise_diff(matched_diff(ext[, .(focal, stat, value)], "stat"), "stat")

cat("\n[04] matched-set differences (anchor - mean of its controls), anchor-chromosome bootstrap:\n")
print(res[, .(marker_set, stat, bin, anchor = round(anchor, 4), control = round(control, 4),
              diff = sprintf("%.4f [%.4f,%.4f]", diff, lo, hi), n_anchors)])
cat("\n[04] within-unit haplotype extent (matched):\n")
print(res_ext[, .(stat, anchor = round(anchor, 3), control = round(control, 3), diff = sprintf("%.3f [%.3f,%.3f]", diff, lo, hi))])

saveRDS(list(result = res, extent = res_ext, balance = balance, matches = matches, unmatched = unmatched,
             n_sorted = nrow(anc_all), n_matched = nrow(anc), K_CTRL = K_CTRL, CALIPER = CALIPER, B = B),
        file.path(OUT_DATA, "04_neighbourhood.rds"))

## ---- figures --------------------------------------------------------------------------------
lab <- c(similarity = "ancestry-profile similarity", r_w_poly = "within-population LD (both segregating)")
p_abs <- ggplot(melt(res, id.vars = c("marker_set", "bin", "stat"), measure.vars = c("anchor", "control"),
                     variable.name = "set"), aes(bin, value, colour = set, group = set)) +
  geom_line() + geom_point(size = 1.3) +
  facet_grid(stat ~ marker_set, scales = "free_y", labeller = labeller(stat = lab)) +
  scale_colour_manual(values = c(anchor = "#d95f02", control = "grey40"),
                      labels = c(anchor = "sorted loci", control = "matched unsorted loci"), name = NULL) +
  labs(x = NULL, y = "mean over neighbours") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
p_diff <- ggplot(res, aes(bin, diff, group = 1)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line() + geom_point(size = 1.3) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  facet_grid(stat ~ marker_set, scales = "free_y", labeller = labeller(stat = lab)) +
  labs(x = "distance to neighbouring marker", y = "sorted - matched (95% CI)") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
p <- patchwork::wrap_plots(p_abs + ggtitle("a  neighbourhood means"), p_diff + ggtitle("b  matched-set difference"), nrow = 1)
ggsave(file.path(OUT_FIG, "04_neighbourhood.png"), p, width = 12, height = 6.5, dpi = 200)
ggsave(file.path(OUT_FIG, "04_neighbourhood.pdf"), p, width = 12, height = 6.5)

bl <- melt(balance, id.vars = c("variable", "role"), measure.vars = c("smd_before", "smd_after"), variable.name = "when")
bl[, when := factor(when, levels = c("smd_before", "smd_after"), labels = c("vs all unsorted", "vs matched controls"))]
pb <- ggplot(bl, aes(value, variable, colour = when)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_vline(xintercept = c(-0.1, 0.1), colour = "grey70", linetype = 2) +
  geom_point(size = 2) + facet_grid(role ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c("grey50", "#d95f02"), name = NULL) +
  labs(x = "standardised mean difference (sorted - unsorted)", y = NULL) + theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG, "04_balance.png"), pb, width = 7, height = 4, dpi = 200)
ggsave(file.path(OUT_FIG, "04_balance.pdf"), pb, width = 7, height = 4)
cat("[04] done\n")
