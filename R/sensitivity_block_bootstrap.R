## =========================================================
## SHARED SENSITIVITY: block bootstrap of the main A/B/C coefficients
## =========================================================
## WHY. The headline A/B/C coefficients are fitted across very large numbers of
## observations, and their model-based z-statistics are correspondingly enormous
## (|z| up to ~180). Even at the LD-cluster / representative level, where marker
## pseudo-replication is already removed, residual spatial dependence remains
## (neighbouring units share ancestry tracts, local recombination and diversity).
## The model-based SE assumes the units are independent and so overstates the
## precision. A non-parametric BLOCK bootstrap -- resampling whole genomic blocks
## with replacement and refitting -- gives confidence intervals that respect that
## dependence, without changing the point estimates.
##
## This is ONE shared sensitivity procedure over the existing coefficients, not a
## new analysis. It leaves every module's point estimate untouched; it only
## re-quantifies the uncertainty.
##
## BLOCKS (both reported).
##   - chromosome            : 27 blocks. Most conservative; assumption-light about
##                             the dependence length-scale.
##   - contiguous BLOCK_CM cM : finer; more blocks, tighter CIs. 10 cM exceeds the
##                             ~10-15 cM admixture-LD decay in these populations.
##
## LEVEL. Unit level only, matching the pipeline's unit-of-analysis principle:
##   - sorting coefficients : LD-pruning representatives (one per cluster), i.e.
##                            Module B's units. The representative-level sorting
##                            statistics equal Module A's per-SNP output subset to
##                            the representatives (parallelism_stats is per-marker
##                            deterministic), so we READ them from
##                            data/moduleA_snp.rds rather than recomputing.
##   - climate coefficient  : has_eMLG clusters, i.e. Module C's units.
##
## TARGET COEFFICIENTS (each self-validated against its published estimate before
## bootstrapping -- if the full-data refit does not reproduce the module value the
## reconstruction is wrong and the script stops):
##   A/B  DI -> sorting direction   : glm(is_aqu ~ zDI + zr + zmaf + zcs)  coef zDI   [+1.46]
##   B    recomb -> direction       : same glm                            coef zr    [-0.09]
##   B    recomb -> magnitude       : lm (prop_fixed ~ zr + zDI)          coef zr    [+0.05]
##   C    climate rate -> diagnostic: lm (frac_hi ~ p_sig + log n_loci),  coef p_sig
##          n_loci >= 50, primary config (aland_excluded x withOmega), per PC       [+4.5 / +3.7 pp per +10pp]
##
## Reads : data/moduleA_snp.rds, data/hybrids_only_maf005.Rdata (map only),
##         data/Frufa_DTOL_PR.ref_genome.recmap, data/eMLG_5loci_0025_cM05.rds,
##         ./aland_excluded/PC{1,2}_DIEM_withOmega_summary_betai_reg.out
## Writes: data/sensitivity_block_bootstrap.rds,
##         Figures/sensitivity_bootstrap_forest.{pdf,png}
## Run from the repo root. DESCRIPTIVE: wider CIs qualify precision, not direction;
## the neutral-null licence (Module E) is still what any selection reading needs.

suppressMessages({ library(data.table); library(ggplot2) })
set.seed(1)

## ---- parameters ---------------------------------------------------------
B           <- 10000L               # bootstrap replicates
BLOCK_CM    <- 10                   # genetic-block width (cM) for the finer bootstrap
CI          <- c(0.025, 0.975)      # percentile interval
SIG_THR     <- 15                   # BayPass BF(dB) significance threshold (= Module C)
DI_TH       <- -25                  # high-DI definition (= Module C)
MIN_NLOCI_C <- 50L                  # climate stratum: clusters with >= 50 loci (reported)
PRIM_POPSET <- "aland_excluded"; PRIM_OMEGA <- "withOmega"
CLUSTERING  <- "data/eMLG_5loci_0025_cM05.rds"
RECMAP      <- "data/Frufa_DTOL_PR.ref_genome.recmap"
uni         <- c("aquilonia", "polyctena")

## published point estimates, for the self-validation gate (tolerances are loose;
## this only catches a broken reconstruction, not rounding).
PUBLISHED <- list(
  `A/B DI->direction`      =  1.46,
  `B recomb->direction`    = -0.09,
  `B recomb->magnitude`    =  0.05,
  `C rate->diagnostic PC1` =  4.5,
  `C rate->diagnostic PC2` =  3.7)

## ---- helpers ------------------------------------------------------------
elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))
pct_ci  <- function(x) { x <- x[is.finite(x)]; if (!length(x)) c(NA, NA) else unname(quantile(x, CI)) }

## generic block bootstrap. `row_block` = block id per row of the model frame;
## `est(idx)` refits on the given row indices and returns a named coef vector
## (NA on failure). Blocks are resampled with replacement to their own number;
## a block drawn k times contributes its rows k times.
block_boot <- function(row_block, est, B) {
  by_block <- split(seq_along(row_block), row_block)
  blocks   <- names(by_block)
  probe    <- est(unlist(by_block, use.names = FALSE))          # == full-data estimate
  out <- matrix(NA_real_, B, length(probe), dimnames = list(NULL, names(probe)))
  for (b in seq_len(B)) {
    samp <- sample(blocks, length(blocks), replace = TRUE)
    idx  <- unlist(by_block[samp], use.names = FALSE)
    out[b, ] <- tryCatch(est(idx), error = function(e) rep(NA_real_, length(probe)))
  }
  list(boot = out, n_blocks = length(blocks))
}

## ========================================================================
## 1 -- shared inputs: map, recombination, genetic position, clustering
## ========================================================================
t0 <- Sys.time()
em <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = em)
map <- as.data.table(em$map_hyb_005)[, .(marker, Chr, Pos, DiagnosticIndex)]
rm(em); invisible(gc())

rec <- fread(RECMAP); setnames(rec, 1:4, c("chr", "pos", "cM", "cMMb"))
rec[, Chr := sub("chromosome_", "Chr", chr)]
## per-marker recombination rate (cM/Mb) and genetic position (cM), map-interpolated
map[, `:=`(recomb = NA_real_, cM = NA_real_)]
for (ch in unique(map$Chr)) {
  r <- rec[Chr == ch]; if (nrow(r) < 2) next
  idx <- map[, which(Chr == ch)]
  map[idx, recomb := approx(r$pos, r$cMMb, xout = map$Pos[idx], rule = 2)$y]
  map[idx, cM     := approx(r$pos, r$cM,   xout = map$Pos[idx], rule = 2)$y]
}

groups <- readRDS(CLUSTERING)$groups
reps   <- data.table(marker = groups$representative, n_loci = groups$n_loci)
message(sprintf("[inputs] map + recmap + clustering | %.0fs", elapsed(t0)))

## ========================================================================
## 2 -- SORTING frame (representative units) from Module A per-SNP output
## ========================================================================
snp  <- readRDS("data/moduleA_snp.rds")
keep <- c("marker", "DI", "parent_maf", "prop_fixed", "sort_class", "differentiated")
sortfr <- as.data.table(snp)[, ..keep][reps, on = "marker", nomatch = 0]  # reps adds n_loci
sortfr <- map[, .(marker, Chr, Pos, recomb, cM)][sortfr, on = "marker"]   # add position

## du = differentiated representatives with defined recomb & DI (Module B's frame).
## Standardise on du, exactly as Module B (scale() uses du's mean/sd; the outcome's
## NA rows are dropped afterwards, as lm() does internally).
du <- sortfr[differentiated == TRUE & is.finite(recomb) & is.finite(DI)]
du[, `:=`(zr   = as.numeric(scale(log10(recomb + 0.1))),
          zDI  = as.numeric(scale(DI)),
          zmaf = as.numeric(scale(parent_maf)),
          zcs  = as.numeric(scale(log10(n_loci))),
          blk10 = paste(Chr, floor(cM / BLOCK_CM)))]
du_m <- du[is.finite(prop_fixed)]                                  # magnitude frame
cu   <- du[sort_class %in% uni]                                    # direction frame
cu[, is_aqu := as.integer(sort_class == "aquilonia")]
cat(sprintf("[sorting] du_m = %d (magnitude); cu = %d unidirectional (direction)\n",
            nrow(du_m), nrow(cu)))

## precompute model matrices so each refit is a bare lm.fit / glm.fit
Xmag <- model.matrix(~ zr + zDI, du_m);              ymag <- du_m$prop_fixed
Xdir <- model.matrix(~ zDI + zr + zmaf + zcs, cu);   ydir <- cu$is_aqu
est_mag <- function(idx) {                                          # idx into du_m
  b <- .lm.fit(Xmag[idx, , drop = FALSE], ymag[idx])$coefficients
  c(`B recomb->magnitude` = b[2])                                  # col 2 = zr
}
est_dir <- function(idx) {                                          # idx into cu
  fit <- glm.fit(Xdir[idx, , drop = FALSE], ydir[idx], family = binomial())
  if (!isTRUE(fit$converged)) return(c(`A/B DI->direction` = NA_real_, `B recomb->direction` = NA_real_))
  b <- unname(fit$coefficients)                                    # drop glm.fit's coef names
  c(`A/B DI->direction` = b[2], `B recomb->direction` = b[3])      # 2 = zDI, 3 = zr
}

## ========================================================================
## 3 -- CLIMATE frame (has_eMLG clusters), primary config, per PC
## ========================================================================
he  <- groups[has_eMLG == TRUE]
m2g <- he[, .(marker = unlist(members)), by = group_id]
m2g[, DI := setNames(map$DiagnosticIndex, map$marker)[marker]]
di_cl <- m2g[, .(n_hi = sum(DI > DI_TH, na.rm = TRUE), n_di = sum(!is.na(DI))), by = group_id]

import_bf <- function(pc) {
  f <- sprintf("./%s/%s_DIEM_%s_summary_betai_reg.out", PRIM_POPSET, pc, PRIM_OMEGA)
  r <- fread(f)
  stopifnot(nrow(r) == nrow(map), identical(r$MRK, seq_len(nrow(r))))
  setNames(r$`BF(dB)`, map$marker)
}
psig_by_cluster <- function(pc) {
  m2g[, BF := import_bf(pc)[marker]]
  a <- m2g[, .(n_loci = .N, n_sig = sum(BF >= SIG_THR, na.rm = TRUE)), by = group_id]
  a[, p_sig := n_sig / n_loci][, .(group_id, n_loci, p_sig)]
}
cfr <- di_cl[psig_by_cluster("PC1")[, .(group_id, n_loci, p_sig_PC1 = p_sig)], on = "group_id"]
cfr <- psig_by_cluster("PC2")[, .(group_id, p_sig_PC2 = p_sig)][cfr, on = "group_id"]
## cluster genetic position via its representative marker
reppos <- map[data.table(marker = he$representative, group_id = he$group_id),
              on = "marker"][, .(group_id, Chr, cM)]
cfr <- reppos[cfr, on = "group_id"]
cfr <- cfr[n_di > 0 & n_loci >= MIN_NLOCI_C & is.finite(cM)]
cfr[, `:=`(frac_hi = n_hi / n_di, logn = log(n_loci), blk10 = paste(Chr, floor(cM / BLOCK_CM)))]
cat(sprintf("[climate] %d clusters with n_loci >= %d (primary config)\n", nrow(cfr), MIN_NLOCI_C))

Xc1 <- cbind(1, cfr$p_sig_PC1, cfr$logn)
Xc2 <- cbind(1, cfr$p_sig_PC2, cfr$logn); yc <- cfr$frac_hi
est_clim <- function(idx) {                                         # coef reported per +10pp rate
  f1 <- .lm.fit(Xc1[idx, , drop = FALSE], yc[idx])$coefficients
  f2 <- .lm.fit(Xc2[idx, , drop = FALSE], yc[idx])$coefficients
  c(`C rate->diagnostic PC1` = f1[2] * 10, `C rate->diagnostic PC2` = f2[2] * 10)
}

## ========================================================================
## 4 -- self-validation: full-data refit must reproduce the module estimates
## ========================================================================
full_est <- c(est_mag(seq_len(nrow(du_m))), est_dir(seq_len(nrow(cu))), est_clim(seq_len(nrow(cfr))))
val <- data.table(target = names(PUBLISHED),
                  reconstructed = round(unname(full_est[names(PUBLISHED)]), 3),
                  published = unlist(PUBLISHED))
val[, ok := abs(reconstructed - published) <= pmax(0.05, 0.15 * abs(published))]
cat("\n=== self-validation: reconstructed full-data estimate vs published ===\n"); print(val)
if (!all(val$ok)) stop("reconstruction does not reproduce a published estimate -- check the frames")

## naive model-based 95% CI (for contrast with the block-bootstrap CIs)
z95 <- qnorm(0.975)
gterm <- function(m, term, mult = 1) { s <- unname(summary(m)$coefficients[term, ])
  c(est = s[1] * mult, lo = (s[1] - z95 * s[2]) * mult, hi = (s[1] + z95 * s[2]) * mult) }
m_mag <- lm(prop_fixed ~ zr + zDI, du_m)
m_dir <- glm(is_aqu ~ zDI + zr + zmaf + zcs, cu, family = binomial())
m_c1  <- lm(frac_hi ~ p_sig_PC1 + logn, cfr); m_c2 <- lm(frac_hi ~ p_sig_PC2 + logn, cfr)
naive <- rbind(
  `A/B DI->direction`      = gterm(m_dir, "zDI"),
  `B recomb->direction`    = gterm(m_dir, "zr"),
  `B recomb->magnitude`    = gterm(m_mag, "zr"),
  `C rate->diagnostic PC1` = gterm(m_c1, "p_sig_PC1", 10),
  `C rate->diagnostic PC2` = gterm(m_c2, "p_sig_PC2", 10))

## ========================================================================
## 5 -- run the block bootstrap for both block definitions
## ========================================================================
run_all <- function(block_of) {                # block_of(frame) -> block-key vector
  bm <- block_boot(block_of(du_m), est_mag,  B)
  bd <- block_boot(block_of(cu),   est_dir,  B)
  bc <- block_boot(block_of(cfr),  est_clim, B)
  list(boot = cbind(bm$boot, bd$boot, bc$boot),
       n_blocks = c(sorting_mag = bm$n_blocks, sorting_dir = bd$n_blocks, climate = bc$n_blocks))
}
message(sprintf("[boot] chromosome blocks, B = %d ...", B)); t0 <- Sys.time()
chr  <- run_all(function(d) d$Chr)
message(sprintf("[boot] %.0f cM blocks, B = %d | chromosome done %.0fs ...", BLOCK_CM, B, elapsed(t0)))
cm10 <- run_all(function(d) d$blk10)
message(sprintf("[boot] done | %.0fs", elapsed(t0)))

targets <- names(PUBLISHED)
res <- rbindlist(lapply(targets, function(tg) {
  cc <- pct_ci(chr$boot[, tg]); c10 <- pct_ci(cm10$boot[, tg])
  data.table(target = tg, estimate = round(full_est[tg], 3),
             naive_lo = round(naive[tg, "lo"], 3), naive_hi = round(naive[tg, "hi"], 3),
             chr_lo   = round(cc[1], 3),            chr_hi   = round(cc[2], 3),
             cm10_lo  = round(c10[1], 3),           cm10_hi  = round(c10[2], 3),
             chr_excl0  = !(cc[1] <= 0 & 0 <= cc[2]),
             cm10_excl0 = !(c10[1] <= 0 & 0 <= c10[2]),
             chr_nfail  = sum(!is.finite(chr$boot[, tg])),
             cm10_nfail = sum(!is.finite(cm10$boot[, tg])))
}))
cat("\n=== block-bootstrap sensitivity of the main A/B/C coefficients ===\n")
cat(sprintf("(B = %d; chromosome blocks vs %g cM blocks; CI = central %.0f%%)\n",
            B, BLOCK_CM, 100 * (CI[2] - CI[1])))
print(res)

saveRDS(list(results = res, boot_chr = chr$boot, boot_cm10 = cm10$boot,
             full_est = full_est, naive = naive,
             n_blocks = list(chr = chr$n_blocks, cm10 = cm10$n_blocks),
             params = list(B = B, BLOCK_CM = BLOCK_CM, CI = CI, SIG_THR = SIG_THR,
                           DI_TH = DI_TH, MIN_NLOCI_C = MIN_NLOCI_C,
                           config = c(PRIM_POPSET, PRIM_OMEGA))),
        "data/sensitivity_block_bootstrap.rds")

## ========================================================================
## 6 -- forest figure: model-based vs chromosome vs cM CIs, per target (free x)
## ========================================================================
lvl <- c("model-based", "chromosome", paste0(BLOCK_CM, " cM"))
fig <- melt(res, id.vars = "target", measure.vars = patterns("_lo$", "_hi$"),
            variable.name = "method", value.name = c("lo", "hi"))
fig[, method := factor(lvl[method], levels = lvl)]
fig[, est := full_est[target]]
p <- ggplot(fig, aes(est, method)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey65") +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = method), height = 0.25, linewidth = 0.6) +
  geom_point(size = 1.6) +
  facet_wrap(~ target, scales = "free_x", ncol = 1) +
  scale_colour_manual(values = setNames(c("#999999", "#0072B2", "#D55E00"), lvl), guide = "none") +
  labs(x = "coefficient (published scale)", y = NULL) +
  theme_classic(base_size = 9) +
  theme(strip.text = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 8))
dir.create("Figures", showWarnings = FALSE)
ggsave("Figures/sensitivity_bootstrap_forest.pdf", p, width = 150, height = 165, units = "mm")
ggsave("Figures/sensitivity_bootstrap_forest.png", p, width = 150, height = 165, units = "mm", dpi = 300)
cat("\nSaved: data/sensitivity_block_bootstrap.rds, Figures/sensitivity_bootstrap_forest.pdf/.png\n")
