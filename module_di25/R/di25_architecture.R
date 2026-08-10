## =========================================================
## module_di25 (high-DI analyses) -- is sorting stronger in low-recombination regions?
## =========================================================
## Module A's architecture analysis (Table 3 / Fig S4b / Table 4), restricted to
## the high-DI (DI > -25) loci. Full-data finding to compare against: at the UNIT
## level sorting was LOWEST in low recombination (increased with recombination),
## while the per-SNP level looked inflated there from LD pseudoreplication.
##
## Recombination (cM/Mb) is interpolated per chromosome from the recmap exactly as
## Module A does (approx(..., rule = 2)). Sorting is taken at the primary tau = 0.6
## (phi = 0.85): sort_class is re-derived from the tau-independent counts.
##
## Run from the repo root:  Rscript module_di25/R/di25_architecture.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()

RECMAP <- "data/Frufa_DTOL_PR.ref_genome.recmap"
TAU    <- 0.6
OUTDIR <- "module_di25/data"; FIGDIR <- "module_di25/Figures"

## ---- recombination interpolation (cM/Mb), Module A convention ------------
rec <- fread(RECMAP); setnames(rec, 1:4, c("chr", "pos", "cM", "cMMb"))
rec[, Chr := sub("chromosome_", "Chr", chr)]
assign_recomb <- function(chr_int, pos) {
  out <- rep(NA_real_, length(pos))
  for (ch in unique(chr_int)) {
    r <- rec[Chr == paste0("Chr", ch)]; if (nrow(r) < 2) next
    idx <- which(chr_int == ch)
    out[idx] <- approx(r$pos, r$cMMb, xout = pos[idx], rule = 2)$y
  }
  out
}

## ---- load sorting counts + positions, add recomb & tau=0.6 sort class ----
add_sort <- function(dt) {
  ok <- dt$differentiated & dt$n_obs > 0 & is.finite(dt$uni_score)
  cl <- rep(NA_character_, nrow(dt))
  cl[ok] <- classify_sort(dt$n_aqu[ok], dt$n_pol[ok], dt$n_obs[ok],
                          sort_th = TAU, sort_rule = "binom", alpha = 0.05)
  dt[, `:=`(cls = cl,
            sorted = cl %in% c("aquilonia", "polyctena", "unresolved", "ambiguous"),
            is_aqu = cl == "aquilonia", directional = cl %in% c("aquilonia", "polyctena"))]
  dt
}

snp <- as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_snp.rds")))
snp[, `:=`(chr = as.integer(sub("Chr", "", sub(":.*", "", marker))),
           pos = as.integer(sub(".*:", "", marker)))]
snp[, recomb := assign_recomb(chr, pos)]
add_sort(snp)

emlg <- as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_emlg.rds")))
g <- readRDS(file.path(OUTDIR, "di25_clustering_cM5.rds"))$groups
emlg[, rep_mk := g$representative[match(group_id, g$group_id)]]
emlg[, `:=`(chr = as.integer(sub("Chr", "", sub(":.*", "", rep_mk))),
            pos = as.integer(sub(".*:", "", rep_mk)))]
emlg[, recomb := assign_recomb(chr, pos)]
add_sort(emlg)

## ---- fraction sorted by recombination decile (breaks from the SNP set) ---
base_s <- snp[differentiated == TRUE & n_obs > 0 & is.finite(recomb)]
base_e <- emlg[differentiated == TRUE & n_obs > 0 & is.finite(recomb)]
brk <- quantile(base_s$recomb, 0:10 / 10, na.rm = TRUE)
wilson <- function(k, n) {                          # Wilson 95% CI for a proportion
  z <- 1.959964; p <- k / n; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  cbind(lo = pmax(0, ctr - hw), hi = pmin(1, ctr + hw))
}
decile_tab <- function(x, lab) {
  x[, rbin := cut(recomb, brk, include.lowest = TRUE, labels = FALSE)]
  t <- x[!is.na(rbin), .(level = lab, med_recomb = round(median(recomb), 2), n = .N,
                         n_sorted = sum(sorted), frac_sorted = mean(sorted),
                         frac_aqu = sum(is_aqu) / sum(directional)), by = rbin][order(rbin)]
  ci <- wilson(t$n_sorted, t$n); t[, `:=`(lo = ci[, 1], hi = ci[, 2])][]
}
dec <- rbind(decile_tab(base_s, "SNP"), decile_tab(base_e, "eMLG"))
saveRDS(dec, file.path(OUTDIR, "di25_architecture_deciles.rds"))

## ---- STATISTICAL TEST: does P(sorted) depend on recombination? ----------
## eMLG units are ~independent (that is what the LD reduction buys), so one
## logistic observation per unit is a valid, non-pseudoreplicated test.
eu <- base_e[is.finite(recomb)][, `:=`(y = as.integer(sorted), lr = log10(recomb + 0.1))]
N_MIN   <- 100L                                     # "reasonable n" per decile
good_bins <- dec[level == "eMLG" & n >= N_MIN, rbin]
r_lo <- min(base_e[rbin %in% good_bins, recomb])    # recomb floor of the reliable range
glm_all  <- glm(y ~ lr,             data = eu,                family = binomial())
glm_allDI<- glm(y ~ lr + scale(DI),  data = eu[is.finite(DI)], family = binomial())
glm_rr   <- glm(y ~ lr,             data = eu[recomb >= r_lo], family = binomial())
mcf <- function(m) { m0 <- update(m, . ~ 1); round(1 - as.numeric(logLik(m) / logLik(m0)), 4) }
sm  <- function(m) { c <- summary(m)$coef["lr", ]; sprintf("%+.3f (p=%.1e)", c[1], c[4]) }
## effect size in probability terms: predicted P(sorted) across the reliable range
pr_at <- function(x) predict(glm_rr, data.frame(lr = log10(x + 0.1)), type = "response")
p_lo <- pr_at(r_lo); p_hi <- pr_at(quantile(eu[recomb >= r_lo]$recomb, 0.95))
test_txt <- sprintf("P(sorted eMLG)~log10(recomb): all %s (R2=%.3f) ; reliable(recomb>=%.1f,n=%d) %s (R2=%.3f) | P(sorted) %.1f->%.1f%% across the reliable range",
                    sm(glm_all), mcf(glm_all), r_lo, nrow(eu[recomb >= r_lo]), sm(glm_rr), mcf(glm_rr),
                    100 * p_lo, 100 * p_hi)

## ---- magnitude & direction regressions (unit level, Module A form) -------
du <- base_e[is.finite(DI)]
du[, `:=`(zr = as.numeric(scale(log10(recomb + 0.1))), zDI = as.numeric(scale(DI)))]
lm_mag <- lm(prop_fixed ~ zr + zDI, data = du)

cu <- base_e[directional == TRUE & is.finite(DI) & is.finite(parent_maf)]
cu[, `:=`(zr = as.numeric(scale(log10(recomb + 0.1))), zDI = as.numeric(scale(DI)),
          zmaf = as.numeric(scale(parent_maf)), zcs = as.numeric(scale(log10(n_loci))))]
glm_dir <- glm(is_aqu ~ zDI + zr + zmaf + zcs, data = cu, family = binomial())

## also the SNP-level magnitude slope, for the pseudoreplication contrast
ds <- base_s[is.finite(DI)]
ds[, `:=`(zr = as.numeric(scale(log10(recomb + 0.1))), zDI = as.numeric(scale(DI)))]
lm_mag_snp <- lm(prop_fixed ~ zr + zDI, data = ds)

sp <- function(a, b) round(cor(a, b, method = "spearman", use = "pairwise.complete.obs"), 3)

## ---- report -------------------------------------------------------------
cat("\n===== DI25: sorting vs recombination (tau = 0.6, phi = 0.85) =====\n")
cat(sprintf("Spearman DI vs recomb (SNP): %s  |  DI vs recomb (eMLG): %s\n",
            sp(base_s$DI, base_s$recomb), sp(base_e$DI, base_e$recomb)))
cat("\nFraction sorted by recombination decile (1 = lowest recomb):\n")
print(dcast(dec, rbin + med_recomb ~ level, value.var = "frac_sorted")[order(rbin)])
cat("\nFraction toward F. aquilonia by recombination decile:\n")
print(dcast(dec, rbin ~ level, value.var = "frac_aqu")[order(rbin)])
cat(sprintf("\nMagnitude  prop_fixed ~ log10(recomb)+DI :\n  eMLG  recomb %+.3f (t=%.1f), DI %+.3f (t=%.1f)\n  SNP   recomb %+.3f (t=%.1f), DI %+.3f (t=%.1f)\n",
            coef(lm_mag)["zr"], summary(lm_mag)$coef["zr","t value"],
            coef(lm_mag)["zDI"], summary(lm_mag)$coef["zDI","t value"],
            coef(lm_mag_snp)["zr"], summary(lm_mag_snp)$coef["zr","t value"],
            coef(lm_mag_snp)["zDI"], summary(lm_mag_snp)$coef["zDI","t value"]))
cat(sprintf("\nDirection  P(aquilonia) ~ DI+recomb+MAF+logNloci (eMLG, directional):\n  DI %+.3f (z=%.1f) ; recomb %+.3f (z=%.1f)\n",
            coef(glm_dir)["zDI"], summary(glm_dir)$coef["zDI","z value"],
            coef(glm_dir)["zr"], summary(glm_dir)$coef["zr","z value"]))
cat("\nSTATISTICAL TEST -- ", test_txt, "\n", sep = "")
cat("n independent eMLG units per recomb decile: ",
    paste(dec[level == "eMLG", n], collapse = ", "), "\n")

saveRDS(list(deciles = dec, lm_mag = coef(summary(lm_mag)),
             lm_mag_snp = coef(summary(lm_mag_snp)), glm_dir = coef(summary(glm_dir)),
             glm_sorted_all = coef(summary(glm_all)), glm_sorted_allDI = coef(summary(glm_allDI)),
             glm_sorted_reliable = coef(summary(glm_rr)), n_min = N_MIN, recomb_floor = r_lo,
             di_recomb_rho = c(SNP = sp(base_s$DI, base_s$recomb), eMLG = sp(base_e$DI, base_e$recomb))),
        file.path(OUTDIR, "di25_architecture.rds"))

## ---- figure: eMLG fraction sorted vs recombination decile, with n bars ---
## x = recombination decile (median cM/Mb labelled). Background grey bars = the
## number of INDEPENDENT eMLG units in each decile (scaled to the right axis) --
## they collapse to a handful in low recombination. eMLG points carry Wilson 95%
## CIs and are sized by n; SNP is a faint reference. Deciles below the reliable-n
## floor are shaded.
e <- dec[level == "eMLG"]; s <- dec[level == "SNP"]
ymax <- max(e$hi, s$frac_sorted, na.rm = TRUE) * 1.02
bs   <- ymax / max(e$n)                              # bar scale to the fraction axis
lo_bins <- e[n < N_MIN, rbin]
pA <- ggplot(e, aes(rbin)) +
  { if (length(lo_bins)) annotate("rect", xmin = min(lo_bins) - 0.5, xmax = max(lo_bins) + 0.5,
                                  ymin = -Inf, ymax = Inf, fill = "grey95") } +
  geom_col(aes(y = n * bs), fill = "grey85", width = 0.85) +
  geom_line(data = s, aes(y = frac_sorted), colour = "grey65", linewidth = 0.6) +
  geom_point(data = s, aes(y = frac_sorted), colour = "grey65", size = 1.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, colour = "#315B7D") +
  geom_line(aes(y = frac_sorted), colour = "#315B7D", linewidth = 0.8) +
  geom_point(aes(y = frac_sorted, size = n), colour = "#315B7D") +
  scale_size_area(max_size = 5, guide = "none") +
  scale_x_continuous(breaks = 1:10, labels = round(e$med_recomb, 1)) +
  scale_y_continuous(name = "fraction sorted (eMLG; SNP faint)",
                     sec.axis = sec_axis(~ . / bs, name = "n independent eMLG units (bars)")) +
  labs(x = "recombination decile (median cM/Mb)",
       title = "Sorting is weakly stronger at lower recombination (high-DI loci)",
       subtitle = sprintf("real & highly significant but SMALL: reliable range %s, McFadden R2=%.3f; P(sorted) %.0f%%->%.0f%% across the range",
                          sm(glm_rr), mcf(glm_rr), 100 * p_lo, 100 * p_hi)) +
  theme_bw(base_size = 10)
ggsave(file.path(FIGDIR, "di25_sorting_vs_recomb.png"), pA, width = 9, height = 5.2, dpi = 200)
cat("\n[di25-architecture] wrote", file.path(FIGDIR, "di25_sorting_vs_recomb.png"), "\n")
