## =========================================================
## module_di25 (high-DI analyses) -- sorting vs recombination, swept over tau
## =========================================================
## Is high-DI sorting stronger at low recombination, and does that hold across
## the sorting threshold tau? phi = 0.85 fixed; tau in {0.5,0.6,0.7,0.8}.
##
## Inference note: LD-reduced units are NOT independent -- they are much more
## independent than SNPs (residual LD is reduced, not removed). So the per-unit
## logistic model SE is anticonservative. We therefore lead with EFFECT SIZE
## (McFadden R2; predicted P(sorted) across the range) and take CIs on the
## recombination coefficient from a GENOMIC BLOCK BOOTSTRAP that resamples whole
## chromosomes (preserving within-chromosome residual dependence). The per-SNP
## level is shown only as the pseudoreplicated foil.
##
## Run from the repo root:  Rscript module_di25/R/di25_recomb_tau_sweep.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("moduleA_sorting/R/parallelism_stats.R")     # classify_sort()

TAU_GRID <- c(0.5, 0.6, 0.7, 0.8)
N_BOOT   <- 400L
N_MIN    <- 100L
SEED     <- 1
OUTDIR   <- "module_di25/data"; FIGDIR <- "module_di25/Figures"
set.seed(SEED)

## ---- recombination ------------------------------------------------------
rec <- fread("data/Frufa_DTOL_PR.ref_genome.recmap"); setnames(rec, 1:4, c("chr","p","cM","cMMb"))
rec[, Chr := sub("chromosome_", "Chr", chr)]
rc <- function(ch, pos) { o <- rep(NA_real_, length(pos))
  for (cc in unique(ch)) { r <- rec[Chr == paste0("Chr", cc)]; if (nrow(r) < 2) next
    i <- which(ch == cc); o[i] <- approx(r$p, r$cMMb, xout = pos[i], rule = 2)$y }; o }

## ---- load counts + positions (tau-independent) --------------------------
e <- as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_emlg.rds")))
g <- readRDS(file.path(OUTDIR, "di25_clustering_cM5.rds"))$groups
e[, rmk := g$representative[match(group_id, g$group_id)]]
e[, `:=`(chr = as.integer(sub("Chr","",sub(":.*","",rmk))), pos = as.integer(sub(".*:","",rmk)))]
e[, recomb := rc(chr, pos)]; e <- e[differentiated == TRUE & n_obs > 0 & is.finite(uni_score) & is.finite(recomb)]
e[, lr := log10(recomb + 0.1)]

s <- as.data.table(readRDS(file.path(OUTDIR, "di25_sorting_snp.rds")))
s[, `:=`(chr = as.integer(sub("Chr","",sub(":.*","",marker))), pos = as.integer(sub(".*:","",marker)))]
s[, recomb := rc(chr, pos)]; s <- s[differentiated == TRUE & n_obs > 0 & is.finite(uni_score) & is.finite(recomb)]

## decile breaks from the SNP recomb distribution (shared across levels & tau)
brk <- quantile(s$recomb, 0:10/10, na.rm = TRUE)
e[, rbin := cut(recomb, brk, include.lowest = TRUE, labels = FALSE)]
s[, rbin := cut(recomb, brk, include.lowest = TRUE, labels = FALSE)]
good <- e[!is.na(rbin), .N, by = rbin][N >= N_MIN, rbin]; r_lo <- min(e[rbin %in% good, recomb])
chr_idx <- split(seq_len(nrow(e[recomb >= r_lo])), e[recomb >= r_lo]$chr)   # for block bootstrap

sorted_at <- function(dt, tau) {
  cl <- classify_sort(dt$n_aqu, dt$n_pol, dt$n_obs, sort_th = tau, sort_rule = "binom", alpha = 0.05)
  as.integer(cl %in% c("aquilonia","polyctena","unresolved","ambiguous"))
}
mcf <- function(m) { m0 <- update(m, . ~ 1); as.numeric(1 - logLik(m)/logLik(m0)) }

## ---- sweep --------------------------------------------------------------
dec_all <- list(); summ <- list()
for (tau in TAU_GRID) {
  e[, y := sorted_at(e, tau)]; s[, y := sorted_at(s, tau)]
  de <- e[!is.na(rbin), .(level="eMLG", tau=tau, n=.N, med_recomb=round(median(recomb),1),
                          frac_sorted=mean(y)), by=rbin][order(rbin)]
  ds <- s[!is.na(rbin), .(level="SNP",  tau=tau, n=.N, med_recomb=round(median(recomb),1),
                          frac_sorted=mean(y)), by=rbin][order(rbin)]
  dec_all[[length(dec_all)+1]] <- rbind(de, ds)

  er <- e[recomb >= r_lo]
  m  <- glm(y ~ lr, data = er, family = binomial())
  ## genomic block bootstrap: resample whole chromosomes (residual LD is within-chr)
  chs <- names(chr_idx)
  bc <- vapply(seq_len(N_BOOT), function(b) {
    idx <- unlist(chr_idx[sample(chs, length(chs), replace = TRUE)], use.names = FALSE)
    coef(glm.fit(cbind(1, er$lr[idx]), er$y[idx], family = binomial()))[2]
  }, numeric(1))
  pr <- predict(m, data.frame(lr = log10(c(r_lo, quantile(er$recomb, 0.95)) + 0.1)), type = "response")
  summ[[length(summ)+1]] <- data.table(
    tau = tau, frac_sorted_overall = round(mean(e$y), 3),
    coef_recomb = round(coef(m)["lr"], 3),
    boot_lo = round(quantile(bc, 0.025), 3), boot_hi = round(quantile(bc, 0.975), 3),
    McFadden_R2 = round(mcf(m), 4),
    P_sorted_loRecomb = round(100*pr[1], 1), P_sorted_hiRecomb = round(100*pr[2], 1))
}
dec <- rbindlist(dec_all); sweep <- rbindlist(summ)
saveRDS(list(deciles = dec, sweep = sweep, r_lo = r_lo,
             n_units_per_decile = e[!is.na(rbin), .N, by=rbin][order(rbin)]),
        file.path(OUTDIR, "di25_recomb_tau_sweep.rds"))

cat("\n===== sorting vs recombination, tau sweep (eMLG, reliable range recomb>=", round(r_lo,1),
    "cM/Mb; block-bootstrap CI over chromosomes) =====\n", sep="")
print(sweep)
cat("\nn independent eMLG units per recomb decile (tau-independent): ",
    paste(e[!is.na(rbin), .N, by=rbin][order(rbin)]$N, collapse=", "), "\n")

## ---- figure -------------------------------------------------------------
ne <- e[!is.na(rbin), .N, by=rbin][order(rbin)]
med_lab <- e[!is.na(rbin), round(median(recomb), 1), by=rbin][order(rbin)]$V1
de <- dec[level == "eMLG"]; de[, tau := factor(tau)]
ymax <- max(de$frac_sorted); bs <- ymax*0.98 / max(ne$N)
lo_bins <- ne[N < N_MIN, rbin]
pA <- ggplot(de, aes(rbin)) +
  { if (length(lo_bins)) annotate("rect", xmin=min(lo_bins)-0.5, xmax=max(lo_bins)+0.5, ymin=-Inf, ymax=Inf, fill="grey95") } +
  geom_col(data = ne, aes(rbin, N*bs), fill = "grey86", width = 0.85, inherit.aes = FALSE) +
  geom_line(aes(y = frac_sorted, colour = tau), linewidth = 0.8) +
  geom_point(aes(y = frac_sorted, colour = tau), size = 1.6) +
  scale_colour_viridis_d(end = 0.9, name = expression(tau)) +
  scale_x_continuous(breaks = 1:10, labels = med_lab) +
  scale_y_continuous(name = "fraction sorted (eMLG units)",
                     sec.axis = sec_axis(~ ./bs, name = "n independent units (bars)")) +
  labs(x = "recombination decile (median cM/Mb)",
       title = "Sorting vs recombination across tau (high-DI eMLG units)",
       subtitle = "sorting is weakly higher at low recombination at every tau; grey = low-n deciles (not reliably estimable)") +
  theme_bw(base_size = 10)
pB <- ggplot(sweep, aes(factor(tau), coef_recomb)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_errorbar(aes(ymin = boot_lo, ymax = boot_hi), width = 0.18, colour = "#315B7D") +
  geom_point(size = 2.4, colour = "#315B7D") +
  labs(x = expression(tau), y = "recomb coef (logit / log10 cM/Mb)",
       title = "b  Effect + block-bootstrap 95% CI",
       subtitle = "negative = sorting stronger at low recomb") +
  theme_bw(base_size = 10)
ggsave(file.path(FIGDIR, "di25_recomb_tau_sweep.png"), pA + pB + plot_layout(widths = c(2, 1)),
       width = 13, height = 5.2, dpi = 200)
cat("\n[di25-recomb-sweep] wrote", file.path(FIGDIR, "di25_recomb_tau_sweep.png"), "\n")
