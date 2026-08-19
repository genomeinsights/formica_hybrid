## =============================================================================
## module_di25 -- windowed local LD (ld_w_0.95) vs BDMI coverage, on its own
## =============================================================================
## Standalone analysis of the ld_w <-> BDMI relationship (distinct from the
## residual test in di25_ldw_residual_vs_sorting_bdmi.R): does raw windowed local
## LD coincide with BDMI candidate density?
##
## Central control: the SAME correlation is computed for the empirical ld_w AND
## for the neutral 1000-replicate simulated median. Both share the marker panel
## and recombination landscape, so:
##   * if BDMI regions are simply low-recombination / high-LD blocks, empirical
##     and simulated correlations will be EQUAL (a landscape coincidence);
##   * an empirical EXCESS (emp rho > sim rho) is LD at BDMI regions beyond what
##     neutral drift on the same map makes -- the epistasis/selection-relevant part.
##
## Swept over every BDMI X^2 cutoff. For the primary cutoff a per-chromosome
## circular-rotation null (phase-randomise BDMI coverage vs the fixed ld_w track,
## wrapping within chromosome) gives a p-value that respects spatial autocorrelation.
##
## Run from the formica_hybrid repo root:
##   Rscript module_di25/R/di25_ldw_vs_bdmi.R [primary_cutoff_tag]   (default 17_002)
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

args    <- commandArgs(trailingOnly = TRUE)
PRIMARY <- if (length(args) >= 1) args[1] else "17_002"
ENV     <- "module_di25/data/di25_ldw_envelope.rds"
BEDDIR  <- "data/liftoff_Frufa_DTOL_PR"
OUTPNG  <- "module_di25/Figures/di25_ldw_vs_bdmi.png"
OUTPDF  <- sub("\\.png$", ".pdf", OUTPNG)
WIN_BP  <- 500000L
N_PERM  <- 2000L
set.seed(1)

## ---- windows with empirical + simulated-median ld_w ------------------------
pl <- as.data.table(readRDS(ENV)$pl)[is.finite(emp) & is.finite(sim_med)]
pl[, `:=`(win_start = win * WIN_BP, win_end = win * WIN_BP + WIN_BP)]

## ---- per-window BDMI bp-coverage for a given bed ---------------------------
bdmi_coverage <- function(bed) {
  b <- fread(bed, header = FALSE, col.names = c("chr", "start", "end"))
  b[, chr_id := as.integer(sub("chromosome_", "", chr))]
  w <- pl[, .(chr_id, win, ws = win_start, we = win_end)]
  setkey(b, chr_id, start, end)
  ov <- foverlaps(w, b[, .(chr_id, start, end)], by.x = c("chr_id", "ws", "we"),
                  type = "any", nomatch = 0L)
  ov[, olen := pmin(we, end) - pmax(ws, start)]
  cov <- ov[, .(bdmi_cov = sum(pmax(0, olen)) / WIN_BP), by = .(chr_id, win)]
  merge(pl[, .(chr_id, win)], cov, by = c("chr_id", "win"), all.x = TRUE)[
    is.na(bdmi_cov), bdmi_cov := 0][]
}

## ---- sweep over every BDMI cutoff ------------------------------------------
beds <- list.files(BEDDIR, "^bdmi_candidates\\..*\\.bed$", full.names = TRUE)
parse_tag <- function(f) sub(".*cutoff_([0-9]+_[0-9]+)\\.liftoff.*", "\\1", basename(f))
sp <- function(x, y) suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))

sweep <- rbindlist(lapply(beds, function(f) {
  tag <- parse_tag(f); x2 <- as.numeric(paste0("0.", sub(".*_", "", tag)))
  cov <- bdmi_coverage(f)
  d   <- merge(pl, cov, by = c("chr_id", "win"))
  data.table(tag = tag, X2 = x2, n_int = nrow(fread(f, header = FALSE)),
             mean_cov = mean(d$bdmi_cov),
             rho_emp = sp(d$emp, d$bdmi_cov), rho_sim = sp(d$sim_med, d$bdmi_cov))
}))
setorder(sweep, -X2)
cat("\n=== ld_w (windowed median) vs BDMI coverage -- across X^2 cutoffs ===\n")
cat("   (rho_emp = empirical Spearman; rho_sim = neutral-sim Spearman; excess = emp-sim)\n")
print(sweep[, .(tag, X2, n_int, mean_cov = round(mean_cov, 3),
                rho_emp = round(rho_emp, 3), rho_sim = round(rho_sim, 3),
                excess = round(rho_emp - rho_sim, 3))])

## ---- primary cutoff: rotation null for the empirical correlation ------------
covP <- bdmi_coverage(file.path(BEDDIR, sprintf("bdmi_candidates.cutoff_%s.liftoff.Frufa_DTOL_PR.bed", PRIMARY)))
dP   <- merge(pl, covP, by = c("chr_id", "win")); setorder(dP, chr_id, win)
obs_emp <- sp(dP$emp, dP$bdmi_cov); obs_sim <- sp(dP$sim_med, dP$bdmi_cov)

## circularly rotate the BDMI track within each chromosome; rotate ONCE per perm
## and score both LD tracks against the SAME rotated track (so the emp-sim excess
## null holds phase fixed across the two tracks). Two-sided p-values.
chr_idx <- split(seq_len(nrow(dP)), dP$chr_id)
rot_both <- function() {
  cov_r <- dP$bdmi_cov
  for (ix in chr_idx) {
    k <- length(ix); if (k > 1) cov_r[ix] <- dP$bdmi_cov[ix][ (seq_len(k) - 1 + sample.int(k, 1)) %% k + 1 ]
  }
  c(sp(dP$emp, cov_r), sp(dP$sim_med, cov_r))
}
null <- replicate(N_PERM, rot_both())          # 2 x N_PERM
null_emp <- null[1, ]; null_sim <- null[2, ]
obs_exc  <- obs_emp - obs_sim; null_exc <- null_emp - null_sim
two_sided <- function(null, obs, centre = 0)
  (1 + sum(abs(null - centre) >= abs(obs - centre))) / (length(null) + 1)
p_emp <- two_sided(null_emp, obs_emp)
p_sim <- two_sided(null_sim, obs_sim)
p_exc <- two_sided(null_exc, obs_exc, centre = mean(null_exc))
cat(sprintf("\n=== primary cutoff %s: rotation null (%d perms, within-chr phase, two-sided) ===\n", PRIMARY, N_PERM))
cat(sprintf("  empirical rho = %+.3f   (rot p = %.4f; null mean %+.3f)\n", obs_emp, p_emp, mean(null_emp)))
cat(sprintf("  neutral  rho = %+.3f   (rot p = %.4f; null mean %+.3f)\n", obs_sim, p_sim, mean(null_sim)))
cat(sprintf("  excess (emp - sim) = %+.3f   (rot p = %.4f; null mean %+.3f)\n", obs_exc, p_exc, mean(null_exc)))

## ---- figures ---------------------------------------------------------------
pA <- ggplot(dP, aes(bdmi_cov, emp)) +
  geom_point(aes(colour = sim_med), size = 1.1, alpha = 0.85, stroke = 0) +
  geom_smooth(method = "loess", se = FALSE, span = 1, colour = "#d95f02", linewidth = 0.7) +
  scale_colour_viridis_c(name = expression(sim~ld[w]), option = "D") +
  labs(x = sprintf("BDMI coverage per window (cutoff %s)", PRIMARY),
       y = expression("empirical windowed  " * ld[w][0.95]),
       title = sprintf("empirical ld_w vs BDMI coverage  (Spearman %+.2f, rot p=%.3f)", obs_emp, p_emp)) +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())

sw <- melt(sweep, id.vars = c("tag", "X2"), measure.vars = c("rho_emp", "rho_sim"),
           variable.name = "src", value.name = "rho")
sw[, src := factor(src, c("rho_emp", "rho_sim"), c("empirical", "neutral sim"))]
pB <- ggplot(sw, aes(X2, rho, colour = src)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  scale_x_reverse() +
  scale_colour_manual(values = c("empirical" = "#d95f02", "neutral sim" = "#1b9e77"), name = NULL) +
  labs(x = expression("BDMI "*X^2*" cutoff  (left = stringent/small regions, right = permissive/large coverage)"),
       y = expression("Spearman( windowed "*ld[w]*", BDMI coverage )"),
       title = "ld_w-BDMI correlation across cutoffs: empirical vs neutral") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "top")

p <- pA + pB + plot_annotation(
  title = "Windowed local LD (ld_w_0.95) vs BDMI candidate coverage",
  subtitle = "empirical vs neutral-sim correlation controls for the shared recombination landscape; excess (emp>sim) = LD at BDMIs beyond drift")
ggsave(OUTPDF, p, width = 12, height = 5.3)
ggsave(OUTPNG, p, width = 12, height = 5.3, dpi = 150)
saveRDS(list(sweep = sweep, primary = PRIMARY, obs_emp = obs_emp, obs_sim = obs_sim,
             p_emp = p_emp, p_sim = p_sim, obs_exc = obs_exc, p_exc = p_exc),
        "module_di25/data/di25_ldw_vs_bdmi.rds")
cat("saved:", OUTPNG, "\n")
