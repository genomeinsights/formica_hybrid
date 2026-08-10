## =========================================================================
## DI25 among-region -- SCHEMATIC: the detectability window for pairwise DMIs.
## =========================================================================
## Conceptual figure (illustrative curves, not fitted) for the power argument:
## an equal-fitness DMI builds WITHIN-population LD while the two loci are still
## segregating -> EMMAX detects it; once both loci FIX within every deme the within-pop
## variation (hence the detectable LD) collapses to zero, and only the AMONG-population
## co-sorting direction remains -- which is confounded by the shared founding LD common to
## all 20 replicates. So we are powered for IN-PROGRESS DMIs and BLIND to completed ones.
## Top: within-population genotype cartoons (individuals x 2 loci) at three sorting stages,
## showing the variance -- and thus the detectable LD -- vanishing at fixation.
## Bottom: the two LD components vs sorting progress, with the detection window and the
## blind spot. Writes Figures/di25D_detectability_schematic.png.

gcol <- c("#B2182B", "#F0F0F0", "#2166AC")            # 0 pol / 1 het / 2 aqu
set.seed(1)
## a within-population sample of n individuals at ancestry freq p_aqu with coupling r
gen_stage <- function(p_aqu, r, n = 18) {
  A <- rbinom(n, 2, p_aqu)
  Bi <- rbinom(n, 2, p_aqu)
  B <- ifelse(rbinom(n, 1, r) == 1, A, Bi)            # r = prob B tracks A (coupling)
  cbind(A, B)[order(-A - B), ]
}
stages <- list(
  list(s = 0.20, ttl = "admixed",  M = gen_stage(0.50, 0.25)),
  list(s = 0.50, ttl = "sorting",  M = gen_stage(0.62, 0.85)),
  list(s = 0.88, ttl = "fixed",    M = gen_stage(0.97, 0.90)))

## schematic LD trajectories over sorting progress s in [0,1]
s  <- seq(0, 1, length.out = 300)
W  <- (1 - exp(-5 * s)) * (1 - s)^1.4; W <- W / max(W)   # within-pop LD: build-up x remaining variance = hump
A  <- s^1.9                                              # among-pop LD (R_st / D'2st): rises to fixation
thr <- 0.33                                              # detection threshold
win <- range(s[W >= thr]); s2 <- win[2]                  # detection window [s1,s2]; blind spot = s > s2

png("module_di25/among_region_association/Figures/di25D_detectability_schematic.png",
    width = 1800, height = 1500, res = 200)
layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE), heights = c(1, 1.7))

## ---- top: within-population genotype cartoons ----------------------------
par(mar = c(2.2, 1.5, 2.6, 1.5))
for (st in stages) {
  M <- st$M; n <- nrow(M)
  plot(NA, xlim = c(0, 2), ylim = c(0, n), axes = FALSE, xlab = "", ylab = "",
       main = sprintf("%s  (s=%.2f)", st$ttl, st$s), cex.main = 1.0)
  for (j in 1:2) for (i in 1:n) rect(j - 1, n - i, j, n - i + 1, col = gcol[M[i, j] + 1], border = "white", lwd = 0.4)
  text(c(0.5, 1.5), -0.06 * n, c("locus A", "locus B"), cex = 0.75, xpd = NA)
  sub <- c(admixed = "weak LD, high variation", sorting = "strong within-pop LD",
           fixed = "~monomorphic: LD gone")[st$ttl]
  mtext(sub, side = 1, line = 1.1, cex = 0.62, font = 3,
        col = if (st$ttl == "sorting") "#2166AC" else "grey35")
}

## ---- bottom: LD components vs sorting progress ---------------------------
par(mar = c(4.2, 4.4, 1.2, 1.5))
plot(NA, xlim = c(0, 1), ylim = c(0, 1.08), xlab = "ancestry sorting progress   (admixed  ->  fixed in all demes)",
     ylab = "pairwise LD signal (schematic)", xaxs = "i", yaxs = "i", cex.lab = 1.0)
## blind spot (fixed) + detection window shading
rect(s2, 0, 1, 1.08, col = "#8888880F", border = NA)
rect(win[1], 0, s2, 1.08, col = "#1B7A3712", border = NA)
lines(s, A, col = "#D55E00", lwd = 3)                    # among-pop LD
lines(s, W, col = "#2166AC", lwd = 3)                    # within-pop LD (detectable)
abline(h = thr, lty = 3, col = "grey55")
## connectors to the three cartoon stages
for (st in stages) segments(st$s, 0, st$s, 1.08, col = "grey75", lty = 2)
text(sapply(stages, `[[`, "s"), 1.04, sapply(stages, `[[`, "ttl"), cex = 0.75, col = "grey35")
## labels
text(mean(c(win[1], s2)), 0.93, "DETECTION WINDOW\nwithin-pop LD -> EMMAX", col = "#1B5E20", cex = 0.82, font = 2)
text((s2 + 1) / 2, 0.55, "BLIND SPOT\nboth loci fixed:\nno within-pop LD;\nonly among-pop signal,\nconfounded by shared\nfounding LD", col = "grey30", cex = 0.75, font = 1)
text(0.30, W[which.min(abs(s - 0.30))] + 0.07, "within-pop LD\n(EMMAX / within_pop_r)", col = "#2166AC", cex = 0.8, font = 2, adj = 0)
text(0.86, 0.82, "among-pop LD\n(R_st / D'2st)", col = "#D55E00", cex = 0.8, font = 2, adj = 1)
text(0.02, thr + 0.03, "detection threshold", col = "grey45", cex = 0.68, adj = 0)
## where the empirical result sits
points(0.45, W[which.min(abs(s - 0.45))], pch = 21, bg = "#2166AC", col = "white", cex = 1.5, lwd = 1.2)
text(0.45, W[which.min(abs(s - 0.45))] - 0.10, "our lead\nF9614-F9879\n(segregating, 11 pops)", col = "#2166AC", cex = 0.68)
text(0.985, A[length(A)] - 0.10, "data: D'2st ~ 0\n(no systematic\namong-pop signal)", col = "#D55E00", cex = 0.66, adj = 1, font = 3)
dev.off()
cat("[done] Figures/di25D_detectability_schematic.png\n")
