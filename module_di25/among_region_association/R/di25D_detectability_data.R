## =========================================================================
## DI25 among-region -- DATA-GROUNDED companion to the detectability schematic.
## =========================================================================
## (A) within-population segregation (wp_maf) collapses as loci sort (prop_fixed -> 1): the
##     empirical basis for the within-pop LD "hump" falling to zero at fixation -- and the
##     "blind spot" (near-fixed loci, untestable for in-progress DMI) is a real, large fraction.
## (B) LD DECAY vs genetic distance for SAME-chromosome candidate pairs: AMONG-population LD
##     (R_st^2) decays to the inter-chromosomal baseline by ~10 cM (admixture-LD scale, as the
##     Ohta arm validated), but WITHIN-population LD (within_pop_r^2) does NOT -- it persists
##     chromosome-wide. This is why ALL EMMAX leads are same-chromosome (>10 cM but still
##     within-chromosome long-range LD), and why the 10 cM "unlinked" rule (valid for the
##     among-pop axis) does NOT make a same-chr pair trans on the within-pop axis.
## Writes Figures/di25D_detectability_data.png.

suppressPackageStartupMessages(library(data.table))
u <- readRDS("module_di25/among_region_association/data/di25D_units.rds")
N <- readRDS("module_di25/among_region_association/data/di25D_network.rds")
D <- u$dosage; pops <- as.character(u$pops); chr_of <- u$chr_of; cm_of <- u$cm_of
m <- u$metrics; g <- u$gate
pop_idx <- split(seq_len(nrow(D)), pops)
lead_ids <- unique(c(N$meta_edges$ma, N$meta_edges$mb))
cand_ids <- intersect(colnames(D), g[bi_score >= 0.2, group_id])

## ---- (A) within-pop segregation vs sortedness ---------------------------
dA <- merge(m[, .(group_id, wp_maf)], g[differentiated == TRUE, .(group_id, prop_fixed, bi_score)], by = "group_id")
dA[, pf_bin := cut(prop_fixed, seq(0, 1, 0.1), include.lowest = TRUE)]
aggA <- dA[, .(wp = median(wp_maf, na.rm = TRUE), q1 = quantile(wp_maf, .25, na.rm = TRUE),
               q3 = quantile(wp_maf, .75, na.rm = TRUE), mid = mean(as.numeric(sub("[^,]*,([0-9.]+).*", "\\1", pf_bin))) - 0.05,
               n = .N), by = pf_bin][order(pf_bin)]
aggA[, mid := seq(0.05, 0.95, 0.1)[.I]]
blind_frac <- dA[prop_fixed > 0.8, .N] / nrow(dA)

## ---- (B) same-chromosome LD decay: within-pop vs among-pop ---------------
decay <- rbindlist(lapply(unique(chr_of[cand_ids]), function(ch) {
  ids <- cand_ids[chr_of[cand_ids] == ch]; k <- length(ids); if (k < 3) return(NULL)
  cm <- cm_of[ids]
  Warr <- simplify2array(lapply(pop_idx, function(ix)
    suppressWarnings(cor(D[ix, ids, drop = FALSE], use = "pairwise.complete.obs"))))
  wpr <- apply(Warr, c(1, 2), median, na.rm = TRUE)                      # median within-pop r
  PF  <- vapply(pop_idx, function(ix) colMeans(D[ix, ids, drop = FALSE], na.rm = TRUE), numeric(k))  # k x nPop
  rst <- suppressWarnings(cor(t(PF)))                                    # among-pop r (over pops)
  dcm <- as.matrix(dist(cm)); ut <- upper.tri(dcm)
  data.table(dcM = dcm[ut], wpr2 = wpr[ut]^2, rst2 = rst[ut]^2)
}))
decay <- decay[is.finite(dcM) & is.finite(wpr2) & is.finite(rst2)]
decay[, cmbin := cut(dcM, c(seq(0, 60, 5), Inf), include.lowest = TRUE, right = FALSE)]
aggB <- decay[, .(cm = mean(dcM), wpr2 = mean(wpr2), rst2 = mean(rst2), n = .N), by = cmbin][order(cm)][n >= 30]
base_rst <- readRDS("module_di25/among_region_association/data/di25D_ohta.rds")$baseline_rst2

png("module_di25/among_region_association/Figures/di25D_detectability_data.png",
    width = 2000, height = 900, res = 190)
par(mfrow = c(1, 2), mar = c(4.3, 4.4, 3, 1.2), mgp = c(2.4, 0.7, 0))

## Panel A
plot(aggA$mid, aggA$wp, type = "n", xlim = c(0, 1), ylim = c(0, max(aggA$q3) * 1.05),
     xlab = "sortedness  (prop_fixed: fraction of demes near-fixed)", ylab = "within-population MAF (deme segregation)",
     main = "a  Segregation collapses as loci sort (blind spot small here)")
rect(0.8, 0, 1, par("usr")[4], col = "#8888881A", border = NA)
text(0.9, max(aggA$q3) * 0.55, sprintf("blind spot\n~%.0f%% of loci\n(sorting mostly\nincomplete)", 100 * blind_frac),
     cex = 0.66, col = "grey35")
arrows(aggA$mid, aggA$q1, aggA$mid, aggA$q3, angle = 90, code = 3, length = 0.02, col = "grey70")
## candidate + lead loci positions (drawn first, medians on top)
cf <- g[group_id %in% cand_ids, prop_fixed]; cw <- m[match(g[group_id %in% cand_ids, group_id], group_id), wp_maf]
points(cf, cw, pch = 16, col = "#0072B233", cex = 0.4)
lines(aggA$mid, aggA$wp, col = "grey25", lwd = 2.2); points(aggA$mid, aggA$wp, pch = 16, col = "grey25")
lf <- g[group_id %in% lead_ids, prop_fixed]; lw <- m[match(g[group_id %in% lead_ids, group_id], group_id), wp_maf]
points(lf, lw, pch = 21, bg = "#D55E00", col = "white", cex = 1.4, lwd = 1.1)
legend("bottomleft", bty = "n", cex = 0.72, pch = c(16, 16, 21), col = c("grey25", "#0072B2", "black"),
       pt.bg = c(NA, NA, "#D55E00"), legend = c("median (all diff loci)", "bidirectional candidates", "EMMAX lead loci"))

## Panel B -- log y so the low decay AND the extreme leads both fit
le <- copy(N$meta_edges)[, `:=`(dcm = abs(cm_of[ma] - cm_of[mb]), wpr2 = within_pop_r^2)]
yl <- range(c(aggB$wpr2, aggB$rst2, le$wpr2))
matplot(aggB$cm, cbind(aggB$wpr2, aggB$rst2), type = "b", pch = 16, lty = 1, lwd = 2, log = "y", ylim = yl,
        col = c("#2166AC", "#D55E00"), xlab = "genetic distance between same-chr loci (cM)",
        ylab = expression("mean LD  ("*r^2*", log)"), main = "b  Among-pop LD flattens by 10 cM; within-pop extends further")
abline(v = 10, lty = 2, col = "grey55"); abline(h = base_rst, lty = 3, col = "#D55E00")
points(le$dcm, le$wpr2, pch = 8, col = "#2166AC", cex = 1.4, lwd = 1.8)   # the EMMAX leads = long-range outliers
text(12, base_rst * 1.9, "10 cM\n(unlinked rule)", cex = 0.62, col = "grey40", pos = 4)
text(max(aggB$cm) * 0.7, base_rst, "inter-chr baseline", cex = 0.62, col = "#D55E00", pos = 1)
text(max(le$dcm), max(le$wpr2) * 0.82, "EMMAX leads =\nlong-range within-chr", cex = 0.66, col = "#2166AC", pos = 2)
legend("right", bty = "n", cex = 0.74, lwd = c(2, 2, NA), pch = c(16, 16, 8), col = c("#2166AC", "#D55E00", "#2166AC"),
       legend = c("within-pop (within_pop_r2)", "among-pop (R_st2)", "leads (within_pop_r2)"))
dev.off()
cat(sprintf("[done] Figures/di25D_detectability_data.png | blind-spot (prop_fixed>0.8) = %.0f%% of loci; %d same-chr pairs\n",
            100 * blind_frac, nrow(decay)))
