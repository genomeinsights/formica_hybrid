## =========================================================
## How much LD at the DI markers is within-species vs admixture?
## =========================================================
##
## Question: if we seed SLiM founders by independent binomial draws at
## empirical per-species frequencies, we reproduce ADMIXTURE LD (between-
## species Delta p) but NOT within-species haplotype LD. How big is the
## part we would lose? At the 100%-fixed-difference extreme it is zero
## (no within-species variation), so this quantifies where the real
## DI > -25 set sits on that continuum.
##
## Approach: composite (dosage-correlation) r^2 decay computed on the SAME
## DI marker set in matched panels, all subsampled to n = 15 diploids so
## the finite-sample r^2 bias (~1/2n) is identical across panels:
##   - aqu-only     : within F. aquilonia parents        (within-species LD)
##   - pol-only     : within F. polyctena parents         (within-species LD)
##   - pooled 50/50 : 8 aqu + 7 pol                        (~admixture LD at founding)
##   - hybrid       : LangholmenW                          (empirical hybrid LD, reference)
## Stratified by local recombination rate (low vs high), because within-
## species LD and the sorting signal both concentrate in low-recomb regions.

suppressMessages({library(data.table); library(ggplot2)})

set.seed(1)

## ---- data ----
e <- new.env(); load("./data/hybrids_and_parents_maf005.Rdata", envir = e)
GTs <- e$GTs_with_parents            # samples x markers, 0/1/2 dosage
sd  <- e$sample_data_with_parents
map <- copy(e$map_hyb_005)

DI_TH   <- -25
MAXDIST <- 300000                    # bp, max pair distance
MAXPART <- 60                        # cap downstream partners per marker
N_IND   <- 15                        # common panel size
MAF_MIN <- 0.05                      # polymorphism threshold within a panel

di <- map[DiagnosticIndex > DI_TH, which = TRUE]
cat(sprintf("DI > %d markers: %d\n", DI_TH, length(di)))

## ---- per-marker recombination rate (cM/Mb) from the recmap ----
rec <- fread("./data/Frufa_DTOL_PR.ref_genome.recmap")
setnames(rec, c("chr", "pos", "cM", "cMMb"))
rec[, Chr := sub("chromosome_", "Chr", chr)]
map[, recomb := NA_real_]
for (ch in unique(map$Chr)) {
  r <- rec[Chr == ch]
  if (nrow(r) < 2) next
  idx <- map[, which(Chr == ch)]
  map[idx, recomb := approx(r$pos, r$cMMb, xout = map$Pos[idx], rule = 2)$y]
}
rec_split <- median(map$recomb[di], na.rm = TRUE)
cat(sprintf("median recomb at DI set: %.2f cM/Mb (low/high split)\n", rec_split))

## ---- panel membership (row indices into GTs) ----
rows_of <- function(pops, n) {
  idx <- which(sd$Population %in% pops)
  if (length(idx) > n) idx <- sample(idx, n)
  idx
}
panels <- list(
  aqu    = rows_of("aquilonia_parent", N_IND),
  pol    = rows_of("polyctena_parent", N_IND),
  pooled = c(sample(which(sd$Population == "aquilonia_parent"), 8),
             sample(which(sd$Population == "polyctena_parent"), 7)),
  hybrid = rows_of("LangholmenW", N_IND)
)
cat("panel sizes:", paste(names(panels), sapply(panels, length), sep="="), "\n")

## ---- within-species polymorphism at the DI set (first-order answer) ----
maf_in <- function(rows) {
  p <- colMeans(GTs[rows, di, drop = FALSE], na.rm = TRUE) / 2
  pmin(p, 1 - p)
}
maf_aqu <- maf_in(panels$aqu)
maf_pol <- maf_in(panels$pol)
cat(sprintf("\nDI markers polymorphic (MAF > %.2f) within aquilonia: %.1f%%\n",
            MAF_MIN, 100 * mean(maf_aqu > MAF_MIN, na.rm = TRUE)))
cat(sprintf("DI markers polymorphic (MAF > %.2f) within polyctena: %.1f%%\n",
            MAF_MIN, 100 * mean(maf_pol > MAF_MIN, na.rm = TRUE)))
cat("within-species MAF quantiles (aqu):\n"); print(round(quantile(maf_aqu, c(.5,.75,.9,.95,.99), na.rm=TRUE),3))
cat("within-species MAF quantiles (pol):\n"); print(round(quantile(maf_pol, c(.5,.75,.9,.95,.99), na.rm=TRUE),3))

## ---- composite r^2 decay for one panel ----
dist_breaks <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 300) * 1000

ld_decay_panel <- function(rows) {
  G <- GTs[rows, di, drop = FALSE]                 # inds x DI markers
  info <- map[di, .(Chr, Pos, recomb)]
  ## drop markers monomorphic in this panel
  poly <- apply(G, 2, function(x) {
    v <- x[!is.na(x)]; length(v) >= 6 && length(unique(v)) > 1
  })
  out <- vector("list", 0)
  for (ch in unique(info$Chr)) {
    j <- which(info$Chr == ch & poly)
    if (length(j) < 2) next
    o  <- order(info$Pos[j]); j <- j[o]
    pos <- info$Pos[j]; rr <- info$recomb[j]
    Gc <- G[, j, drop = FALSE]
    for (a in seq_len(length(j) - 1)) {
      b_hi <- min(a + MAXPART, length(j))
      b <- (a + 1):b_hi
      d <- pos[b] - pos[a]
      keep <- d <= MAXDIST
      if (!any(keep)) next
      b <- b[keep]; d <- d[keep]
      r <- suppressWarnings(cor(Gc[, a], Gc[, b, drop = FALSE],
                                use = "pairwise.complete.obs"))
      out[[length(out) + 1]] <- data.table(dist = d, r2 = as.numeric(r)^2,
                                            recomb = rr[a])
    }
  }
  res <- rbindlist(out)
  res[is.finite(r2)]
}

cat("\ncomputing r^2 decay per panel ...\n")
decay <- rbindlist(lapply(names(panels), function(nm) {
  d <- ld_decay_panel(panels[[nm]])
  d[, panel := nm][]
}))

decay[, dbin := cut(dist, dist_breaks, labels = FALSE)]
decay[, dmid := (dist_breaks[dbin] + dist_breaks[dbin + 1]) / 2]
decay[, rec_stratum := fifelse(recomb <= rec_split, "low recomb", "high recomb")]

## ---- summaries ----
summ <- decay[, .(mean_r2 = mean(r2), median_r2 = median(r2), n_pairs = .N),
              by = .(panel, rec_stratum, dmid)][order(panel, rec_stratum, dmid)]

overall <- decay[, .(mean_r2 = round(mean(r2), 4), median_r2 = round(median(r2), 4),
                     n_pairs = .N), by = .(panel, rec_stratum)][order(rec_stratum, panel)]
cat("\n=== overall mean composite r^2 by panel x recomb stratum ===\n")
print(overall)

## short-range only (<= 20 kb), where within-species LD is strongest
sr <- decay[dist <= 20000, .(mean_r2 = round(mean(r2),4), n=.N),
            by=.(panel, rec_stratum)][order(rec_stratum, panel)]
cat("\n=== short-range (<=20kb) mean composite r^2 ===\n"); print(sr)

## ---- plot ----
lab <- c(aqu="aquilonia (within-sp)", pol="polyctena (within-sp)",
         pooled="pooled 50/50 (admixture)", hybrid="LangholmenW (hybrid)")
summ[, panel_lab := factor(lab[panel], levels = lab)]

p <- ggplot(summ[n_pairs >= 20],
            aes(dmid/1000, mean_r2, color = panel_lab)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.3) +
  facet_wrap(~ rec_stratum) +
  scale_x_log10() +
  scale_color_manual(values = c("#1b9e77","#7570b3","#d95f02","grey40"), name=NULL) +
  labs(x = "distance (kb, log scale)", y = expression("mean composite " * r^2),
       title = sprintf("LD decay at DI > %d markers (n=%d per panel)", DI_TH, N_IND),
       subtitle = "within-species vs admixture vs hybrid") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

dir.create("Figures", showWarnings = FALSE)
ggsave("Figures/parent_LD_diagnostic.png", p, width = 11, height = 5, dpi = 110)
cat("\nsaved Figures/parent_LD_diagnostic.png\n")
