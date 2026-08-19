## =========================================================
## MODULE A (clean _PK) -- ancestry-sorting landscape along the genome
## =========================================================
## "Interpolated Manhattan" plots of sorting along the genome, faceted by
## chromosome. Two outputs:
##   (1) moduleA_sorting_manhattan.png -- fraction SORTED per window (any class),
##       one line per sorting threshold tau (nested bands; sorted-at-higher-tau
##       is a subset of lower-tau).
##   (2) moduleA_sorting_manhattan_directional.png -- stacked rows: fraction
##       fixing toward AQUILONIA, toward POLYCTENA, the DIRECTION-UNRESOLVED fraction,
##       the NET direction (aqu - pol), and the diagnostic-index (DI) landscape.
##
## Classification is recomputed at each tau from the stored per-SNP counts in
## moduleA_sorting/data/moduleA_snp.rds, using the sort_rule = "binom" rule:
##   sorted (any)  iff   prop_fixed >= tau                        (magnitude gate)
##   aquilonia     iff   sorted & p_binom < ALPHA & uni_score > 0  (predictable, aqu)
##   polyctena     iff   sorted & p_binom < ALPHA & uni_score < 0
##   unresolved   iff   sorted & p_binom >= ALPHA & n_fixed >= N_POW  (direction unresolved)
## i.e. DIRECTION is the binomial random-direction test (p_binom; tau-independent),
## not the fixed 1/4 split -- a majority-but-not-significant lean is direction-unresolved,
## not unidirectional. No parallelism_stats re-run is needed. Coordinates are
## parsed from the marker id "Chr:Pos".
##
## (A random-direction null overlay was tried here but removed: conditioning on
## the observed amount of near-fixation is easy to misread, and the appropriate
## reference is the simulated neutral null from Module E. Compare against that
## when it is available.)
##
## Reads  moduleA_sorting/data/moduleA_snp.rds
## Writes moduleA_sorting/Figures/moduleA_sorting_manhattan.png,
##        moduleA_sorting/Figures/moduleA_sorting_manhattan_directional.png
## Run from the repo root.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(wesanderson); library(patchwork)
})

## ---- knobs --------------------------------------------------------------
SORT_TH_GRID <- c(0.5, 0.6, 0.7, 0.8)   # magnitude thresholds tau to overlay
ALPHA        <- 0.05                     # binom direction-test significance (matches Module A)
N_POW        <- ceiling(log(ALPHA/2)/log(0.5))  # smallest n_fixed testable (=6 at 0.05)
WIN          <- 100000L                 # window width (bp)
MIN_DENS     <- 50L                     # drop windows with < this many SNPs/100kb (low density)
SMOOTH_K     <- 5L                      # rolling-mean width (in windows) for the trend

## ---- data ---------------------------------------------------------------
snp <- as.data.table(readRDS("moduleA_sorting/data/moduleA_snp.rds"))
d <- snp[differentiated == TRUE & is.finite(uni_score),
         .(marker, uni_score, bi_score, prop_fixed, p_binom, n_fixed, DI)]
d[, c("Chr", "Pos") := tstrsplit(marker, ":", fixed = TRUE)]
d[, Pos := as.integer(Pos)]
d[, `:=`(uni_mag = abs(uni_score), win = Pos %/% WIN)]

## low-density filter: drop 100-kb windows with < MIN_DENS differentiated SNPs,
## uniformly across every row (their per-window sorting/net estimates are noisy)
d[, dens := .N, by = .(Chr, win)]
n_win0 <- uniqueN(d[, .(Chr, win)])
d <- d[dens >= MIN_DENS]
cat(sprintf("SNP-density filter (>= %d SNPs/100kb): dropped %d of %d windows\n",
            MIN_DENS, n_win0 - uniqueN(d[, .(Chr, win)]), n_win0))

## marker-level local LD (ld_w_095) from the map, joined on marker.
## Read from hybrids_only_maf005.Rdata: module0 attaches ld_w_095 to map_hyb_005
## only when re-saving that file (after the hybrids_and_parents save), so the
## parents file's map_hyb_005 lacks the column. Same markers in both.
e <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e)
mp <- as.data.table(e$map_hyb_005)[, .(marker, ld_w_095)]
d[mp, on = "marker", ld_w_095 := i.ld_w_095]; rm(e)

num    <- suppressWarnings(as.integer(sub("Chr", "", unique(d$Chr))))
chr_lv <- unique(d$Chr)[order(num)]

## ---- per-window fractions at each tau (any / aquilonia / polyctena) ------
land <- rbindlist(lapply(SORT_TH_GRID, function(th)
  d[, .(n   = .N,
        any = mean(prop_fixed >= th),
        aqu = mean(prop_fixed >= th & p_binom <  ALPHA & uni_score > 0),
        pol = mean(prop_fixed >= th & p_binom <  ALPHA & uni_score < 0),
        bi  = mean(prop_fixed >= th & p_binom >= ALPHA & n_fixed >= N_POW)),
    by = .(Chr, win)][, sort_th := th][]))
valcols <- c("any", "aqu", "pol", "bi")
land[n < MIN_DENS, (valcols) := NA_real_]
land[, pos_mb := (win * WIN + WIN / 2) / 1e6]
setorder(land, sort_th, Chr, win)
land[, (valcols) := lapply(.SD, frollmean, n = SMOOTH_K, align = "center", na.rm = TRUE),
     by = .(sort_th, Chr), .SDcols = valcols]
land[, net := aqu - pol]                 # net direction (+ = aquilonia, - = polyctena)
land[, Chr := factor(Chr, levels = chr_lv)]

## ---- DI landscape (tau-independent) ------------------------------------
di_land <- d[, .(n = .N, DI = mean(DI, na.rm = TRUE)), by = .(Chr, win)]
di_land[n < MIN_DENS, DI := NA_real_]
di_land[, pos_mb := (win * WIN + WIN / 2) / 1e6]
setorder(di_land, Chr, win)
di_land[, DI_s := frollmean(DI, SMOOTH_K, align = "center", na.rm = TRUE), by = Chr]
di_land[, Chr := factor(Chr, levels = chr_lv)]

## ---- raw marker-level values (points) for the DI and ld_w rows ----------
di_pts <- d[is.finite(DI),
            .(Chr = factor(Chr, levels = chr_lv), pos_mb = Pos / 1e6, DI)]
ld_pts <- d[is.finite(ld_w_095),
            .(Chr = factor(Chr, levels = chr_lv), pos_mb = Pos / 1e6, ld_w_095)]

## windowed rolling-mean line for the ld_w row (matches the DI row's line)
ld_land <- d[, .(n = .N, ldw = mean(ld_w_095, na.rm = TRUE)), by = .(Chr, win)]
ld_land[n < MIN_DENS, ldw := NA_real_]
ld_land[, pos_mb := (win * WIN + WIN / 2) / 1e6]
setorder(ld_land, Chr, win)
ld_land[, ldw_s := frollmean(ldw, SMOOTH_K, align = "center", na.rm = TRUE), by = Chr]
ld_land[, Chr := factor(Chr, levels = chr_lv)]

## ---- shared styling -----------------------------------------------------
pal <- wes_palette("Zissou1", length(SORT_TH_GRID), type = "continuous")
strip_theme <- theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(1, "pt"),
        strip.background = element_blank(), strip.text = element_text(size = 6),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank())
tau_layers <- list(
  facet_grid(. ~ Chr, scales = "free_x", space = "free_x"),
  scale_colour_manual(values = pal),
  scale_x_continuous(expand = expansion(mult = 0.02)),
  labs(colour = expression(tau)), strip_theme)

## ---- (1) any-sorted landscape ------------------------------------------
p_any <- ggplot(land, aes(pos_mb, any, colour = factor(sort_th), group = sort_th)) +
  geom_line(linewidth = 0.35) + tau_layers +
  labs(x = "position (Mbp)", y = "fraction sorted (100 kb window)",
       title = "Ancestry-sorting landscape across sorting thresholds",
       subtitle = sprintf("fraction of differentiated SNPs sorted (prop_fixed >= tau) per %d-kb window (windows < %d SNPs removed as low-density); %d-window rolling mean; fix_major = 0.85; sort_rule = binom (alpha = %.2f)",
                          WIN %/% 1000L, MIN_DENS, SMOOTH_K, ALPHA)) +
  theme(legend.position = "top")
ggsave("moduleA_sorting/Figures/moduleA_sorting_manhattan.png", p_any, width = 16, height = 4, dpi = 300)

## ---- (2) directional (aqu / pol) + net direction + DI -------------------
frac_max <- max(c(land$aqu, land$pol), na.rm = TRUE)         # shared y for aqu vs pol
dir_long <- melt(land, id.vars = c("Chr", "pos_mb", "sort_th"),
                 measure.vars = c("aqu", "pol"),
                 variable.name = "metric", value.name = "value")
dir_long[, metric := factor(metric, levels = c("aqu", "pol"),
         labels = c("toward aquilonia", "toward polyctena"))]
dir_long[, Chr := factor(Chr, levels = chr_lv)]

p_dir_sort <- ggplot(dir_long, aes(pos_mb, value, colour = factor(sort_th), group = sort_th)) +
  geom_line(linewidth = 0.35) +
  facet_grid(metric ~ Chr, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = pal) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(0, frac_max)) +
  labs(x = NULL, y = "fraction fixing toward parent", colour = expression(tau)) +
  strip_theme

p_bi <- ggplot(land, aes(pos_mb, bi, colour = factor(sort_th), group = sort_th)) +
  geom_line(linewidth = 0.35) +
  facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = pal) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  labs(x = NULL, y = "fraction\ndirection-\nunresolved", colour = expression(tau)) +
  strip_theme + theme(strip.text = element_blank())

p_net <- ggplot(land, aes(pos_mb, net, colour = factor(sort_th), group = sort_th)) +

  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.35) +
  facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = pal) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  labs(x = NULL, y = "net direction\n(aqu - pol)", colour = expression(tau)) +
  strip_theme + theme(strip.text = element_blank())

p_di <- ggplot() +
  geom_point(data = di_pts, aes(pos_mb, DI), colour = "grey50", alpha = 0.5,
             size = 0.5, stroke = 0) +
  geom_hline(yintercept = -25, linetype = 2, colour = "#B22222", linewidth = 0.3) +
  geom_line(data = di_land, aes(pos_mb, DI_s, group = Chr),
            linewidth = 0.35, colour = "grey20") +
  facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  coord_cartesian(ylim = c(-100, 0)) +
  labs(x = NULL, y = "DI") +
  strip_theme + theme(strip.text = element_blank())

p_ldw <- ggplot() +
  geom_point(data = ld_pts, aes(pos_mb, ld_w_095), colour = "grey50", alpha = 0.5,
             size = 0.5, stroke = 0) +
  geom_line(data = ld_land, aes(pos_mb, ldw_s, group = Chr),
            linewidth = 0.35, colour = "grey20") +
  facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  labs(x = "position (Mbp)", y = "ld_w(0.95)") +
  strip_theme + theme(strip.text = element_blank())

## no plot title/subtitle: the supplementary caption describes the figure
p_dir <- (p_dir_sort / p_bi / p_net / p_di / p_ldw) +
  plot_layout(heights = c(2, 1, 1, 1, 1), guides = "collect") &
  theme(legend.position = "top")
ggsave("moduleA_sorting/Figures/moduleA_sorting_manhattan_directional.png", p_dir,
       width = 16, height = 7, dpi = 300)

cat("Saved: moduleA_sorting/Figures/moduleA_sorting_manhattan.png,",
    "moduleA_sorting/Figures/moduleA_sorting_manhattan_directional.png\n")
