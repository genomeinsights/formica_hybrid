## =====================================================================
## Figures.R  --  manuscript-ready "plain" figures for Modules 0 and A
## =====================================================================
## Curated, deliberately manual figure script. For each main/supplementary
## figure of Modules 0 (LD decay & complexity reduction) and A (ancestry
## sorting) it READS THE SAVED TIDY DATA produced by the pipeline and rebuilds
## the plot WITHOUT the in-figure descriptive title/subtitle (the manuscript
## supplies those via LaTeX \caption) and with PARSED, descriptive axis / legend
## / facet labels (no bare "ld_w", "prop_fixed", "aqu", "Fst", ... ).
##
## Outputs are collected in a single repo-root Figures/ directory, each file
## clearly prefixed by its module (module0_* / moduleA_*), leaving the titled
## pipeline versions in the module folders untouched:
##   Figures/module0_<name>.{png,pdf}
##   Figures/moduleA_<name>.{png,pdf}
##
## NOTE: the Stage-1-vs-combined figure (figS05) is produced by the LDscnR
## package plotter plot_pruning_comparison() and is fixed AT SOURCE (clean panel
## titles there) rather than duplicated here -- it is intentionally NOT in this file.
##
## Run from the formica_hybrid repo root:  Rscript Figures.R
## (or source specific sections interactively).

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(wesanderson); library(pROC)
})

FIGDIR <- "Figures"                       # single repo-root output directory
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

## ---- shared descriptive / parsed labels ---------------------------------
## species names in italics; statistics with proper subscripts; no underscores.
AQU   <- "italic('F. aquilonia')"     # plotmath string (use with label_parsed / parse())
POL   <- "italic('F. polyctena')"
lab_aqu    <- expression(italic("F. aquilonia"))
lab_pol    <- expression(italic("F. polyctena"))
lab_recomb <- "recombination (cM/Mb, log)"
lab_DI     <- "diagnostic index (DI)"
lab_ldw    <- expression("local LD support " * ld[w] * " (" * rho == 0.95 * ")")
lab_tau    <- expression("sorting threshold " * tau)
lab_phi    <- expression("near-fixation floor " * phi)
lab_pctdiff<- "% of differentiated SNPs"

## a clean base theme; nothing sets plot.title because we never add one.
theme_plain <- function(base = 9)
  theme_classic(base_size = base) +
  theme(axis.title = element_text(size = base),
        axis.text  = element_text(size = base - 2),
        legend.title = element_text(size = base - 1),
        legend.text  = element_text(size = base - 2),
        strip.background = element_blank(),
        strip.text = element_text(size = base - 1))

save_fig <- function(plot, stem, width, height, dpi = 300, pdf = FALSE) {
  png <- file.path(FIGDIR, paste0(stem, ".png"))
  ggsave(png, plot, width = width, height = height, units = "mm", dpi = dpi)
  if (pdf) ggsave(file.path(FIGDIR, paste0(stem, ".pdf")), plot,
                  width = width, height = height, units = "mm")
  cat("Saved:", png, if (pdf) "(+ .pdf)" else "", "\n")
}

wes_grad <- function(n) wes_palette("Zissou1", n, type = "continuous")

## =====================================================================
## MODULE 0
## =====================================================================

## ---- [Fig S: ROC] discrimination of low-recombination regions -----------
## reads the saved ld_w / decay-rate vs map-recombination comparison table and
## rebuilds the ROC curves (lowest 5/10/25% of windows). Plain: no super-title;
## panels labelled by the low-recombination definition; legend series spelled out.
fig_roc <- function() {
  cmp <- readRDS("module0_ld_pruning/data/ldw_a_recmap_comparison.rds")$comp_dt
  q_grid <- c(0.05, 0.10, 0.25)
  roc_plots <- lapply(q_grid, function(q) {
    label <- as.integer(cmp$cM_Mb_med <= quantile(cmp$cM_Mb_med, q, na.rm = TRUE))
    r_ldw <- roc(label, cmp$ld_w_med, direction = "auto", quiet = TRUE)
    r_a   <- roc(label, cmp$a,        direction = "auto", quiet = TRUE)
    ggroc(list("local LD support" = r_ldw, "decay rate a" = r_a),
          linewidth = 0.8) +
      geom_abline(slope = 1, intercept = 1, linetype = 2, colour = "grey60") +
      scale_colour_manual(values = c("#0072B2", "#D55E00"), name = NULL) +
      labs(x = "specificity", y = "sensitivity",
           subtitle = sprintf("lowest %d%% recombination", q * 100)) +
      theme_plain(11) + theme(plot.subtitle = element_text(size = 9, hjust = 0))
  })
  p <- wrap_plots(roc_plots, ncol = 3) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  save_fig(p, "module0_roc_low_recombination", width = 200, height = 85)
}

## ---- [Fig S: LD tracks] ld_w co-localises with the recomb environment ----
## 2x2: rows = example chromosomes, cols = ld_w vs decay-rate a / vs recomb rate.
## Plain: drop the per-panel "ChrNN: ld_w and ..." titles (keep only panel tags a-d).
fig_ld_tracks <- function() {
  CHRS <- c("Chr26", "Chr10")
  aw  <- readRDS("module0_ld_pruning/data/ld_tracks_a_windows.rds")[
           Chr %in% CHRS & regime == "structured" & is.finite(a)]
  ldw <- readRDS("module0_ld_pruning/data/ld_tracks_ldw_persnp.rds")[
           Chr %in% CHRS & is.finite(ld_w_095)]
  rec <- fread("data/Frufa_DTOL_PR.ref_genome.recmap")
  rec[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
  rec <- rec[Chr %in% CHRS, .(Chr, pos, cMMb = `cM/Mb`)]
  sc01 <- function(v) { v <- as.numeric(v); (v - min(v, na.rm = TRUE)) / diff(range(v, na.rm = TRUE)) }
  col_a <- "#D55E00"; col_r <- "#0072B2"; col_ldw <- "grey55"
  panel <- function(ch, which) {
    lw <- ldw[Chr == ch]; lw[, `:=`(x = Pos / 1e6, y = sc01(ld_w_095))]
    if (which == "a") { tr <- aw[Chr == ch][order(mid)]; tr[, `:=`(x = mid / 1e6, y = sc01(a))]; col <- col_a }
    else               { tr <- rec[Chr == ch][order(pos)]; tr[, `:=`(x = pos / 1e6, y = sc01(cMMb))]; col <- col_r }
    ggplot() +
      geom_point(data = lw, aes(x, y), colour = col_ldw, size = 0.22, alpha = 0.35) +
      geom_line(data = tr, aes(x, y), colour = col, linewidth = 0.45) +
      scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
      labs(x = paste0(ch, " position (Mbp)"), y = "scaled to 0-1") +
      theme_plain(8) + theme(plot.tag = element_text(face = "bold", size = 10))
  }
  fig <- (panel("Chr26", "a") | panel("Chr26", "rec")) /
         (panel("Chr10", "a") | panel("Chr10", "rec")) +
    plot_annotation(tag_levels = "a")
  save_fig(fig, "module0_ld_tracks_chr26_chr10", width = 180, height = 110, pdf = TRUE)
}

## ---- [Fig S: fidelity] consensus round-trip fidelity histogram ----------
## reads the >=5-marker construction sweep; plain: parse facet labels (local LD
## support instead of ld_w; cluster size instead of n_loci).
fig_fidelity <- function() {
  Accent <- "#315B7D"; Grey <- "#C7CDD2"; Pol <- "#D08A45"; DROP_TH <- 0.2
  res <- readRDS("module0_ld_pruning/data/results_min_loci5.rds")
  fl <- rbindlist(lapply(res, function(r) {
    g <- r$groups[startsWith(r$groups$group_id, "F")]; g[, .(score, th, n_loci)] }))
  fl <- fl[th != DROP_TH]
  fl[, th_label := sprintf("local LD support threshold = %s  (n = %s, mean = %.3f)",
                           th, format(.N, big.mark = ","), mean(score)), by = th]
  fl[, th_label := factor(th_label, levels = unique(th_label[order(th)]))]
  fl[, cluster_type := fifelse(n_loci == 1L, "single-marker cluster (no merge possible)",
                               "multi-marker cluster (merged)")]
  hb <- fl[, { br <- seq(0.80, 1, by = 0.01); h <- hist(score, breaks = br, plot = FALSE)
    .(bin_mid = h$mids, count = h$counts) }, by = .(th_label, cluster_type)]
  hb <- fl[, .(N = .N), by = th_label][hb, on = "th_label"][, pct := 100 * count / N]
  p <- ggplot(hb, aes(bin_mid, pct, fill = cluster_type)) +
    geom_col(colour = "white", linewidth = 0.1, width = 0.01) +
    geom_vline(xintercept = 0.80, linetype = 2, colour = Pol, linewidth = 0.5) +
    scale_fill_manual(values = c("single-marker cluster (no merge possible)" = Grey,
                                 "multi-marker cluster (merged)" = Accent), name = NULL) +
    facet_wrap(~ th_label, ncol = 1) + coord_cartesian(xlim = c(0.80, 1)) +
    labs(x = expression("consensus fidelity  " * score[eMLG] == cor(round(x), x)^2),
         y = "% of that run's flagged clusters") +
    theme_plain(9) + theme(legend.position = "top")
  save_fig(p, "module0_fidelity_hist", width = 120, height = 150, pdf = TRUE)
}

## =====================================================================
## MODULE A
## =====================================================================

## ---- [Fig S1] per-SNP sorting-classification threshold sweep ------------
## Two directionally-resolved class panels (free_y) + a dedicated panel for the
## direction-unresolved share of sorted loci (the "any sorted" total is omitted).
fig_sorting_sweep <- function() {
  sw <- as.data.table(readRDS("moduleA_sorting/data/moduleA_sortth_fixth_sweep.rds"))
  taus <- sort(unique(sw$sort_th)); pal <- wes_grad(uniqueN(sw$fix_major))
  ## (left) directionally-resolved classes in ONE panel: phi = colour, direction = linetype.
  ## The aquilonia/polyctena split here is POOLED over all DI; their balance is DI-dependent
  ## -- that caveat belongs in the manuscript \caption, not baked into the figure.
  ## NB: the saved sweep still uses the pre-relabel column names.
  long <- melt(sw, id.vars = c("sort_th", "fix_major"),
               measure.vars = c("pct_aquilonia", "pct_polyctena"),
               variable.name = "metric", value.name = "pct")
  long[, direction := factor(metric, levels = c("pct_aquilonia", "pct_polyctena"),
                             labels = c(AQU, POL))]
  p_class <- ggplot(long, aes(sort_th, pct, colour = factor(fix_major),
                              linetype = direction, group = interaction(fix_major, direction))) +
    geom_line() + geom_point(size = 1.0) +
    scale_x_continuous(breaks = taus) +
    scale_colour_manual(values = pal, name = lab_phi) +
    scale_linetype_manual(values = c("solid", "longdash"),
                          labels = function(l) parse(text = l), name = "sorting direction") +
    labs(x = lab_tau, y = lab_pctdiff) + theme_plain(9)
  ## (right) direction-unresolved as a fraction of SORTED loci (unresolved / sorted)
  sw[, prop_unres := unresolved / sorted]
  p_prop <- ggplot(sw, aes(sort_th, prop_unres, colour = factor(fix_major), group = fix_major)) +
    geom_line() + geom_point(size = 1.0) +
    scale_x_continuous(breaks = taus) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_colour_manual(values = pal, name = lab_phi) +
    labs(x = lab_tau, y = "direction unresolved\n(% of sorted loci)") + theme_plain(9)
  p <- (p_class | p_prop) + plot_layout(widths = c(1.25, 1), guides = "collect") &
    theme(legend.position = "right")
  save_fig(p, "moduleA_sorting_sweep", width = 195, height = 62)
}

## ---- [Fig S2] DI-decile sorting direction, threshold sweep --------------
fig_panelB_sweep <- function() {
  pb <- readRDS("moduleA_sorting/data/moduleA_panelB_sweep.rds")
  p <- ggplot(pb, aes(DI_decile, frac_aqu_of_unidir, colour = factor(sort_th), group = sort_th)) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey55") +
    geom_line() + geom_point(size = 1.1) +
    facet_wrap(~ fix_major, nrow = 1,
               labeller = as_labeller(function(x) paste0("phi == '",
                 formatC(as.numeric(x), format = "f", digits = 2), "'"), default = label_parsed)) +
    scale_x_continuous(breaks = 1:10) + scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = wes_grad(uniqueN(pb$sort_th)), name = lab_tau) +
    labs(x = "diagnostic-index decile (low to high)",
         y = expression("fraction fixing toward " * italic("F. aquilonia"))) +
    theme_plain(9)
  save_fig(p, "moduleA_panelB_sweep", width = 200, height = 62)
}

## ---- [Fig S3] direction-model coefficients vs sorting threshold ---------
fig_direction_sweep <- function() {
  ds <- as.data.table(readRDS("moduleA_sorting/data/moduleA_architecture.rds")$direction_sweep)
  ## term is already stored in descriptive form; fix the panel order only.
  ds[, term := factor(term, levels = c("diagnostic index", "recombination", "parental MAF", "log block size"))]
  p <- ggplot(ds, aes(sort_th, estimate)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_line(colour = "#315B7D") + geom_point(colour = "#315B7D", size = 1.8) +
    geom_errorbar(aes(ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se),
                  width = 0.02, colour = "#315B7D") +
    facet_wrap(~ term, nrow = 1, scales = "free_y") +
    scale_x_continuous(breaks = sort(unique(ds$sort_th))) +
    labs(x = lab_tau, y = "standardised coefficient (95% CI)") +
    theme_plain(9)
  save_fig(p, "moduleA_direction_sweep", width = 200, height = 70)
}

## ---- [Fig 1, main text] genomic architecture of sorting (a/b/c) ---------
fig_architecture <- function() {
  a  <- readRDS("moduleA_sorting/data/moduleA_architecture.rds")
  zsc <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  tt <- theme_plain(9) + theme(plot.tag = element_text(face = "bold", size = 11),
                               legend.position = "bottom", legend.title = element_blank(),
                               legend.key.size = unit(3, "mm"))
  ## (a) standardised Fst / dxy / within-parent H_E across recombination deciles
  cha <- melt(a$arch_tab[, .(med_recomb, Fst = zsc(Fst), dxy = zsc(dxy), pi_within = zsc(pi_within))],
              id.vars = "med_recomb", variable.name = "metric", value.name = "z")
  cha[, metric := factor(metric, levels = c("Fst", "dxy", "pi_within"),
                         labels = c("F[ST]", "d[xy]", "'within-parent'~H[E]"))]
  p_a <- ggplot(cha, aes(med_recomb, z, colour = metric)) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.6) + scale_x_log10() +
    scale_colour_manual(values = c("F[ST]" = "#D55E00", "d[xy]" = "#009E73",
                                   "'within-parent'~H[E]" = "#0072B2"),
                        labels = function(l) parse(text = l)) +
    labs(x = lab_recomb, y = "standardised value (z)") + tt
  ## (b) fraction sorted vs recombination: LD-reduced unit vs SNP (stored separately)
  cmp <- rbind(a$unit_by_recomb, a$snp_by_recomb)
  cmp[, level := factor(level, levels = c("SNP", "LD-reduced unit"))]
  p_b <- ggplot(cmp, aes(med_r, frac_sorted, colour = level)) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.6) + scale_x_log10() +
    scale_colour_manual(values = c("LD-reduced unit" = "#315B7D", "SNP" = "#C7CDD2")) +
    labs(x = lab_recomb, y = "fraction sorted") + tt
  ## (c) direction vs recombination: LD-reduced unit vs SNP (both levels, as in panel b).
  ## The SNP series spikes toward F. aquilonia in the lowest-recombination decile
  ## (pseudoreplication in a few large blocks); independent units are flat. Legend is
  ## suppressed here -- panel b's SNP/unit legend applies.
  cmp_c <- rbind(a$unit_by_recomb, a$snp_by_recomb)
  cmp_c[, level := factor(level, levels = c("SNP", "LD-reduced unit"))]
  p_c <- ggplot(cmp_c, aes(med_r, frac_aqu_of_unidir, colour = level)) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
    geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
    scale_x_log10() + scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = c("LD-reduced unit" = "#315B7D", "SNP" = "#C7CDD2")) +
    labs(x = lab_recomb, y = expression("fraction fixing toward " * italic("F. aquilonia"))) +
    tt + theme(legend.position = "none")
  fig <- p_a + p_b + p_c + plot_annotation(tag_levels = "a")
  save_fig(fig, "moduleA_architecture_fig", width = 180, height = 78, pdf = TRUE)
}

## ---- [Fig S4] directional ancestry-sorting landscape (stacked rows) -----
## Reclassifies per-SNP counts at each tau (binom rule) and builds the 5-row
## stack. Already title-free in the pipeline; here we parse the row/axis labels
## (local LD support, diagnostic index, species names, "near-fixed").
fig_manhattan_directional <- function() {
  SORT_TH_GRID <- c(0.5, 0.6, 0.7, 0.8); ALPHA <- 0.05
  N_POW <- ceiling(log(ALPHA / 2) / log(0.5)); WIN <- 100000L; MIN_DENS <- 50L; SMOOTH_K <- 5L
  snp <- as.data.table(readRDS("moduleA_sorting/data/moduleA_snp.rds"))
  d <- snp[differentiated == TRUE & is.finite(uni_score),
           .(marker, uni_score, bi_score, prop_fixed, p_binom, n_fixed, DI)]
  d[, c("Chr", "Pos") := tstrsplit(marker, ":", fixed = TRUE)]; d[, Pos := as.integer(Pos)]
  d[, win := Pos %/% WIN]
  d[, dens := .N, by = .(Chr, win)]; d <- d[dens >= MIN_DENS]
  e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
  mp <- as.data.table(e$map_hyb_005)[, .(marker, ld_w_095)]
  d[mp, on = "marker", ld_w_095 := i.ld_w_095]; rm(e)
  num <- suppressWarnings(as.integer(sub("Chr", "", unique(d$Chr)))); chr_lv <- unique(d$Chr)[order(num)]

  land <- rbindlist(lapply(SORT_TH_GRID, function(th)
    d[, .(n = .N, any = mean(prop_fixed >= th),
          aqu = mean(prop_fixed >= th & p_binom <  ALPHA & uni_score > 0),
          pol = mean(prop_fixed >= th & p_binom <  ALPHA & uni_score < 0),
          bi  = mean(prop_fixed >= th & p_binom >= ALPHA & n_fixed >= N_POW)),
      by = .(Chr, win)][, sort_th := th][]))
  vc <- c("any", "aqu", "pol", "bi"); land[n < MIN_DENS, (vc) := NA_real_]
  land[, pos_mb := (win * WIN + WIN / 2) / 1e6]; setorder(land, sort_th, Chr, win)
  land[, (vc) := lapply(.SD, frollmean, n = SMOOTH_K, align = "center", na.rm = TRUE),
       by = .(sort_th, Chr), .SDcols = vc]
  land[, net := aqu - pol]; land[, Chr := factor(Chr, levels = chr_lv)]

  di_land <- d[, .(n = .N, DI = mean(DI, na.rm = TRUE)), by = .(Chr, win)]
  di_land[n < MIN_DENS, DI := NA_real_]; di_land[, pos_mb := (win * WIN + WIN / 2) / 1e6]
  setorder(di_land, Chr, win); di_land[, DI_s := frollmean(DI, SMOOTH_K, align = "center", na.rm = TRUE), by = Chr]
  di_land[, Chr := factor(Chr, levels = chr_lv)]
  di_pts <- d[is.finite(DI), .(Chr = factor(Chr, levels = chr_lv), pos_mb = Pos / 1e6, DI)]
  ld_pts <- d[is.finite(ld_w_095), .(Chr = factor(Chr, levels = chr_lv), pos_mb = Pos / 1e6, ld_w_095)]
  ld_land <- d[, .(n = .N, ldw = mean(ld_w_095, na.rm = TRUE)), by = .(Chr, win)]
  ld_land[n < MIN_DENS, ldw := NA_real_]; ld_land[, pos_mb := (win * WIN + WIN / 2) / 1e6]
  setorder(ld_land, Chr, win); ld_land[, ldw_s := frollmean(ldw, SMOOTH_K, align = "center", na.rm = TRUE), by = Chr]
  ld_land[, Chr := factor(Chr, levels = chr_lv)]

  pal <- wes_grad(length(SORT_TH_GRID))
  strip_theme <- theme_bw(base_size = 8) +
    theme(panel.spacing.x = unit(1, "pt"), strip.background = element_blank(),
          strip.text = element_text(size = 6), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), panel.grid.minor = element_blank())
  frac_max <- max(c(land$aqu, land$pol), na.rm = TRUE)
  dir_long <- melt(land, id.vars = c("Chr", "pos_mb", "sort_th"), measure.vars = c("aqu", "pol"),
                   variable.name = "metric", value.name = "value")
  dir_long[, metric := factor(metric, levels = c("aqu", "pol"),
           labels = c("toward~italic('F. aquilonia')", "toward~italic('F. polyctena')"))]
  dir_long[, Chr := factor(Chr, levels = chr_lv)]

  p_dir_sort <- ggplot(dir_long, aes(pos_mb, value, colour = factor(sort_th), group = sort_th)) +
    geom_line(linewidth = 0.35) +
    facet_grid(metric ~ Chr, scales = "free_x", space = "free_x", labeller = labeller(metric = label_parsed)) +
    scale_colour_manual(values = pal) + scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, frac_max)) +
    labs(x = NULL, y = "fraction near-fixed toward parent", colour = expression(tau)) + strip_theme
  p_bi <- ggplot(land, aes(pos_mb, bi, colour = factor(sort_th), group = sort_th)) +
    geom_line(linewidth = 0.35) + facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
    scale_colour_manual(values = pal) + scale_x_continuous(expand = expansion(mult = 0.02)) +
    labs(x = NULL, y = "fraction\ndirection\nunresolved", colour = expression(tau)) +
    strip_theme + theme(strip.text = element_blank())
  p_net <- ggplot(land, aes(pos_mb, net, colour = factor(sort_th), group = sort_th)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) + geom_line(linewidth = 0.35) +
    facet_grid(. ~ Chr, scales = "free_x", space = "free_x") + scale_colour_manual(values = pal) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    labs(x = NULL, y = "net direction\n(aquilonia - polyctena)", colour = expression(tau)) +
    strip_theme + theme(strip.text = element_blank())
  p_di <- ggplot() +
    geom_point(data = di_pts, aes(pos_mb, DI), colour = "grey50", alpha = 0.5, size = 0.5, stroke = 0) +
    geom_hline(yintercept = -25, linetype = 2, colour = "#B22222", linewidth = 0.3) +
    geom_line(data = di_land, aes(pos_mb, DI_s, group = Chr), linewidth = 0.35, colour = "grey20") +
    facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
    scale_x_continuous(expand = expansion(mult = 0.02)) + coord_cartesian(ylim = c(-100, 0)) +
    labs(x = NULL, y = "diagnostic\nindex (DI)") + strip_theme + theme(strip.text = element_blank())
  p_ldw <- ggplot() +
    geom_point(data = ld_pts, aes(pos_mb, ld_w_095), colour = "grey50", alpha = 0.5, size = 0.5, stroke = 0) +
    geom_line(data = ld_land, aes(pos_mb, ldw_s, group = Chr), linewidth = 0.35, colour = "grey20") +
    facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    labs(x = "position (Mbp)", y = expression(atop("local LD", "support " * ld[w]))) +
    strip_theme + theme(strip.text = element_blank())
  p <- (p_dir_sort / p_bi / p_net / p_di / p_ldw) +
    plot_layout(heights = c(2, 1, 1, 1, 1), guides = "collect") & theme(legend.position = "top")
  save_fig(p, "moduleA_sorting_manhattan_directional", width = 400, height = 180)
}

## =====================================================================
## SUPPLEMENTARY TABLES
## =====================================================================

## ---- [Table S] genomic architecture across recombination deciles --------
## Compact decile table matching Fig 1a: median recombination, DI, Fst, dxy,
## within-parent H_E, mean cluster size. Writes a booktabs LaTeX table (needs
## \usepackage{booktabs}) to Tables/moduleA_architecture_deciles.tex.
make_arch_table <- function() {
  dir.create("Tables", showWarnings = FALSE, recursive = TRUE)
  at <- as.data.table(readRDS("moduleA_sorting/data/moduleA_architecture.rds")$arch_tab)
  setorder(at, rbin)
  rows <- at[, sprintf("%d & %.2f & %.1f & %.3f & %.3f & %.3f & %.1f \\\\",
                       rbin, med_recomb, DI, Fst, dxy, pi_within, cluster_size)]
  tex <- c(
    "\\begin{table}[h]",
    "\\centering",
    "\\caption{\\textbf{Genomic architecture across recombination-rate deciles.}",
    "Per decile of map recombination rate (decile~1 = lowest): median recombination",
    "rate, mean marker-level diagnostic index (DI), parental relative differentiation",
    "$F_{\\mathrm{ST}}$, between-parent allelic difference $d_{xy}$, within-parent",
    "expected heterozygosity $H_E$, and mean LD-cluster size. Allele-frequency",
    "summaries are conditional on the retained SNPs.}",
    "\\label{tabS:architecture_deciles}",
    "\\small",
    "\\begin{tabular}{ccccccc}",
    "\\toprule",
    "Decile & Median recomb. & DI & $F_{\\mathrm{ST}}$ & $d_{xy}$ & Within-parent $H_E$ & Mean cluster \\\\",
    "       & (cM/Mb)        &    &                   &          &                     & size \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}")
  writeLines(tex, "Tables/moduleA_architecture_deciles.tex")
  cat("Saved: Tables/moduleA_architecture_deciles.tex\n")
}

## ---- [Table S] sorting across recombination deciles, by analytical level -
## Backs Fig 1b/c: fraction sorted and fraction toward F. aquilonia per decile,
## for LD-reduced units vs the random SNP sample. The decile-1 SNP inflation is
## the spatial-pseudoreplication signature. Uses \Faq (manuscript macro), booktabs
## + \cmidrule. Writes Tables/moduleA_sorting_by_level.tex.
make_sorting_level_table <- function() {
  dir.create("Tables", showWarnings = FALSE, recursive = TRUE)
  a  <- readRDS("moduleA_sorting/data/moduleA_architecture.rds")
  u  <- as.data.table(a$unit_by_recomb)[order(rbin)]
  s  <- as.data.table(a$snp_by_recomb)[order(rbin)]
  at <- as.data.table(a$arch_tab)[order(rbin)]
  stopifnot("decile rows misaligned across objects" =
              all(u$rbin == s$rbin) && all(u$rbin == at$rbin))
  rows <- sprintf("%d & %.2f & %.3f & %.3f & %.3f & %.3f \\\\",
                  at$rbin, at$med_recomb, u$frac_sorted, s$frac_sorted,
                  u$frac_aqu_of_unidir, s$frac_aqu_of_unidir)
  tex <- c(
    "\\begin{table}[h]",
    "\\centering",
    "\\caption{\\textbf{Ancestry sorting across recombination-rate deciles, by",
    "analytical level.} Per decile (decile~1 = lowest recombination): median",
    "recombination rate, the fraction of eligible loci sorted, and the fraction of",
    "directionally resolved loci fixing toward \\Faq{}, for LD-reduced units and for a",
    "random sample of individual markers (SNP). The inflated SNP values in decile~1",
    "reflect spatial pseudoreplication---redundant markers within a few large",
    "low-recombination blocks.}",
    "\\label{tabS:sorting_by_level}",
    "\\begin{tabular}{cccccc}",
    "\\toprule",
    " & & \\multicolumn{2}{c}{Fraction sorted} & \\multicolumn{2}{c}{Fraction toward \\Faq{}} \\\\",
    "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
    "Decile & Median recomb. (cM/Mb) & Unit & SNP & Unit & SNP \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}")
  writeLines(tex, "Tables/moduleA_sorting_by_level.tex")
  cat("Saved: Tables/moduleA_sorting_by_level.tex\n")
}

## ---- [Table] architecture regression models (magnitude + direction) -----
## Combines the two unit-level models so their coefficients need not sit in prose;
## supersedes/merges with the manuscript's direction-model table. booktabs.
make_models_table <- function() {
  dir.create("Tables", showWarnings = FALSE, recursive = TRUE)
  a   <- readRDS("moduleA_sorting/data/moduleA_architecture.rds")
  mag <- coef(summary(a$lm_magnitude)); dr <- coef(summary(a$glm_direction))
  ## SE at 4 dp so the tiny magnitude SEs (~0.0004) do not display as 0.000.
  row <- function(label, est, se, stat, sy)
    sprintf("\\quad %s & $%s$ & $%s$ & $%s=%s$ \\\\", label,
            formatC(est, format = "f", digits = 3), formatC(se, format = "f", digits = 4),
            sy, formatC(stat, format = "f", digits = 1))
  tex <- c(
    "\\begin{table}[!ht]", "\\centering",
    "\\caption{\\textbf{Regression models of ancestry sorting among LD-reduced units.}",
    "Standardised coefficients. \\emph{Magnitude}: linear model of the sorted fraction",
    "$p_{\\mathrm{fixed}}$ on recombination and DI. \\emph{Direction}: logistic model of",
    "the probability of fixing toward \\Faq{} among directionally resolved units, on DI,",
    "recombination, parental MAF and log cluster size. All predictors standardised over",
    "the differentiated-unit set; the statistic is a Wald $t$ (linear) or $z$ (logistic),",
    "and all terms are significant (weakest, parental MAF: $P=0.009$).}",
    "\\label{tab:architecture-models}",
    "\\begin{tabular}{lrrr}", "\\toprule",
    "Predictor & Estimate & Std.\\ error & Statistic \\\\", "\\midrule",
    "\\multicolumn{4}{l}{\\emph{Magnitude} ($p_{\\mathrm{fixed}}$; linear)} \\\\",
    row("$\\log_{10}$(recombination)", mag["zr", "Estimate"],  mag["zr", "Std. Error"],  mag["zr", "t value"],  "t"),
    row("Diagnostic index",            mag["zDI", "Estimate"], mag["zDI", "Std. Error"], mag["zDI", "t value"], "t"),
    "\\addlinespace",
    "\\multicolumn{4}{l}{\\emph{Direction} ($P(\\text{fix}\\rightarrow\\Faq)$; logistic)} \\\\",
    row("Diagnostic index",            dr["zDI", "Estimate"],  dr["zDI", "Std. Error"],  dr["zDI", "z value"],  "z"),
    row("$\\log_{10}$(recombination)", dr["zr", "Estimate"],   dr["zr", "Std. Error"],   dr["zr", "z value"],   "z"),
    row("Pooled parental MAF",         dr["zmaf", "Estimate"], dr["zmaf", "Std. Error"], dr["zmaf", "z value"], "z"),
    row("$\\log_{10}$(cluster size)",  dr["zcs", "Estimate"],  dr["zcs", "Std. Error"],  dr["zcs", "z value"],  "z"),
    "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex, "Tables/moduleA_architecture_models.tex")
  cat("Saved: Tables/moduleA_architecture_models.tex\n")
}

## ---- [Table S] admixture-geometry exposure by DI decile x phi -----------
## Reads the exposure diagnostic (moduleA_admixture_exposure.R must have run) and
## writes the DI-decile x phi table: fraction of observed near-fixed calls whose
## admixture expectation q_tilde = h*p_aqu + (1-h)*p_pol also clears phi (i.e. could
## be geometry, not sorting). booktabs. Skipped with a note if the .rds is absent.
make_exposure_table <- function() {
  f <- "moduleA_sorting/data/moduleA_admixture_exposure.rds"
  if (!file.exists(f)) { cat("[skip] exposure table:", f, "not found (run moduleA_admixture_exposure.R)\n"); return(invisible()) }
  dir.create("Tables", showWarnings = FALSE, recursive = TRUE)
  d <- as.data.table(readRDS(f)$tab_dec)
  phis <- sort(unique(d$phi))
  w  <- dcast(d, dec ~ phi, value.var = "frac_exp")[order(dec)]
  ov <- d[, .(f = sum(frac_exp * n_call) / sum(n_call)), by = phi][order(phi)]   # call-weighted overall
  fmt <- function(x) ifelse(x == 0, "0", formatC(x, format = "f", digits = 3))
  declab <- function(i) ifelse(i == 1L, "1 (low)", ifelse(i == 10L, "10 (high)", as.character(i)))
  body <- vapply(seq_len(nrow(w)), function(r)
    sprintf("%s & %s \\\\", declab(w$dec[r]),
            paste(vapply(as.character(phis), function(p) fmt(w[[p]][r]), character(1)), collapse = " & ")),
    character(1))
  overall <- sprintf("\\textbf{overall} & %s \\\\",
    paste(vapply(ov$f, function(x) if (x == 0) "\\textbf{0}" else sprintf("\\textbf{%.3f}", x), character(1)), collapse = " & "))
  tex <- c(
    "\\begin{table}[h]", "\\centering",
    "\\caption{\\textbf{Admixture-geometry exposure by diagnostic-index decile.} Fraction of",
    "observed near-fixed calls whose per-population admixture expectation $\\tilde q = h\\,p_{\\mathrm{aqu}}",
    "+ (1-h)\\,p_{\\mathrm{pol}}$ also clears the near-fixation floor $\\phi$---i.e.\\ that could arise from",
    "admixture geometry rather than sorting---by DI decile (1 = lowest) and $\\phi$. Exposure is zero at",
    "$\\phi \\ge 0.90$ and, at $\\phi = 0.85$, concentrated in the intermediate deciles rather than the",
    "low/high deciles that drive the DI-dependent sorting direction.}",
    "\\label{tabS:admixture_exposure}",
    paste0("\\begin{tabular}{c", strrep("c", length(phis)), "}"), "\\toprule",
    paste0("DI decile & ", paste(sprintf("$\\phi=%.2f$", phis), collapse = " & "), " \\\\"), "\\midrule",
    body, "\\midrule", overall, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex, "Tables/moduleA_admixture_exposure.tex")
  cat("Saved: Tables/moduleA_admixture_exposure.tex\n")
}

## =====================================================================
## MODULE B — climate-GEA Manhattans + ancestry confound
## =====================================================================
## Title-free reproductions of the Module B figures (counts that were in the titles
## belong in the LaTeX caption). The SNP-level Manhattans read a slim per-SNP BF
## cache built ONCE from the 141 MB BayPass .out files.
MB_D       <- "moduleB_climate_GEA"
MB_ASSOC   <- file.path(MB_D, "data/moduleB_eMLG_association.rds")
MB_BFCACHE <- file.path(MB_D, "data/moduleB_snp_bf.rds")
MB_CLUST   <- "module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds"
## canonical best-SNP marker per cluster (group_id -> best_marker); used to anchor
## each eMLG at its real representative SNP rather than the centrality representative.
mB_bestmk  <- { bs <- as.data.table(readRDS("module0_ld_pruning/data/eMLG_5loci_0025_cM05_bestsnp.rds")$rep_snp_all)
                setNames(bs$rep_snp, bs$group_id) }
MB_SIG <- 15; MB_EMLG <- 15; MB_MIN5 <- 5L
mB_clpal <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33",
              "#A65628","#F781BF","#1B9E77","#D95F02","#7570B3","#66A61E")
mB_theme <- theme_bw(base_size = 9) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

## slim per-SNP BF cache (marker/Chr/Pos/BF1/BF2), built once from raw BayPass .out
mB_snp_bf <- function() {
  if (file.exists(MB_BFCACHE)) return(as.data.table(readRDS(MB_BFCACHE)))
  cat("[build] slim SNP-BF cache from raw BayPass .out (one-time, ~300 MB read)...\n")
  ## map from hybrids_and_parents (hybrids_only_maf005.Rdata is corrupt); same canonical map_hyb_005.
  e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
  map <- as.data.table(e$map_hyb_005)[, .(marker, Chr, Pos)]; rm(e); invisible(gc())
  rd <- function(pc) { r <- fread(sprintf("aland_excluded/%s_DIEM_withOmega_summary_betai_reg.out", pc))
    stopifnot(nrow(r) == nrow(map), identical(r$MRK, seq_len(nrow(r)))); r$`BF(dB)` }
  map[, `:=`(BF1 = rd("PC1"), BF2 = rd("PC2"))]
  saveRDS(map, MB_BFCACHE); cat("[cache] wrote", MB_BFCACHE, "\n"); map[]
}

## shared genome x-coordinate transform (from all markers, so every B figure aligns)
mB_coords <- function(snp) {
  cl <- snp[, .(len = max(Pos)), by = Chr]
  cl[, chr_num := suppressWarnings(as.integer(gsub("[^0-9]", "", Chr)))]; setorder(cl, chr_num)
  sp <- 0.01 * sum(cl$len); cl[, start := c(0, head(cumsum(len + sp), -1))]
  sx <- cl[, .(Chr, start)][snp, on = "Chr"][, x := Pos + start]
  list(clen = cl, xkey = setNames(sx$x, sx$marker),
       shade = cl[chr_num %% 2 == 0, .(xmin = start, xmax = start + len)],
       cmid  = cl[, .(mid = start + len / 2, lab = gsub("^Chr", "", Chr))])
}
mB_xscale <- function(co) scale_x_continuous(breaks = co$cmid$mid, labels = co$cmid$lab, expand = c(0.01, 0.01))

## ---- [Fig] eMLG-level climate Manhattan (categories + FDR floor-survivor triangles)
fig_moduleB_eMLG_manhattan <- function() {
  dt <- as.data.table(readRDS(MB_ASSOC)); co <- mB_coords(mB_snp_bf())
  dt <- co$clen[, .(Chr, start)][dt, on = "Chr"][, x := Pos + start]
  lev <- c("ns", "1-4 sig SNPs", "eMLG-only (0 sig SNPs)", ">=5 sig SNPs (candidate)")
  pal <- c("ns"="grey78", "eMLG-only (0 sig SNPs)"="#D81B60", "1-4 sig SNPs"="#1E88E5", ">=5 sig SNPs (candidate)"="#000000")
  panel <- function(bf, ct, fl, lab) {
    d <- dt[, .(x, BF = get(bf), cat = factor(get(ct), levels = lev), fl = get(fl))]; setorder(d, cat)
    ggplot() +
      geom_rect(data = co$shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey93") +
      geom_point(data = d[cat == "ns"], aes(x, BF), colour = "grey78", size = 0.4) +
      geom_point(data = d[cat != "ns" & !fl], aes(x, BF, colour = cat), size = 1.5) +
      geom_point(data = d[fl == TRUE], aes(x, BF, colour = cat), shape = 17, size = 3.2, show.legend = FALSE) +
      geom_point(data = d[fl == TRUE], aes(x, BF), shape = 2, size = 3.2, colour = "grey15", stroke = 0.5, show.legend = FALSE) +
      geom_hline(yintercept = MB_EMLG, linetype = 2, colour = "red", linewidth = 0.4) +
      scale_colour_manual(values = pal, name = NULL, drop = FALSE) + mB_xscale(co) +
      labs(x = if (lab == "PC2") "chromosome" else NULL, y = sprintf("eMLG BF (dB) - %s", lab)) +
      mB_theme + theme(legend.position = "top")
  }
  p <- panel("eBF1","cat1","floor1","PC1") / panel("eBF2","cat2","floor2","PC2") +
    plot_layout(guides = "collect") & theme(legend.position = "top")
  save_fig(p, "moduleB_eMLG_manhattan", width = 380, height = 175, dpi = 150)
}

## ---- [Fig] all-SNP Manhattan; member SNPs of FDR floor-survivor clusters coloured
fig_moduleB_fdr_snp_manhattan <- function() {
  snp <- mB_snp_bf(); co <- mB_coords(snp)
  dt  <- as.data.table(readRDS(MB_ASSOC))
  m2g <- readRDS(MB_CLUST)$groups[has_eMLG == TRUE, .(marker = unlist(members)), by = group_id]
  fpal <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#1B9E77","#D95F02","#7570B3")
  panel <- function(bfcol, fdr_groups, lab) {
    bfk <- setNames(snp[[bfcol]], snp$marker)
    bg  <- data.table(x = co$xkey[snp$marker], BF = snp[[bfcol]])
    fg  <- m2g[group_id %in% fdr_groups, .(marker, group_id)][
             , `:=`(x = co$xkey[marker], BF = bfk[marker], grp = factor(group_id, levels = fdr_groups))]
    ggplot() +
      geom_rect(data = co$shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey93") +
      geom_point(data = bg, aes(x, BF), colour = "grey80", size = 0.25, alpha = 0.5) +
      geom_point(data = fg, aes(x, BF, colour = grp), size = 1.1) +
      geom_hline(yintercept = MB_SIG, linetype = 2, colour = "red", linewidth = 0.4) +
      scale_colour_manual(values = setNames(fpal[seq_along(fdr_groups)], fdr_groups), name = "FDR eMLG", drop = FALSE) +
      mB_xscale(co) +
      labs(x = if (lab == "PC2") "chromosome" else NULL, y = sprintf("SNP BF (dB) - %s", lab)) +
      mB_theme + theme(legend.position = "right")
  }
  p <- panel("BF1", dt[floor1 == TRUE, group_id], "PC1") /
       panel("BF2", dt[floor2 == TRUE, group_id], "PC2")
  save_fig(p, "moduleB_fdr_snp_manhattan", width = 380, height = 175, dpi = 150)
}

## ---- [Fig] eMLG-filtered SNP Manhattan; clusters with >=5 loci BF>=15 coloured (PC1 & PC2)
fig_moduleB_snp_manhattan_clustered <- function() {
  snp <- mB_snp_bf(); co <- mB_coords(snp)
  he  <- readRDS(MB_CLUST)$groups[has_eMLG == TRUE]
  m2g <- he[, .(marker = unlist(members)), by = group_id]
  gp  <- snp[he[, .(group_id, marker = mB_bestmk[group_id])], on = "marker"][, .(group_id, Chr, Pos)]
  gp[, chr_num := suppressWarnings(as.integer(gsub("[^0-9]", "", Chr)))]; setorder(gp, chr_num, Pos)
  gp[, col := mB_clpal[(seq_len(.N) - 1) %% length(mB_clpal) + 1]]           # global, position-ordered
  m2g <- gp[, .(group_id, col)][m2g, on = "group_id"]
  emlg <- snp[, .(marker, x = co$xkey[marker], BF1, BF2)][m2g, on = "marker"]  # members only
  one <- function(bfcol, lab) {
    d   <- emlg[, .(group_id, col, x, BF = get(bfcol))]
    sig <- d[, .(n = sum(BF >= MB_SIG, na.rm = TRUE)), by = group_id][n >= MB_MIN5, group_id]
    bg  <- d[!group_id %in% sig]; fg <- d[group_id %in% sig][order(x)]
    p <- ggplot() +
      geom_rect(data = co$shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey85", alpha = 0.5) +
      geom_point(data = bg, aes(x, BF), size = 0.4, colour = "grey60", alpha = 0.5) +
      geom_point(data = fg, aes(x, BF, colour = col), size = 0.6, alpha = 0.9) + scale_colour_identity() +
      geom_hline(yintercept = MB_SIG, linetype = 2, colour = "red", linewidth = 0.4) +
      mB_xscale(co) + labs(x = "chromosome", y = "BF (dB)") + mB_theme
    save_fig(p, sprintf("moduleB_snp_manhattan_%s", lab), width = 380, height = 110, dpi = 150)
  }
  one("BF1", "PC1"); one("BF2", "PC2")
}

## ---- [Fig] climate PC vs genome-wide aquilonia ancestry (ancestry confound)
fig_moduleB_ancestry <- function() {
  m <- as.data.table(readRDS(file.path(MB_D, "data/moduleB_ancestry_confound.rds")))
  long <- melt(m, id.vars = c("Population", "ancestry"), measure.vars = c("PC1", "PC2"),
               variable.name = "PC", value.name = "score")
  long[, extreme := Population %in% c("Aland", "Sielva")]
  p <- ggplot(long, aes(ancestry, score)) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey70", linewidth = 0.5) +
    geom_point(aes(colour = extreme), size = 1.8) +
    geom_text(data = long[extreme == TRUE], aes(label = Population), vjust = -0.8, size = 2.6) +
    scale_colour_manual(values = c("FALSE" = "grey40", "TRUE" = "#D55E00"), guide = "none") +
    facet_wrap(~ PC) +
    labs(x = expression("genome-wide " * italic("F. aquilonia") * " ancestry (per population)"),
         y = "climate PC score") + theme_plain(9)
  save_fig(p, "moduleB_ancestry_confound", width = 180, height = 80)
}

## ---- [Fig] all six Manhattan panels stacked (patchwork) -----------------
## a,b eMLG-filtered SNP, clusters of >=5 loci BF>=15 coloured; c,d all-SNP with
## FDR floor-survivor cluster members coloured; e,f eMLG-level (categories + FDR
## floor-survivor triangles). Shared x; eMLG colour + shape legend collected at bottom.
fig_moduleB_manhattans <- function() {
  snp <- mB_snp_bf(); co <- mB_coords(snp)
  dt  <- as.data.table(readRDS(MB_ASSOC))
  he  <- readRDS(MB_CLUST)$groups[has_eMLG == TRUE]
  m2g <- he[, .(marker = unlist(members)), by = group_id]
  dtx <- co$clen[, .(Chr, start)][dt, on = "Chr"][, x := Pos + start]
  strip_x <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  shade_l <- function(fill) geom_rect(data = co$shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = fill)
  hthr    <- geom_hline(yintercept = MB_SIG, linetype = 2, colour = "red", linewidth = 0.35)

  ## a,b -- eMLG-level
  lev <- c("ns", "1-4 sig SNPs", "eMLG-only (0 sig SNPs)", ">=5 sig SNPs (candidate)")
  epal <- c("ns"="grey60", "eMLG-only (0 sig SNPs)"="#D81B60", "1-4 sig SNPs"="#1E88E5", ">=5 sig SNPs (candidate)"="#000000")
  emlg_p <- function(bf, ct, fl, lab, bottom = FALSE, showleg = TRUE) {
    d <- dtx[, .(x, BF = get(bf), cat = factor(get(ct), levels = lev), fl = get(fl))]; setorder(d, cat)
    p <- ggplot() + shade_l("grey85") +
      geom_point(data = d[cat == "ns"], aes(x, BF), colour = "grey60", size = 0.3) +
      geom_point(data = d[cat != "ns" & !fl], aes(x, BF, colour = cat), size = 1.1) +
      geom_point(data = d[fl == TRUE], aes(x, BF, colour = cat, shape = "FDR floor-survivor"), size = 2.4) +
      geom_point(data = d[fl == TRUE], aes(x, BF), shape = 2, size = 2.4, colour = "grey15", stroke = 0.4) +
      hthr + scale_colour_manual(values = epal, name = NULL, drop = FALSE) +
      scale_shape_manual(values = c("FDR floor-survivor" = 17), name = NULL) +
      (if (showleg) guides(colour = guide_legend(override.aes = list(shape = 16)),
                           shape  = guide_legend(override.aes = list(colour = "grey20")))
       else guides(colour = "none", shape = "none")) +
      mB_xscale(co) + labs(x = "chromosome", y = sprintf("eMLG BF (dB) - %s", lab)) + mB_theme
    if (!bottom) p + strip_x else p
  }
  ## c,d -- all-SNP, FDR floor-survivor cluster members coloured
  fpal <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#1B9E77","#D95F02","#7570B3")
  fdr_p <- function(bfcol, fg_groups, lab, bottom = FALSE) {
    bfk <- setNames(snp[[bfcol]], snp$marker)
    bg  <- data.table(x = co$xkey[snp$marker], BF = snp[[bfcol]])
    fg  <- m2g[group_id %in% fg_groups, .(marker, group_id)][
             , `:=`(x = co$xkey[marker], BF = bfk[marker], grp = factor(group_id, levels = fg_groups))]
    p <- ggplot() + shade_l("grey85") +
      geom_point(data = bg, aes(x, BF), colour = "grey60", size = 0.2, alpha = 0.5) +
      geom_point(data = fg, aes(x, BF, colour = grp), size = 0.9) +
      hthr + scale_colour_manual(values = setNames(fpal[seq_along(fg_groups)], fg_groups), guide = "none") +
      mB_xscale(co) + labs(x = "chromosome", y = sprintf("SNP BF (dB) - %s", lab)) + mB_theme
    if (!bottom) p + strip_x else p
  }
  ## e,f -- eMLG-filtered SNP, clusters of >=5 loci BF>=15 coloured
  gp <- snp[he[, .(group_id, marker = mB_bestmk[group_id])], on = "marker"][, .(group_id, Chr, Pos)]
  gp[, chr_num := suppressWarnings(as.integer(gsub("[^0-9]", "", Chr)))]; setorder(gp, chr_num, Pos)
  gp[, col := mB_clpal[(seq_len(.N) - 1) %% length(mB_clpal) + 1]]
  emlg <- snp[, .(marker, x = co$xkey[marker], BF1, BF2)][gp[, .(group_id, col)][m2g, on = "group_id"], on = "marker"]
  clu_p <- function(bfcol, lab, bottom = FALSE) {
    d   <- emlg[, .(group_id, col, x, BF = get(bfcol))]
    sig <- d[, .(n = sum(BF >= MB_SIG, na.rm = TRUE)), by = group_id][n >= MB_MIN5, group_id]
    p <- ggplot() + shade_l("grey85") +
      geom_point(data = d[!group_id %in% sig], aes(x, BF), size = 0.3, colour = "grey60", alpha = 0.5) +
      geom_point(data = d[group_id %in% sig][order(x)], aes(x, BF, colour = col), size = 0.5, alpha = 0.9) +
      scale_colour_identity() + hthr + mB_xscale(co) +
      labs(x = "chromosome", y = sprintf("SNP BF (dB) - %s", lab)) + mB_theme
    if (!bottom) p + strip_x else p
  }
  ## SNP-level analyses first (a-d), eMLG-level last (e,f); bottom panel keeps the x-axis
  panels <- list(clu_p("BF1","PC1"), clu_p("BF2","PC2"),
                 fdr_p("BF1", dt[floor1 == TRUE, group_id], "PC1"), fdr_p("BF2", dt[floor2 == TRUE, group_id], "PC2"),
                 emlg_p("eBF1","cat1","floor1","PC1"), emlg_p("eBF2","cat2","floor2","PC2", bottom = TRUE, showleg = FALSE))
  comb <- wrap_plots(panels, ncol = 1) + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
  save_fig(comb, "moduleB_manhattans", width = 400, height = 330, dpi = 150)
}

## =====================================================================
## MODULE C — climate-association strength vs the Omega-structured null
## =====================================================================
## Genome-wide test of whether the STRENGTH of climate association (BayPass BF
## on PC1/PC2) is organised by ancestry sorting, DI, or recombination beyond
## what population structure alone produces. Both figures read the single saved
## results object; observed climate axes are drawn over the 10,000-covariate
## Omega-structured null. Titles/subtitles are dropped (LaTeX caption supplies
## the "shaded = null 95% interval; dashed = null median" note).
MC_RES   <- "moduleC_climate_vs_sorting/data/moduleC_results.rds"
mC_pcpal <- c(PC1 = "#0072B2", PC2 = "#D55E00")   # observed climate axes
## primary genome-wide statistics, with underscore-free plain labels
mC_LAB <- c(rho_DI                  = "DI (Spearman correlation)",
            rho_rec                 = "recombination (Spearman correlation)",
            sort_gap_differentiated = "directional sorting (BF percentile gap)",
            rho_sort_magnitude      = "sorting magnitude (Spearman correlation)")

## ---- [Fig] structured-null calibration of the primary genome-wide stats --
fig_moduleC_calibration <- function() {
  r <- readRDS(MC_RES); obs <- r$observed; null <- r$null_stats
  fig_stats <- names(mC_LAB)
  long <- rbindlist(lapply(fig_stats, function(s)
    data.table(label = mC_LAB[[s]], value = null[, s])))
  long[, label := factor(label, levels = mC_LAB[fig_stats])]
  mk <- rbindlist(lapply(fig_stats, function(s) data.table(
    label = mC_LAB[[s]], median = median(null[, s]),
    lo = quantile(null[, s], 0.025), hi = quantile(null[, s], 0.975),
    PC1 = obs["PC1", s], PC2 = obs["PC2", s])))
  mk[, label := factor(label, levels = mC_LAB[fig_stats])]
  p <- ggplot(long, aes(value)) +
    geom_histogram(bins = 60, fill = "grey80", colour = NA) +
    geom_rect(data = mk, inherit.aes = FALSE,
              aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf), fill = "#999999", alpha = 0.18) +
    geom_vline(data = mk, aes(xintercept = median), linetype = 2, colour = "grey40") +
    geom_vline(data = mk, aes(xintercept = PC1, colour = "PC1"), linewidth = 0.7) +
    geom_vline(data = mk, aes(xintercept = PC2, colour = "PC2"), linewidth = 0.7) +
    facet_wrap(~ label, scales = "free", ncol = 2) +
    scale_colour_manual(values = mC_pcpal, name = NULL) +
    labs(x = "genome-wide statistic", y = "count (of 10,000 null covariates)") +
    theme_plain(9) + theme(legend.position = "bottom")
  save_fig(p, "moduleC_null_calibration", width = 175, height = 140, pdf = TRUE)
}

## ---- [Fig S] rank-threshold sensitivity (top 0.1% / 0.5% / 1%) -----------
fig_moduleC_threshold_sensitivity <- function() {
  r <- readRDS(MC_RES); null <- r$null_stats; th <- as.data.table(r$threshold)
  lev  <- c("top 0.1%", "top 0.5%", "top 1%")
  mlev <- c("median DI", "median recomb", "prop. differentiated", "prop. directional | diff.")
  setfac <- function(d) { d[, frac := factor(frac, levels = lev)]
                          d[, metric := factor(metric, levels = mlev)]; d }
  th_meta <- unique(th[, .(stat, frac, metric)])
  thl <- rbindlist(lapply(th_meta$stat, function(s) data.table(stat = s, value = null[, s])))
  thl <- setfac(th_meta[thl, on = "stat"])
  th_stat_meta <- unique(th[, .(stat, frac, metric, median = null_median, lo = null_lo, hi = null_hi)])
  th_obs <- dcast(th, stat ~ axis, value.var = "observed")
  thm <- setfac(th_obs[th_stat_meta, on = "stat"])
  p <- ggplot(thl, aes(value)) +
    geom_histogram(bins = 45, fill = "grey80", colour = NA) +
    geom_rect(data = thm, inherit.aes = FALSE,
              aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf), fill = "#999999", alpha = 0.18) +
    geom_vline(data = thm, aes(xintercept = median), linetype = 2, colour = "grey40") +
    geom_vline(data = thm, aes(xintercept = PC1, colour = "PC1"), linewidth = 0.6) +
    geom_vline(data = thm, aes(xintercept = PC2, colour = "PC2"), linewidth = 0.6) +
    facet_wrap(vars(metric, frac), scales = "free", ncol = 3,
               labeller = label_wrap_gen(multi_line = FALSE)) +
    scale_colour_manual(values = mC_pcpal, name = NULL) +
    labs(x = "value within the top-ranked fraction of eMLGs", y = "count (null covariates)") +
    theme_plain(8) + theme(legend.position = "bottom")
  save_fig(p, "moduleC_threshold_sensitivity", width = 185, height = 165, pdf = TRUE)
}

## =====================================================================
## build all
## =====================================================================
cat("== Module 0 ==\n")
fig_roc(); fig_ld_tracks(); fig_fidelity()
cat("== Module A ==\n")
fig_sorting_sweep(); fig_panelB_sweep(); fig_direction_sweep()
fig_architecture(); fig_manhattan_directional()
cat("== Module B ==\n")
fig_moduleB_ancestry(); fig_moduleB_manhattans()
## individual manhattans (fig_moduleB_eMLG_manhattan / _fdr_snp_manhattan /
## _snp_manhattan_clustered) remain defined above; call them if you want them split.
cat("== Module C ==\n")
fig_moduleC_calibration(); fig_moduleC_threshold_sensitivity()
cat("== Tables ==\n")
make_arch_table(); make_sorting_level_table(); make_models_table()
## make_exposure_table() -- admixture-exposure is an internal check, left OUT of the
## manuscript for now (2026-08); re-enable this call to rebuild the table if it returns.
cat("\nAll figures and tables written.\n")
