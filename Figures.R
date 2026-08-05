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
  sw[, prop_unres := bidirectional / sorted]
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
  row <- function(label, est, stat, sy)
    sprintf("\\quad %s & $%s$ & $%s=%s$ \\\\", label,
            formatC(est, format = "f", digits = 3), sy, formatC(stat, format = "f", digits = 1))
  tex <- c(
    "\\begin{table}[h]", "\\centering",
    "\\caption{\\textbf{Architecture regression models of ancestry sorting among",
    "LD-reduced units.} Standardised coefficients. \\emph{Magnitude}: linear model of",
    "the sorted fraction $p_{\\mathrm{fixed}}$ on recombination and DI.",
    "\\emph{Direction}: logistic model of the probability of fixing toward \\Faq{} among",
    "directionally resolved units, on DI, recombination, parental MAF and log cluster",
    "size. All predictors standardised over the differentiated-unit set; $t$/$z$ are",
    "Wald statistics. All terms significant (weakest: parental MAF, $P=0.009$).}",
    "\\label{tab:architecture-models}",
    "\\begin{tabular}{lrr}", "\\toprule",
    "Predictor & Coefficient & Statistic \\\\", "\\midrule",
    "\\multicolumn{3}{l}{\\emph{Magnitude} ($p_{\\mathrm{fixed}}$, linear)} \\\\",
    row("Recombination", mag["zr", "Estimate"],  mag["zr", "t value"],  "t"),
    row("DI",            mag["zDI", "Estimate"], mag["zDI", "t value"], "t"),
    "\\addlinespace",
    "\\multicolumn{3}{l}{\\emph{Direction} ($P(\\text{fix}\\rightarrow\\Faq)$, logistic)} \\\\",
    row("DI",               dr["zDI", "Estimate"],  dr["zDI", "z value"],  "z"),
    row("Recombination",    dr["zr", "Estimate"],   dr["zr", "z value"],   "z"),
    row("Parental MAF",     dr["zmaf", "Estimate"], dr["zmaf", "z value"], "z"),
    row("Log cluster size", dr["zcs", "Estimate"],  dr["zcs", "z value"],  "z"),
    "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex, "Tables/moduleA_architecture_models.tex")
  cat("Saved: Tables/moduleA_architecture_models.tex\n")
}

## =====================================================================
## build all
## =====================================================================
cat("== Module 0 ==\n")
fig_roc(); fig_ld_tracks(); fig_fidelity()
cat("== Module A ==\n")
fig_sorting_sweep(); fig_panelB_sweep(); fig_direction_sweep()
fig_architecture(); fig_manhattan_directional()
cat("== Tables ==\n")
make_arch_table(); make_sorting_level_table(); make_models_table()
cat("\nAll figures and tables written.\n")
