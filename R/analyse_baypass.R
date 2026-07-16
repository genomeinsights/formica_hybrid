library(data.table)
library(ggplot2)

# ------------------------------------------------------------
# Analysing BayPass results
# ------------------------------------------------------------
# Moved out of ld_pruning_DIEM.R (formerly baypass.R), which is now scoped
# to LD-pruning/eMLG generation only. Input-file generation for BayPass
# itself lives in write_baypass_inputs.R / prepare_with_aland.R /
# prepare_aland_excluded.R; the actual BayPass runs happen via
# with_aland/run_baypass.sh and aland_excluded/run_baypass.sh, each
# producing PC1_DIEM_noOmega_*, PC2_DIEM_noOmega_*, PC1_DIEM_withOmega_*,
# PC2_DIEM_withOmega_* summary files. This file analyses whichever of
# those result sets you point it at.

message("=== Loading data ===")
load("./data/hybrids_only_maf005.Rdata")

## --- import a BayPass covariate association result into map_hyb_005 ---
##
## summary_betai_reg.out has no marker ID column at all -- MRK is just
## BayPass's row index into whatever -countdatafile it was given (the FULL
## unpruned set). This is therefore a POSITIONAL join (MRK==row number),
## not a name-based one, and only holds because u_DIEM.geno was written
## from GTs_hybrids_005 in map_hyb_005's row order and never reordered
## afterward -- if that ever changes, this breaks silently (wrong BF
## values on wrong markers, no error thrown). The stopifnot below is the
## guard against exactly that (confirmed to catch a real stale-data
## mismatch directly: 1,122,864 vs 1,114,340 markers after the DIEM
## bi-allelic-filter rerun).
import_baypass_covariate <- function(file, n_expected) {
  res <- fread(file)
  stopifnot(
    "summary_betai_reg.out row count doesn't match map_hyb_005 -- MRK-based positional join would silently misalign" =
      nrow(res) == n_expected,
    "MRK isn't 1:N in order -- positional join assumption doesn't hold for this file" =
      identical(res$MRK, seq_len(nrow(res)))
  )
  data.table(BF = res$`BF(dB)`, eBPis = res$eBPis, Beta = res$Beta_is)
}

## --- Manhattan plot: cumulative bp on x-axis, spacer between
## chromosomes, alternating background shading every second chromosome ---
##
## chr_order defaults to natural numeric sort (Chr1, Chr2, ..., Chr10, ...)
## rather than lexicographic (Chr1, Chr10, Chr11, ..., Chr2, ...), which is
## what a plain factor()/sort() would give since Chr names are strings.
manhattan_plot <- function(dt, value_col, chr_col = "Chr", pos_col = "Pos",
                            chr_order = NULL, spacer_frac = 0.01,
                            threshold = NULL, point_color = "grey20",
                            shade_color = "grey85", point_size = 0.5,
                            title = NULL, ylab = value_col) {

  dt <- data.table::copy(dt)
  data.table::setnames(dt, c(chr_col, pos_col, value_col), c("Chr", "Pos", "value"))
  dt <- dt[!is.na(value) & !is.na(Pos)]

  if (is.null(chr_order)) {
    u <- unique(dt$Chr)
    chr_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", u)))
    chr_order <- u[order(chr_num, u)]
  }
  dt <- dt[Chr %in% chr_order]
  dt[, Chr := factor(Chr, levels = chr_order)]
  data.table::setorder(dt, Chr, Pos)

  chr_len <- dt[, .(len = max(Pos)), by = Chr]
  chr_len <- chr_len[match(chr_order, Chr)]
  spacer <- spacer_frac * sum(chr_len$len)
  chr_len[, chr_start := c(0, head(cumsum(len + spacer), -1))]

  dt <- chr_len[, .(Chr, chr_start)][dt, on = "Chr"]
  dt[, x := Pos + chr_start]

  chr_mid <- chr_len[, .(Chr, mid = chr_start + len / 2)]

  shade_chrs <- chr_order[seq(2, length(chr_order), by = 2)]
  shade_rects <- chr_len[Chr %in% shade_chrs, .(xmin = chr_start, xmax = chr_start + len)]

  p <- ggplot(dt, aes(x, value)) +
    geom_point(size = point_size, color = point_color) +
    scale_x_continuous(breaks = chr_mid$mid, labels = gsub("^Chr", "", chr_mid$Chr), expand = c(0.01, 0.01)) +
    labs(x = "Chromosome", y = ylab, title = title) +
    theme_bw(base_size = 13) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

  if (nrow(shade_rects)) {
    p <- p + geom_rect(
      data = shade_rects, aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf,
      inherit.aes = FALSE, fill = shade_color, alpha = 0.5
    )
    ## shading was added as the LAST layer (drawn on top) -- move it first
    ## so points render over the shading, not under it
    p$layers <- rev(p$layers)
  }

  if (!is.null(threshold)) {
    p <- p + geom_hline(yintercept = threshold, linetype = 2, color = "red", linewidth = 0.4)
  }
  p
}

## --- import + plot one PC's association result, for a given results
## folder/outprefix combination (out_folder/{prefix}_summary_betai_reg.out)
## -- call this once per (population set) x (PC) x (with/without omega)
## combination you want a Manhattan plot for.
analyse_baypass_pc <- function(out_folder, prefix, map, fig_prefix, ylab = "BF (dB)") {
  res <- import_baypass_covariate(
    paste0(out_folder, prefix, "_summary_betai_reg.out"), nrow(map)
  )
  map <- copy(map)
  map[, BF := res$BF]
  map[, eBPis := res$eBPis]

  p <- manhattan_plot(
    map, value_col = "BF", ylab = ylab,
    title = paste0(prefix, " association -- BayPass BF(dB)")
  )
  ggsave(paste0("./Figures/manhattan_", fig_prefix, ".png"), p, width = 16, height = 5, dpi = 150)

  list(map_snp = map, plot = p)
}

## --- example usage ---
## out_baypass_2/ is what's currently available (pre-Aland-exclusion,
## pre-poolsize-fix -- see project memory / earlier conversation for why
## these are stale). Once with_aland/run_baypass.sh and
## aland_excluded/run_baypass.sh have actually been run, point out_folder
## at those instead, e.g.:
##   analyse_baypass_pc("./with_aland/", "PC1_DIEM_noOmega", map_hyb_005, "PC1_with_aland_noOmega")
##   analyse_baypass_pc("./aland_excluded/", "PC1_DIEM_withOmega", map_hyb_005, "PC1_aland_excluded_withOmega")
## and so on for PC2 / the other omega setting.
pc1 <- analyse_baypass_pc("./out_baypass_2/", "PC1_DIEM", map_hyb_005, "PC1_BF")
pc2 <- analyse_baypass_pc("./out_baypass_2/", "PC2_DIEM", map_hyb_005, "PC2_BF")
