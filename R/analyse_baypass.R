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
#
# Only output kept: one Manhattan plot per PC, background (all
# eMLG-filtered markers, grey) + foreground (clusters with n_loci >= 10,
# colored by LD cluster) in a single plot. Earlier exploratory variants
# (plain Manhattan, noOmega-vs-withOmega comparison, has_eMLG-only,
# min_n_sig_loci-based filtering) were dropped -- see git history if
# any of those are wanted again.

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

load_baypass_bf <- function(out_folder, prefix, map) {
  res <- import_baypass_covariate(paste0(out_folder, prefix, "_summary_betai_reg.out"), nrow(map))
  map <- copy(map)
  map[, BF := res$BF]
  map[, eBPis := res$eBPis]
  list(map_snp = map)
}

## with_aland/ (all populations, including Aland) and aland_excluded/
## (Aland's samples excluded) results, withOmega setting -- each
## -countdatafile u_DIEM.geno (full, unpruned set, matching map_hyb_005's
## row count -- verified directly: 1,114,340 rows in every case). Both
## population sets share the SAME marker/eMLG-cluster structure below,
## since that was computed once genome-wide from the full sample set,
## independent of which BayPass run's poolsizes are used. BF(dB) >= 20
## ("strong" evidence on Jeffreys' scale) is the significance threshold.
sig_threshold <- 20
pc1_with_aland_withOmega     <- load_baypass_bf("./with_aland/",     "PC1_DIEM_withOmega", map_hyb_005)
pc2_with_aland_withOmega     <- load_baypass_bf("./with_aland/",     "PC2_DIEM_withOmega", map_hyb_005)
pc1_aland_excluded_withOmega <- load_baypass_bf("./aland_excluded/", "PC1_DIEM_withOmega", map_hyb_005)
pc2_aland_excluded_withOmega <- load_baypass_bf("./aland_excluded/", "PC2_DIEM_withOmega", map_hyb_005)

## noOmega: same PC association runs, but WITHOUT the Omega pre-estimated
## from the pruned genotype set (step 1 of run_baypass.sh) -- BayPass
## estimates its own Omega internally from the full, unpruned genotype set
## during the run instead. Comparing against withOmega checks how much the
## result depends on which Omega estimate is used.
pc1_with_aland_noOmega     <- load_baypass_bf("./with_aland/",     "PC1_DIEM_noOmega", map_hyb_005)
pc2_with_aland_noOmega     <- load_baypass_bf("./with_aland/",     "PC2_DIEM_noOmega", map_hyb_005)
pc1_aland_excluded_noOmega <- load_baypass_bf("./aland_excluded/", "PC1_DIEM_noOmega", map_hyb_005)
pc2_aland_excluded_noOmega <- load_baypass_bf("./aland_excluded/", "PC2_DIEM_noOmega", map_hyb_005)

## --- eMLG-filtered marker set ---
##
## Restricts to markers belonging to a cluster that actually got an eMLG
## (has_eMLG==TRUE, i.e. n_loci >= min_n_loci_eMLG under this run's
## settings -- see ./data/eMLG_5loci_0025_cM05.rds, ld_w_threshold=0.025,
## min_n_loci_flag=5, min_n_loci_eMLG=5).
eMLG_result <- readRDS("./data/eMLG_5loci_0025_cM05.rds")
has_eMLG_groups <- eMLG_result$groups[has_eMLG == TRUE]
marker_group <- has_eMLG_groups[TRUE, .(marker = unlist(members), n_loci = n_loci), by = group_id]
message(
  "eMLG-filtered marker set: ", nrow(marker_group), "/", nrow(map_hyb_005), " markers, ",
  uniqueN(marker_group$group_id), " clusters"
)

## --- global, stable color assignment per LD cluster ---
##
## Clusters are the SAME objects across every PC/population-set figure
## (built once, genome-wide, independent of any particular BayPass run --
## only WHICH clusters pass min_n_sig_loci differs per plot), so a color
## assigned per-figure (by that figure's own first-appearance order) gives
## the same physical region a different color in each plot purely by
## chance. Assigning color ONCE here, in genome position order, means the
## same cluster is always the same color -- e.g. a Chr24 region flagged as
## an outlier in both PC1 and PC2 will visibly match between figures
## instead of only being findable by comparing chromosome/position by eye.
pal_cluster <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33",
                  "#A65628","#F781BF","#1B9E77","#D95F02","#7570B3","#66A61E")
group_pos <- map_hyb_005[, .(marker, Chr, Pos)][
  has_eMLG_groups[TRUE, .(group_id, representative)], on = c(marker = "representative")
]
group_pos[TRUE, chr_num := suppressWarnings(as.integer(gsub("[^0-9]", "", Chr)))]
data.table::setorder(group_pos, chr_num, Pos)
group_pos[TRUE, group_color := pal_cluster[(seq_len(.N) - 1) %% length(pal_cluster) + 1]]
marker_group <- group_pos[TRUE, .(group_id, group_color)][marker_group, on = "group_id"]

## --- background (all eMLG-filtered markers, light grey) + foreground
## (clusters with >= min_n_sig_loci INDIVIDUALLY significant members,
## colored, drawn on top) in one plot, so the amount of noise removed is
## directly visible. has_eMLG==TRUE alone only guarantees a cluster is
## big enough (n_loci>=5), not that it's actually associated -- filtering
## on cluster SIZE alone (tried directly: n_loci>=10 regardless of BF)
## doesn't declutter at all, since it's independent of association
## strength (confirmed: identical background/foreground counts for PC1
## and PC2 that way). What matters is how many members are INDIVIDUALLY
## significant: a cluster where only a handful of loci cross the
## threshold within an otherwise-unremarkable LD-connected block is more
## likely a false positive than a real signal (shown directly in the
## LDscnR paper) -- min_n_sig_loci requires the signal to be broadly
## supported across the cluster, not just a few individually-lucky
## markers.
##
## Both layers share ONE coordinate transform (computed from the full,
## unfiltered marker set) so background and foreground points align
## correctly -- computing chr_start independently per layer risks
## misalignment if the two subsets don't span identical per-chromosome
## position ranges.
plot_clustered_manhattan_with_background <- function(baypass_result, fig_prefix, ylab = "BF (dB)",
                                                       threshold = sig_threshold, min_n_sig_loci = 10,
                                                       spacer_frac = 0.01) {
  full <- baypass_result$map_snp[, .(marker, Chr, Pos, BF)]

  chr_order <- {
    u <- unique(full$Chr)
    chr_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", u)))
    u[order(chr_num, u)]
  }
  chr_len <- full[, .(len = max(Pos)), by = Chr][match(chr_order, Chr)]
  spacer <- spacer_frac * sum(chr_len$len)
  chr_len[, chr_start := c(0, head(cumsum(len + spacer), -1))]
  chr_mid <- chr_len[, .(Chr, mid = chr_start + len / 2)]
  shade_chrs <- chr_order[seq(2, length(chr_order), by = 2)]
  shade_rects <- chr_len[Chr %in% shade_chrs, .(xmin = chr_start, xmax = chr_start + len)]
  add_x <- function(d) {
    d <- chr_len[, .(Chr, chr_start)][d, on = "Chr"]
    d[, x := Pos + chr_start]
    d
  }

  eMLG_dt <- full[marker_group, on = "marker"]
  sig_ids <- eMLG_dt[TRUE, .(n_sig = sum(BF >= threshold, na.rm = TRUE)), by = group_id][n_sig >= min_n_sig_loci, group_id]

  ## bg is exactly the markers EXCLUDED by min_n_sig_loci (not the full
  ## has_eMLG superset) -- keeps the grey/colored split exact rather than
  ## relying on draw-order overplotting to hide the difference
  bg <- add_x(eMLG_dt[!group_id %in% sig_ids, .(marker, Chr, Pos, BF)])
  ## group_color already carried through from marker_group's global,
  ## position-ordered assignment (see above) -- NOT recomputed per-call,
  ## so the same cluster is always the same color across every figure
  fg <- add_x(eMLG_dt[group_id %in% sig_ids])
  data.table::setorder(fg, x)

  p <- ggplot() +
    { if (nrow(shade_rects)) geom_rect(
        data = shade_rects, aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf,
        fill = "grey85", alpha = 0.5
      ) } +
    geom_point(data = bg, aes(x, BF), size = 0.5, color = "grey60", alpha = 0.5) +
    geom_point(data = fg, aes(x, BF, color = group_color), size = 0.7, alpha = 0.9) +
    scale_color_identity() +
    geom_hline(yintercept = threshold, linetype = 2, color = "red", linewidth = 0.4) +
    scale_x_continuous(breaks = chr_mid$mid, labels = gsub("^Chr", "", chr_mid$Chr), expand = c(0.01, 0.01)) +
    labs(
      x = "Chromosome", y = ylab,
      title = sprintf(
        "%s -- grey: eMLG-filtered background (%s markers); color: clusters with >=%d loci BF>=%s (%s markers, %d clusters)",
        fig_prefix, format(nrow(bg), big.mark = ","), min_n_sig_loci, threshold,
        format(nrow(fg), big.mark = ","), uniqueN(fg$group_id)
      )
    ) +
    theme_bw(base_size = 13) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

  fname <- paste0(
    "./Figures/manhattan_", fig_prefix, "_eMLG_clustered_min", min_n_sig_loci, "SigLoci_withBackground.png"
  )
  ggsave(fname, p, width = 16, height = 5, dpi = 150)
  message(
    fig_prefix, ": background ", nrow(bg), " markers, foreground ", nrow(fg),
    " markers (", uniqueN(fg$group_id), " clusters) -> ", fname
  )
  list(background = bg, foreground = fg, plot = p)
}

pc1_with_aland_withOmega_layered     <- plot_clustered_manhattan_with_background(pc1_with_aland_withOmega,     "PC1_with_aland_withOmega",     min_n_sig_loci = 10)
pc2_with_aland_withOmega_layered     <- plot_clustered_manhattan_with_background(pc2_with_aland_withOmega,     "PC2_with_aland_withOmega",     min_n_sig_loci = 10)
pc1_aland_excluded_withOmega_layered <- plot_clustered_manhattan_with_background(pc1_aland_excluded_withOmega, "PC1_aland_excluded_withOmega", min_n_sig_loci = 10)
pc2_aland_excluded_withOmega_layered <- plot_clustered_manhattan_with_background(pc2_aland_excluded_withOmega, "PC2_aland_excluded_withOmega", min_n_sig_loci = 10)

pc1_with_aland_noOmega_layered     <- plot_clustered_manhattan_with_background(pc1_with_aland_noOmega,     "PC1_with_aland_noOmega",     min_n_sig_loci = 10)
pc2_with_aland_noOmega_layered     <- plot_clustered_manhattan_with_background(pc2_with_aland_noOmega,     "PC2_with_aland_noOmega",     min_n_sig_loci = 10)
pc1_aland_excluded_noOmega_layered <- plot_clustered_manhattan_with_background(pc1_aland_excluded_noOmega, "PC1_aland_excluded_noOmega", min_n_sig_loci = 10)
pc2_aland_excluded_noOmega_layered <- plot_clustered_manhattan_with_background(pc2_aland_excluded_noOmega, "PC2_aland_excluded_noOmega", min_n_sig_loci = 10)


