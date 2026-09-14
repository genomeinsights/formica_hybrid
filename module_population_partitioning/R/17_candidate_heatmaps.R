## =========================================================================
## module_population_partitioning -- candidate-locus follow-up 17: item 1,
## locus x population heatmap of candidate best-SNP allele frequencies.
##
## Frequencies are STANDARDIZED WITHIN EACH LOCUS (z-score across the 19
## hybrid populations) for visualization only -- unstandardized (raw,
## oriented [0,1]) values are retained in the saved data for analysis.
## Every locus is annotated by: BayPass statistic, raw-candidate status,
## floor-survivor status, chromosome, position, DI, recombination rate,
## sorting direction, and BDMI overlap.
##
## Loops over the 4 targets (PC1, PC2, bio_winter, mitoC2) x {raw, floor,
## top-10} status, skipping a status with zero loci (e.g. mitoC2's floor
## set, which is empty by design -- see 15_candidate_import.R). PC1/PC2 raw
## sets and the mitoC2 set are never labelled "discovery"/"discoveries" in
## any generated text; bio_winter's floor set is the only one highlighted
## as such, per TARGETS$<tgt>$is_discovery/no_discovery_set from script 15.
##
## Run from the formica_hybrid repo root, after 15_candidate_import.R:
##   Rscript module_population_partitioning/R/17_candidate_heatmaps.R
## Reads : module_population_partitioning/data/followup/15_candidate_data.rds
## Writes: module_population_partitioning/data/followup/17_heatmap_data.rds
##         module_population_partitioning/Figures/followup/17_heatmap_<target>_<status>.png
## =========================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
t_start <- Sys.time()
OUTDIR <- "module_population_partitioning/data/followup"
FIGDIR <- "module_population_partitioning/Figures/followup"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
r15 <- readRDS(file.path(OUTDIR, "15_candidate_data.rds"))
u <- r15$u; Fmat_all <- r15$Fmat_all
stopifnot("candidate unit table must have exactly 18,361 rows -- rerun 15_candidate_import.R" = nrow(u) == 18361L,
         "Fmat_all columns must match u$group_id in order" = identical(colnames(Fmat_all), u$group_id))

AQU <- "#21918C"; POL <- "#D3C93B"
theme_ms <- theme_bw(base_size = 11) + theme(strip.background = element_blank(), panel.grid.minor = element_blank())
STATUS_LABEL <- c(raw = "raw candidates", floor = "floor survivors", topN = sprintf("top-%d ranked", r15$N_TOP))

## ---------------------------------------------------------------------
## build one heatmap (+ annotation strips) for one target/status
## ---------------------------------------------------------------------
build_heatmap <- function(tgt, status, ids) {
  cfg <- r15$TARGETS[[tgt]]
  sub <- u[match(ids, group_id)]
  setorder(sub, ChrNum, Pos)
  ord_ids <- sub$group_id

  Fsub <- Fmat_all[, ord_ids, drop = FALSE]           # pops x n_loci, RAW (unstandardized) oriented frequency
  Zsub <- scale(Fsub, center = TRUE, scale = TRUE)     # per-LOCUS (column) z-score across the 19 pops -- viz only

  dm_raw <- as.data.table(as.table(Fsub)); setnames(dm_raw, c("pop", "group_id", "f_aqu_raw"))
  dm_z   <- as.data.table(as.table(Zsub)); setnames(dm_z,   c("pop", "group_id", "f_aqu_z"))
  dm <- merge(dm_raw, dm_z, by = c("pop", "group_id"))
  dm[, group_id := factor(group_id, levels = ord_ids)]

  ann <- sub[, .(group_id, Chr, Pos, DI, recomb, sort_class, bdmi_overlap,
                stat = get(cfg$stat_col), is_floor = get(cfg$floor_col))]
  ann[, group_id := factor(group_id, levels = ord_ids)]
  ann[, chr_stripe := as.integer(sub("Chr", "", Chr)) %% 2L]

  p_stat <- ggplot(ann, aes(group_id, 1, fill = stat)) + geom_tile() +
    scale_fill_gradient(low = "grey90", high = "firebrick", name = cfg$stat_col) +
    theme_void() + theme(legend.position = "top", legend.key.height = unit(3, "mm"))
  p_floor <- ggplot(ann, aes(group_id, 1, fill = is_floor)) + geom_tile() +
    scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "grey90"), name = "floor survivor") +
    theme_void() + theme(legend.position = "top", legend.key.height = unit(3, "mm"))
  p_bdmi <- ggplot(ann, aes(group_id, 1, fill = bdmi_overlap)) + geom_tile() +
    scale_fill_manual(values = c(`TRUE` = "#7570b3", `FALSE` = "grey90"), name = "BDMI overlap") +
    theme_void() + theme(legend.position = "top", legend.key.height = unit(3, "mm"))
  p_sort <- ggplot(ann, aes(group_id, 1, fill = sort_class)) + geom_tile() +
    scale_fill_manual(values = c(aquilonia = AQU, polyctena = POL, unresolved = "grey50",
                                 unsorted = "grey90", ambiguous = "grey70"), na.value = "white", name = "sort_class") +
    theme_void() + theme(legend.position = "top", legend.key.height = unit(3, "mm"))
  p_chr <- ggplot(ann, aes(group_id, 1, fill = chr_stripe)) + geom_tile() +
    scale_fill_gradient(low = "grey85", high = "grey55", guide = "none") +
    theme_void()

  p_hm <- ggplot(dm, aes(group_id, pop, fill = f_aqu_z)) + geom_tile() +
    scale_fill_gradient2(low = POL, mid = "grey95", high = AQU, midpoint = 0, name = "standardized\naquilonia freq.\n(z, within-locus)") +
    labs(x = sprintf("%s: %s (%s, n=%d loci, chromosome-ordered)", tgt, STATUS_LABEL[status],
                     if (cfg$is_discovery) "highlighted discovery set" else if (isTRUE(cfg$no_discovery_set) && status != "raw") "NOT a discovery set" else "not described as a discovery",
                     length(ord_ids)),
        y = NULL) +
    theme_ms + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())

  fig <- p_chr + p_stat + p_floor + p_bdmi + p_sort + p_hm +
    plot_layout(ncol = 1, heights = c(0.4, 1, 1, 1, 1, 12))
  fname <- sprintf("17_heatmap_%s_%s.png", tgt, status)
  ggsave(file.path(FIGDIR, fname), fig, width = max(7, min(20, length(ord_ids) * 0.15 + 4)), height = 8, dpi = 200, limitsize = FALSE)
  cat(sprintf("[heatmap] %-10s %-8s: %3d loci -> %s\n", tgt, status, length(ord_ids), fname))

  list(target = tgt, status = status, ids = ord_ids, dm = dm, ann = ann, figure = fname)
}

## ---------------------------------------------------------------------
## loop over targets x status, skipping empty sets
## ---------------------------------------------------------------------
heatmap_data <- list()
for (tgt in names(r15$TARGETS)) {
  sets <- list(raw = r15$cand_ids[[tgt]], floor = r15$floor_ids[[tgt]], topN = r15$topN_ids[[tgt]])
  for (status in names(sets)) {
    ids <- sets[[status]]
    if (!length(ids)) { cat(sprintf("[heatmap] %-10s %-8s: 0 loci, skipped\n", tgt, status)); next }
    heatmap_data[[paste(tgt, status, sep = "_")]] <- build_heatmap(tgt, status, ids)
  }
}

result <- list(heatmap_data = heatmap_data, STATUS_LABEL = STATUS_LABEL,
              session_info = sessionInfo(), run_time = Sys.time(),
              elapsed_secs = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
saveRDS(result, file.path(OUTDIR, "17_heatmap_data.rds"))
cat(sprintf("\n[heatmap] saved -> %s (elapsed %.1fs)\n", file.path(OUTDIR, "17_heatmap_data.rds"), result$elapsed_secs))
