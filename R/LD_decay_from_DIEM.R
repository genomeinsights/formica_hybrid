library(ggplot2)
library(data.table)
library(patchwork)
library(SNPRelate)
devtools::load_all("~/gitlab/LDscnR/")
# ------------------------------------------------------------
# Parse diem data to GTs/map 
# ------------------------------------------------------------
if(!file.exists("./data/diem_parsed.rds")){
  sample_data <- fread("data/Sample_covariate_info_outlier_analysis_20.txt")
  DIEM <- fread("data/Formica_hybrids_filtered_diem_output.bed.gz")
  DIEM_samples <- colnames(DIEM)[10]
  DIEM_samples <- strsplit(DIEM_samples,"|",fixed = TRUE)[[1]]
  
  DIEM <- DIEM[DIEM$nVNTs==2]
  ## map
  map <- DIEM[,.(Chr=gsub("chromosome_","Chr",`#Chrom`),Pos=End)]
  map[,marker := paste(Chr,Pos,sep=":")]
  map <- cbind(map,DIEM[,.(Polarity,DiagnosticIndex)])
  
  
  # Pre-extract only the column you need so workers do not carry full DIEM
  track_strings <- as.character(DIEM[[10]])
  
  #track_strings
  
  # Process a chunk of rows per worker, not one row per future
  idx_chunks <- split(seq_len(nrow(DIEM)),ceiling(seq_len(nrow(DIEM)) / 1000))
  
  pb <- txtProgressBar(min = 0, max = nrow(map)-1, style = 3)
  
  #idx <- idx_chunks[[1]]
  DIEM_data <- rbindlist(lapply(idx_chunks,function(idx){
    setTxtProgressBar(pb, max(idx))
    out <- vector("list", length(idx))
    #j <- 1
    for (j in seq_along(idx)) {
      x <- idx[j]
      
      parts <- strsplit(track_strings[x], "|", fixed = FALSE)[[1]]
      track <- suppressWarnings(as.numeric(parts[-1]))
      
      
      out[[j]] <- cbind(map[x], track=list(track))
    }
    rbindlist(out[!vapply(out, is.null, logical(1))])
  }))
  
  GTs <- do.call(cbind,DIEM_data$track)
  rownames(GTs) <- DIEM_samples
  colnames(GTs) <-  map$marker
  
  keep_inds <- c(sample_data[,Sample_ID],DIEM_samples[grepl("Fpol",DIEM_samples) | grepl("Faqu",DIEM_samples)])
  
  GTs <- GTs[keep_inds,]
 
  ## ploarize
  GTs[,map$Polarity==1] <- 2-GTs[,map$Polarity==1]
  
  saveRDS(list(GTs=GTs,map=map),"./data/diem_parsed.rds")
  
}else{
  tmp <- readRDS("./data/diem_parsed.rds")
  GTs <- tmp$GTs
  map <- tmp$map
  sample_data <- fread("data/Sample_covariate_info_outlier_analysis_20.txt")
  rm(tmp)
  gc()
}

# ------------------------------------------------------------
# compute LD-decay
# ------------------------------------------------------------
# keep only hybrids
GTs_hybrids <- GTs[sample_data$Sample_ID,]

## filter by maf easies through SNP relate
gds_hyb <- create_gds_from_geno(geno=GTs_hybrids, map, "gds_hybrids.gds") #wrapper from LDscnR
#maf <- snpgdsSNPRateFreq(gds_hyb)$MinorFreq
map[,maf_hyb:= snpgdsSNPRateFreq(gds_hyb)$MinorFreq]
map_hyb_005 <- map[maf_hyb>=0.05]
GTs_hybrids_005 <- GTs_hybrids[TRUE,map_hyb_005$marker]

snpgdsClose(gds_hyb)
rm(GTs,GTs_hybrids)
gc()

#save(GTs_hybrids_005,map_hyb_005,sample_data,file = "./data/hybrids_only_maf005.Rdata")

if(!file.exists("./data/ld_decay_DIEM_100w.rds")){
  
  gds_hyb <- create_gds_from_geno(geno=GTs_hybrids_005, map_hyb_005, "gds_hybrids.gds")
  
  ld_decay_DIEM_100w <- compute_LD_decay(
    gds_hyb,n_win_decay = 100,
    el_data_folder = "./el_diem/", # too large to keep in memory
    keep_el = TRUE,
    slide=400, 
    cores = 10,ld_method = "corr"
  )
  saveRDS(ld_decay_DIEM_100w,"./data/ld_decay_DIEM_100w.rds")
  
  # keep only diagnostic markers
  gds_hyb_DI <- create_gds_from_geno(geno=GTs_hybrids_005[TRUE,map_hyb_005[DiagnosticIndex > (-25),marker]], map_hyb_005[DiagnosticIndex > (-25)], "gds_hybridsID.gds")
  
  ld_decay_DI <- compute_LD_decay(
    gds_hyb_DI,
    keep_el = TRUE,
    slide=100, ## slide~100 is enough
    cores = 10,ld_method = "corr"
  )
  
  saveRDS(ld_decay_DI,"./data/datald_decay_DI.rds")
  
}else{
  ld_decay <- readRDS("./data/ld_decay_DIEM_100w.rds")  
  ld_decay_DI <- readRDS("./data/datald_decay_DI.rds")  
}

## compare
# plot(ld_decay_DI)
# plot(ld_decay_DI,type="chr",chr="Chr18")
# 
# plot(ld_decay_DIEM)
# plot(ld_decay_DIEM,type="chr",chr="Chr18")
# 
# ld_w_095_DI <- compute_ld_w(0.95,ld_decay = ld_decay_DI)
# ld_w_095_sim <- compute_ld_w(0.95,ld_decay = ld_decay_sim)

# ------------------------------------------------------------
# LD-decay and ld_w compared to recombination rate
# ------------------------------------------------------------
# work on full data set
ld_decay <- readRDS("./data/ld_decay_DIEM_100w.rds")


ld_w_095 <- as.vector(compute_ld_w(0.95,ld_decay = ld_decay))
map_hyb_005[,ld_w_095:=ld_w_095]
save(GTs_hybrids_005,map_hyb_005,sample_data,ld_decay,file = "./data/hybrids_only_maf005.Rdata")

map <- copy(map_hyb_005)

rec_map <- fread("./data/Frufa_DTOL_PR.ref_genome.recmap")
rec_map[, Chr := paste0("Chr", sub("chromosome_", "", chr))]

## rec_rate for the interval ENDING at each row must be computed within each
## chromosome (a plain diff() bleeds across chromosome boundaries) and mid
## must average this row's pos with the PREVIOUS row's pos (the interval the
## rate actually describes) -- using the next row's pos (shift(type="lead"))
## silently shifts every rec_rate one marker out of register with its true
## physical location.
rec_map[, rec_rate := c(NA_real_, diff(cM) / diff(pos)), by = Chr]
rec_map[TRUE, mid := (pos + shift(pos, type = "lag")) / 2, by = Chr]

plot_ld_decay_tracks <- function(
    ld_decay,
    ld_w_095 = NULL,
    rec_dt = NULL,
    map = NULL,
    outliers = NULL,
    plot_ld_w = TRUE,
    plot_rec_rate = FALSE,
    highlight_outliers = TRUE,
    a_threshold = 0.0025,
    ncol = 5
) {
  
  library(data.table)
  library(ggplot2)
  
  ## Build LD-decay window table
  plot_dt <- rbindlist(lapply(names(ld_decay$by_chr), function(ch) {
    
    chr_obj <- ld_decay$by_chr[[ch]]
    decay <- copy(chr_obj$decay)
    
    decay[, Chr := ch]
    decay[, mid := rowMeans(.SD, na.rm = TRUE),
          .SDcols = c("start", "end")]
    decay[, chr_mean := chr_obj$decay_sum$a[1]]
    
    decay
  }), fill = TRUE)
  
  plot_dt[, Chr := factor(Chr, levels = names(ld_decay$by_chr))]

  ## Windows with too few LD pairs are dropped entirely upstream (not kept as
  ## NA rows), so large assembly gaps / unmappable regions show up as a jump
  ## in `mid` between consecutive *existing* rows rather than an NA to break
  ## on. Flag those jumps here so every track's line can be split there
  ## instead of geom_line() silently bridging straight across the hole.
  setorder(plot_dt, Chr, start)
  plot_dt[, contig_run := {
    gaps <- diff(mid)
    brk <- c(TRUE, gaps > 1.5 * stats::median(gaps, na.rm = TRUE))
    cumsum(brk)
  }, by = Chr]

  ## Add median ld_w per window
  if (plot_ld_w && !is.null(map)) {

    ## built from map (which carries ld_w_095 positionally aligned to its own
    ## rows) rather than names(ld_w_095) -- as.vector() upstream strips names,
    ## so that round-trip silently produced an empty `snp` column.
    ld_dt <- map[, .(Chr, Pos, ld_w = ld_w_095)]

    plot_dt[, ld_w_med := NA_real_]
    
    for (ch in names(ld_decay$by_chr)) {
      
      ld_ch <- ld_dt[Chr == ch]
      
      plot_dt[Chr == ch, ld_w_med := sapply(seq_len(.N), function(i) {
        median(
          ld_ch[Pos >= start[i] & Pos <= end[i], ld_w],
          na.rm = TRUE
        )
      })]
    }
    
    ## scale ld_w to comparable axis as a
    plot_dt[, scale_fac_ld := max(a, na.rm = TRUE) /
              max(ld_w_med, na.rm = TRUE),
            by = Chr]
    
    plot_dt[, ld_w_scaled := ld_w_med * scale_fac_ld]

    ## break the line both at position gaps (contig_run) and at windows with
    ## too few SNPs for a median (NA) -- either should produce a real gap
    ## instead of geom_line() bridging straight across it.
    plot_dt[, ld_w_run := rleid(contig_run, is.na(ld_w_scaled)), by = Chr]
  }
  
  ## Add recombination rate per LD-decay window
  if (plot_rec_rate && !is.null(rec_dt)) {
    
    rec_dt <- copy(rec_dt)
    
    if (!"Chr" %in% names(rec_dt)) {
      rec_dt[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
    }
    
    setDT(rec_dt)
    
    plot_dt[, win_id := .I]
    
    rec_win <- rec_dt[
      plot_dt,
      on = .(Chr, pos >= start, pos <= end),
      allow.cartesian = TRUE
    ][
      ,
      .(
        rec_rate = mean(`cM/Mb`, na.rm = TRUE),
        rec_cM   = mean(cM, na.rm = TRUE)
      ),
      by = win_id
    ]
    
    plot_dt <- rec_win[plot_dt, on = "win_id"]
    
    plot_dt[, scale_fac_rec := max(a, na.rm = TRUE) /
              max(rec_rate, na.rm = TRUE),
            by = Chr]
    
    plot_dt[, rec_scaled := rec_rate * scale_fac_rec]

    ## same gap-preserving fix as ld_w_run above
    plot_dt[, rec_run := rleid(contig_run, is.na(rec_scaled)), by = Chr]
  }
  
  ## Base plot
  plot_dt[is.na(a),a:=0]
  p <- ggplot(plot_dt, aes(mid / 1e6)) +
    geom_line(aes(y = a, color = "a (windowed decay rate)", group = interaction(Chr, contig_run)), linewidth = 0.6) +
    facet_wrap(~ Chr, scales = "free", ncol = ncol) +
    labs(
      x = "Chromosome position (Mbp)",
      y = "Scaled tracks"
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  ## Add threshold
  if (!is.null(a_threshold)) {
    p <- p +
      geom_hline(
        yintercept = a_threshold,
        linetype = 2,
        linewidth = 0.5
      )
  }
  
  ## Add LD-cluster SNPs
  if (!is.null(map)) {
    
    
    p <- p +
      geom_point(
        data = map,
        aes(x = Pos / 1e6, y = ld_w_095 / 100),
        inherit.aes = FALSE,
        size = 0.25,
        alpha = 0.5,
        col = "grey80"
      )
  }
  
  ## Add window median ld_w
  if (plot_ld_w && "ld_w_scaled" %in% names(plot_dt)) {
    p <- p +
      geom_line(
        aes(y = ld_w_scaled, color = "ld_w (window median)", group = interaction(Chr, ld_w_run)),
        linewidth = 0.4,
        alpha = 0.9
      )
  }

  ## Add recombination rate
  if (plot_rec_rate && "rec_scaled" %in% names(plot_dt)) {
    p <- p +
      geom_line(
        aes(y = rec_scaled, color = "recombination rate (recmap)", group = interaction(Chr, rec_run)),
        linewidth = 0.5,
        alpha = 0.8
      )
  }

  p <- p +
    scale_color_manual(
      name = NULL,
      values = c(
        "a (windowed decay rate)"      = "salmon",
        "ld_w (window median)"         = "forestgreen",
        "recombination rate (recmap)"  = "steelblue"
      )
    ) +
    theme(legend.position = "bottom")
  
  ## Highlight outliers
  if (highlight_outliers && !is.null(outliers)) {
    
    outliers <- copy(outliers)
    
    if (!"mid" %in% names(outliers)) {
      outliers[, mid := rowMeans(.SD, na.rm = TRUE),
               .SDcols = c("start", "end")]
    }
    
    p <- p +
      geom_point(
        data = outliers,
        aes(x = mid / 1e6, y = a),
        inherit.aes = FALSE,
        size = 1,
        col = "red"
      )
  }
  
  p
}

p <- plot_ld_decay_tracks(
  ld_decay = ld_decay,
  ld_w_095 = ld_w_095,
  rec_dt = rec_map,
  map = map,
  plot_ld_w = TRUE,
  plot_rec_rate = TRUE,
  outliers = NULL,
  highlight_outliers = TRUE
)

ggsave("./Figures/ld_tracks.png",p, width = 20, height = 20)

# ------------------------------------------------------------
# comparing ld_w, a and rec_rate
# ------------------------------------------------------------
library(pROC)
library(patchwork)

## --- attach recombination-map estimates to each LD-decay window --- ##
## (same start/end window resolution as plot_dt / plot_ld_decay_tracks())
win_dt <- copy(plot_dt)
win_dt[, win_id := .I]

rec_win <- rec_map[
  win_dt,
  on = .(Chr, pos >= start, pos <= end),
  allow.cartesian = TRUE
][
  ,
  .(
    cM_Mb_med    = median(`cM/Mb`, na.rm = TRUE),   # native, smoothed recmap rate
    rec_rate_med = median(rec_rate, na.rm = TRUE),  # recomputed adjacent-marker rate
    n_map_pts    = .N
  ),
  by = win_id
]

comp_dt <- rec_win[win_dt, on = "win_id"]
comp_dt <- comp_dt[n_map_pts >= 5 & !is.na(ld_w_med) & !is.na(a)]

## Scope note: this comparison validates ld_w/a against the recmap, i.e. against
## the NEUTRAL recombination-rate landscape only -- the recmap comes from a lab
## F1 cross, where selection isn't acting, so it cannot validate (or refute)
## ld_w's sensitivity to selection-driven excess LD (sweeps, differentiation
## islands) in the actual wild hybrid population. A SNP-level test stratified by
## window-level `a` was tried and showed ld_w tracks little beyond what `a`
## already captures once conditioned on the recmap -- expected, since any
## selection signal wouldn't be present in this ground truth to begin with, not
## evidence against ld_w. Testing the selection-sensitivity claim needs a
## different validation (simulation with selection, or cross-population
## comparison), which is a separate LDscnR-methods question, not pursued here.

## --- correlations (Spearman: relationship is monotonic, not linear) --- ##
cor_spearman <- function(x, y) cor(x, y, method = "spearman", use = "pairwise.complete.obs")

cor_pairs <- list(
  c("a", "ld_w_med"),
  c("a", "cM_Mb_med"),
  c("a", "rec_rate_med"),
  c("ld_w_med", "cM_Mb_med"),
  c("ld_w_med", "rec_rate_med"),
  c("cM_Mb_med", "rec_rate_med")
)

cor_tbl <- rbindlist(lapply(cor_pairs, function(p) {
  data.table(
    x = p[1], y = p[2],
    rho_pooled     = cor_spearman(comp_dt[[p[1]]], comp_dt[[p[2]]]),
    rho_within_chr = comp_dt[TRUE, cor_spearman(get(p[1]), get(p[2])), by = Chr][TRUE, mean(V1, na.rm = TRUE)]
  )
}))
print(cor_tbl)

## rec_rate (diff(cM)/diff(pos) between adjacent markers) is ~65% exact zero
## genome-wide: with ~300 F1s, most adjacent-marker intervals see zero observed
## recombinants, so its signal is only informative in the upper part of its
## distribution. cM/Mb is smoothed and doesn't have this problem, so it is used
## below as the ground truth for ROC/AUC; rec_rate is kept in cor_tbl above for
## transparency only.
cat("Fraction of windows with rec_rate_med == 0:",
    comp_dt[, mean(rec_rate_med == 0, na.rm = TRUE)], "\n")

## --- scatter plots --- ##
p_ldw_recmap <- ggplot(comp_dt, aes(cM_Mb_med, ld_w_med)) +
  geom_point(alpha = 0.4, size = 0.8) +
  theme_bw(base_size = 14) +
  labs(x = "Recombination rate (cM/Mb, recmap)", y = expression(ld["w"]~"(window median)"))

p_a_recmap <- ggplot(comp_dt, aes(cM_Mb_med, a)) +
  geom_point(alpha = 0.4, size = 0.8) +
  theme_bw(base_size = 14) +
  labs(x = "Recombination rate (cM/Mb, recmap)", y = "Windowed decay rate (a)")

p_gt_agree <- ggplot(comp_dt, aes(cM_Mb_med, rec_rate_med)) +
  geom_point(alpha = 0.3, size = 0.8) +
  theme_bw(base_size = 14) +
  labs(x = "cM/Mb (recmap, smoothed)", y = "rec_rate (adjacent-marker diff(cM)/diff(pos))")

p_a_ldw <- ggplot(comp_dt, aes(a, ld_w_med)) +
  geom_point(alpha = 0.4, size = 0.8) +
  theme_bw(base_size = 14) +
  labs(x = "Windowed decay rate (a)", y = expression(ld["w"]~"(window median)"))

p_scatter <- (p_ldw_recmap + p_a_recmap + p_gt_agree + p_a_ldw) + plot_layout(nrow = 1)
ggsave("./Figures/p_ldw_a_recmap_scatter.png", p_scatter, width = 20, height = 5)

## --- ROC/AUC: can ld_w / a detect low-recombination windows? --- ##
## "low recombination" = bottom q of cM/Mb (centromeres, inversions). Repeated
## across a few thresholds so the conclusion doesn't hinge on one arbitrary cutoff.
q_grid <- c(0.05, 0.10, 0.25)

roc_list <- lapply(q_grid, function(q) {
  label <- as.integer(comp_dt$cM_Mb_med <= quantile(comp_dt$cM_Mb_med, q, na.rm = TRUE))
  list(
    q = q,
    roc_ld_w = roc(label, comp_dt$ld_w_med, direction = "auto", quiet = TRUE, ci = TRUE),
    roc_a    = roc(label, comp_dt$a,        direction = "auto", quiet = TRUE, ci = TRUE)
  )
})

roc_summary <- rbindlist(lapply(roc_list, function(x) {
  ci_ldw <- ci.auc(x$roc_ld_w)
  ci_a   <- ci.auc(x$roc_a)
  data.table(
    quantile = x$q,
    auc_ld_w = as.numeric(auc(x$roc_ld_w)), auc_ld_w_lo = ci_ldw[1], auc_ld_w_hi = ci_ldw[3],
    auc_a    = as.numeric(auc(x$roc_a)),    auc_a_lo    = ci_a[1],   auc_a_hi    = ci_a[3]
  )
}))
print(roc_summary)

roc_plots <- lapply(roc_list, function(x) {
  ggroc(list(ld_w = x$roc_ld_w, a = x$roc_a), linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 1, linetype = 2, col = "grey60") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("bottom ", x$q * 100, "% cM/Mb"),
      subtitle = sprintf("AUC ld_w=%.3f | a=%.3f", auc(x$roc_ld_w), auc(x$roc_a)),
      col = NULL
    )
})

p_roc <- wrap_plots(roc_plots, ncol = 3) + plot_layout(guides = "collect")
ggsave("./Figures/p_roc_low_recombination.png", p_roc, width = 13, height = 5)

saveRDS(
  list(comp_dt = comp_dt, cor_tbl = cor_tbl, roc_summary = roc_summary),
  "./data/ldw_a_recmap_comparison.rds"
)
