## =============================================================================
## Module E -- local-LD ("ld_w") Manhattan: empirical vs neutral null over generations
## =============================================================================
## Per-SNP local LD (median composite r^2 among neighbours within WIN_BP, MAF>=0.1)
## along the genome, for the empirical hybrids and for the pooled neutral demes at
## several generations (fixed K). POOLED across all populations/demes on purpose:
## between-population differentiation contributes to LD and is a useful diagnostic.
## Same 40k DI-stratified markers, matched N, same recombination map on both sides.
##
## Negative values are clamped to 0 (r^2 medians are already >=0; this also clamps
## the loess trend, which can dip below 0 as a smoothing artefact).
##
## Output: data/moduleE_sim/moduleE_ldw_manhattan.rds
##         Figures/moduleE_ldw_manhattan.{pdf,png}
## =============================================================================

suppressMessages({ library(data.table); library(ggplot2) })

FORMICA  <- "/Users/petrikem/gitlab/formica_hybrid"
SIMDIR   <- file.path(FORMICA, "data/moduleE_sim/distrat")
FVCF_DIR <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/SLiM/founders_maf015_DIstrat4000"
EMP_RD   <- file.path(FORMICA, "data/hybrids_only_maf005.Rdata")
OUTRDS   <- file.path(FORMICA, "data/moduleE_sim/moduleE_ldw_manhattan.rds")
FIGDIR   <- file.path(FORMICA, "Figures"); dir.create(FIGDIR, showWarnings = FALSE)

TAG <- "Naq30_Npol13"; K <- 6250; GENS <- c(20, 60, 100, 160)   # earlier -> later
WIN_BP <- 100000L; MAXPART <- 100L; MAF_MIN <- 0.1

## per-SNP local LD (median r^2 with downstream partners within WIN_BP); >=0
local_ld <- function(geno, Chr, Pos) {
  res <- rep(NA_real_, ncol(geno))
  for (ch in unique(Chr)) {
    j <- which(Chr == ch); o <- order(Pos[j]); j <- j[o]
    pos <- Pos[j]; G <- geno[, j, drop=FALSE]
    p <- colMeans(G, na.rm=TRUE)/2; poly <- is.finite(p) & pmin(p,1-p) >= MAF_MIN
    for (a in seq_along(j)) {
      if (!poly[a]) next
      b <- which(pos <= pos[a] + WIN_BP & seq_along(pos) > a & poly)
      if (length(b) > MAXPART) b <- b[seq_len(MAXPART)]
      if (!length(b)) next
      r <- suppressWarnings(cor(G[,a], G[, b, drop=FALSE], use="pairwise.complete.obs"))
      res[j[a]] <- median(as.numeric(r)^2, na.rm=TRUE)
    }
  }
  pmax(res, 0)                                          # clamp negatives -> 0
}

## ---- canonical marker set ---------------------------------------------------
fm <- rbindlist(lapply(list.files(FVCF_DIR, "founders_ch.*[.]vcf$", full.names=TRUE),
  function(f) fread(f, skip="#CHROM", select=c(1,2,3), header=TRUE, col.names=c("chs","Pos","marker"))))
fm[, `:=`(Chr=sub("ch","Chr",chs), chr_id=as.integer(sub("ch","",chs)))]
setorder(fm, chr_id, Pos); markers <- fm$marker
message(sprintf("markers: %d over %d chromosomes", length(markers), uniqueN(fm$chr_id)))

read_dose <- function(f) {
  dt <- fread(f, skip="#CHROM", header=TRUE, sep="\t", showProgress=FALSE)
  mk <- paste(sub("ch","Chr",dt[[1]]), dt[[2]], sep=":")
  G <- as.matrix(dt[,10:ncol(dt)]); gt <- sub(":.*$","",G)
  d <- matrix(suppressWarnings(as.integer(sub("[|/].*$","",gt))) +
              suppressWarnings(as.integer(sub("^.*[|/]","",gt))), nrow=nrow(G))
  m <- match(markers, mk); keep <- !is.na(m)
  D <- matrix(0L, length(markers), ncol(d)); D[keep,] <- d[m[keep],,drop=FALSE]; t(D)
}

## ---- empirical (all hybrids pooled) ----------------------------------------
message("[1] empirical ...")
e <- new.env(); load(EMP_RD, envir=e)
GTh <- e$GTs_hybrids_005; sdh <- as.data.table(e$sample_data); maph <- as.data.table(e$map_hyb_005)
hyb <- which(!sdh$Population %in% c("aquilonia_parent","polyctena_parent"))
emp_gt <- GTh[hyb, match(markers, maph$marker), drop=FALSE]; emp_gt[is.na(emp_gt)] <- 0L
N <- nrow(emp_gt)
message(sprintf("    empirical N=%d", N))
emp_w <- local_ld(emp_gt, fm$Chr, fm$Pos)

## ---- simulated cells: pool 20 demes, match N=164 ---------------------------
sim_list <- lapply(GENS, function(g) {
  message("[2] sim K", K, " gen", g, " ...")
  fs <- file.path(SIMDIR, sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf", TAG, K, 1:20, g))
  fs <- fs[file.exists(fs)]
  gt <- do.call(rbind, lapply(fs, read_dose))
  set.seed(1); if (nrow(gt) > N) gt <- gt[sample(nrow(gt), N), ]
  data.table(marker=markers, ld_w=local_ld(gt, fm$Chr, fm$Pos), src=sprintf("sim gen%d", g))
})
dt <- rbind(data.table(marker=markers, ld_w=emp_w, src="empirical"), rbindlist(sim_list))
dt <- fm[, .(marker, Chr, Pos, chr_id)][dt, on="marker"]
saveRDS(list(dt=dt, cell=sprintf("%s K%d", TAG, K), gens=GENS, N=N), OUTRDS)

## ---- clamped loess trend per (chromosome, source) --------------------------
src_lv <- c("empirical", sprintf("sim gen%d", GENS))
dt[, src := factor(src, levels=src_lv)]
trend <- dt[is.finite(ld_w), {
  if (.N >= 10) { fit <- pmax(0, predict(loess(ld_w ~ Pos, span=0.2))); .(Pos=Pos, y=fit) }
  else .(Pos=Pos, y=ld_w)
}, by=.(chr_id, Chr, src)]
trend[, Chr := factor(Chr, levels=paste0("Chr", sort(unique(fm$chr_id))))]

cols <- c(empirical="#d95f02",
          setNames(colorRampPalette(c("#bce4d8","#1b9e77","#0b3d2e"))(length(GENS)),
                   sprintf("sim gen%d", GENS)))
p <- ggplot(trend, aes(Pos/1e6, y, colour=src)) +
  geom_line(linewidth=0.5) +
  facet_grid(. ~ Chr, scales="free_x", space="free_x") +
  scale_colour_manual(values=cols, name=NULL) +
  labs(x="position (Mbp)", y="local LD (median r^2, 100kb window)",
       title=sprintf("Local-LD landscape: empirical vs neutral null over generations (%s K%d, pooled)", TAG, K),
       subtitle=sprintf("same %d DI-stratified markers, matched N=%d; loess trend (clamped >=0)", length(markers), N)) +
  theme_bw(base_size=8) +
  theme(panel.grid.minor=element_blank(), panel.spacing=unit(1,"pt"),
        strip.text.x=element_text(size=5, angle=90), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position="top")
ggsave(file.path(FIGDIR,"moduleE_ldw_manhattan.pdf"), p, width=14, height=4.2)
ggsave(file.path(FIGDIR,"moduleE_ldw_manhattan.png"), p, width=14, height=4.2, dpi=130)

## ---- summary: mean level + landscape correlation vs empirical --------------
cat("\n=== mean local-LD and per-chromosome landscape correlation vs empirical ===\n")
pc <- dcast(dt[, .(m=mean(ld_w,na.rm=TRUE)), by=.(chr_id, src)], chr_id ~ src, value.var="m")
for (s in src_lv) {
  ov <- dt[src==s, mean(ld_w, na.rm=TRUE)]
  rho <- if (s=="empirical") NA else cor(pc[["empirical"]], pc[[s]], method="spearman", use="complete.obs")
  cat(sprintf("  %-12s mean=%.3f  per-chr Spearman vs emp=%s\n", s, ov, ifelse(is.na(rho),"-",sprintf("%.3f",rho))))
}
cat("\nsaved: ", OUTRDS, " ; Figures/moduleE_ldw_manhattan.{pdf,png}\n")
