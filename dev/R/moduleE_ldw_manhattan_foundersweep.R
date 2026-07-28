## =============================================================================
## Module E -- local-LD ("ld_w") Manhattan: empirical vs the FOUNDER-NUMBER sweep
## =============================================================================
## Per-SNP local LD (median composite r^2 among neighbours within WIN_BP, MAF>=0.1)
## along the genome, POOLED across all 20 demes (between-population differentiation
## contributes to LD -- the intended diagnostic). Compares the empirical hybrids to
## the neutral null at four founding depths (30/13 deep -> 1200/520 whole-pool) at a
## fixed generation. Tests whether a shallower founding bottleneck lowers the pooled
## LD level toward empirical (the Wahlund prediction: pooled LD tracks neutral F_ST),
## and whether the landscape SHAPE (recombination-driven peaks) is reproduced
## (per-chromosome Spearman vs empirical).
##
## Sim demes are subsampled to the empirical per-population sizes (rep i -> emp_ns[i])
## before pooling, so both sides pool the SAME size distribution. Same 15k DI-strat
## markers, same recombination map on both sides. Negatives clamped to 0.
##
## Output: data/moduleE_sim/moduleE_ldw_manhattan_foundersweep.rds
##         Figures/moduleE_ldw_manhattan_foundersweep.{pdf,png}
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2) })

FORMICA  <- "/Users/petrikem/gitlab/formica_hybrid"
SIMDIR   <- file.path(FORMICA, "moduleE_slim/output_foundersweep")
FVCF_DIR <- file.path(FORMICA, "moduleE_slim/founders/maf015_DIstrat1500_x40")
BUNDLE   <- file.path(FORMICA, "moduleE_slim/inputs/empirical_bundle.rds")
EMP_RD   <- file.path(FORMICA, "data/hybrids_only_maf005.Rdata")
OUTRDS   <- file.path(FORMICA, "data/moduleE_sim/moduleE_ldw_manhattan_foundersweep.rds")
FIGDIR   <- file.path(FORMICA, "Figures"); dir.create(FIGDIR, showWarnings = FALSE)

K <- 6250; GEN <- 60L
FOUNDINGS <- c("Naq30_Npol13","Naq150_Npol65","Naq450_Npol195","Naq1200_Npol520")
FLAB <- c("Naq30_Npol13"="found 30/13 (deep)", "Naq150_Npol65"="found 150/65",
          "Naq450_Npol195"="found 450/195", "Naq1200_Npol520"="found 1200/520 (whole pool)")
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
  pmax(res, 0)
}
read_dose <- function(f) {
  dt <- fread(f, skip="#CHROM", header=TRUE, sep="\t", showProgress=FALSE)
  mk <- paste(sub("ch","Chr",dt[[1]]), dt[[2]], sep=":")
  G <- as.matrix(dt[,10:ncol(dt)]); gt <- sub(":.*$","",G)
  d <- matrix(suppressWarnings(as.integer(sub("[|/].*$","",gt))) +
              suppressWarnings(as.integer(sub("^.*[|/]","",gt))), nrow=nrow(G))
  m <- match(markers, mk); keep <- !is.na(m)
  D <- matrix(0L, length(markers), ncol(d)); D[keep,] <- d[m[keep],,drop=FALSE]; t(D)
}

## ---- canonical marker set (15k) --------------------------------------------
fm <- rbindlist(lapply(list.files(FVCF_DIR, "founders_ch.*[.]vcf$", full.names=TRUE),
  function(f) fread(f, skip="#CHROM", select=c(1,2,3), header=TRUE, col.names=c("chs","Pos","marker"))))
fm[, `:=`(Chr=sub("ch","Chr",chs), chr_id=as.integer(sub("ch","",chs)))]
setorder(fm, chr_id, Pos); markers <- fm$marker
message(sprintf("markers: %d over %d chromosomes", length(markers), uniqueN(fm$chr_id)))
emp_ns <- readRDS(BUNDLE)$emp_ns                              # per-pop diploid sizes (3-20)

## ---- empirical (all hybrids pooled) ----------------------------------------
message("[1] empirical ...")
e <- new.env(); load(EMP_RD, envir=e)
GTh <- e$GTs_hybrids_005; sdh <- as.data.table(e$sample_data); maph <- as.data.table(e$map_hyb_005)
hyb <- which(!sdh$Population %in% c("aquilonia_parent","polyctena_parent"))
emp_gt <- GTh[hyb, match(markers, maph$marker), drop=FALSE]; emp_gt[is.na(emp_gt)] <- 0L
N <- nrow(emp_gt); rm(e); gc()
message(sprintf("    empirical N=%d", N))
emp_w <- local_ld(emp_gt, fm$Chr, fm$Pos)

## ---- each founding: pool 20 demes, rep i subsampled to emp_ns[i] ------------
sim_list <- lapply(FOUNDINGS, function(TAG) {
  message("[2] ", TAG, " gen", GEN, " ...")
  parts <- lapply(1:20, function(i) {
    f <- file.path(SIMDIR, sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf", TAG, K, i, GEN))
    if (!file.exists(f)) return(NULL)
    g <- read_dose(f); tgt <- emp_ns[(i-1)%%length(emp_ns)+1]
    set.seed(100+i); if (nrow(g) > tgt) g <- g[sample(nrow(g), tgt), , drop=FALSE]; g
  })
  gt <- do.call(rbind, parts)
  data.table(marker=markers, ld_w=local_ld(gt, fm$Chr, fm$Pos), src=unname(FLAB[TAG]), n_pooled=nrow(gt))
})
dt <- rbind(data.table(marker=markers, ld_w=emp_w, src="empirical", n_pooled=N), rbindlist(sim_list))
dt <- fm[, .(marker, Chr, Pos, chr_id)][dt, on="marker"]
saveRDS(list(dt=dt, gen=GEN, K=K, emp_N=N), OUTRDS)

## ---- loess trend per (chromosome, source) ----------------------------------
src_lv <- c("empirical", unname(FLAB[FOUNDINGS]))
dt[, src := factor(src, levels=src_lv)]
trend <- dt[is.finite(ld_w), {
  if (.N >= 10) { fit <- pmax(0, predict(loess(ld_w ~ Pos, span=0.2))); .(Pos=Pos, y=fit) }
  else .(Pos=Pos, y=ld_w)
}, by=.(chr_id, Chr, src)]
trend[, Chr := factor(Chr, levels=paste0("Chr", sort(unique(fm$chr_id))))]

cols <- c("empirical"="#d95f02",
          setNames(colorRampPalette(c("#0b3d2e","#1b9e77","#7fcdb8"))(length(FOUNDINGS)),
                   unname(FLAB[FOUNDINGS])))
p <- ggplot(trend, aes(Pos/1e6, y, colour=src)) +
  geom_line(linewidth=0.5) +
  facet_grid(. ~ Chr, scales="free_x", space="free_x") +
  scale_colour_manual(values=cols, name=NULL) +
  labs(x="position (Mbp)", y="local LD (median r^2, 100kb window)",
       title=sprintf("Local-LD (ld_w) landscape: empirical vs founder-number sweep (K%d, gen %d, pooled)", K, GEN),
       subtitle=sprintf("same %d DI-stratified markers, per-deme N matched to empirical (pooled N~%d); loess trend", length(markers), N)) +
  theme_bw(base_size=8) +
  theme(panel.grid.minor=element_blank(), panel.spacing=unit(1,"pt"),
        strip.text.x=element_text(size=5, angle=90), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position="top")
ggsave(file.path(FIGDIR,"moduleE_ldw_manhattan_foundersweep.pdf"), p, width=14, height=4.2)
ggsave(file.path(FIGDIR,"moduleE_ldw_manhattan_foundersweep.png"), p, width=14, height=4.2, dpi=130)

## ---- summary: mean level + per-chromosome landscape correlation vs empirical -
cat("\n=== mean local-LD (ld_w) and per-chromosome landscape Spearman vs empirical ===\n")
pc <- dcast(dt[, .(m=mean(ld_w,na.rm=TRUE)), by=.(chr_id, src)], chr_id ~ src, value.var="m")
for (s in src_lv) {
  ov  <- dt[src==s, mean(ld_w, na.rm=TRUE)]
  rho <- if (s=="empirical") NA else cor(pc[["empirical"]], pc[[s]], method="spearman", use="complete.obs")
  cat(sprintf("  %-28s mean ld_w=%.3f  per-chr Spearman vs emp=%s\n",
              s, ov, ifelse(is.na(rho),"-",sprintf("%.3f",rho))))
}
cat("\n(empirical mean is the target level; Spearman ~ how well the null reproduces the LD landscape shape)\n")
cat("saved: ", OUTRDS, " ; Figures/moduleE_ldw_manhattan_foundersweep.{pdf,png}\n")
