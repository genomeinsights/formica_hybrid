## =============================================================================
## Module E -- final figures for the write-up
## =============================================================================
## Regenerates the figures that support the reported results:
##   Fig 1  F_ST by DI bin, observed vs neutral null           (the dose-response)
##   Fig 2  sorted-locus fraction by DI bin (log scale)
##   Fig 3  LD-decay on the low-DI anchor markers, observed vs null
## and prints the numbers quoted in the text.
## =============================================================================

suppressMessages({library(data.table); library(ggplot2); library(patchwork)})

FORMICA  <- "/Users/petrikem/gitlab/formica_hybrid"
RDS      <- file.path(FORMICA, "data/moduleE_sim/moduleE_di_stratified.rds")
DISTRAT  <- file.path(FORMICA, "data/moduleE_sim/distrat")
POOL     <- file.path(FORMICA, "data/moduleE_founders/moduleE_founder_haplotypes.rds")
FVCF_DIR <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/SLiM/founders_maf015_DIstrat4000"
EMP_RD   <- file.path(FORMICA, "data/hybrids_only_maf005.Rdata")
FIGDIR   <- file.path(FORMICA, "Figures")
dir.create(FIGDIR, showWarnings = FALSE)

TAG <- "Naq30_Npol13"; BEST_K <- 1500; BEST_GEN <- 160
DI_BREAKS <- c(-Inf,-90,-75,-60,-50,-40,-30,-25,-20,-15,Inf)
ANCHOR_MAXDI <- -60
MAXDIST <- 250000L; MAXPART <- 100L
DIST_BREAKS <- c(0,1,2,5,10,20,50,100,150,250)*1000L; MIN_MAF <- 0.05
MATCH_SEED <- 1L

x <- readRDS(RDS)
emp_bin <- as.data.table(x$emp_bin)
best <- x$cells[[which.min(x$summ$ld_rmse)]]
sim_bin <- as.data.table(best$bins)

cmp <- merge(emp_bin, sim_bin, by="dibin", suffixes=c("_emp","_sim"))
cmp[, dibin := factor(dibin, levels=levels(emp_bin$dibin))]

## ---------------- Fig 1: F_ST by DI bin ----------------
d1 <- melt(cmp[, .(dibin, Observed=fst_emp, `Neutral null`=fst_sim)], id.vars="dibin",
           variable.name="source", value.name="fst")
p1 <- ggplot(d1, aes(dibin, fst, colour=source, group=source)) +
  geom_line(linewidth=0.9) + geom_point(size=2) +
  scale_colour_manual(values=c("Observed"="#d95f02","Neutral null"="#1b9e77"), name=NULL) +
  labs(x=NULL, y=expression(F[ST]~"among populations"),
       title="Differentiation scales with ancestry-informativeness in the data, not under drift") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="top",
        panel.grid.minor=element_blank())

## ---------------- Fig 2: sorted fraction by DI bin ----------------
d2 <- melt(cmp[, .(dibin, Observed=sorted_emp, `Neutral null`=sorted_sim)], id.vars="dibin",
           variable.name="source", value.name="sorted")
d2[sorted <= 0, sorted := NA]
p2 <- ggplot(d2, aes(dibin, sorted, colour=source, group=source)) +
  geom_line(linewidth=0.9) + geom_point(size=2) + scale_y_log10() +
  scale_colour_manual(values=c("Observed"="#d95f02","Neutral null"="#1b9e77"), name=NULL) +
  labs(x="DiagnosticIndex bin (low = neutral background, high = ancestry-informative)",
       y="fraction of loci sorted (log)") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none",
        panel.grid.minor=element_blank())

ggsave(file.path(FIGDIR,"moduleE_fig1_dose_response.pdf"), p1/p2, width=8, height=7)

## ---------------- Fig 3: LD-decay on the low-DI anchor ----------------
sim_markers <- unique(unlist(lapply(
  list.files(FVCF_DIR,"founders_ch.*[.]vcf$",full.names=TRUE),
  function(f) fread(f, skip="#CHROM", select=3, header=TRUE)[[1]])))
ph <- readRDS(POOL); pmap <- as.data.table(ph$map)
pj <- match(sim_markers, pmap$marker); ok <- !is.na(pj)
markers <- sim_markers[ok]; pj <- pj[ok]
DI <- pmap$DI[pj]; mk_chr <- pmap$Chr[pj]; mk_pos <- pmap$Pos[pj]
anchor <- which(DI <= ANCHOR_MAXDI); chr_idx <- split(anchor, mk_chr[anchor])

ld_anchor <- function(dose) {
  out <- list()
  for (ch in names(chr_idx)) {
    i <- chr_idx[[ch]]; o <- order(mk_pos[i]); i <- i[o]
    pos <- mk_pos[i]; G <- dose[i,,drop=FALSE]
    p <- rowMeans(G,na.rm=TRUE)/2
    jj <- which(is.finite(p) & pmin(p,1-p) >= MIN_MAF); if (length(jj)<2L) next
    pos <- pos[jj]; G <- G[jj,,drop=FALSE]
    for (a in seq_len(length(pos)-1L)) {
      b <- (a+1L):min(a+MAXPART,length(pos)); d <- pos[b]-pos[a]
      k <- d<=MAXDIST; if(!any(k)) next; b <- b[k]; d <- d[k]
      r <- suppressWarnings(cor(G[a,], t(G[b,,drop=FALSE]), use="pairwise.complete.obs"))
      out[[length(out)+1L]] <- data.table(dist=d, r2=as.numeric(r)^2)
    }
  }
  res <- rbindlist(out)[is.finite(r2)]
  res[, db := cut(dist, DIST_BREAKS, labels=FALSE)]
  res[!is.na(db), .(r2=mean(r2), dmid=(DIST_BREAKS[db]+DIST_BREAKS[db+1L])/2), by=db][order(db)]
}
e <- new.env(); load(EMP_RD, envir=e)
GTh <- e$GTs_hybrids_005; sdh <- as.data.table(e$sample_data); maph <- as.data.table(e$map_hyb_005)
ci <- match(markers, maph$marker)
pops <- setdiff(unique(sdh$Population), c("aquilonia_parent","polyctena_parent"))
emp_ns <- sapply(pops, function(p) sum(sdh$Population==p))
emp_ld <- rbindlist(lapply(pops, function(p)
  ld_anchor(t(GTh[sdh$Population==p, ci, drop=FALSE]))))[, .(r2=mean(r2), dmid=dmid[1]), by=db][order(db)]

read_dosage <- function(vcf) {
  dt <- fread(vcf, skip="#CHROM", header=TRUE, sep="\t", showProgress=FALSE)
  mk <- paste(sub("ch","Chr",dt[[1]]), dt[[2]], sep=":")
  G <- as.matrix(dt[,10:ncol(dt)]); gt <- sub(":.*$","",G)
  d <- matrix(suppressWarnings(as.integer(sub("[|/].*$","",gt))) +
              suppressWarnings(as.integer(sub("^.*[|/]","",gt))), nrow=nrow(G))
  m <- match(markers, mk); keep <- !is.na(m)
  D <- matrix(0L, length(markers), ncol(d)); D[keep,] <- d[m[keep],,drop=FALSE]; D
}
sim_curves <- rbindlist(lapply(1:20, function(i){
  f <- file.path(DISTRAT, sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf", TAG, BEST_K, i, BEST_GEN))
  if (!file.exists(f)) return(NULL)
  D <- read_dosage(f); tgt <- emp_ns[(i-1)%%length(emp_ns)+1]
  if (ncol(D) > tgt) { set.seed(MATCH_SEED+i); D <- D[, sample(ncol(D), tgt), drop=FALSE] }
  ld_anchor(D)
}))[, .(r2=mean(r2), dmid=dmid[1]), by=db][order(db)]

d3 <- rbind(emp_ld[, .(dmid, r2, source="Observed")],
            sim_curves[, .(dmid, r2, source="Neutral null")])
p3 <- ggplot(d3, aes(dmid/1000, r2, colour=source)) +
  geom_line(linewidth=0.9) + geom_point(size=1.8) + scale_x_log10() +
  scale_colour_manual(values=c("Observed"="#d95f02","Neutral null"="#1b9e77"), name=NULL) +
  labs(x="distance (kb, log scale)", y=expression("composite "*r^2),
       title=expression("LD decay on low-DI anchor markers (DI"<=-60*")")) +
  theme_bw(base_size=11) + theme(legend.position="top", panel.grid.minor=element_blank())
ggsave(file.path(FIGDIR,"moduleE_fig2_ld_anchor.pdf"), p3, width=7, height=4.5)

## ---------------- numbers for the text ----------------
cat("\n=== numbers quoted in the write-up ===\n")
print(cmp[, .(dibin, fst_emp=round(fst_emp,3), fst_sim=round(fst_sim,3),
              pi_emp=round(pi_emp,3), pi_sim=round(pi_sim,3),
              srt_emp=round(sorted_emp,3), srt_sim=round(sorted_sim,3),
              ratio=round(ifelse(sorted_sim>0, sorted_emp/sorted_sim, NA),1))])
cat(sprintf("\nanchor-matched cell: K=%d gen=%d (ld_rmse=%.4f)\n", best$K, best$gen, best$ld_rmse))
cat(sprintf("neutral background (DI<=-90): observed F_ST=%.3f vs null %.3f\n",
            cmp[1,fst_emp], cmp[1,fst_sim]))
cat(sprintf("implied Nm from observed background F_ST: %.1f\n",
            (1/cmp[1,fst_emp] - 1)/4))
cat("\nfigures written:\n  Figures/moduleE_fig1_dose_response.pdf\n  Figures/moduleE_fig2_ld_anchor.pdf\n")
