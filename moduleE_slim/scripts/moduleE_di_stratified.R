## =============================================================================
## Module E -- DI-STRATIFIED comparison: internal anchor + dose-response
## =============================================================================
## Universe = parent_maf >= 0.15 (the locked Module A gate), 40,000 markers drawn
## equally from 10 DI bins spanning the full DI range (-168 .. -4).
##
## THE DESIGN
##   - LOW-DI bins cannot be driven by ancestry sorting, so they are an INTERNAL
##     demographic anchor -- same runs, same markers, no separate calibration.
##   - The gap vs DI is the test:
##        gap RISES with DI  -> sorting excess is real (dose-response)
##        gap FLAT across DI -> demographic misspecification, not sorting
## This replaces the single fragile LD-decay anchor everything previously hung on.
##
## Both datasets are stratified by the SAME DI and the SAME pooled-parental MAF
## (taken from the phased founder pool, which is what the sim was seeded from),
## so the bins are identical on both sides. Sim demes are subsampled to the
## empirical per-population sizes (sort_class and r^2 are both sample-size
## sensitive).
##
## Output: data/moduleE_sim/moduleE_di_stratified.rds
##         Figures/moduleE_di_stratified.pdf
## =============================================================================

suppressMessages({library(data.table); library(ggplot2); library(parallel)})
source("dev/R/parallelism_stats.R")

FORMICA  <- "/Users/petrikem/gitlab/formica_hybrid"
DISTRAT  <- file.path(FORMICA, "data/moduleE_sim/distrat")
POOL     <- file.path(FORMICA, "data/moduleE_founders/moduleE_founder_haplotypes.rds")
FVCF_DIR <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/SLiM/founders_maf015_DIstrat4000"
EMP_RD   <- file.path(FORMICA, "data/hybrids_only_maf005.Rdata")
EMPTGT   <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/moduleE/inputs/empirical_targets.rds"
OUTRDS   <- file.path(FORMICA, "data/moduleE_sim/moduleE_di_stratified.rds")
FIGDIR   <- file.path(FORMICA, "Figures")

TAG  <- "Naq30_Npol13"
KS   <- c(1500, 6250)
GENS <- seq(20, 160, by = 20)
DI_BREAKS <- c(-Inf,-90,-75,-60,-50,-40,-30,-25,-20,-15,Inf)
ANCHOR_MAXDI <- -60          # bins at or below this DI pool into the LD anchor
FIX_TH <- 0.15; SORT_TH <- 0.5; MIN_PARENT_MAF <- 0.15; NULL_PROB <- 0.5
MAXDIST <- 250000L; MAXPART <- 100L
DIST_BREAKS <- c(0,1,2,5,10,20,50,100,150,250)*1000L; MIN_MAF <- 0.05
CORES <- 6L; MATCH_SEED <- 1L
dir.create(FIGDIR, showWarnings = FALSE)

## ---- markers, DI bins, parental reference (one source for both datasets) ----
message("[1] markers / DI bins ...")
sim_markers <- unique(unlist(lapply(
  list.files(FVCF_DIR, "founders_ch.*[.]vcf$", full.names=TRUE),
  function(f) fread(f, skip="#CHROM", select=3, header=TRUE)[[1]])))
ph <- readRDS(POOL); pmap <- as.data.table(ph$map)
pj <- match(sim_markers, pmap$marker); ok <- !is.na(pj)
markers <- sim_markers[ok]; pj <- pj[ok]
DI <- pmap$DI[pj]; mk_chr <- pmap$Chr[pj]; mk_pos <- pmap$Pos[pj]
dibin <- cut(DI, DI_BREAKS)
f_aq <- colMeans(ph$aqu[,pj,drop=FALSE]); f_pol <- colMeans(ph$pol[,pj,drop=FALSE])
p0 <- (colSums(ph$aqu[,pj,drop=FALSE])+colSums(ph$pol[,pj,drop=FALSE]))/(nrow(ph$aqu)+nrow(ph$pol))
parent_maf <- pmin(p0,1-p0); names(parent_maf) <- markers
anchor <- which(DI <= ANCHOR_MAXDI)
chr_idx_anchor <- split(anchor, mk_chr[anchor])
message(sprintf("    %d markers, %d DI bins; LD anchor = %d markers (DI <= %d)",
                length(markers), nlevels(dibin), length(anchor), ANCHOR_MAXDI))

## ---- empirical ----
message("[2] empirical ...")
e <- new.env(); load(EMP_RD, envir=e)
GTh <- e$GTs_hybrids_005; sdh <- as.data.table(e$sample_data); maph <- as.data.table(e$map_hyb_005)
ci <- match(markers, maph$marker)
pops <- setdiff(unique(sdh$Population), c("aquilonia_parent","polyctena_parent"))
emp_ns <- sapply(pops, function(p) sum(sdh$Population==p))
emp_M <- do.call(rbind, lapply(pops, function(p) colMeans(GTh[sdh$Population==p, ci, drop=FALSE], na.rm=TRUE)))
rownames(emp_M) <- pops; colnames(emp_M) <- markers

## ---- helpers ----
read_dosage <- function(vcf) {
  dt <- fread(vcf, skip="#CHROM", header=TRUE, sep="\t", showProgress=FALSE)
  mk <- paste(sub("ch","Chr",dt[[1]]), dt[[2]], sep=":")
  G <- as.matrix(dt[,10:ncol(dt)]); gt <- sub(":.*$","",G)
  d <- matrix(suppressWarnings(as.integer(sub("[|/].*$","",gt))) +
              suppressWarnings(as.integer(sub("^.*[|/]","",gt))), nrow=nrow(G))
  m <- match(markers, mk); keep <- !is.na(m)
  D <- matrix(0L, length(markers), ncol(d)); D[keep,] <- d[m[keep],,drop=FALSE]; D
}
ld_anchor <- function(dose) {          # LD-decay on the low-DI anchor markers only
  out <- list()
  for (ch in names(chr_idx_anchor)) {
    i <- chr_idx_anchor[[ch]]; o <- order(mk_pos[i]); i <- i[o]
    pos <- mk_pos[i]; G <- dose[i,,drop=FALSE]
    p <- rowMeans(G, na.rm=TRUE)/2
    jj <- which(is.finite(p) & pmin(p,1-p) >= MIN_MAF); if (length(jj)<2L) next
    pos <- pos[jj]; G <- G[jj,,drop=FALSE]
    for (a in seq_len(length(pos)-1L)) {
      b <- (a+1L):min(a+MAXPART,length(pos)); d <- pos[b]-pos[a]
      k <- d<=MAXDIST; if (!any(k)) next; b <- b[k]; d <- d[k]
      r <- suppressWarnings(cor(G[a,], t(G[b,,drop=FALSE]), use="pairwise.complete.obs"))
      out[[length(out)+1L]] <- data.table(dist=d, r2=as.numeric(r)^2)
    }
  }
  if (!length(out)) return(NULL)
  res <- rbindlist(out)[is.finite(r2)]
  res[, db := cut(dist, DIST_BREAKS, labels=FALSE)]
  res[!is.na(db), .(r2=mean(r2)), by=db][order(db)]
}
run_sort <- function(mm) {
  pm <- rbind(mm, aquilonia_parent=f_aq*2, polyctena_parent=f_pol*2); colnames(pm) <- markers
  parallelism_stats(list(pop_means=pm), hybrid_pops=rownames(mm),
    aqu_pops="aquilonia_parent", pol_pops="polyctena_parent", fix_th=FIX_TH,
    sort_th=SORT_TH, null_prob=NULL_PROB, parent_maf=parent_maf, min_parent_maf=MIN_PARENT_MAF)
}
by_bin <- function(mm, s) {
  p <- mm/2; pb <- colMeans(p,na.rm=TRUE)
  Hs <- colMeans(2*p*(1-p),na.rm=TRUE); Ht <- 2*pb*(1-pb)
  sc <- s[match(markers, marker), sort_class]
  data.table(dibin=dibin, Hs=Hs, Ht=Ht, sorted=!is.na(sc) & sc!="unsorted",
             tested=!is.na(sc))[, .(pi=mean(Hs,na.rm=TRUE),
      fst=sum(Ht-Hs,na.rm=TRUE)/sum(Ht,na.rm=TRUE),
      sorted=sum(sorted)/sum(tested)), by=dibin][order(dibin)]
}

emp_sort <- run_sort(emp_M); emp_bin <- by_bin(emp_M, emp_sort)
## LD anchor must be computed WITHIN each population and averaged -- pooling all
## individuals across differentiated populations manufactures LD (Wahlund) and
## makes the anchor meaningless (it inflates empirical r^2 so no cell can match).
emp_ld_anchor <- rbindlist(lapply(pops, function(p)
    ld_anchor(t(GTh[sdh$Population == p, ci, drop = FALSE]))
  ))[, .(r2 = mean(r2)), by = db][order(db)]

## ---- simulated cells ----
message("[3] simulated cells ...")
do_cell <- function(K, gen) {
  pm <- matrix(NA_real_, 20, length(markers)); curves <- list()
  for (i in 1:20) {
    f <- file.path(DISTRAT, sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf", TAG, K, i, gen))
    if (!file.exists(f)) next
    D <- read_dosage(f)
    tgt <- emp_ns[(i-1) %% length(emp_ns) + 1]; ni <- ncol(D)
    if (ni > tgt) { set.seed(MATCH_SEED+i); D <- D[, sample(ni, tgt), drop=FALSE] }
    pm[i,] <- rowMeans(D, na.rm=TRUE)
    cv <- ld_anchor(D); if (!is.null(cv)) curves[[length(curves)+1L]] <- cv
  }
  keep <- !is.na(pm[,1]); pm <- pm[keep,,drop=FALSE]
  if (nrow(pm) < 10) return(NULL)
  rownames(pm) <- paste0("deme",seq_len(nrow(pm))); colnames(pm) <- markers
  ld <- rbindlist(curves)[, .(r2=mean(r2)), by=db][order(db)]
  ev <- setNames(emp_ld_anchor$r2, emp_ld_anchor$db); m <- match(ld$db, as.integer(names(ev)))
  b <- by_bin(pm, run_sort(pm))
  list(K=K, gen=gen, ld_rmse=sqrt(mean((ld$r2-ev[m])^2, na.rm=TRUE)), bins=b)
}
cells <- CJ(K=KS, gen=GENS, sorted=FALSE)
out <- mcMap(function(K,g) do_cell(K,g), cells$K, cells$gen, mc.cores=CORES)
out <- out[!vapply(out, is.null, logical(1))]
summ <- rbindlist(lapply(out, function(z) data.table(K=z$K, gen=z$gen, ld_rmse=z$ld_rmse)))
saveRDS(list(cells=out, summ=summ, emp_bin=emp_bin, emp_ld_anchor=emp_ld_anchor), OUTRDS)

cat("\n=== LD-decay anchor (LOW-DI markers only): sim vs empirical ===\n")
print(summ[order(ld_rmse)][1:8, .(K, gen, ld_rmse=round(ld_rmse,4))])
best <- out[[which.min(summ$ld_rmse)]]
cat(sprintf("\n>>> anchor-matched cell: K=%d gen=%d (ld_rmse=%.4f)\n", best$K, best$gen, best$ld_rmse))

cmp <- merge(emp_bin, best$bins, by="dibin", suffixes=c("_emp","_sim"))
cmp[, `:=`(pi_gap=pi_emp-pi_sim, fst_gap=fst_emp-fst_sim,
           sorted_ratio=ifelse(sorted_sim>0, sorted_emp/sorted_sim, NA_real_))]
cat("\n=== DOSE-RESPONSE: gap by DI bin at the anchor-matched cell ===\n")
print(cmp[, .(dibin, pi_emp=round(pi_emp,3), pi_sim=round(pi_sim,3), pi_gap=round(pi_gap,3),
              fst_emp=round(fst_emp,3), fst_sim=round(fst_sim,3), fst_gap=round(fst_gap,3),
              srt_emp=round(sorted_emp,3), srt_sim=round(sorted_sim,3),
              ratio=round(sorted_ratio,1))])
cat("\nREAD: if gaps are ~0 in the LOW-DI bins and grow toward high DI, the demography\n")
cat("is right and the sorting excess is real. If gaps are FLAT across DI, the\n")
cat("mismatch is demographic and the sorting excess is not interpretable.\n")

pl <- melt(cmp[, .(dibin, pi_gap, fst_gap)], id.vars="dibin")
p <- ggplot(pl, aes(dibin, value, group=variable, colour=variable)) +
  geom_line() + geom_point() + geom_hline(yintercept=0, linetype=2) +
  labs(x="DI bin (low -> high ancestry-informativeness)", y="empirical - simulated",
       title="Module E: gap vs DI at the low-DI-anchored neutral demography") +
  theme_bw(base_size=11) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(FIGDIR,"moduleE_di_stratified.pdf"), p, width=9, height=5)
cat("\nsaved: ", OUTRDS, "\n")
