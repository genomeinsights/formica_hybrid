## =============================================================================
## Module E -- founder-number (shallow-bottleneck) sweep: cross-founding trend
## =============================================================================
## Reads the run_foundersweep.sh output and reports, for each founding setting and
## generation, the NEUTRAL-bin (DI<=-90) F_ST and pi and the HIGH-DI (>-15) F_ST /
## sorted fraction, all sample-size matched to the empirical per-pop sizes and using
## the locked sorting conventions -- so we can see whether a shallower founding
## bottleneck moves neutral F_ST -> ~0 and pi -> empirical without touching the
## diagnostic signal. Same estimators/logic as analyze_di_stratified.R.
##
## Usage: Rscript compare_foundersweep.R [sim_output_dir] [founder_vcf_dir]
##   defaults: ../output_foundersweep  ../founders/maf015_DIstrat1500_x40
## =============================================================================
suppressMessages({ library(data.table); library(parallel) })
.args <- commandArgs(trailingOnly = FALSE); .self <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT <- normalizePath(file.path(dirname(.self), ".."))
pos  <- commandArgs(trailingOnly = TRUE)
SIMDIR   <- if (length(pos)>=1) pos[1] else file.path(ROOT, "output_foundersweep")
FVCF_DIR <- if (length(pos)>=2) pos[2] else file.path(ROOT, "founders/maf015_DIstrat1500_x40")
BUNDLE   <- file.path(ROOT, "inputs/empirical_bundle.rds")
source(file.path(ROOT, "scripts/parallelism_stats.R"))
DI_BREAKS <- c(-Inf,-90,-75,-60,-50,-40,-30,-25,-20,-15,Inf)
FIX_TH<-0.15; SORT_TH<-0.5; MIN_PARENT_MAF<-0.15; NULL_PROB<-0.5; MATCH_SEED<-1L
CORES <- max(1L, detectCores()-1L)

b <- readRDS(BUNDLE); uni <- as.data.table(b$universe); setkey(uni, marker)
emp_ns <- b$emp_ns; pops <- b$pops
sim_markers <- unique(unlist(lapply(list.files(FVCF_DIR, "founders_ch.*[.]vcf$", full.names=TRUE),
  function(f) fread(f, skip="#CHROM", select=3, header=TRUE)[[1]])))
mk <- intersect(sim_markers, uni$marker); um <- uni[.(mk)]; setorder(um, Chr, Pos)
markers <- um$marker; DI <- um$DI; dibin <- cut(DI, DI_BREAKS)
f_aq <- b$f_aq_par[markers]; f_pol <- b$f_pol_par[markers]
parent_maf <- setNames(um$parent_maf, markers)
emp_M <- (b$emp_mean[, markers, drop=FALSE]) / 1000
neu <- levels(dibin)[1]; hi <- levels(dibin)[nlevels(dibin)]
message(sprintf("markers: %d ; neutral(DI<=-90): %d ; high(DI>-15): %d",
                length(markers), sum(dibin==neu), sum(dibin==hi)))

read_dosage <- function(vcf) {
  dt <- fread(vcf, skip="#CHROM", header=TRUE, sep="\t", showProgress=FALSE)
  m0 <- paste(sub("ch","Chr",dt[[1]]), dt[[2]], sep=":")
  G <- as.matrix(dt[,10:ncol(dt)]); gt <- sub(":.*$","",G)
  d <- matrix(suppressWarnings(as.integer(sub("[|/].*$","",gt))) +
              suppressWarnings(as.integer(sub("^.*[|/]","",gt))), nrow=nrow(G))
  mm <- match(markers, m0); keep <- !is.na(mm)
  D <- matrix(0L, length(markers), ncol(d)); D[keep,] <- d[mm[keep],,drop=FALSE]; D
}
run_sort <- function(mm) {
  pm <- rbind(mm, aquilonia_parent=f_aq*2, polyctena_parent=f_pol*2); colnames(pm) <- markers
  parallelism_stats(list(pop_means=pm), hybrid_pops=rownames(mm),
    aqu_pops="aquilonia_parent", pol_pops="polyctena_parent", fix_th=FIX_TH,
    sort_th=SORT_TH, null_prob=NULL_PROB, parent_maf=parent_maf, min_parent_maf=MIN_PARENT_MAF)
}
by_bin <- function(mm) {
  s <- run_sort(mm); sc <- s[match(markers, marker), sort_class]
  p <- mm/2; pb <- colMeans(p,na.rm=TRUE); Hs <- colMeans(2*p*(1-p),na.rm=TRUE); Ht <- 2*pb*(1-pb)
  data.table(dibin=dibin, Hs=Hs, Ht=Ht, sorted=!is.na(sc)&sc!="unsorted", tested=!is.na(sc))[
    , .(pi=mean(Hs,na.rm=TRUE), fst=sum(Ht-Hs,na.rm=TRUE)/sum(Ht,na.rm=TRUE),
        sorted=sum(sorted)/sum(tested)), by=dibin][order(dibin)]
}
emp_bin <- by_bin(emp_M)

## founding tags present, sorted by founder number
tags <- unique(na.omit(sub(".*(Naq[0-9]+_Npol[0-9]+)_K.*", "\\1", list.files(SIMDIR, "\\.vcf$"))))
tags <- tags[order(as.integer(sub("Naq([0-9]+)_.*","\\1",tags)))]
gens <- sort(as.integer(unique(na.omit(sub(".*_gen([0-9]+)\\.vcf$","\\1", list.files(SIMDIR,"\\.vcf$"))))))
Ks   <- sort(as.integer(unique(na.omit(sub(".*_K([0-9]+)_rep.*","\\1", list.files(SIMDIR,"\\.vcf$"))))))
message("foundings: ", paste(tags,collapse=", "), " ; K: ", paste(Ks,collapse=","), " ; gens: ", paste(gens,collapse=","))

do_cell <- function(TAG, K, gen) {
  pm <- matrix(NA_real_, 20, length(markers))
  for (i in 1:20) {
    f <- file.path(SIMDIR, sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf", TAG, K, i, gen))
    if (!file.exists(f)) next
    D <- read_dosage(f); tgt <- emp_ns[(i-1)%%length(emp_ns)+1]
    if (ncol(D) > tgt) { set.seed(MATCH_SEED+i); D <- D[, sample(ncol(D), tgt), drop=FALSE] }
    pm[i,] <- rowMeans(D, na.rm=TRUE)
  }
  keep <- !is.na(pm[,1]); pm <- pm[keep,,drop=FALSE]
  if (nrow(pm) < 10) return(NULL)
  rownames(pm) <- paste0("deme",seq_len(nrow(pm))); colnames(pm) <- markers
  bb <- by_bin(pm)
  data.table(founding=TAG, N_AQ=as.integer(sub("Naq([0-9]+)_.*","\\1",TAG)), K=K, gen=gen,
             neu_fst=bb[dibin==neu,fst], neu_pi=bb[dibin==neu,pi], neu_sorted=bb[dibin==neu,sorted],
             hi_fst=bb[dibin==hi,fst],  hi_sorted=bb[dibin==hi,sorted])
}
grid <- CJ(TAG=tags, K=Ks, gen=gens, sorted=FALSE)
res <- rbindlist(mcMap(function(t,k,g) do_cell(t,k,g), grid$TAG, grid$K, grid$gen, mc.cores=CORES))
setorder(res, N_AQ, gen)
saveRDS(list(res=res, emp_bin=emp_bin, neu=neu, hi=hi), file.path(SIMDIR, "foundersweep_trend.rds"))

ef <- emp_bin[dibin==neu, fst]; ep <- emp_bin[dibin==neu, pi]; ehf <- emp_bin[dibin==hi, fst]
cat(sprintf("\n=== EMPIRICAL targets: neutral F_ST=%.3f pi=%.3f | high-DI F_ST=%.3f sorted=%.3f ===\n",
            ef, ep, ehf, emp_bin[dibin==hi, sorted]))
cat("\n=== founder-number sweep (K=", Ks[1], ", per generation) ===\n", sep="")
pr <- res[, .(founding, gen, neu_fst=round(neu_fst,3), neu_pi=round(neu_pi,3),
              neu_sorted=round(neu_sorted,3), hi_fst=round(hi_fst,3), hi_sorted=round(hi_sorted,3))]
print(pr)
cat(sprintf("\n(neutral targets to hit: F_ST->%.3f, pi->%.3f. Shallower founding = larger N_AQ.)\n", ef, ep))
