## =============================================================================
## Module E -- DI-stratified analysis (PORTABLE)
## =============================================================================
## Self-contained version for the migrated machine. Resolves paths relative to
## the moduleE_slim/ bundle and reads the compact empirical_bundle.rds instead of
## the 1 GB hybrids_only_maf005.Rdata.
##
## Computes observed vs neutral-null pi, F_ST and sorted-locus fraction across DI
## bins, and picks the neutral cell that best matches the LOW-DI (neutral) bin --
## the internal control that cannot be driven by ancestry sorting. The high-DI
## excess is read at that cell. See ../HANDOFF.md for design, current result and
## the next step (sweep K upward + founder number).
##
## The empirical LD-decay anchor (optional, for a secondary check) uses the
## per-individual genotypes stored in the bundle ($emp_geno_anchor); it is NOT
## used to choose the cell here because F_ST/pi on the neutral bin are the
## trustworthy anchors (LD favours the opposite, over-drifted, direction).
##
## Usage:  Rscript analyze_di_stratified.R [sim_output_dir] [founding_tag]
##   sim_output_dir : where the SLiM VCFs are (default ../output)
##   founding_tag   : e.g. Naq30_Npol13 (default). K values auto-detected.
## =============================================================================

suppressMessages({library(data.table); library(parallel)})

.args <- commandArgs(trailingOnly = FALSE)
.self <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT  <- normalizePath(file.path(dirname(.self), ".."))
pos   <- commandArgs(trailingOnly = TRUE)
SIMDIR   <- if (length(pos) >= 1) pos[1] else file.path(ROOT, "output")
TAG      <- if (length(pos) >= 2) pos[2] else "Naq30_Npol13"

BUNDLE   <- file.path(ROOT, "inputs/empirical_bundle.rds")
FVCF_DIR <- file.path(ROOT, "founders/maf015_DIstrat4000")
source(file.path(ROOT, "scripts/parallelism_stats.R"))

DI_BREAKS <- c(-Inf,-90,-75,-60,-50,-40,-30,-25,-20,-15,Inf)
FIX_TH <- 0.15; SORT_TH <- 0.5; MIN_PARENT_MAF <- 0.15; NULL_PROB <- 0.5
CORES <- max(1L, detectCores() - 1L); MATCH_SEED <- 1L

## ---- empirical from the compact bundle -------------------------------------
b <- readRDS(BUNDLE)
uni <- as.data.table(b$universe); setkey(uni, marker)
emp_ns <- b$emp_ns; pops <- b$pops

sim_markers <- unique(unlist(lapply(
  list.files(FVCF_DIR, "founders_ch.*[.]vcf$", full.names=TRUE),
  function(f) fread(f, skip="#CHROM", select=3, header=TRUE)[[1]])))
mk <- intersect(sim_markers, uni$marker)
um <- uni[.(mk)]; setorder(um, Chr, Pos)
markers <- um$marker; DI <- um$DI
dibin <- cut(DI, DI_BREAKS)
f_aq <- b$f_aq_par[markers]; f_pol <- b$f_pol_par[markers]
parent_maf <- setNames(um$parent_maf, markers)
emp_M <- (b$emp_mean[, markers, drop=FALSE]) / 1000            # pops x markers dosage
message(sprintf("markers analysed: %d ; DI bins: %d ; pops: %d",
                length(markers), nlevels(dibin), length(pops)))

## ---- statistics ------------------------------------------------------------
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
neu <- levels(dibin)[1]; hi <- levels(dibin)[length(levels(dibin))]

## ---- simulated cells: (K, generation), sample-size matched -----------------
Ks   <- as.integer(unique(na.omit(sub(sprintf(".*_%s_K([0-9]+)_.*", TAG), "\\1",
          grep(sprintf("_%s_K", TAG), list.files(SIMDIR, "\\.vcf$"), value=TRUE)))))
gens <- as.integer(unique(na.omit(sub(".*_gen([0-9]+)\\.vcf$", "\\1",
          list.files(SIMDIR, "\\.vcf$")))))
if (!length(Ks)) stop("no VCFs matching tag ", TAG, " in ", SIMDIR)
message("K present: ", paste(sort(Ks), collapse=","), " ; gens: ", paste(sort(gens), collapse=","))

do_cell <- function(K, gen) {
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
  bins <- by_bin(pm)
  list(K=K, gen=gen, bins=bins,
       neu_fst=bins[dibin==neu,fst], neu_pi=bins[dibin==neu,pi])
}
cells <- CJ(K=Ks, gen=gens, sorted=FALSE)
out <- mcMap(function(K,g) do_cell(K,g), cells$K, cells$gen, mc.cores=CORES)
out <- out[!vapply(out, is.null, logical(1))]

## anchor on the NEUTRAL bin: minimise |F_ST gap| + |pi gap| there
ef <- emp_bin[dibin==neu, fst]; ep <- emp_bin[dibin==neu, pi]
score <- sapply(out, function(z) abs(z$neu_fst-ef) + abs(z$neu_pi-ep))
best <- out[[which.min(score)]]
saveRDS(list(cells=out, emp_bin=emp_bin, best=best, TAG=TAG),
        file.path(SIMDIR, "di_stratified_results.rds"))

cat(sprintf("\n=== empirical: neutral bin F_ST=%.3f pi=%.3f | high bin F_ST=%.3f ===\n",
            ef, ep, emp_bin[dibin==hi, fst]))
cat("\n=== neutral-bin fit across cells (best anchors the demography) ===\n")
r <- rbindlist(lapply(out, function(z) data.table(K=z$K, gen=z$gen,
      neu_fst=z$neu_fst, neu_pi=z$neu_pi)))[order(abs(neu_fst-ef)+abs(neu_pi-ep))]
print(r[1:min(10,.N), .(K, gen, neu_fst=round(neu_fst,3), neu_fst_gap=round(neu_fst-ef,3),
                        neu_pi=round(neu_pi,3), neu_pi_gap=round(neu_pi-ep,3))])

cmp <- merge(emp_bin, best$bins, by="dibin", suffixes=c("_emp","_sim"))[order(dibin)]
cmp[, ratio := ifelse(sorted_sim>0, sorted_emp/sorted_sim, NA)]
cat(sprintf("\n>>> anchor-matched cell: K=%d gen=%d\n", best$K, best$gen))
cat("=== DOSE-RESPONSE by DI bin at that cell ===\n")
print(cmp[, .(dibin, fst_emp=round(fst_emp,3), fst_sim=round(fst_sim,3),
              pi_emp=round(pi_emp,3), pi_sim=round(pi_sim,3),
              srt_emp=round(sorted_emp,3), srt_sim=round(sorted_sim,3), ratio=round(ratio,1))])
cat("\nsaved: ", file.path(SIMDIR, "di_stratified_results.rds"), "\n")
