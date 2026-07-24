## =============================================================================
## Module E -- MIXED-AGE null: do populations of differing age reconcile
##             LD-decay with pi / F_ST?
## =============================================================================
## No single-age neutral cell matches all three summaries: LD-decay wants little
## drift (K >= 1500, gen 40-60), while pi and F_ST want K ~ 100-200. Small Ne
## generates long-range LD by drift and flattens the decay curve, so the two
## demands are irreconcilable at ONE age.
##
## The empirical populations are NOT one age: within-population heterozygosity
## spans 0.225 (Parikkala) to 0.419 (Sielva), which a single-age model cannot
## produce. A MIXTURE may satisfy both demands at once:
##    young demes  -> retain high, shallow LD  (matches the LD-decay target)
##    old demes    -> supply drift             (low pi, high F_ST)
##
## Tested with the runs we ALREADY have: every deme was sampled at generations
## 20,40,...,160, so a mixed-age null is just a different choice of which
## generation to read per deme. No new simulations.
##
## Sampling: random, subsampled to the empirical per-population sizes (the nest
## effect is small -- measured within-nest excess similarity is only 0.047,
## consistent with polygyny -- so random sampling is the appropriate baseline).
##
## Output: data/moduleE_sim/moduleE_mixed_age.rds
## =============================================================================

suppressMessages({library(data.table); library(parallel)})
source("dev/R/parallelism_stats.R")

FORMICA  <- "/Users/petrikem/gitlab/formica_hybrid"
SCREEN   <- file.path(FORMICA, "data/moduleE_sim/screen")
KSWEEP   <- file.path(FORMICA, "data/moduleE_sim/ksweep")
POOL     <- file.path(FORMICA, "data/moduleE_founders/moduleE_founder_haplotypes.rds")
FVCF_DIR <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/SLiM/founders_DIm25_thin15000"
EMPTGT   <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/moduleE/inputs/empirical_targets.rds"
OUTRDS   <- file.path(FORMICA, "data/moduleE_sim/moduleE_mixed_age.rds")

FOUNDINGS <- list(list(lab="12/12", tag="Naq12_Npol12"), list(lab="30/13", tag="Naq30_Npol13"))
KS <- c(375, 750, 1500)
FIX_TH <- 0.15; SORT_TH <- 0.5; MIN_PARENT_MAF <- 0.15; NULL_PROB <- 0.5
CORES <- 8L; MATCH_SEED <- 1L; EMP_SORTED <- 0.136

## age scenarios: generation assigned to each of the 20 demes
SCEN <- list(
  "single gen40"       = rep(40, 20),
  "single gen120"      = rep(120, 20),
  "uniform 20-160"     = rep(c(20,40,60,80,100,120,140,160), length.out = 20),
  "bimodal 20 / 160"   = c(rep(20,10), rep(160,10)),
  "young-skew 20-60"   = c(rep(20,7), rep(40,7), rep(60,6)),
  "old-skew 120-160"   = c(rep(120,7), rep(140,7), rep(160,6))
)

## ---- empirical targets ----
emp <- readRDS(EMPTGT)
emp_ld <- as.data.table(emp$emp_ld); emp_pi <- emp$emp_pi; emp_fst <- emp$emp_fst
emp_ns <- emp$emp_ns
MAXDIST <- emp$params$MAXDIST; MAXPART <- emp$params$MAXPART
DIST_BREAKS <- emp$params$DIST_BREAKS; MIN_MAF <- emp$params$MIN_MAF
emp_vec <- setNames(emp_ld$r2, emp_ld$dbin)

## ---- markers / parents ----
sim_markers <- unique(unlist(lapply(
  list.files(FVCF_DIR, "founders_ch.*[.]vcf$", full.names=TRUE),
  function(f) fread(f, skip="#CHROM", select=3, header=TRUE)[[1]])))
ph <- readRDS(POOL); pmap <- as.data.table(ph$map)
pj <- match(sim_markers, pmap$marker); ok <- !is.na(pj)
markers <- sim_markers[ok]; pj <- pj[ok]
mk_chr <- pmap$Chr[pj]; mk_pos <- pmap$Pos[pj]
chr_idx <- split(seq_along(markers), mk_chr)
f_aq <- colMeans(ph$aqu[,pj,drop=FALSE]); f_pol <- colMeans(ph$pol[,pj,drop=FALSE])
p0 <- (colSums(ph$aqu[,pj,drop=FALSE])+colSums(ph$pol[,pj,drop=FALSE]))/(nrow(ph$aqu)+nrow(ph$pol))
parent_maf <- pmin(p0,1-p0); names(parent_maf) <- markers

## ---- helpers ----
read_dosage <- function(vcf) {
  dt <- fread(vcf, skip="#CHROM", header=TRUE, sep="\t", showProgress=FALSE)
  mk <- paste(sub("ch","Chr",dt[[1]]), dt[[2]], sep=":")
  G <- as.matrix(dt[, 10:ncol(dt)]); gt <- sub(":.*$","",G)
  d <- matrix(suppressWarnings(as.integer(sub("[|/].*$","",gt))) +
              suppressWarnings(as.integer(sub("^.*[|/]","",gt))), nrow=nrow(G))
  m <- match(markers, mk); keep <- !is.na(m)
  D <- matrix(0L, length(markers), ncol(d)); D[keep,] <- d[m[keep],,drop=FALSE]; D
}
ld_curve <- function(dose) {
  out <- list()
  for (ch in names(chr_idx)) {
    i <- chr_idx[[ch]]; o <- order(mk_pos[i]); i <- i[o]
    pos <- mk_pos[i]; G <- dose[i,,drop=FALSE]
    p <- rowMeans(G, na.rm=TRUE)/2
    jj <- which(is.finite(p) & pmin(p,1-p) >= MIN_MAF); if (length(jj)<2L) next
    pos <- pos[jj]; G <- G[jj,,drop=FALSE]
    for (a in seq_len(length(pos)-1L)) {
      b <- (a+1L):min(a+MAXPART,length(pos)); d <- pos[b]-pos[a]
      k <- d <= MAXDIST; if (!any(k)) next; b <- b[k]; d <- d[k]
      r <- suppressWarnings(cor(G[a,], t(G[b,,drop=FALSE]), use="pairwise.complete.obs"))
      out[[length(out)+1L]] <- data.table(dist=d, r2=as.numeric(r)^2)
    }
  }
  if (!length(out)) return(NULL)
  res <- rbindlist(out)[is.finite(r2)]
  res[, dbin := cut(dist, DIST_BREAKS, labels=FALSE)]
  res[!is.na(dbin), .(r2=mean(r2)), by=dbin]
}
run_sort <- function(mm) {
  pm <- rbind(mm, aquilonia_parent=f_aq*2, polyctena_parent=f_pol*2); colnames(pm) <- markers
  parallelism_stats(list(pop_means=pm), hybrid_pops=rownames(mm),
    aqu_pops="aquilonia_parent", pol_pops="polyctena_parent", fix_th=FIX_TH,
    sort_th=SORT_TH, null_prob=NULL_PROB, parent_maf=parent_maf, min_parent_maf=MIN_PARENT_MAF)
}

do_cell <- function(lab, tag, K, sname) {
  gens <- SCEN[[sname]]
  dir  <- if (K == 6250) SCREEN else KSWEEP
  ttag <- if (K == 6250) tag else paste0(tag, "_K", K)
  pm <- matrix(NA_real_, 20, length(markers)); curves <- list()
  for (i in 1:20) {
    f <- file.path(dir, sprintf("hyb_neutral_realfounders_%s_rep%d_gen%d.vcf", ttag, i, gens[i]))
    if (!file.exists(f)) next
    D <- read_dosage(f)
    tgt <- emp_ns[(i-1) %% length(emp_ns) + 1]; ni <- ncol(D)
    if (ni > tgt) { set.seed(MATCH_SEED+i); D <- D[, sample(ni, tgt), drop=FALSE] }
    pm[i,] <- rowMeans(D, na.rm=TRUE)
    cv <- ld_curve(D); if (!is.null(cv)) curves[[length(curves)+1L]] <- cv
  }
  keep <- !is.na(pm[,1]); pm <- pm[keep,,drop=FALSE]
  if (nrow(pm) < 10 || !length(curves)) return(NULL)
  rownames(pm) <- paste0("deme", seq_len(nrow(pm))); colnames(pm) <- markers
  ld <- rbindlist(curves)[, .(r2=mean(r2)), by=dbin][order(dbin)]
  m  <- match(ld$dbin, as.integer(names(emp_vec)))
  s  <- run_sort(pm); tt <- s[differentiated==TRUE, .N, by=sort_class]
  p  <- pm/2; pb <- colMeans(p,na.rm=TRUE); Hs <- colMeans(2*p*(1-p),na.rm=TRUE); Ht <- 2*pb*(1-pb)
  data.table(founding=lab, K=K, scenario=sname, n_demes=nrow(pm),
    ld_rmse=sqrt(mean((ld$r2-emp_vec[m])^2, na.rm=TRUE)),
    pi=mean(rowMeans(2*p*(1-p), na.rm=TRUE), na.rm=TRUE),
    fst=sum(Ht-Hs,na.rm=TRUE)/sum(Ht,na.rm=TRUE),
    sorted=1-tt[sort_class=="unsorted",N]/sum(tt$N))
}

grid <- CJ(fi=seq_along(FOUNDINGS), K=KS, sn=names(SCEN), sorted=FALSE)
message(sprintf("[run] %d cells ...", nrow(grid)))
res <- rbindlist(mcMap(function(fi,K,sn) do_cell(FOUNDINGS[[fi]]$lab, FOUNDINGS[[fi]]$tag, K, sn),
                       grid$fi, grid$K, grid$sn, mc.cores=CORES), fill=TRUE)
res[, `:=`(pi_gap=pi-emp_pi, fst_gap=fst-emp_fst, excess=EMP_SORTED/sorted)]
saveRDS(list(res=res, emp_pi=emp_pi, emp_fst=emp_fst, emp_sorted=EMP_SORTED), OUTRDS)

cat(sprintf("\n=== EMPIRICAL: pi=%.3f  F_ST=%.3f  sorted=%.3f ===\n", emp_pi, emp_fst, EMP_SORTED))
cat("\n=== all cells, ranked by LD-decay RMSE ===\n")
print(res[order(ld_rmse), .(founding, K, scenario, ld=round(ld_rmse,4), pi=round(pi,3),
      pi_gap=round(pi_gap,3), fst=round(fst,3), fst_gap=round(fst_gap,3),
      sorted=round(sorted,3), excess=round(excess,1))])
cat("\n=== cells matching LD (<0.04) AND pi (|gap|<0.04) AND F_ST (|gap|<0.05) ===\n")
g <- res[ld_rmse<0.04 & abs(pi_gap)<0.04 & abs(fst_gap)<0.05][order(ld_rmse)]
if (nrow(g)) { print(g[, .(founding,K,scenario,ld=round(ld_rmse,4),pi=round(pi,3),
                           fst=round(fst,3),sorted=round(sorted,3),excess=round(excess,1))])
} else { cat("(none)\n") }
cat("\nsaved: ", OUTRDS, "\n")
