## =============================================================================
## Module E -- LD-decay / pi / F_ST / sorting under COLONY (nest) sampling
## =============================================================================
## Completes the colony-sampling analysis. Previously LD-decay was compared with
## the sim sampled as RANDOM unrelated females while the empirical data are
## workers from a few colonies. Sibling structure inflates composite r^2, so the
## sim's LD must be generated under the same sampling scheme before any
## LD-matched demography can be identified.
##
## Scans founding x K x generation, colony-sampling each deme with its matched
## empirical nest structure, and reports the cell where LD-decay, pi and F_ST all
## match the empirical values simultaneously. The sorting excess read off THAT
## cell is the number worth reporting.
##
## Output: data/moduleE_sim/moduleE_colony_ld.rds
##         Figures/moduleE_colony_ld.pdf
## =============================================================================

suppressMessages({library(data.table); library(ggplot2); library(parallel)})
source("dev/R/parallelism_stats.R")

FORMICA  <- "/Users/petrikem/gitlab/formica_hybrid"
SCREEN   <- file.path(FORMICA, "data/moduleE_sim/screen")
KSWEEP   <- file.path(FORMICA, "data/moduleE_sim/ksweep")
POOL     <- file.path(FORMICA, "data/moduleE_founders/moduleE_founder_haplotypes.rds")
FVCF_DIR <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/SLiM/founders_DIm25_thin15000"
NESTCSV  <- file.path(FORMICA, "data/moduleE_sim/empirical_nest_structure.csv")
EMPTGT   <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/moduleE/inputs/empirical_targets.rds"
OUTRDS   <- file.path(FORMICA, "data/moduleE_sim/moduleE_colony_ld.rds")
FIGDIR   <- file.path(FORMICA, "Figures")

FOUNDINGS <- list(list(lab="12/12", tag="Naq12_Npol12"), list(lab="30/13", tag="Naq30_Npol13"))
KS   <- c(200, 375, 750, 1500)
GENS <- c(40, 80, 120, 160)
FIX_TH <- 0.15; SORT_TH <- 0.5; MIN_PARENT_MAF <- 0.15; NULL_PROB <- 0.5
MATE_W <- c(49,38,15,7,3,1)
CORES  <- 8L; SEED <- 1L
dir.create(FIGDIR, showWarnings = FALSE)

## ---- empirical targets (already nest-structured by nature) -----------------
emp <- readRDS(EMPTGT)
emp_ld <- as.data.table(emp$emp_ld); emp_pi <- emp$emp_pi; emp_fst <- emp$emp_fst
MAXDIST <- emp$params$MAXDIST; MAXPART <- emp$params$MAXPART
DIST_BREAKS <- emp$params$DIST_BREAKS; MIN_MAF <- emp$params$MIN_MAF
emp_vec <- setNames(emp_ld$r2, emp_ld$dbin)
EMP_SORTED <- 0.136                      # from the matched-conventions run
message(sprintf("empirical: pi=%.3f F_ST=%.3f sorted=%.3f", emp_pi, emp_fst, EMP_SORTED))

## ---- markers, cM map, parental reference -----------------------------------
sim_markers <- unique(unlist(lapply(
  list.files(FVCF_DIR, "founders_ch.*[.]vcf$", full.names=TRUE),
  function(f) fread(f, skip="#CHROM", select=3, header=TRUE)[[1]])))
ph <- readRDS(POOL); pmap <- as.data.table(ph$map)
pj <- match(sim_markers, pmap$marker); ok <- !is.na(pj)
markers <- sim_markers[ok]; pj <- pj[ok]
mk_chr <- pmap$Chr[pj]; mk_cM <- pmap$cM[pj]; mk_pos <- pmap$Pos[pj]
chr_idx <- split(seq_along(markers), mk_chr)
chr_len <- sapply(chr_idx, function(i) diff(range(mk_cM[i])))
f_aq_par <- colMeans(ph$aqu[, pj, drop=FALSE]); f_pol_par <- colMeans(ph$pol[, pj, drop=FALSE])
p_pooled <- (colSums(ph$aqu[, pj, drop=FALSE]) + colSums(ph$pol[, pj, drop=FALSE])) /
            (nrow(ph$aqu) + nrow(ph$pol))
parent_maf <- pmin(p_pooled, 1-p_pooled); names(parent_maf) <- markers
nest <- fread(NESTCSV)
nsz  <- lapply(nest$nest_sizes, function(s) as.integer(strsplit(s, ",")[[1]]))

## ---- helpers ---------------------------------------------------------------
meiosis <- function(h1, h2) {
  out <- h1
  for (ch in names(chr_idx)) {
    i <- chr_idx[[ch]]; cm <- mk_cM[i]; L <- chr_len[[ch]]
    if (!is.finite(L) || L <= 0) L <- 0
    nx <- rpois(1, L/100)
    xs <- if (nx > 0) sort(runif(nx, min(cm), max(cm))) else numeric(0)
    par <- (findInterval(cm, xs) + sample(0:1, 1)) %% 2
    out[i] <- ifelse(par == 0, h1[i], h2[i])
  }
  out
}
read_queens <- function(vcf) {
  dt <- fread(vcf, skip="#CHROM", header=TRUE, sep="\t", showProgress=FALSE)
  mk <- paste(sub("ch","Chr",dt[[1]]), dt[[2]], sep=":")
  G <- as.matrix(dt[, 10:ncol(dt)]); gt <- sub(":.*$","",G)
  H1 <- matrix(suppressWarnings(as.integer(sub("[|/].*$","",gt))), nrow=nrow(G))
  H2 <- matrix(suppressWarnings(as.integer(sub("^.*[|/]","",gt))), nrow=nrow(G))
  m <- match(markers, mk); keep <- !is.na(m)
  f <- function(H){ M <- matrix(0L, length(markers), ncol(H)); M[keep,] <- H[m[keep],,drop=FALSE]; M }
  list(h1=f(H1), h2=f(H2))
}
colony_sample <- function(q, nest_sizes) {
  nq <- ncol(q$h1); if (nq < 2) return(NULL)
  qsel <- sample(nq, min(length(nest_sizes), nq)); cols <- list()
  for (j in seq_along(qsel)) {
    qi <- qsel[j]; k <- nest_sizes[j]
    nm  <- sample(seq_along(MATE_W), 1, prob=MATE_W)
    src <- sample(setdiff(seq_len(nq), qi), min(nm, nq-1))
    sires <- lapply(src, function(s) meiosis(q$h1[,s], q$h2[,s]))
    for (w in seq_len(k))
      cols[[length(cols)+1L]] <- meiosis(q$h1[,qi], q$h2[,qi]) + sires[[sample(length(sires),1)]]
  }
  do.call(cbind, cols)
}
ld_decay_curve <- function(dose) {
  out <- list()
  for (ch in names(chr_idx)) {
    i <- chr_idx[[ch]]; o <- order(mk_pos[i]); i <- i[o]
    pos <- mk_pos[i]; G <- dose[i,,drop=FALSE]
    p <- rowMeans(G, na.rm=TRUE)/2
    jj <- which(is.finite(p) & pmin(p,1-p) >= MIN_MAF); if (length(jj) < 2L) next
    pos <- pos[jj]; G <- G[jj,,drop=FALSE]
    for (a in seq_len(length(pos)-1L)) {
      b <- (a+1L):min(a+MAXPART, length(pos)); d <- pos[b]-pos[a]
      keep <- d <= MAXDIST; if (!any(keep)) next; b <- b[keep]; d <- d[keep]
      r <- suppressWarnings(cor(G[a,], t(G[b,,drop=FALSE]), use="pairwise.complete.obs"))
      out[[length(out)+1L]] <- data.table(dist=d, r2=as.numeric(r)^2)
    }
  }
  if (!length(out)) return(NULL)
  res <- rbindlist(out)[is.finite(r2)]
  res[, dbin := cut(dist, DIST_BREAKS, labels=FALSE)]
  res[!is.na(dbin), .(r2=mean(r2)), by=dbin]
}
run_sort <- function(pop_means) {
  pm <- rbind(pop_means, aquilonia_parent=f_aq_par*2, polyctena_parent=f_pol_par*2)
  colnames(pm) <- markers
  parallelism_stats(list(pop_means=pm), hybrid_pops=rownames(pop_means),
    aqu_pops="aquilonia_parent", pol_pops="polyctena_parent",
    fix_th=FIX_TH, sort_th=SORT_TH, null_prob=NULL_PROB,
    parent_maf=parent_maf, min_parent_maf=MIN_PARENT_MAF)
}

## ---- one cell --------------------------------------------------------------
do_cell <- function(lab, tag, K, gen) {
  dir <- if (K == 6250) SCREEN else KSWEEP
  ttag <- if (K == 6250) tag else paste0(tag, "_K", K)
  vcfs <- list.files(dir, sprintf("hyb_neutral_realfounders_%s_rep[0-9]+_gen%d[.]vcf$", ttag, gen),
                     full.names=TRUE)
  if (length(vcfs) < 5) return(NULL)
  pm <- matrix(NA_real_, length(vcfs), length(markers))
  curves <- list()
  for (i in seq_along(vcfs)) {
    set.seed(SEED + 1000*K + gen + i)
    q <- read_queens(vcfs[i])
    W <- colony_sample(q, nsz[[(i-1) %% length(nsz) + 1]])
    if (is.null(W)) next
    pm[i, ] <- rowMeans(W)
    cv <- ld_decay_curve(W); if (!is.null(cv)) curves[[length(curves)+1L]] <- cv
  }
  keep <- !is.na(pm[,1]); pm <- pm[keep,,drop=FALSE]
  if (nrow(pm) < 5 || !length(curves)) return(NULL)
  rownames(pm) <- paste0("deme", seq_len(nrow(pm))); colnames(pm) <- markers
  ld <- rbindlist(curves)[, .(r2=mean(r2)), by=dbin][order(dbin)]
  m  <- match(ld$dbin, as.integer(names(emp_vec)))
  s  <- run_sort(pm); tt <- s[differentiated==TRUE, .N, by=sort_class]
  p  <- pm/2
  data.table(founding=lab, K=K, gen=gen, n_demes=nrow(pm),
    ld_rmse = sqrt(mean((ld$r2 - emp_vec[m])^2, na.rm=TRUE)),
    pi  = mean(rowMeans(2*p*(1-p), na.rm=TRUE), na.rm=TRUE),
    fst = { pb<-colMeans(p,na.rm=TRUE); Hs<-colMeans(2*p*(1-p),na.rm=TRUE); Ht<-2*pb*(1-pb)
            sum(Ht-Hs,na.rm=TRUE)/sum(Ht,na.rm=TRUE) },
    sorted = 1 - tt[sort_class=="unsorted", N]/sum(tt$N))
}

cells <- CJ(fi=seq_along(FOUNDINGS), K=KS, gen=GENS, sorted=FALSE)
message(sprintf("[run] %d cells, colony sampling + LD ...", nrow(cells)))
res <- rbindlist(mcMap(function(fi,K,gen) do_cell(FOUNDINGS[[fi]]$lab, FOUNDINGS[[fi]]$tag, K, gen),
                       cells$fi, cells$K, cells$gen, mc.cores=CORES), fill=TRUE)
res[, `:=`(pi_gap = pi - emp_pi, fst_gap = fst - emp_fst, sorted_ratio = EMP_SORTED/sorted)]
saveRDS(list(res=res, emp_pi=emp_pi, emp_fst=emp_fst, emp_sorted=EMP_SORTED, emp_ld=emp_ld), OUTRDS)

cat(sprintf("\n=== EMPIRICAL: pi=%.3f  F_ST=%.3f  sorted=%.3f ===\n", emp_pi, emp_fst, EMP_SORTED))
cat("\n=== colony-sampled cells, ranked by LD-decay RMSE ===\n")
print(res[order(ld_rmse), .(founding, K, gen, n=n_demes, ld_rmse=round(ld_rmse,4),
     pi=round(pi,3), pi_gap=round(pi_gap,3), fst=round(fst,3), fst_gap=round(fst_gap,3),
     sorted=round(sorted,3), excess=round(sorted_ratio,1))][1:15])
cat("\n=== cells matching LD (rmse<0.04) AND pi (|gap|<0.04) AND F_ST (|gap|<0.05) ===\n")
good <- res[ld_rmse < 0.04 & abs(pi_gap) < 0.04 & abs(fst_gap) < 0.05][order(ld_rmse)]
if (nrow(good)) {
  print(good[, .(founding, K, gen, ld_rmse=round(ld_rmse,4), pi=round(pi,3),
                 fst=round(fst,3), sorted=round(sorted,3), excess=round(sorted_ratio,1))])
} else {
  cat("(none -- no cell matches all three simultaneously)\n")
}
cat("\nsaved: ", OUTRDS, "\n")
