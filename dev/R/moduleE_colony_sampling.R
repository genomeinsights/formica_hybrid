## =============================================================================
## Module E -- COLONY (nest) structured sampling of the neutral null
## =============================================================================
##
## WHY
## ---
## The empirical samples are workers taken from a few COLONIES, not random
## individuals: 164 samples in only 39 nests, median 2 nests/population and 4
## workers/nest, and 9 of 20 populations are a SINGLE colony. Nestmates are
## siblings. We had been sampling the simulation as random unrelated females, so
## the two are not comparable. Sibling sampling
##   - lowers within-sample diversity (pi)        -> empirical pi looked too low
##   - makes each population's allele frequencies family-biased -> F_ST inflated
##   - makes loci look near-fixed                 -> sort_class inflated
## i.e. it can produce all three "excess" signals with no selection at all.
##
## WHY THIS IS DONE IN R RATHER THAN IN SLiM
## -----------------------------------------
## A simulated female IS a queen; the model has no within-colony worker pool, so a
## queen's lifetime surviving daughters number only ~3-4 -- too few to draw 4-10
## sisters from -- and SLiM cannot generate offspring outside a reproduction()
## callback. But `outputIndividualsToVCF` writes PHASED diploid genotypes, so each
## simulated queen hands us her two haplotypes. A worker is exactly
##      (recombinant of the queen's two haplotypes)  +  (a sire haplotype),
## which we can generate faithfully by meiosis along the empirical cM map. This
## reuses every existing run -- no re-simulation -- and reproduces the sibling
## structure exactly.
##
## Sires: haploid Formica males arise from unfertilised eggs, i.e. a male's genome
## is itself a recombinant of some queen's two haplotypes. Sires are therefore
## drawn as recombinants of other queens in the same deme. Queens mate multiply
## (Pamilo 1993); workers are assigned among a queen's sires.
##
## Output: data/moduleE_sim/moduleE_colony_sampling.rds
## =============================================================================

suppressMessages({library(data.table); library(parallel)})
source("dev/R/parallelism_stats.R")

## ------------------------------- CONFIG -------------------------------------
FORMICA  <- "/Users/petrikem/gitlab/formica_hybrid"
SCREEN   <- file.path(FORMICA, "data/moduleE_sim/screen")
KSWEEP   <- file.path(FORMICA, "data/moduleE_sim/ksweep")
POOL     <- file.path(FORMICA, "data/moduleE_founders/moduleE_founder_haplotypes.rds")
FVCF_DIR <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/SLiM/founders_DIm25_thin15000"
NESTCSV  <- file.path(FORMICA, "data/moduleE_sim/empirical_nest_structure.csv")
EMP_RD   <- file.path(FORMICA, "data/hybrids_only_maf005.Rdata")
OUTRDS   <- file.path(FORMICA, "data/moduleE_sim/moduleE_colony_sampling.rds")

## cells to re-sample with colony structure (label, dir, tag, gen)
CELLS <- list(
  list(label="30/13 K6250 gen40", dir=SCREEN, tag="Naq30_Npol13",     gen=40),
  list(label="12/12 K375 gen100", dir=KSWEEP, tag="Naq12_Npol12_K375",gen=100),
  list(label="12/12 K100 gen100", dir=KSWEEP, tag="Naq12_Npol12_K100",gen=100)
)
FIX_TH <- 0.15; SORT_TH <- 0.5; MIN_PARENT_MAF <- 0.15; NULL_PROB <- 0.5
MATE_W <- c(49,38,15,7,3,1)      # Pamilo 1993, conditioned on >=1 mating
SEED   <- 1L
## ---------------------------------------------------------------------------
set.seed(SEED)

## ------------------- markers, cM map, parental reference --------------------
message("[1] markers / cM map / parents ...")
sim_markers <- unique(unlist(lapply(
  list.files(FVCF_DIR, "founders_ch.*[.]vcf$", full.names=TRUE),
  function(f) fread(f, skip="#CHROM", select=3, header=TRUE)[[1]])))
ph <- readRDS(POOL); pmap <- as.data.table(ph$map)
pj <- match(sim_markers, pmap$marker); ok <- !is.na(pj)
markers <- sim_markers[ok]; pj <- pj[ok]
mk_chr <- pmap$Chr[pj]; mk_cM <- pmap$cM[pj]
ordm <- order(mk_chr, mk_cM)                     # meiosis walks markers in cM order
f_aq_par <- colMeans(ph$aqu[, pj, drop=FALSE]); f_pol_par <- colMeans(ph$pol[, pj, drop=FALSE])
p_pooled <- (colSums(ph$aqu[, pj, drop=FALSE]) + colSums(ph$pol[, pj, drop=FALSE])) /
            (nrow(ph$aqu) + nrow(ph$pol))
parent_maf <- pmin(p_pooled, 1-p_pooled); names(parent_maf) <- markers

## per-chromosome marker index + genetic length (for crossover simulation)
chr_idx <- split(seq_along(markers), mk_chr)
chr_len <- sapply(chr_idx, function(i) diff(range(mk_cM[i])))

## ------------------- meiosis: one recombinant haplotype ---------------------
## crossovers ~ Poisson(L_cM/100) per chromosome, positions uniform in cM; the
## parental haplotype alternates at each crossover (start chosen at random).
meiosis <- function(h1, h2) {
  out <- h1
  for (ch in names(chr_idx)) {
    i  <- chr_idx[[ch]]; cm <- mk_cM[i]
    L  <- chr_len[[ch]]; if (!is.finite(L) || L <= 0) L <- 0
    nx <- rpois(1, L/100)
    xs <- if (nx > 0) sort(runif(nx, min(cm), max(cm))) else numeric(0)
    par <- (findInterval(cm, xs) + sample(0:1, 1)) %% 2      # 0 -> h1, 1 -> h2
    seg <- ifelse(par == 0, h1[i], h2[i])
    out[i] <- seg
  }
  out
}

## ------------------- read a deme VCF -> queens' phased haplotypes ------------
read_queens <- function(vcf) {
  dt <- fread(vcf, skip="#CHROM", header=TRUE, sep="\t", showProgress=FALSE)
  mk <- paste(sub("ch","Chr",dt[[1]]), dt[[2]], sep=":")
  G  <- as.matrix(dt[, 10:ncol(dt)]); gt <- sub(":.*$","",G)
  a1 <- suppressWarnings(as.integer(sub("[|/].*$","",gt)))
  a2 <- suppressWarnings(as.integer(sub("^.*[|/]","",gt)))
  H1 <- matrix(a1, nrow=nrow(G)); H2 <- matrix(a2, nrow=nrow(G))
  m  <- match(markers, mk)                     # markers lost to drift -> 0
  f  <- function(H) { M <- matrix(0L, length(markers), ncol(H))
                      keep <- !is.na(m); M[keep, ] <- H[m[keep], , drop=FALSE]; M }
  list(h1 = f(H1), h2 = f(H2))                 # markers x queens
}

## ------------------- colony-structured sample for one deme ------------------
## nest_sizes: workers per colony, e.g. c(4,3,3). Returns markers x workers dosage.
colony_sample <- function(q, nest_sizes) {
  nq <- ncol(q$h1); if (nq < 2) return(NULL)
  qsel <- sample(nq, min(length(nest_sizes), nq))
  cols <- vector("list", 0)
  for (j in seq_along(qsel)) {
    k <- nest_sizes[j]; qi <- qsel[j]
    ## sires: haploid males = recombinants of OTHER queens in the deme
    nm  <- sample(seq_along(MATE_W), 1, prob=MATE_W)
    src <- sample(setdiff(seq_len(nq), qi), min(nm, nq-1), replace=FALSE)
    sires <- lapply(src, function(s) meiosis(q$h1[, s], q$h2[, s]))
    for (w in seq_len(k)) {
      mat <- meiosis(q$h1[, qi], q$h2[, qi])          # maternal recombinant
      pat <- sires[[sample(length(sires), 1)]]        # one of her mates
      cols[[length(cols)+1L]] <- mat + pat            # dosage 0/1/2
    }
  }
  do.call(cbind, cols)
}

## ------------------- empirical reference ------------------------------------
message("[2] empirical ...")
nest <- fread(NESTCSV)
e <- new.env(); load(EMP_RD, envir=e)
GTh <- e$GTs_hybrids_005; sdh <- as.data.table(e$sample_data); maph <- as.data.table(e$map_hyb_005)
ci  <- match(markers, maph$marker)
hyb_pops <- nest$Population
emp_means <- do.call(rbind, lapply(hyb_pops, function(pp)
  colMeans(GTh[sdh$Population==pp, ci, drop=FALSE], na.rm=TRUE)))
rownames(emp_means) <- hyb_pops; colnames(emp_means) <- markers

run_sort <- function(pop_means) {
  pm <- rbind(pop_means, aquilonia_parent=f_aq_par*2, polyctena_parent=f_pol_par*2)
  colnames(pm) <- markers
  parallelism_stats(list(pop_means=pm), hybrid_pops=rownames(pop_means),
    aqu_pops="aquilonia_parent", pol_pops="polyctena_parent",
    fix_th=FIX_TH, sort_th=SORT_TH, null_prob=NULL_PROB,
    parent_maf=parent_maf, min_parent_maf=MIN_PARENT_MAF)
}
pi_of  <- function(pm){ p<-pm/2; mean(rowMeans(2*p*(1-p), na.rm=TRUE), na.rm=TRUE) }
fst_of <- function(pm){ f<-pm/2; pb<-colMeans(f,na.rm=TRUE); Hs<-colMeans(2*f*(1-f),na.rm=TRUE)
                        Ht<-2*pb*(1-pb); sum(Ht-Hs,na.rm=TRUE)/sum(Ht,na.rm=TRUE) }
emp_sort <- run_sort(emp_means)
emp_stats <- list(pi=pi_of(emp_means), fst=fst_of(emp_means),
                  sc=emp_sort[differentiated==TRUE, .N, by=sort_class])

## ------------------- simulate colony sampling per cell ----------------------
res <- list()
for (cell in CELLS) {
  message("[3] cell ", cell$label)
  vcfs <- list.files(cell$dir,
    sprintf("hyb_neutral_realfounders_%s_rep[0-9]+_gen%d[.]vcf$", cell$tag, cell$gen),
    full.names=TRUE)
  if (!length(vcfs)) { message("   (no runs)"); next }
  nsz <- lapply(nest$nest_sizes, function(s) as.integer(strsplit(s, ",")[[1]]))

  ## RANDOM sampling (what we did before) vs COLONY sampling, same demes
  pm_rand <- matrix(NA_real_, length(vcfs), length(markers))
  pm_col  <- matrix(NA_real_, length(vcfs), length(markers))
  for (i in seq_along(vcfs)) {
    q  <- read_queens(vcfs[i])
    ns <- nsz[[(i-1) %% length(nsz) + 1]]
    ## random: draw sum(ns) unrelated queens (matched n, no family structure)
    nq <- ncol(q$h1); take <- min(sum(ns), nq)
    rs <- sample(nq, take)
    pm_rand[i, ] <- rowMeans(q$h1[, rs, drop=FALSE] + q$h2[, rs, drop=FALSE])
    ## colony: workers from length(ns) colonies
    W <- colony_sample(q, ns)
    if (!is.null(W)) pm_col[i, ] <- rowMeans(W)
  }
  dimnames(pm_rand) <- dimnames(pm_col) <- list(paste0("deme",seq_along(vcfs)), markers)
  s_rand <- run_sort(pm_rand); s_col <- run_sort(pm_col)
  sc <- function(s) { t <- s[differentiated==TRUE, .N, by=sort_class]
                      setNames(t$N/sum(t$N), t$sort_class) }
  res[[length(res)+1L]] <- data.table(
    cell = cell$label,
    sampling = c("random (unrelated)", "colony (nest-structured)"),
    pi  = c(pi_of(pm_rand),  pi_of(pm_col)),
    fst = c(fst_of(pm_rand), fst_of(pm_col)),
    sorted = c(1 - sc(s_rand)[["unsorted"]], 1 - sc(s_col)[["unsorted"]]))
}
out <- rbindlist(res)
emp_sorted <- 1 - emp_stats$sc[sort_class=="unsorted", N]/sum(emp_stats$sc$N)
saveRDS(list(sim=out, emp_pi=emp_stats$pi, emp_fst=emp_stats$fst,
             emp_sorted=emp_sorted, nest=nest), OUTRDS)

cat(sprintf("\n=== EMPIRICAL: pi=%.3f  F_ST=%.3f  sorted=%.3f ===\n",
            emp_stats$pi, emp_stats$fst, emp_sorted))
cat("\n=== simulated, random vs COLONY sampling (same demes) ===\n")
print(out[, .(cell, sampling, pi=round(pi,3), fst=round(fst,3), sorted=round(sorted,3))])
cat("\nREAD: if colony sampling moves sim pi DOWN, F_ST UP and sorted UP toward the\n")
cat("empirical values, the 'excess' was substantially a family-sampling artefact.\n")
cat("\nsaved: ", OUTRDS, "\n")
