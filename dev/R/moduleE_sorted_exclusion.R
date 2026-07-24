## =============================================================================
## Module E, E2h -- exclude SORTED loci, then re-compare empirical vs neutral null
## =============================================================================
##
## THE QUESTION
## ------------
## The K-sweep and the founder-heterogeneity ceiling both failed to explain the
## empirical excess (pi 0.299 vs sim ~0.42; among-pop F_ST 0.313 vs sim ~0.12-0.20).
## But those statistics were measured at DI>-25 markers, which INCLUDE the sorted
## loci -- so they cannot serve as both the demographic calibration and the test.
##
## Fix: remove the sorted loci and ask whether the remaining BACKGROUND is
## consistent with the neutral null.
##   - background MATCHES  -> the neutral demography is validated on sorting-free
##     loci and the excess is CONFINED to sorted loci (supports localised selection)
##   - background STILL MISMATCHES -> the discrepancy is genome-wide, i.e.
##     demographic misspecification (weakens the selection interpretation)
##
## TWO THINGS THAT MAKE OR BREAK THIS TEST
## ---------------------------------------
## 1. SYMMETRIC EXCLUSION. "Sorted" == near-fixation in different directions across
##    populations == high among-population differentiation BY DEFINITION. Dropping
##    sorted loci from the empirical alone mechanically lowers its F_ST and raises
##    its pi, moving it toward the sim whether or not selection is involved. So
##    sort_class is computed on the SIMULATED demes too, with the same code path and
##    thresholds, and loci sorted in EITHER dataset are dropped -> a COMMON marker
##    set for the comparison.
## 2. SAMPLE-SIZE MATCHING. sort_class depends on n: a 3-diploid population can
##    near-fix trivially, 50 diploids essentially cannot. Sim demes are therefore
##    subsampled to the empirical per-population sizes before anything is computed.
##    Without this the sim looks artificially unsorted.
##
## Both datasets are oriented with the SAME parental reference (the real phased
## founder pool), and use the locked Module A conventions
## (parent_maf >= 0.15, fix_th = 0.15, sort_th = 0.5, null_prob = 0.5).
##
## BY-PRODUCT: this yields sort_class proportions in the neutral null vs observed
## -- i.e. the step (b) sorting statistic -- reported below.
##
## Output: data/moduleE_sim/moduleE_sorted_exclusion.rds
##         Figures/moduleE_sorted_exclusion.pdf
## =============================================================================

suppressMessages({library(data.table); library(ggplot2); library(parallel)})
source("dev/R/parallelism_stats.R")          # the locked Module A statistic

## ------------------------------- CONFIG -------------------------------------
FORMICA   <- "/Users/petrikem/gitlab/formica_hybrid"
SCREEN    <- file.path(FORMICA, "data/moduleE_sim/screen")
KSWEEP    <- file.path(FORMICA, "data/moduleE_sim/ksweep")
POOL      <- file.path(FORMICA, "data/moduleE_founders/moduleE_founder_haplotypes.rds")
FVCF_DIR  <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/SLiM/founders_DIm25_thin15000"
EMP_RDATA <- file.path(FORMICA, "data/hybrids_only_maf005.Rdata")
EMLG      <- file.path(FORMICA, "data/eMLG_5loci_0025_cM05.rds")
OUTRDS    <- file.path(FORMICA, "data/moduleE_sim/moduleE_sorted_exclusion.rds")
FIGDIR    <- file.path(FORMICA, "Figures")

## LD-decay-matched neutral cells to test (label, dir, tag, gen)
CELLS <- list(
  list(label = "30/13 K6250 gen40",  dir = SCREEN, tag = "Naq30_Npol13",      gen = 40),
  list(label = "12/12 K6250 gen60",  dir = SCREEN, tag = "Naq12_Npol12",      gen = 60),
  list(label = "12/12 K375 gen100",  dir = KSWEEP, tag = "Naq12_Npol12_K375", gen = 100)
)

## locked Module A conventions
FIX_TH <- 0.15; SORT_TH <- 0.5; MIN_PARENT_MAF <- 0.15; NULL_PROB <- 0.5
MATCH_SEED <- 1L
dir.create(FIGDIR, showWarnings = FALSE)

## ------------------- 1. markers, parental reference -------------------------
message("[1] markers + parental reference ...")
sim_markers <- unique(unlist(lapply(
  list.files(FVCF_DIR, "founders_ch.*[.]vcf$", full.names = TRUE),
  function(f) fread(f, skip = "#CHROM", select = 3, header = TRUE)[[1]])))

ph  <- readRDS(POOL); pmap <- as.data.table(ph$map)
pj  <- match(sim_markers, pmap$marker); keep <- !is.na(pj); pj <- pj[keep]
markers <- sim_markers[keep]
## parental allele frequencies from the REAL phased founder pool (same reference
## used to orient BOTH datasets). aq = 30 haplotypes; pol = 13 diploid females.
f_aq_par  <- colMeans(ph$aqu[, pj, drop = FALSE])
Hpol      <- ph$pol[, pj, drop = FALSE]
f_pol_par <- colMeans(Hpol)
## pooled parental MAF (the primary gate)
p_pooled   <- (colSums(ph$aqu[, pj, drop=FALSE]) + colSums(Hpol)) /
              (nrow(ph$aqu) + nrow(Hpol))
parent_maf <- pmin(p_pooled, 1 - p_pooled); names(parent_maf) <- markers
message(sprintf("    markers: %d ; parental |dp| mean = %.3f",
                length(markers), mean(abs(f_aq_par - f_pol_par))))

## ------------------- 2. empirical: per-population dosage --------------------
message("[2] empirical hybrids ...")
e <- new.env(); load(EMP_RDATA, envir = e)
GTh <- e$GTs_hybrids_005; sdh <- as.data.table(e$sample_data); maph <- as.data.table(e$map_hyb_005)
ci  <- match(markers, maph$marker)
hyb_pops <- setdiff(unique(sdh$Population), c("aquilonia_parent","polyctena_parent"))
emp_ns   <- sapply(hyb_pops, function(pp) sum(sdh$Population == pp))
emp_means <- do.call(rbind, lapply(hyb_pops, function(pp)
  colMeans(GTh[sdh$Population == pp, ci, drop = FALSE], na.rm = TRUE)))   # dosage 0-2
rownames(emp_means) <- hyb_pops; colnames(emp_means) <- markers
message(sprintf("    %d pops, sizes %s", length(hyb_pops), paste(sort(emp_ns), collapse=",")))

## ------------------- 3. sort_class via the locked statistic -----------------
## build a prep-like object: rows = populations + the two parental references
run_sort <- function(pop_means) {
  pm <- rbind(pop_means,
              aquilonia_parent = f_aq_par * 2,   # freq -> dosage scale
              polyctena_parent = f_pol_par * 2)
  colnames(pm) <- markers
  parallelism_stats(list(pop_means = pm),
                    hybrid_pops = rownames(pop_means),
                    aqu_pops = "aquilonia_parent", pol_pops = "polyctena_parent",
                    fix_th = FIX_TH, sort_th = SORT_TH, null_prob = NULL_PROB,
                    parent_maf = parent_maf, min_parent_maf = MIN_PARENT_MAF)
}
emp_sort <- run_sort(emp_means)
message("    empirical sort_class: ",
        paste(sprintf("%s=%d", emp_sort[differentiated==TRUE, .N, by=sort_class]$sort_class,
                      emp_sort[differentiated==TRUE, .N, by=sort_class]$N), collapse=" "))

## ------------------- 4. simulated cells -------------------------------------
read_vcf_dosage <- function(vcf) {
  dt <- fread(vcf, skip="#CHROM", header=TRUE, sep="\t", showProgress=FALSE)
  mk <- paste(sub("ch","Chr",dt[[1]]), dt[[2]], sep=":")
  G <- as.matrix(dt[, 10:ncol(dt)]); gt <- sub(":.*$","",G)
  d <- matrix(suppressWarnings(as.integer(sub("[|/].*$","",gt))) +
              suppressWarnings(as.integer(sub("^.*[|/]","",gt))), nrow=nrow(G))
  list(dose = d, marker = mk)
}
sim_pop_means <- function(cell) {
  pat  <- sprintf("hyb_%s%s_rep[0-9]+_gen%d[.]vcf$",
                  if (grepl("_K", cell$tag)) "neutral_realfounders_" else "neutral_realfounders_",
                  cell$tag, cell$gen)
  vcfs <- list.files(cell$dir, pat, full.names = TRUE)
  if (!length(vcfs)) stop("no VCFs for cell ", cell$label)
  M <- matrix(NA_real_, nrow = length(vcfs), ncol = length(markers),
              dimnames = list(paste0("deme", seq_along(vcfs)), markers))
  for (i in seq_along(vcfs)) {
    x <- read_vcf_dosage(vcfs[i])
    ## SAMPLE-SIZE MATCH to an empirical population before anything is computed
    tgt <- emp_ns[(i-1L) %% length(emp_ns) + 1L]; ni <- ncol(x$dose)
    if (ni > tgt) { set.seed(MATCH_SEED + i); x$dose <- x$dose[, sample(ni, tgt), drop=FALSE] }
    mm <- match(markers, x$marker)
    M[i, ] <- ifelse(is.na(mm), NA_real_, rowMeans(x$dose, na.rm=TRUE)[mm])
  }
  ## markers absent from a deme's VCF were lost to drift -> fixed for the
  ## reference allele in that deme (dosage 0), not missing
  M[is.na(M)] <- 0
  M
}

## ------------------- 5. eMLG cluster map (for cluster-level exclusion) ------
message("[3] eMLG cluster membership ...")
eml <- readRDS(EMLG)
grp <- as.data.table(eml$groups)
mk2cl <- rbindlist(lapply(seq_len(nrow(grp)), function(i)
  data.table(marker = grp$members[[i]], cluster = grp$group_id[i])))
setkey(mk2cl, marker)
cl_of <- mk2cl[.(markers), cluster]           # NA if not in any cluster
message(sprintf("    markers with a cluster: %d / %d", sum(!is.na(cl_of)), length(markers)))

## ------------------- 6. statistics on marker subsets ------------------------
pi_of  <- function(pop_means, idx) {
  p <- pop_means[, idx, drop=FALSE] / 2
  mean(rowMeans(2 * p * (1 - p), na.rm = TRUE), na.rm = TRUE)   # mean within-pop He
}
fst_of <- function(pop_means, idx) {
  f <- pop_means[, idx, drop=FALSE] / 2
  pbar <- colMeans(f, na.rm=TRUE); Hs <- colMeans(2*f*(1-f), na.rm=TRUE)
  Ht <- 2*pbar*(1-pbar); sum(Ht-Hs, na.rm=TRUE)/sum(Ht, na.rm=TRUE)
}
fst_per_locus <- function(pop_means, idx) {
  f <- pop_means[, idx, drop=FALSE] / 2
  pbar <- colMeans(f, na.rm=TRUE); Hs <- colMeans(2*f*(1-f), na.rm=TRUE)
  Ht <- 2*pbar*(1-pbar); ifelse(Ht > 0, (Ht-Hs)/Ht, NA_real_)
}

results <- list(); dists <- list()
for (cell in CELLS) {
  message("[4] cell: ", cell$label)
  sm <- sim_pop_means(cell)
  sim_sort <- run_sort(sm)

  ## align sort_class to the marker vector
  es <- emp_sort[match(markers, marker)]; ss <- sim_sort[match(markers, marker)]
  gated <- es$differentiated & ss$differentiated & !is.na(es$sort_class) & !is.na(ss$sort_class)
  emp_sorted <- gated & es$sort_class != "unsorted"
  sim_sorted <- gated & ss$sort_class != "unsorted"

  ## --- exclusion sets (COMMON marker set: drop if sorted in EITHER dataset) ---
  idx_all   <- which(gated)
  idx_locus <- which(gated & !emp_sorted & !sim_sorted)
  bad_cl    <- unique(cl_of[which(emp_sorted | sim_sorted)]); bad_cl <- bad_cl[!is.na(bad_cl)]
  idx_clust <- which(gated & !emp_sorted & !sim_sorted & !(cl_of %in% bad_cl))

  tab <- rbindlist(lapply(
    list(list("all gated", idx_all), list("exclude sorted loci", idx_locus),
         list("exclude sorted clusters", idx_clust)),
    function(z) data.table(
      cell = cell$label, subset = z[[1]], n_markers = length(z[[2]]),
      emp_pi = pi_of(emp_means, z[[2]]), sim_pi = pi_of(sm, z[[2]]),
      emp_fst = fst_of(emp_means, z[[2]]), sim_fst = fst_of(sm, z[[2]]))))
  tab[, `:=`(pi_gap = emp_pi - sim_pi, fst_gap = emp_fst - sim_fst)]
  results[[length(results)+1L]] <- tab

  ## threshold-free: per-locus F_ST distributions on all gated markers
  dists[[length(dists)+1L]] <- data.table(
    cell = cell$label,
    fst = c(fst_per_locus(emp_means, idx_all), fst_per_locus(sm, idx_all)),
    src = rep(c("empirical","simulated"), each = length(idx_all)))

  ## (b) by-product: sort_class proportions, observed vs neutral null
  sc <- merge(es[gated, .N, by=sort_class][, .(sort_class, empirical = N/sum(N))],
              ss[gated, .N, by=sort_class][, .(sort_class, simulated = N/sum(N))],
              by = "sort_class", all = TRUE)
  results[[length(results)]] <- tab
  attr(results[[length(results)]], "sort_class") <- sc
  message("    sort_class (obs vs null): ",
          paste(sprintf("%s %.3f/%.3f", sc$sort_class, sc$empirical, sc$simulated), collapse="  "))
}

res  <- rbindlist(results)
dist <- rbindlist(dists)
saveRDS(list(results = res, fst_dist = dist,
             sort_class = lapply(results, attr, "sort_class"), cells = CELLS,
             conventions = list(FIX_TH=FIX_TH, SORT_TH=SORT_TH,
                                MIN_PARENT_MAF=MIN_PARENT_MAF)), OUTRDS)

## ------------------- 7. report + figure -------------------------------------
cat("\n=== pi / F_ST before and after SYMMETRIC exclusion of sorted loci ===\n")
print(res[, .(cell, subset, n_markers,
              emp_pi = round(emp_pi,3), sim_pi = round(sim_pi,3), pi_gap = round(pi_gap,3),
              emp_fst = round(emp_fst,3), sim_fst = round(sim_fst,3), fst_gap = round(fst_gap,3))])
cat("\nREAD: if pi_gap and fst_gap collapse toward 0 in the excluded subsets, the\n")
cat("neutral null fits the sorting-free background and the excess is confined to\n")
cat("the sorted loci. If the gaps persist, the mismatch is genome-wide.\n")

p <- ggplot(dist, aes(fst, colour = src)) + stat_ecdf(linewidth = 0.8) +
  facet_wrap(~ cell) + coord_cartesian(xlim = c(0, 1)) +
  labs(x = "per-locus F_ST", y = "ECDF", colour = NULL,
       title = "Per-locus F_ST: observed vs neutral null (all gated markers)",
       subtitle = "threshold-free view: background should overlay; excess appears as an upper tail") +
  theme_bw(base_size = 11)
ggsave(file.path(FIGDIR, "moduleE_sorted_exclusion.pdf"), p, width = 11, height = 4.5)
cat("\nsaved: ", OUTRDS, "\n       ", file.path(FIGDIR,"moduleE_sorted_exclusion.pdf"), "\n")
