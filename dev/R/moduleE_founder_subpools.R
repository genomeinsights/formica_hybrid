## =============================================================================
## Module E, E2g-1 -- founder SUB-POOLS for the founder-heterogeneity experiment
## =============================================================================
##
## QUESTION THIS SERVES
## --------------------
## The K-sweep showed no neutral (K, generation) reproduces the empirical
## among-population F_ST (0.313) while holding the LD-decay match -- the best
## acceptable-LD cell reached ~0.20. But our demes all draw founders from the SAME
## rangewide pool, so founders are near-identical across demes and among-deme F_ST
## is low BY CONSTRUCTION. If the 20 real hybrid populations were founded from
## DIFFERENT local source populations, their among-pop F_ST would be inflated
## NEUTRALLY by founder heterogeneity. This script sets up that test.
##
## CHEAP UPPER BOUND (read this before running any simulation)
## ----------------------------------------------------------
## Founder heterogeneity can only inflate among-deme F_ST up to the differentiation
## that actually exists AMONG the parental sub-pools, measured at the SAME markers.
## So we report, per k:
##   - within-species F_ST among the k aq sub-pools and among the k pol sub-pools
##   - the AT-FOUNDING among-deme F_ST: F_ST among the k admixed founding frequency
##     vectors (0.5*aq_subpool + 0.5*pol_subpool), i.e. the head start heterogeneity
##     gives before any drift.
## If that at-founding F_ST is far below the 0.11 gap we need to close, founder
## heterogeneity CANNOT rescue the neutral explanation -- and no simulation is needed.
##
## CONTROLLING THE CONFOUND
## ------------------------
## A local draw necessarily uses fewer founders, and fewer founders alone raises
## F_ST. So N_AQ/N_POL are held FIXED across all heterogeneity levels (k=1 rangewide
## control vs k>1 local), and must fit inside the smallest sub-pool. The script
## validates this and reports the max feasible N per k.
##
## OUTPUT
##   data/moduleE_sim/founder_subpools.rds        -- assignments + diagnostics
##   data/moduleE_sim/founder_assignments.csv     -- harness input:
##        k, deme, subpool, aq_keep, pol_keep   (1-based indices, comma-separated)
##
## Review, then run (fast; pure data prep + diagnostics). NOT yet the simulation.
## =============================================================================

suppressMessages({library(data.table)})

## ------------------------------- CONFIG -------------------------------------
FORMICA <- "/Users/petrikem/gitlab/formica_hybrid"
POOL    <- file.path(FORMICA, "data/moduleE_founders/moduleE_founder_haplotypes.rds")
FVCF    <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/SLiM/founders_DIm25_thin15000"
OUTRDS  <- file.path(FORMICA, "data/moduleE_sim/founder_subpools.rds")
OUTCSV  <- file.path(FORMICA, "data/moduleE_sim/founder_assignments.csv")

KS      <- c(2, 3, 4)      # heterogeneity levels (number of sub-pools); k=1 = control
NDEMES  <- 20
## founder numbers, held FIXED across k (must fit the smallest sub-pool at max k)
N_AQ    <- 7L
N_POL   <- 3L
SEED    <- 1L
## ---------------------------------------------------------------------------

set.seed(SEED)
dir.create(dirname(OUTRDS), recursive = TRUE, showWarnings = FALSE)

## ---- founder pool, restricted to the markers the sim actually carries ----
ph  <- readRDS(POOL)
map <- as.data.table(ph$map)
sim_markers <- unique(unlist(lapply(
  list.files(FVCF, "founders_ch.*\\.vcf$", full.names = TRUE),
  function(f) fread(f, skip = "#CHROM", select = 3, header = TRUE)[[1]])))
ci <- match(sim_markers, map$marker); ci <- ci[!is.na(ci)]
cat(sprintf("markers used (= the sim's thinned DI>-25 set): %d\n", length(ci)))

Haq <- ph$aqu[, ci, drop = FALSE]                       # 30 haplotypes x markers (0/1)
Hpol<- ph$pol[, ci, drop = FALSE]                       # 26 haplotypes
## pol as 13 diploid females (rows 2i-1,2i are one female)
Dpol <- Hpol[seq(1, nrow(Hpol), 2), , drop = FALSE] + Hpol[seq(2, nrow(Hpol), 2), , drop = FALSE]
cat(sprintf("pool: %d aq haplotypes, %d pol diploid females\n", nrow(Haq), nrow(Dpol)))

## ---- helper: F_ST among groups given a group x marker frequency matrix ----
fst_hat <- function(freqs) {
  pbar <- colMeans(freqs, na.rm = TRUE)
  Hs   <- colMeans(2 * freqs * (1 - freqs), na.rm = TRUE)
  Ht   <- 2 * pbar * (1 - pbar)
  sum(Ht - Hs, na.rm = TRUE) / sum(Ht, na.rm = TRUE)
}
grp_freq <- function(M, grp, ploidy) {   # M: units x markers; returns groups x markers freq
  do.call(rbind, lapply(sort(unique(grp)),
    function(g) colMeans(M[grp == g, , drop = FALSE], na.rm = TRUE) / ploidy))
}

## ---- cluster each species into k sub-pools (structure that really exists) ----
## distance = Manhattan on allele state, so clusters follow real within-species
## structure at these markers (no locality labels are available for the parents).
d_aq  <- dist(Haq,  method = "manhattan")
d_pol <- dist(Dpol, method = "manhattan")
h_aq  <- hclust(d_aq,  method = "ward.D2")
h_pol <- hclust(d_pol, method = "ward.D2")

diag_rows <- list(); assign_rows <- list()

for (k in KS) {
  g_aq  <- cutree(h_aq,  k)
  g_pol <- cutree(h_pol, k)
  sz_aq <- table(g_aq); sz_pol <- table(g_pol)

  ## --- CEILING DIAGNOSTICS at the sim's markers ---
  f_aq  <- grp_freq(Haq,  g_aq,  1)      # haplotypes -> freq
  f_pol <- grp_freq(Dpol, g_pol, 2)      # diploid dosage -> freq
  fst_aq  <- fst_hat(f_aq)
  fst_pol <- fst_hat(f_pol)
  ## at-founding among-deme F_ST if demes are founded from DIFFERENT sub-pools
  ## (50/50 admixture of sub-pool s of each species)
  n_s <- min(k, nrow(f_aq), nrow(f_pol))
  f_found <- 0.5 * f_aq[seq_len(n_s), , drop = FALSE] + 0.5 * f_pol[seq_len(n_s), , drop = FALSE]
  fst_founding <- fst_hat(f_found)

  diag_rows[[length(diag_rows) + 1L]] <- data.table(
    k = k, min_aq = min(sz_aq), min_pol = min(sz_pol),
    max_feasible_N_AQ = min(sz_aq), max_feasible_N_POL = min(sz_pol),
    fst_among_aq_subpools = fst_aq, fst_among_pol_subpools = fst_pol,
    fst_at_founding = fst_founding)

  ## --- per-deme founder assignment: deme d uses sub-pool (d mod k) ---
  feasible <- (N_AQ <= min(sz_aq)) && (N_POL <= min(sz_pol))
  for (d in seq_len(NDEMES)) {
    s <- ((d - 1L) %% k) + 1L
    aq_pool  <- which(g_aq  == s); pol_pool <- which(g_pol == s)
    if (!feasible) next
    assign_rows[[length(assign_rows) + 1L]] <- data.table(
      k = k, deme = d, subpool = s,
      aq_keep  = paste(sort(sample(aq_pool,  N_AQ)),  collapse = ","),
      pol_keep = paste(sort(sample(pol_pool, N_POL)), collapse = ","))
  }
  if (!feasible)
    cat(sprintf("  !! k=%d SKIPPED: N_AQ=%d/N_POL=%d exceed smallest sub-pool (%d aq, %d pol)\n",
                k, N_AQ, N_POL, min(sz_aq), min(sz_pol)))
}

## ---- k = 1 CONTROL: same N, drawn RANGEWIDE from the whole pool ----
for (d in seq_len(NDEMES)) {
  assign_rows[[length(assign_rows) + 1L]] <- data.table(
    k = 1L, deme = d, subpool = 1L,
    aq_keep  = paste(sort(sample(nrow(Haq),  N_AQ)),  collapse = ","),
    pol_keep = paste(sort(sample(nrow(Dpol), N_POL)), collapse = ","))
}

diags   <- rbindlist(diag_rows)
assigns <- rbindlist(assign_rows)[order(k, deme)]

fwrite(assigns, OUTCSV)
saveRDS(list(diagnostics = diags, assignments = assigns,
             N_AQ = N_AQ, N_POL = N_POL, markers = length(ci),
             clust = list(aq = h_aq, pol = h_pol)), OUTRDS)

## ------------------------------ REPORT --------------------------------------
cat("\n=== sub-pool sizes and CEILING diagnostics (at the sim's markers) ===\n")
print(diags[, .(k, min_aq, min_pol,
                fst_aq = round(fst_among_aq_subpools, 4),
                fst_pol = round(fst_among_pol_subpools, 4),
                fst_at_founding = round(fst_at_founding, 4))])
cat(sprintf("\nfounder numbers held fixed across k: N_AQ=%d, N_POL=%d\n", N_AQ, N_POL))
cat(sprintf("assignments written for k = %s (+ k=1 rangewide control), %d demes each\n",
            paste(sort(unique(assigns$k[assigns$k > 1])), collapse = ","), NDEMES))
cat("\nINTERPRETATION: 'fst_at_founding' is the among-deme F_ST that founder\n")
cat("heterogeneity alone supplies at generation 0 (drift then adds on top).\n")
cat("The K-sweep left a gap of ~0.11 to empirical F_ST (0.313 vs ~0.20 reached).\n")
cat("If fst_at_founding is well below that gap, founder heterogeneity CANNOT\n")
cat("rescue a neutral explanation and the simulations are not needed.\n")
cat("\nsaved: ", OUTCSV, "\n       ", OUTRDS, "\n")
