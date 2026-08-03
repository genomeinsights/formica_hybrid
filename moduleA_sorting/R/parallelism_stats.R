## =========================================================
## Binomial parallelism test for ancestry sorting
## =========================================================
##
## Companion to sorting_stats(): where sorting_stats() measures how MUCH
## each locus has sorted (magnitude), parallelism_stats() asks whether the
## DIRECTION of that sorting is repeatable across the replicate hybrid
## populations.
##
## Biology / the two scenarios this separates
## -------------------------------------------
## Under a Dobzhansky-Muller incompatibility (or any locus under
## consistent selection) we expect near-fixation of a parental allele in
## each population. Two sub-cases:
##
##   (1) One allelic combination outperforms the other
##         -> the SAME parental allele fixes in (nearly) all populations
##         -> direction is predictable / parallel
##         -> binomial test rejects the random-direction null.
##
##   (2) The two combinations are equally fit
##         -> populations still near-fix, but in RANDOM directions
##         -> direction is not predictable
##         -> binomial test does NOT reject the null.
##
## Note that plain neutral drift in small, young populations also produces
## case (2) (near-fixation in random directions). So a SIGNIFICANT
## departure toward one parent -- case (1) -- is the signal of predictable,
## selection-like sorting; it is the direct, per-locus analogue of the
## cross-population Spearman correlations used in Nouhaud et al. (2022,
## PLOS Biology), but far more powerful here because we have ~20 replicate
## populations rather than 3.
##
##
## Input is a `prep` object from ohta_fast_prepare() (LDscnR package), the
## same object sorting_stats() consumes, so the two results join on
## `marker`.

library(data.table)

## ---------------------------------------------------------
## parallelism_stats()
##
## prep         : list from ohta_fast_prepare() (needs $pop_means, rows =
##                populations, cols = markers).
## hybrid_pops  : character vector of hybrid population names to test.
## aqu_pops     : parental population(s) used to define the aquilonia allele.
## pol_pops     : parental population(s) used to define the polyctena allele.
## fix_th       : near-fixation tolerance. A population counts toward
##                aquilonia if its (oriented) aquilonia-allele freq
##                >= 1 - fix_th, toward polyctena if <= fix_th. Note this
##                interacts with population sample size (a 3-diploid pop can
##                only reach freq 0 or 1); loosen for small pops.
## DI           : DiagnosticIndex per marker, as a marker-NAMED numeric
##                vector (aligned by name), or an unnamed vector in
##                colnames(prep$pop_means) order. Reported in the output and
##                used as a COVARIATE downstream (DI vs recombination / cluster
##                size / PC association). It is NOT the sorting gate: gating on
##                DI would truncate its range and make "does DI predict
##                sorting?" circular. Pass e.g. setNames(map$DiagnosticIndex, map$marker).
## parent_maf   : pooled-parental minor-allele frequency per marker, marker-
##                NAMED (or unnamed in prep order). From POOLED parental allele
##                counts, e.g. p <- colMeans(GTs_parents, na.rm=TRUE)/2;
##                pmin(p, 1-p). The PRIMARY sorting gate (see min_parent_maf):
##                keeps only loci polymorphic in the parents, so a high sorting
##                index reflects real sorting rather than an allele that was
##                near-monomorphic in the founding pool -- while leaving DI free
##                to vary. Also reported in the output.
## min_parent_maf : PRIMARY gate. Loci with parent_maf < min_parent_maf are not
##                tested / not classified (sort_class = NA). NULL = off.
## min_DI       : optional SECONDARY ancestry-diagnostic gate (loci with
##                DI < min_DI dropped). Default NULL (off) -- available but not
##                the primary gate; prefer min_parent_maf. Requires DI.
## min_parent_diff : optional secondary safeguard, gating on empirical parental
##                differentiation |f_aqu_parent - f_pol_parent| >= min_parent_diff
##                (only when orient = "parents"). Default 0 (off).
## admix_prop   : assumed aquilonia admixture proportion at founding; used
##                only to report founding_f (the expected founding aquilonia
##                frequency, admix_prop*f_aqu_parent + (1-admix_prop)*f_pol_parent).
## sort_th      : THE UNIFIED THRESHOLD. A fraction of the hybrid-population
##                panel. prop_fixed decomposes exactly as
##                    prop_fixed = |uni_score| + bi_score
##                with uni_score = (n_aqu-n_pol)/n_obs (signed, +=aquilonia)
##                and bi_score = 2*min(n_aqu,n_pol)/n_obs, both fractions of
##                populations in [0,1]. A unit is classed by whichever
##                pattern wins, if the winner reaches sort_th:
##                  |uni_score| >= bi_score & |uni_score| >= sort_th
##                     -> "aquilonia" / "polyctena" (unidirectional)
##                  bi_score > |uni_score| & bi_score >= sort_th
##                     -> "bidirectional"
##                  else "unsorted".
## sort_rule    : how sort_class is called from the fixation counts.
##                "prop_fixed" (DEFAULT): magnitude gate = total near-fixation
##                (prop_fixed >= sort_th), then the uni/bi split (bidirectional iff
##                minority > 1/4 of fixed pops). "component" (ORIGINAL): the larger
##                of |uni_score|/bi_score must itself reach sort_th (under-calls
##                bidirectional). "binom": magnitude gate = prop_fixed >= sort_th,
##                but DIRECTION is the binomial random-direction test -- unidirectional
##                only when the split is SIGNIFICANTLY biased toward one parent
##                (p_binom < alpha); non-significant sorted loci are "bidirectional"
##                if n_fixed had power (n_fixed >= smallest n with 2*0.5^n < alpha),
##                else "ambiguous". See the inline note at the sort_class block.
## alpha        : significance level for the "binom" direction test (default 0.05).
##                Because all three scores are population fractions, the SAME
##                sort_th is comparable across SNPs and eMLGs (run this
##                function on eMLG consensus genotypes) and genome
##                proportions are weighted tallies of sort_class. (fix_th, by
##                contrast, is the per-population near-fixation convention and
##                should be held constant across all analyses -- it is not the
##                cross-unit comparison knob.)
## orient       : "parents" (default) polarises each locus from aqu_pops /
##                pol_pops allele frequencies; "dosage" trusts the dosage
##                coding directly (see dosage_is_aqu).
## dosage_is_aqu: only used when orient = "dosage"; TRUE if dosage counts
##                the aquilonia-associated allele.
##
## Returns a data.table keyed on `marker`:
##   DI           DiagnosticIndex passed in (NA if none), for convenience
##   n_obs        hybrid pops with data and a defined orientation
##   n_aqu, n_pol pops near-fixed toward aquilonia / polyctena
##   n_fixed      n_aqu + n_pol
##   n_unsorted   n_obs - n_fixed
##   prop_fixed   n_fixed / n_obs  (magnitude of sorting; bidirectional axis)
##   f_aqu_pooled mean hybrid aquilonia-allele frequency
##   f_aqu_parent aquilonia-allele freq in the aquilonia parents (oriented)
##   f_pol_parent aquilonia-allele freq in the polyctena parents (oriented)
##   parent_diff  f_aqu_parent - f_pol_parent (between-species differentiation)
##   parent_maf   pooled-parental minor-allele frequency (the primary gate input)
##   founding_f   expected founding aquilonia-allele frequency
##   differentiated  passed the ancestry-diagnostic gate(s) (min_DI /
##                min_parent_diff) -- i.e. was this locus tested/classified
##   dir_bias     (n_aqu - n_pol) / n_fixed  in [-1, 1]
##   direction    "aquilonia" / "polyctena" / "tie"
##   uni_score    (n_aqu - n_pol) / n_obs  (signed unidirectional fraction)
##   bi_score     2*min(n_aqu,n_pol) / n_obs  (bidirectional fraction)
##   sort_class   "aquilonia"/"polyctena"/"bidirectional"/"unsorted" at sort_th
## ---------------------------------------------------------
parallelism_stats <- function(prep,
                              hybrid_pops,
                              aqu_pops = NULL,
                              pol_pops = NULL,
                              fix_th = 0.1,
                              DI = NULL,
                              min_DI = NULL,
                              parent_maf = NULL,
                              min_parent_maf = NULL,
                              min_parent_diff = 0,
                              admix_prop = 0.5,
                              sort_th = 0.5,
                              sort_rule = c("prop_fixed","component","binom"),
                              alpha = 0.05,
                              orient = c("parents", "dosage"),
                              dosage_is_aqu = TRUE) {

  orient    <- match.arg(orient)
  sort_rule <- match.arg(sort_rule)

  P <- prep$pop_means / 2                 # populations x markers allele freq
  pops_avail <- rownames(P)
  markers <- colnames(P)

  ## DiagnosticIndex aligned to markers (by NAME if named, else assumed in
  ## marker order). Aligning by name avoids position/key-order mismatches.
  DI_vec <- rep(NA_real_, length(markers))
  if (!is.null(DI)) {
    if (!is.null(names(DI))) {
      DI_vec <- as.numeric(DI[markers])
    } else {
      if (length(DI) != length(markers)) {
        stop("unnamed DI must have length ncol(prep$pop_means); prefer a marker-named vector")
      }
      DI_vec <- as.numeric(DI)
    }
  }

  ## pooled-parental MAF aligned to markers (by NAME if named, else marker order)
  parent_maf_vec <- rep(NA_real_, length(markers))
  if (!is.null(parent_maf)) {
    if (!is.null(names(parent_maf))) {
      parent_maf_vec <- as.numeric(parent_maf[markers])
    } else {
      if (length(parent_maf) != length(markers)) {
        stop("unnamed parent_maf must have length ncol(prep$pop_means); prefer a marker-named vector")
      }
      parent_maf_vec <- as.numeric(parent_maf)
    }
  }

  hybrid_pops <- intersect(hybrid_pops, pops_avail)
  if (length(hybrid_pops) == 0L) stop("no hybrid_pops found in prep$pop_means")

  #orient = "parents"
  ## ---- per-locus orientation into aquilonia-allele frequency ----
  if (orient == "parents") {
    aqu_pops <- intersect(aqu_pops, pops_avail)
    pol_pops <- intersect(pol_pops, pops_avail)
    if (length(aqu_pops) == 0L || length(pol_pops) == 0L) {
      stop("orient = 'parents' requires aqu_pops and pol_pops present in prep")
    }
    f_aqu_par <- colMeans(P[aqu_pops, , drop = FALSE], na.rm = TRUE)
    f_pol_par <- colMeans(P[pol_pops, , drop = FALSE], na.rm = TRUE)
    parent_delta <- f_aqu_par - f_pol_par         # + => dosage allele is aqu
    sign_aqu <- sign(parent_delta)                # +1 / -1 / 0 / NA
    ## parental frequencies of the *aquilonia* allele (after orientation)
    f_aqu_parent <- ifelse(sign_aqu > 0, f_aqu_par, 1 - f_aqu_par)
    f_pol_parent <- ifelse(sign_aqu > 0, f_pol_par, 1 - f_pol_par)
  } else {
    parent_delta <- rep(NA_real_, length(markers))
    sign_aqu <- rep(if (dosage_is_aqu) 1 else -1, length(markers))
    f_aqu_parent <- rep(NA_real_, length(markers))
    f_pol_parent <- rep(NA_real_, length(markers))
  }
  parent_diff <- abs(parent_delta)                # between-species differentiation
  ## expected aquilonia-allele frequency in the founding admixture
  founding_f <- admix_prop * f_aqu_parent + (1 - admix_prop) * f_pol_parent

  ## oriented aquilonia-allele frequency in each hybrid population
  P_hyb <- P[hybrid_pops, , drop = FALSE]         # pops x markers
  flip <- which(sign_aqu < 0)
  undef <- which(is.na(sign_aqu) | sign_aqu == 0)

  f_aqu_hyb <- P_hyb
  if (length(flip)) f_aqu_hyb[, flip] <- 1 - P_hyb[, flip]
  if (length(undef)) f_aqu_hyb[, undef] <- NA_real_

  ## ---- per-population near-fixation calls ----
  is_aqu <- f_aqu_hyb >= (1 - fix_th)
  is_pol <- f_aqu_hyb <= fix_th

  n_obs   <- colSums(!is.na(f_aqu_hyb))
  n_aqu   <- colSums(is_aqu, na.rm = TRUE)
  n_pol   <- colSums(is_pol, na.rm = TRUE)
  n_fixed <- n_aqu + n_pol
  n_unsorted <- n_obs - n_fixed

  f_aqu_pooled <- colMeans(f_aqu_hyb, na.rm = TRUE)

  # ## ---- null probability ----
  # if (identical(null_prob, "pooled")) {
  #   p0 <- f_aqu_pooled
  # } else if (length(null_prob) == 1L) {
  #   p0 <- rep(as.numeric(null_prob), length(markers))
  # } else {
  #   if (length(null_prob) != length(markers)) {
  #     stop("null_prob must be a scalar, 'pooled', or length ncol(prep$pop_means)")
  #   }
  #   p0 <- as.numeric(null_prob)
  # }

  ## ---- ancestry-diagnostic gate ----
  ## Near-fixation is only evidence of SORTING if the locus was segregating
  ## in the founding admixture, i.e. the parents differ. Loci that are near-
  ## monomorphic everywhere (low DI / low parent_diff) would otherwise let
  ## every hybrid pop trivially "near-fix" the same allele -> spurious 20/20
  ## parallelism. Primary gate is min_DI (matches the DIEM-based pipeline and
  ## the DI exploration axis); min_parent_diff is an optional mechanistic
  ## safeguard (default off). Both are applied if both are set.
  differentiated <- rep(TRUE, length(markers))
  if (!is.null(min_parent_maf)) {
    if (is.null(parent_maf)) stop("min_parent_maf requires parent_maf (a marker-named pooled-parental MAF vector)")
    differentiated <- differentiated & !is.na(parent_maf_vec) & parent_maf_vec >= min_parent_maf
  }
  if (!is.null(min_DI)) {
    if (is.null(DI)) stop("min_DI requires DI (a marker-named DiagnosticIndex vector)")
    differentiated <- differentiated & !is.na(DI_vec) & DI_vec >= min_DI
  }
  if (min_parent_diff > 0 && orient == "parents") {
    differentiated <- differentiated & parent_diff >= min_parent_diff
  }
  differentiated[is.na(differentiated)] <- FALSE
  if (is.null(min_parent_maf) && is.null(min_DI) && !(min_parent_diff > 0)) {
    warning("no polymorphism/ancestry gate active (min_parent_maf, min_DI, min_parent_diff all unset): ",
            "near-monomorphic loci can produce spurious parallelism")
  }
  # 
  # ## ---- binomial test ----
  # run <- n_fixed >= min_fixed & !is.na(p0) & p0 > 0 & p0 < 1 & differentiated
  # p_binom <- rep(NA_real_, length(markers))
  # p_binom[run] <- .binom_two_sided_p(n_aqu[run], n_fixed[run], p0[run])
  # 
  # q_binom <- rep(NA_real_, length(markers))
  # q_binom[run] <- p.adjust(p_binom[run], method = "BH")

  ## ---- direction / magnitude summaries ----
  dir_bias <- ifelse(n_fixed > 0, (n_aqu - n_pol) / n_fixed, NA_real_)
  direction <- ifelse(
    is.na(dir_bias), NA_character_,
    ifelse(dir_bias > 0, "aquilonia",
           ifelse(dir_bias < 0, "polyctena", "tie"))
  )
  prop_fixed <- ifelse(n_obs > 0, n_fixed / n_obs, NA_real_)

  ## ---- population-fraction decomposition (all on one [0,1] scale) ----
  ## prop_fixed = uni_score + bi_score, so a single threshold sort_th on
  ## these commensurable fractions classifies every unit (SNP or eMLG)
  ## identically and lets genome proportions be simple weighted tallies.
  uni_score <- ifelse(n_obs > 0, (n_aqu - n_pol) / n_obs, NA_real_)   # signed, + = aquilonia
  bi_score  <- ifelse(n_obs > 0, 2 * pmin(n_aqu, n_pol) / n_obs, NA_real_)

  uni_mag <- abs(uni_score)
  ## two-sided binomial p that the fixation DIRECTION is unbiased under the
  ## random-direction null k_aqu ~ Binomial(n_fixed, 0.5). Vectorised exact form
  ## (symmetric null): twice the lower tail at the minority count, capped at 1;
  ## NA where nothing fixed. Used only by sort_rule = "binom" but always reported.
  p_binom <- ifelse(n_fixed > 0,
                    pmin(1, 2 * stats::pbinom(pmin(n_aqu, n_pol), n_fixed, 0.5)),
                    NA_real_)

  sort_class <- rep(NA_character_, length(markers))
  ok <- differentiated & n_obs > 0 & !is.na(uni_score)
  sort_class[ok] <- "unsorted"
  is_uni <- is_bi <- is_amb <- logical(length(markers))
  ## Magnitude gate + direction call. All three share the fixation counts:
  ##   "component" (ORIGINAL): the larger of |uni_score|/bi_score must itself reach
  ##     sort_th. A both-directions split shrinks BOTH, so near-balanced loci need
  ##     prop_fixed ~ 1 to pass -> "bidirectional" is systematically under-called.
  ##   "prop_fixed" (default): magnitude gate = TOTAL near-fixation
  ##     (prop_fixed = uni_mag + bi_score) >= sort_th; direction = the same uni/bi
  ##     split (bidirectional iff the minority > 1/4 of fixed pops). Removes the
  ##     "valley" but the 1/4 split is a fixed ratio, not sample-size aware.
  ##   "binom": magnitude gate = prop_fixed >= sort_th; DIRECTION is decided by the
  ##     binomial random-direction test, so "unidirectional" means the split is
  ##     SIGNIFICANTLY biased toward one parent (predictable), not merely a majority.
  ##     Sorted-but-not-significant loci are "bidirectional" when n_fixed had power
  ##     to detect a one-way bias (n_fixed >= n_pow, the smallest n at which an
  ##     all-one-way split reaches p < alpha), else "ambiguous" (underpowered).
  if (sort_rule == "component") {
    is_uni <- ok & uni_mag >= bi_score & uni_mag >= sort_th
    is_bi  <- ok & bi_score >  uni_mag & bi_score >= sort_th
  } else if (sort_rule == "prop_fixed") {
    sorted <- ok & prop_fixed >= sort_th
    is_bi  <- sorted & bi_score > uni_mag
    is_uni <- sorted & !is_bi
  } else {                                   # "binom"
    n_pow  <- ceiling(log(alpha / 2) / log(0.5))     # smallest n_fixed with 2*0.5^n < alpha
    sorted <- ok & prop_fixed >= sort_th
    sig    <- sorted & !is.na(p_binom) & p_binom < alpha
    is_uni <- sig
    is_bi  <- sorted & !sig & n_fixed >= n_pow
    is_amb <- sorted & !sig & n_fixed <  n_pow
  }
  sort_class[is_uni] <- ifelse(uni_score[is_uni] > 0, "aquilonia", "polyctena")
  sort_class[is_bi]  <- "bidirectional"
  sort_class[is_amb] <- "ambiguous"

  out <- data.table(
    marker       = markers,
    DI           = DI_vec,
    n_obs        = n_obs,
    n_aqu        = n_aqu,
    n_pol        = n_pol,
    n_fixed      = n_fixed,
    n_unsorted   = n_unsorted,
    prop_fixed   = prop_fixed,
    f_aqu_pooled = f_aqu_pooled,
    f_aqu_parent = f_aqu_parent,
    f_pol_parent = f_pol_parent,
    parent_diff  = parent_diff,
    parent_maf   = parent_maf_vec,
    founding_f   = founding_f,
    differentiated = differentiated,
    dir_bias     = dir_bias,
    direction    = direction,
    uni_score    = uni_score,
    bi_score     = bi_score,
    p_binom      = p_binom,
    sort_class   = sort_class
  )
  setkey(out, marker)
  out[]
}


## ---------------------------------------------------------
## Example usage (not run)
## ---------------------------------------------------------
if (FALSE) {

  library(data.table)
  devtools::load_all("~/gitlab/LDscnR/")  # ohta_fast_prepare() (LDscnR package)
  source("dev/R/parallelism_stats.R")

  e <- new.env()
  load("./data/hybrids_and_parents_maf005.Rdata", envir = e)
  GTs         <- e$GTs_with_parents            # samples x markers, 0/1/2
  sample_data <- e$sample_data_with_parents    # Population, Sample_ID, ...
  map         <- e$map_hyb_005                 # marker, DiagnosticIndex, ...

  prep <- ohta_fast_prepare(GTs, pops = sample_data$Population)

  aqu_pops <- "aquilonia_parent"
  pol_pops <- "polyctena_parent"
  hybrid_pops <- setdiff(unique(sample_data$Population),c(aqu_pops, pol_pops))
  ## drop Sielva to mirror earlier analyses, if desired:
  ## hybrid_pops <- setdiff(hybrid_pops, "Sielva")
  # DIvec <- setNames(map$DiagnosticIndex, map$marker)
  # r <- parallelism_stats(prep, hyb, aqu, pol, DI = DIvec, min_DI = -15, sort_th = 0.5)
  # 
  
  par_res <- parallelism_stats(
    prep,
    hybrid_pops = hybrid_pops,
    aqu_pops    = aqu_pops,
    pol_pops    = pol_pops,
    DI          = setNames(map$DiagnosticIndex, map$marker),  # aligned by name
    min_DI      = -50,          # primary ancestry-diagnostic gate
    min_fixed = 5,
    sort_th  = 0.5,
    fix_th      = 0.15,
    null_prob   = 0.5           # try "pooled" to correct for admixture skew
  )

  ## join with the marker map / sorting_stats() output on `marker`
  map[,indx := .I]
  par_res <- map[par_res, on = "marker"]
  setorder(par_res,indx)
  par_res[differentiated==TRUE, .N, by=sort_class]
  tab <- par_res[differentiated==TRUE , .N, by=sort_class][order(-N)]
  tab[TRUE, pct:=round(100*N/sum(N),1)]
  print(tab)

  markers <- par_res[differentiated==TRUE & sort_class=="aquilonia" ,marker]
  image(t(GTs[,markers]))

  markers <- par_res[differentiated==TRUE & sort_class=="polyctena" ,marker]
  image(t(GTs[,markers]))
  
  markers <- par_res[differentiated==TRUE & sort_class=="bidirectional",marker]
  image(t(GTs[,markers]))
 
  ## example code for heatmaps
  # suppressPackageStartupMessages({
  #   library(data.table)
  #   library(magick)
  #   devtools::load_all("~/gitlab/LDscnR/", quiet = TRUE)
  # })
  # eMLG_5loci_0025 <- readRDS("./data/eMLG_5loci_0025_cM05.rds")
  # pruned_stage1   <- readRDS("/tmp/pruned_stage1_genome.rds")  # or wherever your Stage 1 result lives
  # load("./data/hybrids_only_maf005.Rdata")  # GTs_hybrids_005, map_hyb_005, sample_data
  # 
  # grp <- eMLG_5loci_0025$groups[group_id == "F3968"]
  # mk  <- grp$members[[1]]
  
  # pop_order    <- sample_data[order(Population), Sample_ID]
  # cl_by_marker <- setNames(pruned_stage1$map_snp$CL_id, pruned_stage1$map_snp$marker)
  # pop_by_ind   <- setNames(sample_data$Population, sample_data$Sample_ID)
  # ## custom colors: named vector, one entry per level actually present in
  # ## THIS group -- not every population/cluster genome-wide
  # pop_levels <- unique(pop_by_ind[pop_order])
  # pop_colors <- c(
  #   Aland = "#E41A1C", Bunkkeri = "#377EB8", Grundsund = "#4DAF4A",
  #   Heinamaki = "#984EA3", Hiivola = "#FF7F00", Jarvenpaa = "#FFFF33",
  #   Karsikas = "#A65628", Katiskoski = "#F781BF", Kummunmaki = "#999999",
  #   LangholmenR = "#66C2A5", LangholmenW = "#FC8D62", Nyrhispera74 = "#8DA0CB",
  #   Nyrhispera75 = "#E78AC3", Parikkala = "#A6D854", Pikkala = "#FFD92F",
  #   Sielva = "#E5C494", Svanvik1 = "#B3B3B3", Svanvik2 = "#1B9E77",
  #   Tvarminne = "#D95F02", Vuosaari = "#7570B3"
  # )[pop_levels]  # subset down to just this group's actual populations
  # 
  # cl_levels <- as.character(unique(cl_by_marker[mk]))
  # cl_colors <- setNames(grDevices::hcl.colors(length(cl_levels), "Spectral"), cl_levels)
  # 
  # ht <- plot_genotype_heatmap(
  #   GTs_hybrids_005[, mk],
  #   row_order = pop_order,
  #   col_annotation = cl_by_marker[mk], col_annotation_name = "Stage 1 cluster",
  #   col_annotation_colors = cl_colors, col_annotation_legend = FALSE,
  #   row_annotation = pop_by_ind[pop_order], row_annotation_name = "Population",
  #   row_annotation_colors = pop_colors,
  #   title = "Chr2: F3968",
  #   out_file = "./Figures/Chr2_F3968_heatmap.png"
  # )
  # 
}



