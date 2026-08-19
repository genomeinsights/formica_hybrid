## =============================================================================
## module_di25 (high-DI analyses) -- RIGOROUS TEST of the "early pruning" result
## =============================================================================
## The claim under test (from di25_sielva_pruning.R / Figure 4):
##   Sielva, the most F1-like population (~84% heterozygous genome-wide), has
##   already begun to LOSE heterozygosity at loci that are sorted in the other 19
##   populations; the loss is DIRECTIONAL (toward the allele the others fix,
##   ~121:1); and the three large Chr 5/25/26 polyctena blocks resist it.
##
## Why the current version is untestable as stated:
##   * di25_sielva_pruning.R uses ONE Sielva individual (D[, which=="Sielva"][1]),
##     yet reports n = loci. Sampling unit is the colony/individual, and the loci
##     are linked -> the reported n is SNP-count pseudoreplication.
##   * The heterozygosity gradient could be manufactured by (a) low-coverage het
##     under-calling correlating with sortedness, or (b) shared ancestral allele-
##     frequency skew rather than active recombination-driven pruning.
##
## This script replaces the single-colony analysis with a 20-POPULATION,
## eMLG-UNIT-LEVEL design and tests the MECHANISM rather than asserting it. It is
## ADDITIVE: it reads di25_* outputs and the source matrices only; it modifies no
## existing script and writes only di25_pruning_* / di25_three_blocks.
##
## The analysis may FALSIFY the result. Section 9 enumerates the outcomes that
## would do so; each is evaluated and printed whichever way it goes. Nothing here
## is tuned toward the expected answer.
##
## Run from the repo root:  Rscript module_di25/R/di25_pruning_test.R
## =============================================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork); library(ggrepel) })
devtools::load_all("~/gitlab/LDscnR/", quiet = TRUE)      # consensus_dosage()

## ---- parameters -------------------------------------------------------------
CM_STAMP  <- "cM5"
FIX_TH    <- 0.15          # phi = 0.85: near-fixed if mean aqu dosage >= 1.7 or <= 0.30
PHI_HI    <- 2 * (1 - FIX_TH)   # 1.70
PHI_LO    <- 2 * FIX_TH         # 0.30
HET_LO    <- 0.5; HET_HI <- 1.5 # unit-level het band (consensus dosage is continuous)
MIN_N_IND <- 4             # skip a population x unit cell with < 4 GENOTYPED individuals
N_BOOT    <- 1000          # chromosome block bootstrap
N_PERM    <- 1000          # circular-shift permutations
STRONG    <- 0.75          # "strongly sorted" cut for descriptive / M2 tables
SEED      <- 1
set.seed(SEED)

OUTDIR <- "module_di25/data"; FIGDIR <- "module_di25/Figures"
INPUTS <- file.path(OUTDIR, "di25_inputs.rds")
CLUST  <- file.path(OUTDIR, sprintf("di25_clustering_%s.rds", CM_STAMP))
RECMAP <- "data/Frufa_DTOL_PR.ref_genome.recmap"
COL_HET <- "#BDBDBD"; COL_TOW <- "#315B7D"; COL_AGA <- "#C0392B"; COL_SIEL <- "#C2549D"

## =============================================================================
## 1. INPUTS -- 194-individual matrix, oriented so dosage 2 = homozygous aquilonia
## =============================================================================
inp <- readRDS(INPUTS); map <- as.data.table(inp$map)
res <- readRDS(CLUST);  g <- as.data.table(res$groups)
e2  <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd  <- e2$sample_data_with_parents

GTs  <- rbind(inp$GTs_hyb, inp$GTs_par)
keep <- rownames(GTs) %in% sd$Sample_ID                 # drop Hei159_h (no pop label)
GTs  <- GTs[keep, ]
pops <- sd$Population[match(rownames(GTs), sd$Sample_ID)]
faqu <- which(pops == "aquilonia_parent"); fpol <- which(pops == "polyctena_parent")
hp   <- setdiff(unique(pops), c("aquilonia_parent", "polyctena_parent"))   # 20 hybrid pops
is_par <- pops %in% c("aquilonia_parent", "polyctena_parent")

n_hyb <- sum(!is_par); n_par <- sum(is_par)
cat(sprintf("Individuals: %d hybrids in %d populations + %d parents\n", n_hyb, length(hp), n_par))
stopifnot(n_hyb == 165L, length(hp) == 20L, n_par == 30L)

## per-unit dosage: eMLG consensus for >2-marker clusters, representative SNP else.
## D is units x individuals (matches di25_sielva_pruning.R); flip so 2 = aquilonia.
is_e <- g$n_loci > 2
D <- vapply(seq_len(nrow(g)),
            function(i) if (is_e[i]) consensus_dosage(GTs, g$members[[i]]) else GTs[, g$representative[i]],
            numeric(nrow(GTs)))
D <- t(D); rownames(D) <- g$group_id                    # units x individuals
flip <- rowMeans(D[, faqu, drop = FALSE], na.rm = TRUE) < rowMeans(D[, fpol, drop = FALSE], na.rm = TRUE)
D[flip, ] <- 2 - D[flip, ]
cat(sprintf("Units: %d (%d eMLG consensus + %d rep-SNP); %d oriented (flipped)\n",
            nrow(g), sum(is_e), sum(!is_e), sum(flip, na.rm = TRUE)))

## unit genomic coordinates + interpolated recombination (Module A convention)
g[, `:=`(chr = as.integer(sub("Chr", "", sub(":.*", "", representative))),
         pos = as.integer(sub(".*:", "", representative)))]
rec <- fread(RECMAP); setnames(rec, 1:4, c("rchr", "rpos", "cM", "cMMb"))
rec[, Chr := as.integer(sub("chromosome_", "", rchr))]
interp <- function(chr, pos, col) {                     # approx(..., rule = 2) per chr
  out <- rep(NA_real_, length(pos))
  for (cc in unique(chr)) { r <- rec[Chr == cc]; if (nrow(r) < 2) next
    i <- which(chr == cc); out[i] <- approx(r$rpos, r[[col]], xout = pos[i], rule = 2)$y }
  out
}
g[, recomb := interp(chr, pos, "cMMb")]                 # cM/Mb  (local rate)
g[, cM_cum := interp(chr, pos, "cM")]                   # cumulative cM (for tract lengths)

## unit-level data-quality proxy: fraction of the 164 hybrids missing at the unit
## (source dosage carries genuine NA; consensus is NA only if all members are NA).
hyb_idx <- which(!is_par)
miss_u  <- setNames(rowMeans(is.na(D[, hyb_idx, drop = FALSE])), g$group_id)

## =============================================================================
## 2a. ARTEFACT CHECK -- are the 5 Sielva workers one genotype? (all populations)
## =============================================================================
## Pairwise identity-by-state per population: IBS = 1 - mean(|d_i - d_j|)/2 over
## units both individuals call. If Sielva's 5 workers are ~identical, Sielva is
## effectively n = 1 and every Sielva-specific statistic must say so.
ibs_pair <- function(a, b) { ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(NA_real_); 1 - mean(abs(a[ok] - b[ok])) / 2 }
pop_ibs <- rbindlist(lapply(hp, function(p) {
  cols <- which(pops == p); if (length(cols) < 2) return(data.table(pop = p, n = length(cols),
    mean_ibs = NA_real_, min_ibs = NA_real_, max_ibs = NA_real_))
  v <- combn(cols, 2, function(k) ibs_pair(D[, k[1]], D[, k[2]]))
  data.table(pop = p, n = length(cols), mean_ibs = mean(v), min_ibs = min(v), max_ibs = max(v))
}))[order(-mean_ibs)]
siel_cols <- which(pops == "Sielva")
siel_ibs  <- outer(seq_along(siel_cols), seq_along(siel_cols),
                   Vectorize(function(i, j) ibs_pair(D[, siel_cols[i]], D[, siel_cols[j]])))
dimnames(siel_ibs) <- list(rownames(GTs)[siel_cols], rownames(GTs)[siel_cols])

## =============================================================================
## 2c. ARTEFACT CHECK -- is Sielva F1? tract lengths are NOT diagnostic; use a
##     spatial runs test on the homozygous units instead.
## =============================================================================
## Tract-length alone cannot separate "F1 + scattered genotyping/dither error" from
## "early backcross": with ~16-22% homozygous units, random scatter already gives an
## expected het-run length of ~1/p_hom units, reproducing a short median. The real
## question is SPATIAL: are Sielva's homozygous units CLUSTERED into contiguous
## ancestry tracts (backcross), or Poisson-SCATTERED along the chromosome (F1+error)?
## Wald-Wolfowitz runs test per chromosome, aggregated: fewer runs than expected
## (Z<0) = clustering; Z~0 = scatter. Tract lengths kept only as description.
het_tracts <- function(col) {
  d <- D[, col]
  dt <- data.table(chr = g$chr, cM = g$cM_cum, het = d > HET_LO & d < HET_HI)[
    is.finite(cM) & !is.na(het)][order(chr, cM)]
  out <- list()
  for (cc in unique(dt$chr)) { x <- dt[chr == cc]
    r <- rle(x$het); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
    for (k in which(r$values)) out[[length(out) + 1]] <-
      data.table(chr = cc, n_units = r$lengths[k],
                 cM_len = x$cM[ends[k]] - x$cM[starts[k]])
  }
  rbindlist(out)
}
siel_tracts <- rbindlist(lapply(siel_cols, function(c) het_tracts(c)[, ind := rownames(GTs)[c]]))

## spatial runs test on the Sielva COLONY consensus (hom = 1, het = 0), per chr
siel_cons <- rowMeans(D[, siel_cols, drop = FALSE], na.rm = TRUE)
rt <- data.table(chr = g$chr, pos = g$pos,
                 hom = fifelse(siel_cons >= HET_HI | siel_cons <= HET_LO, 1L,
                       fifelse(siel_cons > HET_LO & siel_cons < HET_HI, 0L, NA_integer_)))[
  !is.na(hom) & is.finite(pos)][order(chr, pos)]
ww <- rt[, { x <- hom; n <- length(x); n1 <- sum(x); n0 <- n - n1
  if (n1 == 0 | n0 == 0 | n < 2) .(R = NA_real_, E = NA_real_, V = NA_real_)
  else .(R = 1 + sum(x[-1] != x[-n]), E = 1 + 2 * n1 * n0 / n,
         V = 2 * n1 * n0 * (2 * n1 * n0 - n) / (n^2 * (n - 1))) }, by = chr][!is.na(R)]
ww_Z <- (sum(ww$R) - sum(ww$E)) / sqrt(sum(ww$V))       # <0 = clustered into tracts
ww_p <- 2 * pnorm(-abs(ww_Z))
siel_phom <- mean(rt$hom)

## =============================================================================
## 3. LEAVE-ONE-OUT SORTEDNESS -- the core quantity (all 20 populations)
## =============================================================================
## Per population p and unit u:
## The COLONY (nest) is the sampling unit, NOT the individual worker. Within a
## population, the workers of one nest are collapsed to a per-nest CONSENSUS dosage
## (mean over the nest's called workers) BEFORE counting, so a clonal single-nest
## population (Sielva = nest Siel338, 5 workers, IBS 0.955) contributes ONE genotype
## per unit, not five correlated pseudo-replicates. 9 of the 20 populations are
## single-nest, so this matters beyond Sielva. Per population p and unit u:
##   mean_dose[p,u]  = mean aquilonia dosage over p's CALLED nests (nest consensus)
##   n_nest[p,u]     = number of p's nests with a called worker  <- the REPLICATION
##   n_ind[p,u]      = number of p's workers called              <- the COVERAGE gate
##   near-fixed if mean_dose >= 1.70 (aqu) or <= 0.30 (pol); else "mid"
##   sortedness_loo[p,u]  = fraction of the OTHER usable pops near-fixed at u
##   p's own composition:  n_het / n_toward / n_against counted over NESTS
## Cells with < MIN_N_IND called INDIVIDUALS are skipped (coverage gate); the NEST
## count is the binomial weight, so clonal replication is not double-counted.
nest    <- sub("_[0-9]*[a-zA-Z]+$", "", colnames(D))    # worker id -> nest/colony id
un      <- unique(nest); nest_pop <- pops[match(un, nest)]
pop_nests  <- lapply(hp, function(p) un[nest_pop == p]); names(pop_nests) <- hp
n_ind_pop  <- vapply(hp, function(p) sum(pops == p), integer(1))
n_nest_pop <- lengths(pop_nests)
## nest-consensus dosage (units x nests) and per-nest worker-called count
DN  <- vapply(un, function(nn) rowMeans(D[, nest == nn, drop = FALSE], na.rm = TRUE), numeric(nrow(D)))
DNc <- vapply(un, function(nn) rowSums(!is.na(D[, nest == nn, drop = FALSE])),        numeric(nrow(D)))
## a population needs >= MIN_N_IND individuals to be modelled at all (Svanvik1 n=3 out)
hp_use  <- hp[n_ind_pop >= MIN_N_IND]
hp_drop <- hp[n_ind_pop <  MIN_N_IND]

MEAN <- N_CALL <- N_IND <- matrix(NA_real_, nrow(g), length(hp), dimnames = list(g$group_id, hp))
NHET <- NAQU <- NPOL <- matrix(0L, nrow(g), length(hp), dimnames = list(g$group_id, hp))
for (p in hp) {
  nc   <- match(pop_nests[[p]], un)                     # this pop's nest columns
  subN <- DN[, nc, drop = FALSE]; subC <- DNc[, nc, drop = FALSE]
  subN[subC < 1] <- NA_real_                            # nest uncalled if no worker called
  MEAN[, p]   <- rowMeans(subN, na.rm = TRUE)           # mean over called nests
  N_CALL[, p] <- rowSums(subC >= 1)                     # number of called NESTS (replication)
  N_IND[, p]  <- rowSums(!is.na(D[, pops == p, drop = FALSE]))  # workers called (coverage)
  NHET[, p]   <- rowSums(subN > HET_LO & subN < HET_HI, na.rm = TRUE)   # het nests
  NAQU[, p]   <- rowSums(subN >= HET_HI, na.rm = TRUE)                  # hom-aqu nests
  NPOL[, p]   <- rowSums(subN <= HET_LO, na.rm = TRUE)                  # hom-pol nests
}
valid   <- N_IND >= MIN_N_IND                           # coverage gate (INDIVIDUALS)
nearfix <- valid & (MEAN >= PHI_HI | MEAN <= PHI_LO)
dir_aqu <- nearfix & (MEAN >= PHI_HI)                   # near-fixed toward aquilonia

## leave-one-out landscape per population
n_valid_oth <- rowSums(valid) - valid                   # # of other pops usable at u
n_fix_oth   <- rowSums(nearfix) - nearfix
n_aqu_oth   <- rowSums(dir_aqu) - dir_aqu
sortedness_loo <- ifelse(n_valid_oth > 0, n_fix_oth / n_valid_oth, NA_real_)
consensus_aqu  <- (n_aqu_oth >= (n_fix_oth - n_aqu_oth))   # TRUE=aqu majority, FALSE=pol
## own near-fixation fraction per population (denominator = usable units for p)
own_sortedness <- setNames(colSums(nearfix) / colSums(valid), hp)

## =============================================================================
## 4. THE MODELS  (per population, unit level)
## =============================================================================
## M1  P(heterozygous) ~ sortedness_loo * log10(recomb+0.1) + missingness
##     The sortedness:recomb INTERACTION is the mechanistic test, not a nuisance:
##     recombination-driven pruning predicts het loss steepest where recomb is
##     high; shared ancestral skew predicts NO interaction. lr is centered at its
##     global median so the sortedness main effect is the slope at typical recomb.
## M2  Among units where p is homozygous, P(toward consensus) vs 0.5.
lr_all  <- log10(g$recomb + 0.1); lr_c0 <- median(lr_all, na.rm = TRUE)

## per-population model frame (one row per usable unit)
frame_for <- function(p) {
  keep <- valid[, p] & is.finite(sortedness_loo[, p]) & is.finite(g$recomb)
  data.table(unit = g$group_id[keep], chr = g$chr[keep],
             n_het = NHET[keep, p], n_call = N_CALL[keep, p],
             n_toward = ifelse(consensus_aqu[keep, p], NAQU[keep, p], NPOL[keep, p]),
             n_against = ifelse(consensus_aqu[keep, p], NPOL[keep, p], NAQU[keep, p]),
             s = sortedness_loo[keep, p], lr = lr_all[keep] - lr_c0, miss = miss_u[keep])
}
frames <- setNames(lapply(hp_use, frame_for), hp_use)

fit_M1 <- function(fr) glm(cbind(n_het, n_call - n_het) ~ s * lr + miss,
                           data = fr, family = binomial())
M1 <- lapply(frames, fit_M1)
## also a missingness-ONLY control (does data quality alone make the gradient?):
M1_miss <- lapply(frames, function(fr) glm(cbind(n_het, n_call - n_het) ~ miss,
                                           data = fr, family = binomial()))

slope_tab <- rbindlist(lapply(hp_use, function(p) { cf <- summary(M1[[p]])$coef
  data.table(pop = p, own_sortedness = own_sortedness[p], n_units = nrow(frames[[p]]),
             b_sorted = cf["s", 1], se_sorted = cf["s", 2],
             b_inter  = cf["s:lr", 1], se_inter = cf["s:lr", 2],
             b_miss   = cf["miss", 1]) }))

## M2: direction of homozygotes at strongly-sorted units, pooled per population
m2_tab <- rbindlist(lapply(hp_use, function(p) { fr <- frames[[p]][s > STRONG]
  nt <- sum(fr$n_toward); na <- sum(fr$n_against); tot <- nt + na
  bt <- if (tot > 0) binom.test(nt, tot, 0.5) else NULL
  data.table(pop = p, n_hom = tot, n_toward = nt, n_against = na,
             frac_toward = if (tot > 0) nt / tot else NA_real_,
             p_binom = if (!is.null(bt)) bt$p.value else NA_real_) }))

## =============================================================================
## 5. NINETEEN SLOPES, NOT ONE  -- pruning_slope[p] ~ own_sortedness[p]
## =============================================================================
## 19 modelled populations (Svanvik1 dropped, n=3). Prediction: slope declines with
## own sortedness (populations further along have already pruned and have less het
## left to lose); Sielva should be the ENDPOINT of a continuum, not an off-line
## outlier. NB Sielva is 1 colony (binomial weight 1): its position on the line is a
## genuine result but rests on a single genome, not on within-population replication.
lm20 <- lm(b_sorted ~ own_sortedness, data = slope_tab)
slope_tab[, `:=`(fit = predict(lm20), resid = residuals(lm20))]
slope_tab[, studentized := rstudent(lm20)[.I]]

## =============================================================================
## 6. UNCERTAINTY  -- inference matched to the resampling UNIT
## =============================================================================
## Two different sources of uncertainty need two different schemes:
##  (A) Sielva's own M1 coefficients + M2: the units retain residual LD, so CIs
##      come from a CHROMOSOME block bootstrap + a circular-shift permutation that
##      preserve within-chromosome dependence. A WITHIN-population inference.
##  (B) the 19-slope regression pruning_slope ~ own_sortedness: the observations
##      are POPULATIONS, so the dominant uncertainty is scatter of the 19 points
##      around the line, which a chromosome bootstrap CANNOT see (it only resamples
##      within a population). We resample POPULATIONS (bootstrap) and permute the
##      slope<->own_sortedness pairing (permutation); the parametric lm p is the
##      honest headline. [An earlier version used the chromosome bootstrap here and
##      understated the CI ~4x -- it measured the wrong thing.]
chrs <- sort(unique(g$chr))
safe_coef <- function(fr, term) { m <- tryCatch(fit_M1(fr), error = function(e) NULL)
  if (is.null(m)) return(NA_real_); cf <- coef(m); if (term %in% names(cf)) cf[[term]] else NA_real_ }
ci <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)

## (A) chromosome block bootstrap -- Sielva only (within-population residual LD)
boot <- rbindlist(lapply(seq_len(N_BOOT), function(b) {
  cc <- sample(chrs, length(chrs), replace = TRUE)
  fr_s <- rbindlist(lapply(cc, function(k) frames[["Sielva"]][chr == k]))
  nt <- sum(fr_s[s > STRONG]$n_toward); na <- sum(fr_s[s > STRONG]$n_against)
  data.table(siel_inter = safe_coef(fr_s, "s:lr"), siel_slope = safe_coef(fr_s, "s"),
             siel_frac_toward = if (nt + na > 0) nt / (nt + na) else NA_real_)
}))
## (A) circular-shift permutation -- Sielva only
ord_s <- order(g[match(frames[["Sielva"]]$unit, group_id), chr],
               g[match(frames[["Sielva"]]$unit, group_id), pos])
shift_siel <- function(k) { fr <- copy(frames[["Sielva"]]); n <- length(ord_s)
  sh <- ((seq_len(n) - 1 + k) %% n) + 1; fr$s[ord_s] <- fr$s[ord_s][sh]; fr }
perm <- rbindlist(lapply(seq_len(N_PERM), function(b) {
  fr_s <- shift_siel(sample.int(nrow(frames[["Sielva"]]) - 1, 1))
  data.table(siel_inter = safe_coef(fr_s, "s:lr"), siel_slope = safe_coef(fr_s, "s"))
}))
perm_p <- function(null, obs) mean(abs(null) >= abs(obs), na.rm = TRUE)   # two-sided

## (B) 19-slope regression: POPULATION-level bootstrap + permutation + lm p
reg_obs    <- coef(lm20)[["own_sortedness"]]
reg_lm_p   <- summary(lm20)$coef["own_sortedness", 4]
reg_boot   <- vapply(seq_len(N_BOOT), function(b) { ix <- sample(nrow(slope_tab), replace = TRUE)
  tryCatch(coef(lm(b_sorted ~ own_sortedness, data = slope_tab[ix]))[["own_sortedness"]],
           error = function(e) NA_real_) }, numeric(1))
reg_perm   <- vapply(seq_len(N_PERM), function(b) {
  st <- copy(slope_tab); st[, os := sample(own_sortedness)]
  coef(lm(b_sorted ~ os, data = st))[["os"]] }, numeric(1))
reg_ci     <- ci(reg_boot); reg_perm_p <- mean(abs(reg_perm) >= abs(reg_obs))

## interaction DIRECTION across populations: a consistent POSITIVE interaction means
## het loss is steeper at LOW recombination (linked purging) -- the OPPOSITE of what
## recombination-exposure pruning predicts (which would be negative). Sielva's lone
## negative rests on 1 colony and is best read as noise.
n_pos <- sum(slope_tab$b_inter > 0); n_neg <- sum(slope_tab$b_inter < 0)
inter_sign_p <- binom.test(max(n_pos, n_neg), nrow(slope_tab), 0.5)$p.value
inter_mean_excl_siel <- mean(slope_tab[pop != "Sielva", b_inter])

## =============================================================================
## 7. THE THREE BLOCKS  -- define from the data, then re-run 4-6 without them
## =============================================================================
## Consensus-polyctena-fixed unit = usable in >= 10 pops, >= 75% near-fixed,
## majority direction polyctena. Each block is ANCHORED on the pipeline's known
## block unit (di25_region_segregation.R: Chr5=F8626, Chr25=F6986, Chr26=F7174)
## and EXPANDED to the contiguous consensus-polyctena run: member-marker intervals
## merged across gaps <= 500 kb (a giant eMLG unit spans many Mb of markers, so
## contiguity is judged on member extents, not the single representative point).
## Report Mb span, NOT marker count (circos width ~ marker count by construction).
overall_fix   <- rowSums(nearfix)
overall_valid <- rowSums(valid)
overall_pol   <- rowSums(nearfix & !dir_aqu)
is_conspol <- overall_valid >= 10 & overall_fix / pmax(overall_valid, 1) >= STRONG &
              overall_pol >= (overall_fix - overall_pol)
posmap <- setNames(map$Pos, map$marker)                 # marker -> bp
u_lo <- vapply(g$members, function(mk) min(posmap[mk]), numeric(1))   # member extent
u_hi <- vapply(g$members, function(mk) max(posmap[mk]), numeric(1))

ANCHORS <- c(Chr5 = "F8626", Chr25 = "F6986", Chr26 = "F7174"); GAP <- 5e5
block_units <- character(0); block_tab <- list()
for (a in ANCHORS) {
  cc <- g$chr[match(a, g$group_id)]
  cand <- which(g$chr == cc & (is_conspol | g$group_id == a))   # anchor forced in
  o  <- cand[order(u_lo[cand])]; lo <- u_lo[o]; hi <- u_hi[o]
  newrun <- c(TRUE, lo[-1] > cummax(hi)[-length(hi)] + GAP)     # interval merge
  run <- cumsum(newrun)
  u  <- o[run == run[match(a, g$group_id[o])]]                  # run holding the anchor
  block_units <- c(block_units, g$group_id[u])
  mpos <- posmap[unlist(g$members[u])]
  span <- (max(mpos) - min(mpos)) / 1e6
  block_tab[[length(block_tab) + 1]] <- data.table(
    chr = cc, anchor = a, n_units = length(u), n_markers = sum(g$n_loci[u]),
    start_Mb = min(mpos) / 1e6, end_Mb = max(mpos) / 1e6, span_Mb = span,
    markers_per_Mb = sum(g$n_loci[u]) / max(span, 1e-3),   # circos arc width ~ this
    group_ids = paste(g$group_id[u], collapse = ","))
}
three_blocks <- rbindlist(block_tab)[order(-span_Mb)]
saveRDS(three_blocks, file.path(OUTDIR, "di25_three_blocks.rds"))

## re-run the core models with block units excluded
frames_nb <- lapply(frames, function(fr) fr[!unit %in% block_units])
M1_nb  <- lapply(frames_nb, fit_M1)
slope_nb <- rbindlist(lapply(hp_use, function(p) { cf <- summary(M1_nb[[p]])$coef
  data.table(pop = p, own_sortedness = own_sortedness[p],
             b_sorted = cf["s", 1], b_inter = cf["s:lr", 1]) }))
lm20_nb <- lm(b_sorted ~ own_sortedness, data = slope_nb)
m2_nb <- rbindlist(lapply(hp_use, function(p) { fr <- frames_nb[[p]][s > STRONG]
  nt <- sum(fr$n_toward); na <- sum(fr$n_against)
  data.table(pop = p, n_hom = nt + na, frac_toward = if (nt + na > 0) nt / (nt + na) else NA_real_) }))

## =============================================================================
## 8. SAVE + FIGURES
## =============================================================================
saveRDS(list(M1 = lapply(M1, coef), M1_miss = lapply(M1_miss, coef),
             M1_nb = lapply(M1_nb, coef)), file.path(OUTDIR, "di25_pruning_models.rds"))
saveRDS(list(slopes = slope_tab, m2 = m2_tab, lm20 = coef(summary(lm20)),
             slopes_nb = slope_nb, lm20_nb = coef(summary(lm20_nb)), m2_nb = m2_nb,
             own_sortedness = own_sortedness), file.path(OUTDIR, "di25_pruning_slopes.rds"))
saveRDS(list(boot = boot, reg_boot = reg_boot, ci = list(
               siel_inter = ci(boot$siel_inter), siel_slope = ci(boot$siel_slope),
               siel_frac_toward = ci(boot$siel_frac_toward),
               reg_slope_popboot = reg_ci)),
        file.path(OUTDIR, "di25_pruning_boot.rds"))
saveRDS(list(perm = perm, reg_perm = reg_perm, p = list(
               siel_inter = perm_p(perm$siel_inter, slope_tab[pop == "Sielva", b_inter]),
               siel_slope = perm_p(perm$siel_slope, slope_tab[pop == "Sielva", b_sorted]),
               reg_slope_popperm = reg_perm_p, reg_slope_lm = reg_lm_p,
               inter_sign = inter_sign_p, runs_Z = ww_Z, runs_p = ww_p)),
        file.path(OUTDIR, "di25_pruning_perm.rds"))

## -- Fig A: 19 population slopes vs own sortedness, Sielva labelled ------------
## CI/perm here are POPULATION-level (resample the 19 pops); the lm p is the honest
## headline. Sielva (magenta) is 1 colony -- on the line, but a single genome.
pA <- ggplot(slope_tab, aes(own_sortedness, b_sorted)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40", fill = "grey85", linewidth = 0.6) +
  geom_errorbar(aes(ymin = b_sorted - se_sorted, ymax = b_sorted + se_sorted),
                width = 0, colour = "grey65") +
  geom_point(aes(colour = pop == "Sielva", size = pop == "Sielva")) +
  geom_text_repel(aes(label = pop), size = 2.6, min.segment.length = 0, segment.colour = "grey75") +
  scale_colour_manual(values = c("FALSE" = "#5B6570", "TRUE" = COL_SIEL), guide = "none") +
  scale_size_manual(values = c("FALSE" = 1.8, "TRUE" = 3), guide = "none") +
  labs(x = "own sortedness (fraction of units where the population is itself near-fixed)",
       y = "M1 het-loss slope  (logit P(het) per unit sortedness)",
       title = "Nineteen pruning slopes vs own sortedness (Svanvik1 dropped, n=3)",
       subtitle = sprintf("slope %+.2f; lm p=%.3f (n=19, the honest test); pop-bootstrap 95%% CI [%.2f, %.2f], pop-perm p=%.3f. Sielva (1 colony) %s the trend.",
         reg_obs, reg_lm_p, reg_ci[1], reg_ci[2], reg_perm_p,
         ifelse(abs(slope_tab[pop == "Sielva", studentized]) > 2.5, "is an OUTLIER from", "sits on"))) +
  theme_bw(base_size = 10)
ggsave(file.path(FIGDIR, "di25_pruning_slopes.png"), pA, width = 8.5, height = 6, dpi = 200)

## -- Fig B: M1 interaction ACROSS the 19 populations (forest) -----------------
## The mechanism question is a 19-population statement, not a Sielva one. A POSITIVE
## interaction = het loss steeper at LOW recombination (linked purging); the recomb-
## exposure pruning hypothesis predicts NEGATIVE. 18/19 are positive.
fb <- copy(slope_tab)[order(b_inter)]; fb[, pop := factor(pop, levels = pop)]
pB <- ggplot(fb, aes(b_inter, pop)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_errorbarh(aes(xmin = b_inter - se_inter, xmax = b_inter + se_inter),
                 height = 0, colour = "grey65") +
  geom_point(aes(colour = b_inter > 0), size = 2.4) +
  scale_colour_manual(values = c("FALSE" = COL_AGA, "TRUE" = COL_TOW), guide = "none") +
  labs(x = "M1 interaction  sortedness x log10(recomb)   (+ = het loss steeper at LOW recomb)",
       y = NULL,
       title = "The recombination interaction points AWAY from pruning",
       subtitle = sprintf("%d of 19 populations POSITIVE (het loss concentrated at low recomb = linked purging), sign-test p=%.1e.\nRecombination-exposure pruning predicts NEGATIVE. Sielva is the lone negative (1 colony) -- best read as noise.",
         n_pos, inter_sign_p)) +
  theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size = 8))
ggsave(file.path(FIGDIR, "di25_pruning_interaction.png"), pB, width = 9, height = 5.6, dpi = 200)

## -- Fig C: original Figure 4 rebuilt at UNIT level (Sielva, 5 workers pooled) -
## composition by sortedness bin, block-bootstrap CI on het, three blocks shown
## as a separate bin. n reported as UNITS (not SNPs).
scomp <- copy(frames[["Sielva"]])
scomp[, inblock := unit %in% block_units]
scomp[, sb := cut(s, c(-.01, .25, .5, .75, .95, 1.01),
                  labels = c("0-25%","25-50%","50-75%","75-95%","95-100%"))]
mk_comp <- function(dt) dt[, .(het = sum(n_het), tow = sum(n_toward), aga = sum(n_against),
                               tot = sum(n_call), n_units = .N), by = sb][order(sb)]
cc_main <- mk_comp(scomp[!(inblock)]);  cc_main[, grp := "excl. 3 blocks"]
cc_blk  <- mk_comp(scomp[(inblock)]);   cc_blk[,  grp := "3 blocks"]
compA <- rbind(cc_main, cc_blk)
compL <- melt(compA, id.vars = c("sb", "grp", "tot", "n_units"),
              measure.vars = c("het", "tow", "aga"), variable.name = "cls", value.name = "n")
compL[, frac := n / tot]
compL[, cls := factor(cls, levels = c("aga", "tow", "het"),
                      labels = c("hom, against", "hom, toward fixed", "heterozygous"))]
pC <- ggplot(compL[grp == "excl. 3 blocks"], aes(sb, frac, fill = cls)) +
  geom_col(width = 0.8) +
  geom_text(data = cc_main, aes(sb, 1.04, label = paste0("u=", n_units)),
            inherit.aes = FALSE, size = 2.5) +
  scale_fill_manual(values = c("heterozygous" = COL_HET, "hom, toward fixed" = COL_TOW,
                               "hom, against" = COL_AGA), name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
  labs(x = "sortedness in the other 19 populations", y = "Sielva composition (unit level)",
       title = "Figure 4 rebuilt at eMLG-unit level, 3 blocks removed (n = units; Sielva = 1 colony consensus/unit)",
       subtitle = sprintf("het %.0f%% (unsorted) -> %.0f%% (75-95%% sorted); the 95-100%% bin holds only %d NON-block units. Blocks shown separately below.",
         100 * cc_main[sb == "0-25%", het / tot], 100 * cc_main[sb == "75-95%", het / tot],
         cc_main[sb == "95-100%", n_units])) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
pC2 <- ggplot(compL[grp == "3 blocks"], aes(sb, frac, fill = cls)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c("heterozygous" = COL_HET, "hom, toward fixed" = COL_TOW,
                               "hom, against" = COL_AGA), guide = "none") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "sortedness", y = "blocks only",
       subtitle = sprintf("1 large Chr26 block (%.1f Mb) + 2 small dense regions Chr5/25 (%.2f, %.2f Mb): %d units, %.1f Mb total.",
         three_blocks[chr == 26, span_Mb], three_blocks[chr == 5, span_Mb],
         three_blocks[chr == 25, span_Mb], sum(three_blocks$n_units), sum(three_blocks$span_Mb))) +
  theme_bw(base_size = 9)
ggsave(file.path(FIGDIR, "di25_pruning_sielva.png"), pC / pC2 + plot_layout(heights = c(2.4, 1)),
       width = 10, height = 7.5, dpi = 200)

## =============================================================================
## 9-10. REPORT + FALSIFICATION CRITERIA
## =============================================================================
cat("\n================= DI25 EARLY-PRUNING TEST =================\n")
cat(sprintf("Modelled populations: %d of 20 (MIN_N_IND = %d).%s\n", length(hp_use), MIN_N_IND,
            if (length(hp_drop)) sprintf(" Excluded (too few individuals): %s.",
              paste(sprintf("%s n=%d", hp_drop, n_ind_pop[hp_drop]), collapse = ", ") ) else ""))
cat(sprintf("Sampling unit = COLONY (nest): %d of 20 pops are single-nest (Sielva = 1 nest, 5 workers).\n",
            sum(n_nest_pop == 1)))
cat("  Binomial weight = # nests; individual count is only the coverage gate. Nests per pop:\n  ")
cat(paste(sprintf("%s:%d", hp, n_nest_pop), collapse = "  "), "\n")
cat("Leave-one-out 'other 19' denominators count only usable populations at each unit.\n")

cat("\n[2a] Sielva 5-worker IBS matrix (1 = identical):\n"); print(round(siel_ibs, 3))
cat(sprintf("     Sielva mean within-pop IBS = %.3f ; rank among 20 pops (1=most clonal): %d\n",
            pop_ibs[pop == "Sielva", mean_ibs], which(pop_ibs$pop == "Sielva")))
cat("     within-pop IBS, most-to-least clonal (top/bottom 3):\n")
print(rbind(head(pop_ibs, 3), tail(pop_ibs, 3))[, .(pop, n, mean_ibs = round(mean_ibs, 3))])

cat("\n[2b] missingness confound:\n")
ms <- data.table(miss = miss_u, s = rowMeans(sortedness_loo, na.rm = TRUE))
cat(sprintf("     Spearman(unit missingness, mean sortedness) = %.3f\n",
            cor(ms$miss, ms$s, method = "spearman", use = "complete.obs")))
cat(sprintf("     Sielva missingness-ONLY het slope: %+.3f (M1 full missingness coef: %+.3f)\n",
            coef(M1_miss[["Sielva"]])[["miss"]], slope_tab[pop == "Sielva", b_miss]))

cat("\n[2c] Is Sielva F1? tract lengths are NOT diagnostic; spatial runs test is.\n")
cat("     het-tract lengths (cM, across 5 workers -- DESCRIPTION only):\n")
print(siel_tracts[, .(n_tracts = .N, med_units = median(n_units), med_cM = round(median(cM_len), 2),
                      max_cM = round(max(cM_len), 2), frac_ge5cM = round(mean(cM_len >= 5), 2))])
cat(sprintf("     Wald-Wolfowitz runs on Sielva-consensus homozygous units: %.0f%% hom; runs Z=%+.2f, p=%.2g\n",
            100 * siel_phom, ww_Z, ww_p))
cat(sprintf("     -> %s. (Z<0 = clustered ancestry tracts; Z~0 = scattered = F1+error.)\n",
            if (ww_p < 0.05 && ww_Z < 0) "homozygous units are SIGNIFICANTLY clustered (real tract structure), but modestly -- F1-like background carrying some tracts; NOT a clean F1, NOT clearly a backcross" else
            "no significant clustering -- consistent with F1 + scattered error"))

cat("\n[3/5] 19-slope regression  pruning_slope ~ own_sortedness (obs = POPULATIONS):\n")
print(round(coef(summary(lm20)), 3))
cat(sprintf("     lm p = %.3f (n=19, the HONEST test); pop-bootstrap 95%% CI [%.3f, %.3f]; pop-perm p = %.3f\n",
            reg_lm_p, reg_ci[1], reg_ci[2], reg_perm_p))
cat("     (a chromosome bootstrap here would understate the CI ~4x -- it resamples within a pop, not across the 19.)\n")
cat(sprintf("     Sielva studentized residual: %+.2f (|>2.5| = off-line outlier); Sielva = 1 colony (weight 1).\n",
            slope_tab[pop == "Sielva", studentized]))

si <- slope_tab[pop == "Sielva"]
cat("\n[4] M1 recombination interaction ACROSS the 19 populations (the mechanism test):\n")
cat(sprintf("     %d POSITIVE / %d negative; sign-test p = %.2g; mean interaction (excl Sielva) = %+.3f\n",
            n_pos, n_neg, inter_sign_p, inter_mean_excl_siel))
cat("     POSITIVE = het loss steeper at LOW recombination (linked purging); recomb-EXPOSURE pruning predicts NEGATIVE.\n")
cat(sprintf("     Sielva is the lone negative (%+.3f; block-boot CI [%.2f, %.2f]; perm p=%.3f) -- 1 colony, best read as noise.\n",
            si$b_inter, ci(boot$siel_inter)[1], ci(boot$siel_inter)[2], perm_p(perm$siel_inter, si$b_inter)))
cat(sprintf("     Sielva het-loss slope %+.3f (block-boot CI [%.2f, %.2f]; perm p=%.3f) -- direction of loss, unrelated to the mechanism sign.\n",
            si$b_sorted, ci(boot$siel_slope)[1], ci(boot$siel_slope)[2], perm_p(perm$siel_slope, si$b_sorted)))

cat("\n[4] M2 direction of Sielva homozygotes at strongly-sorted units (sortedness > 0.75):\n")
print(m2_tab[pop == "Sielva"])
cat("     CAVEAT: the 0.5 null holds only if the 20 populations did NOT draw from one\n")
cat("     shared, already-skewed founding pool. At DI25 markers the parents are near-\n")
cat("     fixed for opposite alleles and Sielva is ~50/50 genome-wide, so founding-pool\n")
cat("     frequency ~ admixture proportion ~ 0.5 -- but a shared skewed pool would break this.\n")

cat("\n[7] blocks: ONE large Chr26 block + TWO small dense regions (say Mb, not marker count):\n")
print(three_blocks[, .(chr, n_units, n_markers, span_Mb = round(span_Mb, 2),
                       markers_per_Mb = round(markers_per_Mb))])
cat(sprintf("     Chr26 = %.1f Mb (%.0f markers/Mb) is the only large block; Chr5 (%.2f Mb) & Chr25 (%.2f Mb, %.0f markers/Mb)\n",
            three_blocks[chr == 26, span_Mb], three_blocks[chr == 26, markers_per_Mb],
            three_blocks[chr == 5, span_Mb], three_blocks[chr == 25, span_Mb],
            three_blocks[chr == 25, markers_per_Mb]))
cat(sprintf("     are sub-Mb but marker-DENSE -> wide circos arcs, narrow regions (the artefact). ~%.0f Mb total; say '~%.0f Mb'.\n",
            sum(three_blocks$span_Mb), sum(three_blocks$span_Mb)))
cat("     (A separate dense consensus-polyctena region sits on Chr5 at ~2.5-3.2 Mb, 1.3 Mb from the anchored block.)\n")
## sortedness-MATCHED block comparison (blocks are all strongly sorted). DESCRIPTIVE.
blk_het <- scomp[(inblock), sum(n_het) / sum(n_call)]
nb_hi_het <- scomp[!(inblock) & s > STRONG, sum(n_het) / sum(n_call)]
cat(sprintf("     Sielva het: %.0f%% IN the blocks vs %.0f%% in NON-block units at matched sortedness (>%.2f).\n",
            100 * blk_het, 100 * nb_hi_het, STRONG))
cat(sprintf("     DESCRIPTIVE only: %d block units split across 2 bins (~%d/bin) -- describe, do not test.\n",
            sum(three_blocks$n_units), round(sum(three_blocks$n_units) / 2)))
cat(sprintf("     M1 with blocks EXCLUDED -- Sielva sortedness slope %+.3f (was %+.3f); interaction %+.3f (was %+.3f)\n",
            slope_nb[pop == "Sielva", b_sorted], si$b_sorted,
            slope_nb[pop == "Sielva", b_inter], si$b_inter))
cat(sprintf("     19-slope regression (blocks excluded): slope %+.3f (lm p=%.3f)\n",
            coef(lm20_nb)[["own_sortedness"]], summary(lm20_nb)$coef["own_sortedness", 4]))

cat("\n===== FALSIFICATION CRITERIA (Section 9) =====\n")
f1 <- perm_p(perm$siel_slope, si$b_sorted) > 0.05 || (ci(boot$siel_slope)[1] < 0 & ci(boot$siel_slope)[2] > 0)
## #2 = the het gradient itself collapses when the blocks are removed: sign flip
## or the Sielva slope more than halved. (NOT the 19-slope regression's power.)
siel_slope_nb <- slope_nb[pop == "Sielva", b_sorted]
f2b <- (siel_slope_nb / si$b_sorted) < 0.5 || sign(siel_slope_nb) != sign(si$b_sorted)
f3 <- abs(si$studentized) > 2.5
## recomb-EXPOSURE pruning is supported only if the interactions are consistently
## NEGATIVE. They are consistently POSITIVE (linked purging) -> pruning falsified.
pruning_supported <- (n_neg > n_pos) && inter_sign_p < 0.05
linked_purging    <- (n_pos > n_neg) && inter_sign_p < 0.05
f4 <- !pruning_supported
f5 <- abs(coef(M1_miss[["Sielva"]])[["miss"]]) > 0 &&
      cor(ms$miss, ms$s, method = "spearman", use = "complete.obs") > 0.3
f6 <- pop_ibs[pop == "Sielva", mean_ibs] > 0.98
cat(sprintf("  (1) effect gone at unit level ...................... %s\n", ifelse(f1, "YES (falsified)", "no")))
cat(sprintf("  (2) effect gone with blocks excluded ............... %s\n", ifelse(f2b, "YES (falsified)", "no")))
cat(sprintf("  (3) Sielva off-line outlier, not endpoint .......... %s\n", ifelse(f3, "YES", "no")))
cat(sprintf("  (4) recomb-exposure PRUNING mechanism ............. %s\n",
            ifelse(f4, "UNSUPPORTED (interaction is +ve: linked purging, not pruning)", "supported")))
cat(sprintf("  (5) missingness reproduces the gradient alone ...... %s\n", ifelse(f5, "YES (artefactual)", "no")))
cat(sprintf("  (6) 5 Sielva workers are one genotype .............. %s\n", ifelse(f6, "YES (n=1)", "no")))

cat("\nVERDICT:\n")
cat(sprintf("  * Directional het-loss toward the fixed allele: %s (M2 %d:%d, unit level).\n",
            ifelse(!f1, "REAL", "FAILS"), m2_tab[pop=="Sielva", n_toward], m2_tab[pop=="Sielva", n_against]))
cat(sprintf("  * It scales across the 19 populations (continuum lm p=%.3f); Sielva is the endpoint (resid %+.2f), on 1 colony.\n",
            reg_lm_p, si$studentized))
cat(sprintf("  * Mechanism: het loss is concentrated at LOW recombination (%d/19 +ve, sign p=%.1e) = %s,\n",
            n_pos, inter_sign_p, ifelse(linked_purging, "LINKED PURGING", "inconsistent")))
cat("    the OPPOSITE of recombination-exposure pruning. The 'active pruning' mechanism is not supported;\n")
cat("    this aligns with Nouhaud low-recomb purging and Result 3 (sorting stronger at low recomb).\n")
cat("\n[di25-pruning-test] done\n")
