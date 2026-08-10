## =========================================================================
## DI25 among-region association -- STEP 0: assemble the unit dosage matrix + gate.
## =========================================================================
## Module D's among-region screen consumes (a) a continuous unit dosage matrix
## (individuals x LD-reduced units, oriented so 2 = F. aquilonia), (b) a per-unit
## gate table (differentiated / DI / sort_class / prop_fixed / uni_score), and
## (c) a genetic-map cM per unit. Module D read these from the CANONICAL all-marker
## clustering (eMLG_5loci_0025_cM05.rds) + the moduleC gate (moduleC_C3_cl.rds).
##
## Here we rebuild the SAME contract from the HIGH-DI-ONLY pipeline, so every unit
## summarises high-DI variation exclusively:
##   * units          <- di25_clustering_cM5.rds  (from-scratch DI25 LD clustering)
##   * gate           <- di25_sorting_emlg.rds     (a SUPERSET of the moduleC gate)
##   * dosages         : eMLG consensus (>2 markers) or representative SNP (1-2), as in
##                       diem_circos_di25_eMLG.R, oriented to F. aquilonia via parents.
##
## Individuals: rbind(GTs_hyb, GTs_par) = 195. sample_data_with_parents lists only 194
## (one di25 hybrid, Hei159_h, is absent). di25_sorting.R DROPPED it; here we KEEP the
## full 165-hybrid set and RECOVER its population from the colony prefix in the name
## (Hei159_h -> colony Hei159 -> Heinamaki). PC1/PC2 are not needed (the EMMAX arm builds
## its own genome PCs from the dosage matrix), so the recovered row carries only a label.
## -> 165 hybrids (20 pops) + 30 parents. Parents orient each unit and are then set aside;
## the analysis matrix is the 165 hybrids.
##
## Writes ONE object, data/di25D_units.rds, consumed by every downstream di25D script
## (di25D_ohta_dmi.R, di25D_emmax.R, di25D_network_build.R, ...). Run from the repo
## root: Rscript module_di25/among_region_association/R/di25D_assemble_units.R
## =========================================================================

suppressMessages({ library(data.table) })
devtools::load_all("~/gitlab/LDscnR/")            # consensus_dosage() -- same as diem_circos_di25_eMLG.R
set.seed(1)

## ---- PARAMETERS ---------------------------------------------------------
CM_STAMP <- "cM5"                                 # 5 cM = the saturating DI25 cap (locked)
CLUST    <- sprintf("module_di25/data/di25_clustering_%s.rds", CM_STAMP)
INPUTS   <- "module_di25/data/di25_inputs.rds"
GATE     <- "module_di25/data/di25_sorting_emlg.rds"
SAMPLES  <- "data/hybrids_and_parents_maf005.Rdata"   # sample_data_with_parents (Sample_ID, Population)
RECMAP   <- "data/Frufa_DTOL_PR.ref_genome.recmap"
OUT      <- "module_di25/among_region_association/data/di25D_units.rds"
EMLG_MIN <- 3                                     # >2 markers -> eMLG consensus (matches di25 clustering)

## ---- inputs -------------------------------------------------------------
res <- readRDS(CLUST); g <- res$groups            # 11,052 units: group_id, Chr, representative, n_loci, members
inp <- readRDS(INPUTS)                            # $GTs_hyb (165 x m, 012), $GTs_par (30 x m), $map
gate <- as.data.table(readRDS(GATE))             # 11,052-row gate (superset of moduleC_C3_cl)
e2 <- new.env(); load(SAMPLES, envir = e2)
sd <- as.data.table(e2$sample_data_with_parents)[, .(Sample_ID, Population)]

GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)        # 195 x markers (012, NA = missing), marker colnames
## Recover any di25 individual absent from sample_data (Hei159_h) via its colony prefix:
## population is shared within a colony (<Prefix><colony>_<suffix>), so Hei159_h inherits
## Heinamaki from its Hei159_* nestmates. PC1/PC2 are not needed here, so a label suffices.
colony_of <- function(x) sub("_[^_]*$", "", x)    # Hei159_h -> Hei159
missing   <- setdiff(rownames(GTs_all), sd$Sample_ID)
if (length(missing)) {
  pop_by_colony <- sd[, .(Population = Population[1]), by = .(colony = colony_of(Sample_ID))]
  add <- data.table(Sample_ID = missing,
                    Population = pop_by_colony$Population[match(colony_of(missing), pop_by_colony$colony)])
  stopifnot(!anyNA(add$Population))               # every recovered sample must resolve to a known colony
  message(sprintf("[assemble] recovered %d hybrid(s) absent from sample_data via colony prefix: %s",
                  nrow(add), paste(sprintf("%s->%s", add$Sample_ID, add$Population), collapse = ", ")))
  sd <- rbind(sd, add)
}
pops <- sd$Population[match(rownames(GTs_all), sd$Sample_ID)]
stopifnot(!anyNA(pops))                           # now every di25 individual (incl. Hei159_h) is labelled
is_par   <- grepl("_parent$", pops)
aqu_rows <- which(pops == "aquilonia_parent")
pol_rows <- which(pops == "polyctena_parent")
stopifnot(length(aqu_rows) > 0, length(pol_rows) > 0)
message(sprintf("[assemble] %d individuals (%d hybrids in %d pops + %d parents); %d units",
                nrow(GTs_all), sum(!is_par), length(unique(pops[!is_par])), sum(is_par), nrow(g)))

## ---- one dosage vector per unit (eMLG consensus or representative SNP) ---
## Consensus computed over ALL kept individuals in one pass (single polarization per
## unit, so hybrid & parent sides share a sign convention) -- exactly diem_circos_di25_eMLG.R.
is_emlg <- g$n_loci >= EMLG_MIN
D <- vapply(seq_len(nrow(g)), function(i) {
  if (is_emlg[i]) consensus_dosage(GTs_all, g$members[[i]])   # 0..2 polarized consensus
  else            GTs_all[, g$representative[i]]               # 0/1/2 representative SNP
}, numeric(nrow(GTs_all)))                                     # individuals x units
colnames(D) <- g$group_id

## ---- orient each unit to F. aquilonia via the parent species rows -------
maqu <- colMeans(D[aqu_rows, , drop = FALSE], na.rm = TRUE)
mpol <- colMeans(D[pol_rows, , drop = FALSE], na.rm = TRUE)
flip <- which(maqu < mpol)
D[, flip] <- 2 - D[, flip]                                     # 2 = aquilonia everywhere
message(sprintf("[assemble] %d eMLG-consensus + %d rep-SNP units; %d units flipped to aqu-orientation",
                sum(is_emlg), sum(!is_emlg), length(flip)))

## ---- split: analysis matrix = hybrids; parents kept aside ---------------
U <- D[!is_par, , drop = FALSE]                                # 165 hybrids x 11,052 units, aqu-oriented 0..2
P <- D[is_par,  , drop = FALSE]                                # 30 parents x units (orientation provenance)
pops_hyb <- pops[!is_par]
stopifnot(nrow(U) == length(pops_hyb))

## ---- per-unit chromosome + genetic-map cM (median member cM) ------------
## markers are "Chr<n>:<pos>" strings, so Chr/Pos parse directly; cM by recmap approx.
chr_of <- setNames(as.character(g$Chr), g$group_id)
rec <- fread(RECMAP); setnames(rec, 1:4, c("chr", "pos", "cM", "cMMb"))
rec[, Chr := sub("chromosome_", "Chr", chr)]
memL <- data.table(group_id = rep(g$group_id, lengths(g$members)), marker = unlist(g$members))
memL[, `:=`(Chr = sub(":.*", "", marker), Pos = as.integer(sub(".*:", "", marker)))]
memL[, cM := NA_real_]
for (ch in unique(memL$Chr)) { r <- rec[Chr == ch]; if (nrow(r) < 2) next
  ix <- memL[, which(Chr == ch)]; memL[ix, cM := approx(r$pos, r$cM, xout = Pos, rule = 2)$y] }
cm_of  <- memL[, .(cM = median(cM, na.rm = TRUE)), by = group_id][, setNames(cM, group_id)]
cm_of  <- cm_of[g$group_id]                                    # align to unit order
rate_of <- { rr <- memL[, .(Pos = median(Pos), Chr = Chr[1]), by = group_id]
             rate <- vapply(seq_len(nrow(rr)), function(i) { r <- rec[Chr == rr$Chr[i]]
               if (nrow(r) < 2) NA_real_ else approx(r$pos, r$cMMb, xout = rr$Pos[i], rule = 2)$y }, numeric(1))
             setNames(rate, rr$group_id)[g$group_id] }        # local cM/Mb (recomb annotation, Module-E null)

## ---- per-marker observed heterozygosity over hybrids (paralogy het) -----
GT_hyb_raw <- GTs_all[!is_par, , drop = FALSE]
marker_Ho  <- colMeans(GT_hyb_raw == 1, na.rm = TRUE)          # marker-named; di25D_paralogy het_of uses this

## ---- align the gate to the unit order + sanity ---------------------------
stopifnot(setequal(gate$group_id, g$group_id))
setkey(gate, group_id); gate <- gate[g$group_id]
stopifnot(identical(gate$group_id, g$group_id))

## ---- per-unit HYBRID-BASED metrics (the gate only has parent-based ones) --------
## For studying the loci: within-hybrid MAF (overall + within-population segregation), het,
## expected het, per-locus Fst/Fis, segregating-population count, missingness. Aqu-oriented U.
uid <- colnames(U); is_emlg_u <- is_emlg[match(uid, g$group_id)]
pop_idx <- split(seq_len(nrow(U)), pops_hyb)
pf   <- vapply(pop_idx, function(ix) colMeans(U[ix, , drop = FALSE], na.rm = TRUE) / 2, numeric(ncol(U)))   # units x pops (aqu freq)
hetp <- vapply(pop_idx, function(ix) colMeans(round(U[ix, , drop = FALSE]) == 1, na.rm = TRUE), numeric(ncol(U)))
aqu_freq <- colMeans(U, na.rm = TRUE) / 2
He_pool  <- 2 * aqu_freq * (1 - aqu_freq)
wp_He    <- rowMeans(2 * pf * (1 - pf), na.rm = TRUE)                  # mean within-pop expected het
wp_het   <- rowMeans(hetp, na.rm = TRUE)                               # mean within-pop observed het
metrics <- data.table(
  group_id  = uid, Chr = chr_of[uid], cM = cm_of[uid], rate = rate_of[uid],
  n_loci    = g$n_loci[match(uid, g$group_id)], is_emlg = is_emlg_u,
  aqu_freq  = aqu_freq,
  maf       = pmin(aqu_freq, 1 - aqu_freq),                            # overall within-hybrid MAF
  het_frac  = colMeans(round(U) == 1, na.rm = TRUE),                   # pooled observed het (Ho)
  He_pool   = He_pool,                                                 # pooled expected het
  wp_maf    = rowMeans(pmin(pf, 1 - pf), na.rm = TRUE),                # mean within-population MAF (deme segregation)
  wp_het    = wp_het, wp_He = wp_He,                                   # mean within-population Ho / He
  n_seg_pops = rowSums(pmin(pf, 1 - pf) >= 0.05, na.rm = TRUE),        # # populations segregating (within-pop MAF >= 0.05)
  fst_locus = ifelse(He_pool > 0, (He_pool - wp_He) / He_pool, NA_real_),  # per-locus among-pop differentiation
  fis_locus = ifelse(wp_He > 0, 1 - wp_het / wp_He, NA_real_),         # within-pop HWE deviation (Ho vs He)
  miss_frac = colMeans(is.na(U)))
stopifnot(identical(metrics$group_id, g$group_id))
message(sprintf("[metrics] per-unit hybrid metrics: MAF median %.2f, wp_maf median %.2f, het_frac median %.2f",
                median(metrics$maf), median(metrics$wp_maf), median(metrics$het_frac)))

saveRDS(list(
  dosage    = U,                 # 165 hybrids x 11,052 units, aqu-oriented 0..2 (NA = missing)
  parents   = P,                 # 30 parents x units (orientation provenance)
  pops      = pops_hyb,          # length-165 population labels for the hybrids
  gate      = gate,              # di25_sorting_emlg superset of the moduleC gate, unit-aligned (PARENT-based metrics)
  metrics   = metrics,           # per-unit HYBRID-based metrics (MAF/het/Fst/Fis/segregation) -- unit-aligned
  groups    = g,                 # clustering groups (members, representative, n_loci) -- paralogy/heatmaps
  chr_of    = chr_of, cm_of = cm_of, rate_of = rate_of,   # per-unit Chr / median cM / local cM-Mb
  marker_Ho = marker_Ho,         # per-marker Ho over hybrids (paralogy excess-Ho corroboration)
  params    = list(CM_STAMP = CM_STAMP, EMLG_MIN = EMLG_MIN, clustering = CLUST, gate = GATE)
), OUT)
message(sprintf("[out] %s  |  %d hybrids x %d units  (%d differentiated, %d eMLG-consensus)",
                OUT, nrow(U), ncol(U), sum(gate$differentiated), sum(is_emlg)))
cat(sprintf("\n[done] di25D_units.rds: %d hybrids x %d high-DI units (%s), aqu-oriented; gate + map attached\n",
            nrow(U), ncol(U), CM_STAMP))
