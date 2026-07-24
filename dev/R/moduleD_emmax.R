## =========================================================
## MODULE D (EMMAX arm) -- structure-corrected two-locus LD via a linear mixed model.
## =========================================================
## The Ohta screen (moduleD_ohta_dmi.R) showed the raw among-population LD between
## unlinked clusters is dominated by the shared ancestry axis / drift (D2st >> D2is,
## D'2st ~ 0). This script removes that confound the principled way: treat one eMLG's
## consensus dosage as a quantitative trait and test its association against every
## other eMLG in a mixed model with a genomic-relationship random effect K (EMMAX;
## Li, Kemppainen, Rastas & Merila -- the method is built for LD-clustered data).
##
## WHY THIS SUBSUMES THE ANCESTRY-AXIS AND WAHLUND CONTROLS. K's leading eigenvectors
## are the ancestry/structure PCs, so conditioning on K removes the global ancestry
## gradient (the confound that swamped the Ohta screen); and K also carries the
## nestmate/sibling relatedness (9 of 20 "populations" are a single colony), so the
## family-structure Wahlund inflation of within-population LD is regressed out too.
## A residual peak between two loci -- beyond K -- is pair-specific association with
## both the ancestry gradient and the colony structure already removed. Continuous
## eMLGs are used directly (no hard-calling, unlike the Ohta gametic tabulation).
##
## READING A MANHATTAN (one focal eMLG vs the genome):
##   - the self-peak at the focal locus, and its local LD block, are the built-in
##     positive control that the pipeline works; they are expected and ignored.
##   - an ADDITIONAL peak that is UNLINKED from the focal locus (different chromosome,
##     or same chromosome > LINK_CM cM) is the signal: pair-specific coupling not
##     explained by structure -- a candidate epistasis / DMI edge. Its SIGN
##     (coupling vs repulsion) distinguishes cis- from trans-compatible combinations.
##   - PATTERN MAP: unidirectional/fully-sorted loci are near-constant traits with no
##     power (MAF-LD bound; owned by Module A); BIDIRECTIONAL loci sort orthogonally
##     to the ancestry axis, so K does NOT absorb them and the LMM is MAXIMALLY
##     powered exactly there (the symmetric-DMI case); within-deme coupling (Ohta_D>0)
##     becomes a peak with the Wahlund confound already removed.
##
## SIGNIFICANCE. emmax()'s B option draws Ynew ~ MVN(0, vg*K + ve*I) and builds the
## max-F null -> a genome-wide-corrected p under the STRUCTURE null, per focal locus.
## This gives calibrated significance now; the neutral SIM still validates that the
## LMM null is adequate for admixture LD and licenses any DMI aligned WITH the main
## structure axis (which K, by design, absorbs).
##
## SCOPE (settled with PK): TARGETED focal set first (proof of concept), LOCO GRM.
## Run from repo root. Reads canonical clustering, hybrids_only, moduleC_C3_cl.rds,
## recmap, moduleD_ohta.rds; sources PK's dev/emmax.R (LDscnR). Writes
## data/moduleD_emmax.rds and Figures/moduleD_emmax_manhattans.pdf.

suppressMessages({ library(data.table); library(MASS); library(parallel) })
source("dev/R/emmax.R")        # PK's emmax()/emma.REMLE() (Li, Kemppainen, Rastas & Merila); copied verbatim
source("dev/R/moduleD_paralogy.R")   # cross-chromosome duplication filter (shared with the Ohta arm)
set.seed(1)

## ---- PARAMETERS ---------------------------------------------------------
CLUSTERING <- "data/eMLG_5loci_0025_cM05.rds"
CL_GATE    <- "data/moduleC_C3_cl.rds"
RECMAP     <- "data/Frufa_DTOL_PR.ref_genome.recmap"
OHTA       <- "data/moduleD_ohta.rds"                     # for negative/positive-R_st focal picks
LINK_CM    <- 10          # unlinked = different chr or same chr > LINK_CM (matches Ohta arm)
B_PERM     <- 0           # permutations for the max-F structure null; 0 = F-test p only
                          # (peak-calling then uses a Bonferroni genome-wide threshold)
CORES      <- 8
SIG_CORR   <- 0.05        # corrected-p threshold when B_PERM > 0 (else Bonferroni 0.05/m is used)
PARALOGY_R <- 0.9         # |within-pop r| above which an unlinked hit is flagged as a duplicate
## focal-set sizes (the curated Pattern-2 candidates + controls)
N_BI <- 20    # most bidirectional clusters (bi_score) -- prime symmetric-DMI candidates
N_NEG <- 15   # clusters in the most NEGATIVE-R_st Ohta pairs (trans-compatible candidates)
N_POS <- 8    # clusters in the most POSITIVE-R_st Ohta pairs (should be ABSORBED by K -> control)
N_UNI <- 4    # strongly unidirectional clusters (expect self-peak only -> control)
N_RAND <- 5   # random differentiated clusters (expect flat Manhattan -> control)

## =========================================================================
## Setup: eMLG dosage matrix, per-cluster chr + cM, gate metadata
## =========================================================================
clust  <- readRDS(CLUSTERING); eMLG <- clust$eMLG; groups <- clust$groups
cl     <- readRDS(CL_GATE)
e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1)
sample_data <- as.data.table(e1$sample_data)
pops_all <- sample_data[match(rownames(eMLG), Sample_ID), Population]
stopifnot(!anyNA(pops_all))

chr_of <- setNames(as.character(groups$Chr), groups$group_id)
## per-cluster cM (median member cM), interpolated from the recmap (as in the Ohta arm)
map_dt <- as.data.table(e1$map_hyb_005)[, .(marker, Chr = as.character(Chr), Pos)]
rec <- fread(RECMAP); setnames(rec, 1:4, c("chr", "pos", "cM", "cMMb")); rec[, Chr := sub("chromosome_", "Chr", chr)]
map_dt[, cM := NA_real_]
for (ch in unique(map_dt$Chr)) { r <- rec[Chr == ch]; if (nrow(r) < 2) next
  ix <- map_dt[, which(Chr == ch)]; map_dt[ix, cM := approx(r$pos, r$cM, xout = Pos, rule = 2)$y] }
marker_cM <- setNames(map_dt$cM, map_dt$marker)
memL <- groups[.(colnames(eMLG)), on = "group_id", .(marker = unlist(members)), by = group_id]
memL[, cM := marker_cM[marker]]
cm_of <- memL[, .(cM = median(cM, na.rm = TRUE)), by = group_id][, setNames(cM, group_id)]

## predictor set X = differentiated eMLGs (informative markers), NA-imputed to column mean
scope <- intersect(colnames(eMLG), cl[differentiated == TRUE, group_id])
X <- eMLG[, scope, drop = FALSE]
X <- apply(X, 2, function(v) { v[!is.finite(v)] <- mean(v, na.rm = TRUE); v })  # mean-impute
chr_X <- chr_of[scope]; cm_X <- cm_of[scope]
DI_of <- setNames(cl$DI, cl$group_id); sc_of <- setNames(cl$sort_class, cl$group_id)
## per-cluster observed heterozygosity (raw member genotypes) for the paralogy filter
marker_Ho <- colMeans(e1$GTs_hybrids_005 == 1, na.rm = TRUE)
het_of <- moduleD_cluster_het(groups, scope, marker_Ho)
message(sprintf("[setup] X = %d individuals x %d differentiated eMLGs", nrow(X), ncol(X)))

## =========================================================================
## LOCO genomic-relationship matrices (one per chromosome, VanRaden GRM on eMLGs
## EXCLUDING that chromosome -- so a focal locus is never in its own random effect).
## emmax() re-normalises K internally, so absolute scaling is irrelevant.
## =========================================================================
make_grm <- function(M) {                    # M = individuals x markers (imputed)
  p <- colMeans(M) / 2; Z <- sweep(M, 2, 2 * p)
  tcrossprod(Z) / sum(2 * p * (1 - p))
}
chroms <- unique(chr_X)
K_loco <- setNames(lapply(chroms, function(ch) make_grm(X[, chr_X != ch, drop = FALSE])), chroms)
message(sprintf("[setup] built %d LOCO GRMs (n = %d)", length(K_loco), nrow(X)))

## =========================================================================
## Focal set: curated Pattern-2 candidates + controls
## =========================================================================
cl_sc <- cl[group_id %in% scope]
cl_sc[, bi_score := prop_fixed - abs(uni_score)]         # bidirectional component (Module A decomposition)
foc_bi  <- cl_sc[order(-bi_score)][1:N_BI, group_id]
foc_uni <- cl_sc[order(-abs(uni_score))][1:N_UNI, group_id]
foc_rand <- sample(setdiff(scope, c(foc_bi, foc_uni)), N_RAND)
## Ohta negative / positive-R_st partners (unique clusters from the extreme pairs)
oh <- readRDS(OHTA)$pairs
neg_ids <- unique(unlist(oh[order(R_st)][1:(3 * N_NEG), .(cluster1, cluster2)]))
foc_neg <- intersect(neg_ids, scope)[1:N_NEG]
pos_ids <- unique(unlist(oh[order(-R_st)][1:(3 * N_POS), .(cluster1, cluster2)]))
foc_pos <- intersect(pos_ids, scope)[1:N_POS]
focal <- data.table(group_id = unique(c(foc_bi, foc_neg, foc_pos, foc_uni, foc_rand)))
focal[, set := fifelse(group_id %in% foc_bi, "bidirectional",
                fifelse(group_id %in% foc_neg, "neg_Rst",
                fifelse(group_id %in% foc_pos, "pos_Rst(ctrl)",
                fifelse(group_id %in% foc_uni, "unidir(ctrl)", "random(ctrl)"))))]
focal[, `:=`(Chr = chr_of[group_id], cM = cm_of[group_id],
             DI = DI_of[group_id], sort_class = sc_of[group_id])]
set_of <- setNames(focal$set, focal$group_id)
message(sprintf("[focal] %d focal loci: %s", nrow(focal),
                paste(names(table(focal$set)), table(focal$set), sep = "=", collapse = ", ")))

## =========================================================================
## Scan: one EMMAX Manhattan per focal locus (LOCO K for its chromosome)
## =========================================================================
## structure-corrected SIGN of an association: sign of the GLS coefficient, i.e. the
## sign of the whitened covariance between focal and partner (recomputed only for hits).
whitened_sign <- function(y, Xsub, K) {
  Kn <- (nrow(K) - 1) / sum((diag(nrow(K)) - matrix(1, nrow(K), nrow(K)) / nrow(K)) * K) * K
  nu <- emma.REMLE(y, matrix(1, length(y)), Kn)
  M <- solve(chol(nu$vg * Kn + nu$ve * diag(nrow(K))))
  yt <- crossprod(M, y); Xt <- crossprod(M, Xsub)
  sign(as.numeric(crossprod(Xt - mean(Xt), yt - mean(yt))))   # per column
}

scan_one <- function(gid) {
  ch <- chr_of[gid]; y <- eMLG[, gid]; y[!is.finite(y)] <- mean(y, na.rm = TRUE)
  res <- emmax(Y = y, X = X, K = K_loco[[ch]], B = if (B_PERM > 0) B_PERM else NULL, cores = 1)
  dt <- data.table(partner = colnames(X), Chr = chr_X, cM = cm_X,
                   F = as.numeric(res$F), pval = as.numeric(res$pval), Rsq = as.numeric(res$Rsq))
  if (!is.null(res$pval.corr)) dt[, pval_corr := as.numeric(res$pval.corr)]
  ## unlinked = different chromosome, or same chromosome > LINK_CM from the focal locus
  dt[, unlinked := Chr != ch | (Chr == ch & is.finite(cM) & abs(cM - cm_of[gid]) > LINK_CM)]
  dt[, focal := gid]
  dt
}
t0 <- Sys.time()
res_list <- mclapply(focal$group_id, scan_one, mc.cores = CORES)     # one Manhattan per focal locus
names(res_list) <- focal$group_id
message(sprintf("[scan] %d focal loci scanned | %.0fs", nrow(focal), as.numeric(difftime(Sys.time(), t0, units = "secs"))))

## significance: max-F permutation p if available, else Bonferroni over the m tested markers
pcol    <- if (B_PERM > 0) "pval_corr" else "pval"
sig_thr <- if (B_PERM > 0) SIG_CORR else 0.05 / ncol(X)

## candidate edge list: significant UNLINKED peaks (self-locus and linked block excluded)
hits <- rbindlist(lapply(res_list, function(dt) dt[unlinked == TRUE & get(pcol) < sig_thr]))
if (nrow(hits) > 0) {
  hits[, sign := mapply(function(f, p) whitened_sign(
    { y <- eMLG[, f]; y[!is.finite(y)] <- mean(y, na.rm = TRUE); y },
    matrix(X[, p], ncol = 1), K_loco[[chr_of[f]]]), focal, partner)]
  hits[, `:=`(focal_set = set_of[focal], coupling = ifelse(sign > 0, "coupling", "repulsion"))]
  ## cross-chromosome paralogy / duplication filter (within-pop |r| ~ 1 = same reads)
  hits <- flag_paralogy(hits, "focal", "partner", eMLG, pops_all, het_of = het_of,
                        thr = PARALOGY_R, cores = CORES)
}
message(sprintf("[scan] %d significant unlinked peaks (%s < %.2g, %s)",
                nrow(hits), pcol, sig_thr, if (B_PERM > 0) "max-F perm" else "Bonferroni"))
if (nrow(hits) > 0)
  message(sprintf("[paralogy] %d/%d hits flagged as duplicates (within-pop |r| > %.2f); %d clean candidates remain",
                  sum(hits$paralog), nrow(hits), PARALOGY_R, sum(!hits$paralog)))

## =========================================================================
## Outputs + Manhattan figure (one page per focal locus)
## =========================================================================
saveRDS(list(results = res_list, hits = hits, focal = focal, het_of = het_of,
             params = list(LINK_CM = LINK_CM, B_PERM = B_PERM, SIG_CORR = SIG_CORR,
                           PARALOGY_R = PARALOGY_R, scope_n = ncol(X), clustering = CLUSTERING)),
        "data/moduleD_emmax.rds")
message("[out] saved data/moduleD_emmax.rds")

dir.create("Figures", showWarnings = FALSE)
chr_order <- unique(chr_X[order(as.integer(sub("Chr", "", chr_X)))])
chr_col <- setNames(rep(c("grey55", "grey30"), length.out = length(chr_order)), chr_order)
pdf("Figures/moduleD_emmax_manhattans.pdf", width = 9, height = 3.2)
for (gid in focal$group_id) {
  dt <- copy(res_list[[gid]]); dt[, Chr := factor(Chr, levels = chr_order)]
  setorder(dt, Chr, cM); dt[, x := .I]; dt[, lp := -log10(pmax(get(pcol), .Machine$double.xmin))]
  thr <- -log10(sig_thr)
  plot(dt$x, dt$lp, pch = 20, cex = 0.35, col = chr_col[as.character(dt$Chr)],
       xlab = "", ylab = expression(-log[10]~p), xaxt = "n",
       main = sprintf("%s  [%s]  %s  DI=%.0f  (self=red, hit=orange, paralog=grey x)",
                      gid, focal[group_id == gid, set], focal[group_id == gid, Chr],
                      focal[group_id == gid, DI]))
  abline(h = thr, lty = 2, col = "grey70")
  points(dt[partner == gid, x], dt[partner == gid, lp], col = "red", pch = 20, cex = 0.9)   # self
  par_partners <- if (nrow(hits) > 0) hits[focal == gid & paralog == TRUE, partner] else character(0)
  hh <- dt[unlinked == TRUE & get(pcol) < sig_thr]
  if (nrow(hh) > 0) {
    hh[, par := partner %in% par_partners]
    points(hh[par == FALSE, x], hh[par == FALSE, lp], col = "#E69F00", pch = 1, cex = 1.3, lwd = 1.5)  # clean hit
    points(hh[par == TRUE,  x], hh[par == TRUE,  lp], col = "grey55", pch = 4, cex = 1.1, lwd = 1.2)   # paralog (x)
  }
  cm <- dt[, .(mid = mean(x)), by = Chr]; axis(1, at = cm$mid, labels = sub("Chr", "", cm$Chr), cex.axis = 0.5, las = 2)
}
dev.off()
cat(sprintf("\n[done] data/moduleD_emmax.rds, Figures/moduleD_emmax_manhattans.pdf | %d hits (%d clean, %d paralog) over %d focal loci\n",
            nrow(hits), if (nrow(hits) > 0) sum(!hits$paralog) else 0L,
            if (nrow(hits) > 0) sum(hits$paralog) else 0L, nrow(focal)))
