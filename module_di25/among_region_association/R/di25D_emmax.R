## =========================================================================
## DI25 among-region EMMAX arm -- PAIRWISE CANDIDATE NETWORK (structure-corrected).
## =========================================================================
## DESIGN CHANGE (PK, from focal-vs-all to candidate-vs-candidate). Instead of testing a
## few focal loci against ALL ~11k LD units, we define a marginal CANDIDATE SET and test
## only the pairwise associations AMONG the candidates -- a symmetric, network-native scan
## whose nodes are the candidates and whose edges are the structure-corrected pairwise LD.
##
## CANDIDATE SET (marginal, non-circular -- selection is per-locus, independent of the
## pairwise LD being tested): BIDIRECTIONAL loci, bi_score >= BI_MIN = 0.2. bi_score measures
## strong sorting in INCONSISTENT directions across the replicate populations (the equal-
## fitness / DMI-plausible signature). Because a bidirectional locus MUST segregate to sort
## both ways, this criterion also ENFORCES appreciable MAF automatically (min within-hybrid
## MAF 0.16 at bi_score>=0.2) -- so it fixes the low-MAF inflation that the earlier prop_fixed
## gate suffered (prop_fixed anti-correlates with MAF, cor -0.94 with within-pop MAF; it had
## SELECTED near-fixed low-MAF loci, median MAF 0.15, whose within-pop LD was inflated by a
## thin het layer). K = 2083 -> ~2.17M unlinked pairwise tests.
##   * Trade-off (accepted): a candidate criterion tied to marginal sorting is BLIND to
##     ancestry-INDEPENDENT trans-DMIs where BOTH loci look marginally unsorted (Module D's
##     cleanest lead F11431-F49480 was of that class). This arm targets the sorted-but-mixed
##     class; the both-unsorted class would need the (heavier) all-vs-all.
##
## STRUCTURE CORRECTION unchanged and still full-data: K = double-LOCO VanRaden GRM and the
## top-N_PC genome PCs are built from ALL differentiated units (not just the candidates), so
## each pairwise test carries the full ancestry-axis + relatedness correction. Only the
## TESTED traits/partners shrink to the candidate set.
##
## SIGNIFICANCE = pooled (whole-experiment) BH-FDR over all unlinked pairs (~2.17M) at
## Q_FDR, mapped to one r^2 line via emmax's own F = R^2/(1-R^2)*(n-2). This EMMAX FDR is a
## PRE-FILTER (candidate edge list), NOT the DMI verdict -- the discriminant is Module E's
## recombination-matched neutral null; the FDR just bounds false positives vs the no-
## association null within this (sorting-enriched) design.
##
## Run from repo root. Reads di25D_units.rds; sources the local emmax.R + di25D_paralogy.R.
## Writes data/di25D_emmax.rds (pairs [FDR universe], edges [FDR-sig, signed, paralogy-
## flagged], candidates, params) + Figures/di25D_emmax_network_summary.png.

suppressMessages({ library(data.table); library(MASS); library(parallel) })
source("module_di25/among_region_association/R/emmax.R")
source("module_di25/among_region_association/R/di25D_paralogy.R")
set.seed(1)

## ---- PARAMETERS ---------------------------------------------------------
UNITS      <- "module_di25/among_region_association/data/di25D_units.rds"
OUT        <- "module_di25/among_region_association/data/di25D_emmax.rds"
FIG        <- "module_di25/among_region_association/Figures/di25D_emmax_network_summary.png"
BI_MIN     <- 0.2         # candidate gate: bidirectional (bi_score), which also enforces MAF>=0.16
LINK_CM    <- 10          # unlinked = different chr or same chr > LINK_CM cM
Q_FDR      <- 0.05        # pooled Benjamini-Hochberg level over all unlinked pairs
CORES      <- 8
PARALOGY_R <- 0.9         # |within-pop r| above which an edge is flagged as a duplicate
N_PC       <- 10          # top genome PCs to condition each focal trait on (from ALL diff units)
DOUBLE_LOCO <- TRUE       # exclude focal + partner chromosome from K

## =========================================================================
## Setup: units, GRM/PC basis (ALL differentiated), candidate set (prop_fixed gate)
## =========================================================================
u      <- readRDS(UNITS)
eMLG   <- u$dosage; groups <- u$groups; cl <- u$gate
pops_all <- as.character(u$pops); stopifnot(!anyNA(pops_all), length(pops_all) == nrow(eMLG))
chr_of <- u$chr_of; cm_of <- u$cm_of; marker_Ho <- u$marker_Ho
n_ind  <- nrow(eMLG)

## GRM / PC basis = all differentiated units (K must contain the ancestry axis it removes)
scope <- intersect(colnames(eMLG), cl[differentiated == TRUE, group_id])
X <- eMLG[, scope, drop = FALSE]
X <- apply(X, 2, function(v) { v[!is.finite(v)] <- mean(v, na.rm = TRUE); v })
chr_X <- chr_of[scope]
PCS <- if (N_PC > 0) prcomp(X, center = TRUE, scale. = FALSE)$x[, seq_len(N_PC), drop = FALSE] else NULL

## candidate set C = bidirectional units (bi_score >= BI_MIN; also enforces MAF >= 0.16)
C <- intersect(colnames(eMLG), cl[bi_score >= BI_MIN, group_id])
Cchr <- chr_of[C]; Ccm <- cm_of[C]
XC <- eMLG[, C, drop = FALSE]
XC <- apply(XC, 2, function(v) { v[!is.finite(v)] <- mean(v, na.rm = TRUE); v })  # candidate predictors
message(sprintf("[setup] %d individuals | GRM/PC basis = %d diff units%s | %d candidates (bi_score >= %.2f)",
                n_ind, ncol(X), if (N_PC > 0) sprintf(" + top %d PCs", N_PC) else "", length(C), BI_MIN))

## per-chromosome VanRaden contributions -> double-LOCO K
Nc <- list(); dc <- numeric()
for (cc in unique(chr_X)) { M <- X[, chr_X == cc, drop = FALSE]
  p <- colMeans(M) / 2; Z <- sweep(M, 2, 2 * p); Nc[[cc]] <- tcrossprod(Z); dc[cc] <- sum(2 * p * (1 - p)) }
Ntot <- Reduce(`+`, Nc); dtot <- sum(dc)
Kdl <- function(a, b) { N <- Ntot - Nc[[a]]; d <- dtot - dc[[a]]
  if (DOUBLE_LOCO && !is.na(b) && b != a) { N <- N - Nc[[b]]; d <- d - dc[[b]] }; N / d }

## precompute for the cheap uncorrected estimates (independent of K / partner-chromosome
## blocking, so computed ONCE per focal over all candidates, not per emmax call)
pop_idx <- split(seq_len(nrow(XC)), pops_all)                 # individual indices per population
XCbar   <- do.call(rbind, lapply(pop_idx, function(ix) colMeans(XC[ix, , drop = FALSE], na.rm = TRUE)))  # nPop x K candidates

## =========================================================================
## Pairwise scan: each candidate as a trait vs every candidate partner (K depends on the
## partner chromosome under double LOCO, so partners are tested in chromosome blocks).
## Alongside the structure-corrected EMMAX Rsq we return three CHEAP, UNCORRECTED estimates
## so the correction is transparent: r2_raw = plain pairwise r^2 (all individuals unrelated,
## no K/PC); R_st = among-population correlation (cor of the 20 per-pop mean dosages, signed);
## within_pop_r = median signed within-population correlation (structure removed by
## stratification -- a stratified analogue of the GRM correction).
## =========================================================================
scan_one <- function(gid) {
  chA <- chr_of[gid]; y <- eMLG[, gid]; y[!is.finite(y)] <- mean(y, na.rm = TRUE)
  ## cheap uncorrected estimates over ALL candidates (indexed per block below)
  r2raw <- as.vector(suppressWarnings(cor(y, XC)))^2
  ybar  <- vapply(pop_idx, function(ix) mean(y[ix], na.rm = TRUE), numeric(1))
  rst   <- as.vector(suppressWarnings(cor(ybar, XCbar)))
  Wp <- vapply(pop_idx, function(ix) {                        # within-pop r per pop: K x nPop
    yi <- y[ix]; if (length(ix) < 3 || sd(yi) == 0) return(rep(NA_real_, ncol(XC)))
    as.vector(suppressWarnings(cor(yi, XC[ix, , drop = FALSE]))) }, numeric(ncol(XC)))
  wpr <- apply(Wp, 1, median, na.rm = TRUE)                   # median over populations
  dt <- rbindlist(lapply(unique(Cchr), function(chB) {
    cols <- which(Cchr == chB)
    res <- emmax(Y = y, X = XC[, cols, drop = FALSE], K = Kdl(chA, chB), Covar = PCS, B = NULL, cores = 1)
    data.table(focal = gid, partner = C[cols], Chr = chB, cM = Ccm[cols],
               pval = as.numeric(res$pval), Rsq = as.numeric(res$Rsq),
               r2_raw = r2raw[cols], R_st = rst[cols], within_pop_r = wpr[cols])
  }))
  dt[partner != gid]                                          # drop self
}
t0 <- Sys.time()
res <- rbindlist(mclapply(C, scan_one, mc.cores = CORES))
message(sprintf("[scan] %d candidates scanned pairwise | %.0fs", length(C), as.numeric(difftime(Sys.time(), t0, units = "secs"))))

## unlinked = different chr, or same chr > LINK_CM cM apart on the map
res[, unlinked := Chr != chr_of[focal] | (Chr == chr_of[focal] & is.finite(cM) & abs(cM - cm_of[focal]) > LINK_CM)]
## collapse each unordered pair to one row (the test is ~symmetric under the symmetric K)
res[, key := paste(pmin(focal, partner), pmax(focal, partner))]
setorder(res, pval)
pairs <- res[!duplicated(key)]
pairsU <- pairs[unlinked == TRUE]
message(sprintf("[pairs] %s unique pairs (%s unlinked)", format(nrow(pairs), big.mark = ","),
                format(nrow(pairsU), big.mark = ",")))

## =========================================================================
## Pooled BH-FDR over all unlinked pairs -> one r^2 threshold (emmax's own F/df)
## =========================================================================
## genomic inflation over the unlinked bulk. NB the candidate set is SELECTED for strong,
## direction-variable sorting, so lambda > 1 is EXPECTED here (diffuse residual co-ancestry
## among sorted loci; the full focal-vs-all set gives lambda ~= 1.02, so this is enrichment,
## NOT a failure of the K+PC structure model). We report it and do NOT genomic-control the
## FDR (option a): the FDR is only a pre-filter and the discriminant is Module E's neutral
## null, which recalibrates against exactly this residual co-ancestry. lambda>1 does mean the
## uniform-p null is anti-conservative, so the edge list is a lead set, not a DMI count.
#Q_FDR = 0.05
lambda <- median(qchisq(pairsU$pval, df = 1, lower.tail = FALSE)) / qchisq(0.5, 1)
pairsU[, padj := p.adjust(pval, "BH")]
rej    <- pairsU[padj <= Q_FDR]
pstar  <- if (nrow(rej) > 0) max(rej$pval) else NA_real_
df2    <- n_ind - 2                                           # emmax uses df2 = n - 2 (emmax.R:69)
Fc     <- if (is.finite(pstar)) qf(1 - pstar, 1, df2) else NA_real_
r2crit <- if (is.finite(Fc)) Fc / (Fc + df2) else NA_real_
message(sprintf("[FDR] BH q<%.2g over %s unlinked pairs -> %d edges | p* = %.2e, r2crit = %.3f | lambda = %.3f",
                Q_FDR, format(nrow(pairsU), big.mark = ","), nrow(rej), pstar, r2crit, lambda))

## =========================================================================
## Edge list: FDR-significant unlinked pairs -> sign (coupling/repulsion) + paralogy
## =========================================================================
whitened_sign <- function(y, Xsub, K) {
  if (!is.null(PCS)) y <- residuals(lm(y ~ PCS))
  n <- nrow(K); Kn <- (n - 1) / sum((diag(n) - matrix(1, n, n) / n) * K) * K
  nu <- emma.REMLE(y, matrix(1, n), Kn); M <- solve(chol(nu$vg * Kn + nu$ve * diag(n)))
  yt <- crossprod(M, y); xt <- crossprod(M, Xsub); ot <- crossprod(M, rep(1, n))
  sign(lsfit(cbind(ot, xt), yt, intercept = FALSE)$coefficients[2])
}
edges <- copy(rej)
setnames(edges, c("focal", "partner"), c("cluster1", "cluster2"))
if (nrow(edges) > 0) {
  edges[, sign := mapply(function(a, b) whitened_sign(
    { y <- eMLG[, a]; y[!is.finite(y)] <- mean(y, na.rm = TRUE); y },
    matrix({ x <- eMLG[, b]; x[!is.finite(x)] <- mean(x, na.rm = TRUE); x }, ncol = 1),
    Kdl(chr_of[a], chr_of[b])), cluster1, cluster2)]
  edges[, coupling := ifelse(sign > 0, "coupling", "repulsion")]
  het_of <- moduleD_cluster_het(groups, C, marker_Ho)
  edges <- flag_paralogy(edges, "cluster1", "cluster2", eMLG, pops_all, het_of = het_of,
                         thr = PARALOGY_R, cores = CORES)
  ## annotate nodes with gate metadata
  meta <- setNames(cl$DI, cl$group_id); scv <- setNames(cl$sort_class, cl$group_id)
  dirv <- setNames(cl$direction, cl$group_id); pfv <- setNames(cl$prop_fixed, cl$group_id)
  edges[, `:=`(Chr1 = chr_of[cluster1], Chr2 = chr_of[cluster2],
               DI1 = meta[cluster1], DI2 = meta[cluster2], sort1 = scv[cluster1], sort2 = scv[cluster2],
               dir1 = dirv[cluster1], dir2 = dirv[cluster2], pf1 = pfv[cluster1], pf2 = pfv[cluster2])]
  message(sprintf("[edges] %d FDR edges | %d coupling / %d repulsion | %d paralog-flagged; %d clean",
                  nrow(edges), sum(edges$coupling == "coupling"), sum(edges$coupling == "repulsion"),
                  sum(edges$paralog), sum(!edges$paralog)))
}

## =========================================================================
## Outputs + compact network summary figure
## =========================================================================
candidates <- cl[group_id %in% C, .(group_id, Chr = chr_of[group_id], cM = cm_of[group_id],
                                     DI, sort_class, direction, prop_fixed, uni_score, bi_score)]
saveRDS(list(pairs = pairs, edges = edges, candidates = candidates,
             params = list(BI_MIN = BI_MIN, LINK_CM = LINK_CM, Q_FDR = Q_FDR,
                           PARALOGY_R = PARALOGY_R, N_PC = N_PC, n_ind = n_ind, scope_n = ncol(X),
                           pstar = pstar, r2crit = r2crit, lambda = lambda, units = UNITS)),
        OUT)
message("[out] saved ", OUT)

dir.create(dirname(FIG), showWarnings = FALSE)
png(FIG, width = 2000, height = 800, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.2, 0.6, 0))
## (a) Rsq distribution of unlinked pairs with the pooled-FDR threshold
h <- hist(pmin(pairsU$Rsq, 0.6), breaks = 80, plot = FALSE)
plot(h$mids, h$counts + 1, type = "h", log = "y", col = "grey60", lwd = 2,
     xlab = expression(r^2~"(unlinked candidate pairs)"), ylab = "count (log)",
     main = sprintf("a  %s unlinked pairs; FDR q<%.2g", format(nrow(pairsU), big.mark = ","), Q_FDR))
if (is.finite(r2crit)) { abline(v = r2crit, lty = 2, col = "red")
  text(r2crit, max(h$counts), sprintf("r2crit=%.2f", r2crit), pos = 4, col = "red", cex = 0.8) }
## (b) degree distribution of the FDR network (clean edges)
if (nrow(edges) > 0 && sum(!edges$paralog) > 0) {
  ce <- edges[paralog == FALSE]; deg <- table(c(ce$cluster1, ce$cluster2))
  barplot(table(factor(deg, levels = 1:max(deg))), col = "grey70",
          xlab = "node degree (clean FDR edges)", ylab = "# candidate nodes",
          main = sprintf("b  %d edges, %d nodes", nrow(ce), length(deg)))
} else { plot.new(); title("b  no clean FDR edges") }
dev.off()
cat(sprintf("\n[done] %s + %s | %d candidates, %s unlinked pairs, %d FDR edges (%d clean) at q<%.2g (r2crit=%.3f)\n",
            OUT, FIG, length(C), format(nrow(pairsU), big.mark = ","), nrow(rej),
            if (nrow(edges) > 0) sum(!edges$paralog) else 0L, Q_FDR, r2crit))
