## =========================================================================
## Build the Module D trans/cis association network from the EMMAX scan, on ONE
## honest multiple-testing scale, with same-region clusters collapsed.
## =========================================================================
## Three steps, in order:
##  (1) GLOBAL FDR. The EMMAX arm tested each of the 51 focal clusters against all
##      differentiated partners; the original hit-calling applied a Bonferroni
##      threshold *independently per focal Manhattan*, which does not correct for
##      testing 51 focals. Here we pool ALL unlinked (focal, partner) tests onto a
##      single Benjamini-Hochberg scale and keep q < Q_FDR. (q<0.01 reproduces the
##      old hit count ~ 276 but on a whole-experiment footing.)
##  (2) PARALOGY FILTER. Drop cross-locus duplicates (|within-pop r| > PARALOGY_R);
##      technical, not epistasis (see moduleD_paralogy.R).
##  (3) THIRD-LEVEL (context-only) MERGE. The LD-clustering deliberately does not
##      over-merge (so eMLGs stay clean), so a single large low-recombination region
##      can appear as several correlated clusters -- which then reads as "multiple
##      independent interacting loci". For THIS network only, clusters on the same
##      chromosome within LINK_CM that remain correlated (|consensus r| > MERGE_R)
##      are collapsed into one meta-node (average-linkage clustering, MERGE_METHOD; single
##      linkage / connected components chains uncorrelated ends -- do not use). This is NOT a
##      change to the canonical clustering; it is a display/interpretation grouping.
##  Then: a single-population LEVERAGE filter (leave-one-pop-out |r| >= LEV_LOO, which subsumes
##  any near-fixed/MAF cutoff); a low-recombination ANNOTATION (STRUCT_PCT) that is carried to
##  Module E's local-recomb-matched null, NOT used as a filter; and cross-chromosome MODULES
##  (MODULE_R). See manuscript_notes/moduleD_plan.md for the minimal-pipeline rationale.
##
## Emits data/moduleD_network.rds = list(meta_nodes, meta_edges, params), consumed by
## dev/R/moduleD_network_circos.R and moduleD_module_heatmaps.R. Run from the repo root.

suppressPackageStartupMessages({ library(data.table); library(igraph) })
source("moduleD_among_region_association/R/moduleD_paralogy.R")
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

Q_FDR      <- 0.01   # global Benjamini-Hochberg threshold across all unlinked tests
LINK_CM    <- 10     # same-chromosome merge window (matches the unlinked definition)
MERGE_R    <- 0.5    # |consensus r| above which two within-window clusters are one region
MERGE_METHOD <- "average"   # linkage for the third-level merge; "single" = old connected-components (chains)
MODULE_R   <- 0.4    # cross-chromosome correlation MODULES: average-linkage cut at |consensus r| (display grouping)
PARALOGY_R <- 0.9    # |within-pop r| duplicate threshold
LEV_LOO    <- 0.3    # leverage: keep an edge only if its among-pop |r| survives leave-one-pop-out
STRUCT_PCT <- 0.10   # low-recomb 'structure' ANNOTATION (never a filter): recomb percentile below this.
                     # Low recomb is where coadaptation lives, so these are carried to the null (Module E),
                     # judged against a LOCAL-recomb-matched neutral baseline -- not discarded. See moduleD_plan.md.
EMMAX      <- "moduleD_among_region_association/data/moduleD_emmax.rds"
CLUSTERING <- "data/eMLG_5loci_0025_cM05.rds"
CL_GATE    <- "data/moduleC_C3_cl.rds"
RECMAP     <- "data/Frufa_DTOL_PR.ref_genome.recmap"

## ---- inputs --------------------------------------------------------------
E <- readRDS(EMMAX); clust <- readRDS(CLUSTERING); groups <- clust$groups; eMLG <- clust$eMLG
cl <- readRDS(CL_GATE)
chr_of <- setNames(as.character(groups$Chr), groups$group_id)
e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1)
sample_data <- as.data.table(e1$sample_data)
pops_all <- sample_data[match(rownames(eMLG), Sample_ID), Population]
map <- as.data.table(e1$map_hyb_005)[, .(marker, Chr = as.character(Chr), Pos)]
rec <- fread(RECMAP); setnames(rec, 1:4, c("chr","pos","cM","cMMb")); rec[, Chr := sub("chromosome_","Chr",chr)]
map[, `:=`(cM = NA_real_, rate = NA_real_)]
for (ch in unique(map$Chr)) { r <- rec[Chr == ch]; if (nrow(r) < 2) next
  ix <- map[, which(Chr == ch)]
  map[ix, cM   := approx(r$pos, r$cM,   xout = Pos, rule = 2)$y]
  map[ix, rate := approx(r$pos, r$cMMb, xout = Pos, rule = 2)$y] }
memL <- groups[, .(marker = unlist(members)), by = group_id]
memL <- map[memL, on = "marker"]
cpos <- memL[, .(cM = median(cM, na.rm = TRUE), bp = median(Pos, na.rm = TRUE),
                 start = min(Pos, na.rm = TRUE), end = max(Pos, na.rm = TRUE),
                 rate = median(rate, na.rm = TRUE)), by = group_id]
cmv <- setNames(cpos$cM, cpos$group_id)
## recombination percentile of each cluster vs the differentiated-cluster background (structure label)
crate <- setNames(cpos$rate, cpos$group_id); bg_rate <- crate[cl[differentiated == TRUE, group_id]]
rec_pct <- function(id) mean(bg_rate < crate[id], na.rm = TRUE)
## consensus minor-allele frequency of a cluster (for the meta-node MAF column / diagnostics)
maf_of <- function(id) { p <- mean(eMLG[, id], na.rm = TRUE) / 2; min(p, 1 - p) }

## ---- (1) global FDR over all unlinked tests ------------------------------
all <- rbindlist(lapply(names(E$results), function(g)
  E$results[[g]][unlinked == TRUE, .(focal = g, partner, pval)]))
all[, key := ifelse(focal < partner, paste(focal, partner), paste(partner, focal))]
setorder(all, pval); u <- all[!duplicated(key)]            # unique unordered pair, best direction
m_tests <- nrow(u); u[, q := p.adjust(pval, "BH")]
d <- u[q < Q_FDR]
d[, r := mapply(function(a, b) cor(eMLG[, a], eMLG[, b], use = "pairwise.complete.obs"), focal, partner)]
d[, coupling := ifelse(r < 0, "trans", "cis")]
cat(sprintf("[fdr] %s unlinked tests; q<%.g -> %d pairs (trans %d / cis %d)\n",
            format(m_tests, big.mark = ","), Q_FDR, nrow(d), d[coupling=="trans",.N], d[coupling=="cis",.N]))

## ---- (2) paralogy filter -------------------------------------------------
marker_Ho <- colMeans(e1$GTs_hybrids_005 == 1, na.rm = TRUE)
het_of <- moduleD_cluster_het(groups, unique(c(d$focal, d$partner)), marker_Ho)
d <- flag_paralogy(d, "focal", "partner", eMLG, pops_all, het_of = het_of, thr = PARALOGY_R)
cat(sprintf("[paralogy] dropped %d paralogous pairs (|within-pop r| > %.2f)\n", d[paralog == TRUE, .N], PARALOGY_R))
d <- d[paralog == FALSE]
nodes <- unique(c(d$focal, d$partner))

## ---- (3) third-level merge: cluster within-chromosome nodes by |consensus r| ---
## Connected components of the |r|>MERGE_R graph (the old default) is SINGLE linkage: a
## chain A-B-C merges even when A and C are uncorrelated. Instead cluster each chromosome's
## nodes hierarchically on distance 1-|r| with MERGE_METHOD linkage (default "complete",
## which requires the whole merged group to be mutually correlated) and cut at 1-MERGE_R.
## Pairs farther than LINK_CM cM apart are given an un-mergeable distance (1), so a merged
## region is both correlated AND proximal, and cannot chain across a >LINK_CM gap.
cr <- function(a, b) abs(cor(eMLG[, a], eMLG[, b], use = "pairwise.complete.obs"))
comp <- setNames(integer(length(nodes)), nodes); nid <- 0L
for (ch in unique(chr_of[nodes])) {
  nc <- nodes[chr_of[nodes] == ch]
  if (length(nc) == 1L) { nid <- nid + 1L; comp[nc] <- nid; next }
  D <- matrix(10, length(nc), length(nc), dimnames = list(nc, nc))   # 10 = un-mergeable (>LINK_CM)
  for (i in seq_len(length(nc) - 1L)) for (j in (i + 1L):length(nc))
    D[i, j] <- D[j, i] <- if (abs(cmv[nc[i]] - cmv[nc[j]]) <= LINK_CM) 1 - cr(nc[i], nc[j]) else 10
  diag(D) <- 0
  ct <- cutree(hclust(as.dist(D), method = MERGE_METHOD), h = 1 - MERGE_R)
  comp[nc] <- nid + ct; nid <- nid + max(ct)
}
odeg <- table(c(d$focal, d$partner))                       # pre-merge degree -> pick representative
rep_of <- tapply(names(comp), comp, function(mem) mem[order(-as.integer(odeg[mem]), cmv[mem])][1])
meta_of <- setNames(rep_of[as.character(comp)], names(comp))
cat(sprintf("[merge:%s] %d nodes -> %d meta-nodes (%d clusters absorbed into %d multi-cluster meta-nodes)\n",
            MERGE_METHOD, length(nodes), length(unique(meta_of)),
            length(nodes) - length(unique(meta_of)), sum(tabulate(comp) > 1)))

## ---- collapse edges onto meta-nodes --------------------------------------
d[, `:=`(ma = meta_of[focal], mb = meta_of[partner])]
d[, mkey := ifelse(ma < mb, paste(ma, mb), paste(mb, ma))]
meta_edges <- d[ma != mb][, .(
  q = min(q), coupling = names(which.max(table(coupling))), n_pairs = .N,
  focal_pairs = paste(sprintf("%s-%s", focal, partner), collapse = ";")), by = .(mkey)]
meta_edges[, c("ma", "mb") := tstrsplit(mkey, " ", fixed = TRUE)]

## ---- (3) leverage filter (structure loci are LABELLED, not removed -- see below) ---
memb_of <- tapply(names(meta_of), meta_of, c)              # meta rep -> member cluster ids
## single-population leverage: drop edges whose among-pop |r| does not survive leave-one-pop-out.
## This subsumes any near-fixed / low-MAF cutoff (a rare-allele edge that a single deme drives
## fails here), so no separate MIN_MAF filter is needed on the maf005-gated data.
popm <- function(id) tapply(eMLG[, id], pops_all, mean, na.rm = TRUE)
loo_min <- function(a, b) { pa <- popm(a); pb <- popm(b); ok <- is.finite(pa) & is.finite(pb)
  pa <- pa[ok]; pb <- pb[ok]; min(abs(sapply(seq_along(pa), function(k) cor(pa[-k], pb[-k])))) }
meta_edges[, loo_min := mapply(loo_min, ma, mb)]
n_lev <- meta_edges[loo_min < LEV_LOO, .N]; meta_edges <- meta_edges[loo_min >= LEV_LOO]
cat(sprintf("[leverage] dropped %d edges with leave-one-pop-out |r| < %.2f\n", n_lev, LEV_LOO))

## ---- meta-node table (span = union of member-cluster extents; for Step-4 bands) ----
memb_of <- memb_of[names(memb_of) %in% c(meta_edges$ma, meta_edges$mb)]   # surviving meta-nodes
bi_reps <- tryCatch(readRDS("moduleD_among_region_association/data/moduleD_bidirectional.rds")$reps$group_id, error = function(e) character(0))
sc_of   <- setNames(as.character(cl$sort_class), cl$group_id)
sort_cat <- function(ids) {                               # meta sort category for the colour bar
  if (any(ids %in% bi_reps)) return("bidirectional")
  t <- table(sc_of[ids]); t <- t[!names(t) %in% c("", "NA", "unsorted")]
  if (length(t) == 0) "unsorted" else names(which.max(t))
}
mdeg <- table(c(meta_edges$ma, meta_edges$mb))
meta_nodes <- rbindlist(lapply(names(memb_of), function(rep) {
  ids <- memb_of[[rep]]; ext <- cpos[group_id %in% ids]
  data.table(meta = rep, chr = chr_of[rep], cM = cmv[rep],
             start = min(ext$start, na.rm = TRUE), end = max(ext$end, na.rm = TRUE),
             n_clusters = length(ids), members = paste(ids, collapse = ";"),
             degree = as.integer(mdeg[rep] %||% 0L), sort_class = sort_cat(ids),
             MAF = min(sapply(ids, maf_of)), rec_pct = min(sapply(ids, rec_pct))) }))
meta_nodes[, structure := rec_pct < STRUCT_PCT]  # low-recomb ANNOTATION (kept + carried to E, not a filter)
meta_nodes <- meta_nodes[meta %in% c(meta_edges$ma, meta_edges$mb)]   # keep connected meta-nodes only
## cross-chromosome correlation MODULES: average-linkage clustering of the meta-node
## consensuses (across chromosomes), cut at |r| > MODULE_R. Unlike the within-chromosome
## merge, this groups co-varying regions genome-wide, so each module is one coherent
## co-ancestry block (used to give a single-module heatmap and within-module circos links).
mcons <- sapply(meta_nodes$meta, function(mt) {
  ids <- strsplit(meta_nodes[meta == mt, members], ";")[[1]]; rowMeans(eMLG[, ids, drop = FALSE], na.rm = TRUE) })
mcons <- apply(mcons, 2, function(v) { v[!is.finite(v)] <- mean(v, na.rm = TRUE); v })
hm <- hclust(as.dist(1 - abs(cor(mcons))), method = "average")
meta_nodes[, module := cutree(hm, h = 1 - MODULE_R)[meta]]
setorder(meta_nodes, -degree)
cat(sprintf("[modules] %d cross-chromosome modules at |r|>%.2f (largest %d meta-nodes)\n",
            uniqueN(meta_nodes$module), MODULE_R, max(table(meta_nodes$module))))

saveRDS(list(meta_nodes = meta_nodes, meta_edges = meta_edges,
             params = list(Q_FDR = Q_FDR, LINK_CM = LINK_CM, MERGE_R = MERGE_R, MERGE_METHOD = MERGE_METHOD,
                           MODULE_R = MODULE_R, PARALOGY_R = PARALOGY_R, LEV_LOO = LEV_LOO,
                           STRUCT_PCT = STRUCT_PCT, m_tests = m_tests, clustering = CLUSTERING)),
        "moduleD_among_region_association/data/moduleD_network.rds")
cat(sprintf("[out] %d meta-edges over %d meta-nodes (trans %d / cis %d); %d hubs (deg>=10); %d structure-labelled\n",
            nrow(meta_edges), nrow(meta_nodes), meta_edges[coupling=="trans",.N],
            meta_edges[coupling=="cis",.N], meta_nodes[degree >= 10, .N], meta_nodes[structure == TRUE, .N]))
cat("top meta-hubs:\n")
print(head(meta_nodes[, .(meta, chr, cM = round(cM,1), n_clusters, degree, MAF = round(MAF,2), structure, sort_class)], 8))
