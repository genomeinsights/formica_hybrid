## =========================================================================
## DI25 among-region -- build the candidate network read-out from the pairwise EMMAX FDR.
## =========================================================================
## The pooled FDR + paralogy filter already ran in di25D_emmax.R (pairwise-candidate scan),
## so this takes the CLEAN FDR edges and applies the remaining honest read-out filters, in
## Module D's order and with the same constants:
##   (1) third-level MERGE  -- collapse within-chromosome, <=LINK_CM, still-correlated nodes
##       into one meta-node (average linkage on 1-|r|, cut at 1-MERGE_R). A single large
##       low-recomb region can appear as several correlated units; this is display grouping,
##       NOT a change to the clustering.
##   (2) LEVERAGE (WITHIN-population) -- drop an edge whose WITHIN-population correlation does
##       not survive leave-one-population-out (min over dropped pops of |median within-pop r|
##       >= LEV_LOO). *** This REPLACES Module D's among-population R_st leverage. *** The
##       r2_raw / R_st / within_pop_r decomposition (di25D_emmax.R) showed the EMMAX Rsq is
##       ORTHOGONAL to among-pop R_st (cor 0.13) and tracks the WITHIN-population signal (0.36):
##       raw r^2 is the structure, EMMAX strips it. So an among-pop R_st leverage would select
##       FOR structure and discard the within-population DMI candidates (the whole point of the
##       structure correction). The within-pop LOO is the correct single-deme robustness test
##       for these individual-level edges. (The among-pop R_st + its LOO are still REPORTED.)
##   (3) recomb ANNOTATION  -- flag low-recombination meta-nodes (rec_pct < STRUCT_PCT) as
##       'structure'; CARRIED to Module E's local-recomb-matched null, never a filter (low
##       recombination is where coadaptation lives, not evidence of neutrality).
##   (4) cross-chromosome MODULES (|consensus r| > MODULE_R, average linkage) -- display.
##
## The surviving meta-edges are the DMI/coadaptation LEADS carried to Module E; the FDR here
## is a pre-filter (lambda = 1.34 in the enriched set), so E's neutral null makes the call.
##
## Reads di25D_emmax.rds (edges, candidates, params) + di25D_units.rds (dosages, pops, map).
## Writes data/di25D_network.rds (meta_nodes, meta_edges, params).

suppressPackageStartupMessages({ library(data.table) })
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a
EMMAX  <- "module_di25/among_region_association/data/di25D_emmax.rds"
UNITS  <- "module_di25/among_region_association/data/di25D_units.rds"
OUT    <- "module_di25/among_region_association/data/di25D_network.rds"
LINK_CM <- 10; MERGE_R <- 0.5; MERGE_METHOD <- "average"; MODULE_R <- 0.4
LEV_LOO <- 0.3; STRUCT_PCT <- 0.10

E <- readRDS(EMMAX); u <- readRDS(UNITS)
eMLG <- u$dosage; pops_all <- as.character(u$pops); groups <- u$groups; gate <- u$gate
chr_of <- u$chr_of; cmv <- u$cm_of; rate_of <- u$rate_of

## per-cluster physical extent (from member markers "Chr:pos") for meta-node spans
memL <- data.table(group_id = rep(groups$group_id, lengths(groups$members)),
                   pos = as.integer(sub(".*:", "", unlist(groups$members))))
cpos <- memL[, .(bp = median(pos), start = min(pos), end = max(pos)), by = group_id]
setkey(cpos, group_id)
bg_rate <- rate_of[gate[differentiated == TRUE, group_id]]
rec_pct <- function(id) mean(bg_rate < rate_of[id], na.rm = TRUE)
maf_of  <- function(id) { p <- mean(eMLG[, id], na.rm = TRUE) / 2; min(p, 1 - p) }
cr      <- function(a, b) abs(cor(eMLG[, a], eMLG[, b], use = "pairwise.complete.obs"))

## ---- start from the clean FDR edges (FDR + paralogy done in di25D_emmax.R) ----
d <- E$edges[paralog == FALSE, .(focal = cluster1, partner = cluster2, r2 = Rsq, pval, coupling)]
d[, r := mapply(function(a, b) cor(eMLG[, a], eMLG[, b], use = "pairwise.complete.obs"), focal, partner)]
nodes <- unique(c(d$focal, d$partner))
cat(sprintf("[in] %d clean FDR edges over %d nodes (r2crit=%.3f, lambda=%.2f)\n",
            nrow(d), length(nodes), E$params$r2crit, E$params$lambda %||% NA_real_))

## ---- (1) within-chromosome average-linkage merge -> meta-nodes -----------
comp <- setNames(integer(length(nodes)), nodes); nid <- 0L
for (ch in unique(chr_of[nodes])) {
  nc <- nodes[chr_of[nodes] == ch]
  if (length(nc) == 1L) { nid <- nid + 1L; comp[nc] <- nid; next }
  D <- matrix(10, length(nc), length(nc), dimnames = list(nc, nc))     # 10 = un-mergeable (>LINK_CM)
  for (i in seq_len(length(nc) - 1L)) for (j in (i + 1L):length(nc))
    D[i, j] <- D[j, i] <- if (abs(cmv[nc[i]] - cmv[nc[j]]) <= LINK_CM) 1 - cr(nc[i], nc[j]) else 10
  diag(D) <- 0
  ct <- cutree(hclust(as.dist(D), method = MERGE_METHOD), h = 1 - MERGE_R)
  comp[nc] <- nid + ct; nid <- nid + max(ct)
}
odeg   <- table(c(d$focal, d$partner))
rep_of <- tapply(names(comp), comp, function(mem) mem[order(-as.integer(odeg[mem]), cmv[mem])][1])
meta_of <- setNames(rep_of[as.character(comp)], names(comp))
cat(sprintf("[merge:%s] %d nodes -> %d meta-nodes\n", MERGE_METHOD, length(nodes), length(unique(meta_of))))

## ---- collapse edges onto meta-nodes --------------------------------------
d[, `:=`(ma = meta_of[focal], mb = meta_of[partner])]
d[, mkey := ifelse(ma < mb, paste(ma, mb), paste(mb, ma))]
meta_edges <- d[ma != mb][, .(coupling = names(which.max(table(coupling))), r2 = max(r2),
  pval = min(pval), n_pairs = .N, pairs = paste(sprintf("%s-%s", focal, partner), collapse = ";")), by = mkey]
meta_edges[, c("ma", "mb") := tstrsplit(mkey, " ", fixed = TRUE)]

## ---- (2) leverage: WITHIN-population leave-one-pop-out |r| ----------------
pop_idx <- split(seq_len(nrow(eMLG)), pops_all)
popm <- function(id) vapply(pop_idx, function(ix) mean(eMLG[ix, id], na.rm = TRUE), numeric(1))
## among-population (REPORTED for context, not the filter)
rst_of  <- function(a, b) cor(popm(a), popm(b), use = "complete.obs")
rst_loo <- function(a, b) { pa <- popm(a); pb <- popm(b); ok <- is.finite(pa) & is.finite(pb)
  pa <- pa[ok]; pb <- pb[ok]; if (length(pa) < 4) return(NA_real_)
  min(abs(sapply(seq_along(pa), function(k) cor(pa[-k], pb[-k])))) }
## per-population signed within-pop correlation (NA where a pop is near-fixed at either locus)
wp_vec <- function(a, b) vapply(pop_idx, function(ix) {
  xa <- eMLG[ix, a]; xb <- eMLG[ix, b]; ok <- is.finite(xa) & is.finite(xb)
  if (sum(ok) < 3 || sd(xa[ok]) == 0 || sd(xb[ok]) == 0) return(NA_real_)
  cor(xa[ok], xb[ok]) }, numeric(1))
## within_pop_r = median over informative pops; wpr_loo = min over dropped pops of |median of rest|
wpr_stat <- function(a, b) { rp <- wp_vec(a, b); rp <- rp[is.finite(rp)]; n <- length(rp)
  if (n < 2) return(c(within_pop_r = if (n == 1) rp else NA_real_, wpr_loo = NA_real_, n_pop = n))
  c(within_pop_r = median(rp), wpr_loo = min(abs(vapply(seq_len(n), function(k) median(rp[-k]), numeric(1)))), n_pop = n) }
ws <- t(mapply(wpr_stat, meta_edges$ma, meta_edges$mb))
meta_edges[, `:=`(rst = mapply(rst_of, ma, mb), rst_loo = mapply(rst_loo, ma, mb),
                  within_pop_r = ws[, "within_pop_r"], wpr_loo = ws[, "wpr_loo"], n_pop = as.integer(ws[, "n_pop"]))]
n_lev <- meta_edges[!(wpr_loo >= LEV_LOO), .N]
kept  <- meta_edges[wpr_loo >= LEV_LOO]
cat(sprintf("[leverage:within-pop] dropped %d/%d edges with leave-one-pop-out |within-pop median r| < %.2f -> %d surviving\n",
            n_lev, nrow(meta_edges), LEV_LOO, nrow(kept)))

## ---- meta-node table + recomb annotation ---------------------------------
memb_of <- tapply(names(meta_of), meta_of, c)
memb_of <- memb_of[names(memb_of) %in% c(kept$ma, kept$mb)]
sc_of  <- setNames(as.character(gate$sort_class), gate$group_id)
dir_of <- setNames(as.character(gate$direction), gate$group_id)
DI_of  <- setNames(gate$DI, gate$group_id)
sort_cat <- function(ids) { t <- table(sc_of[ids]); t <- t[!names(t) %in% c("", "NA", "unsorted")]
  if (length(t) == 0) "unsorted" else names(which.max(t)) }
mdeg <- table(c(kept$ma, kept$mb))
meta_nodes <- if (length(memb_of) == 0) data.table() else rbindlist(lapply(names(memb_of), function(rp) {
  ids <- memb_of[[rp]]; ext <- cpos[ids]
  data.table(meta = rp, chr = chr_of[rp], cM = cmv[rp],
             start = min(ext$start, na.rm = TRUE), end = max(ext$end, na.rm = TRUE),
             n_clusters = length(ids), members = paste(ids, collapse = ";"),
             degree = as.integer(mdeg[rp] %||% 0L), sort_class = sort_cat(ids), direction = dir_of[rp],
             DI = DI_of[rp], MAF = min(sapply(ids, maf_of)), rec_pct = min(sapply(ids, rec_pct))) }))
if (nrow(meta_nodes)) meta_nodes[, structure := rec_pct < STRUCT_PCT]

## ---- (4) cross-chromosome modules ----------------------------------------
if (nrow(meta_nodes) >= 2) {
  mcons <- sapply(meta_nodes$meta, function(mt) {
    ids <- strsplit(meta_nodes[meta == mt, members], ";")[[1]]; rowMeans(eMLG[, ids, drop = FALSE], na.rm = TRUE) })
  mcons <- apply(mcons, 2, function(v) { v[!is.finite(v)] <- mean(v, na.rm = TRUE); v })
  hm <- hclust(as.dist(1 - abs(cor(mcons))), method = "average")
  meta_nodes[, module := cutree(hm, h = 1 - MODULE_R)[meta]]
} else if (nrow(meta_nodes)) meta_nodes[, module := 1L]

setorder(kept, pval)
saveRDS(list(meta_nodes = meta_nodes, meta_edges = kept, all_edges = meta_edges,
             params = list(LINK_CM = LINK_CM, MERGE_R = MERGE_R, MERGE_METHOD = MERGE_METHOD,
                           MODULE_R = MODULE_R, LEV_LOO = LEV_LOO, STRUCT_PCT = STRUCT_PCT,
                           r2crit = E$params$r2crit, lambda = E$params$lambda %||% NA_real_,
                           n_candidates = nrow(E$candidates))),
        OUT)
cat(sprintf("[out] %s | %d surviving meta-edges over %d meta-nodes (%d structure-labelled)\n",
            OUT, nrow(kept), nrow(meta_nodes), if (nrow(meta_nodes)) sum(meta_nodes$structure) else 0L))

## ---- the candidate DMI lead table (what goes to Module E) ----------------
if (nrow(kept)) {
  lead <- copy(kept)
  lead[, `:=`(Chr1 = chr_of[ma], Chr2 = chr_of[mb], dir1 = dir_of[ma], dir2 = dir_of[mb],
              sort1 = sc_of[ma], sort2 = sc_of[mb], DI1 = round(DI_of[ma]), DI2 = round(DI_of[mb]),
              recpct1 = round(sapply(ma, rec_pct), 2), recpct2 = round(sapply(mb, rec_pct), 2))]
  cat("\n[leads] candidate pairwise DMI/coadaptation leads carried to Module E (within-pop leverage):\n")
  print(lead[, .(ma, mb, Chr1, Chr2, coupling, Rsq = round(r2, 3),
                 within_pop_r = round(within_pop_r, 2), wpr_loo = round(wpr_loo, 2), n_pop,
                 rst = round(rst, 2), dir1, dir2, DI1, DI2, recpct1, recpct2)])
}
