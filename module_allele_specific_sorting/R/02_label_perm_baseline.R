## =========================================================================
## module_allele_specific_sorting -- 02: population-label permutation baseline
##
## The realised fraction of the among-population LD ceiling (sum c^2 G_i G_j /
## sum G_i G_j) and mean concordance are compared with a baseline in which the
## loci share NO population partition: for every unit independently, the 20
## population labels are permuted. This keeps each unit's own profile and
## differentiation (G_i, hence the ceiling, is unchanged) and removes every
## correspondence between units. For two unrelated profiles E[c^2] = 1/19, but
## the permutation gives the baseline, with its interval, for exactly the pairs,
## missingness and G-weights of each distance bin.
##
## Within-chromosome bins: explicit pairs (9.6 million), binned by bp and cM.
## Unlinked (207 million cross-chromosome pairs): exact algebraic shortcut. With
## unit-norm standardised profiles z_i (c_ij = z_i . z_j) and per-chromosome
##   M_a = sum_{i in a} g_i z_i z_i'   and   s_a = sum_{i in a} z_i,
##   sum_{a<b} sum_{i in a, j in b} c_ij^2 g_i g_j = (||sum_a M_a||_F^2 - sum_a ||M_a||_F^2) / 2
##   sum_{a<b} sum_{i in a, j in b} c_ij           = (||sum_a s_a||^2   - sum_a ||s_a||^2) / 2
## The observed (unpermuted) values are asserted to reproduce 01_agg_*.rds.
##
## Inputs : module_population_partitioning/data (via 00_utils.R), data/01_units.rds,
##          data/01_agg_within.rds, data/01_agg_cross.rds
## Output : data/02_baseline.rds
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/02_label_perm_baseline.R [N_PERM] [N_CORES]
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
args <- commandArgs(trailingOnly = TRUE)
N_PERM  <- if (length(args) >= 1) as.integer(args[1]) else 500L
N_CORES <- if (length(args) >= 2) as.integer(args[2]) else 1L
SEED <- 1L

U  <- load_units()
u  <- readRDS(file.path(OUT_DATA, "01_units.rds"))
stopifnot(identical(u$group_id, U$u$group_id))
aw <- readRDS(file.path(OUT_DATA, "01_agg_within.rds"))
ac <- readRDS(file.path(OUT_DATA, "01_agg_cross.rds"))

Za <- among_pop_Z(U$Fmat); Zr <- among_pop_Z(U$Resid)
npop <- nrow(Za)
g <- u$G_nei
chrs <- unique(u$Chr)
cols <- split(u$idx, factor(u$Chr, levels = chrs))

## ---- fixed pair structure per chromosome (bins, ceiling weights) -------------
pairs <- lapply(cols, function(ix) {
  m <- length(ix); ut <- which(upper.tri(diag(m)))
  r <- row(diag(m))[ut]; cc <- col(diag(m))[ut]
  i <- ix[r]; j <- ix[cc]
  list(ut = ut, m = m,
       bp = cut(abs(u$Pos[j] - u$Pos[i]), BP_BREAKS, labels = FALSE, right = FALSE),
       cm = cut(abs(u$cM_pos[j] - u$cM_pos[i]), CM_BREAKS, labels = FALSE, right = FALSE),
       gg = g[i] * g[j])
})

## sums for one (possibly permuted) pair of standardised matrices
stats_for <- function(A, R) {
  w <- lapply(seq_along(cols), function(k) {
    p <- pairs[[k]]
    ca <- crossprod(A[, cols[[k]], drop = FALSE])[p$ut]
    cr <- crossprod(R[, cols[[k]], drop = FALSE])[p$ut]
    num <- ca^2 * p$gg
    v <- cbind(num = ifelse(is.na(num), 0, num), den = ifelse(is.na(num), 0, p$gg),
               sc = ifelse(is.na(ca), 0, ca), nc = !is.na(ca),
               sr = ifelse(is.na(cr), 0, cr), nr = !is.na(cr))
    list(bp = rowsum(v, p$bp), cm = rowsum(v[!is.na(p$cm), , drop = FALSE], p$cm[!is.na(p$cm)]))
  })
  agg <- function(type, labels) {
    x <- Reduce(function(a, b) { k <- union(rownames(a), rownames(b))
                                 out <- matrix(0, length(k), ncol(a), dimnames = list(k, colnames(a)))
                                 out[rownames(a), ] <- out[rownames(a), ] + a
                                 out[rownames(b), ] <- out[rownames(b), ] + b; out },
                lapply(w, `[[`, type))
    x <- x[order(as.integer(rownames(x))), , drop = FALSE]
    data.table(bin = labels[as.integer(rownames(x))], realised = x[, "num"] / x[, "den"],
               conc = x[, "sc"] / x[, "nc"], conc_resid = x[, "sr"] / x[, "nr"],
               num = x[, "num"], sc = x[, "sc"], sr = x[, "sr"], den = x[, "den"], nc = x[, "nc"], nr = x[, "nr"])
  }
  ## unlinked: algebraic shortcut over chromosomes
  A0 <- A; A0[is.na(A0)] <- 0; R0 <- R; R0[is.na(R0)] <- 0
  okA <- colSums(is.na(A)) == 0; okR <- colSums(is.na(R)) == 0
  gw <- ifelse(is.na(g) | !okA, 0, g)
  M <- lapply(cols, function(ix) A0[, ix, drop = FALSE] %*% (gw[ix] * t(A0[, ix, drop = FALSE])))
  sa <- lapply(cols, function(ix) rowSums(A0[, ix, drop = FALSE]))
  sr <- lapply(cols, function(ix) rowSums(R0[, ix, drop = FALSE]))
  gsum <- sapply(cols, function(ix) sum(gw[ix])); na <- sapply(cols, function(ix) sum(okA[ix])); nr <- sapply(cols, function(ix) sum(okR[ix]))
  Msum <- Reduce(`+`, M); Ssa <- Reduce(`+`, sa); Ssr <- Reduce(`+`, sr)
  half <- function(tot, parts) (tot - parts) / 2
  num_u <- half(sum(Msum^2), sum(sapply(M, function(x) sum(x^2))))
  den_u <- half(sum(gsum)^2, sum(gsum^2))
  sc_u  <- half(sum(Ssa^2), sum(sapply(sa, function(x) sum(x^2))))
  nc_u  <- half(sum(na)^2, sum(na^2))
  sr_u  <- half(sum(Ssr^2), sum(sapply(sr, function(x) sum(x^2))))
  nr_u  <- half(sum(nr)^2, sum(nr^2))
  unl <- data.table(bin = "unlinked", realised = num_u / den_u, conc = sc_u / nc_u, conc_resid = sr_u / nr_u,
                    num = num_u, sc = sc_u, sr = sr_u, den = den_u, nc = nc_u, nr = nr_u)
  list(bp = rbind(agg("bp", BP_LABELS), unl), cm = rbind(agg("cm", CM_LABELS), unl))
}

## ---- observed, validated against 01 -------------------------------------------
obs <- stats_for(Za, Zr)
chk <- function(o, a, col) {
  ref <- a[, .(num = sum(s_rST2), den = sum(s_ceil2), sc = sum(s_conc), nc = sum(n_conc),
               sr = sum(s_conc_resid), nr = sum(n_conc_resid)), by = c(bin = col)]
  m <- merge(o, ref, by = "bin", suffixes = c("", "_01"))
  for (v in c("num", "den", "sc", "nc", "sr", "nr")) {
    rel <- abs(m[[v]] - m[[paste0(v, "_01")]]) / pmax(abs(m[[paste0(v, "_01")]]), 1e-12)
    if (!all(rel < 1e-6)) stop(sprintf("observed %s does not reproduce 01", v))
  }
}
chk(obs$bp[bin != "unlinked"], copy(aw$bp)[, bin := as.character(bp_bin)], "bin")
chk(obs$cm[bin != "unlinked"], copy(aw$cm)[, bin := as.character(cm_bin)], "bin")
cross_ref <- ac[, .(bin = "unlinked", num = sum(s_rST2), den = sum(s_ceil2), sc = sum(s_conc), nc = sum(n_conc),
                    sr = sum(s_conc_resid), nr = sum(n_conc_resid))]
for (v in c("num", "den", "sc", "nc", "sr", "nr"))
  if (abs(obs$bp[bin == "unlinked"][[v]] - cross_ref[[v]]) / abs(cross_ref[[v]]) >= 1e-6)
    stop(sprintf("unlinked shortcut %s does not reproduce 01", v))
cat("[02] observed sums reproduce 01 exactly (within-chromosome bins and the unlinked shortcut)\n")

## ---- permutations ----------------------------------------------------------------
perm_cols <- function(Z, P) matrix(Z[cbind(as.vector(P), rep(seq_len(ncol(Z)), each = npop))], npop,
                                   dimnames = dimnames(Z))
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
t0 <- Sys.time()
nulls <- parallel::mclapply(seq_len(N_PERM), mc.cores = N_CORES, mc.set.seed = TRUE, FUN = function(b) {
  P <- replicate(ncol(Za), sample.int(npop))            # one permutation per unit, shared by raw and residual
  s <- stats_for(perm_cols(Za, P), perm_cols(Zr, P))
  list(bp = s$bp[, .(bin, realised, conc, conc_resid)], cm = s$cm[, .(bin, realised, conc, conc_resid)])
})
bad <- vapply(nulls, function(x) inherits(x, "try-error") || is.null(x), logical(1))
stopifnot("permutation workers failed" = !any(bad))
cat(sprintf("[02] %d permutations in %.1f min\n", N_PERM, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

summ <- function(type, labels) {
  nl <- rbindlist(lapply(nulls, `[[`, type), idcol = "rep")
  s <- nl[, .(realised_null = mean(realised), realised_lo = quantile(realised, 0.025), realised_hi = quantile(realised, 0.975),
              conc_null = mean(conc), conc_lo = quantile(conc, 0.025), conc_hi = quantile(conc, 0.975),
              conc_resid_null = mean(conc_resid), conc_resid_lo = quantile(conc_resid, 0.025),
              conc_resid_hi = quantile(conc_resid, 0.975)), by = bin]
  o <- merge(obs[[type]][, .(bin, realised_obs = realised, conc_obs = conc, conc_resid_obs = conc_resid)], s, by = "bin")
  o[, realised_excess := realised_obs - realised_null]
  o[, bin := factor(bin, levels = c(labels, "unlinked"))]
  setorder(o, bin); o
}
res <- list(bp = summ("bp", BP_LABELS), cm = summ("cm", CM_LABELS), n_perm = N_PERM, seed = SEED)
cat("\n[02] realised fraction: observed vs population-label permutation baseline (cM bins):\n")
print(res$cm[, .(bin, realised_obs, realised_null, realised_lo, realised_hi, realised_excess)], digits = 4)
saveRDS(res, file.path(OUT_DATA, "02_baseline.rds"))
cat("[02] done\n")
