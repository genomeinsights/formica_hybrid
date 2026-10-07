## =========================================================================
## module_allele_specific_sorting -- EXPLORATORY (not in the document, not committed):
## can founders built as GENOTYPE MOSAICS of the empirical parents reproduce the
## empirical parental LD, without phasing?
##
## A synthetic founder is built chromosome by chromosome: start with a random donor
## among the 15 empirical parents of the species, copy the donor's genotypes, and
## switch to a new random donor at breakpoints placed along the genetic map as a
## Poisson process with rate LAMBDA switches per cM. Allele frequencies and
## heterozygosity therefore come from the empirical parents directly.
## Trade-off: LAMBDA -> 0 gives exact copies of donors (LD preserved, no new
## diversity); larger LAMBDA mixes donors but breaks LD between markers separated
## by a switch, so mosaic LD is always <= donor LD.
## Evaluation (same as explore_parent_ld.R): within-species r^2 between the 20,807
## unit SNPs in a synthetic sample of n = 15 (same finite-sample floor as the
## empirical parents), by cM bin (zero-map-distance pairs excluded) and unlinked;
## close pairs (< 5 kb) split by same / different simulation (min_r2 = 0.2) cluster;
## diversity = mean pairwise genotype identity among founders.
## Output: data/explore_mosaic_founders.rds, Figures/explore_mosaic_founders.png
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/exploratory/explore_mosaic_founders.R [N_SETS] [N_CORES]
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
args <- commandArgs(trailingOnly = TRUE)
N_SETS  <- if (length(args) >= 1) as.integer(args[1]) else 10L
N_CORES <- if (length(args) >= 2) as.integer(args[2]) else 8L
LAMBDAS <- c(0, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2)
N_OUT   <- 15L

u0 <- readRDS(file.path(OUT_DATA, "01_units.rds"))[, .(Chr, Pos, cM_pos, unit_marker)]
gi <- as.data.table(readRDS("group_info_new.rds"))
gi[, mk := paste0("Chr", sub("chromosome_", "", chromosome), ":", as.integer(position))]
u0[, simcl := gi$group_id[match(unit_marker, gi$mk)]]
chrs <- unique(u0$Chr); cols <- split(seq_len(nrow(u0)), factor(u0$Chr, levels = chrs))
## cM for mosaic breakpoints: units outside the map range have NA cM_pos -> nearest mapped neighbour
cmp_pos <- u0$cM_pos
for (ix in cols) { v <- cmp_pos[ix]; if (anyNA(v)) cmp_pos[ix] <- approx(seq_along(v)[!is.na(v)], v[!is.na(v)], seq_along(v), rule = 2)$y }

pairs <- lapply(cols, function(ix) {
  m <- length(ix); ut <- which(upper.tri(diag(m)))
  i <- ix[row(diag(m))[ut]]; j <- ix[col(diag(m))[ut]]
  dcm <- abs(u0$cM_pos[j] - u0$cM_pos[i])
  cm <- cut(dcm, CM_BREAKS, labels = FALSE, right = FALSE); cm[which(dcm == 0)] <- NA
  close <- abs(u0$Pos[j] - u0$Pos[i]) < 5e3
  list(ut = ut, cm = cm, close = close, same = u0$simcl[i] == u0$simcl[j])
})

ld_eval <- function(G) {
  X <- sweep(G, 2, colMeans(G, na.rm = TRUE)); X[is.na(X)] <- 0
  nrm <- sqrt(colSums(X^2)); poly <- nrm > 0
  Z <- sweep(X, 2, ifelse(poly, nrm, NA), "/")
  w <- rbindlist(lapply(seq_along(cols), function(k) {
    p <- pairs[[k]]; r2 <- crossprod(Z[, cols[[k]], drop = FALSE])[p$ut]^2
    rbind(data.table(bin = p$cm, r2 = r2)[!is.na(bin) & !is.na(r2), .(s = sum(r2), n = .N), by = bin],
          data.table(bin = ifelse(p$same, 101L, 102L)[p$close], r2 = r2[p$close])[!is.na(r2), .(s = sum(r2), n = .N), by = bin])
  }))[, .(s = sum(s), n = sum(n)), by = bin]
  Z0 <- Z; Z0[is.na(Z0)] <- 0
  M <- lapply(cols, function(ix) tcrossprod(Z0[, ix, drop = FALSE])); np <- sapply(cols, function(ix) sum(poly[ix]))
  half <- function(tot, parts) (tot - parts) / 2
  rbind(w, data.table(bin = 99L, s = half(sum(Reduce(`+`, M)^2), sum(sapply(M, function(m) sum(m^2)))),
                      n = half(sum(np)^2, sum(np^2))))[, .(bin, r2 = s / n)]
}
identity_mean <- function(G) {
  n <- nrow(G); v <- numeric(0)
  for (a in 1:(n - 1)) for (b in (a + 1):n) v <- c(v, mean(G[a, ] == G[b, ], na.rm = TRUE))
  mean(v)
}

mosaic <- function(D, lambda, n_out) {
  out <- matrix(NA_integer_, n_out, ncol(D))
  for (f in seq_len(n_out)) for (ix in cols) {
    cmv <- cmp_pos[ix]; len <- max(cmv) - min(cmv)
    nb <- if (lambda > 0) rpois(1, lambda * len) else 0L
    brk <- sort(runif(nb, min(cmv), max(cmv)))
    seg <- findInterval(cmv, brk) + 1L
    donor <- sample.int(nrow(D), length(brk) + 1L, replace = TRUE)
    out[f, ix] <- D[cbind(donor[seg], ix)]
  }
  out
}

inp <- readRDS("module_di25/data/di25_inputs.rds")
P <- inp$GTs_par[, u0$unit_marker]
donors <- list(aquilonia = P[grepl("^Faqu", rownames(P)), ], polyctena = P[grepl("^Fpol", rownames(P)), ])

emp <- rbindlist(lapply(names(donors), function(sp) ld_eval(donors[[sp]])[, `:=`(species = sp, lambda = NA_real_, set = 0L)]))
emp_id <- sapply(donors, identity_mean)

grid <- CJ(species = names(donors), lambda = LAMBDAS, set = seq_len(N_SETS))
RNGkind("L'Ecuyer-CMRG"); set.seed(1)
res <- parallel::mclapply(seq_len(nrow(grid)), mc.cores = N_CORES, mc.set.seed = TRUE, FUN = function(k) {
  g <- grid[k]; G <- mosaic(donors[[g$species]], g$lambda, N_OUT)
  list(ld = ld_eval(G)[, `:=`(species = g$species, lambda = g$lambda, set = g$set)],
       id = data.table(species = g$species, lambda = g$lambda, set = g$set, identity = identity_mean(G)))
})
mos <- rbindlist(lapply(res, `[[`, "ld")); ids <- rbindlist(lapply(res, `[[`, "id"))

lab_bin <- function(b) ifelse(b == 99L, "unlinked", ifelse(b == 101L, "<5kb same sim cluster",
                       ifelse(b == 102L, "<5kb different sim cluster", CM_LABELS[pmin(b, length(CM_LABELS))])))
mos_s <- mos[, .(r2 = mean(r2)), by = .(species, lambda, bin)]
tab <- dcast(rbind(emp[, .(species, lambda = -1, bin, r2)], mos_s), species + bin ~ lambda, value.var = "r2")
setnames(tab, c("-1"), "empirical")
tab[, bin := lab_bin(bin)]
cat("\n[mosaic] mean within-species r^2: empirical parents vs mosaic founders (columns = switches per cM; n = 15):\n")
print(tab, digits = 3)
cat("\n[mosaic] mean pairwise genotype identity among founders (empirical parents:",
    paste(sprintf("%s %.3f", names(emp_id), emp_id), collapse = ", "), ")\n")
print(dcast(ids[, .(identity = mean(identity)), by = .(species, lambda)], species ~ lambda, value.var = "identity"), digits = 3)
saveRDS(list(table = tab, mosaic = mos, identity = ids, emp_identity = emp_id, lambdas = LAMBDAS),
        file.path(OUT_DATA, "explore_mosaic_founders.rds"))

pd <- rbind(emp[, .(species, lambda = "empirical parents", bin, r2)],
            mos_s[lambda %in% c(0, 0.1, 0.5, 2), .(species, lambda = sprintf("mosaic, %g switches/cM", lambda), bin, r2)])
pd <- pd[bin <= 99L][, bin_lab := factor(lab_bin(bin), levels = c(CM_LABELS, "unlinked"))]
p <- ggplot(pd, aes(bin_lab, r2, colour = lambda, group = lambda)) + geom_line() + geom_point(size = 1.2) +
  facet_wrap(~ species) + scale_y_log10() +
  scale_colour_manual(values = c("black", "#fdbe85", "#fd8d3c", "#e6550d", "#a63603"), name = NULL) +
  labs(x = "genetic distance between unit SNPs", y = expression("mean within-species "*r^2*" (log), n = 15")) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
ggsave(file.path(OUT_FIG, "explore_mosaic_founders.png"), p, width = 10, height = 4.8, dpi = 200)
cat("[mosaic] done\n")
