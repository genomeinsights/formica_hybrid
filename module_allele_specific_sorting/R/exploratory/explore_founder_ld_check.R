## =========================================================================
## module_allele_specific_sorting -- EXPLORATORY: where does the excess short-range
## within-population LD of the neutral simulations come from? (chromosomes 1-6)
##
## (1) Same within-population LD statistic (sim_stats_lib.R::profile, ancestry-adjusted,
##     both loci segregating) for: empirical hybrids; empirical parents; a phased mosaic
##     founder pool (pure species, and synthetic first-generation hybrids); simulated
##     hybrids at cycle 60 with the standard founder switch rate (lambda = 1 per cM) and
##     with lambda = 20 per cM (more founder haplotype diversity).
##     Result (2026-10-07): DI25 founders match the empirical parents/hybrids
##     (0.16 at 0.001-0.01 cM); the excess (0.56) appears during the hybrid phase and
##     does not change with lambda = 20 -> admixture LD, not founder diversity.
## (2) F. polyctena only, near-neutral units: empirical parents, Beagle-phased parents,
##     a resample of them, and the founder queens -> founders reproduce parental LD.
## (3) Within-population ancestry variability at DI25 units: mean p(1-p), % unit x
##     population fixed / intermediate, empirical vs simulated (K 6,250, 1,000 founders).
## Inputs (git-ignored, sim_founder_fix/out/): founder_pool_phased_seed3/ (make_mosaic_founders.R,
## lambda 1, 50+50 founders, seed 3, phased), lam20_cycle60/ (4 runs, K 6,250, 1,000
## founders, lambda 20, chromosomes 1-6, 60 cycles), calib/ (calibration grid).
## Output: printed tables, data/explore_founder_ld_check.rds
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/exploratory/explore_founder_ld_check.R
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
POOL <- "sim_founder_fix/out/founder_pool_phased_seed3"; CAL <- "sim_founder_fix/out/calib"; LAM20 <- "sim_founder_fix/out/lam20_cycle60"
CHRS <- paste0("Chr", 1:6)
U <- load_units(); u25 <- readRDS(file.path(OUT_DATA, "01_units.rds")); GG <- load_oriented_genotypes(U)
inp <- readRDS("module_di25/data/di25_inputs.rds")
sgn_all <- sign(colMeans(inp$GTs_par[grepl("^Faqu", rownames(inp$GTs_par)), u25$unit_marker], na.rm = TRUE) -
                colMeans(inp$GTs_par[grepl("^Fpol", rownames(inp$GTs_par)), u25$unit_marker], na.rm = TRUE))
k <- u25$Chr %in% CHRS; u25 <- u25[k]; sgn <- sgn_all[k]
neu <- readRDS(file.path(OUT_DATA, "06_lowDI.rds"))$units[Chr %in% CHRS, .(marker = rep_snp, Chr, Pos, cM_pos)]
e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
sets <- list(DI25 = u25[, .(marker = unit_marker, Chr, Pos, cM_pos)], neutral = neu)
all_markers <- c(u25$unit_marker, neu$marker)
source("module_allele_specific_sorting/R/sim_stats_lib.R")
orient <- function(G) { G[, which(sgn < 0)] <- 2L - G[, which(sgn < 0)]; G[, which(is.na(sgn) | sgn == 0)] <- NA; G }
wld <- function(G, pop, label) {
  Gor <- orient(G[, u25$unit_marker])
  rbind(profile(Gor, pop, sets$DI25, Gor)[, set := "DI25"],
        profile(G[, neu$marker], pop, sets$neutral, Gor)[, set := "neutral"])[, .(set, bin, within_LD, data = label)]
}
set.seed(1)

## ---- (1) within-population LD by data source ----
out <- list()
Gh <- cbind(inp$GTs_hyb[rownames(GG$G), u25$unit_marker], e$GTs_with_parents[rownames(GG$G), neu$marker])
out$emp_hyb <- wld(Gh, GG$pop, "empirical hybrids (165, 20 pops)")
par <- rownames(inp$GTs_par)
Gp <- cbind(inp$GTs_par[, u25$unit_marker], e$GTs_with_parents[par, neu$marker])
out$emp_par <- wld(Gp, ifelse(grepl("^Faqu", par), "aq", "pol"), "empirical parents (15+15, 2 species)")
fs <- list.files(POOL, "^founders_ch[0-9]+\\.vcf$", full.names = TRUE)
v <- rbindlist(lapply(fs, fread, skip = "#CHROM", colClasses = "character"))[ID %in% all_markers]
aq <- grep("^aq_hap", names(v), value = TRUE); pq <- grep("^pol_fem", names(v), value = TRUE)
mk <- function(M) { G <- matrix(0L, nrow(M), length(all_markers), dimnames = list(NULL, all_markers)); G[, v$ID] <- M; G }
A <- t(sapply(aq, function(s) as.integer(v[[s]])))
Q1 <- t(sapply(pq, function(s) as.integer(substr(v[[s]], 1, 1)))); Q2 <- t(sapply(pq, function(s) as.integer(substr(v[[s]], 3, 3))))
na <- length(aq); np <- length(pq)
out$pool_pure <- wld(mk(rbind(A[seq(1, na, 2), ] + A[seq(2, na, 2), ], Q1 + Q2)), rep(c("aq", "pol"), c(na / 2, np)),
                     "founder pool, pure species (aq male pairs + pol queens)")
n <- 40; ia <- sample(na, n, TRUE); iq <- sample(np, n, TRUE); w <- rbinom(n, 1, 0.5)
out$pool_F1 <- wld(mk(A[ia, ] + w * Q1[iq, ] + (1 - w) * Q2[iq, ]), rep(sprintf("p%d", 1:4), each = 10),
                   "founder pool, synthetic F1 (40, 4 pops)")
simset <- function(files, label) {
  parts <- lapply(files, read_sim, n = 10); pop <- rep(sprintf("r%d", seq_along(files)), sapply(parts, nrow))
  wld(do.call(rbind, parts), pop, label)
}
out$sim1 <- simset(head(list.files(CAL, "^females_ckl60_K6k_N1000_[0-9]+\\.vcf\\.gz$", full.names = TRUE), 4), "sim cycle 60, lambda 1 (4 x 10)")
out$sim20 <- simset(list.files(LAM20, "^females_ckl60_lam20_.*\\.vcf\\.gz$", full.names = TRUE), "sim cycle 60, lambda 20 (4 x 10)")
ld <- rbindlist(out)[, bin := ifelse(bin == 99L, "unlinked", CM_LABELS[bin])]
options(width = 200)
cat("[founder-ld] (1) ancestry-adjusted within-population LD (both loci segregating)\n")
print(dcast(ld, set + data ~ factor(bin, levels = c(CM_LABELS, "unlinked")), value.var = "within_LD"), digits = 2)

## ---- (2) F. polyctena only, near-neutral ----
neuonly <- function(G, pop) { Gor <- orient(G[, u25$unit_marker]); profile(G[, neu$marker], pop, sets$neutral, Gor)$within_LD[1:3] }
two <- function(n) rep(c("a", "b"), length.out = n)                 # profile() needs >= 2 groups
pol <- par[grepl("^Fpol", par)]
Gpol <- cbind(inp$GTs_par[pol, u25$unit_marker], e$GTs_with_parents[pol, neu$marker])
ph <- readRDS("sim_founder_fix/out/phased/parents_phased.rds")$H
Hp <- ph[grepl("^Fpol", rownames(ph)), all_markers]; Gph <- Hp[c(TRUE, FALSE), ] + Hp[c(FALSE, TRUE), ]
Gq <- mk(Q1 + Q2); rs <- sample(nrow(Gph), np, TRUE)
pol_tab <- rbind(`empirical pol parents (15)` = neuonly(Gpol, two(15)),
                 `Beagle-phased pol parents (15)` = neuonly(Gph, two(15)),
                 `resample of phased parents` = neuonly(Gph[rs, ], two(np)),
                 `founder pol queens (mosaic)` = neuonly(Gq, two(np)))
colnames(pol_tab) <- CM_LABELS[1:3]
cat("\n[founder-ld] (2) F. polyctena only, near-neutral within-population LD\n"); print(round(pol_tab, 3))

## ---- (3) within-population ancestry variability at DI25 units ----
het <- function(G, pop) { Gor <- orient(G[, u25$unit_marker]) / 2
  Fm <- t(sapply(unique(pop), function(p) colMeans(Gor[pop == p, , drop = FALSE], na.rm = TRUE)))
  c(mean_within_pop_pq = mean(Fm * (1 - Fm), na.rm = TRUE), pct_fixed = 100 * mean(Fm <= 0.05 | Fm >= 0.95, na.rm = TRUE),
    pct_intermediate = 100 * mean(Fm > 0.2 & Fm < 0.8, na.rm = TRUE)) }
anc <- rbind(empirical = het(inp$GTs_hyb[rownames(GG$G), ], GG$pop))
for (cy in c(60, 500, 1000)) {
  f <- head(list.files(CAL, sprintf("^females_ckl%d_K6k_N1000_[0-9]+\\.vcf\\.gz$", cy), full.names = TRUE), 6)
  parts <- lapply(f, read_sim, n = 10)
  anc <- rbind(anc, het(do.call(rbind, parts), rep(seq_along(f), sapply(parts, nrow))))
  rownames(anc)[nrow(anc)] <- sprintf("sim K 6,250, 1,000 founders, cycle %d", cy)
}
cat("\n[founder-ld] (3) within-population ancestry variability at DI25 units\n"); print(round(anc, 3))
saveRDS(list(within_ld = ld, polyctena_neutral = pol_tab, ancestry_variability = anc),
        file.path(OUT_DATA, "explore_founder_ld_check.rds"))
cat("[founder-ld] done\n")
