## =========================================================================
## sim_founder_fix -- generate founders as GENOTYPE MOSAICS of the empirical parents
## (no phasing needed)
##
## Each founder is built chromosome by chromosome: start with a random donor among
## the 15 empirical parents of its species, copy the donor's genotypes, and switch to
## a new random donor at breakpoints placed along the genetic map as a Poisson
## process with rate LAMBDA switches per cM. Allele frequencies, heterozygosity and
## graded short-range LD therefore come from the empirical parents. LAMBDA trades LD
## fidelity (low) against new haplotype diversity (high); 0.5-1 per cM reproduced
## the empirical parental LD decay best (see NOTE_parent_ld.pdf).
##
## Ploidy, matching the SLiM founder set-up (p1 = F. aquilonia males, haploid;
## p2 = F. polyctena queens, diploid):
##   aquilonia males : one allele from the copied diploid genotype (0 -> 0, 2 -> 1,
##                     heterozygote -> random allele)
##   polyctena queens: the copied diploid genotype; heterozygotes split at random
##                     between the two haplosomes (phase among neighbouring
##                     heterozygous sites is random; heterozygosity at DI25 sites is
##                     0.05-0.12, so this affects few pairs)
## Missing donor genotypes (2-3%) are filled, per site, from a random donor of the
## same species with an observed genotype.
## Allele "1" = the coded allele of the empirical genotype matrices (dosage counts it).
##
## Output: <OUT_DIR>/founders_ch<id>.vcf in the format of the earlier real-founder
## SLiM script (SpecIAnt_rufa_neutral_realfounders.slim, moduleE_slim/), read with readHaplosomesFromVCF
## (columns aq_hapNN haploid, then pol_femNN phased diploid; POS = genome position),
## plus <OUT_DIR>/provenance.txt.
##
## Run from the formica_hybrid repo root:
##   Rscript sim_founder_fix/make_mosaic_founders.R <OUT_DIR> [LAMBDA=1] [N_AQ=50] [N_POL=50] [SEED=1] [PANEL=DI25]
## PANEL = "DI25" (51,612 ancestry-informative SNPs) or "DI25+neutral" (adds ~14,100
## near-neutral SNPs, DI <= -90 and pooled parental MAF >= 0.15, as a calibration
## anchor; their IDs are listed in <OUT_DIR>/neutral_markers.txt).
## =========================================================================
source("sim_founder_fix/parent_ld_lib.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript sim_founder_fix/make_mosaic_founders.R <OUT_DIR> [LAMBDA] [N_AQ] [N_POL] [SEED]")
OUT_DIR <- args[1]
LAMBDA  <- if (length(args) >= 2) as.numeric(args[2]) else 1
N_AQ    <- if (length(args) >= 3) as.integer(args[3]) else 50L
N_POL   <- if (length(args) >= 4) as.integer(args[4]) else 50L
SEED    <- if (length(args) >= 5) as.integer(args[5]) else 1L
PANEL   <- if (length(args) >= 6) args[6] else "DI25"
stopifnot(PANEL %in% c("DI25", "DI25+neutral"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

panel <- load_panel(include_neutral = PANEL == "DI25+neutral"); map <- panel$map; P <- panel$parents
stopifnot(identical(colnames(P), map$marker))
chrs <- unique(map$Chr)
cols <- split(seq_len(nrow(map)), factor(map$Chr, levels = chrs))
## breakpoint coordinate: cM, with off-map SNPs given their nearest mapped neighbour's value
bcm <- map$cM
for (ix in cols) { v <- bcm[ix]; if (all(is.na(v))) { bcm[ix] <- 0; next }
                   if (anyNA(v)) bcm[ix] <- approx(which(!is.na(v)), v[!is.na(v)], seq_along(v), rule = 2)$y }

fill_missing <- function(D) {
  for (j in which(colSums(is.na(D)) > 0)) {
    obs <- which(!is.na(D[, j])); miss <- which(is.na(D[, j]))
    D[miss, j] <- if (length(obs)) D[obs[sample.int(length(obs), length(miss), replace = TRUE)], j] else 0L
  }
  D
}
donors <- list(aquilonia = fill_missing(P[grepl("^Faqu", rownames(P)), ]),
               polyctena = fill_missing(P[grepl("^Fpol", rownames(P)), ]))

mosaic <- function(D, n_out) {
  out <- matrix(NA_integer_, n_out, ncol(D))
  for (f in seq_len(n_out)) for (ix in cols) {
    v <- bcm[ix]; lo <- min(v); hi <- max(v)
    nb <- if (LAMBDA > 0 && hi > lo) rpois(1, LAMBDA * (hi - lo)) else 0L
    brk <- sort(runif(nb, lo, hi)); seg <- findInterval(v, brk) + 1L
    donor <- sample.int(nrow(D), nb + 1L, replace = TRUE)
    out[f, ix] <- D[cbind(donor[seg], ix)]
  }
  out
}
A <- mosaic(donors$aquilonia, N_AQ)       # diploid dosages, to be reduced to haploid
Q <- mosaic(donors$polyctena, N_POL)

aq_hap <- ifelse(A == 0L, 0L, ifelse(A == 2L, 1L, rbinom(length(A), 1, 0.5)))
dim(aq_hap) <- dim(A)
q1 <- ifelse(Q == 2L, 1L, ifelse(Q == 0L, 0L, rbinom(length(Q), 1, 0.5))); dim(q1) <- dim(Q)
q2 <- Q - q1                                    # second haplosome: the remaining allele
stopifnot(all(q2 %in% 0:1))

aq_names  <- sprintf("aq_hap%02d", seq_len(N_AQ)); pol_names <- sprintf("pol_fem%02d", seq_len(N_POL))
for (k in seq_along(chrs)) {
  ix <- cols[[k]]; chs <- sub("Chr", "ch", chrs[k])
  gt <- cbind(t(aq_hap[, ix, drop = FALSE]), matrix(paste0(t(q1[, ix, drop = FALSE]), "|", t(q2[, ix, drop = FALSE])), length(ix)))
  body <- paste(sprintf("%s\t%d\t%s\tA\tC\t.\tPASS\t.\tGT", chs, map$Pos[ix], map$marker[ix]),
                apply(gt, 1, paste, collapse = "\t"), sep = "\t")
  writeLines(c("##fileformat=VCFv4.2", sprintf("##contig=<ID=%s>", chs),
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", aq_names, pol_names), collapse = "\t"),
               body), file.path(OUT_DIR, sprintf("founders_%s.vcf", chs)))
}
writeLines(c(sprintf("generated: %s", format(Sys.time())), sprintf("LAMBDA (switches per cM): %g", LAMBDA),
             sprintf("N_AQ haploid males: %d; N_POL diploid queens: %d; SEED: %d", N_AQ, N_POL, SEED),
             sprintf("donors: %d aquilonia, %d polyctena empirical parents; %d SNPs on %d chromosomes",
                     nrow(donors$aquilonia), nrow(donors$polyctena), nrow(map), length(chrs)),
             sprintf("panel: %s (%d DI25 + %d near-neutral SNPs)", PANEL, sum(map$class == "DI25"), sum(map$class == "neutral")),
             "allele 1 = coded allele of the empirical genotype matrices"),
           file.path(OUT_DIR, "provenance.txt"))
if (PANEL == "DI25+neutral") writeLines(map[class == "neutral", marker], file.path(OUT_DIR, "neutral_markers.txt"))
cat(sprintf("[mosaic] wrote %d chromosome VCFs (%d aq males, %d pol queens, lambda %g) to %s\n",
            length(chrs), N_AQ, N_POL, LAMBDA, OUT_DIR))
