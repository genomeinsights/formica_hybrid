## =========================================================
## module_di25 -- best single-SNP genotype, missing data filled from the consensus
## =========================================================
## For every eMLG block, take the member SNP most correlated with the consensus
## (di25_emlg_best_member.rds) and replace ONLY its missing calls with the block
## consensus. Result: a SNP-level genotype (keeps SNP-specific signal) with the
## block eMLG repairing its missingness.
##
## Handling:
##   * ORIENTATION -- the consensus is polarized to an arbitrary within-cluster
##     reference, so the best SNP may be negatively correlated with it. We orient
##     the consensus to the SNP (2 - consensus when r < 0) before filling.
##   * SCALE -- consensus is on the same 0-2 dosage scale; fills are rounded to
##     0/1/2 so the output stays an integer genotype.
##   * Observed SNP calls are never overwritten; only NA cells are filled.
##
## Run from the repo root:  Rscript module_di25/R/di25_bestSNP_consensus_filled.R
## =========================================================

suppressMessages({ library(data.table) })

OUTDIR <- "module_di25/data"
CM     <- "cM5"

inp    <- readRDS(file.path(OUTDIR, "di25_inputs.rds"))
GTs    <- inp$GTs_hyb                                        # hybrids x markers (012, NA)
res    <- readRDS(file.path(OUTDIR, sprintf("di25_clustering_%s.rds", CM)))
emlg   <- res$eMLG                                           # hybrids x blocks (consensus, 0-2)
dt     <- as.data.table(readRDS(file.path(OUTDIR, "di25_emlg_best_member.rds")))

G  <- GTs[rownames(emlg), , drop = FALSE]                   # align hybrids to consensus rows
gi <- dt[match(colnames(emlg), group_id)]                  # best_marker per block, same order

n_ind   <- nrow(emlg)
n_block <- ncol(emlg)
filled  <- matrix(NA_real_, n_ind, n_block,
                  dimnames = list(rownames(emlg), gi$best_marker))

flipped     <- logical(n_block)
n_obs       <- integer(n_block)   # calls taken from the SNP
n_fill      <- integer(n_block)   # NA calls repaired from consensus
n_resid_na  <- integer(n_block)   # still NA (consensus also missing)

for (j in seq_len(n_block)) {
  snp  <- G[, gi$best_marker[j]]
  cons <- emlg[, j]

  ## orient consensus to the SNP's allele coding
  r <- suppressWarnings(cor(snp, cons, use = "pairwise.complete.obs"))
  if (is.finite(r) && r < 0) { cons <- 2 - cons; flipped[j] <- TRUE }

  miss <- is.na(snp)
  out  <- snp
  out[miss] <- pmin(pmax(round(cons[miss]), 0), 2)          # rounded, clamped 0-2

  filled[, j]  <- out
  n_obs[j]     <- sum(!miss)
  n_fill[j]    <- sum(miss & !is.na(cons))
  n_resid_na[j]<- sum(is.na(out))
}

## integer genotype matrix
storage.mode(filled) <- "integer"

map <- data.table(
  group_id    = colnames(emlg),
  best_marker = gi$best_marker,
  n_loci      = gi$n_loci,
  best_abs_r  = gi$best_abs_r,
  flipped     = flipped,
  n_obs       = n_obs,
  n_filled    = n_fill,
  n_resid_na  = n_resid_na
)

saveRDS(list(geno = filled, map = map, cM = CM),
        file.path(OUTDIR, "di25_bestSNP_consensus_filled.rds"))

## ---- report ---------------------------------------------------------------
tot <- n_ind * n_block
cat(sprintf("\n===== [%s] best-SNP genotype, consensus-filled (%d hybrids x %d blocks) =====\n",
            CM, n_ind, n_block))
cat(sprintf("Cells total          : %d\n", tot))
cat(sprintf("  observed (from SNP): %d  (%.2f%%)\n", sum(n_obs), 100*sum(n_obs)/tot))
cat(sprintf("  filled (consensus) : %d  (%.2f%%)\n", sum(n_fill), 100*sum(n_fill)/tot))
cat(sprintf("  still missing      : %d  (%.3f%%)\n", sum(n_resid_na), 100*sum(n_resid_na)/tot))
cat(sprintf("\nMissingness: best-SNP raw %.2f%%  ->  after fill %.3f%%\n",
            100*sum(is.na(G[, gi$best_marker]))/tot, 100*sum(n_resid_na)/tot))
cat(sprintf("Blocks with orientation flip (r<0): %d (%.1f%%)\n",
            sum(flipped), 100*mean(flipped)))
cat(sprintf("Blocks fully complete after fill  : %d (%.1f%%)\n",
            sum(n_resid_na == 0), 100*mean(n_resid_na == 0)))
cat(sprintf("\nWrote %s\n  $geno = %d x %d integer matrix (colnames = best SNP marker)\n  $map  = per-block provenance table\n",
            file.path(OUTDIR, "di25_bestSNP_consensus_filled.rds"), n_ind, n_block))
