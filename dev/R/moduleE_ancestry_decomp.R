## =============================================================================
## Module E -- decompose diagnostic differentiation: ancestry axis vs sorting
## =============================================================================
## Neutral F_ST ~ 0 but diagnostic F_ST ~ 0.31, and drift is DI-blind, so the
## diagnostic differentiation is NOT drift. It is one of two non-drift sources:
##   (a) variable FOUNDING ADMIXTURE proportion -- a single ancestry number per
##       population, so p_ij = f_pol_i + a_j * dp_i  (rank-1; demographic), or
##   (b) locus-specific ancestry SORTING/selection -- structure beyond that axis.
##
## Decomposition of the among-population variance at diagnostic markers (oriented
## to the aquilonia allele):
##   V_total    = observed among-pop variance of freqs
##   V_sampling = binomial sampling variance p(1-p)/(2n)  (subtracted out)
##   V_true     = V_total - V_sampling
##   V_ancestry = variance explained by the single per-pop ancestry axis a_j
##   V_residual = V_true - V_ancestry  (genuine locus-specific = sorting)
## A residual F_ST and a PCA of the residuals show whether what's left is
## population-structured (repeatable sorting) or noise.
##
## Runs from the COMPACT bundle (no 1 GB dataset).
## Output: data/moduleE_sim/moduleE_ancestry_decomp.rds
##         Figures/moduleE_ancestry_decomp.pdf
## =============================================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

FORMICA <- "/Users/petrikem/gitlab/formica_hybrid"
BUNDLE  <- "/Users/petrikem/gitlab/formica_hybrid/moduleE_slim/inputs/empirical_bundle.rds"
OUTRDS  <- file.path(FORMICA, "data/moduleE_sim/moduleE_ancestry_decomp.rds")
FIGDIR  <- file.path(FORMICA, "Figures"); dir.create(FIGDIR, showWarnings=FALSE)

b <- readRDS(BUNDLE); uni <- as.data.table(b$universe)
pops <- b$pops; n2 <- b$emp_ns * 2                       # allele counts per pop

## ---- diagnostic markers, oriented to the aquilonia allele ------------------
sel <- which(uni$DI > -25)
faq <- b$f_aq_par[uni$marker[sel]]; fpol <- b$f_pol_par[uni$marker[sel]]
P   <- (b$emp_mean[, uni$marker[sel], drop=FALSE]) / 2000   # pops x markers freq
flip <- (faq - fpol) < 0
faq [flip] <- 1 - faq[flip]; fpol[flip] <- 1 - fpol[flip]
P[, flip] <- 1 - P[, flip]
dp <- faq - fpol                                          # > 0 (aquilonia allele)
keep <- is.finite(dp) & dp > 0.05 & colSums(!is.finite(P)) == 0   # informative + complete across pops
P <- P[, keep, drop=FALSE]; dp <- dp[keep]; fpol <- fpol[keep]; faq <- faq[keep]
cat(sprintf("diagnostic markers used: %d ; %d populations\n", length(dp), length(pops)))

## ---- best-fit ancestry per population: p_ij - fpol_i = a_j * dp_i -----------
a <- sapply(seq_along(pops), function(j) {
  y <- P[j, ] - fpol; sum(y * dp) / sum(dp^2) })          # LS ancestry
names(a) <- pops
pred <- outer(a, dp) + matrix(fpol, nrow=length(pops), ncol=length(dp), byrow=TRUE)
resid <- P - pred                                         # locus-specific deviation

## ---- variance decomposition (per marker, then summed) ----------------------
J <- length(pops)
V_total_i  <- apply(P, 2, function(x) mean((x-mean(x))^2))            # among-pop var (/J)
V_pred_i   <- apply(pred, 2, function(x) mean((x-mean(x))^2))         # ancestry-explained
V_resid_i  <- apply(resid, 2, function(x) mean((x-mean(x))^2))        # residual (obs)
V_samp_i   <- colMeans(P*(1-P) / matrix(n2, J, length(dp)))           # binomial sampling
pbar <- colMeans(P); Ht_i <- 2*pbar*(1-pbar)

V_total <- sum(V_total_i); V_pred <- sum(V_pred_i)
V_resid <- sum(V_resid_i); V_samp <- sum(V_samp_i); Ht <- sum(Ht_i)
V_true  <- V_total - V_samp
V_resid_true <- max(V_resid - V_samp, 0)

fst_total    <- 2*V_total      / Ht                # F_ST-analog (var = He/2 scale)
fst_true     <- 2*V_true       / Ht
fst_ancestry <- 2*V_pred       / Ht
fst_residual <- 2*V_resid_true / Ht

## ---- residual structure: F_ST on residuals + PCA ---------------------------
## PCA of the residual matrix (pops x markers) -> do pops cluster beyond ancestry?
rr <- resid; rr[!is.finite(rr)] <- 0
pr <- prcomp(rr, center=TRUE, scale.=FALSE)
ve <- round(100*pr$sdev^2/sum(pr$sdev^2), 1)
pc <- data.table(Population=pops, ancestry=a, PC1=pr$x[,1], PC2=pr$x[,2])

saveRDS(list(a=a, decomp=list(V_total=V_total, V_true=V_true, V_ancestry=V_pred,
             V_residual_true=V_resid_true, V_sampling=V_samp),
             fst=c(total=fst_total, true=fst_true, ancestry=fst_ancestry, residual=fst_residual),
             pca=pc, ve=ve), OUTRDS)

## ---- report ----------------------------------------------------------------
cat(sprintf("\nper-population best-fit ancestry (aq): range %.3f - %.3f, sd %.3f\n",
            min(a), max(a), sd(a)))
cat("\n=== among-population DIAGNOSTIC variance decomposition ===\n")
cat(sprintf("  total (observed)          : %.5f\n", V_total))
cat(sprintf("  - sampling (binomial)     : %.5f  (%.0f%% of observed)\n", V_samp, 100*V_samp/V_total))
cat(sprintf("  = true among-pop variance : %.5f\n", V_true))
cat(sprintf("      of which ancestry axis: %.5f  (%.0f%% of true)\n", V_pred, 100*V_pred/V_true))
cat(sprintf("      of which locus-specific: %.5f  (%.0f%% of true)\n", V_resid_true, 100*V_resid_true/V_true))
cat("\n=== F_ST-analog decomposition ===\n")
cat(sprintf("  total F_ST     : %.3f\n", fst_total))
cat(sprintf("  ancestry-axis  : %.3f  (variable admixture, demographic)\n", fst_ancestry))
cat(sprintf("  residual F_ST  : %.3f  (locus-specific = sorting/selection)\n", fst_residual))
cat(sprintf("\nresidual PCA variance explained: PC1=%.1f%% PC2=%.1f%% PC3=%.1f%%\n", ve[1], ve[2], ve[3]))
cat("\nREAD: residual F_ST ~ 0 => variable admixture explains the differentiation (no selection).\n")
cat("      residual F_ST >> 0 with structured residual PCA => locus-specific parallel sorting.\n")

## figure: ancestry axis (top) + residual PCA (bottom)
p1 <- ggplot(data.table(Population=pops, a=a), aes(reorder(Population,a), a)) +
  geom_col(fill="#1b9e77") + coord_flip() +
  labs(x=NULL, y="best-fit aquilonia ancestry", title="Per-population ancestry (the single admixture axis)") +
  theme_bw(base_size=8)
p2 <- ggplot(pc, aes(PC1, PC2, colour=ancestry, label=Population)) +
  geom_point(size=2) + geom_text(size=2, vjust=-0.8) +
  scale_colour_viridis_c(option="plasma") +
  labs(title=sprintf("Residual PCA (after removing the ancestry axis): PC1=%.0f%% PC2=%.0f%%", ve[1], ve[2]),
       subtitle="structure here = locus-specific sorting beyond variable admixture") +
  theme_bw(base_size=8)
ggsave(file.path(FIGDIR,"moduleE_ancestry_decomp.pdf"), p1/p2, width=9, height=8)
ggsave(file.path(FIGDIR,"moduleE_ancestry_decomp.png"), p1/p2, width=9, height=8, dpi=130)
cat("\nsaved: figures + ", OUTRDS, "\n")
