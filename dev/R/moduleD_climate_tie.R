## Tie the EMMAX trans-hub (F33028) to Module C (the extrinsic climate arm): is the
## hub's co-sorting module a climate signal, or not? Tests (a) at the LOCUS level whether
## the dominant genome ancestry axes align with the BayPass climate association, and
## (b) whether F33028 + its module are climate outliers. Reads the primary BayPass config
## (aland_excluded x withOmega), the clustering, moduleC_C3_cl.rds, moduleD_emmax.rds.
## Run from repo root. RESULT (2026-07-25): NEGATIVE -- F33028 and its module are among the
## LEAST climate-associated loci, and the genome axes barely track climate at the locus
## level; the individual-level genome-PC ~ climate-PC correlation is geographic entanglement,
## not locus-specific selection. => the hub belongs to neither arm (candidate rearrangement).
suppressMessages(library(data.table))

e2 <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2); map <- e2$map_hyb_005
imp <- function(f) { r <- fread(f); stopifnot(nrow(r) == nrow(map)); setNames(r$`BF(dB)`, map$marker) }
bf1 <- imp("aland_excluded/PC1_DIEM_withOmega_summary_betai_reg.out")   # climate PC1, primary config
bf2 <- imp("aland_excluded/PC2_DIEM_withOmega_summary_betai_reg.out")   # climate PC2
clust <- readRDS("data/eMLG_5loci_0025_cM05.rds"); groups <- clust$groups; eMLG <- clust$eMLG
cl <- readRDS("data/moduleC_C3_cl.rds"); E <- readRDS("data/moduleD_emmax.rds")

## per-cluster mean climate Bayes factor (continuous; the binary BF>=15 rate is too sparse)
m2g <- groups[, .(marker = unlist(members)), by = group_id][, `:=`(b1 = bf1[marker], b2 = bf2[marker])]
cs <- m2g[, .(mb1 = mean(b1, na.rm = TRUE), mb2 = mean(b2, na.rm = TRUE)), by = group_id]; setkey(cs, group_id)
diff_ids <- cl[differentiated == TRUE, group_id]

## (a) locus level: do the dominant genome ancestry axes align with the climate association?
Xg <- eMLG[, diff_ids]; Xg <- apply(Xg, 2, function(v) { v[!is.finite(v)] <- mean(v, na.rm = TRUE); v })
pc <- prcomp(Xg, scale. = FALSE)
ld <- data.table(group_id = diff_ids,
                 l1 = abs(as.numeric(cor(Xg, pc$x[, 1]))),
                 l2 = abs(as.numeric(cor(Xg, pc$x[, 2])))) [cs, on = "group_id", nomatch = 0]
cat("=== (a) locus level: dominant genome axes vs BayPass climate association ===\n")
cat(sprintf("cor(|genome-PC1 loading|, climate-PC1 mean BF) = %.2f\n", cor(ld$l1, ld$mb1, use = "pairwise")))
cat(sprintf("cor(|genome-PC2 loading|, climate-PC2 mean BF) = %.2f  (weak => geographic entanglement, not selection)\n",
            cor(ld$l2, ld$mb2, use = "pairwise")))

## (b) is the F33028 module a climate outlier set?
sig <- 0.05 / E$params$scope_n
mod <- c("F33028", E$results[["F33028"]][unlinked == TRUE & pval < sig, partner])
cat(sprintf("\n=== (b) F33028 module (n=%d) vs genome-wide differentiated background ===\n", length(mod)))
for (nm in c("mb1", "mb2")) cat(sprintf("climate-%s mean BF: module median %.2f vs background %.2f (Wilcoxon p=%.2g)\n",
    ifelse(nm == "mb1", "PC1", "PC2"), median(cs[.(mod), get(nm)], na.rm = TRUE),
    median(cs[.(diff_ids), get(nm)], na.rm = TRUE), wilcox.test(cs[.(mod), get(nm)], cs[.(diff_ids), get(nm)])$p.value))
cat(sprintf("F33028 itself: PC1 BF pctile %.2f, PC2 BF pctile %.2f (low => NOT a climate outlier)\n",
    mean(cs[.(diff_ids), mb1] < cs[.("F33028"), mb1], na.rm = TRUE),
    mean(cs[.(diff_ids), mb2] < cs[.("F33028"), mb2], na.rm = TRUE)))
