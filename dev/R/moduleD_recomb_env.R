## Recombination environment of the EMMAX trans-hub (F33028 + module) and the Table-1
## pair-specific candidates: do they map to low-recombination (centromeric) regions?
## If the hub blocks are centromeric, their large LD is the reference recombination
## landscape (a fixed feature), not a segregating inversion. Reads the clustering,
## moduleC_C3_cl.rds, moduleD_emmax.rds, moduleD_nonhub_trans.rds, the recmap. Run from repo root.
## RESULT (2026-07-25): F33028 at 0.4 cM/Mb (1st percentile); 33/37 module clusters in the
## low-recomb tail (centromeric) => NOT an inversion, but low-recomb regions retaining founding
## admixture co-ancestry. Table-1 candidates are normal-to-high recombination (real signal).
suppressMessages(library(data.table))

clust <- readRDS("data/eMLG_5loci_0025_cM05.rds"); groups <- clust$groups
cl <- readRDS("data/moduleC_C3_cl.rds"); E <- readRDS("data/moduleD_emmax.rds")
e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1)
map <- as.data.table(e1$map_hyb_005)[, .(marker, Chr = as.character(Chr), Pos)]
rec <- fread("data/Frufa_DTOL_PR.ref_genome.recmap"); setnames(rec, 1:4, c("chr","pos","cM","cMMb")); rec[, Chr := sub("chromosome_","Chr",chr)]
## recombination RATE (cM/Mb) per marker -> per-cluster median
map[, rate := NA_real_]
for (ch in unique(map$Chr)) { r <- rec[Chr == ch]; if (nrow(r) < 2) next
  ix <- map[, which(Chr == ch)]; map[ix, rate := approx(r$pos, r$cMMb, xout = Pos, rule = 2)$y] }
rr <- setNames(map$rate, map$marker)
crate <- groups[, .(marker = unlist(members)), by = group_id][, rate := rr[marker]][
  , .(rate = median(rate, na.rm = TRUE)), by = group_id][, setNames(rate, group_id)]
diff_ids <- cl[differentiated == TRUE, group_id]; bg <- crate[diff_ids]
pctile <- function(g) mean(bg < crate[g], na.rm = TRUE)
cat(sprintf("genome-wide recomb rate (differentiated clusters): median %.1f cM/Mb (10th %.1f, 90th %.1f)\n",
    median(bg, na.rm = TRUE), quantile(bg, 0.1, na.rm = TRUE), quantile(bg, 0.9, na.rm = TRUE)))

sig <- 0.05 / E$params$scope_n
mod <- c("F33028", E$results[["F33028"]][unlinked == TRUE & pval < sig, partner])
cat(sprintf("\n=== hub F33028 + module (n=%d) ===\n", length(mod)))
cat(sprintf("F33028 (Chr10): %.1f cM/Mb, percentile %.2f\n", crate["F33028"], pctile("F33028")))
cat(sprintf("module clusters in the low-recomb tail (pctile<0.1): %d/%d; module median percentile %.2f\n",
    sum(sapply(mod, pctile) < 0.1), length(mod), median(sapply(mod, pctile))))

cat("\n=== Table-1 pair-specific candidates (recombination percentiles) ===\n")
nh <- readRDS("data/moduleD_nonhub_trans.rds")
print(nh[, .(pair = paste(focal, partner),
             rate = sprintf("%.1f/%.1f", crate[focal], crate[partner]),
             pctile = sprintf("%.2f/%.2f", sapply(focal, pctile), sapply(partner, pctile)),
             disposition)])
