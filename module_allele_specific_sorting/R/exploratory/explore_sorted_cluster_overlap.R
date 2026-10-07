## =========================================================================
## module_allele_specific_sorting -- EXPLORATORY: do the large sorted clusters
## overlap other candidate regions?
##
## Clusters: explore_sorted_cluster_sizes.R, gap <= 50 kb. "Large" = >= MIN_SORTED
## sorted units. Compared against
##   (1) BDMI candidate regions (colleague's merged intervals, 23 X2 cutoffs;
##       data/liftoff_Frufa_DTOL_PR/bdmi_candidates.cutoff_*.bed, Frufa_DTOL_PR coords);
##   (2) BayPass climate candidate regions (module_BayPass/tables/
##       formica_region_concordance.tsv; covariate, primary/secondary tier);
##   (3) the Chr7 among-region lead (F9614-F9879, Chr7:6,646,791-6,647,856).
## Overlap = interval intersection after padding clusters by PAD_KB on each side.
## Null for (1)-(2): each large cluster circularly shifted along its chromosome
## (N_ROT rotations; cluster sizes and target regions fixed) -> expected number of
## large clusters overlapping, and permutation p (one-sided, excess).
## Output: printed tables, data/explore_sorted_cluster_overlap.rds
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/exploratory/explore_sorted_cluster_overlap.R [MIN_SORTED=4] [PAD_KB=10] [N_ROT=2000]
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
args <- commandArgs(trailingOnly = TRUE)
MIN_SORTED <- if (length(args) >= 1) as.integer(args[1]) else 4L
PAD_KB <- if (length(args) >= 2) as.numeric(args[2]) else 10
N_ROT <- if (length(args) >= 3) as.integer(args[3]) else 2000L
set.seed(1)
cl <- readRDS(file.path(OUT_DATA, "explore_sorted_cluster_sizes.rds"))$clusters[["50"]]$observed
big <- cl[n_sorted >= MIN_SORTED][order(-n_sorted)]
big[, `:=`(s = start - PAD_KB * 1e3, e = end + PAD_KB * 1e3)]
## chromosome lengths (bp) for rotations, as in the SLiM model
CH_ID <- c(1:22, 24:27)
CH_LEN <- c(16449812, 18558483, 15805178, 16666556, 14957794, 10925039, 13372598, 13028660, 10584038, 13832892, 11830133,
            11260715, 11449430, 7730741, 10965718, 11699404, 11021207, 8644910, 7297899, 9820556, 10314430, 7671390,
            6269374, 13114028, 9969264, 7047176)
chlen <- setNames(CH_LEN, paste0("Chr", CH_ID))
cat(sprintf("[overlap] %d large sorted clusters (>= %d sorted units, gap <= 50 kb), padded +-%g kb\n", nrow(big), MIN_SORTED, PAD_KB))

ov_any <- function(ch, s, e, R) { if (!nrow(R)) return(rep(FALSE, length(ch)))
  mapply(function(c, a, b) any(R$Chr == c & R$start <= b & R$end >= a), ch, s, e) }
rot_null <- function(R) {
  obs <- sum(ov_any(big$Chr, big$s, big$e, R))
  w <- big$e - big$s; L <- chlen[big$Chr]
  nul <- replicate(N_ROT, { s <- (big$s + runif(nrow(big), 0, L)) %% L
    s <- pmin(s, L - w); sum(ov_any(big$Chr, s, s + w, R)) })
  c(observed = obs, expected = mean(nul), fold = obs / mean(nul), p = (1 + sum(nul >= obs)) / (N_ROT + 1),
    target_cover_pct = 100 * sum(R$end - R$start) / sum(chlen))
}

## ---- (1) BDMI candidate regions ----
bf <- list.files("data/liftoff_Frufa_DTOL_PR", "^bdmi_candidates\\.cutoff_.*\\.bed$", full.names = TRUE)
tok <- sub("^bdmi_candidates\\.cutoff_[0-9]+_([0-9]+)\\..*", "\\1", basename(bf))
x2 <- as.numeric(sub("^(0)(\\d+)$", "0.\\2", tok))
bed <- rbindlist(lapply(seq_along(bf), function(i) fread(bf[i], header = FALSE, col.names = c("chr", "start", "end"))[, X2 := x2[i]]))
bed[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
bd <- rbindlist(lapply(sort(unique(bed$X2), decreasing = TRUE), function(k) data.table(X2 = k, t(rot_null(bed[X2 == k])))))
cat("\n[overlap] (1) BDMI candidate regions: large sorted clusters overlapping, vs rotation null (higher X2 = more stringent)\n")
print(bd, digits = 3)
## per cluster: most stringent cutoff at which it overlaps a BDMI region
big[, bdmi_max_X2 := sapply(seq_len(.N), function(i) { k <- bed[Chr == big$Chr[i] & start <= big$e[i] & end >= big$s[i], X2]
  if (length(k)) max(k) else NA_real_ })]

## ---- (2) BayPass climate regions ----
bp <- fread("module_manuscript_rho05/module_BayPass/tables/formica_region_concordance.tsv")
bp <- bp[, .(covariate, region_id, Chr = chromosome, start = start_bp, end = end_bp, tier = candidate_tier,
             ld_reduced = ld_reduced != "")]
cat(sprintf("\n[overlap] (2) BayPass regions: %d (%s)\n", nrow(bp), paste(bp[, .N, by = covariate][, paste(covariate, N)], collapse = ", ")))
bpt <- rbind(data.table(set = "all regions", t(rot_null(bp))),
             data.table(set = "primary tier", t(rot_null(bp[tier == "primary"]))),
             rbindlist(lapply(unique(bp$covariate), function(v) data.table(set = v, t(rot_null(bp[covariate == v]))))))
print(bpt, digits = 3)
big[, baypass := sapply(seq_len(.N), function(i) { k <- bp[Chr == big$Chr[i] & start <= big$e[i] & end >= big$s[i]]
  if (nrow(k)) paste(unique(paste0(k$region_id, ifelse(k$tier == "", "", paste0("(", k$tier, ")")))), collapse = ",") else "" })]

## ---- (3) Chr7 among-region lead ----
lead <- data.table(Chr = "Chr7", start = 6646791, end = 6647856)
big[, dist_chr7_lead_kb := ifelse(Chr == "Chr7", pmax(0, pmax(lead$start - end, start - lead$end)) / 1e3, NA_real_)]

options(width = 220)
cat("\n[overlap] large sorted clusters and their overlaps\n")
print(big[, .(Chr, start, end, n_sorted, purity = round(purity, 2), span_kb = round(span_kb), share_aqu = round(share_aqu, 2),
              bdmi_max_X2, baypass, dist_chr7_lead_kb)], nrows = 200)
saveRDS(list(big = big, bdmi = bd, baypass = bpt, min_sorted = MIN_SORTED, pad_kb = PAD_KB, n_rot = N_ROT),
        file.path(OUT_DATA, "explore_sorted_cluster_overlap.rds"))
cat("[overlap] done\n")
