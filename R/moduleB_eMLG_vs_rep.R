## =========================================================
## Validation: eMLG consensus vs representative SNP
## =========================================================
## Do we gain anything (for the SORTING analysis) by using the eMLG consensus
## instead of the cheaper LD-pruning representative (tag SNP)? Two questions:
##
##   Q1 CALL concordance (decision-relevant): does the cluster-level sort_class
##      from the eMLG CONSENSUS differ from the sort_class of the cluster's
##      REPRESENTATIVE SNP? Cross-tabbed overall and by cluster size. If they
##      agree, the representative suffices for sorting; if they disagree --
##      especially large clusters where the consensus says "unsorted" but the
##      representative sorts -- the eMLG earns its keep by removing
##      pseudo-replication.
##
##   Q2 DOSAGE fidelity (mechanistic): cor(consensus, representative) across the
##      164 hybrids, by cluster size. Sorting is a THRESHOLD call, so a high
##      correlation does not guarantee an equal sort_class -- this just says how
##      much of the consensus a single tag SNP recovers.
##
## NB: eMLGs remain REQUIRED for the Ohta / among-region LD analyses (Module D
## needs consensus dosages). This asks only about the sorting-unit choice.
##
## Reads moduleA_clusters.rds, moduleA_snp.rds, the clustering (Q1, cheap) and
## eMLG_sorted_cM05.rds + GTs_hybrids (Q2; ~1.6 GB peak). Run from repo root.
##
## ---------------------------------------------------------------------------
## CONCLUSIONS (run 2026-07-21; parent_maf>=0.15 gate, canonical _cM05 clustering)
## ---------------------------------------------------------------------------
## Do we gain anything (for SORTING) from the eMLG consensus vs the tag SNP?
##
## Q1 CALL concordance: exact sort_class match 94.3% overall, but that is carried
##    by singletons (60% of units, 100% match by construction). Among multi-marker
##    clusters exact match is 71-84%; the disagreement is almost entirely
##    sorted-vs-unsorted (~10-21% "representative sorts, consensus washes out"),
##    NOT direction -- opposite-direction (aqu<->pol) calls are <1%.
## Q2 DOSAGE fidelity |cor|(consensus, representative), 164 hybrids: 2-4 loci
##    ~0.96 (near-perfect proxy); 5-49 loci ~0.83 (moderate, <5% below 0.7);
##    50+ ~0.91 (large blocks are coherent). A 9.5-22.7% NEGATIVE-cor tail is
##    arbitrary allele polarity (tag SNP coded opposite to the consensus
##    reference); harmless to the sorting call (parents re-orient both) but
##    corrupting for any RAW-dosage LD analysis -> use the consensus there.
##
## => DIRECTION (DI->direction; Module B B3): NO gain -- rep == consensus ~99%,
##    robust to unit choice.
## => MAGNITUDE (sorted vs unsorted, >=5-loci blocks): REAL gain -- the consensus
##    gives the principled block-integrated call; a single tag SNP is one
##    arbitrary member and mis-calls 15-30% of multi-marker clusters.
## => Module D (Ohta/LD): eMLGs REQUIRED -- consistent orientation (the sign-flip
##    rate above) + block integration.
##
## PIPELINE DECISION (validates the existing A/B split; no redesign):
##  - Module A definitive per-cluster calls: eMLG consensus (magnitude matters).
##  - Module B genome-wide sorting-vs-recomb + direction: representatives OK, but
##    read the recombination TRENDS and the SNP-vs-unit FLATTENING (robust), not
##    absolute frac_sorted magnitudes for large blocks (unit-definition dependent).
##  - Sorted-cluster set feeding Module C (climate) and Module D (Ohta): consensus.
## ---------------------------------------------------------------------------

suppressMessages({ library(data.table); library(ggplot2) })

cl  <- readRDS("data/moduleA_clusters.rds")   # consensus calls; unit_id = group_id
snp <- readRDS("data/moduleA_snp.rds")        # per-SNP calls
g   <- readRDS("data/eMLG_5loci_0025_cM05.rds")$groups

uni <- c("aquilonia", "polyctena")
size_brk <- c(0, 1, 4, 9, 49, Inf); size_lab <- c("1", "2-4", "5-9", "10-49", "50+")

## representative marker per unit, and its SNP-level call
cl[, rep_marker := g[.(unit_id), on = "group_id", representative]]
cl[snp, on = .(rep_marker = marker), `:=`(rep_sort = i.sort_class, rep_diff = i.differentiated)]
cl[, size_class := cut(n_loci, breaks = size_brk, labels = size_lab)]

## ========================================================================
## Q1 -- call concordance (consensus vs representative), where both classified
## ========================================================================
q1 <- cl[differentiated == TRUE & !is.na(sort_class) & rep_diff == TRUE & !is.na(rep_sort)]
cat("Units with both consensus and representative classified:", nrow(q1), "\n")

cat("\n=== Q1: exact sort_class agreement (consensus vs representative) ===\n")
cat("overall exact-match:", round(100 * mean(q1$sort_class == q1$rep_sort), 1), "%\n")
cat("\nby cluster size:\n")
print(q1[, .(n = .N,
             exact_match_pct  = round(100 * mean(sort_class == rep_sort), 1),
             sorted_agree_pct = round(100 * mean((sort_class != "unsorted") ==
                                                  (rep_sort != "unsorted")), 1),
             ## the pseudo-replication case: representative sorts but consensus does not
             rep_sorts_consensus_unsorted_pct =
               round(100 * mean(rep_sort != "unsorted" & sort_class == "unsorted"), 1)
            ), by = size_class][order(size_class)])

cat("\n=== Q1: full cross-tab (rows = consensus, cols = representative) ===\n")
print(dcast(q1, sort_class ~ rep_sort, value.var = "unit_id", fun.aggregate = length))

## ========================================================================
## Q2 -- dosage fidelity cor(consensus, representative) across hybrids
## ========================================================================
cat("\n[Q2] loading consensus matrix + hybrid genotypes for the dosage correlation ...\n")
se <- readRDS("data/eMLG_sorted_cM05.rds")            # $eMLG = hybrids x units
E  <- se$eMLG
e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1)
GTh <- e1$GTs_hybrids_005                             # hybrids x markers

## multi-marker units only (a singleton's consensus IS its representative -> cor 1)
mm <- cl[n_loci > 1 & unit_id %in% colnames(E) & rep_marker %in% colnames(GTh)]
ind <- rownames(E)
cor_rep <- vapply(seq_len(nrow(mm)), function(i) {
  cons <- E[ind, mm$unit_id[i]]
  rep  <- GTh[ind, mm$rep_marker[i]]
  suppressWarnings(stats::cor(cons, rep, use = "pairwise.complete.obs"))
}, numeric(1))
mm[, `:=`(cor_rep = cor_rep, abscor = abs(cor_rep))]

## Dosage sign is arbitrary: a representative coded on the opposite allele from
## the consensus's internal reference gives a NEGATIVE cor despite high shared
## information. parallelism_stats re-orients both by the parents, so the sign is
## irrelevant to the sorting call (cf. Q1's ~99% direction agreement). The honest
## fidelity metric is |cor|; the sign-flipped share is reported because it matters
## for raw-dosage LD (Module D), where the consistently-oriented consensus -- not
## the tag SNP -- must be used.
cat("\n=== Q2: |cor|(consensus, representative) across 164 hybrids, multi-marker units ===\n")
print(mm[, .(n = .N,
             median_abscor     = round(median(abscor, na.rm = TRUE), 3),
             pct_abscor_ge_0.9 = round(100 * mean(abscor >= 0.9, na.rm = TRUE), 1),
             pct_abscor_lt_0.7 = round(100 * mean(abscor <  0.7, na.rm = TRUE), 1),
             pct_sign_flipped  = round(100 * mean(cor_rep <  0,   na.rm = TRUE), 1)
            ), by = size_class][order(size_class)])

## small figure: cor distribution by size class
dir.create("Figures", showWarnings = FALSE)
p <- ggplot(mm[!is.na(cor_rep)], aes(size_class, cor_rep)) +
  geom_violin(fill = "#0072B2", alpha = 0.25, colour = NA) +
  geom_boxplot(width = 0.25, outlier.size = 0.3, linewidth = 0.3) +
  labs(x = "cluster size (n loci)", y = "cor(consensus, representative)") +
  theme_classic(base_size = 9)
ggsave("Figures/eMLG_vs_rep_cor.png", p, width = 90, height = 70, units = "mm", dpi = 300)

saveRDS(list(q1 = q1[, .(unit_id, n_loci, size_class, sort_class, rep_sort)],
             q2 = mm[, .(unit_id, n_loci, size_class, sort_class, rep_sort, cor_rep, abscor)]),
        "data/moduleB_eMLG_vs_rep.rds")
cat("\nSaved: data/moduleB_eMLG_vs_rep.rds, Figures/eMLG_vs_rep_cor.png\n")
