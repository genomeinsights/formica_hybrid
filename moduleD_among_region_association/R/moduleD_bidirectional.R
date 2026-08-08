## =========================================================
## MODULE D -- bidirectional-locus ANNOTATION (not a DMI screen).
## =========================================================
## Produces the set of LD clusters that sort BIDIRECTIONALLY: strong sorting in RANDOM
## directions across replicate populations (each deme independently fixes one of two
## equally-fit compatible combinations) -- the equal-fitness-DMI trace, for which the
## single-locus parallelism test is null by construction (n_aqu ~ n_pol). Used ONLY as a
## descriptive annotation: the "bidirectional" category in the sorting ring of the network
## circos and per-module heatmaps reads data/moduleD_bidirectional.rds$reps$group_id.
##
## NOTE: the SNP-level symmetric-DMI PAIR test (pairwise co-sorting R_st + per-locus
## permutation null + candidate network) was explored and SET ASIDE -- a weak, diffuse
## excess whose permutation null overstates significance by not modelling shared drift
## (manuscript_notes/moduleD_plan.md). Only the bidirectional-locus identification remains.
##
## WHY SNP-LEVEL: bidirectional sorting is real at the SNP level but is averaged below the
## sort threshold at the eMLG consensus, so per LD cluster we take its MOST bidirectional
## member SNP as the representative. Run from repo root. Reads moduleA_snp.rds + clustering.

suppressMessages(library(data.table))

snp   <- readRDS("data/moduleA_snp.rds")
clust <- readRDS("data/eMLG_5loci_0025_cM05.rds"); groups <- clust$groups

## bidirectional reps: one (most bidirectional) SNP per LD cluster
bi <- snp[sort_class == "bidirectional"]
bi[, c("Chr", "Pos") := tstrsplit(marker, ":", fixed = TRUE)][, Pos := as.integer(Pos)]
m2g <- groups[, .(marker = unlist(members)), by = group_id]
bi <- m2g[bi, on = "marker"]
reps <- bi[order(-bi_score, -prop_fixed)][, .SD[1], by = group_id][
  , .(group_id, marker, Chr, Pos, DI, bi_score, prop_fixed, f_aqu_pooled)]
message(sprintf("[bi] %d bidirectional reps (1 per cluster) over %d chromosomes",
                nrow(reps), reps[, uniqueN(Chr)]))

saveRDS(list(reps = reps,
             params = list(source = "moduleA_snp sort_class==bidirectional, 1 rep/cluster",
                           role = "descriptive annotation only, not a DMI screen")),
        "moduleD_among_region_association/data/moduleD_bidirectional.rds")
message("[bi] saved data/moduleD_bidirectional.rds (bidirectional-locus annotation)")
