## =============================================================================
## Module E -- export the empirical NEST (colony) structure per population
## =============================================================================
## The empirical samples are workers from a few colonies, not random individuals:
## 164 samples across only 39 nests, median 2 nests/population and 4 workers/nest,
## and 9 of 20 populations are a SINGLE colony. Workers within a colony are
## siblings, so a nest sample has far less diversity than an equal number of
## random individuals, and its allele frequencies are family-biased. The
## simulation must be sampled the same way or pi / F_ST / sort_class are not
## comparable. This writes the per-population nest sizes for the sim harness.
## =============================================================================
suppressMessages(library(data.table))
e <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e)
sd <- as.data.table(e$sample_data)[!Population %in% c("aquilonia_parent","polyctena_parent")]
sd[, nest := sub("_[^_]+$", "", Sample_ID)]        # strip the individual suffix
ns <- sd[, .N, by = .(Population, nest)][order(Population, -N)]
out <- ns[, .(n_nests = .N, n_samples = sum(N),
              nest_sizes = paste(N, collapse = ",")), by = Population][order(Population)]
fwrite(out, "data/moduleE_sim/empirical_nest_structure.csv")
print(out)
cat(sprintf("\ntotal %d samples in %d nests across %d populations\n",
            sum(out$n_samples), sum(out$n_nests), nrow(out)))
