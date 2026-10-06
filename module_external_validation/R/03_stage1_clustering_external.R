## module_external_validation / 03: Stage-1 LD clustering (ld_complexity_reduction,
## rho = 0.5) of the external hybrid samples, MAF > 0.1 SNPs, using the LD-decay fit
## from 02_ld_decay_external.R (slide = 100: physical window ~128 kb, comparable to
## the original ~105 kb). Edge lists are rebuilt on the fly from the GDS.
##
## The cluster representatives (`pruned`) are the LD-pruned marker set intended for
## the EMMAX relatedness matrix.
##
## Run from the repo root:  Rscript module_external_validation/R/03_stage1_clustering_external.R [slide]
## Reads : data/ext_geno_maf10.rds, data/ld_decay_ext_maf10_slide<slide>.rds, data/ext_maf10_slide<slide>.gds
## Writes: data/ext_stage1_maf10_rho05.rds   (ld_complexity_reduction object)
##         data/ext_pruned_markers_maf10.txt (one marker per line)
##         results/stage1_external_summary.txt

suppressPackageStartupMessages({ library(data.table); library(SNPRelate); library(parallel) })
devtools::load_all("~/gitlab/LDscnR/", quiet = TRUE)
OUT <- "module_external_validation"
SLIDE <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) as.integer(a[1]) else 100L }
TAG <- if (SLIDE == 400L) "" else paste0("_slide", SLIDE)

d <- readRDS(file.path(OUT, "data/ext_geno_maf10.rds")); map <- d$map
ld_decay <- readRDS(file.path(OUT, paste0("data/ld_decay_ext_maf10", TAG, ".rds")))
gds_path <- file.path(OUT, paste0("data/ext_maf10", TAG, ".gds"))
stopifnot(file.exists(gds_path), identical(ld_decay$params$slide, SLIDE))

t0 <- Sys.time()
s1 <- ld_complexity_reduction(map = map, LD_decay = ld_decay, rho = 0.5, cores = 8, gds = gds_path)
message("Stage 1 elapsed: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
saveRDS(s1, file.path(OUT, "data/ext_stage1_maf10_rho05.rds"))
writeLines(s1$pruned, file.path(OUT, "data/ext_pruned_markers_maf10.txt"))

cl <- as.data.table(s1$clusters)
orig <- as.data.table(readRDS("module0_ld_pruning/data/pruned_stage1.rds")$clusters)
fmt <- function(x, n_in) sprintf("markers in = %d; clusters = %d (%.1f%% of markers); singletons = %d; >=5 SNPs = %d; max cluster = %d; median size = %.0f; mean size = %.2f",
                                 n_in, nrow(x), 100 * nrow(x) / n_in, sum(x$n_snps == 1), sum(x$n_snps >= 5), max(x$n_snps), median(x$n_snps), mean(x$n_snps))
sink(file.path(OUT, "results/stage1_external_summary.txt"), split = TRUE)
cat("EXTERNAL (21 inds, MAF>0.1, rho=0.5, slide", SLIDE, "):\n ", fmt(cl, nrow(map)), "\n")
cat("ORIGINAL (165 inds, MAF>=0.05, rho=0.5):\n ", fmt(orig, sum(orig$n_snps)), "\n\n")
cat("Clusters per chromosome (external):\n"); print(cl[, .(n_clusters = .N, n_markers = sum(n_snps), max_size = max(n_snps)), by = Chr][order(as.integer(sub("Chr", "", Chr)))])
cat("\nCluster-size distribution (external):\n"); print(table(cut(cl$n_snps, c(0, 1, 2, 4, 9, 19, 49, Inf), labels = c("1", "2", "3-4", "5-9", "10-19", "20-49", "50+"))))
cat("\nPhysical span of clusters with >=2 SNPs, kb (external vs original):\n")
span <- function(x) { m <- x[n_snps >= 2, members]; sapply(m, function(z) { p <- as.integer(sub(".*:", "", z)); (max(p) - min(p)) / 1e3 }) }
print(rbind(external = quantile(span(cl), c(.5, .9, .99)), original = quantile(span(orig), c(.5, .9, .99))))
sink()
