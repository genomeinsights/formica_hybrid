## module_external_validation / 02: LD-decay estimation for the external hybrid
## samples (21 individuals), SNPs with MAF > 0.1 only.
##
## Purpose: the GRM for EMMAX needs an LD-pruned marker set built with the same
## Stage-1 clustering (ld_complexity_reduction) as the original hybrid data, and
## that clustering needs a decay fit made on the SAME marker set (LDscnR docs:
## decay parameters depend on marker density because `slide` is counted in SNPs).
##
## Settings mirror R/LD_decay_from_DIEM.R (n_win_decay = 100, slide = 400,
## ld_method = "corr"); seed fixed so the fit is reproducible. Edge lists are
## NOT kept (keep_el = FALSE): ld_complexity_reduction() can rebuild them on the
## fly from the GDS, which is faster than writing/reading them.
##
## CAVEAT built into the data: n = 21 individuals, so r^2 between unlinked SNPs
## has expectation ~1/n ~ 0.05 and a 95th percentile near 0.15-0.2 -- the
## background-LD estimate and the decay asymptote will be much higher than for
## the 165-individual original (see results/ld_decay_external_summary.txt).
##
## Run from the repo root:  Rscript module_external_validation/R/02_ld_decay_external.R
## Needs bcftools. Reads : data/HybridSamples_SNPQ30...FiSeDe.vcf.gz,
##                         module_external_validation/data/vcf_sites_maf21.tsv.gz (from 00_)
## Writes: data/ext_geno_maf10.rds   (list: geno [21 x M ALT-dosage], map)
##         data/ext_maf10.gds        (GDS used by the decay fit and, later, clustering)
##         data/ld_decay_ext_maf10.rds
##         results/ld_decay_external_summary.txt, results/ld_decay_external.png

suppressPackageStartupMessages({ library(data.table); library(SNPRelate); library(parallel); library(ggplot2) })
devtools::load_all("~/gitlab/LDscnR/", quiet = TRUE)
OUT <- "module_external_validation"
VCF <- "data/HybridSamples_SNPQ30.biall.fixedHeader.indDP.minDP8.hwe.376inds.AN10percMiss_FiSeDe.vcf.gz"
MAF_MIN <- 0.1
CORES <- 8L
## optional arg: slide (SNPs). Default 400 mirrors the original; slide_bp differs because the marker set is sparser.
SLIDE <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) as.integer(a[1]) else 400L }
TAG <- if (SLIDE == 400L) "" else paste0("_slide", SLIDE)
geno_f <- file.path(OUT, "data/ext_geno_maf10.rds")

## ---- 1. genotypes for MAF > 0.1 sites on the 26 chromosomes --------------------
if (!file.exists(geno_f)) {
  v <- fread(file.path(OUT, "data/vcf_sites_maf21.tsv.gz"), col.names = c("Chr", "Pos", "REF", "ALT", "MAF", "FMISS"))
  ## positions occurring in >1 record (split multiallelic sites) are dropped entirely:
  ## marker IDs must be unique and the two records are not independent SNPs
  v[, dupPos := duplicated(paste(Chr, Pos)) | duplicated(paste(Chr, Pos), fromLast = TRUE)]
  message("dropping ", sum(v$dupPos & v$MAF > MAF_MIN), " MAF>", MAF_MIN, " records at duplicated positions")
  v <- v[MAF > MAF_MIN & !dupPos & grepl("^Chr[0-9]+$", Chr)]
  tg <- tempfile(fileext = ".tsv"); fwrite(v[, .(chr = sub("^Chr", "chromosome_", Chr), Pos)], tg, sep = "\t", col.names = FALSE)
  gt <- tempfile(fileext = ".tsv")
  stopifnot(system(sprintf("bcftools query -T %s -f '%%CHROM\\t%%POS\\t[%%GT\\t]\\n' %s 2>/dev/null > %s", tg, VCF, gt)) == 0)
  g <- fread(gt, header = FALSE, sep = "\t"); g[, V1 := sub("^chromosome_", "Chr", V1)]
  samp <- system(sprintf("bcftools query -l %s 2>/dev/null", VCF), intern = TRUE)
  ns <- length(samp); g <- g[, 1:(ns + 2)]
  stopifnot(nrow(g) == nrow(v), all(g$V1 == v$Chr), all(g$V2 == v$Pos))
  gm <- as.matrix(g[, 3:(ns + 2)])
  dos <- matrix(NA_integer_, nrow(gm), ns)
  dos[gm %in% c("0/0", "0|0")] <- 0L; dos[gm %in% c("0/1", "1/0", "0|1", "1|0")] <- 1L; dos[gm %in% c("1/1", "1|1")] <- 2L
  geno <- t(dos); colnames(geno) <- paste0(v$Chr, ":", v$Pos); rownames(geno) <- samp
  saveRDS(list(geno = geno, map = v[, .(Chr, Pos, marker = paste0(Chr, ":", Pos), REF, ALT, maf21 = MAF, fmiss21 = FMISS)]), geno_f)
}
d <- readRDS(geno_f); geno <- d$geno; map <- d$map
stopifnot(identical(colnames(geno), map$marker))
## SNPRelate needs chromosomes in order and positions sorted within chromosome
map[, chr_n := as.integer(sub("Chr", "", Chr))]; o <- order(map$chr_n, map$Pos); stopifnot(identical(o, seq_len(nrow(map))))

## ---- 2. GDS + LD decay -----------------------------------------------------------
gds_f <- file.path(OUT, paste0("data/ext_maf10", TAG, ".gds"))
if (file.exists(gds_f)) file.remove(gds_f)
gds <- create_gds_from_geno(geno, map, gds_f)
t0 <- Sys.time()
ld_decay <- compute_LD_decay(gds, n_win_decay = 100, slide = SLIDE, ld_method = "corr",
                             keep_el = FALSE, cores = CORES, min_maf_decay = MAF_MIN, seed = 1)
message("LD decay elapsed: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
saveRDS(ld_decay, file.path(OUT, paste0("data/ld_decay_ext_maf10", TAG, ".rds")))

## ---- 3. summary ----------------------------------------------------------------------
sink(file.path(OUT, paste0("results/ld_decay_external_summary", TAG, ".txt")), split = TRUE)
cat("External hybrids: ", nrow(geno), "individuals,", ncol(geno), "SNPs with MAF >", MAF_MIN, "on", uniqueN(map$Chr), "chromosomes\n")
cat("missing genotype fraction:", round(mean(is.na(geno)), 4), "\n\n")
cat("Background LD (r^2) estimate:\n"); print(ld_decay$params$background_ld %||% ld_decay$params)
cat("\nPer-chromosome decay summary:\n"); print(as.data.frame(ld_decay$decay_sum), digits = 4)
cat("\nRecommendation:\n"); print(ld_decay$recommendation)
orig <- readRDS("module0_ld_pruning/data/ld_decay_DIEM_100w.rds")
cat("\n--- original hybrid data (165 inds, 1.1M SNPs) for comparison ---\n"); print(as.data.frame(orig$decay_sum), digits = 4)
sink()
p <- plot(ld_decay)
if (inherits(p, "ggplot")) ggsave(file.path(OUT, paste0("results/ld_decay_external", TAG, ".png")), p, width = 9, height = 6, dpi = 200)
snpgdsClose(gds)
