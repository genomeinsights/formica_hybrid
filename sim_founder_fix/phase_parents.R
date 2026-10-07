## =========================================================================
## sim_founder_fix -- statistically phase the 30 empirical parents (Beagle 5)
##
## Why: the genotype-mosaic founders (make_mosaic_founders.R) split heterozygous
## sites at random between haplosomes, which destroys within-species LD among
## polymorphisms shared by the two species (near-neutral loci lost ~95% of their
## short-range LD in the first simulations). Phased parental haplotypes keep it.
##
## Phasing set: all 195 samples (30 parents + 165 hybrids; the hybrids add haplotypes
## and are related within colonies, which helps phasing) on all SNPs of the full
## MAF >= 0.05 matrix (dense markers phase better than the panel alone), with the DI25
## panel SNPs taken from the DI25 inputs. Sporadic missing genotypes are imputed by Beagle.
## Genetic map: the reference recombination map (Frufa_DTOL_PR.ref_genome.recmap).
##
## Output: <OUT_DIR>/parents_phased.rds = list(map = panel map (load_panel order,
## DI25 + near-neutral), H = 60 x n_panel haplotype matrix (0/1, allele 1 = coded
## allele of the empirical genotype matrices), rows "<parent>_h1"/"_h2")
## plus the Beagle input/output files.
##
## Run from the formica_hybrid repo root:
##   Rscript sim_founder_fix/phase_parents.R <OUT_DIR> <BEAGLE_JAR> [JAVA=java] [THREADS=8]
## =========================================================================
source("sim_founder_fix/parent_ld_lib.R")
args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- args[1]; JAR <- args[2]
JAVA <- if (length(args) >= 3) args[3] else "java"
THREADS <- if (length(args) >= 4) as.integer(args[4]) else 8L
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
options(scipen = 999)

panel <- load_panel(include_neutral = TRUE)
inp <- readRDS("module_di25/data/di25_inputs.rds")
e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
Gall <- e$GTs_with_parents; m <- as.data.table(e$map_hyb_005)[, .(Chr, Pos = as.integer(Pos), marker)]; rm(e); invisible(gc())
smp <- c(rownames(inp$GTs_par), rownames(inp$GTs_hyb))
stopifnot(setequal(smp, rownames(Gall)))
Gall <- Gall[smp, , drop = FALSE]

## Panel SNPs keep the coding the founders use: DI25 SNPs from the DI25 inputs, near-
## neutral SNPs from the full matrix (as in load_panel()). The two sources code the
## allele in opposite directions at most shared SNPs, so DI25 SNPs are NOT taken from
## the full matrix; its remaining SNPs only add phasing information.
G25 <- rbind(inp$GTs_par, inp$GTs_hyb)[smp, , drop = FALSE]
m25 <- as.data.table(inp$map)[, .(Chr, Pos = as.integer(Pos), marker)]
stopifnot(identical(m25$marker, colnames(G25)))
sup <- which(!(m$marker %in% m25$marker))
G <- cbind(G25, Gall[, sup, drop = FALSE]); rm(Gall, G25); invisible(gc())
M <- rbind(m25, m[sup]); stopifnot(identical(M$marker, colnames(G)))
M[, chrn := as.integer(sub("Chr", "", Chr))]
keep <- M[, .I[!duplicated(paste(Chr, Pos))]]   # unique positions for Beagle; DI25 SNPs come first, so they win
o <- keep[order(M$chrn[keep], M$Pos[keep])]
M <- M[o]; G <- G[, o, drop = FALSE]
stopifnot(all(panel$map$marker %in% M$marker))
cat(sprintf("[phase] %d samples x %d SNPs (%d panel SNPs)\n", nrow(G), ncol(G), nrow(panel$map)))

## ---- VCF (allele 1 = coded allele -> ALT) ----
vcf <- file.path(OUT_DIR, "phase_in.vcf.gz")
gt <- matrix(c("0/0", "0/1", "1/1")[G + 1L], nrow(G)); gt[is.na(gt)] <- "./."
con <- gzfile(vcf, "w")
writeLines(c("##fileformat=VCFv4.2", sprintf("##contig=<ID=%s>", unique(M$Chr)),
             "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
             paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", smp), collapse = "\t")), con)
step <- 20000L
for (s in seq(1L, ncol(G), by = step)) {
  ix <- s:min(ncol(G), s + step - 1L)
  writeLines(paste(sprintf("%s\t%d\t%s\tA\tC\t.\tPASS\t.\tGT", M$Chr[ix], M$Pos[ix], M$marker[ix]),
                   apply(gt[, ix, drop = FALSE], 2, paste, collapse = "\t"), sep = "\t"), con)
}
close(con); rm(gt); invisible(gc())

## ---- PLINK-format genetic map from the reference recombination map ----
gm <- fread("data/Frufa_DTOL_PR.ref_genome.recmap"); gm[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
gm <- gm[Chr %in% M$Chr & !is.na(cM)][order(as.integer(sub("Chr", "", Chr)), pos)]
gm <- gm[, .SD[!duplicated(pos)][order(pos)][, cM := cummax(cM)], by = Chr]  # map must be non-decreasing
fwrite(gm[, .(Chr, ".", cM, pos)], file.path(OUT_DIR, "phase.map"), sep = " ", col.names = FALSE)

## ---- Beagle ----
outp <- file.path(OUT_DIR, "phase_out")
cmd <- sprintf("%s -Xmx24g -jar %s gt=%s map=%s out=%s nthreads=%d seed=1", JAVA, JAR, vcf,
               file.path(OUT_DIR, "phase.map"), outp, THREADS)
cat("[phase] ", cmd, "\n", sep = ""); stopifnot(system(cmd) == 0)

## ---- parental haplotypes at the panel SNPs ----
ph <- fread(cmd = sprintf("gzip -dc %s.vcf.gz | grep -v '^##'", outp), sep = "\t", colClasses = "character")
ph <- ph[match(panel$map$marker, ph$ID)]; stopifnot(!anyNA(ph$ID))
par <- rownames(inp$GTs_par)
H <- do.call(rbind, lapply(par, function(s) { x <- ph[[s]]
  rbind(as.integer(substr(x, 1, 1)), as.integer(substr(x, 3, 3))) }))
rownames(H) <- as.vector(rbind(paste0(par, "_h1"), paste0(par, "_h2"))); colnames(H) <- panel$map$marker
## the phased dosage must reproduce every observed parental genotype
Dp <- H[c(TRUE, FALSE), ] + H[c(FALSE, TRUE), ]
Pobs <- panel$parents[par, panel$map$marker]
cat(sprintf("[phase] phased vs observed parental genotypes: %.5f identical; %.2f%% imputed\n",
            mean(Dp == Pobs, na.rm = TRUE), 100 * mean(is.na(Pobs))))
saveRDS(list(map = panel$map, H = H), file.path(OUT_DIR, "parents_phased.rds"))
cat("[phase] wrote", file.path(OUT_DIR, "parents_phased.rds"), "\n")
