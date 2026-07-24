## =============================================================================
## Module E, Script E2b -- Build SLiM-ready founder VCFs from the phased pool
## =============================================================================
##
## PURPOSE
## -------
## Turn the E1 phased haplotype pool (moduleE_founder_haplotypes.rds) into
## per-chromosome VCFs that the neutral SLiM model reads as founders via
## readHaplosomesFromVCF(..., mutationType=m10). One-time prep; the SLiM model
## does the per-deme DRAW of N founders internally (so we don't explode into
## reps x N files here).
##
## ONE COMBINED, UNCOMPRESSED VCF PER CHROMOSOME (validated design):
##   - 30 HAPLOID aquilonia samples (GT "0"/"1"), then 13 DIPLOID polyctena
##     samples (GT "a|b"). aq haplotypes seed founder MALES (one haplotype each);
##     pol diploids seed founder FEMALES (rows 2i-1,2i of the pool are a female's
##     two haplotypes, from E1's interleaving -> preserves real within-female phase
##     and thus within-species LD).
##   - Combined (not separate aq/pol) because SLiM's readHaplosomesFromVCF creates
##     a NEW mutation object per read call; reading aq and pol separately makes TWO
##     distinct mutations at each position. One combined read -> one shared
##     mutation per position (validated: freqs match the pool exactly).
##   - Uncompressed .vcf: readHaplosomesFromVCF cannot read .vcf.gz.
## The SLiM model reads all 30+13 founders, then samples N per deme internally.
##
## Restricted to DI > -25 markers (== the regenerated AIMs, E2a): AIM-AIM LD
## depends only on the recmap distance between the AIMs, so seeding only these is
## sufficient and keeps SLiM light. CHROM renamed Chr<c> -> ch<c> to match the
## SLiM chromosome symbols. POS kept 1-based (SLiM VCF I/O is 1-based; it places
## the mutation at POS-1, i.e. the 0-based AIM position).
##
## OUTPUT (into the branch repo):
##   <BRANCH>/SLiM/founders_DIm25/aq_founders_ch<c>.vcf(.gz)   30 haploid samples
##   <BRANCH>/SLiM/founders_DIm25/pol_founders_ch<c>.vcf(.gz)  13 diploid samples
##   founders_DIm25_manifest.txt
##
## Review, then run.
## =============================================================================

suppressMessages(library(data.table))

## ------------------------------- CONFIG -------------------------------------
FORMICA <- "/Users/petrikem/gitlab/formica_hybrid"
POOL    <- file.path(FORMICA, "data/moduleE_founders/moduleE_founder_haplotypes.rds")
BRANCH  <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution"

DI_THRESHOLD <- -25

## THINNING (for the memory-light founder-number SCREEN). Set MODULEE_THIN to a
## marker count (e.g. 15000) to write a random genome-wide subset of the DI>-25
## markers to founders_DIm25_thin<N>/. Empty = full 50,072 set (production null).
## LD-decay/pi/FST are well estimated from ~15k markers; the production sorting
## null uses the full set. The thinned set is a strict subset (seeded, reproducible).
THIN_N    <- suppressWarnings(as.integer(Sys.getenv("MODULEE_THIN", "")))
THIN_SEED <- 1L

## DI-STRATIFIED MODE (MODULEE_STRAT=<markers per DI bin>). Universe becomes the
## parent_maf >= 0.15 gate (650,950 markers, the locked Module A gate) instead of
## DI > -25, and markers are drawn EQUALLY from fixed DI bins spanning the whole
## DI range. Equal allocation is deliberate: uniform thinning would be swamped by
## the mid-DI bulk and starve both tails. Low-DI bins then serve as a demographic
## anchor (they cannot be driven by ancestry sorting) while the DI gradient gives
## a dose-response test.
STRAT_N   <- suppressWarnings(as.integer(Sys.getenv("MODULEE_STRAT", "")))
MAF_GATE  <- 0.15
DI_BREAKS <- c(-Inf,-90,-75,-60,-50,-40,-30,-25,-20,-15,Inf)

OUTDIR <- if (!is.na(STRAT_N)) file.path(BRANCH, sprintf("SLiM/founders_maf015_DIstrat%d", STRAT_N)) else
          if (!is.na(THIN_N))  file.path(BRANCH, sprintf("SLiM/founders_DIm25_thin%d", THIN_N)) else
                               file.path(BRANCH, "SLiM/founders_DIm25")

## smoke-test hook: e.g. MODULEE_TEST_CHROMS="Chr27" to build just those
TEST_CHROMS <- Sys.getenv("MODULEE_TEST_CHROMS", "")
## ---------------------------------------------------------------------------

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

ph  <- readRDS(POOL)
Haq <- ph$aqu                 # 30 hap x M (0/1)
Hpol<- ph$pol                 # 26 hap x M (0/1); rows 2i-1,2i = female i
map <- as.data.table(ph$map)  # marker, Chr, Pos, cM, DI, maf_hyb
stopifnot(nrow(Hpol) %% 2L == 0L)

cat(sprintf("pool: aq %d haplotypes, pol %d haplotypes (%d diploid females), %d markers\n",
            nrow(Haq), nrow(Hpol), nrow(Hpol) %/% 2L, ncol(Haq)))

if (!is.na(STRAT_N)) {
  ## ---- DI-STRATIFIED on the parent_maf >= 0.15 universe --------------------
  p_pool <- (colSums(Haq) + colSums(Hpol)) / (nrow(Haq) + nrow(Hpol))
  map[, parent_maf := pmin(p_pool, 1 - p_pool)]
  uni <- which(map$parent_maf >= MAF_GATE & is.finite(map$DI))
  cat(sprintf("parent_maf >= %.2f universe: %d markers\n", MAF_GATE, length(uni)))
  bins <- cut(map$DI[uni], DI_BREAKS)
  set.seed(THIN_SEED)
  sel <- sort(unlist(lapply(split(uni, bins), function(ix)
    if (length(ix) <= STRAT_N) ix else sample(ix, STRAT_N))))
  cat(sprintf("DI-STRATIFIED: %d markers from %d bins (<=%d per bin) -> %s\n",
              length(sel), nlevels(bins), STRAT_N, OUTDIR))
  print(data.table(DI_bin = levels(bins),
                   n = as.integer(table(cut(map$DI[sel], DI_BREAKS)))))
} else {
  ## ---- original: DI > -25, optionally randomly thinned ---------------------
  map[, keep := DI > DI_THRESHOLD]
  sel <- which(map$keep)
  cat(sprintf("DI > %d markers: %d / %d\n", DI_THRESHOLD, length(sel), nrow(map)))
  if (!is.na(THIN_N)) {
    if (THIN_N >= length(sel)) stop("MODULEE_THIN >= number of DI>-25 markers")
    set.seed(THIN_SEED)
    sel <- sort(sample(sel, THIN_N))
    cat(sprintf("THINNED to %d markers (seed %d) -> %s\n", length(sel), THIN_SEED, OUTDIR))
  }
}

chroms <- unique(map$Chr[sel])
if (nzchar(TEST_CHROMS)) chroms <- intersect(chroms, strsplit(TEST_CHROMS, ",")[[1]])

## ---- combined VCF writer (aq haploid samples first, then pol diploid) ----
N_AQ  <- nrow(Haq)              # 30 aq haplotypes -> 30 haploid samples (males)
N_POL <- nrow(Hpol) %/% 2L      # 13 pol diploid females
aq_samples  <- sprintf("aq_hap%02d",  seq_len(N_AQ))
pol_samples <- sprintf("pol_fem%02d", seq_len(N_POL))
all_samples <- c(aq_samples, pol_samples)

write_combined_vcf <- function(ch_slim, pos, ids, gaq, gpol, out_vcf) {
  ## gaq: markers x N_AQ (haploid "0"/"1"); gpol: markers x N_POL ("a|b")
  hdr <- c("##fileformat=VCFv4.2",
           sprintf("##contig=<ID=%s>", ch_slim),
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
           paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
                   all_samples), collapse = "\t"))
  fixed <- sprintf("%s\t%d\t%s\tA\tC\t.\tPASS\t.\tGT", ch_slim, pos, ids)
  gt    <- cbind(gaq, gpol)                                        # markers x 43
  body  <- paste(fixed, apply(gt, 1L, paste, collapse = "\t"), sep = "\t")
  writeLines(c(hdr, body), out_vcf)                                # UNCOMPRESSED
  out_vcf
}

manifest <- rbindlist(lapply(chroms, function(ch) {
  ci  <- sel[map$Chr[sel] == ch]
  ci  <- ci[order(map$Pos[ci])]
  pos <- map$Pos[ci]; ids <- map$marker[ci]
  ch_slim <- sub("Chr", "ch", ch)

  ## aquilonia haploid genotypes (single allele)
  gaq <- matrix(as.character(t(Haq[, ci, drop = FALSE])), nrow = length(ci))  # mk x 30
  ## polyctena diploid phased genotypes, pairing pool rows (2i-1, 2i)
  H1 <- Hpol[seq(1, nrow(Hpol), by = 2), ci, drop = FALSE]
  H2 <- Hpol[seq(2, nrow(Hpol), by = 2), ci, drop = FALSE]
  gpol <- matrix(paste(t(H1), t(H2), sep = "|"), nrow = length(ci))           # mk x 13

  f <- write_combined_vcf(ch_slim, pos, ids, gaq, gpol,
                          file.path(OUTDIR, sprintf("founders_%s.vcf", ch_slim)))
  data.table(Chr = ch, ch_slim = ch_slim, n_markers = length(ci),
             n_aq_hap = N_AQ, n_pol_dip = N_POL, vcf = basename(f))
}))

fwrite(manifest, file.path(OUTDIR, "founders_DIm25_manifest.txt"), sep = "\t")
cat("\ncombined founder VCFs written to ", OUTDIR, ":\n", sep = "")
print(manifest)
cat(sprintf("\nlayout per VCF: %d haploid aq samples (males) then %d diploid pol samples (females)\n",
            N_AQ, N_POL))
cat("SLiM: read all into c(p1_male_haplosomes, p2_female_haplosomes), then sample N per deme.\n")
