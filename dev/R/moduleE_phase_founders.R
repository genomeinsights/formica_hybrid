## =============================================================================
## Module E, Script E1 -- Phase + impute the parental founder haplotypes
## =============================================================================
##
## PURPOSE
## -------
## Turn the 15 (diploid, female) parental samples per species into a pool of
## real, phased, gap-filled haplotypes that will SEED the SLiM neutral null
## (Module E). Conditioning founders on real haplotypes -- rather than
## independent binomial draws at empirical frequencies -- is what preserves
## within-species LD, which the parent-LD diagnostic showed is 60-85% of
## admixture-scale LD in low-recombination regions and is threshold-proof.
##
## We phase + impute EACH SPECIES SEPARATELY. Pooling the two species would let
## Beagle borrow haplotype information across the species boundary and
## manufacture artificial intermediate haplotypes -- corrupting exactly the
## within-species structure we are trying to keep. Imputation is LD-based, so it
## mildly reinforces the common-haplotype structure; that is acceptable (and far
## better than fixed differences), and the downstream LD-decay calibration is
## where any over-cleaning would show up.
##
## Phasing is informed by the empirical genetic (cM) map, so phase confidence is
## highest exactly where surviving (low-recombination) LD matters most.
##
## SCOPE: all 1,114,340 markers are phased/imputed here, so the marker-set
## decision (50k DI>-25 vs the maf>=0.15 set vs a thinned genome) is deferred to
## the pre-SLiM thinning step (Script E2/E3) and the SAME phased pool serves
## every choice.
##
## OUTPUTS (under data/moduleE_founders/)
##   vcf_phased/<species>_<Chr>.phased.vcf.gz(+.tbi)  -- per-chrom SLiM inputs
##   moduleE_founder_haplotypes.rds                   -- genome-wide phased pool:
##        list(aqu = <2*n_aqu x M 0/1 haplotype matrix>,
##             pol = <2*n_pol x M 0/1 haplotype matrix>,
##             map = <data.table: marker, Chr, Pos, cM, DiagnosticIndex, ...>,
##             dropped_samples = <char>, params = <list>)
##   moduleE_phase_QC.txt                             -- human-readable QC report
##
## NB. NOTHING here is calibrated on the sorting statistics -- this is pure
## data preparation. Run it, review the QC, THEN proceed to founder-number
## screening (E2).
##
## Author: (Module E build) -- review before running.
## =============================================================================

suppressMessages({
  library(data.table)
})

## ------------------------------- CONFIG -------------------------------------
FORMICA   <- "/Users/petrikem/gitlab/formica_hybrid"
RDATA     <- file.path(FORMICA, "data/hybrids_and_parents_maf005.Rdata")
RECMAP    <- file.path(FORMICA, "data/Frufa_DTOL_PR.ref_genome.recmap")

JAVA      <- "/opt/homebrew/opt/openjdk/bin/java"
BEAGLE    <- path.expand("~/software/beagle/beagle.27Feb25.75f.jar")
BCFTOOLS  <- "bcftools"      # on PATH (/opt/homebrew/bin)
BGZIP     <- "bgzip"
TABIX     <- "tabix"

OUTDIR    <- file.path(FORMICA, "data/moduleE_founders")
VCF_IN    <- file.path(OUTDIR, "vcf_in")       # unphased input VCFs
VCF_OUT   <- file.path(OUTDIR, "vcf_phased")   # Beagle output
MAPDIR    <- file.path(OUTDIR, "beagle_maps")  # per-chrom cM maps for Beagle

SPECIES <- c(aqu = "aquilonia_parent", pol = "polyctena_parent")

## smoke-test hook: export MODULEE_TEST_CHROMS="Chr27" (comma-sep) to run the
## full pipeline on just those chromosomes. Empty = all chromosomes.
TEST_CHROMS <- Sys.getenv("MODULEE_TEST_CHROMS", "")

## an individual with more missingness than this is dropped before phasing
MAX_IND_MISS   <- 0.40
## a marker all-missing (or, optionally, too sparse) within a species is dropped
MAX_MARK_MISS  <- 0.90     # per-species, applied per marker; 1.0 = only drop all-NA

## Beagle requires STRICTLY increasing cM in its map file, but interpolation
## gives many markers identical cM in flat/zero-recombination stretches and in
## the chromosome tail beyond the last recmap control point. We break ties with
## a negligible per-marker ramp (total added across a whole chromosome ~= N*EPS,
## e.g. 40k * 1e-7 = 0.004 cM -- genetically nil).
MAP_EPS_CM     <- 1e-7

## compute budget (tuned for an 8-core / 30 GB envelope on an 11-core M4)
## Each job is tiny (<=15 samples, ~40k markers/chrom), so throughput is best
## with many single-threaded jobs rather than few multi-threaded ones.
N_JOBS         <- 8L       # how many Beagle jobs (species x chrom) at once -> 8 cores
BEAGLE_THREADS <- 1L       # threads per Beagle job (N_JOBS*BEAGLE_THREADS = 8)
BEAGLE_XMX     <- "3g"     # JVM heap per job; 8 * 3g = 24 GB peak (< 30 GB)

## Beagle: n=15 statistical phasing, no external reference panel. Defaults are
## fine for phasing; we pass the cM map. 'ne' is left at Beagle's default.
## ---------------------------------------------------------------------------

dir.create(VCF_IN,  recursive = TRUE, showWarnings = FALSE)
dir.create(VCF_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(MAPDIR,  recursive = TRUE, showWarnings = FALSE)

QC <- c()
qc <- function(...) { line <- paste0(...); message(line); QC[[length(QC) + 1L]] <<- line }

qc("== Module E / E1 phase+impute founders ==  ", as.character(Sys.time()))

## ------------------------------ 1. LOAD -------------------------------------
qc("\n[1] loading ", basename(RDATA))
e <- new.env(); load(RDATA, envir = e)
GTs <- e$GTs_with_parents          # samples x markers, dosage 0/1/2, NA missing
sd  <- as.data.table(e$sample_data_with_parents)
map <- as.data.table(e$map_hyb_005)
stopifnot(nrow(GTs) == nrow(sd), ncol(GTs) == nrow(map))
qc("    GTs: ", nrow(GTs), " samples x ", ncol(GTs), " markers")

## marker table we carry through (order == columns of GTs)
map[, idx := .I]
setnames(map, old = intersect(c("DiagnosticIndex"), names(map)), new = "DI", skip_absent = TRUE)

## ------------------------- 2. cM PER MARKER ---------------------------------
## Interpolate cM at each marker position from the recombination map (same
## linear interpolation used across the pipeline; rule=2 extrapolates flat at
## chromosome ends). This drives both the Beagle map and the marker table.
qc("\n[2] interpolating cM per marker from ", basename(RECMAP))
rec <- fread(RECMAP)                                  # chr,pos,cM,cM/Mb
rec[, Chr := sub("chromosome_", "Chr", chr)]
map[, cM := NA_real_]
for (ch in unique(map$Chr)) {
  r <- rec[Chr == ch]
  if (nrow(r) < 2L) { qc("    WARNING: <2 recmap rows for ", ch, " -- cM left NA"); next }
  j <- map[Chr == ch, which = TRUE]
  map[j, cM := approx(r$pos, r$cM, xout = map$Pos[j], rule = 2)$y]
}
qc("    markers with cM: ", map[!is.na(cM), .N], " / ", nrow(map))

## ------------------- 3. helper: dosage -> VCF (one chrom, one species) -------
## Writes a minimal but valid biallelic VCF. Dosage counts the ALT allele:
##   0 -> 0/0, 1 -> 0/1, 2 -> 1/1, NA -> ./.  (REF=A, ALT=C, arbitrary but
## consistent, so phased output maps back to the same dosage coding).
dose_to_gt <- function(d) {
  ## d is an integer vector 0/1/2/NA -> character GT
  out <- c("0/0", "0/1", "1/1")[d + 1L]
  out[is.na(d)] <- "./."
  out
}

write_species_chrom_vcf <- function(sp_key, ch, rows, mark_idx) {
  ## rows      : sample row indices (this species, after dropping high-miss inds)
  ## mark_idx  : marker column indices on this chromosome (already filtered)
  samp_ids <- sd$Sample_ID[rows]
  sub <- GTs[rows, mark_idx, drop = FALSE]            # inds x markers (dosage)
  pos <- map$Pos[mark_idx]
  o   <- order(pos)                                    # VCF must be position-sorted
  sub <- sub[, o, drop = FALSE]; pos <- pos[o]; mk <- map$marker[mark_idx][o]

  ## GT matrix (markers x inds). t(sub) is markers x inds; as.integer() flattens
  ## it column-major, so reshape with the DEFAULT byrow=FALSE to reconstruct the
  ## same layout (byrow=TRUE would transpose genotypes onto the wrong markers).
  gt <- matrix(dose_to_gt(as.integer(t(sub))), nrow = length(pos))

  vcf <- data.table(
    `#CHROM` = ch, POS = pos, ID = mk, REF = "A", ALT = "C",
    QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT"
  )
  vcf <- cbind(vcf, as.data.table(gt))
  setnames(vcf, (ncol(vcf) - length(samp_ids) + 1L):ncol(vcf), samp_ids)

  out_vcf <- file.path(VCF_IN, sprintf("%s_%s.vcf", sp_key, ch))
  hdr <- c(
    "##fileformat=VCFv4.2",
    sprintf("##contig=<ID=%s>", ch),
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
  )
  writeLines(hdr, out_vcf)
  fwrite(vcf, out_vcf, sep = "\t", quote = FALSE, append = TRUE, col.names = TRUE)
  ## bgzip + index
  system2(BGZIP, c("-f", shQuote(out_vcf)))
  system2(TABIX, c("-f", "-p", "vcf", shQuote(paste0(out_vcf, ".gz"))))
  paste0(out_vcf, ".gz")
}

## ---------------------- 4. per-chrom Beagle cM map --------------------------
## Beagle map format: 4 cols, no header: CHROM  ID  cM  bp
## NB. Written PER SPECIES (not just per chrom): each species filters markers by
## its own missingness, so the maps differ; and a shared per-chrom path would be
## written concurrently by the two species jobs under mclapply -> race condition
## that corrupts the file (a truncated/duplicated line -> Beagle map error).
write_beagle_map <- function(sp_key, ch, mark_idx) {
  m <- map[mark_idx][order(Pos)]
  bm <- data.table(Chr = ch, ID = m$marker, cM = m$cM, bp = m$Pos)
  ## Beagle requires finite, non-decreasing cM; guard against any NA/dupes
  bm <- bm[is.finite(cM)]
  ## strictly increasing: cummax removes any decrease, ramp breaks ties.
  ## (non-decreasing x) + (strictly increasing ramp) is strictly increasing.
  bm[, cM := cummax(cM) + (seq_len(.N) - 1L) * MAP_EPS_CM]
  f <- file.path(MAPDIR, sprintf("%s_%s.map", sp_key, ch))
  fwrite(bm, f, sep = "\t", col.names = FALSE)
  f
}

## ------------------- 5. drop high-missing individuals -----------------------
qc("\n[5] per-species individuals (dropping missingness > ", MAX_IND_MISS, ")")
sp_rows <- list(); dropped <- c()
for (k in names(SPECIES)) {
  rows <- which(sd$Population == SPECIES[[k]])
  miss <- rowMeans(is.na(GTs[rows, , drop = FALSE]))
  keep <- rows[miss <= MAX_IND_MISS]
  drp  <- setdiff(rows, keep)
  if (length(drp))
    dropped <- c(dropped, sd$Sample_ID[drp])
  sp_rows[[k]] <- keep
  qc("    ", k, " (", SPECIES[[k]], "): kept ", length(keep), " / ", length(rows),
     if (length(drp)) paste0("  DROPPED: ", paste(sd$Sample_ID[drp], collapse = ", ")) else "")
}

## ------------------- 6. build inputs + run Beagle ---------------------------
## One job per (species, chromosome). Marker filtering is per-species: drop
## markers all-/mostly-missing WITHIN that species (Beagle needs data present).
chroms <- unique(map$Chr)
if (nzchar(TEST_CHROMS)) {
  chroms <- intersect(chroms, strsplit(TEST_CHROMS, ",")[[1]])
  qc("    [TEST MODE] restricting to chromosomes: ", paste(chroms, collapse = ", "))
}
jobs <- CJ(sp = names(SPECIES), ch = chroms, sorted = FALSE)
qc("\n[6] building VCFs + maps and phasing (", nrow(jobs), " species x chrom jobs)")

run_one <- function(i) {
  sp <- jobs$sp[i]; ch <- jobs$ch[i]
  rows <- sp_rows[[sp]]
  mk_all <- map[Chr == ch, which = TRUE]
  ## per-species marker missingness filter
  mmiss <- colMeans(is.na(GTs[rows, mk_all, drop = FALSE]))
  mk <- mk_all[mmiss <= MAX_MARK_MISS]
  if (length(mk) < 2L) return(data.table(sp = sp, ch = ch, n_in = length(mk),
                                          n_out = 0L, status = "skipped_too_few"))
  vcf_in  <- write_species_chrom_vcf(sp, ch, rows, mk)
  bmap    <- write_beagle_map(sp, ch, mk)
  out_pre <- file.path(VCF_OUT, sprintf("%s_%s.phased", sp, ch))
  ## delete any stale output first, so a Beagle failure can't be misread as
  ## success by the file.exists() check below.
  unlink(paste0(out_pre, ".vcf.gz")); unlink(paste0(out_pre, ".vcf.gz.tbi"))
  args <- c(sprintf("-Xmx%s", BEAGLE_XMX), "-jar", BEAGLE,
            sprintf("gt=%s", vcf_in),
            sprintf("map=%s", bmap),
            sprintf("out=%s", out_pre),
            sprintf("nthreads=%d", BEAGLE_THREADS))
  log <- suppressWarnings(system2(JAVA, args, stdout = TRUE, stderr = TRUE))
  writeLines(log, paste0(out_pre, ".log"))
  ok <- file.exists(paste0(out_pre, ".vcf.gz"))
  if (ok) system2(TABIX, c("-f", "-p", "vcf", shQuote(paste0(out_pre, ".vcf.gz"))))
  data.table(sp = sp, ch = ch, n_in = length(mk),
             n_out = if (ok) length(mk) else 0L,
             status = if (ok) "ok" else "FAILED")
}

if (requireNamespace("parallel", quietly = TRUE) && N_JOBS > 1L) {
  res <- data.table::rbindlist(parallel::mclapply(seq_len(nrow(jobs)), run_one,
                                                  mc.cores = N_JOBS))
} else {
  res <- data.table::rbindlist(lapply(seq_len(nrow(jobs)), run_one))
}
qc("\n[6] phasing job status:")
qc(paste(capture.output(print(res[, .N, by = status])), collapse = "\n"))
if (any(res$status == "FAILED")) qc("    !! some jobs FAILED -- inspect .log files in ", VCF_OUT)

## ------------- 7. read phased haplotypes back into a genome-wide pool --------
## Extract phased GTs per chrom/species with bcftools, split "a|b" into two
## haplotype rows per individual. Assemble a 0/1 haplotype matrix per species,
## aligned to a single genome-wide marker order.
qc("\n[7] reassembling genome-wide phased haplotype pool")

read_phased <- function(sp, ch) {
  f <- file.path(VCF_OUT, sprintf("%s_%s.phased.vcf.gz", sp, ch))
  if (!file.exists(f)) return(NULL)
  ## marker | ind1_h1 ind1_h2 ind2_h1 ... via bcftools
  txt <- system2(BCFTOOLS, c("query", "-f", shQuote("%ID[\\t%GT]\\n"), shQuote(f)),
                 stdout = TRUE)
  if (!length(txt)) return(NULL)
  ## force tab sep: GT values contain '|', which fread would otherwise
  ## auto-detect as the field separator and split on the wrong character.
  dt <- fread(text = paste(txt, collapse = "\n"), header = FALSE, sep = "\t")
  ids <- dt[[1]]; g <- as.matrix(dt[, -1])       # markers x inds ("a|b")
  ## split each "a|b" into two integer haplotypes -> markers x (2*inds)
  h1 <- matrix(as.integer(sub("\\|.*$", "", g)), nrow = nrow(g))
  h2 <- matrix(as.integer(sub("^.*\\|", "", g)), nrow = nrow(g))
  H <- matrix(0L, nrow = nrow(g), ncol = 2L * ncol(g))
  H[, seq(1, by = 2, length.out = ncol(g))] <- h1
  H[, seq(2, by = 2, length.out = ncol(g))] <- h2
  list(ids = ids, H = H)                          # H: markers x haplotypes
}

assemble_species <- function(sp) {
  pieces <- lapply(chroms, function(ch) read_phased(sp, ch))
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  ids <- unlist(lapply(pieces, `[[`, "ids"))
  Hm  <- do.call(rbind, lapply(pieces, `[[`, "H"))   # markers x haplotypes
  list(ids = ids, H = t(Hm))                          # -> haplotypes x markers
}

pool <- lapply(names(SPECIES), assemble_species)
names(pool) <- names(SPECIES)

## Align both species to a common marker order (intersection, in map order)
common <- Reduce(intersect, lapply(pool, `[[`, "ids"))
common <- map$marker[map$marker %in% common]          # keep genome order
qc("    common phased markers across species: ", length(common), " / ", nrow(map))

align <- function(p) {
  colnames(p$H) <- p$ids
  p$H[, common, drop = FALSE]
}
Haqu <- align(pool$aqu); Hpol <- align(pool$pol)
map_common <- map[match(common, marker),
                  .(marker, Chr, Pos, cM,
                    DI = if ("DI" %in% names(map)) DI else NA_real_,
                    maf_hyb = if ("maf_hyb" %in% names(map)) maf_hyb else NA_real_)]

out <- list(
  aqu = Haqu, pol = Hpol, map = map_common,
  dropped_samples = dropped,
  params = list(MAX_IND_MISS = MAX_IND_MISS, MAX_MARK_MISS = MAX_MARK_MISS,
                beagle = basename(BEAGLE), when = as.character(Sys.time()))
)
saveRDS(out, file.path(OUTDIR, "moduleE_founder_haplotypes.rds"))
qc("    saved ", file.path(OUTDIR, "moduleE_founder_haplotypes.rds"),
   "  (aqu ", nrow(Haqu), " hap x ", ncol(Haqu), " mk ; pol ", nrow(Hpol), " hap)")

## ------------------------- 8. founding-pool QC ------------------------------
## Quick diversity check on the PHASED pool (pre-SLiM): within-species allele
## freq / heterozygosity and between-species delta, to confirm the pool looks
## like the raw dosage characterisation and imputation did not distort it.
qc("\n[8] founding-pool QC (phased)")
p_aqu <- colMeans(Haqu); p_pol <- colMeans(Hpol)
He <- function(p) mean(2 * p * (1 - p))
qc("    within-aqu He = ", round(He(p_aqu), 4), " ; within-pol He = ", round(He(p_pol), 4))
qc("    mean |dp| between species = ", round(mean(abs(p_aqu - p_pol)), 3))
qc("    frac markers |dp|>0.5 = ", round(mean(abs(p_aqu - p_pol) > 0.5), 3))

writeLines(unlist(QC), file.path(OUTDIR, "moduleE_phase_QC.txt"))
qc("\n[done] QC written to ", file.path(OUTDIR, "moduleE_phase_QC.txt"))
