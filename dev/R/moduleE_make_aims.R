## =============================================================================
## Module E, Script E2a -- Regenerate SLiM AIMs from our DI > -25 markers
## =============================================================================
##
## PURPOSE
## -------
## The SLiM model (SpecIAnt_rufa_neutral_realfounders.slim) plants ancestry
## markers at "AIM" positions read from per-chromosome files. The AIMs shipped in
## the repo (AIMs/, from species_diagnostic_markers.DTOL_PS.10FaquFpol.txt,
## 52,819 markers) may be stale; we want the AIMs to be EXACTLY the SNPs with
## DiagnosticIndex > -25 in OUR current data (map_hyb_005), i.e. the 50,072-marker
## set that Modules A/B/C and the phased founder pool (E1) are built on. This keeps
## the null's marker set identical to the observed analysis (apples-to-apples).
##
## Coordinates are already verified consistent with the SLiM genome: all DI>-25
## markers fall within the per-chromosome lengths Ls, Chr23 is absent in our data
## (matching the model's exclusion of the social-supergene chromosome), and
## positions are 1-based -> written 0-based for SLiM (as AIMs_for_SLiM.R does).
##
## OUTPUT (into the branch repo, new dir so the originals are untouched):
##   <BRANCH>/SLiM/AIMs_DIm25/AIMs_ch<id>.txt   -- two cols, tab-sep, no header:
##        chr_id   pos_0based
##   + AIMs_DIm25_manifest.txt (per-chrom counts + provenance)
##
## Review, then run (fast, pure data prep).
## =============================================================================

suppressMessages(library(data.table))

## ------------------------------- CONFIG -------------------------------------
FORMICA   <- "/Users/petrikem/gitlab/formica_hybrid"
RDATA     <- file.path(FORMICA, "data/hybrids_and_parents_maf005.Rdata")

## output into the SLiM branch (persistent clone); new dir, does NOT clobber AIMs/
BRANCH    <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution"
OUTDIR    <- file.path(BRANCH, "SLiM/AIMs_DIm25")

DI_THRESHOLD <- -25

## SLiM genome: chromosome lengths indexed by chr id 1..27 (NA at 23 = excluded).
## Used only as a safety check that no marker sits past a chromosome end.
Ls <- c(16449812,18558483,15805178,16666556,14957794,10925039,13372598,13028660,
        10584038,13832892,11830133,11260715,11449430,7730741,10965718,11699404,
        11021207,8644910,7297899,9820556,10314430,7671390,NA,6269374,13114028,
        9969264,7047176)
## ---------------------------------------------------------------------------

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ---- load marker map ----
e <- new.env(); load(RDATA, envir = e)
map <- as.data.table(e$map_hyb_005)[, .(Chr, Pos, marker, DiagnosticIndex)]
map[, chr_id := as.integer(sub("Chr", "", Chr))]

## ---- select DI > threshold ----
aims <- map[DiagnosticIndex > DI_THRESHOLD]
setorder(aims, chr_id, Pos)
cat(sprintf("DI > %d markers: %d across %d chromosomes\n",
            DI_THRESHOLD, nrow(aims), length(unique(aims$chr_id))))

## ---- safety checks against the SLiM genome ----
aims[, L := Ls[chr_id]]
if (aims[, any(is.na(L))])
  stop("AIM on a chromosome with no SLiM length (Chr23 or unknown): ",
       paste(unique(aims[is.na(L), Chr]), collapse = ", "))
if (aims[, any(Pos > L)])
  stop("AIM position beyond chromosome end -- would be rejected/trimmed by SLiM")
if (aims[, min(Pos)] < 1L)
  stop("expected 1-based positions (min Pos < 1)")

## ---- write per-chromosome 0-based AIM files (repo format) ----
manifest <- aims[, {
  out <- data.table(chr_id = chr_id, pos_0based = Pos - 1L)   # 1-based -> 0-based
  f <- file.path(OUTDIR, sprintf("AIMs_ch%d.txt", chr_id[1]))
  fwrite(out, f, sep = "\t", col.names = FALSE)
  .(n = .N, min_pos0 = min(out$pos_0based), max_pos0 = max(out$pos_0based),
    L = L[1], file = basename(f))
}, by = chr_id][order(chr_id)]

fwrite(manifest, file.path(OUTDIR, "AIMs_DIm25_manifest.txt"), sep = "\t")
cat("\nper-chromosome AIM files written to ", OUTDIR, ":\n", sep = "")
print(manifest)
cat(sprintf("\nTOTAL AIMs written: %d\n", sum(manifest$n)))
cat("provenance: map_hyb_005 (hybrids_and_parents_maf005.Rdata), DI > ",
    DI_THRESHOLD, "\n", sep = "")
