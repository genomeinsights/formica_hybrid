## =============================================================================
## Module E -- INFLATED DI-stratified founder VCFs (for the founder-number / shallow-
## bottleneck sweep)
## =============================================================================
## Same DI-stratified marker draw as make_founders.R (identical DI_BREAKS/MAF_GATE/
## SEED, so the marker SET matches the K-sweep founders for a given n_per_bin), but
## each real founder haplotype is written `copies` times. The pool becomes
## 30*copies aq (haploid) + 13*copies pol (diploid). The SLiM model then subsamples
## N_AQ/N_POL WITHOUT replacement from this inflated pool == a resample WITH
## replacement from the 30/13 real founders. So N_AQ/N_POL can exceed 30/13 to make
## the founding bottleneck shallow (the F_ST / drift-LD knob); genetic diversity is
## still capped at the real 30/13 haplotypes (copies add census, not variants).
##
## Usage:  Rscript make_founders_inflated.R <n_per_bin> <copies>
##   e.g.  Rscript make_founders_inflated.R 1500 40   -> ~15k markers, pool 1200/520
##         -> ../founders/maf015_DIstrat1500_x40/
## =============================================================================
suppressMessages(library(data.table))
.args <- commandArgs(trailingOnly = FALSE)
.self <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT  <- normalizePath(file.path(dirname(.self), ".."))
ca <- commandArgs(trailingOnly = TRUE)
N <- as.integer(ca[1]); C <- as.integer(ca[2])
if (is.na(N)) stop("give n_per_bin, e.g. 1500")
if (is.na(C) || C < 1) stop("give copies (>=1), e.g. 40")

DI_BREAKS <- c(-Inf,-90,-75,-60,-50,-40,-30,-25,-20,-15,Inf); MAF_GATE <- 0.15; SEED <- 1L
ph  <- readRDS(file.path(ROOT, "inputs/moduleE_founder_haplotypes.rds"))
map <- as.data.table(ph$map)
Haq <- ph$aqu; Hpol <- ph$pol
p <- (colSums(Haq)+colSums(Hpol))/(nrow(Haq)+nrow(Hpol)); map[, parent_maf := pmin(p,1-p)]
uni <- which(map$parent_maf >= MAF_GATE & is.finite(map$DI))
set.seed(SEED)                                             # identical selection to make_founders.R
sel <- sort(unlist(lapply(split(uni, cut(map$DI[uni], DI_BREAKS)),
                          function(ix) if (length(ix)<=N) ix else sample(ix, N))))

nAq  <- nrow(Haq); nPol <- nrow(Hpol) %/% 2L               # 30 aq haplotypes, 13 pol females
OUT <- file.path(ROOT, sprintf("founders/maf015_DIstrat%d_x%d", N, C))
dir.create(OUT, TRUE, showWarnings=FALSE)
## inflated sample names (unique) and the copy index maps
aqIdx  <- rep(seq_len(nAq),  times=C)                      # length 30*C, values in 1..30
polIdx <- rep(seq_len(nPol), times=C)                      # length 13*C, values in 1..13
aq_s  <- sprintf("aq_hap%04d",  seq_along(aqIdx))
pol_s <- sprintf("pol_fem%04d", seq_along(polIdx))

for (ch in unique(map$Chr[sel])) {
  ci <- sel[map$Chr[sel]==ch]; ci <- ci[order(map$Pos[ci])]
  pos <- map$Pos[ci]; ids <- map$marker[ci]; chs <- sub("Chr","ch",ch)
  gaq0  <- matrix(as.character(t(Haq[,ci,drop=FALSE])), nrow=length(ci))          # markers x 30
  H1 <- Hpol[seq(1,nrow(Hpol),2), ci, drop=FALSE]; H2 <- Hpol[seq(2,nrow(Hpol),2), ci, drop=FALSE]
  gpol0 <- matrix(paste(t(H1), t(H2), sep="|"), nrow=length(ci))                  # markers x 13
  gaq  <- gaq0[,  aqIdx,  drop=FALSE]                       # markers x (30*C) copies
  gpol <- gpol0[, polIdx, drop=FALSE]                       # markers x (13*C) copies
  hdr <- c("##fileformat=VCFv4.2", sprintf("##contig=<ID=%s>", chs),
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
           paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",aq_s,pol_s), collapse="\t"))
  body <- paste(sprintf("%s\t%d\t%s\tA\tC\t.\tPASS\t.\tGT", chs, pos, ids),
                apply(cbind(gaq,gpol), 1, paste, collapse="\t"), sep="\t")
  writeLines(c(hdr, body), file.path(OUT, sprintf("founders_%s.vcf", chs)))
}
cat(sprintf("wrote %d markers x (%d aq / %d pol) inflated founders to %s (%d chr)\n",
            length(sel), length(aqIdx), length(polIdx), OUT, length(unique(map$Chr[sel]))))
cat(sprintf("  -> run SLiM with  N_AQ_POOL=%d N_POL_POOL=%d  and sweep N_AQ/N_POL up to those.\n",
            length(aqIdx), length(polIdx)))
