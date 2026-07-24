## =============================================================================
## Module E -- regenerate DI-stratified founder VCFs (PORTABLE)
## =============================================================================
## Draws `n_per_bin` markers from each of 10 fixed DI bins over the parent_maf>=0.15
## universe, and writes one combined per-chromosome founder VCF (30 haploid aq +
## 13 diploid pol samples) that the SLiM model reads. Use a smaller n_per_bin to
## cut per-deme memory when pushing K higher (memory scales ~linearly with markers).
##
## Usage:  Rscript make_founders.R <n_per_bin>   e.g. 2500  -> ../founders/maf015_DIstrat2500/
## =============================================================================
suppressMessages(library(data.table))
.args <- commandArgs(trailingOnly = FALSE)
.self <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT  <- normalizePath(file.path(dirname(.self), ".."))
N <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(N)) stop("give n_per_bin, e.g. 2500")

DI_BREAKS <- c(-Inf,-90,-75,-60,-50,-40,-30,-25,-20,-15,Inf); MAF_GATE <- 0.15; SEED <- 1L
ph  <- readRDS(file.path(ROOT, "inputs/moduleE_founder_haplotypes.rds"))
map <- as.data.table(ph$map)
Haq <- ph$aqu; Hpol <- ph$pol
p <- (colSums(Haq)+colSums(Hpol))/(nrow(Haq)+nrow(Hpol)); map[, parent_maf := pmin(p,1-p)]
uni <- which(map$parent_maf >= MAF_GATE & is.finite(map$DI))
set.seed(SEED)
sel <- sort(unlist(lapply(split(uni, cut(map$DI[uni], DI_BREAKS)),
                          function(ix) if (length(ix)<=N) ix else sample(ix, N))))
OUT <- file.path(ROOT, sprintf("founders/maf015_DIstrat%d", N)); dir.create(OUT, TRUE, showWarnings=FALSE)
aq_s  <- sprintf("aq_hap%02d", seq_len(nrow(Haq)))
pol_s <- sprintf("pol_fem%02d", seq_len(nrow(Hpol)%/%2L))
for (ch in unique(map$Chr[sel])) {
  ci <- sel[map$Chr[sel]==ch]; ci <- ci[order(map$Pos[ci])]
  pos <- map$Pos[ci]; ids <- map$marker[ci]; chs <- sub("Chr","ch",ch)
  gaq  <- matrix(as.character(t(Haq[,ci,drop=FALSE])), nrow=length(ci))
  H1 <- Hpol[seq(1,nrow(Hpol),2), ci, drop=FALSE]; H2 <- Hpol[seq(2,nrow(Hpol),2), ci, drop=FALSE]
  gpol <- matrix(paste(t(H1), t(H2), sep="|"), nrow=length(ci))
  hdr <- c("##fileformat=VCFv4.2", sprintf("##contig=<ID=%s>", chs),
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
           paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",aq_s,pol_s), collapse="\t"))
  body <- paste(sprintf("%s\t%d\t%s\tA\tC\t.\tPASS\t.\tGT", chs, pos, ids),
                apply(cbind(gaq,gpol), 1, paste, collapse="\t"), sep="\t")
  writeLines(c(hdr, body), file.path(OUT, sprintf("founders_%s.vcf", chs)))
}
cat(sprintf("wrote %d markers to %s (%d chromosomes)\n", length(sel), OUT, length(unique(map$Chr[sel]))))
