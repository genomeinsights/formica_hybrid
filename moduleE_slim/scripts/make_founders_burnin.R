## =============================================================================
## Module E -- DIVERSIFIED founders via a per-species RECOMBINATION burn-in
## =============================================================================
## Breaks the finite-founder short-range LD before the hybrid phase. Per chromosome,
## per species (aquilonia / polyctena SEPARATELY, so the diagnostic between-species
## differences are untouched), seed a large panmictic pool from the real haplotypes
## and apply recombination on the empirical cM map SCALED UP by RECOMB_SCALE with NO
## mutation. Recombination is frequency-neutral, so allele frequencies (hence DI) are
## preserved up to a little drift; scaling by S reaches short-range linkage
## equilibrium in 1/S the generations (validated in dev/R/moduleE_burnin_prototype.R).
## High-recombination regions relax to LE (killing the finite-panel excess); low-
## recombination regions keep their recombination-limited LD (the real landscape).
## Output founders are read by the hybrid model exactly like the inflated pool.
##
## Usage: Rscript make_founders_burnin.R [S] [T] [N]   (defaults 100 8 1000)
##   -> ../founders/maf015_DIstrat1500_burnin/   (1200 aq haploid / 520 pol diploid)
## =============================================================================
suppressMessages(library(data.table))
.args <- commandArgs(trailingOnly = FALSE); .self <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT  <- normalizePath(file.path(dirname(.self), ".."))
ca <- commandArgs(trailingOnly = TRUE)
S <- as.numeric(ca[1]); if (is.na(S)) S <- 100          # recombination scale
Tg<- as.integer(ca[2]); if (is.na(Tg)) Tg <- 8L         # burn-in generations (at scale S)
N <- as.integer(ca[3]); if (is.na(N)) N <- 1000L        # burn-in pop size (diploids)
N_AQ_OUT <- 1200L; N_POL_OUT <- 520L                    # output pool (matches x40 inflated sizes)
SRC <- file.path(ROOT, "founders/maf015_DIstrat1500_x40")  # marker set to match
OUT <- file.path(ROOT, "founders/maf015_DIstrat1500_burnin"); dir.create(OUT, TRUE, showWarnings=FALSE)

ph   <- readRDS(file.path(ROOT, "inputs/moduleE_founder_haplotypes.rds"))
pmap <- as.data.table(ph$map)[, idx := .I]; setkey(pmap, marker)
Haq  <- ph$aqu; Hpol <- ph$pol                          # 30 x M_all ; 26 x M_all
cat(sprintf("burn-in: S=%g  T=%d gens (U=%g)  N=%d ; out pool %d aq / %d pol\n", S, Tg, S*Tg, N, N_AQ_OUT, N_POL_OUT))

## one WF recombination step (panmictic, no mutation): 2N gametes, each a mosaic of
## two random current haplotypes with crossovers per the S-scaled Haldane distances
recomb_step <- function(pool, r_switch) {
  n <- nrow(pool); M <- ncol(pool)
  Ha <- pool[sample.int(n, n, replace=TRUE),,drop=FALSE]; Hb <- pool[sample.int(n, n, replace=TRUE),,drop=FALSE]
  sw <- matrix(as.integer(runif(n*(M-1)) < r_switch), nrow=n)
  ch <- matrix(0L, n, M); ch[,1] <- sample.int(2, n, replace=TRUE)-1L
  for (j in 2:M) ch[,j] <- bitwXor(ch[,j-1], sw[,j-1])
  Ha*(1L-ch) + Hb*ch
}
burn <- function(H0, r_switch) {                        # H0: h x M real haplotypes -> 2N burned
  pool <- H0[sample.int(nrow(H0), 2*N, replace=TRUE),,drop=FALSE]
  for (t in seq_len(Tg)) pool <- recomb_step(pool, r_switch)
  pool
}

vcf_files <- list.files(SRC, "founders_ch.*\\.vcf$", full.names=TRUE)
dfreq_all <- c()
for (f in vcf_files) {
  hdr <- fread(f, skip="#CHROM", nrows=0)                # marker table for this chromosome
  mk  <- fread(f, skip="#CHROM", select=c(1,2,3), header=TRUE, col.names=c("chs","Pos","marker"))
  chs <- mk$chs[1]; id <- as.integer(sub("ch","",chs))
  sel <- pmap[mk$marker]; ord <- order(mk$Pos); mk <- mk[ord]; sel <- sel[ord]
  Mc  <- nrow(mk); Pos <- mk$Pos
  Haq0 <- Haq[, sel$idx, drop=FALSE]; Hpol0 <- Hpol[, sel$idx, drop=FALSE]
  ## cM map -> inter-marker Haldane switch prob, scaled by S
  rc <- fread(file.path(ROOT,"rec_maps",sprintf("ch_%d.recmap",id)), header=FALSE, col.names=c("pos","rate"))
  rc <- rc[is.finite(pos)&is.finite(rate)][order(pos)]
  cM <- approxfun(rc$pos, c(0, cumsum(rc$rate[-nrow(rc)]*diff(rc$pos)/1e6)), rule=2)(Pos)
  dM <- pmax(diff(cM)/100, 0); r_switch <- 0.5*(1 - exp(-2*dM*S))
  ## burn each species, take the output pools
  set.seed(1000+id)
  aqp  <- burn(Haq0,  r_switch)[seq_len(N_AQ_OUT),,drop=FALSE]        # 1200 x Mc haploid
  polp <- burn(Hpol0, r_switch)[seq_len(2*N_POL_OUT),,drop=FALSE]     # 1040 x Mc -> 520 diploids
  dfreq_all <- c(dfreq_all,
    mean(abs(colMeans(aqp)-colMeans(Haq0))), mean(abs(colMeans(polp)-colMeans(Hpol0))))
  ## assemble VCF (same format as make_founders_inflated.R)
  aq_s  <- sprintf("aq_hap%04d",  seq_len(N_AQ_OUT)); pol_s <- sprintf("pol_fem%04d", seq_len(N_POL_OUT))
  gaq  <- t(aqp)                                                       # Mc x 1200 alleles
  H1 <- polp[seq(1,nrow(polp),2),,drop=FALSE]; H2 <- polp[seq(2,nrow(polp),2),,drop=FALSE]
  gpol <- matrix(paste(t(H1), t(H2), sep="|"), nrow=Mc)                # Mc x 520 "a|b"
  head <- c("##fileformat=VCFv4.2", sprintf("##contig=<ID=%s>", chs),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",aq_s,pol_s), collapse="\t"))
  body <- paste(sprintf("%s\t%d\t%s\tA\tC\t.\tPASS\t.\tGT", chs, Pos, mk$marker),
                apply(cbind(gaq,gpol), 1, paste, collapse="\t"), sep="\t")
  writeLines(c(head, body), file.path(OUT, sprintf("founders_%s.vcf", chs)))
  cat(sprintf("  %s: %d markers burned\n", chs, Mc))
}
cat(sprintf("\nwrote diversified founders to %s\n", OUT))
cat(sprintf("mean |delta allele-freq| from burn-in drift: %.4f  (frequencies/DI ~preserved)\n", mean(dfreq_all)))
cat(sprintf("-> run hybrid model with FDIR=%s/  N_AQ_POOL=%d N_POL_POOL=%d\n", OUT, N_AQ_OUT, N_POL_OUT))
