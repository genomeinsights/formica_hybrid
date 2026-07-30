## =============================================================================
## Module E -- PROTOTYPE: within-species recombination burn-in to break the
## founder-haplotype short-range LD, with recombination-rate rescaling
## =============================================================================
## Tests option-2 (recombination burn-in) + the rescaling idea: seed a large
## panmictic population from the real within-species founder haplotypes, apply
## recombination (empirical cM map, scaled up by S), NO mutation, and let short-range
## LD relax toward the large-population equilibrium implied by the sampled allele
## frequencies. Recombination is frequency-neutral, so at large N allele frequencies
## (hence DI) are preserved; scaling recombination by S reaches a given amount of LD
## breakdown in 1/S the generations AND with 1/S the drift.
##
## Claims validated here (aquilonia founders, chromosome 1, real cM map):
##   (1) MECHANISM: burn-in lowers short-range ld_w from the 30-haplotype value.
##   (2) RESCALING INVARIANCE: ld_w as a function of cumulative recombination U = S*t
##       coincides for S=20 and S=100 -> the speedup is free.
##   (3) FREQUENCY PRESERVATION: mean |delta allele-freq| is tiny and smaller at higher S.
## Genome-wide production would run this in SLiM per species; this R version validates
## the physics and calibration cheaply.
## Output: Figures/moduleE_burnin_prototype.{pdf,png} + printed summary.
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2) })
FORMICA <- "/Users/petrikem/gitlab/formica_hybrid"
FIGDIR  <- file.path(FORMICA, "Figures")

## ---- data: aquilonia founder haplotypes + cM map on chr1 -------------------
ph  <- readRDS(file.path(FORMICA, "data/moduleE_founders/moduleE_founder_haplotypes.rds"))
pm  <- as.data.table(ph$map)[, .(Chr, Pos, marker, idx=.I)]
mk1 <- fread(file.path(FORMICA,"moduleE_slim/founders/maf015_DIstrat1500_x40/founders_ch1.vcf"),
             skip="#CHROM", select=3, header=TRUE)[[1]]                 # the 15k-set chr1 markers
sel <- pm[marker %in% mk1][order(Pos)]
H0  <- ph$aqu[, sel$idx, drop=FALSE]                                    # 30 haplotypes x M
M   <- ncol(H0); Pos <- sel$Pos
cat(sprintf("aquilonia founders: %d haplotypes x %d chr1 markers\n", nrow(H0), M))

## cumulative cM at each marker from the empirical recombination map (cM/Mb) --
rc <- fread(file.path(FORMICA,"moduleE_slim/rec_maps/ch_1.recmap"), header=FALSE, col.names=c("pos","rate"))
rc <- rc[is.finite(pos) & is.finite(rate)][order(pos)]
cumcM_at <- approxfun(rc$pos, c(0, cumsum(rc$rate[-nrow(rc)] * diff(rc$pos)/1e6)), rule=2)
cM  <- cumcM_at(Pos)                                                    # cumulative cM per marker
dM  <- pmax(diff(cM)/100, 0)                                            # inter-marker Morgans (M-1)
cat(sprintf("chr1 spans %.0f cM; median inter-marker gap %.1f kb / %.4f cM\n",
            max(cM), median(diff(Pos))/1e3, median(diff(cM))))

## ---- WF recombination step (panmictic, freq-neutral, no mutation) ----------
## pool: (2N) x M haplotype matrix. Each new gamete recombines two random current
## haplotypes; crossovers via Haldane map function on the S-scaled distances.
recomb_step <- function(pool, r_switch) {
  n <- nrow(pool)
  pa <- sample.int(n, n, replace=TRUE); pb <- sample.int(n, n, replace=TRUE)
  Ha <- pool[pa,,drop=FALSE]; Hb <- pool[pb,,drop=FALSE]
  sw <- matrix(as.integer(runif(n*(M-1)) < r_switch), nrow=n)          # switch between marker j,j+1
  ch <- matrix(0L, n, M); ch[,1] <- sample.int(2, n, replace=TRUE)-1L   # random starting parent
  for (j in 2:M) ch[,j] <- bitwXor(ch[,j-1], sw[,j-1])                  # running parity
  Ha*(1L-ch) + Hb*ch
}
## short-range ld_w on DIPLOID dosages (median r^2 among neighbours <=100kb, MAF>=0.1)
WIN <- 1e5L; MAXP <- 60L
ldw <- function(pool) {
  n <- nrow(pool); D <- pool[seq(1,n,2),,drop=FALSE] + pool[seq(2,n,2),,drop=FALSE]  # dosage
  p <- colMeans(D)/2; poly <- which(pmin(p,1-p) >= 0.1)
  res <- rep(NA_real_, M)
  for (ai in seq_along(poly)) { a <- poly[ai]
    b <- poly[poly>a & Pos[poly]<=Pos[a]+WIN]; if(length(b)>MAXP) b<-b[seq_len(MAXP)]; if(!length(b)) next
    res[a] <- median(suppressWarnings(cor(D[,a], D[,b,drop=FALSE]))^2, na.rm=TRUE) }
  c(mean=mean(res,na.rm=TRUE), median=median(res,na.rm=TRUE))
}

run_burnin <- function(S, snaps, N=800, seed=1) {
  set.seed(seed)
  r_switch <- 0.5*(1 - exp(-2 * dM * S))                               # Haldane, S-scaled
  pool <- H0[sample.int(nrow(H0), 2*N, replace=TRUE),,drop=FALSE]      # seed from the real haplotypes
  f0 <- colMeans(pool)/1                                               # per-haplotype allele freq
  out <- list(); tmax <- max(snaps)
  for (t in seq_len(tmax)) { pool <- recomb_step(pool, r_switch)
    if (t %in% snaps) { l <- ldw(pool)
      out[[length(out)+1]] <- data.table(S=S, t=t, U=S*t, ldw_mean=l["mean"], ldw_med=l["median"],
                                          dfreq=mean(abs(colMeans(pool)-f0))) } }
  rbindlist(out)
}

## ---- baseline (U=0) + two S regimes at matched U = S*t ----------------------
base <- { l<-ldw(H0[sample.int(nrow(H0),1600,replace=TRUE),,drop=FALSE])
          data.table(S=NA,t=0,U=0,ldw_mean=l["mean"],ldw_med=l["median"],dfreq=0) }
cat(sprintf("\nbaseline (30-haplotype founders): mean ld_w=%.3f  median=%.3f\n", base$ldw_mean, base$ldw_med))

U <- c(100,200,400,800,1600)
r20  <- run_burnin(S=20,  snaps=U/20)     # t = 5,10,20,40,80
r100 <- run_burnin(S=100, snaps=U/100)    # t = 1,2,4,8,16   (5x fewer generations)
res <- rbind(base, r20, r100, fill=TRUE)

cat("\n=== rescaling-invariance check (ld_w at matched cumulative recombination U=S*t) ===\n")
cmp <- merge(r20[,.(U,S20_med=round(ldw_med,3),S20_gens=t,S20_dfreq=round(dfreq,4))],
             r100[,.(U,S100_med=round(ldw_med,3),S100_gens=t,S100_dfreq=round(dfreq,4))], by="U")
print(cmp)
cat(sprintf("\nmedian |ld_w(S=20) - ld_w(S=100)| across matched U: %.4f  (should be ~0)\n",
            median(abs(cmp$S20_med-cmp$S100_med))))
cat(sprintf("frequency drift at U=800: S=20 -> %.4f (t=40 gens) vs S=100 -> %.4f (t=8 gens)\n",
            cmp[U==800,S20_dfreq], cmp[U==800,S100_dfreq]))

## empirical hybrid chr1 ld_w as an illustrative target line -------------------
emp_line <- tryCatch({
  e<-new.env(); load(file.path(FORMICA,"data/hybrids_only_maf005.Rdata"),envir=e)
  sdh<-as.data.table(e$sample_data); maph<-as.data.table(e$map_hyb_005)
  hyb<-which(!sdh$Population%in%c("aquilonia_parent","polyctena_parent"))
  g<-e$GTs_hybrids_005[hyb, match(sel$marker, maph$marker),drop=FALSE]; g[is.na(g)]<-0L
  p<-colMeans(g)/2; poly<-which(pmin(p,1-p)>=0.1); res<-rep(NA_real_,M)
  for(ai in seq_along(poly)){a<-poly[ai]; b<-poly[poly>a & Pos[poly]<=Pos[a]+WIN]; if(length(b)>MAXP)b<-b[seq_len(MAXP)]; if(!length(b))next
    res[a]<-median(suppressWarnings(cor(g[,a],g[,b,drop=FALSE]))^2,na.rm=TRUE)}
  median(res,na.rm=TRUE) }, error=function(e) NA)
cat(sprintf("\nempirical hybrid chr1 median ld_w (illustrative target): %.3f\n", emp_line))

## ---- figure ----------------------------------------------------------------
pl <- rbind(r20[,.(U,ldw_med,S=factor("S=20 (t=U/20)"))],
            r100[,.(U,ldw_med,S=factor("S=100 (t=U/100)"))])
pl <- rbind(data.table(U=0,ldw_med=base$ldw_med,S="S=20 (t=U/20)"),
            data.table(U=0,ldw_med=base$ldw_med,S="S=100 (t=U/100)"), pl)
p <- ggplot(pl, aes(U, ldw_med, colour=S)) +
  geom_line(linewidth=.7) + geom_point(size=2) +
  { if(is.finite(emp_line)) geom_hline(yintercept=emp_line, linetype=2, colour="#d95f02") } +
  { if(is.finite(emp_line)) annotate("text",x=0,y=emp_line,label="empirical hybrid chr1",
       hjust=-0.02,vjust=-0.6,size=3,colour="#d95f02") } +
  scale_colour_manual(values=c("S=20 (t=U/20)"="#1b9e77","S=100 (t=U/100)"="#7570b3"),name=NULL)+
  labs(x="cumulative recombination  U = S × t  (scaled generations)",
       y="short-range median ld_w",
       title="Recombination burn-in: short-range LD vs cumulative recombination",
       subtitle="two recombination scalings collapse onto one curve → higher S reaches the same LD in 1/S the generations") +
  theme_bw(base_size=10)
ggsave(file.path(FIGDIR,"moduleE_burnin_prototype.pdf"), p, width=7.5, height=4.5)
ggsave(file.path(FIGDIR,"moduleE_burnin_prototype.png"), p, width=7.5, height=4.5, dpi=130)
cat("\nsaved: Figures/moduleE_burnin_prototype.{pdf,png}\n")
