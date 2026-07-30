## =============================================================================
## Module E -- empirical vs neutral-null summary panels: the F_ST / LD dissociation
## =============================================================================
## One self-contained figure that shows the empirical-vs-null discrepancy at two
## scales, using ONE LD metric throughout (LD between a SNP and its directly
## adjacent SNP(s) = ld_w at a one-SNP window):
##   TOP  aggregate dose-response vs DiagnosticIndex:
##        (A) among-population F_ST  -- empirical RISES, every neutral demography flat
##        (B) adjacent-SNP LD        -- the neutral null RISES, empirical stays flat
##   BOTTOM per-SNP genomic tracks in one high-recombination window (Chr14 2-4 Mb):
##        (C) per-SNP F_ST           (D) adjacent-SNP LD      (E) sorting direction (uni_score)
## The two nulls (single panmictic pool, 20 independent demes) are both DI-blind for
## F_ST; for pooled LD they coincide (same individuals), so LD uses one null line.
##
## Inputs : moduleE_slim/inputs/empirical_bundle.rds, data/hybrids_only_maf005.Rdata,
##          moduleE_slim/output_burnin/ (450/195, K6250, gen60), moduleE_slim/rec_maps/
## Outputs: Figures/moduleE_emp_vs_sim_summary.{pdf,png} + data/moduleE_sim/moduleE_emp_vs_sim_panels.rds
## =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

F        <- "/Users/petrikem/gitlab/formica_hybrid"
SIMDIR   <- file.path(F, "moduleE_slim/output_burnin"); TAG <- "Naq450_Npol195"; K <- 6250; GEN <- 60L
FVCF     <- file.path(F, "moduleE_slim/founders/maf015_DIstrat1500_burnin")
BUNDLE   <- file.path(F, "moduleE_slim/inputs/empirical_bundle.rds")
EMP_RD   <- file.path(F, "data/hybrids_only_maf005.Rdata")
FIG      <- file.path(F, "Figures"); OUT <- file.path(F, "data/moduleE_sim")
WIN_CHR  <- "Chr14"; WIN_MB <- c(2, 4)                       # the per-SNP zoom window
COL      <- c(empirical = "#d95f02", "neutral null" = "#1b9e77", "null: panmictic pool" = "#7570b3")

b   <- readRDS(BUNDLE); uni <- as.data.table(b$universe); setkey(uni, marker); emp_ns <- b$emp_ns
fm  <- rbindlist(lapply(list.files(FVCF, "founders_ch.*[.]vcf$", full.names = TRUE),
        function(f) fread(f, skip = "#CHROM", select = c(1,2,3), header = TRUE, col.names = c("chs","Pos","marker"))))
fm[, `:=`(Chr = sub("ch","Chr",chs), chr_id = as.integer(sub("ch","",chs)))]; setorder(fm, chr_id, Pos)
markers <- fm$marker; DI <- uni[markers, DI]
faq <- b$f_aq_par[markers]; fpol <- b$f_pol_par[markers]; flip <- faq < fpol

read_dose <- function(f) { dt <- fread(f, skip="#CHROM", header=TRUE, sep="\t", showProgress=FALSE)
  mk <- paste(sub("ch","Chr",dt[[1]]), dt[[2]], sep=":"); G <- as.matrix(dt[,10:ncol(dt)]); gt <- sub(":.*$","",G)
  d <- matrix(as.integer(sub("[|/].*$","",gt)) + as.integer(sub("^.*[|/]","",gt)), nrow=nrow(G)); rownames(d) <- mk
  t(d)[, markers, drop=FALSE] }

## ---- per-population / per-deme allele frequencies ---------------------------
Fe <- b$emp_mean[, markers, drop=FALSE] / 2000                                   # empirical, 20 pops (n = 3-20)
demes <- lapply(1:20, function(i) read_dose(file.path(SIMDIR, sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf", TAG, K, i, GEN))))
## null = 20 INDEPENDENT demes, each downsampled to the matching empirical population size
## (one sample per deme; NOT pooled into a single population, NOT full n) -- bias-matched to Fe
subs <- lapply(1:20, function(i){ set.seed(400+i); D <- demes[[i]]; D[sample(nrow(D), min(emp_ns[i], nrow(D))),,drop=FALSE] })
Fn <- t(vapply(subs, function(D) colMeans(D)/2, numeric(length(markers))))         # null: 20 independent demes, sample-matched
Gn <- do.call(rbind, subs)                                                          # null pooled genotypes

## ---- pooled empirical genotypes --------------------------------------------
e <- new.env(); load(EMP_RD, envir=e); sdh <- as.data.table(e$sample_data); maph <- as.data.table(e$map_hyb_005)
hyb <- sdh[!Population %in% c("aquilonia_parent","polyctena_parent"), Sample_ID]
Ge <- e$GTs_hybrids_005[match(hyb, rownames(e$GTs_hybrids_005)), match(markers, maph$marker), drop=FALSE]; Ge[is.na(Ge)] <- 0L; rm(e); gc()

## ---- per-SNP statistics -----------------------------------------------------
fst_snp  <- function(Fr){ pb<-colMeans(Fr,na.rm=TRUE); Ht<-2*pb*(1-pb); Hs<-colMeans(2*Fr*(1-Fr),na.rm=TRUE); ifelse(Ht>0,(Ht-Hs)/Ht,NA_real_) }
uni_snp  <- function(Fr){ Fo<-Fr; Fo[,flip]<-1-Fo[,flip]; (colSums(Fo>=0.85,na.rm=TRUE)-colSums(Fo<=0.15,na.rm=TRUE))/nrow(Fr) }  # +=aq
adjld    <- function(G){ res<-rep(NA_real_,length(markers))
  for (cid in unique(fm$chr_id)) { j<-which(fm$chr_id==cid); if(length(j)<2) next
    for (a in seq_along(j)) { nb<-j[c(a-1,a+1)]; nb<-nb[nb>=min(j)&nb<=max(j)&nb!=j[a]]; nb<-nb[!is.na(nb)]; if(!length(nb)) next
      res[j[a]] <- mean(suppressWarnings(cor(G[,j[a]], G[,nb,drop=FALSE]))^2, na.rm=TRUE) } }
  res }
d <- data.table(marker=markers, Chr=fm$Chr, chr_id=fm$chr_id, Pos=fm$Pos, DI=DI,
                fst_emp=fst_snp(Fe), fst_sim=fst_snp(Fn), uni_emp=uni_snp(Fe), uni_sim=uni_snp(Fn),
                ld_emp=adjld(Ge), ld_sim=adjld(Gn))
## per-marker recombination (cM/Mb), map-interpolated
d[, recomb := NA_real_]
for (cid in unique(d$chr_id)) { rc <- fread(file.path(F,"moduleE_slim/rec_maps",sprintf("ch_%d.recmap",cid)), header=FALSE, col.names=c("pos","r"))
  rc <- rc[is.finite(pos)&is.finite(r)][order(pos)]; d[chr_id==cid, recomb := approx(rc$pos, rc$r, xout=Pos, rule=2)$y] }

## ---- (A,B) dose-response vs DI ---------------------------------------------
BR <- c(-Inf,-90,-75,-60,-50,-40,-30,-25,-20,-15,Inf); d[, dib := cut(DI, BR)]
## bin midpoint and F_ST computed TOGETHER per bin (keyby=dib) so they can't misalign
fst_bin <- function(Fr){ pb<-colMeans(Fr,na.rm=TRUE); Ht<-2*pb*(1-pb); Hs<-colMeans(2*Fr*(1-Fr),na.rm=TRUE)
  data.table(dib=d$dib, di=DI, Ht=Ht, Hs=Hs)[, .(mid=mean(di,na.rm=TRUE), v=sum(Ht-Hs,na.rm=TRUE)/sum(Ht,na.rm=TRUE)), keyby=dib] }
doseA <- rbindlist(list(cbind(fst_bin(Fe), src="empirical"), cbind(fst_bin(Fn), src="neutral null")))
ldbin <- d[, .(emp=median(ld_emp,na.rm=TRUE), sim=median(ld_sim,na.rm=TRUE), mid=mean(DI,na.rm=TRUE)), keyby=dib]
doseB <- rbindlist(list(data.table(mid=ldbin$mid, v=ldbin$emp, src="empirical"),
                        data.table(mid=ldbin$mid, v=ldbin$sim, src="neutral null")))

pA <- ggplot(doseA, aes(mid, v, colour=src)) + geom_line(linewidth=.6) + geom_point(size=1.6) +
  scale_colour_manual(values=COL, name=NULL) +
  labs(x="DiagnosticIndex", y=expression(F[ST]), title="A  Differentiation vs DI") +
  theme_bw(base_size=9) + theme(legend.position="top", legend.text=element_text(size=7))
pB <- ggplot(doseB, aes(mid, v, colour=src)) + geom_line(linewidth=.6) + geom_point(size=1.6) +
  scale_colour_manual(values=COL, name=NULL) +
  labs(x="DiagnosticIndex", y="adjacent-SNP LD", title="B  Adjacent-SNP LD vs DI") +
  theme_bw(base_size=9) + theme(legend.position="top", legend.text=element_text(size=7))

## ---- (C-E rows x low/mid/high DI columns) per-SNP DISTRIBUTIONS, all data ---
## restricted to high-recombination SNPs so adjacent-SNP LD reflects ancestry
## (admixture) rather than physical linkage
hi <- d[is.finite(recomb) & recomb >= median(d$recomb, na.rm=TRUE)]
DILAB <- c("low DI: neutral (<= -90)","mid DI (-90..-25)","high DI: diagnostic (> -25)")
hi[, DIclass := factor(cut(DI, c(-Inf,-90,-25,Inf), labels=DILAB), levels=DILAB)]
cat("high-recomb SNPs per DI class: "); print(table(hi$DIclass))
mkv <- function(vE, vS, lab) rbind(hi[, .(DIclass, val=get(vE), src="empirical",   metric=lab)],
                                   hi[, .(DIclass, val=get(vS), src="neutral null", metric=lab)])
vl <- rbind(mkv("fst_emp","fst_sim","C  per-SNP F_ST"),
            mkv("ld_emp","ld_sim","D  adjacent-SNP LD"),
            mkv("uni_emp","uni_sim","E  sorting direction"))
vl[, metric := factor(metric, levels=c("C  per-SNP F_ST","D  adjacent-SNP LD","E  sorting direction"))]
pGrid <- ggplot(vl, aes(src, val, fill=src, colour=src)) +
  geom_hline(yintercept=0, colour="grey85", linewidth=.2) +
  geom_violin(alpha=.30, linewidth=.3, scale="width") +
  geom_boxplot(width=.16, alpha=.65, outlier.size=.1, linewidth=.3) +
  facet_grid(metric ~ DIclass, scales="free_y", switch="y") +
  scale_fill_manual(values=COL, name=NULL) + scale_colour_manual(values=COL, name=NULL) +
  labs(x=NULL, y=NULL, title="C-E  Per-SNP distributions by DI class (all high-recombination SNPs): the dissociation appears only at high DI") +
  theme_bw(base_size=9) + theme(legend.position="top", panel.grid.minor=element_blank(),
        strip.placement="outside", strip.background=element_rect(fill="grey92"),
        legend.text=element_text(size=7), axis.text.x=element_text(size=8))

comp <- (pA | pB) / pGrid + plot_layout(heights=c(1, 2.2)) +
  plot_annotation(title="Empirical vs neutral null: high F_ST but low LD at ancestry-informative loci",
                  subtitle="Aggregate over all data (A,B) and per-SNP distributions by DI class (C-E). LD throughout = mean r^2 between a SNP and its directly adjacent SNP(s).",
                  theme=theme(plot.title=element_text(face="bold", size=12)))
ggsave(file.path(FIG,"moduleE_emp_vs_sim_summary.pdf"), comp, width=11, height=9.5)
ggsave(file.path(FIG,"moduleE_emp_vs_sim_summary.png"), comp, width=11, height=9.5, dpi=140)
byDI <- hi[, .(n=.N, fst_emp=round(median(fst_emp,na.rm=T),3), fst_null=round(median(fst_sim,na.rm=T),3),
               ld_emp=round(median(ld_emp,na.rm=T),3), ld_null=round(median(ld_sim,na.rm=T),3),
               absuni_emp=round(median(abs(uni_emp),na.rm=T),3), absuni_null=round(median(abs(uni_sim),na.rm=T),3)), keyby=DIclass]
saveRDS(list(per_snp=d, doseA=doseA, doseB=doseB, byDI=byDI), file.path(OUT,"moduleE_emp_vs_sim_panels.rds"))
cat("\nper-SNP medians by DI class (high-recombination SNPs):\n"); print(byDI)
cat("saved: Figures/moduleE_emp_vs_sim_summary.{pdf,png}\n")
