## =========================================================
## DI-threshold sweep: how the within-species LD "problem" and marker
## density trade off as the diagnostic-index cutoff is tightened.
## =========================================================
##
## Extends parent_LD_diagnostic.R across DI in {-15,-20,-25,-30}. LD pairs
## are computed ONCE on the loosest (DI > -30) superset with each pair
## tagged by the minimum DI of its two markers; a pair is "present" at a
## stricter threshold iff both markers pass it. Recomb stratum split and
## panel membership (n = 15) are held fixed across thresholds so the only
## thing changing is the marker set.

suppressMessages({library(data.table); library(ggplot2)})
set.seed(1)

e <- new.env(); load("./data/hybrids_and_parents_maf005.Rdata", envir = e)
GTs <- e$GTs_with_parents
sd  <- e$sample_data_with_parents
map <- copy(e$map_hyb_005)

THRESHOLDS <- c(-15, -20, -25, -30)
SUPER_TH <- min(THRESHOLDS)           # -30 superset
MAXDIST  <- 200000
MAXPART  <- 50
N_IND    <- 15
MAF_MIN  <- 0.05

di <- map[DiagnosticIndex > SUPER_TH, which = TRUE]
cat(sprintf("superset DI > %d: %d markers\n", SUPER_TH, length(di)))

## per-marker recombination (cM/Mb)
rec <- fread("./data/Frufa_DTOL_PR.ref_genome.recmap")
setnames(rec, c("chr","pos","cM","cMMb")); rec[, Chr := sub("chromosome_","Chr",chr)]
map[, recomb := NA_real_]
for (ch in unique(map$Chr)) {
  r <- rec[Chr == ch]; if (nrow(r) < 2) next
  idx <- map[, which(Chr == ch)]
  map[idx, recomb := approx(r$pos, r$cMMb, xout = map$Pos[idx], rule = 2)$y]
}
rec_split <- median(map$recomb[di], na.rm = TRUE)      # fixed split for all thresholds
cat(sprintf("fixed low/high recomb split: %.2f cM/Mb\n", rec_split))

## fixed panels (n = 15)
rows_of <- function(pops, n) { idx <- which(sd$Population %in% pops); if (length(idx) > n) idx <- sample(idx, n); idx }
panels <- list(
  aqu    = rows_of("aquilonia_parent", N_IND),
  pol    = rows_of("polyctena_parent", N_IND),
  pooled = c(sample(which(sd$Population=="aquilonia_parent"),8),
             sample(which(sd$Population=="polyctena_parent"),7))
)

## ---- within-species polymorphism vs threshold ----
p_aqu <- colMeans(GTs[panels$aqu, di, drop=FALSE], na.rm=TRUE)/2; maf_aqu <- pmin(p_aqu,1-p_aqu)
p_pol <- colMeans(GTs[panels$pol, di, drop=FALSE], na.rm=TRUE)/2; maf_pol <- pmin(p_pol,1-p_pol)
di_val <- map$DiagnosticIndex[di]

poly_tab <- rbindlist(lapply(THRESHOLDS, function(t) {
  k <- di_val > t
  data.table(threshold=t, n_markers=sum(k),
             poly_aqu=round(100*mean(maf_aqu[k] > MAF_MIN, na.rm=TRUE),1),
             poly_pol=round(100*mean(maf_pol[k] > MAF_MIN, na.rm=TRUE),1))
}))
cat("\n=== within-species polymorphism (MAF>0.05) & marker count by threshold ===\n")
print(poly_tab)

## ---- LD pairs on superset, tagged with pair-min-DI ----
ld_pairs_panel <- function(rows) {
  G <- GTs[rows, di, drop=FALSE]
  info <- map[di, .(Chr, Pos, recomb)]; dv <- di_val
  poly <- apply(G, 2, function(x){v<-x[!is.na(x)]; length(v)>=6 && length(unique(v))>1})
  out <- list()
  for (ch in unique(info$Chr)) {
    j <- which(info$Chr==ch & poly); if (length(j)<2) next
    o <- order(info$Pos[j]); j <- j[o]
    pos <- info$Pos[j]; rr <- info$recomb[j]; dd <- dv[j]; Gc <- G[, j, drop=FALSE]
    for (a in seq_len(length(j)-1)) {
      b <- (a+1):min(a+MAXPART, length(j)); dist <- pos[b]-pos[a]
      keep <- dist <= MAXDIST; if (!any(keep)) next
      b <- b[keep]; dist <- dist[keep]
      r <- suppressWarnings(cor(Gc[,a], Gc[,b,drop=FALSE], use="pairwise.complete.obs"))
      out[[length(out)+1]] <- data.table(dist=dist, r2=as.numeric(r)^2,
                                          recomb=rr[a], minDI=pmin(dd[a], dd[b]))
    }
  }
  rbindlist(out)[is.finite(r2)]
}

cat("\ncomputing LD pairs on superset per panel ...\n")
pairs <- rbindlist(lapply(names(panels), function(nm){
  d <- ld_pairs_panel(panels[[nm]]); d[, panel:=nm][]
}))
pairs[, rec_stratum := fifelse(recomb <= rec_split, "low recomb", "high recomb")]

## ---- summarise r^2 by threshold (short-range <=20kb, where within-sp LD peaks) ----
sweep <- rbindlist(lapply(THRESHOLDS, function(t){
  s <- pairs[minDI > t & dist <= 20000,
             .(mean_r2 = mean(r2), n_pairs=.N), by=.(panel, rec_stratum)]
  s[, threshold := t][]
}))
wide <- dcast(sweep, threshold + rec_stratum ~ panel, value.var="mean_r2")
wide <- poly_tab[wide, on="threshold"]
wide[, within_over_pooled := round(pmax(aqu,pol)/pooled, 2)]
wide[, `:=`(aqu=round(aqu,3), pol=round(pol,3), pooled=round(pooled,3))]
setcolorder(wide, c("threshold","n_markers","poly_aqu","poly_pol","rec_stratum",
                    "aqu","pol","pooled","within_over_pooled"))
cat("\n=== short-range (<=20kb) mean composite r^2 by threshold x stratum ===\n")
print(wide[order(rec_stratum, -threshold)])

## ---- plot: the tradeoff ----
plt <- melt(sweep[rec_stratum=="low recomb"],
            id.vars=c("threshold","panel"), measure.vars="mean_r2")
plt <- plt[panel %in% c("aqu","pol","pooled")]
plt[, panel := factor(panel, levels=c("pooled","pol","aqu"),
                      labels=c("pooled 50/50 (admixture)","polyctena (within-sp)","aquilonia (within-sp)"))]
mk <- poly_tab[, .(threshold, n_markers)]

p1 <- ggplot(plt, aes(factor(threshold), value, fill=panel)) +
  geom_col(position="dodge") +
  geom_text(data=mk, aes(factor(threshold), y=1.02, label=paste0(n_markers,"\nmarkers")),
            inherit.aes=FALSE, size=3, vjust=0) +
  scale_fill_manual(values=c("#d95f02","#7570b3","#1b9e77"), name=NULL) +
  ylim(0,1.15) +
  labs(x="DI threshold", y=expression("short-range mean "*r^2*" (low recomb, <=20kb)"),
       title="DI-threshold tradeoff: within-species LD vs marker density",
       subtitle="stricter DI shrinks within-species LD (the null 'problem') but costs markers") +
  theme_bw(base_size=13) + theme(legend.position="bottom", panel.grid.minor=element_blank())

dir.create("Figures", showWarnings=FALSE)
ggsave("Figures/parent_LD_DI_sweep.png", p1, width=9, height=5.5, dpi=110)
cat("\nsaved Figures/parent_LD_DI_sweep.png\n")
