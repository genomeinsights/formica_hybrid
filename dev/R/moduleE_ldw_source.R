## =============================================================================
## Module E -- is the excess local LD founder-inherited or drift/Wahlund?
## =============================================================================
## Compares WITHIN-group local LD (removes the between-population / Wahlund
## component) across:
##   founder pool (aq, pol)  -- the LD the founders carry
##   empirical  within-population (avg over pops)
##   simulated  within-deme (avg over demes), gen20 and gen160
## plus the POOLED sim/empirical for reference. All N-matched to 8 diploids
## (16 gametes) so the finite-sample r^2 bias is identical. Reported per DI bin.
##
## Read:
##  - if within-deme sim ~ within-pop empirical, the pooled excess was Wahlund
##    (between-deme differentiation -> Ne too small -> higher K);
##  - if within-deme sim stays well above empirical, the excess is founder-inherited
##    (rangewide founders carry too much within-species LD).
## =============================================================================

suppressMessages({ library(data.table) })

FORMICA  <- "/Users/petrikem/gitlab/formica_hybrid"
SIMDIR   <- file.path(FORMICA, "data/moduleE_sim/distrat")
POOL     <- file.path(FORMICA, "data/moduleE_founders/moduleE_founder_haplotypes.rds")
FVCF_DIR <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/SLiM/founders_maf015_DIstrat4000"
EMP_RD   <- file.path(FORMICA, "data/hybrids_only_maf005.Rdata")
OUTRDS   <- file.path(FORMICA, "data/moduleE_sim/moduleE_ldw_source.rds")

TAG <- "Naq30_Npol13"; K <- 6250; GENS <- c(20, 160)
NGAM <- 16L                                   # common gamete count = 8 diploids
DI_BREAKS <- c(-Inf,-90,-75,-60,-50,-40,-30,-25,-20,-15,Inf)
WIN_BP <- 100000L; MAXPART <- 100L; MAF_MIN <- 0.1
NPOP_MAX <- 12L; NDEME_MAX <- 12L             # cap #groups for runtime

## per-SNP local LD (dosage; >=0). geno: gametes-or-diploids x markers dosage
local_ld <- function(geno, Chr, Pos, ploidy=2) {
  res <- rep(NA_real_, ncol(geno))
  for (ch in unique(Chr)) {
    j <- which(Chr == ch); o <- order(Pos[j]); j <- j[o]
    pos <- Pos[j]; G <- geno[, j, drop=FALSE]
    p <- colMeans(G, na.rm=TRUE)/ploidy; poly <- is.finite(p) & pmin(p,1-p) >= MAF_MIN
    for (a in seq_along(j)) {
      if (!poly[a]) next
      b <- which(pos <= pos[a] + WIN_BP & seq_along(pos) > a & poly)
      if (length(b) > MAXPART) b <- b[seq_len(MAXPART)]
      if (!length(b)) next
      r <- suppressWarnings(cor(G[,a], G[, b, drop=FALSE], use="pairwise.complete.obs"))
      res[j[a]] <- median(as.numeric(r)^2, na.rm=TRUE)
    }
  }
  pmax(res, 0)
}

## ---- markers + DI ----------------------------------------------------------
fm <- rbindlist(lapply(list.files(FVCF_DIR, "founders_ch.*[.]vcf$", full.names=TRUE),
  function(f) fread(f, skip="#CHROM", select=c(1,2,3), header=TRUE, col.names=c("chs","Pos","marker"))))
fm[, `:=`(Chr=sub("ch","Chr",chs), chr_id=as.integer(sub("ch","",chs)))]
setorder(fm, chr_id, Pos); markers <- fm$marker
ph <- readRDS(POOL); pmap <- as.data.table(ph$map)
DI <- pmap$DI[match(markers, pmap$marker)]; dibin <- cut(DI, DI_BREAKS)

## helper: subsample a diploid-dosage matrix (inds x markers) to n gametes/2 inds, local_ld
ldw_group <- function(dose) {                  # dose: inds x markers (diploid 0/1/2)
  n <- NGAM %/% 2L
  if (nrow(dose) > n) { set.seed(1); dose <- dose[sample(nrow(dose), n), , drop=FALSE] }
  local_ld(dose, fm$Chr, fm$Pos, ploidy=2)
}
read_dose <- function(f) {
  dt <- fread(f, skip="#CHROM", header=TRUE, sep="\t", showProgress=FALSE)
  mk <- paste(sub("ch","Chr",dt[[1]]), dt[[2]], sep=":")
  G <- as.matrix(dt[,10:ncol(dt)]); gt <- sub(":.*$","",G)
  d <- matrix(suppressWarnings(as.integer(sub("[|/].*$","",gt))) +
              suppressWarnings(as.integer(sub("^.*[|/]","",gt))), nrow=nrow(G))
  m <- match(markers, mk); keep <- !is.na(m)
  D <- matrix(0L, length(markers), ncol(d)); D[keep,] <- d[m[keep],,drop=FALSE]; t(D)
}
mean_by_bin <- function(w) data.table(dibin=dibin, w=w)[is.finite(w), .(ld=mean(w)), by=dibin][order(dibin)]

## ---- founder pool within-species -------------------------------------------
message("[1] founders ...")
ci <- match(markers, pmap$marker)
Haq <- ph$aqu[, ci, drop=FALSE]; Hpol <- ph$pol[, ci, drop=FALSE]
## pair haplotypes -> pseudo-diploid dosage (random pairing preserves LD)
set.seed(1)
aq_pair  <- { o <- sample(nrow(Haq)); Haq[o[c(TRUE,FALSE)],,drop=FALSE] + Haq[o[c(FALSE,TRUE)],,drop=FALSE] }
pol_pair <- Hpol[seq(1,nrow(Hpol),2),,drop=FALSE] + Hpol[seq(2,nrow(Hpol),2),,drop=FALSE]
f_aq  <- mean_by_bin(ldw_group(aq_pair))[,  src:="founder aq"]
f_pol <- mean_by_bin(ldw_group(pol_pair))[, src:="founder pol"]

## ---- empirical within-population + pooled ----------------------------------
message("[2] empirical ...")
e <- new.env(); load(EMP_RD, envir=e)
GTh <- e$GTs_hybrids_005; sdh <- as.data.table(e$sample_data); maph <- as.data.table(e$map_hyb_005)
eci <- match(markers, maph$marker)
pops <- setdiff(unique(sdh$Population), c("aquilonia_parent","polyctena_parent"))
psz  <- sapply(pops, function(p) sum(sdh$Population==p))
usepops <- names(sort(psz, decreasing=TRUE))[1:min(NPOP_MAX, length(pops))]
emp_within <- rbindlist(lapply(usepops, function(p){
  g <- GTh[sdh$Population==p, eci, drop=FALSE]; g[is.na(g)] <- 0L
  data.table(marker=markers, w=ldw_group(g)) }))[, .(w=mean(w,na.rm=TRUE)), by=marker]
emp_w <- mean_by_bin(emp_within$w[match(markers, emp_within$marker)])[, src:="emp within-pop"]

## ---- simulated within-deme (gen20, gen160) ---------------------------------
sim_w <- rbindlist(lapply(GENS, function(g){
  message("[3] sim within-deme gen", g, " ...")
  fs <- file.path(SIMDIR, sprintf("hyb_neutral_realfounders_%s_K%d_rep%d_gen%d.vcf", TAG, K, 1:NDEME_MAX, g))
  fs <- fs[file.exists(fs)]
  wd <- rbindlist(lapply(fs, function(f) data.table(marker=markers, w=ldw_group(read_dose(f)))))[
        , .(w=mean(w,na.rm=TRUE)), by=marker]
  mean_by_bin(wd$w[match(markers, wd$marker)])[, src:=sprintf("sim within-deme gen%d", g)]
}))

out <- rbind(f_aq, f_pol, emp_w, sim_w)
saveRDS(out, OUTRDS)

## ---- report ----------------------------------------------------------------
w <- dcast(out, dibin ~ src, value.var="ld")
setcolorder(w, c("dibin","founder aq","founder pol","emp within-pop",
                 sprintf("sim within-deme gen%d", GENS)))
cat(sprintf("\n=== WITHIN-group mean local LD by DI bin (N-matched to %d gametes) ===\n", NGAM))
print(w[, lapply(.SD, function(x) if(is.numeric(x)) round(x,3) else x)])
neu <- w[1]
cat(sprintf("\nNEUTRAL bin (DI<=-90):  founder_aq=%.3f  founder_pol=%.3f  emp_within=%.3f  sim_within_gen20=%.3f  sim_within_gen160=%.3f\n",
   neu[["founder aq"]], neu[["founder pol"]], neu[["emp within-pop"]],
   neu[["sim within-deme gen20"]], neu[["sim within-deme gen160"]]))
cat("\nREAD: sim_within ~ emp_within  => pooled excess was Wahlund/drift (Ne).\n")
cat("      sim_within >> emp_within   => founder-inherited (rangewide founders too structured).\n")
cat("      founder_* vs emp_within tells you how much LD the founders start with.\n")
cat("\nsaved: ", OUTRDS, "\n")
