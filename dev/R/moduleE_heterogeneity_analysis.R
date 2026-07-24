## =============================================================================
## Module E, E2g-3 -- founder-heterogeneity analysis
## =============================================================================
## Does founding demes from DIFFERENT local sub-pools (k>1) raise among-deme F_ST
## toward the empirical 0.313 -- while holding the LD-decay match and with founder
## NUMBER held fixed (so the effect is source heterogeneity, not bottleneck size)?
##
##   k = 1 : rangewide control (all demes from the same pool)  <- current behaviour
##   k > 1 : each deme founded from one of k sub-pools
##
## VERDICT LOGIC
##   - if F_ST reaches ~0.313 at some (k, K, gen) with LD_rmse < 0.04 and pi near
##     0.299 -> founder heterogeneity gives a NEUTRAL explanation; the excess
##     ancestry sorting is not evidence of selection.
##   - if F_ST plateaus well below 0.313 (and pi stays high) -> neutral drift +
##     founder heterogeneity together still cannot reproduce the data, and the
##     excess-sorting (selection/DMI) conclusion stands.
##
## Run AFTER run_founder_heterogeneity.sh. Self-contained (helpers as in E2e/E2f).
## Output: data/moduleE_sim/moduleE_heterogeneity_results.rds
##         Figures/moduleE_heterogeneity.pdf
## =============================================================================

suppressMessages({library(data.table); library(ggplot2); library(parallel)})

FORMICA   <- "/Users/petrikem/gitlab/formica_hybrid"
HETDIR    <- file.path(FORMICA, "data/moduleE_sim/founder_het")
FVCF_DIR  <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/SLiM/founders_DIm25_thin15000"
EMP_RDATA <- file.path(FORMICA, "data/hybrids_only_maf005.Rdata")
OUTRDS    <- file.path(FORMICA, "data/moduleE_sim/moduleE_heterogeneity_results.rds")
FIGDIR    <- file.path(FORMICA, "Figures")

HET_KS <- c(1, 2, 3, 4)          # heterogeneity levels present in the runs
KS     <- c(375, 750)
GENS   <- seq(20, 160, by = 20)
CORES  <- 8L; MATCH_SEED <- 1L
MAXDIST <- 250000L; MAXPART <- 100L
DIST_BREAKS <- c(0,1,2,5,10,20,50,100,150,250) * 1000L
MIN_MAF <- 0.05
dir.create(FIGDIR, showWarnings = FALSE)

## ------------------------------ helpers -------------------------------------
read_vcf_dosage <- function(vcf) {
  dt <- fread(vcf, skip = "#CHROM", header = TRUE, sep = "\t", showProgress = FALSE)
  chr <- sub("ch", "Chr", dt[[1]]); pos <- as.integer(dt[[2]])
  G <- as.matrix(dt[, 10:ncol(dt)]); gt <- sub(":.*$", "", G)
  a1 <- suppressWarnings(as.integer(sub("[|/].*$", "", gt)))
  a2 <- suppressWarnings(as.integer(sub("^.*[|/]", "", gt)))
  list(dose = matrix(a1 + a2, nrow = nrow(G)), Chr = chr, Pos = pos,
       marker = paste(chr, pos, sep = ":"))
}
ld_decay_curve <- function(dose, Chr, Pos) {
  out <- vector("list", 0)
  for (ch in unique(Chr)) {
    j <- which(Chr == ch); if (length(j) < 2L) next
    o <- order(Pos[j]); j <- j[o]; pos <- Pos[j]; G <- dose[j, , drop = FALSE]
    p <- rowMeans(G, na.rm = TRUE)/2
    jj <- which(is.finite(p) & pmin(p, 1-p) >= MIN_MAF); if (length(jj) < 2L) next
    pos <- pos[jj]; G <- G[jj, , drop = FALSE]
    for (a in seq_len(length(pos)-1L)) {
      b <- (a+1L):min(a+MAXPART, length(pos)); d <- pos[b]-pos[a]
      keep <- d <= MAXDIST; if (!any(keep)) next; b <- b[keep]; d <- d[keep]
      r <- suppressWarnings(cor(G[a,], t(G[b,,drop=FALSE]), use="pairwise.complete.obs"))
      out[[length(out)+1L]] <- data.table(dist=d, r2=as.numeric(r)^2)
    }
  }
  if (!length(out)) return(NULL)
  res <- rbindlist(out)[is.finite(r2)]
  res[, dbin := cut(dist, DIST_BREAKS, labels=FALSE)]
  res[!is.na(dbin), .(r2=mean(r2), dmid=(DIST_BREAKS[dbin]+DIST_BREAKS[dbin+1L])/2), by=dbin]
}
pi_hat  <- function(dose) { p <- rowMeans(dose, na.rm=TRUE)/2; mean(2*p*(1-p), na.rm=TRUE) }
fst_hat <- function(freqs) {
  pbar <- colMeans(freqs, na.rm=TRUE); Hs <- colMeans(2*freqs*(1-freqs), na.rm=TRUE)
  Ht <- 2*pbar*(1-pbar); sum(Ht-Hs, na.rm=TRUE)/sum(Ht, na.rm=TRUE)
}

## ------------------- 1. empirical targets (as E2e/E2f) ----------------------
message("[1] empirical ...")
e <- new.env(); load(EMP_RDATA, envir = e)
GTh <- e$GTs_hybrids_005; sdh <- as.data.table(e$sample_data); maph <- as.data.table(e$map_hyb_005)
fmark <- unique(rbindlist(lapply(list.files(FVCF_DIR, "founders_ch.*\\.vcf$", full.names=TRUE),
  function(f) fread(f, skip="#CHROM", select=1:3, header=TRUE))))
ci <- match(fmark$ID, maph$marker); ci <- ci[!is.na(ci)]
mk_chr <- maph$Chr[ci]; mk_pos <- maph$Pos[ci]
hyb_pops <- setdiff(unique(sdh$Population), c("aquilonia_parent","polyctena_parent"))
emp_ld <- rbindlist(lapply(hyb_pops, function(pp)
  ld_decay_curve(t(GTh[sdh$Population==pp, ci, drop=FALSE]), mk_chr, mk_pos)))[
  , .(r2=mean(r2), dmid=dmid[1]), by=dbin][order(dbin)]
emp_pi <- mean(sapply(hyb_pops, function(pp) pi_hat(t(GTh[sdh$Population==pp, ci, drop=FALSE]))))
emp_fst <- fst_hat(do.call(rbind, lapply(hyb_pops, function(pp)
  rowMeans(t(GTh[sdh$Population==pp, ci, drop=FALSE]), na.rm=TRUE)/2)))
emp_ns <- sapply(hyb_pops, function(pp) sum(sdh$Population==pp))
emp_vec <- setNames(emp_ld$r2, emp_ld$dbin)
message(sprintf("    empirical pi=%.4f  F_ST=%.4f", emp_pi, emp_fst))

## ------------------- 2. cells: (heterogeneity k, K, generation) -------------
message("[2] heterogeneity cells ...")
cell <- function(hk, K, gen) {
  pat <- sprintf("hyb_hetk%d_Naq[0-9]+_Npol[0-9]+_K%d_rep[0-9]+_gen%d\\.vcf$", hk, K, gen)
  vcfs <- list.files(HETDIR, pat, full.names = TRUE); if (!length(vcfs)) return(NULL)
  per <- lapply(seq_along(vcfs), function(i) {
    x <- read_vcf_dosage(vcfs[i]); dose <- x$dose
    tgt <- emp_ns[(i-1L) %% length(emp_ns) + 1L]; ni <- ncol(dose)
    if (ni > tgt) { set.seed(MATCH_SEED + i); dose <- dose[, sample(ni, tgt), drop=FALSE] }
    list(curve=ld_decay_curve(dose, x$Chr, x$Pos), pi=pi_hat(dose),
         freq=setNames(rowMeans(dose, na.rm=TRUE)/2, x$marker))
  })
  curves <- rbindlist(lapply(seq_along(per), function(i)
    if (!is.null(per[[i]]$curve)) per[[i]]$curve[, deme := i][]))
  ld <- curves[, .(r2=mean(r2), dmid=dmid[1]), by=dbin][order(dbin)]
  common <- Reduce(intersect, lapply(per, function(z) names(z$freq)))
  fmat <- do.call(rbind, lapply(per, function(z) z$freq[common]))
  m <- match(ld$dbin, as.integer(names(emp_vec)))
  data.table(het_k=hk, K=K, gen=gen, n_demes=length(vcfs),
             ld_rmse=sqrt(mean((ld$r2-emp_vec[m])^2, na.rm=TRUE)),
             pi=mean(sapply(per,`[[`,"pi")), fst=fst_hat(fmat), ld=list(ld))
}
cells <- CJ(hk=HET_KS, K=KS, gen=GENS, sorted=FALSE)
sim <- rbindlist(mcMap(cell, cells$hk, cells$K, cells$gen, mc.cores=CORES), fill=TRUE)
if (!nrow(sim)) stop("no heterogeneity runs found in ", HETDIR)

saveRDS(list(sim=sim, emp_pi=emp_pi, emp_fst=emp_fst, emp_ld=emp_ld), OUTRDS)

## ------------------- 3. figure + verdict ------------------------------------
mdt <- melt(sim[, .(het_k, K, gen, ld_rmse, pi, fst)], id.vars=c("het_k","K","gen"))
ref <- data.table(variable=c("pi","fst"), y=c(emp_pi, emp_fst))
p <- ggplot(mdt, aes(gen, value, colour=factor(het_k))) +
  geom_line() + geom_point(size=1.2) +
  geom_hline(data=ref, aes(yintercept=y), linetype=2) +
  facet_grid(variable ~ K, scales="free_y") +
  scale_colour_viridis_d(name="sub-pools (k)\n1 = rangewide", option="viridis", end=0.85) +
  labs(x="generation", y=NULL,
       title="Module E: founder heterogeneity -- does local founding reach empirical F_ST?",
       subtitle="founder number held fixed across k; dashed = empirical pi / F_ST") +
  theme_bw(base_size=11)
ggsave(file.path(FIGDIR, "moduleE_heterogeneity.pdf"), p, width=10, height=6)

cat(sprintf("\n=== empirical: pi=%.3f  F_ST=%.3f ===\n", emp_pi, emp_fst))
cat("\n=== max F_ST reached per heterogeneity level (LD_rmse < 0.04) ===\n")
ok <- sim[ld_rmse < 0.04]
if (nrow(ok)) {
  print(ok[, .SD[which.max(fst)], by=het_k][order(het_k),
        .(het_k, K, gen, ld_rmse=round(ld_rmse,4), pi=round(pi,3),
          fst=round(fst,3), fst_gap=round(fst-emp_fst,3))])
  cat("\nVERDICT: if fst_gap approaches 0 as k rises, founder heterogeneity gives a\n")
  cat("neutral explanation; if it plateaus well short, the excess sorting stands.\n")
} else cat("(no cells with LD_rmse < 0.04)\n")
cat("\nsaved: ", OUTRDS, " ; ", file.path(FIGDIR, "moduleE_heterogeneity.pdf"), "\n")
