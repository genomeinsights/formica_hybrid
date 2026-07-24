## =============================================================================
## Module E, E2f -- K-sweep analysis: can stronger-drift NEUTRAL demography match
##                 empirical pi/F_ST while holding the LD-decay match?
## =============================================================================
## For foundings that matched empirical LD-decay (30/13, 12/12), compute LD-decay
## RMSE, pi and among-pop F_ST as a function of K over {375,750,1500,3000,6250}
## (K<6250 from the ksweep; K=6250 from the screen). Sample-size matched to the
## empirical pops (as in E2e). The question:
##   - does F_ST rise toward empirical (0.31) and pi fall toward empirical (0.30)
##     as K decreases, AND does the LD-decay match hold?
##   - if some K reproduces LD + pi + F_ST jointly -> neutral drift explains the
##     data; if F_ST never reaches ~0.31 at an acceptable LD-RMSE -> the excess
##     ancestry sorting is not explained by neutral drift (selection lead).
##
## Self-contained (helpers + empirical recomputed identically to E2e).
## Output: data/moduleE_sim/moduleE_ksweep_results.rds
##         Figures/moduleE_ksweep.pdf
## =============================================================================

suppressMessages({library(data.table); library(ggplot2); library(parallel)})

## ------------------------------- CONFIG -------------------------------------
FORMICA   <- "/Users/petrikem/gitlab/formica_hybrid"
SCREEN    <- file.path(FORMICA, "data/moduleE_sim/screen")     # K=6250 cells (K-less names)
KSWEEP    <- file.path(FORMICA, "data/moduleE_sim/ksweep")     # K<6250 (K-inclusive names)
FVCF_DIR  <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/SLiM/founders_DIm25_thin15000"
EMP_RDATA <- file.path(FORMICA, "data/hybrids_only_maf005.Rdata")
OUTRDS    <- file.path(FORMICA, "data/moduleE_sim/moduleE_ksweep_results.rds")
FIGDIR    <- file.path(FORMICA, "Figures")

FOUNDINGS <- data.table(setting = c("30/13","12/12"), tag = c("Naq30_Npol13","Naq12_Npol12"))
KS   <- c(50, 100, 200, 375, 750, 1500, 3000, 6250)   # extended DOWNWARD: the
## earlier sweep stopped at 375 with F_ST still climbing (0.12 -> 0.20) and not
## saturating. Reaching the empirical F_ST (0.313) implies Ne ~ 50-135 from
## F_ST ~ 1-exp(-t/2Ne), i.e. below anything previously simulated.
GENS <- seq(20, 160, by = 20)
CORES <- 8L
MATCH_SEED <- 1L

MAXDIST <- 250000L; MAXPART <- 100L
DIST_BREAKS <- c(0,1,2,5,10,20,50,100,150,250) * 1000L
MIN_MAF <- 0.05
dir.create(FIGDIR, showWarnings = FALSE)

## ------------------------------ helpers (as E2e) ----------------------------
read_vcf_dosage <- function(vcf) {
  dt <- fread(vcf, skip = "#CHROM", header = TRUE, sep = "\t", showProgress = FALSE)
  chr <- sub("ch", "Chr", dt[[1]]); pos <- as.integer(dt[[2]])
  G   <- as.matrix(dt[, 10:ncol(dt)])
  gt  <- sub(":.*$", "", G)
  a1  <- suppressWarnings(as.integer(sub("[|/].*$", "", gt)))
  a2  <- suppressWarnings(as.integer(sub("^.*[|/]", "", gt)))
  list(dose = matrix(a1 + a2, nrow = nrow(G)), Chr = chr, Pos = pos,
       marker = paste(chr, pos, sep = ":"))
}
ld_decay_curve <- function(dose, Chr, Pos) {
  out <- vector("list", 0)
  for (ch in unique(Chr)) {
    j <- which(Chr == ch); if (length(j) < 2L) next
    o <- order(Pos[j]); j <- j[o]; pos <- Pos[j]; G <- dose[j, , drop = FALSE]
    p <- rowMeans(G, na.rm = TRUE) / 2
    jj <- which(is.finite(p) & pmin(p, 1 - p) >= MIN_MAF); if (length(jj) < 2L) next
    pos <- pos[jj]; G <- G[jj, , drop = FALSE]
    for (a in seq_len(length(pos) - 1L)) {
      b <- (a + 1L):min(a + MAXPART, length(pos)); d <- pos[b] - pos[a]
      keep <- d <= MAXDIST; if (!any(keep)) next; b <- b[keep]; d <- d[keep]
      r <- suppressWarnings(cor(G[a, ], t(G[b, , drop = FALSE]), use = "pairwise.complete.obs"))
      out[[length(out) + 1L]] <- data.table(dist = d, r2 = as.numeric(r)^2)
    }
  }
  if (!length(out)) return(NULL)
  res <- rbindlist(out)[is.finite(r2)]
  res[, dbin := cut(dist, DIST_BREAKS, labels = FALSE)]
  res[!is.na(dbin), .(r2 = mean(r2), dmid = (DIST_BREAKS[dbin] + DIST_BREAKS[dbin+1L])/2), by = dbin]
}
pi_hat  <- function(dose) { p <- rowMeans(dose, na.rm = TRUE)/2; mean(2*p*(1-p), na.rm = TRUE) }
fst_hat <- function(freqs) {
  pbar <- colMeans(freqs, na.rm = TRUE); Hs <- colMeans(2*freqs*(1-freqs), na.rm = TRUE)
  Ht <- 2*pbar*(1-pbar); sum(Ht - Hs, na.rm = TRUE)/sum(Ht, na.rm = TRUE)
}

## ------------------- 1. empirical targets (identical to E2e) ----------------
message("[1] empirical ...")
e <- new.env(); load(EMP_RDATA, envir = e)
GTh <- e$GTs_hybrids_005; sdh <- as.data.table(e$sample_data); maph <- as.data.table(e$map_hyb_005)
fmark <- unique(rbindlist(lapply(list.files(FVCF_DIR, "founders_ch.*\\.vcf$", full.names = TRUE),
  function(f) fread(f, skip = "#CHROM", select = 1:3, header = TRUE))))
ci <- match(fmark$ID, maph$marker); ci <- ci[!is.na(ci)]
mk_chr <- maph$Chr[ci]; mk_pos <- maph$Pos[ci]
hyb_pops <- setdiff(unique(sdh$Population), c("aquilonia_parent","polyctena_parent"))
emp_curves <- rbindlist(lapply(hyb_pops, function(pp) {
  cur <- ld_decay_curve(t(GTh[sdh$Population == pp, ci, drop = FALSE]), mk_chr, mk_pos)
  if (!is.null(cur)) cur[, pop := pp][] }))
emp_ld <- emp_curves[, .(r2 = mean(r2), dmid = dmid[1]), by = dbin][order(dbin)]
emp_pi <- mean(sapply(hyb_pops, function(pp) pi_hat(t(GTh[sdh$Population == pp, ci, drop = FALSE]))))
emp_fst <- fst_hat(do.call(rbind, lapply(hyb_pops, function(pp)
  rowMeans(t(GTh[sdh$Population == pp, ci, drop = FALSE]), na.rm = TRUE)/2)))
emp_ns <- sapply(hyb_pops, function(pp) sum(sdh$Population == pp))
emp_vec <- setNames(emp_ld$r2, emp_ld$dbin)
message(sprintf("    empirical pi=%.4f  F_ST=%.4f", emp_pi, emp_fst))

## ------------------- 2. sim per (founding, K, gen), sample-size matched -----
message("[2] K cells ...")
sim_cell <- function(tag, K, gen) {
  if (K == 6250) {
    pat <- sprintf("hyb_neutral_realfounders_%s_rep[0-9]+_gen%d\\.vcf$", tag, gen); dir <- SCREEN
  } else {
    pat <- sprintf("hyb_neutral_realfounders_%s_K%d_rep[0-9]+_gen%d\\.vcf$", tag, K, gen); dir <- KSWEEP
  }
  vcfs <- list.files(dir, pat, full.names = TRUE); if (!length(vcfs)) return(NULL)
  per <- lapply(seq_along(vcfs), function(k) {
    x <- read_vcf_dosage(vcfs[k]); dose <- x$dose
    tgt <- emp_ns[(k-1L) %% length(emp_ns) + 1L]; ni <- ncol(dose)
    if (ni > tgt) { set.seed(MATCH_SEED + k); dose <- dose[, sample(ni, tgt), drop = FALSE] }
    list(curve = ld_decay_curve(dose, x$Chr, x$Pos), pi = pi_hat(dose),
         freq = setNames(rowMeans(dose, na.rm = TRUE)/2, x$marker))
  })
  curves <- rbindlist(lapply(seq_along(per), function(i)
    if (!is.null(per[[i]]$curve)) per[[i]]$curve[, deme := i][]))
  ld <- curves[, .(r2 = mean(r2), dmid = dmid[1]), by = dbin][order(dbin)]
  common <- Reduce(intersect, lapply(per, function(z) names(z$freq)))
  fmat <- do.call(rbind, lapply(per, function(z) z$freq[common]))
  m <- match(ld$dbin, as.integer(names(emp_vec)))
  data.table(tag = tag, K = K, gen = gen, n_demes = length(vcfs),
             ld_rmse = sqrt(mean((ld$r2 - emp_vec[m])^2, na.rm = TRUE)),
             pi = mean(sapply(per, `[[`, "pi")), fst = fst_hat(fmat), ld = list(ld))
}
cells <- CJ(tag = FOUNDINGS$tag, K = KS, gen = GENS, sorted = FALSE)
sim <- rbindlist(mcMap(sim_cell, cells$tag, cells$K, cells$gen, mc.cores = CORES), fill = TRUE)
sim <- FOUNDINGS[sim, on = "tag"]

saveRDS(list(sim = sim, emp_ld = emp_ld, emp_pi = emp_pi, emp_fst = emp_fst,
             emp_ns = emp_ns, n_markers = length(ci)), OUTRDS)

## ------------------- 3. figure: pi / F_ST / LD-RMSE vs K --------------------
## use the best-LD generation per (founding,K) for the pi/F_ST vs K summary
bestg <- sim[, .SD[which.min(ld_rmse)], by = .(setting, K)]
mdt <- melt(bestg[, .(setting, K, ld_rmse, pi, fst)], id.vars = c("setting","K"))
ref <- data.table(variable = c("pi","fst","ld_rmse"), y = c(emp_pi, emp_fst, 0))
p <- ggplot(mdt, aes(K, value, colour = setting)) +
  geom_line() + geom_point(size = 2) + scale_x_log10() +
  geom_hline(data = ref[variable != "ld_rmse"], aes(yintercept = y), linetype = 2) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(x = "K (log)", y = NULL,
       title = "Module E K-sweep: pi / among-pop F_ST / LD-decay RMSE vs K",
       subtitle = "at best-LD generation per K; dashed = empirical target (pi, F_ST)") +
  theme_bw(base_size = 11)
ggsave(file.path(FIGDIR, "moduleE_ksweep.pdf"), p, width = 11, height = 4.5)

## ------------------- 4. verdict table --------------------------------------
cat(sprintf("\n=== empirical: pi=%.3f  F_ST=%.3f ===\n", emp_pi, emp_fst))
cat("\n=== best-LD generation per (founding, K) ===\n")
print(bestg[order(setting, K), .(setting, K, best_gen = gen, n = n_demes,
      ld_rmse = round(ld_rmse,4), pi = round(pi,3), fst = round(fst,3))])
## cells that BOTH match LD well AND approach empirical F_ST
cat("\n=== cells with LD_rmse<0.04 ranked by |F_ST - empirical| ===\n")
cand <- sim[ld_rmse < 0.04][order(abs(fst - emp_fst))]
print(head(cand[, .(setting, K, gen, ld_rmse = round(ld_rmse,4),
      pi = round(pi,3), fst = round(fst,3), fst_gap = round(fst - emp_fst,3))], 12))
cat("\nsaved: ", OUTRDS, " ; ", file.path(FIGDIR, "moduleE_ksweep.pdf"), "\n")
