## =============================================================================
## Module E, Script E2e -- founder-number screen: LD-decay / pi / F_ST vs empirical
## =============================================================================
##
## THE GO/NO-GO GATE. For each founder setting (N_AQ,N_POL) and each sampled
## generation (20..160), compute -- on the SAME thinned DI>-25 marker set the sim
## carries -- three demographic summaries and compare them to the EMPIRICAL 20
## hybrid populations:
##   - LD-decay  : composite r^2 vs physical distance, within population (the
##                 primary target; sets how much admixture LD has broken down)
##   - pi        : mean within-population expected heterozygosity
##   - F_ST      : differentiation AMONG the replicate populations/demes
## calibrated on a SORTING-NEUTRAL summary set (LD/pi/FST), with the sorting
## statistics deliberately HELD OUT.
##
## Output: which (founder-N, generation) best reproduces the empirical LD-decay
## (and pi/FST), and whether the conditioned-founder neutral model CAN reproduce
## it at all -- the license for the whole null. Nothing here touches sort_class.
##
## Reads whatever demes are present in the screen dir (robust to the 3 re-running
## small-founder demes not yet being done).
##
## Outputs: data/moduleE_sim/moduleE_screen_results.rds
##          Figures/moduleE_screen_{lddecay,pi_fst}.pdf
## =============================================================================

suppressMessages({library(data.table); library(ggplot2); library(parallel)})

## ------------------------------- CONFIG -------------------------------------
FORMICA   <- "/Users/petrikem/gitlab/formica_hybrid"
SCREEN    <- file.path(FORMICA, "data/moduleE_sim/screen")
FVCF_DIR  <- "/Users/petrikem/gitlab/Replicate-hybrid-evolution/SLiM/founders_DIm25_thin15000"
EMP_RDATA <- file.path(FORMICA, "data/hybrids_only_maf005.Rdata")   # GTs_hybrids_005, map_hyb_005, sample_data
OUTRDS    <- file.path(FORMICA, "data/moduleE_sim/moduleE_screen_results.rds")
FIGDIR    <- file.path(FORMICA, "Figures")

## founder settings (label -> N_AQ,N_POL) and sampled generations
SETTINGS <- data.table(
  setting = c("4/4","6/6","8/8","12/12","20/13","30/13"),
  tag     = c("Naq4_Npol4","Naq6_Npol6","Naq8_Npol8","Naq12_Npol12","Naq20_Npol13","Naq30_Npol13"))
GENS  <- seq(20, 160, by = 20)
CORES <- 8L

## SAMPLE-SIZE MATCHING. Empirical hybrid pops are 3-20 diploids; the sim samples
## 50. Composite r^2 carries a finite-sample bias (~1/2n) that inflates the
## small-n empirical estimates, spuriously matching earlier (higher-LD) sim
## generations. Fix: subsample each sim deme to an empirical pop's actual size so
## the per-unit bias structure is identical on both sides (empirical kept at real
## n; the 20 demes take the 20 empirical sizes as a set). Reproducible seed.
MATCH_N    <- TRUE
MATCH_SEED <- 1L

## smoke-test hook: MODULEE_E2E_TEST=1 -> 2 settings, 2 gens, 3 empirical pops
TEST <- nzchar(Sys.getenv("MODULEE_E2E_TEST"))
if (TEST) { SETTINGS <- SETTINGS[c(1,4)]; GENS <- c(20, 160); message("[TEST MODE]") }

## LD-decay parameters (composite r^2 vs bp distance, within chromosome)
MAXDIST     <- 250000L                                   # max pair distance (bp)
MAXPART     <- 100L                                      # cap downstream partners/marker
DIST_BREAKS <- c(0,1,2,5,10,20,50,100,150,250) * 1000L   # bp bins
MIN_MAF     <- 0.05                                      # drop within-panel monomorphic for LD

dir.create(FIGDIR, showWarnings = FALSE)

## ------------------------------ helpers -------------------------------------

## read a (multi-chrom) VCF -> dosage matrix (markers x inds) + Chr/Pos.
## dosage = count of ALT allele; robust to GT subfields and missing "."
read_vcf_dosage <- function(vcf) {
  dt <- fread(vcf, skip = "#CHROM", header = TRUE, sep = "\t", showProgress = FALSE)
  chr <- sub("ch", "Chr", dt[[1]]); pos <- as.integer(dt[[2]])
  G   <- as.matrix(dt[, 10:ncol(dt)])
  gt  <- sub(":.*$", "", G)
  a1  <- suppressWarnings(as.integer(sub("[|/].*$", "", gt)))
  a2  <- suppressWarnings(as.integer(sub("^.*[|/]", "", gt)))
  d   <- matrix(a1 + a2, nrow = nrow(G))                 # markers x inds, 0/1/2/NA
  list(dose = d, Chr = chr, Pos = pos, marker = paste(chr, pos, sep = ":"))
}

## composite-r^2 LD decay for one population's dosage matrix (markers x inds).
## returns per-distance-bin mean r^2. Mirrors dev/R/parent_LD_diagnostic.R.
ld_decay_curve <- function(dose, Chr, Pos) {
  out <- vector("list", 0)
  for (ch in unique(Chr)) {
    j <- which(Chr == ch)
    if (length(j) < 2L) next
    o <- order(Pos[j]); j <- j[o]
    pos <- Pos[j]; G <- dose[j, , drop = FALSE]
    ## drop markers monomorphic in this panel
    p <- rowMeans(G, na.rm = TRUE) / 2
    poly <- is.finite(p) & pmin(p, 1 - p) >= MIN_MAF
    jj <- which(poly); if (length(jj) < 2L) next
    pos <- pos[jj]; G <- G[jj, , drop = FALSE]
    for (a in seq_len(length(pos) - 1L)) {
      b_hi <- min(a + MAXPART, length(pos)); b <- (a + 1L):b_hi
      d <- pos[b] - pos[a]; keep <- d <= MAXDIST
      if (!any(keep)) next
      b <- b[keep]; d <- d[keep]
      r <- suppressWarnings(cor(G[a, ], t(G[b, , drop = FALSE]),
                                use = "pairwise.complete.obs"))
      out[[length(out) + 1L]] <- data.table(dist = d, r2 = as.numeric(r)^2)
    }
  }
  if (!length(out)) return(NULL)
  res <- rbindlist(out)[is.finite(r2)]
  res[, dbin := cut(dist, DIST_BREAKS, labels = FALSE)]
  res[!is.na(dbin), .(r2 = mean(r2), n = .N,
                      dmid = (DIST_BREAKS[dbin] + DIST_BREAKS[dbin + 1L]) / 2), by = dbin]
}

## mean within-population expected heterozygosity (pi proxy)
pi_hat <- function(dose) {
  p <- rowMeans(dose, na.rm = TRUE) / 2
  mean(2 * p * (1 - p), na.rm = TRUE)
}

## Nei-style F_ST among populations (ratio of averages) at shared markers.
## freqs: pops x markers allele-frequency matrix
fst_hat <- function(freqs) {
  pbar <- colMeans(freqs, na.rm = TRUE)
  Hs   <- colMeans(2 * freqs * (1 - freqs), na.rm = TRUE)
  Ht   <- 2 * pbar * (1 - pbar)
  sum(Ht - Hs, na.rm = TRUE) / sum(Ht, na.rm = TRUE)
}

## ------------------- 1. EMPIRICAL targets (20 hybrid pops) ------------------
message("[1] empirical hybrids ...")
e <- new.env(); load(EMP_RDATA, envir = e)
GTh <- e$GTs_hybrids_005                       # inds x markers dosage
sdh <- as.data.table(e$sample_data)
maph<- as.data.table(e$map_hyb_005)

## the thinned 15k marker set the sim carries (from founder VCF IDs)
fmark <- unique(rbindlist(lapply(list.files(FVCF_DIR, "founders_ch.*\\.vcf$", full.names = TRUE),
  function(f) { d <- fread(f, skip = "#CHROM", select = 1:3, header = TRUE); d })))
sim_markers <- fmark$ID                         # "Chr<id>:Pos"
ci <- match(sim_markers, maph$marker)
ci <- ci[!is.na(ci)]
message(sprintf("    thinned markers matched in empirical map: %d / %d", length(ci), length(sim_markers)))
mk_chr <- maph$Chr[ci]; mk_pos <- maph$Pos[ci]

hyb_pops <- setdiff(unique(sdh$Population), c("aquilonia_parent", "polyctena_parent"))
if (TEST) hyb_pops <- hyb_pops[1:3]
emp_curves <- rbindlist(lapply(hyb_pops, function(pp) {
  rows <- which(sdh$Population == pp)
  dose <- t(GTh[rows, ci, drop = FALSE])        # markers x inds
  cur  <- ld_decay_curve(dose, mk_chr, mk_pos)
  if (is.null(cur)) return(NULL)
  cur[, pop := pp][]
}))
emp_ld <- emp_curves[, .(r2 = mean(r2), sd = sd(r2), dmid = dmid[1]), by = dbin][order(dbin)]
emp_pi <- mean(sapply(hyb_pops, function(pp) pi_hat(t(GTh[sdh$Population == pp, ci, drop = FALSE]))))
emp_freqs <- do.call(rbind, lapply(hyb_pops, function(pp)
  rowMeans(t(GTh[sdh$Population == pp, ci, drop = FALSE]), na.rm = TRUE) / 2))
emp_fst <- fst_hat(emp_freqs)
emp_ns  <- sapply(hyb_pops, function(pp) sum(sdh$Population == pp))  # per-pop sizes
message(sprintf("    empirical pi = %.4f ; among-pop F_ST = %.4f", emp_pi, emp_fst))
message(sprintf("    empirical per-pop sizes: %s (median %d)",
                paste(sort(emp_ns), collapse=","), as.integer(median(emp_ns))))

## ------------------- 2. SIM: per (setting, generation) ----------------------
message("[2] simulated demes ...")
sim_one_cell <- function(tag, gen) {
  vcfs <- list.files(SCREEN, sprintf("hyb_neutral_realfounders_%s_rep[0-9]+_gen%d\\.vcf$", tag, gen),
                     full.names = TRUE)
  if (!length(vcfs)) return(NULL)
  per <- lapply(seq_along(vcfs), function(k) {
    x <- read_vcf_dosage(vcfs[k])
    dose <- x$dose
    ## subsample this deme to an empirical pop size (matched-set), fixing r^2 bias
    if (MATCH_N) {
      tgt <- emp_ns[(k - 1L) %% length(emp_ns) + 1L]
      ni  <- ncol(dose)
      if (ni > tgt) { set.seed(MATCH_SEED + k); dose <- dose[, sample(ni, tgt), drop = FALSE] }
    }
    list(curve = ld_decay_curve(dose, x$Chr, x$Pos),
         pi    = pi_hat(dose),
         freq  = setNames(rowMeans(dose, na.rm = TRUE) / 2, x$marker),
         n     = ncol(dose))
  })
  curves <- rbindlist(lapply(seq_along(per), function(i)
    if (!is.null(per[[i]]$curve)) per[[i]]$curve[, deme := i][]))
  ld <- curves[, .(r2 = mean(r2), dmid = dmid[1]), by = dbin][order(dbin)]
  ## F_ST among the demes at their common markers
  common <- Reduce(intersect, lapply(per, function(z) names(z$freq)))
  fmat <- do.call(rbind, lapply(per, function(z) z$freq[common]))
  data.table(tag = tag, gen = gen, n_demes = length(vcfs),
             mean_n = mean(sapply(per, `[[`, "n")),
             pi = mean(sapply(per, `[[`, "pi")), fst = fst_hat(fmat),
             ld = list(ld))
}

cells <- CJ(tag = SETTINGS$tag, gen = GENS, sorted = FALSE)
sim <- rbindlist(mcMap(sim_one_cell, cells$tag, cells$gen, mc.cores = CORES), fill = TRUE)
sim <- SETTINGS[sim, on = "tag"]

## ------------------- 3. distance-to-empirical (LD-decay RMSE) ---------------
emp_vec <- setNames(emp_ld$r2, emp_ld$dbin)
sim[, ld_rmse := sapply(ld, function(x) {
  m <- match(x$dbin, as.integer(names(emp_vec)))
  sqrt(mean((x$r2 - emp_vec[m])^2, na.rm = TRUE))
})]
sim[, `:=`(pi_diff = pi - emp_pi, fst_diff = fst - emp_fst)]

## ------------------- 4. outputs --------------------------------------------
saveRDS(list(sim = sim, emp_ld = emp_ld, emp_pi = emp_pi, emp_fst = emp_fst,
             params = list(MAXDIST = MAXDIST, DIST_BREAKS = DIST_BREAKS, n_markers = length(ci))),
        OUTRDS)

## LD-decay: facet by founder setting, line per generation, empirical dashed
ld_long <- sim[, rbindlist(lapply(seq_len(.N), function(i)
  ld[[i]][, .(setting = setting[i], gen = gen[i], dmid, r2)]))]
ld_long[, setting := factor(setting, levels = SETTINGS$setting)]
p_ld <- ggplot(ld_long, aes(dmid/1000, r2, colour = factor(gen), group = gen)) +
  geom_line(linewidth = 0.6) +
  geom_line(data = emp_ld, aes(dmid/1000, r2), inherit.aes = FALSE,
            colour = "black", linetype = 2, linewidth = 0.9) +
  facet_wrap(~ setting) + scale_x_log10() +
  scale_colour_viridis_d(name = "generation", option = "plasma", end = 0.9) +
  labs(x = "distance (kb, log)", y = expression("composite " * r^2),
       title = "Module E screen: LD-decay vs empirical (dashed black)",
       subtitle = sprintf("thinned DI>-25 markers (%d); founder setting per facet", length(ci))) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave(file.path(FIGDIR, "moduleE_screen_lddecay.pdf"), p_ld, width = 11, height = 7)

## pi and F_ST vs generation, per founder setting, empirical reference
mdt <- melt(sim[, .(setting = factor(setting, SETTINGS$setting), gen, pi, fst)],
            id.vars = c("setting","gen"))
ref <- data.table(variable = c("pi","fst"), y = c(emp_pi, emp_fst))
p_pf <- ggplot(mdt, aes(gen, value, colour = setting)) +
  geom_line() + geom_point(size = 1) +
  geom_hline(data = ref, aes(yintercept = y), linetype = 2) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(x = "generation", y = NULL,
       title = "Module E screen: pi and among-pop F_ST vs empirical (dashed)") +
  theme_bw(base_size = 11)
ggsave(file.path(FIGDIR, "moduleE_screen_pi_fst.pdf"), p_pf, width = 10, height = 5)

## ------------------- 5. best-matching cell(s) ------------------------------
best <- sim[order(ld_rmse), .(setting, gen, n_demes, ld_rmse = round(ld_rmse,4),
                              pi = round(pi,4), pi_diff = round(pi_diff,4),
                              fst = round(fst,4), fst_diff = round(fst_diff,4))]
cat("\n=== empirical targets ===\n")
cat(sprintf("pi = %.4f ; among-pop F_ST = %.4f\n", emp_pi, emp_fst))
cat("\n=== cells ranked by LD-decay RMSE to empirical (top 12) ===\n")
print(head(best, 12))
cat("\nsaved: ", OUTRDS, "\n  ", file.path(FIGDIR, "moduleE_screen_lddecay.pdf"),
    "\n  ", file.path(FIGDIR, "moduleE_screen_pi_fst.pdf"), "\n")
