## =========================================================
## MODULE D (bidirectional arm) -- the symmetric-DMI test on SNP-level bidirectional loci.
## =========================================================
## An EQUAL-fitness Dobzhansky-Muller incompatibility drives strong sorting in RANDOM
## directions across replicate populations (each deme independently fixes one of two
## equally-fit compatible combinations) -- i.e. BIDIRECTIONAL sorting, not the
## directional sorting Module A's binomial detects (whose single-locus test is null
## for these by construction: n_aqu ~ n_pol). The signature is at the PAIR level: one
## symmetric DMI creates TWO bidirectional loci whose per-population directions are
## strongly CORRELATED across replicates (cis-compatible -> same direction, r>0;
## trans-compatible -> opposite, r<0). This script screens exactly that.
##
## WHY SNP-LEVEL. Bidirectional sorting is real at the SNP level (911 SNPs) but is
## averaged BELOW the sort threshold at the eMLG consensus (0 clusters qualify) --
## the consensus dilutes prop_fixed. So we take, per LD cluster, its MOST bidirectional
## member SNP as the unit (599 reps), de-redundifying without losing the signal.
##
## WHY THE ANCESTRY-AXIS CONFOUND IS ~ABSENT HERE. The confound that swamped the general
## R_st screen is the shared ancestry gradient -- but bidirectional loci by definition do
## NOT track it (they sort both ways), so among-population covariance between two
## bidirectional loci is a comparatively clean DMI signal. A per-locus PERMUTATION null
## (shuffle each locus's per-population values independently -> destroys co-sorting,
## keeps marginals) is therefore a reasonable PRE-SIM calibration here; the neutral sim
## still gives the final null and licenses signal aligned with residual drift structure.
##
## Unlinked = different chromosome or same chromosome > LINK_CM (as the other arms).
## Cross-chromosome paralogy filtered (shared dev/R/moduleD_paralogy.R). Descriptive.
## Run from repo root. Reads moduleA_snp.rds, clustering, hybrids_only, recmap,
## moduleD_ohta.rds (baseline). Writes data/moduleD_bidirectional.rds,
## Figures/moduleD_bidirectional.{pdf,png}.

suppressMessages({ library(data.table); library(ggplot2); library(patchwork); library(parallel) })
source("dev/R/moduleD_paralogy.R")
set.seed(1)

## ---- PARAMETERS ---------------------------------------------------------
LINK_CM <- 10; RST_THR <- 0.7; PARALOGY_R <- 0.9; NPERM <- 1000; CORES <- 8; NET_Q <- 0.99
RECMAP  <- "data/Frufa_DTOL_PR.ref_genome.recmap"

## ---- inputs -------------------------------------------------------------
snp   <- readRDS("data/moduleA_snp.rds")
clust <- readRDS("data/eMLG_5loci_0025_cM05.rds"); groups <- clust$groups
e1 <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e1)
GT <- e1$GTs_hybrids_005                                   # raw ancestry-polarised dosages
pops <- as.data.table(e1$sample_data)[match(rownames(GT), Sample_ID), Population]
pop_lv <- sort(unique(pops)); pop_idx <- split(seq_along(pops), pops)

## ---- bidirectional reps: one (most bidirectional) SNP per LD cluster -----
bi <- snp[sort_class == "bidirectional"]
bi[, c("Chr", "Pos") := tstrsplit(marker, ":", fixed = TRUE)][, Pos := as.integer(Pos)]
m2g <- groups[, .(marker = unlist(members)), by = group_id]
bi <- m2g[bi, on = "marker"]
reps <- bi[order(-bi_score, -prop_fixed)][, .SD[1], by = group_id][
  , .(group_id, marker, Chr, Pos, DI, bi_score, prop_fixed, f_aqu_pooled)]
## cM position of each rep SNP (interpolate from the recmap)
rec <- fread(RECMAP); setnames(rec, 1:4, c("chr", "pos", "cM", "cMMb")); rec[, Chr := sub("chromosome_", "Chr", chr)]
reps[, cM := NA_real_]
for (ch in reps[, unique(Chr)]) { r <- rec[Chr == ch]; if (nrow(r) < 2) next
  reps[Chr == ch, cM := approx(r$pos, r$cM, xout = Pos, rule = 2)$y] }
K <- nrow(reps)
message(sprintf("[bi] %d bidirectional reps (1 per cluster) over %d chromosomes; Chr24 = %d",
                K, reps[, uniqueN(Chr)], reps[Chr == "Chr24", .N]))

## ---- per-population aquilonia-oriented frequency matrix (pops x reps) ----
## orient each rep's dosage to aquilonia using its pooled aqu frequency (flip if the
## raw polarity is the complement), then take per-population means.
GTr <- GT[, reps$marker, drop = FALSE]
for (j in seq_len(K)) {
  d <- GTr[, j]; o <- mean(d, na.rm = TRUE) / 2; fa <- reps$f_aqu_pooled[j]
  if (abs(o - (1 - fa)) < abs(o - fa)) GTr[, j] <- 2 - d       # flip to aqu orientation
}
Fbi <- vapply(pop_lv, function(p) colMeans(GTr[pop_idx[[p]], , drop = FALSE], na.rm = TRUE),
              numeric(K))                                       # K x nPops
Fbi <- t(Fbi)                                                   # nPops x K
colnames(Fbi) <- reps$marker

## ---- unlinked pairwise co-sorting R_st (exhaustive over the reps) --------
Rst <- cor(Fbi, use = "pairwise.complete.obs")
ut  <- which(upper.tri(Rst), arr.ind = TRUE)                    # i<j
pairs <- data.table(i = ut[, 1], j = ut[, 2], R_st = Rst[ut])
pairs[, `:=`(rep1 = reps$marker[i], rep2 = reps$marker[j],
             Chr1 = reps$Chr[i], Chr2 = reps$Chr[j], cM1 = reps$cM[i], cM2 = reps$cM[j])]
pairs[, unlinked := Chr1 != Chr2 | (Chr1 == Chr2 & is.finite(cM1) & is.finite(cM2) & abs(cM1 - cM2) > LINK_CM)]
unl <- pairs[unlinked == TRUE & is.finite(R_st)]
message(sprintf("[bi] %d unlinked bidirectional pairs; |R_st|>=%.2f: %d (%.3f%%)",
                nrow(unl), RST_THR, unl[abs(R_st) >= RST_THR, .N], 100 * unl[, mean(abs(R_st) >= RST_THR)]))

## ---- paralogy filter (cross-chromosome duplicates) ----------------------
unl <- flag_paralogy(unl, "rep1", "rep2", GTr, pops, thr = PARALOGY_R, cores = CORES)
message(sprintf("[bi] paralogy: %d/%d unlinked pairs flagged as duplicates", sum(unl$paralog), nrow(unl)))

## ---- permutation null: shuffle each locus's per-pop values independently -
## => loci keep their marginal sorting but co-sorting is destroyed. Recompute the
## unlinked |R_st|>=thr count NPERM times to calibrate the observed excess.
chr_r <- reps$Chr; cm_r <- reps$cM
unlinked_mask <- outer(seq_len(K), seq_len(K), function(a, b)
  chr_r[a] != chr_r[b] | (chr_r[a] == chr_r[b] & is.finite(cm_r[a]) & is.finite(cm_r[b]) & abs(cm_r[a] - cm_r[b]) > LINK_CM))
um <- which(upper.tri(matrix(0, K, K)) & unlinked_mask)         # linear idx of unlinked upper-tri
obs_n <- sum(abs(Rst[um]) >= RST_THR, na.rm = TRUE)
perm_n <- unlist(mclapply(seq_len(NPERM), function(b) {
  Fp <- apply(Fbi, 2, function(col) col[sample.int(length(col))])   # shuffle rows within each column
  Rp <- cor(Fp, use = "pairwise.complete.obs")
  sum(abs(Rp[um]) >= RST_THR, na.rm = TRUE)
}, mc.cores = CORES))
emp_p <- (1 + sum(perm_n >= obs_n)) / (NPERM + 1)
message(sprintf("[bi] co-sorting excess: observed %d unlinked |R_st|>=%.2f vs permutation null %.1f (95%%<=%d); empirical p=%.4g",
                obs_n, RST_THR, mean(perm_n), quantile(perm_n, 0.95), emp_p))

## general differentiated-pair baseline (ancestry-axis-dominated) for context
base_frac <- tryCatch(readRDS("data/moduleD_ohta.rds")$tail_contrast[R_st_cut == 0.7, frac_unlinked], error = function(e) NA_real_)

## ---- candidate network (top-quantile |R_st|, non-paralog, unlinked) ------
cand <- unl[paralog == FALSE][order(-abs(R_st))]
thr_net <- cand[, quantile(abs(R_st), NET_Q, na.rm = TRUE)]
net_edges <- cand[abs(R_st) >= thr_net, .(rep1, rep2, R_st, Chr1, Chr2,
                                          coupling = ifelse(R_st > 0, "cis", "trans"))]
network <- list(edges = net_edges, threshold = thr_net)
if (requireNamespace("igraph", quietly = TRUE) && nrow(net_edges) > 0) {
  g <- igraph::graph_from_data_frame(net_edges[, .(rep1, rep2)], directed = FALSE)
  comps <- igraph::components(g)
  d <- igraph::degree(g)
  network$hub_degree <- data.table(marker = names(d), degree = as.integer(d),
                                    Chr = reps$Chr[match(names(d), reps$marker)])[order(-degree)]
  message(sprintf("[bi] candidate network: %d loci, %d components (largest %d); Chr24 hubs = %d",
                  igraph::gorder(g), comps$no, max(comps$csize), network$hub_degree[Chr == "Chr24", .N]))
}

## ---- outputs + figure ---------------------------------------------------
saveRDS(list(reps = reps, unlinked = unl, network = network,
             perm = list(obs = obs_n, null = perm_n, emp_p = emp_p, thr = RST_THR),
             baseline_frac = base_frac,
             params = list(LINK_CM = LINK_CM, RST_THR = RST_THR, PARALOGY_R = PARALOGY_R, NPERM = NPERM)),
        "data/moduleD_bidirectional.rds")
message("[bi] saved data/moduleD_bidirectional.rds")

dir.create("Figures", showWarnings = FALSE)
th <- theme_classic(base_size = 8) +
  theme(plot.title = element_text(size = 8.5), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6.5), legend.position = "bottom",
        legend.title = element_blank(), legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"), plot.margin = margin(4, 9, 2, 4))

## (a) observed unlinked-pair R_st vs one permutation-null realisation
one_perm <- { Fp <- apply(Fbi, 2, function(col) col[sample.int(length(col))]); Rp <- cor(Fp, use = "pairwise.complete.obs"); Rp[um] }
dd <- rbind(data.table(R_st = Rst[um], set = "observed"),
            data.table(R_st = one_perm, set = "permuted null"))
p_a <- ggplot(dd, aes(R_st, colour = set)) + geom_vline(xintercept = c(-RST_THR, RST_THR), linetype = 2, colour = "grey70") +
  geom_density() + scale_colour_manual(values = c(observed = "#0072B2", `permuted null` = "grey55")) +
  scale_y_sqrt() + labs(x = expression(R[st]~"(bidirectional co-sorting)"), y = "density (sqrt)",
                        title = "a  Co-sorting vs permutation null") + th

## (b) permutation null of the strong-pair count, with the observed value
p_b <- ggplot(data.table(n = perm_n), aes(n)) + geom_histogram(bins = 40, fill = "grey70") +
  geom_vline(xintercept = obs_n, colour = "#D55E00", linewidth = 0.8) +
  labs(x = sprintf("# unlinked |R_st|>=%.1f pairs", RST_THR), y = "permutations",
       title = sprintf("b  Excess co-sorting (p=%.3g)", emp_p)) + th

## (c) hub degree of the candidate network (Chr24 highlighted)
p_c <- if (!is.null(network$hub_degree) && nrow(network$hub_degree) > 0) {
  hd <- head(network$hub_degree, 20)[, lab := factor(marker, levels = rev(marker))]
  ggplot(hd, aes(degree, lab, fill = Chr == "Chr24")) + geom_col(width = 0.7) +
    scale_fill_manual(values = c(`FALSE` = "grey55", `TRUE` = "#D55E00"), labels = c("other", "Chr24")) +
    labs(x = "co-sorting partners", y = NULL, title = "c  Candidate hubs") + th
} else patchwork::plot_spacer()

fig <- p_a + p_b + p_c + plot_layout(widths = c(1.1, 1, 1))
ggsave("Figures/moduleD_bidirectional.pdf", fig, width = 210, height = 74, units = "mm")
ggsave("Figures/moduleD_bidirectional.png", fig, width = 210, height = 74, units = "mm", dpi = 300)
cat(sprintf("\n[done] observed %d vs null mean %.1f (p=%.3g) | %d candidate edges | baseline frac %.4g\n",
            obs_n, mean(perm_n), emp_p, nrow(net_edges), base_frac))
