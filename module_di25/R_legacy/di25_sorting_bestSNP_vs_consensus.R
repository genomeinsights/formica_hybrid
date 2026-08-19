## =========================================================
## module_di25 -- sorting: filled best-SNP vs eMLG consensus (matched units)
## =========================================================
## Compares ancestry sorting (Module A method, di25_sorting.R conventions) between
## two representations of the SAME 4,010 eMLG blocks (n_loci > 2, where the
## consensus actually differs from any single SNP):
##   (A) CONSENSUS      -- consensus_dosage() over the block members (as di25_sorting.R)
##   (B) BEST-SNP FILLED -- the block's best-consensus-correlated member SNP, with its
##       missing calls filled from the (oriented, rounded) consensus.
## Everything else (phi = 0.85, tau sweep, sort_rule = "binom", ungated DI/pmaf,
## 194 individuals) is identical to di25_sorting.R, so any difference is purely the
## unit representation.
##
## Run from the repo root:  Rscript module_di25/R/di25_sorting_bestSNP_vs_consensus.R
## =========================================================

suppressMessages(library(data.table))
devtools::load_all("~/gitlab/LDscnR/")             # consensus_dosage(), ohta_fast_prepare()
source("moduleA_sorting/R/parallelism_stats.R")    # parallelism_stats(), classify_sort()

CM_STAMP <- "cM5"
OUTDIR   <- "module_di25/data"
FIGDIR   <- "module_di25/Figures"
TAU_GRID <- c(0.5, 0.6, 0.7, 0.8)
FIX_TH   <- 0.15
SORT_RULE <- "binom"; ALPHA <- 0.05

## ---- inputs (identical setup to di25_sorting.R) ---------------------------
inp <- readRDS(file.path(OUTDIR, "di25_inputs.rds"))
res <- readRDS(sprintf("module_di25/data/di25_clustering_%s.rds", CM_STAMP)); g <- res$groups
bm  <- as.data.table(readRDS(file.path(OUTDIR, "di25_emlg_best_member.rds")))

e2 <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2)
sd <- e2$sample_data_with_parents
DI_vec <- setNames(e2$map_hyb_005$DiagnosticIndex, e2$map_hyb_005$marker)

GTs_all <- rbind(inp$GTs_hyb, inp$GTs_par)
keep <- rownames(GTs_all) %in% sd$Sample_ID
GTs_all <- GTs_all[keep, ]
pops <- sd$Population[match(rownames(GTs_all), sd$Sample_ID)]
aqu_pops <- "aquilonia_parent"; pol_pops <- "polyctena_parent"
hybrid_pops <- setdiff(unique(pops), c(aqu_pops, pol_pops))
parent_rows <- grepl("_parent$", pops)

## ---- restrict to the 4,010 eMLG blocks (n_loci > 2) -----------------------
gd  <- as.data.table(g)[n_loci > 2]
gd[, best_marker := bm$best_marker[match(group_id, bm$group_id)]]
stopifnot(!anyNA(gd$best_marker))
cat(sprintf("Comparing %d eMLG blocks (n_loci > 2) on %d individuals\n", nrow(gd), nrow(GTs_all)))

## ---- build the two matched unit matrices (194 individuals x blocks) -------
E_cons <- vapply(seq_len(nrow(gd)), function(i)
  consensus_dosage(GTs_all, gd$members[[i]]), numeric(nrow(GTs_all)))
colnames(E_cons) <- gd$group_id

E_best <- vapply(seq_len(nrow(gd)), function(i) {
  snp  <- GTs_all[, gd$best_marker[i]]
  cons <- E_cons[, i]
  r <- suppressWarnings(cor(snp, cons, use = "pairwise.complete.obs"))
  if (is.finite(r) && r < 0) cons <- 2 - cons          # orient consensus to the SNP
  miss <- is.na(snp)
  snp[miss] <- pmin(pmax(round(cons[miss]), 0), 2)
  snp
}, numeric(nrow(GTs_all)))
colnames(E_best) <- gd$group_id

## ---- per-unit DI (max over members) and pmaf ------------------------------
DI_u <- setNames(vapply(gd$members, function(mk) {
  v <- DI_vec[mk]; if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
}, numeric(1)), gd$group_id)
pmaf_of <- function(M) { f <- colMeans(M[parent_rows, , drop = FALSE], na.rm = TRUE) / 2
                         setNames(pmin(f, 1 - f), colnames(M)) }

## ---- parallelism_stats for each representation ----------------------------
run_ps <- function(M) {
  prep <- ohta_fast_prepare(M, pops = pops)
  parallelism_stats(prep, hybrid_pops = hybrid_pops, aqu_pops = aqu_pops,
                    pol_pops = pol_pops, fix_th = FIX_TH, DI = DI_u, min_DI = NULL,
                    parent_maf = pmaf_of(M), min_parent_maf = NULL,
                    sort_rule = SORT_RULE, alpha = ALPHA)
}
ps_cons <- run_ps(E_cons); ps_best <- run_ps(E_best)

## ---- aggregate tau sweep for each level -----------------------------------
tally_level <- function(ps, level) {
  base <- ps[differentiated == TRUE & n_obs > 0]
  rbindlist(lapply(TAU_GRID, function(tau) {
    cls <- classify_sort(base$n_aqu, base$n_pol, base$n_obs, sort_th = tau,
                         sort_rule = SORT_RULE, alpha = ALPHA)
    n_aqu <- sum(cls == "aquilonia"); n_pol <- sum(cls == "polyctena")
    n_sorted <- sum(cls %in% c("aquilonia","polyctena","unresolved","ambiguous"))
    data.table(level = level, tau = tau, n_units = nrow(base), n_sorted = n_sorted,
               pct_sorted = 100*n_sorted/nrow(base), toward_aqu = n_aqu, toward_pol = n_pol,
               pct_aqu_of_resolved = 100*n_aqu/(n_aqu+n_pol))
  }))
}
sweep <- rbind(tally_level(ps_cons, "consensus"), tally_level(ps_best, "best_SNP"))

cat("\n===== aggregate sorting sweep: consensus vs best-SNP-filled =====\n")
print(sweep[, .(level, tau, n_sorted, pct_sorted = round(pct_sorted,1),
                toward_aqu, toward_pol, pct_aqu_of_resolved = round(pct_aqu_of_resolved,1))])

## ---- per-unit concordance at each tau -------------------------------------
class_vec <- function(ps, tau) {
  cls <- rep("unsorted", nrow(ps))
  ok <- ps$differentiated == TRUE & ps$n_obs > 0
  cls[ok] <- {
    c0 <- classify_sort(ps$n_aqu[ok], ps$n_pol[ok], ps$n_obs[ok], sort_th = tau,
                        sort_rule = SORT_RULE, alpha = ALPHA)
    ifelse(c0 %in% c("aquilonia","polyctena","unresolved","ambiguous"), c0, "unsorted")
  }
  setNames(cls, ps$marker)
}
setkey(ps_cons, marker); setkey(ps_best, marker)
cat("\n===== per-unit concordance (same 4,010 blocks) =====\n")
conc <- rbindlist(lapply(TAU_GRID, function(tau) {
  cc <- class_vec(ps_cons, tau); cb <- class_vec(ps_best, tau)[names(cc)]
  sorted_c <- cc != "unsorted"; sorted_b <- cb != "unsorted"
  data.table(
    tau = tau,
    agree_class      = round(100*mean(cc == cb), 1),
    both_sorted      = sum(sorted_c & sorted_b),
    only_consensus   = sum(sorted_c & !sorted_b),
    only_bestSNP     = sum(!sorted_c & sorted_b),
    dir_disagree     = sum(sorted_c & sorted_b & cc != cb)
  )
}))
print(conc)

saveRDS(list(sweep = sweep, concordance = conc, ps_cons = ps_cons, ps_best = ps_best),
        file.path(OUTDIR, "di25_sorting_bestSNP_vs_consensus.rds"))

## ---- figure ---------------------------------------------------------------
suppressMessages(library(ggplot2)); suppressMessages(library(patchwork))
p1 <- ggplot(sweep, aes(tau, pct_sorted, colour = level)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_colour_manual(values = c(consensus = "#440154", best_SNP = "#21918C")) +
  labs(x = expression(tau), y = "% units sorted", colour = NULL,
       title = "Fraction sorted") + theme_bw(base_size = 11)
p2 <- ggplot(sweep, aes(tau, pct_aqu_of_resolved, colour = level)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) + geom_hline(yintercept = 50, linetype = 3) +
  scale_colour_manual(values = c(consensus = "#440154", best_SNP = "#21918C")) +
  labs(x = expression(tau), y = "% aquilonia of resolved", colour = NULL,
       title = "Direction") + theme_bw(base_size = 11)
ggsave(file.path(FIGDIR, "di25_sorting_bestSNP_vs_consensus.png"),
       (p1 | p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
       width = 9, height = 4.5, dpi = 200)
cat(sprintf("\nWrote %s and the figure.\n", file.path(OUTDIR, "di25_sorting_bestSNP_vs_consensus.rds")))
