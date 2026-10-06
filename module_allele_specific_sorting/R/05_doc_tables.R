## =========================================================================
## module_allele_specific_sorting -- 05: LaTeX tables + figure copies for
## doc_manuscript/ (numerical authority = the saved outputs of 02-04).
## Re-run after any of 02-04 is re-run. 04 is optional: if its output is
## absent, a clearly marked placeholder table is written instead.
##
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/05_doc_tables.R
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
DOC <- "module_allele_specific_sorting/doc_manuscript"
dir.create(file.path(DOC, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DOC, "figures"), showWarnings = FALSE, recursive = TRUE)

ci <- function(m, lo, hi, d = 3) sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), m, lo, hi)
tex_bin <- function(x) gsub("-", "--", gsub(">", "$>$", gsub("<", "$<$", as.character(x))))
write_tex <- function(lines, file) { writeLines(lines, file.path(DOC, "tables", file)); cat("[05] wrote", file, "\n") }

## ---- Table 1: distance decay (02) ------------------------------------------------
d <- readRDS(file.path(OUT_DATA, "02_decay.rds"))$bp
cols <- c(ceiling = "Ceiling $\\overline{G_iG_j}$", among_LD = "Among-pop.\\ LD $\\overline{r_{ST}^2}$",
          realised = "Realised fraction", conc_resid = "Concordance (resid.)",
          within_LD_adj = "Within-pop.\\ LD $\\overline{r_w^2}$")
w <- dcast(d[qty %in% names(cols), .(bin, qty, v = ci(mean, lo, hi, 4))], bin ~ qty, value.var = "v")
setcolorder(w, c("bin", names(cols)))
t1 <- c("\\begin{table}[htbp]", "\\centering", "\\scriptsize",
        "\\caption{Among- and within-population LD of DI25 unit pairs by physical distance. Means per pair with 95\\% chromosome-block bootstrap intervals (2,000 replicates). Ceiling: the among-population LD two loci with their observed $F_{ST}$ would show if they partitioned the populations identically; realised fraction $=\\sum r_{ST}^2/\\sum G_iG_j$, whose chance level for two unrelated 20-population profiles is $1/19=0.053$. Unlinked: all cross-chromosome pairs.}",
        "\\label{tab:decay}",
        paste0("\\resizebox{\\textwidth}{!}{\\begin{tabular}{l", strrep("c", length(cols)), "}"), "\\toprule",
        paste0("Distance & ", paste(cols, collapse = " & "), " \\\\"), "\\midrule",
        w[, paste0(tex_bin(bin), " & ", do.call(paste, c(.SD, sep = " & ")), " \\\\"), .SDcols = names(cols)],
        "\\bottomrule", "\\end{tabular}}", "\\end{table}")
write_tex(t1, "decay_summary.tex")

## ---- Table 2: sorted anchors vs matched controls (03) -----------------------------
a3 <- readRDS(file.path(OUT_DATA, "03_anchor_decay.rds"))
r <- a3$result[dir == "all sorted" & stat %in% c("conc", "conc_resid", "r_w", "r_w_poly")]
slab <- c(conc = "Concordance", conc_resid = "Concordance (resid.)",
          r_w = "Within-pop.\\ LD, pooled", r_w_poly = "Within-pop.\\ LD, both segregating")
r[, stat := factor(stat, levels = names(slab))]; setorder(r, stat, bin)
t2 <- c("\\begin{table}[htbp]", "\\centering", "\\small",
        sprintf("\\caption{Neighbourhood concordance and within-population LD around %d directionally sorted units (anchors) and $K=%d$ unsorted controls per anchor matched on $F_{ST}$, local recombination rate, parental differentiation and local unit density. Difference: anchor $-$ control with 95\\%% chromosome-block bootstrap interval. Pooled within-population LD is diluted for sorted units, which are monomorphic in most populations; the both-segregating version restricts each pair to populations where both loci segregate.}",
                nrow(a3$matches) / a3$K_CTRL, a3$K_CTRL),
        "\\label{tab:anchors}",
        "\\begin{tabular}{llccc}", "\\toprule",
        "Statistic & Distance & Anchors & Controls & Difference [95\\% CI] \\\\", "\\midrule",
        r[, sprintf("%s & %s & %.3f & %.3f & %s \\\\", ifelse(duplicated(stat), "", slab[as.character(stat)]),
                    tex_bin(bin), anchor, control, ci(diff, lo, hi))],
        "\\bottomrule", "\\end{tabular}", "\\end{table}")
write_tex(t2, "anchor_summary.tex")

b <- a3$balance
t2b <- c("\\begin{table}[htbp]", "\\centering", "\\small",
         "\\caption{Covariate balance of the anchor--control matching (means).}", "\\label{tab:balance}",
         "\\begin{tabular}{lcccc}", "\\toprule",
         "Set & $F_{ST}$ & $\\log_{10}$ cM/Mb & Parental differentiation & Units within $\\pm$100 kb \\\\", "\\midrule",
         b[, sprintf("%s & %.3f & %.2f & %.3f & %.1f \\\\", set, FST, log_recomb, parent_diff, density100k)],
         "\\bottomrule", "\\end{tabular}", "\\end{table}")
write_tex(t2b, "anchor_balance.tex")

## ---- Table 3: long-range tails (04), or placeholder --------------------------------
f4 <- file.path(OUT_DATA, "04_longrange_tail.rds")
if (file.exists(f4)) {
  r4 <- readRDS(f4)
  slab4 <- c(conc_resid = "Among-pop.\\ concordance, residualised", conc_raw = "Among-pop.\\ concordance, raw",
             r_w_adj = "Within-pop.\\ LD, hybrid-index adj.")
  rows <- unlist(lapply(names(slab4), function(st) {
    x <- r4[[st]]; tt <- x$tails[abs(thr - x$t_star) < 1e-9]
    an <- x$annotation
    ## tail asymmetry (pos - neg)/(pos + neg): residual population structure / kinship inflate
    ## both tails; same-ancestry epistasis predicts a POSITIVE asymmetry beyond the null's
    asn <- x$tails_null[abs(thr - x$t_star) < 1e-9, (pos - neg) / (pos + neg)]
    sprintf("%s & %.2f & %s (%.1f$\\times$) & %s (%.1f$\\times$) & %.3f & [%.3f, %.3f] & %.3f [%.3f, %.3f] \\\\",
            slab4[st], x$t_star,
            format(tt$pos, big.mark = ","), tt$pos / tt$pos_null_mean,
            format(tt$neg, big.mark = ","), tt$neg / tt$neg_null_mean,
            (tt$pos - tt$neg) / (tt$pos + tt$neg), quantile(asn, 0.025), quantile(asn, 0.975),
            an[annotation == "bdmi"]$frac_tail_endpoints, an[annotation == "bdmi"]$null_lo, an[annotation == "bdmi"]$null_hi)
  }))
  t3 <- c("\\begin{table}[htbp]", "\\centering", "\\scriptsize",
          sprintf("\\caption{Cross-chromosome association tails versus a chromosome-wise permutation null (%d permutations). Positive tail: pairs with $r\\geq T^*$ (same parental ancestry at both loci); negative: $r\\leq -T^*$; in parentheses, the ratio to the null mean. Asymmetry: $(n_+-n_-)/(n_++n_-)$ with the null 95\\%% interval. Both tails are inflated by residual population structure (among populations) and kinship (within populations), which the permutation removes; epistasis between compatible ancestries predicts a positive asymmetry beyond the null's. BDMI: fraction of positive-tail endpoints inside BDMI candidate regions (cutoff %d), with the null 95\\%% interval.}",
                  r4$meta$B_NULL, r4$meta$bdmi_cutoff),
          "\\label{tab:longrange}",
          "\\resizebox{\\textwidth}{!}{\\begin{tabular}{lcccccc}", "\\toprule",
          "Statistic & $T^*$ & $+$ tail & $-$ tail & Asymmetry & Null asymmetry & BDMI endpoints [null] \\\\", "\\midrule",
          rows, "\\bottomrule", "\\end{tabular}}", "\\end{table}")
} else {
  t3 <- c("\\begin{table}[htbp]", "\\centering",
          "\\caption{Cross-chromosome association tails. \\textbf{PENDING}: \\texttt{R/04\\_longrange\\_tail.R} has not finished; re-run \\texttt{R/05\\_doc\\_tables.R} when it has.}",
          "\\label{tab:longrange}", "\\fbox{\\parbox{0.8\\textwidth}{\\centering Results pending.}}", "\\end{table}")
}
write_tex(t3, "longrange_summary.tex")

## ---- figures ---------------------------------------------------------------------------
figs <- c("02_fst_vs_ld_decay.pdf", "02_fst_vs_ld_decay_by_fst_class.png", "02_decay_cM.png",
          "03_sorted_anchor_decay.pdf", "04_longrange_tail.pdf")
for (f in figs) {
  src <- file.path(OUT_FIG, f)
  if (file.exists(src)) { file.copy(src, file.path(DOC, "figures", f), overwrite = TRUE); cat("[05] copied", f, "\n") }
  else cat("[05] MISSING (not yet produced):", f, "\n")
}
