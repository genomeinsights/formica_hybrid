## =========================================================================
## module_allele_specific_sorting -- 05: LaTeX tables + figure copies for
## doc_manuscript/ (numerical authority = the saved outputs of 02-04).
## Re-run after any of 02-04 is re-run; never edit the generated tables by hand.
##
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/05_doc_tables.R
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
DOC <- "module_allele_specific_sorting/doc_manuscript"
dir.create(file.path(DOC, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DOC, "figures"), showWarnings = FALSE, recursive = TRUE)
unlink(list.files(file.path(DOC, "tables"), full.names = TRUE))
unlink(list.files(file.path(DOC, "figures"), full.names = TRUE))

ci <- function(m, lo, hi, d = 3) sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), m, lo, hi)
tex_bin <- function(x) gsub("-", "--", gsub(">", "$>$", gsub("<", "$<$", as.character(x))))
write_tex <- function(lines, file) { writeLines(lines, file.path(DOC, "tables", file)); cat("[05] wrote", file, "\n") }

nb <- readRDS(file.path(OUT_DATA, "04_neighbourhood.rds"))

## ---- balance ----------------------------------------------------------------------------
vlab <- c(FST = "$F_{ST}$", log_recomb = "$\\log_{10}$ recombination rate (cM/Mb)",
          parent_diff = "Parental allele-frequency difference", density100k = "Units within $\\pm$100 kb",
          log2_size = "$\\log_2$ SNPs in own LD unit", span_kb = "Span of own LD unit (kb)")
b <- nb$balance
t_bal <- c("\\begin{table}[htbp]", "\\centering", "\\small",
           sprintf("\\caption{Matching of %d of %d directionally sorted loci to up to %d unsorted loci each (caliper %g SD on every matched covariate). SMD: standardised mean difference, sorted minus unsorted, against all unsorted loci (before) and against the matched controls (after). The last two rows describe the loci's own LD units and were not matched on.}",
                   nb$n_matched, nb$n_sorted, nb$K_CTRL, nb$CALIPER),
           "\\label{tab:balance}",
           "\\resizebox{\\textwidth}{!}{\\begin{tabular}{lrrrrr}", "\\toprule",
           "Variable & Sorted & All unsorted & Matched controls & SMD before & SMD after \\\\", "\\midrule",
           b[, sprintf("%s & %.3f & %.3f & %.3f & %.2f & %.2f \\\\", vlab[variable], mean_anchor, mean_all_unsorted,
                       mean_controls, smd_before, smd_after)],
           "\\bottomrule", "\\end{tabular}}", "\\end{table}")
t_bal <- append(t_bal, "\\midrule", after = which(startsWith(t_bal, "Units within")))
write_tex(t_bal, "balance.tex")

## ---- matched-set differences ----------------------------------------------------------
r <- copy(nb$result)
slab <- c(similarity = "Ancestry-profile similarity", r_w_poly = "Within-population LD")
setorder(r, stat, marker_set, bin)
t_nb <- c("\\begin{table}[htbp]", "\\centering", "\\small",
          "\\caption{Neighbourhoods of sorted loci versus matched unsorted loci. Values are means over all neighbouring markers in each distance bin; difference: mean over sorted loci of (sorted locus minus the mean of its own matched controls), with 95\\% intervals from a bootstrap that resamples the sorted loci's chromosomes and keeps each matched set intact. Within-population LD uses only populations in which both markers segregate. LD-reduced: neighbours are other LD units; unpruned: all DI25 SNPs within 2~Mb, including those in the sorted locus's own LD unit. $n$: sorted loci with at least one neighbour in the bin.}",
          "\\label{tab:neighbourhood}",
          "\\resizebox{\\textwidth}{!}{\\begin{tabular}{lllrrrr}", "\\toprule",
          "Statistic & Markers & Distance & Sorted & Matched & Difference [95\\% CI] & $n$ \\\\", "\\midrule",
          r[, sprintf("%s & %s & %s & %.3f & %.3f & %s & %d \\\\",
                      ifelse(duplicated(stat), "", slab[as.character(stat)]),
                      ifelse(duplicated(paste(stat, marker_set)), "", marker_set),
                      tex_bin(bin), anchor, control, ci(diff, lo, hi), n_anchors)],
          "\\bottomrule", "\\end{tabular}}", "\\end{table}")
write_tex(t_nb, "neighbourhood.tex")

## ---- unmatched sorted loci -------------------------------------------------------------
um <- copy(nb$unmatched)[order(-span_kb)]
t_um <- c("\\begin{table}[htbp]", "\\centering", "\\small",
          sprintf("\\caption{The %d sorted loci for which no unsorted locus fell within the matching caliper, ordered by the physical span of their own LD unit.}", nrow(um)),
          "\\label{tab:unmatched}",
          "\\begin{tabular}{llrlrrrr}", "\\toprule",
          "Unit & Chr. & Position (Mb) & Sorted toward & $F_{ST}$ & cM/Mb & SNPs & Span (kb) \\\\", "\\midrule",
          um[, sprintf("%s & %s & %.2f & %s & %.3f & %.2f & %d & %.0f \\\\", group_id, sub("Chr", "", Chr), Pos / 1e6,
                       ifelse(sort_class == "aquilonia", "\\textit{F. aquilonia}", "\\textit{F. polyctena}"),
                       FST, recomb_rate, n_loci_g, span_kb)],
          "\\bottomrule", "\\end{tabular}", "\\end{table}")
write_tex(t_um, "unmatched.tex")

## ---- figures -------------------------------------------------------------------------------
figs <- c("03_decay_cM.pdf", "03_decay_bp.pdf", "04_neighbourhood.pdf", "04_balance.pdf", "06_lowDI_contrast.pdf")
for (f in figs) {
  src <- file.path(OUT_FIG, f)
  if (!file.exists(src)) stop("missing figure: ", src)
  file.copy(src, file.path(DOC, "figures", f), overwrite = TRUE); cat("[05] copied", f, "\n")
}
