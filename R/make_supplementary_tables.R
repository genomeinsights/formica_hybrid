## =========================================================
## SUPPLEMENTARY TABLES -> LaTeX fragments (manuscript/tables/T*.tex)
## =========================================================
## Generates the supplementary tables from saved analysis outputs as self-contained
## booktabs \table floats, \input by manuscript/supplementary_tables.tex.
## Reads the module .rds; writes manuscript/tables/T1..T8.tex.
## Run from the repo root.

suppressMessages(library(data.table))
dir.create("manuscript/tables", showWarnings = FALSE, recursive = TRUE)
ci  <- function(a, b, d = 2) sprintf("[%.*f, %.*f]", d, a, d, b)
big <- function(x) formatC(x, format = "d", big.mark = ",")

mktab <- function(header, rows, caption, label, colspec, note = NULL) {
  out <- c("\\begin{table}[htbp]", "\\centering", "\\small",
           sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
           sprintf("\\begin{tabular}{%s}", colspec), "\\toprule",
           paste(paste(header, collapse = " & "), "\\\\"), "\\midrule",
           rows, "\\bottomrule", "\\end{tabular}",
           if (!is.null(note)) sprintf("\\par\\smallskip\\footnotesize\\raggedright %s", note),
           "\\end{table}")
  out
}
wr <- function(name, lines) writeLines(lines, sprintf("manuscript/tables/%s.tex", name))

## ---- T1: sort-class counts ----------------------------------------------
snp <- as.data.table(readRDS("data/moduleA_snp.rds"))
g <- snp[differentiated == TRUE]
n_na <- g[is.na(sort_class), .N]; go <- g[!is.na(sort_class)]; tot <- nrow(go)
cnt <- function(cl) go[sort_class == cl, .N]
r1 <- c(
  sprintf("Unsorted & %s & %.2f \\\\", big(cnt("unsorted")), 100*cnt("unsorted")/tot),
  sprintf("Unidirectional (\\emph{polyctena}) & %s & %.2f \\\\", big(cnt("polyctena")), 100*cnt("polyctena")/tot),
  sprintf("Unidirectional (\\emph{aquilonia}) & %s & %.2f \\\\", big(cnt("aquilonia")), 100*cnt("aquilonia")/tot),
  sprintf("Bidirectional & %s & %.2f \\\\", big(cnt("bidirectional")), 100*cnt("bidirectional")/tot),
  "\\midrule",
  sprintf("No defined parental orientation (excluded) & %s & --- \\\\", big(n_na)))
wr("T1", mktab(c("Sort class", "SNPs", "\\% of oriented"), r1,
  sprintf("Genome-wide sort-class counts among parent-polymorphic SNPs (pooled parental MAF $\\geq0.15$; %s SNPs). Percentages are among the %s loci with a defined parental orientation.", big(nrow(g)), big(tot)),
  "tab:sortclass", "l r r"))

## ---- T2: directional-predictability test summary ------------------------
cl <- as.data.table(readRDS("data/moduleA_clusters.rds"))
t2row <- function(x, lab) { d <- x[differentiated == TRUE & n_fixed >= 5]
  sprintf("%s & %s & %s & %.1f \\\\", lab, big(nrow(d)), big(sum(d$sig, na.rm = TRUE)),
          100*sum(d$sig, na.rm = TRUE)/nrow(d)) }
r2 <- c(t2row(snp, "SNP"), t2row(cl, "LD-reduced unit"))
wr("T2", mktab(c("Level", "Tested", "Significant", "\\%"), r2,
  "Directional-predictability test (two-sided exact binomial vs $p=0.5$; units near-fixed in $\\geq5$ populations; Benjamini--Hochberg adjusted, $q<0.05$).",
  "tab:predtest", "l r r r",
  "Tests are conditional and exploratory. The much higher significant fraction at the unit level reflects the parental-polymorphism gate and the larger number of near-fixed populations required, not stronger evidence per locus."))

## ---- T3: cluster-size aggregation + dilution ----------------------------
sz <- as.data.table(readRDS("data/moduleA_di_asymmetry.rds")$size_tab)
r3 <- apply(sz, 1, function(row) sprintf("%s & %s & %.1f \\\\", row[["size_class"]],
              big(as.integer(row[["n"]])), as.numeric(row[["pct_unsorted"]])))
dil <- as.data.table(readRDS("data/moduleA_dilution.rds"))
du <- dil[sort_class == "unsorted" & n_loci >= 2]
wr("T3", mktab(c("Cluster size (markers)", "Clusters", "\\% unsorted after aggregation"), r3,
  "Aggregation of marker-level sorting by cluster size (clusters entering because they contained $\\geq1$ individually sorted marker). The comparison is conditional and does not estimate the genome-wide unsorted fraction.",
  "tab:sizeagg", "l r r",
  sprintf("Dilution check: among the %s multi-marker clusters unsorted after aggregation, consensus fidelity remained high (median $\\mathrm{score}_{\\mathrm{eMLG}}=%.3f$); only %d fell below the 0.80 fidelity floor, so the loss of the aggregate call reflects genuine absence of a joint signal rather than averaging artefact.",
          big(nrow(du)), median(du$score, na.rm = TRUE), sum(du$score < 0.80, na.rm = TRUE))))

## ---- T4: rank correlations (documented final Module B B1 values) --------
r4 <- c(
  "Diagnostic index & Recombination rate & $-0.03$ \\\\",
  "Diagnostic index & Pooled parental MAF & $-0.02$ \\\\",
  "Diagnostic index & Parental expected heterozygosity ($\\pi$) & $-0.46$ \\\\",
  "Diagnostic index & Between-parent allelic difference ($d_{xy}$) & $+0.13$ \\\\")
wr("T4", mktab(c("Variable", "Variable", "Spearman $\\rho$"), r4,
  "Genome-wide rank correlations at analysed SNPs. The diagnostic index is essentially orthogonal to recombination and to pooled parental MAF, and tracks reduced parental expected heterozygosity rather than elevated between-parent difference.",
  "tab:corr", "l l r"))

## ---- T5: architecture regression tables (saved models) ------------------
b <- readRDS("data/moduleB_architecture.rds")
coefrows <- function(m, val = "z value", terms) { s <- coef(summary(m))
  vapply(names(terms), function(k) sprintf("%s & %.3f & %.3f & %.1f & %s \\\\", terms[[k]],
      s[k, 1], s[k, 2], s[k, val], if (s[k, 4] < 1e-3) "$<0.001$" else sprintf("%.3f", s[k, 4])),
    character(1)) }
r5a <- coefrows(b$lm_magnitude, "t value",
  list(zr = "Recombination (std.\\ log)", zDI = "Diagnostic index (std.)"))
r5b <- coefrows(b$glm_direction, "z value",
  list(zDI = "Diagnostic index (std.)", zr = "Recombination (std.\\ log)",
       zmaf = "Pooled parental MAF (std.)", zcs = "Cluster size (std.\\ log)"))
wr("T5", c(
  mktab(c("Predictor", "Coef.", "SE", "$t$", "$P$"), r5a,
    "Sorting magnitude ($\\mathrm{prop\\_fixed}$) at the LD-reduced-unit level (linear model). Cluster size is excluded as a potential collider.",
    "tab:magnitude", "l r r r r"),
  "",
  mktab(c("Predictor", "Coef.", "SE", "$z$", "$P$"), r5b,
    "Direction of sorting (logistic model of \\emph{aquilonia}-directed near-fixation among unidirectional units). Coefficients are associations, interpreted descriptively.",
    "tab:direction", "l r r r r")))

## ---- T6: climate DI-enrichment across strata (primary config) -----------
rb <- readRDS("data/moduleC_rate_based_5_15.rds")
prim <- function(tab, strat_lab, which_test) {
  d <- as.data.table(tab)
  d <- d[popset == "aland_excluded" & omega == "withOmega" &
           predictor == "p_sig" & test == which_test]
  d[, .(pc, strat = strat_lab, effect, p)]
}
di6 <- rbind(prim(rb$full, "all", "DI cluster-level"),
             prim(rb$rob20, "$\\geq20$", "DI cluster-level"),
             prim(rb$rob50, "$\\geq50$", "DI cluster-level"))
r6 <- di6[order(pc), sprintf("%s & %s & %+.2f & %s \\\\", pc, strat, effect,
            ifelse(p < 1e-3, "$<0.001$", sprintf("%.3f", p)))]
so <- rbind(prim(rb$full, "all", "sorting-overlap"),
            prim(rb$rob20, "$\\geq20$", "sorting-overlap"),
            prim(rb$rob50, "$\\geq50$", "sorting-overlap"))
r6b <- so[order(pc), sprintf("%s & %s & %.2f & %s \\\\", pc, strat, effect,
            ifelse(p < 1e-3, "$<0.001$", sprintf("%.3f", p)))]
wr("T6", c(
  mktab(c("Axis", "Stratum", "Effect (pp / +10pp)", "$P$"), r6,
    "Diagnostic-index enrichment: percentage-point change in ancestry-diagnostic content per +10 percentage-point increase in the cluster climate-significance rate (primary configuration: {\\AA}land-excluded, fixed-$\\Omega$; cluster-level linear model with $\\log$ cluster size).",
    "tab:enrich", "l l r r"),
  "",
  mktab(c("Axis", "Stratum", "Odds ratio", "$P$"), r6b,
    "Sorting--climate overlap: cluster-level logistic odds ratio for directional sorting per unit significance rate (primary configuration), showing no association.",
    "tab:overlap", "l l r r")))

# NOTE: the MAF/XtX ladder table and the block-bootstrap CI table were removed from
# the manuscript set to avoid duplicating supplementary figures S11 and S13, whose
# captions now carry the exact values. The generating code stays in
# R/moduleC_maf_power_sensitivity.R and R/sensitivity_block_bootstrap.R.

cat("Wrote manuscript/tables/T1..T6.tex\n")
