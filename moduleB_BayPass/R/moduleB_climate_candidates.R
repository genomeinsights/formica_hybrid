## =========================================================
## MODULE B -- climate-association candidate clusters + eMLG gate
## =========================================================
## Defines the RESTRICTED set of climate-association candidate LD clusters that
## the rest of Module B analyses. It sits AFTER the BayPass runs (moduleB_
## analyse_baypass.R) and BEFORE every downstream climate analysis, because it
## narrows the candidate regions those analyses operate on.
##
## Two-stage candidate definition:
##   (1) member-count entry  -- a cluster is a putative candidate when >= MIN5 of
##       its member markers reach BF(dB) >= SIG_THR on PC1 or PC2 (the set used
##       for visualisation in moduleB_analyse_baypass.R). This is size-dependent:
##       large blocks pass on an absolute count while sitting near the background
##       significant-marker rate, so it over-calls big differentiated regions.
##   (2) eMLG gate           -- of those, keep the clusters whose OWN eMLG (the
##       cluster's consensus genotype, one value per individual) is itself
##       associated with the climate covariate at eMLG BF(dB) >= EMLG_BF_THR.
##       This asks whether the cluster tracks climate AS A UNIT, not just via a
##       handful of individually-significant member SNPs -- and it demotes the
##       "a few lucky SNPs in an LD block" false positives (checked directly:
##       several 30-44% member-%sig Chr2/Chr26 blocks have NEGATIVE eMLG BF).
##       Genome-wide only ~0.16% of eMLGs reach BF(dB) >= 15, so this is a
##       stringent, well-calibrated gate.
##
## The LD-cluster / eMLG association approach (reducing tightly-linked blocks to
## a cluster-level consensus genotype and testing that, rather than individual
## SNPs) follows: Li Z, Kemppainen P, Rastas P, Merila J. 2018. Linkage
## disequilibrium clustering-based approach for association mapping with tightly
## linked genomewide data. Mol Ecol Resour 18:809-824.
##
## SIM-FDR STATUS (R_PK/moduleB_eMLG_null.R, 10,000 Omega-structured null
## covariates; genome-wide sim-FDR = expected floor count / observed, with a
## floor-survivor = eMLG whose real BF beats ALL 10,000 population-structure
## nulls). The 7 member-count candidates below DO NOT survive: each is beaten by
## structure nulls (flagship F63637, BF(dB)=38.4, is exceeded by 12 of 10,000
## structure draws; F63812 by 18). This is EXPECTED and instructive rather than a
## failure of the design -- count-based single-SNP outlier tests are prone to
## differentiation-driven false positives, which is exactly the argument for
## LD-complexity reduction. The signal that best resists structure is instead the
## eMLG-only PC2 set (5 clusters -- e.g. F54421/F53948/F14310 -- whose consensus
## tracks climate although NO member SNP is individually significant; set
## FDR~=0.66, suggestive not significant). It is reachable at all only because the
## eMLG reduction collapses ~1.1M SNP tests to 32,840 cluster tests (~34x lower
## multiple-testing burden, hence higher power per true signal). So the candidate
## set defined here is a member-count/eMLG-BF construct for description; the
## structure-controlled inference lives in moduleB_eMLG_null.rds.
##
## NOTE on eBPis. BayPass's eBPis ("nominal p" = 10^-eBPis) is degenerate here
## -- essentially every one of the 32,840 eMLGs gets p < 0.05 / FDR q < 0.05, so
## it separates nothing. BF(dB) is the discriminating statistic and defines the
## gate; eBPis is retained in the saved object only for completeness.
##
## Primary configuration only: Aland-excluded, fixed LD-pruned Omega (withOmega).
## The eMLG-level association is produced by:
##   R_PK/moduleB_write_eMLG_baypass_inputs.R  -> aland_excluded_eMLG/u_eMLG.geno
##   R_PK/run_baypass_eMLG.sh                  -> aland_excluded_eMLG/PC{1,2}_eMLG_withOmega_summary_betai_reg.out
## This script consumes those .out files (it does not run BayPass itself), the
## same way moduleB_analyse_baypass.R consumes the SNP-level .out files.
##
## Reads : module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds (clustering + eMLG dosages),
##         data/hybrids_only_maf005.Rdata (map: marker order for the positional
##           BF join, Chr/Pos, DiagnosticIndex),
##         {with_aland,aland_excluded}/PC{1,2}_DIEM_{noOmega,withOmega}_summary_{betai_reg,pi_xtx}.out,
##         aland_excluded_eMLG/PC{1,2}_eMLG_withOmega_summary_betai_reg.out,
##         aland_excluded_eMLG/eMLG_group_order.txt
## Writes: moduleB_BayPass/data/moduleB_climate_candidates.rds  (all 35 candidates + pass_eMLG flag;
##           the restricted set is pass_eMLG == TRUE),
##         moduleB_BayPass/doc/tables/climate_candidate_clusters.{txt,tex}  (8-config counts),
##         moduleB_BayPass/doc/tables/climate_candidate_full.tex            (per-cluster + eMLG gate)
## Run from the formica_hybrid repo root.

suppressMessages(library(data.table))

## ---- parameters ---------------------------------------------------------
SIG_THR     <- 15                # marker-level significance: BF(dB) >= 15
MIN5        <- 5; MIN10 <- 10    # member-count entry thresholds
EMLG_BF_THR <- 15                # eMLG gate: cluster-level BF(dB) >= this
DI_TH       <- -25               # ancestry-informative criterion (reported, not a gate)
CLUST       <- "module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds"
OUTDIR      <- "moduleB_BayPass/doc/tables"
EMLG_DIR    <- "aland_excluded_eMLG"
PRIM_POP    <- "aland_excluded"; PRIM_OM <- "withOmega"
## Omega label: noOmega = internal (full data); withOmega = LD-pruned (fixed)
om_label  <- c(noOmega = "Internal (full data)", withOmega = "LD-pruned (fixed)")
pop_label <- c(with_aland = "Full (20 populations)",
               aland_excluded = "\\AA land-excluded (19 populations)")
pop_txt   <- c(with_aland = "Full (20 populations)",
               aland_excluded = "Åland-excluded (19 populations)")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## ---- shared inputs ------------------------------------------------------
groups <- readRDS(CLUST)$groups
he  <- groups[has_eMLG == TRUE]
m2g <- he[, .(marker = unlist(members)), by = group_id]
e <- new.env(); load("data/hybrids_only_maf005.Rdata", envir = e)
map <- as.data.table(e$map_hyb_005)[, .(marker, Chr, Pos, DI = DiagnosticIndex)]
mk  <- map$marker; rm(e); invisible(gc())
m2g <- map[m2g, on = "marker"]
gw_median_size <- median(he$n_loci)

## SNP-level BayPass import: MRK is the row index into the FULL unpruned marker
## set (no marker-ID column), so this is a POSITIONAL join (MRK == row number).
## The stopifnot is the guard against a stale/reordered .out silently
## misaligning BF values onto the wrong markers.
import_col <- function(kind, pop, pc, om, col) {
  f <- sprintf("./%s/%s_DIEM_%s_summary_%s.out", pop, pc, om, kind)
  r <- fread(f)
  stopifnot(nrow(r) == length(mk), identical(r$MRK, seq_len(nrow(r))), col %in% names(r))
  setNames(r[[col]], mk)
}

## eMLG-level BayPass import: one row per eMLG, MRK == row index into
## u_eMLG.geno, whose row order is recorded in eMLG_group_order.txt -> map to
## group_id. Same positional-join guard.
grp_order <- readLines(file.path(EMLG_DIR, "eMLG_group_order.txt"))
import_eMLG <- function(pc) {
  f <- sprintf("%s/%s_eMLG_%s_summary_betai_reg.out", EMLG_DIR, pc, PRIM_OM)
  if (!file.exists(f))
    stop("Missing ", f, " -- run R_PK/moduleB_write_eMLG_baypass_inputs.R then R_PK/run_baypass_eMLG.sh first.")
  r <- fread(f)
  stopifnot(nrow(r) == length(grp_order), identical(r$MRK, seq_len(nrow(r))))
  data.table(group_id = grp_order, BF = r$`BF(dB)`, eBP = r$eBPis, Beta = r$Beta_is)
}

## ========================================================================
## 1 -- putative-candidate summary across the eight configurations
## ========================================================================
## Descriptive context (member-count entry only, no eMLG gate): how many
## clusters pass >=5 / >=10 significant members, per configuration.
cfg <- CJ(pop = c("with_aland", "aland_excluded"), om = c("noOmega", "withOmega"),
          pc = c("PC1", "PC2"), sorted = FALSE)
summ <- rbindlist(lapply(seq_len(nrow(cfg)), function(i) {
  cf <- cfg[i]
  bf <- import_col("betai_reg", cf$pop, cf$pc, cf$om, "BF(dB)")
  ## Hoist the BF lookup OUT of the by-group aggregation: do a single
  ## vectorised name-match over all member markers, then aggregate. Doing
  ## `bf[marker]` inside `by = group_id` instead re-matches against the
  ## 1.11M-name vector once per cluster (x8 configs) -- an O(groups x names)
  ## blowup that runs for hours.
  sig <- bf[m2g$marker] >= SIG_THR
  a <- m2g[, .(n_loci = .N, n_sig = sum(sig[.I], na.rm = TRUE)), by = group_id]
  data.table(pop = cf$pop, om = cf$om, pc = cf$pc,
             out5  = a[n_sig >= MIN5,  .N],
             out10 = a[n_sig >= MIN10, .N],
             med5  = a[n_sig >= MIN5,  as.numeric(median(n_loci))])
}))
cat("=== putative-candidate summary (8 configs, member-count entry) ===\n"); print(summ)

ord <- summ[CJ(pop = c("with_aland", "aland_excluded"), om = c("noOmega", "withOmega"),
               pc = c("PC1", "PC2")), on = c("pop", "om", "pc")]

## ---- .txt (TSV) ---------------------------------------------------------
hdr <- c("Population set", "Covariance matrix (\u03A9)", "Climate axis",
         "Candidate clusters (\u22655 sig. markers)",
         "Candidate clusters (\u226510 sig. markers)",
         "Median candidate cluster size (markers)")
rows_txt <- ord[, sprintf("%s\t%s\t%s\t%d\t%d\t%g",
  pop_txt[pop], ifelse(pop == PRIM_POP & om == PRIM_OM, paste0(om_label[om], " [PRIMARY]"), om_label[om]),
  pc, out5, out10, round(med5))]
writeLines(c("Table [climate-associated cluster candidates]. Candidate climate-associated LD clusters across BayPass configurations.",
  "", paste(hdr, collapse = "\t"), rows_txt, "",
  sprintf("Note: a cluster is a putative candidate when at least 5 (or 10) of its member markers reach BF(dB) >= %d. Genome-wide median cluster size = %g markers (%s clusters with a consensus genotype). The large size of candidate clusters relative to the genome-wide median reflects the size dependence of the member-count criterion; the eMLG gate (see the per-cluster table) restricts this set further.",
          SIG_THR, gw_median_size, format(nrow(he), big.mark = ","))),
  file.path(OUTDIR, "climate_candidate_clusters.txt"))

## ---- shared standalone-LaTeX wrapper ------------------------------------
standalone <- function(caption, tab, note, width = "16.2cm", size = "\\footnotesize") {
  c("\\documentclass[border=12pt]{standalone}", "\\usepackage[T1]{fontenc}", "\\usepackage{lmodern}",
    "\\usepackage{amsmath,booktabs,array,multirow,makecell,xcolor}",
    "\\definecolor{Accent}{HTML}{315B7D}", "\\definecolor{Muted}{HTML}{5B6570}",
    "\\renewcommand{\\arraystretch}{1.18}", "\\begin{document}",
    sprintf("\\begin{minipage}{%s}", width), "\\raggedright",
    sprintf("{\\color{Accent}\\bfseries %s}", caption), "", "\\vspace{6pt}",
    sprintf("\\begin{center}%s", size), tab, "\\end{center}", "",
    sprintf("\\vspace{2pt}{\\scriptsize\\color{Muted} %s\\par}", note),
    "\\end{minipage}", "\\end{document}")
}

## ---- .tex (grouped multirow) --------------------------------------------
sd <- ord
srow <- function(i) {
  r <- sd[i]; omlab <- om_label[r$om]
  if (r$pop == PRIM_POP && r$om == PRIM_OM) omlab <- paste0("\\makecell[l]{", omlab, "\\\\ \\textbf{[primary]}}")
  first_pop <- (i %% 4 == 1); first_om <- (i %% 2 == 1)
  popcell <- if (first_pop) sprintf("\\multirow{4}{*}{%s}", pop_label[r$pop]) else ""
  omcell  <- if (first_om)  sprintf("\\multirow{2}{*}{%s}", omlab) else ""
  line <- sprintf("%s & %s & %s & %d & %d & %g \\\\", popcell, omcell, r$pc, r$out5, r$out10, round(r$med5))
  if (first_om && !first_pop) line <- paste0("\\cmidrule(l){2-6}\n", line)
  if (first_pop && i > 1)     line <- paste0("\\midrule\n", line)
  line
}
tab1 <- c("\\begin{tabular}{@{}l l c r r r@{}}", "\\toprule",
  "Population set & $\\Omega$ (covariance) & Axis & \\makecell{Candidates\\\\($\\geq 5$)} & \\makecell{Candidates\\\\($\\geq 10$)} & \\makecell{Median size\\\\(markers)} \\\\",
  "\\midrule", vapply(seq_len(nrow(sd)), srow, character(1)), "\\bottomrule", "\\end{tabular}")
note1 <- sprintf("A cluster is a putative candidate when at least 5 (or 10) of its member markers reach $\\mathrm{BF(dB)}\\geq %d$. Genome-wide median cluster size $=%g$ markers (%s clusters with a consensus genotype); the excess size of candidate clusters over this median reflects the size dependence of the member-count criterion. The eMLG gate restricts this set to clusters that track climate as a unit (per-cluster table).",
                 SIG_THR, gw_median_size, format(nrow(he), big.mark = ","))
writeLines(standalone(
  "Table [climate-associated cluster candidates]. Putative climate-associated LD clusters across BayPass configurations.",
  tab1, note1), file.path(OUTDIR, "climate_candidate_clusters.tex"))

## ========================================================================
## 2 -- primary per-cluster table + eMLG gate -> restricted candidate set
## ========================================================================
b1 <- import_col("betai_reg", PRIM_POP, "PC1", PRIM_OM, "BF(dB)")
b2 <- import_col("betai_reg", PRIM_POP, "PC2", PRIM_OM, "BF(dB)")
xt <- import_col("pi_xtx",    PRIM_POP, "PC1", PRIM_OM, "M_XtX")
m2g[, `:=`(s1 = b1[marker] >= SIG_THR, s2 = b2[marker] >= SIG_THR,
           XT = xt[marker], hi = DI > DI_TH)]

cl <- m2g[, .(n_loci = .N, Chr = Chr[1], startMb = min(Pos) / 1e6, endMb = max(Pos) / 1e6,
              nsig1 = sum(s1, na.rm = TRUE), nsig2 = sum(s2, na.rm = TRUE),
              XtX = mean(XT, na.rm = TRUE),
              DIpct = 100 * sum(hi, na.rm = TRUE) / sum(!is.na(DI))), by = group_id]
cl[, `:=`(pct1 = 100 * nsig1 / n_loci, pct2 = 100 * nsig2 / n_loci)]

## member-count candidate set (entry), with the axis it entered on
cl[, cand_axis := fifelse(nsig1 >= MIN5 & nsig2 >= MIN5, "both",
                   fifelse(nsig1 >= MIN5, "PC1", fifelse(nsig2 >= MIN5, "PC2", NA_character_)))]
cand <- cl[!is.na(cand_axis)]

## ---- join eMLG-level association (primary config) -----------------------
E1 <- import_eMLG("PC1"); E2 <- import_eMLG("PC2")
cand <- E1[, .(group_id, eBF1 = BF, eBP1 = eBP, eBeta1 = Beta)][cand, on = "group_id"]
cand <- E2[, .(group_id, eBF2 = BF, eBP2 = eBP, eBeta2 = Beta)][cand, on = "group_id"]
## relevant-axis eMLG statistics = the axis the cluster entered on (max over
## both when it entered on both)
cand[, `:=`(
  eMLG_BF   = fifelse(cand_axis == "PC2", eBF2,  fifelse(cand_axis == "PC1", eBF1,  pmax(eBF1, eBF2))),
  eMLG_Beta = fifelse(cand_axis == "PC2", eBeta2, fifelse(cand_axis == "PC1", eBeta1,
                fifelse(eBF1 >= eBF2, eBeta1, eBeta2))),
  memb_pct  = fifelse(cand_axis == "PC2", pct2,  fifelse(cand_axis == "PC1", pct1,  pmax(pct1, pct2))),
  memb_nsig = fifelse(cand_axis == "PC2", nsig2, fifelse(cand_axis == "PC1", nsig1, pmax(nsig1, nsig2)))
)]
cand[, pass_eMLG := eMLG_BF >= EMLG_BF_THR]
setorder(cand, -eMLG_BF)

restricted <- cand[pass_eMLG == TRUE]
cat(sprintf("\n=== eMLG gate (BF(dB) >= %d): %d of %d member-count candidates pass ===\n",
            EMLG_BF_THR, nrow(restricted), nrow(cand)))
print(restricted[, .(group_id, Chr, size = n_loci, cand_axis,
                     memb_pct = round(memb_pct, 1), eMLG_BF = round(eMLG_BF, 1),
                     eMLG_Beta = round(eMLG_Beta, 3))])

## ---- save the candidate object (all candidates + pass flag) -------------
saveRDS(cand, "moduleB_BayPass/data/moduleB_climate_candidates.rds")
cat("\nSaved moduleB_BayPass/data/moduleB_climate_candidates.rds (", nrow(cand),
    " candidates; restricted set = pass_eMLG == TRUE, n=", nrow(restricted), ")\n", sep = "")

## ---- per-cluster LaTeX table (ordered by eMLG BF, gate marked) ----------
esc <- function(x) gsub("_", "\\\\_", x)
frow <- function(r) {
  id <- esc(r$group_id)
  if (r$pass_eMLG) id <- sprintf("\\textbf{%s}", id)          # gate-passing clusters in bold
  axis <- if (r$cand_axis == "both") "PC1+2" else r$cand_axis
  sprintf("%s & %s & %.2f--%.2f & %d & %s & %d & %.1f & %.1f & %+.3f & %.0f & %.1f \\\\",
          id, r$Chr, r$startMb, r$endMb, r$n_loci, axis, r$memb_nsig, r$memb_pct,
          r$eMLG_BF, r$eMLG_Beta, r$XtX, r$DIpct)
}
body <- vapply(seq_len(nrow(cand)), function(i) {
  line <- frow(cand[i])
  ## rule separating gate-passing (top) from failing clusters
  if (i > 1 && cand$pass_eMLG[i - 1] && !cand$pass_eMLG[i]) line <- paste0("\\midrule\n", line)
  line
}, character(1))
tab2 <- c("\\begin{tabular}{@{}l l c r c r r r r r r@{}}", "\\toprule",
  "Cluster & Chr & Position (Mb) & Size & Axis & Sig. & \\% sig. & \\makecell{eMLG\\\\BF(dB)} & \\makecell{eMLG\\\\$\\beta$} & XtX & DI\\% \\\\",
  "\\midrule", body, "\\bottomrule", "\\end{tabular}")
note2 <- sprintf("Ordered by the cluster-level eMLG association on the axis the cluster entered on. Size = member markers; Axis = climate axis of entry (member-count); Sig.\\ = member markers with $\\mathrm{BF(dB)}\\geq %d$; \\%% sig.\\ = Sig./Size; eMLG BF(dB) and $\\beta$ = the cluster's consensus-genotype association from a BayPass run on the eMLGs (primary config); XtX = mean among-population differentiation (a prerequisite for association, uniformly elevated); DI\\%% = proportion of members with $\\mathrm{DI}>%d$. \\textbf{Bold} clusters pass the eMLG gate ($\\mathrm{BF(dB)}\\geq %d$; %d of %d) and constitute the restricted candidate set; the rule separates them from clusters passing on member count alone. Genome-wide median cluster size = %g markers.",
                 SIG_THR, DI_TH, EMLG_BF_THR, nrow(restricted), nrow(cand), gw_median_size)
writeLines(standalone(
  sprintf("Table [climate-associated cluster candidates, full]. Candidate climate-associated LD clusters (primary analysis: {\\AA}land-excluded, LD-pruned $\\Omega$), with the eMLG gate. Bold = restricted set (eMLG $\\mathrm{BF(dB)}\\geq %d$).", EMLG_BF_THR),
  tab2, note2, width = "18.5cm", size = "\\scriptsize"),
  file.path(OUTDIR, "climate_candidate_full.tex"))

## ---- optional: compile the standalone tables to PDF if latexmk present --
if (nzchar(Sys.which("latexmk"))) {
  for (f in c("climate_candidate_clusters", "climate_candidate_full"))
    system(sprintf("cd %s && latexmk -pdf -interaction=nonstopmode %s.tex >/dev/null 2>&1 && latexmk -c %s.tex >/dev/null 2>&1",
                   OUTDIR, f, f))
  invisible(file.remove(Sys.glob(file.path(OUTDIR, "*.aux")),
                        Sys.glob(file.path(OUTDIR, "*.fls")),
                        Sys.glob(file.path(OUTDIR, "*.fdb_latexmk")),
                        Sys.glob(file.path(OUTDIR, "*.log"))))
}
cat("\nWrote moduleB_BayPass/data/moduleB_climate_candidates.rds and moduleB_BayPass/doc/tables/climate_candidate_{clusters.{txt,tex},full.tex}\n")
