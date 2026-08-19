## =========================================================
## module_di25 (high-DI analyses) -- is Sielva already "pruning" alleles at sorted loci?
## =========================================================
## Sielva is an F1-like colony (~84% heterozygous genome-wide). Prediction: even
## so, the strongly-sorted loci (fixed in the other 19 populations) may already
## have lost heterozygosity in Sielva -- toward the SAME allele the others fix.
## Tests Sielva heterozygosity vs how sorted a locus is in the other 19 pops,
## controlling for recombination (sorted loci sit in low-recombination regions,
## where any hybrid is more homozygous). Additive: reads di25_* outputs only.
##
## NB Sielva = one colony (Siel338), so its "genotype" is a single diploid.
##
## Run from the repo root:  Rscript module_di25/R/di25_sielva_pruning.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
devtools::load_all("~/gitlab/LDscnR/")

PHI <- 0.85
FIGDIR <- "module_di25/Figures"; OUTDIR <- "module_di25/data"
COL_HET <- "#BDBDBD"; COL_TOW <- "#315B7D"; COL_AGA <- "#C0392B"

## ---- inputs + oriented per-unit dosage ---------------------------------
inp <- readRDS(file.path(OUTDIR, "di25_inputs.rds"))
g   <- readRDS(file.path(OUTDIR, "di25_clustering_cM5.rds"))$groups
e2  <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2); sd <- e2$sample_data_with_parents
GTs <- rbind(inp$GTs_hyb, inp$GTs_par)
faqu <- grep("^Faqu", rownames(GTs)); fpol <- grep("^Fpol", rownames(GTs))
pops <- sd$Population[match(rownames(GTs), sd$Sample_ID)]
hp <- setdiff(unique(pops[!is.na(pops)]), c("aquilonia_parent", "polyctena_parent"))
is_e <- g$n_loci > 2
D <- vapply(seq_len(nrow(g)), function(i) if (is_e[i]) consensus_dosage(GTs, g$members[[i]]) else GTs[, g$representative[i]],
            numeric(nrow(GTs)))
D <- t(D)
flip <- which(rowMeans(D[, faqu, drop = FALSE], na.rm = TRUE) < rowMeans(D[, fpol, drop = FALSE], na.rm = TRUE))
D[flip, ] <- 2 - D[flip, ]                              # 2 = aquilonia

## ---- sorting strength in the OTHER 19 populations ----------------------
PF <- sapply(hp, function(p) rowMeans(D[, which(pops == p), drop = FALSE], na.rm = TRUE) / 2)
oth <- setdiff(seq_along(hp), which(hp == "Sielva"))
St <- ifelse(PF >= PHI, "aqu", ifelse(PF <= 1 - PHI, "pol", "mid"))
npol <- rowSums(St[, oth] == "pol", na.rm = TRUE); naqu <- rowSums(St[, oth] == "aqu", na.rm = TRUE)
nval <- rowSums(!is.na(St[, oth]))
sortedness <- (npol + naqu) / nval                     # fraction of other 19 near-fixed
majority   <- ifelse(npol >= naqu, "pol", "aqu")       # which allele the others fix

## ---- Sielva (single colony) genotype + recombination -------------------
sd_dose <- D[, which(pops == "Sielva")[1]]
siel <- ifelse(sd_dose >= 1.5, "hom_aqu", ifelse(sd_dose <= 0.5, "hom_pol", "het"))
rec <- fread("data/Frufa_DTOL_PR.ref_genome.recmap"); setnames(rec, 1:4, c("chr","p","cM","cMMb")); rec[, Chr := sub("chromosome_","Chr",chr)]
rchr <- as.integer(sub("Chr","",sub(":.*","",g$representative))); rpos <- as.integer(sub(".*:","",g$representative))
recomb <- rep(NA_real_, nrow(g))
for (cc in unique(rchr)) { r <- rec[Chr == paste0("Chr",cc)]; if (nrow(r) < 2) next
  i <- which(rchr == cc); recomb[i] <- approx(r$p, r$cMMb, xout = rpos[i], rule = 2)$y }

dt <- data.table(sortedness, majority, siel, recomb, chr = rchr, n_loci = g$n_loci)[
  is.finite(sortedness) & is.finite(recomb) & !is.na(siel)]
dt[, `:=`(het = as.integer(siel == "het"),
          dir = fifelse(siel == "het", "het",
                fifelse((majority == "pol" & siel == "hom_pol") | (majority == "aqu" & siel == "hom_aqu"),
                        "hom_toward", "hom_against")))]

## ---- stats -------------------------------------------------------------
m <- glm(het ~ sortedness + log10(recomb + 0.1), data = dt, family = binomial())
cat(sprintf("Genome-wide Sielva heterozygosity: %.1f%%\n", 100 * mean(dt$het)))
cat("Logistic Sielva_het ~ sortedness + log10(recomb):\n"); print(round(summary(m)$coef, 3))
cat("\nDirection of Sielva homozygosity at strongly-sorted loci (sortedness > 0.75):\n")
print(dt[sortedness > 0.75 & siel != "het", table(majority_fixes = majority, Sielva = siel)])
saveRDS(list(data = dt, glm = coef(summary(m))), file.path(OUTDIR, "di25_sielva_pruning.rds"))

## ---- figure ------------------------------------------------------------
dt[, sb := cut(sortedness, c(-.01, .25, .5, .75, .95, 1.01),
               labels = c("0-25%", "25-50%", "50-75%", "75-95%", "95-100%"))]
comp <- dt[, .N, by = .(sb, dir)][, frac := N / sum(N), by = sb]
comp[, dir := factor(dir, levels = c("hom_against", "hom_toward", "het"))]
nlab <- dt[, .(N = .N, het = mean(het)), by = sb][order(sb)]
pA <- ggplot(comp, aes(sb, frac, fill = dir)) +
  geom_col(width = 0.8) +
  geom_text(data = nlab, aes(sb, 1.04, label = paste0("n=", N)), inherit.aes = FALSE, size = 2.7) +
  scale_fill_manual(values = c(het = COL_HET, hom_toward = COL_TOW, hom_against = COL_AGA),
                    labels = c(het = "heterozygous", hom_toward = "homozygous, TOWARD the fixed allele",
                               hom_against = "homozygous, against"),
                    breaks = c("het","hom_toward","hom_against"), name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.06))) +
  labs(x = "sorting strength  (fraction of the other 19 populations near-fixed)",
       y = "Sielva's genotype composition",
       title = "Sielva already prunes heterozygosity at sorted loci -- toward the fixed allele",
       subtitle = sprintf("het falls from ~%.0f%% at unsorted loci to ~%.0f%% at strongly-sorted loci; the homozygous fraction is almost entirely TOWARD the allele the others fix.\nSorting effect survives recombination control (logistic sortedness %.2f, p<1e-16; recomb %.2f). The 95-100%% bin is dominated by the 3 huge blocks Sielva still holds as intact F1 het.",
                          100*mean(dt[sortedness<0.25]$het), 100*mean(dt[sortedness>0.75 & sortedness<0.95]$het),
                          coef(m)["sortedness"], coef(m)["log10(recomb + 0.1)"]) ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 8.5))
ggsave(file.path(FIGDIR, "di25_sielva_pruning.png"), pA, width = 11, height = 5.6, dpi = 200)
cat("\n[sielva-pruning] wrote", file.path(FIGDIR, "di25_sielva_pruning.png"), "\n")
