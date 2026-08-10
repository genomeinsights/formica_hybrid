## =========================================================
## module_di25 (high-DI analyses) -- do the 3 polyctena blocks segregate within populations?
## =========================================================
## Two predictions of the co-adaptation/DMI idea, tested descriptively:
##   (1) the blocks should NOT segregate within populations (recombinants inviable)
##       -> within-population dosage SD ~ 0.
##   (2) Sielva is the exception -- it retains aquilonia at the blocks despite being
##       the most F1-like population (heterozygous genome-wide).
## Additive: reads di25_* outputs only.
##
## Run from the repo root:  Rscript module_di25/R/di25_region_segregation.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(ggrepel) })
devtools::load_all("~/gitlab/LDscnR/")

REG <- c(Chr25 = "F6986", Chr26 = "F7174", Chr5 = "F8626")
FIGDIR <- "module_di25/Figures"; OUTDIR <- "module_di25/data"

inp <- readRDS(file.path(OUTDIR, "di25_inputs.rds"))
g   <- readRDS(file.path(OUTDIR, "di25_clustering_cM5.rds"))$groups
e2  <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2); sd <- e2$sample_data_with_parents
GTs <- rbind(inp$GTs_hyb, inp$GTs_par)
faqu <- grep("^Faqu", rownames(GTs)); fpol <- grep("^Fpol", rownames(GTs))
hyb  <- which(rownames(GTs) %in% sd[!grepl("_parent$", Population), Sample_ID])
hpop <- sd$Population[match(rownames(GTs)[hyb], sd$Sample_ID)]

reg <- sapply(REG, function(gid) { a <- consensus_dosage(GTs, g[group_id == gid, members][[1]])
  if (mean(a[faqu], na.rm = TRUE) < mean(a[fpol], na.rm = TRUE)) a <- 2 - a; a })
H <- reg[hyb, c("Chr5", "Chr25", "Chr26")]
het_gw <- rowMeans(GTs[hyb, ] == 1, na.rm = TRUE)          # genome-wide heterozygosity (F1-ness)

## per-population summary: genome-wide het, and per-block mean dosage + within-pop SD
pp <- rbindlist(lapply(c("Chr5", "Chr25", "Chr26"), function(cc) {
  data.table(pop = names(tapply(H[, cc], hpop, mean)), block = cc,
             gwHet = as.numeric(tapply(het_gw, hpop, mean)),
             dose = as.numeric(tapply(H[, cc], hpop, mean)),
             sdw  = as.numeric(tapply(H[, cc], hpop, sd)),
             n    = as.integer(tapply(H[, cc], hpop, length)))
}))
pp[, block := factor(block, levels = c("Chr5", "Chr25", "Chr26"))]
saveRDS(pp, file.path(OUTDIR, "di25_region_segregation.rds"))

cat("max within-population SD at any block:", round(max(pp$sdw, na.rm = TRUE), 2),
    "(", pp[which.max(sdw), paste(pop, block)], ") -- everything else ~0 = no segregation\n")
cat("Sielva genome-wide het:", round(unique(pp[pop == "Sielva", gwHet]), 2),
    "| next highest:", round(sort(unique(pp$gwHet), decreasing = TRUE)[2], 2), "\n")

## ---- figure: block dosage vs genome-wide F1-ness, per population -------
lab <- pp[pop %in% c("Sielva", "Aland", "Nyrhispera75")]
p <- ggplot(pp, aes(gwHet, dose)) +
  geom_errorbar(aes(ymin = dose - sdw, ymax = dose + sdw), width = 0.006, colour = "grey55") +
  geom_point(aes(colour = pop == "Sielva"), size = 2) +
  geom_text_repel(data = lab, aes(label = pop), size = 3, min.segment.length = 0, segment.colour = "grey70") +
  facet_wrap(~ block, nrow = 1) +
  scale_colour_manual(values = c("FALSE" = "#5B6570", "TRUE" = "#C2549D"), guide = "none") +
  scale_y_continuous(limits = c(-0.15, 2), breaks = c(0, 1, 2),
                     labels = c("0\npolyctena", "1\nheterozygous", "2\naquilonia")) +
  labs(x = "genome-wide heterozygosity (F1-ness) of the population",
       y = "aquilonia dosage at the block  (mean +/- within-population SD)",
       title = "The three polyctena blocks do not segregate within populations",
       subtitle = "error bars ~0 = each population is uniform (no recombinant genotypes). Sielva (F1-like, 81% het) uniquely retains the heterozygous state; Aland segregates only at Chr26.") +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 9), strip.text = element_text(face = "bold"))
ggsave(file.path(FIGDIR, "di25_region_segregation.png"), p, width = 12, height = 4.8, dpi = 200)
cat("[region-segregation] wrote", file.path(FIGDIR, "di25_region_segregation.png"), "\n")
