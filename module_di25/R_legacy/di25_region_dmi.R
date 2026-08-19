## =========================================================
## module_di25 (high-DI analyses) -- LD / missing genotypes between the 3 polyctena blocks
## =========================================================
## Tests the epistasis/DMI idea for the three repeatably-polyctena blocks
## (Chr 5, Chr 25, Chr 26): do their rare F. aquilonia variants co-occur across
## individuals (association / "LD" between regions), and are particular three-
## region genotype combinations missing? Additive: reads di25_* outputs only.
##
## The blocks are near-fixed for polyctena, so only a handful of hybrids carry
## aquilonia at each. The key control is POPULATION: if the co-occurring carriers
## are one population, the "association" is structure, not individual-level
## epistasis. (They are -- see the figure.)
##
## Run from the repo root:  Rscript module_di25/R/di25_region_dmi.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
devtools::load_all("~/gitlab/LDscnR/")             # consensus_dosage()

CARRIER_TH <- 0.5                                   # aquilonia-carrier if oriented dosage >= this
REG <- c(Chr25 = "F6986", Chr26 = "F7174", Chr5 = "F8626")   # top polyctena cluster per chr
FIGDIR <- "module_di25/Figures"; OUTDIR <- "module_di25/data"

## ---- inputs -------------------------------------------------------------
inp <- readRDS(file.path(OUTDIR, "di25_inputs.rds"))
g   <- readRDS(file.path(OUTDIR, "di25_clustering_cM5.rds"))$groups
e2  <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e2); sd <- e2$sample_data_with_parents
GTs <- rbind(inp$GTs_hyb, inp$GTs_par)
faqu <- grep("^Faqu", rownames(GTs)); fpol <- grep("^Fpol", rownames(GTs))
hyb  <- which(rownames(GTs) %in% sd[!grepl("_parent$", Population), Sample_ID])
hpop <- sd$Population[match(rownames(GTs)[hyb], sd$Sample_ID)]

## per-individual oriented aquilonia dosage for each region (2 = aquilonia)
reg <- sapply(REG, function(gid) {
  a <- consensus_dosage(GTs, g[group_id == gid, members][[1]])
  if (mean(a[faqu], na.rm = TRUE) < mean(a[fpol], na.rm = TRUE)) a <- 2 - a
  a
})
colnames(reg) <- names(REG)
H <- reg[hyb, c("Chr5", "Chr25", "Chr26")]          # 164 hybrids x 3 regions
carrier <- (H >= CARRIER_TH) * 1L
hi <- 1 - rowMeans(GTs[hyb, ], na.rm = TRUE) / 2     # genome-wide aquilonia index

## ---- association between regions (raw, partial, and Sielva-excluded) -----
r_raw  <- cor(H, use = "pairwise")
r_part <- cor(apply(H, 2, function(y) resid(lm(y ~ hi))))
noSiel <- hpop != "Sielva"
r_noS  <- suppressWarnings(cor(H[noSiel, ], use = "pairwise"))
## the decisive control: drop the only two populations carrying any aquilonia here
keep2  <- !hpop %in% c("Sielva", "Nyrhispera75")
r_drop2 <- cor(H[keep2, "Chr5"], H[keep2, "Chr25"])

## ---- three-region carrier-combination table (obs vs expected) -----------
combos <- apply(expand.grid(Chr5 = 0:1, Chr25 = 0:1, Chr26 = 0:1), 1, paste, collapse = "")
key <- apply(carrier[, c("Chr5", "Chr25", "Chr26")], 1, paste, collapse = "")
obs <- table(factor(key, levels = combos))
mf  <- colMeans(carrier[, c("Chr5", "Chr25", "Chr26")]); n <- length(hyb)
expd <- vapply(combos, function(k) { b <- as.integer(strsplit(k, "")[[1]])
  n * prod(ifelse(b == 1, mf, 1 - mf)) }, numeric(1))
ctab <- data.table(combo = combos, observed = as.integer(obs), expected = round(expd, 2))

saveRDS(list(r_raw = r_raw, r_partial = r_part, r_noSielva = r_noS, r_drop2 = r_drop2, combos = ctab,
             carrier_pops = lapply(colnames(carrier), function(cc) table(hpop[carrier[,cc]==1]))),
        file.path(OUTDIR, "di25_region_dmi.rds"))

cat("== association (Pearson r) between regional aquilonia dosage ==\n")
cat("raw:\n"); print(round(r_raw, 2))
cat("\ncontrolling for genome-wide ancestry:\n"); print(round(r_part, 2))
cat(sprintf("\nChr5-Chr25 r: all %.2f | drop Sielva+Nyrhispera75 %.2f  <- association is entirely those 2 pops\n",
            r_raw["Chr5","Chr25"], r_drop2))
cat("\n== three-region carrier combinations (obs vs expected under independence) ==\n"); print(ctab)

## ---- figure -------------------------------------------------------------
## (a) the carrier individuals, region ancestry, grouped by population
carr_any <- which(rowSums(carrier) > 0)
dt <- as.data.table(reshape2::melt(H[carr_any, , drop = FALSE]))
setnames(dt, c("ind", "region", "aqu_dose"))
dt[, pop := hpop[carr_any][match(ind, rownames(H)[carr_any])]]
ord_ind <- rownames(H)[carr_any][order(hpop[carr_any], -rowSums(H[carr_any, ]))]
dt[, ind := factor(ind, levels = rev(ord_ind))]
dt[, region := factor(region, levels = c("Chr5", "Chr25", "Chr26"))]
pA <- ggplot(dt, aes(region, ind, fill = aqu_dose)) +
  geom_tile(colour = "grey85") +
  facet_grid(pop ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = "#D3C93B", mid = "#440154", high = "#21918C",
                       midpoint = 1, limits = c(0, 2), name = "aquilonia\ndosage") +
  labs(x = NULL, y = NULL, title = "a  Who carries aquilonia at the three polyctena blocks",
       subtitle = "each row an individual; carriers cluster BY POPULATION -> association is structure") +
  theme_minimal(base_size = 10) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8),
        panel.grid = element_blank(), axis.text.y = element_blank(),
        strip.placement = "outside", legend.position = "right")

## (b) combination counts, observed vs expected
cb <- copy(ctab); cb[, lab := c("000"="pol / pol / pol", "100"="AQU / pol / pol",
  "010"="pol / AQU / pol", "001"="pol / pol / AQU", "110"="AQU / AQU / pol",
  "101"="AQU / pol / AQU", "011"="pol / AQU / AQU", "111"="AQU / AQU / AQU")[combo]]
cb[, lab := factor(lab, levels = lab[order(observed)])]
pB <- ggplot(cb, aes(lab, observed)) +
  geom_col(fill = "#315B7D", width = 0.7) +
  geom_point(aes(y = expected), colour = "#E69F00", size = 2.4) +
  geom_text(aes(label = observed), hjust = -0.5, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  coord_flip() +
  labs(x = "Chr5 / Chr25 / Chr26 carrier combination", y = "number of hybrids",
       title = "b  Which genotype combinations exist",
       subtitle = "bars observed, orange expected if independent; mixed classes absent") +
  theme_bw(base_size = 10) + theme(panel.grid.major.y = element_blank())

fig <- (pA + pB + plot_layout(widths = c(1, 1.15))) +
  plot_annotation(
    title = "The three polyctena blocks: aquilonia variation is population-structured -- between-region LD cannot be tested",
    subtitle = sprintf("Only 3 populations carry any aquilonia (Aland: Chr26; Nyrhispera75: Chr5; Sielva: all three). Chr5-Chr25 r = %.2f -> %.2f without those two pops = population, not epistasis.",
                       r_raw["Chr5","Chr25"], r_drop2),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 9.5)))
ggsave(file.path(FIGDIR, "di25_region_dmi.png"), fig, width = 13.5, height = 6, dpi = 200)
cat("\n[region-dmi] wrote", file.path(FIGDIR, "di25_region_dmi.png"), "\n")
