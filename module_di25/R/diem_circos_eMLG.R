## =========================================================
## module_di25 (high-DI analyses) -- DIEM circos plot, per-eMLG (LD-reduced)
## =========================================================
## The LD-reduced counterpart of diem_circos.R: instead of one ring-cell per
## diagnostic SNP, each cell is an LD-cluster eMLG consensus. This collapses
## linked diagnostic markers into their cluster and shows the same picture at
## the level of (approximately) independent units.
##
## Mapping LD-clusters to the DI25 data:
##   * eMLG_sorted_cM05.rds        -- $eMLG (164 hybrids x 96,461 clusters, dosage
##                                    0..2) + $groups (representative, members).
##   * eMLG_sorted_cM05_parents.rds -- 30 parent individuals x the SAME 96,461
##                                    clusters. rbind -> 194 individuals.
##   * A cluster is kept if any of its member markers is a DI25 diagnostic marker
##     (paste0("Chr", chr, ":", pos) from the SNP file). -> ~5,178 clusters.
##
## Orientation: the eMLG consensus is oriented to an arbitrary within-cluster
##   reference allele, so each cluster is re-oriented to F. aquilonia using the
##   parent side (mean Faqu dosage >= mean Fpol dosage => keep, else 2 - x).
##   Oriented dosage a in [0,2] with 2 = aquilonia; rounded and recoded to the
##   shared DIEM scheme 1=aqu, 2=het, 3=pol via 3 - round(a).
##
## Run from the repo root:  Rscript module_di25/R/diem_circos_eMLG.R
## =========================================================

suppressMessages(library(data.table))
source("module_di25/R/diem_circos_core.R")

SORTED_H <- "moduleA_sorting/data/eMLG_sorted_cM05.rds"
SORTED_P <- "moduleA_sorting/data/eMLG_sorted_cM05_parents.rds"
DI25     <- "data/species_diagnostic_markers_DI25_20pops.tsv.gz"
OUTPNG   <- "module_di25/Figures/diem_circos_eMLG.png"
DI_LABEL <- -25

## ---- load eMLG consensus (hybrid + parent) ------------------------------
message("[diem-circos eMLG] loading eMLG consensus")
s    <- readRDS(SORTED_H)
H    <- s$eMLG                                   # 164 hybrids x clusters
g    <- s$groups
P    <- readRDS(SORTED_P)                        # 30 parents  x clusters
stopifnot(identical(colnames(H), colnames(P)))
E    <- rbind(H, P)                              # 194 individuals x clusters

## ---- keep clusters that carry a DI25 diagnostic marker ------------------
d      <- fread(DI25, select = c("chromosome", "position"))
di_ids <- paste0("Chr", sub("chromosome_", "", d$chromosome), ":", d$position)
mm     <- data.table(group_id = rep(g$group_id, lengths(g$members)),
                     marker   = unlist(g$members))
keep   <- unique(mm[marker %in% di_ids, group_id])
gk     <- g[match(keep, group_id)]
E      <- E[, keep, drop = FALSE]
message(sprintf("  %d clusters carry a DI25 marker (from %d hybrids + %d parents)",
                length(keep), nrow(H), nrow(P)))

## ---- orient each cluster to F. aquilonia via the parent side ------------
faqu <- grep("^Faqu", rownames(E)); fpol <- grep("^Fpol", rownames(E))
maqu <- colMeans(E[faqu, , drop = FALSE], na.rm = TRUE)
mpol <- colMeans(E[fpol, , drop = FALSE], na.rm = TRUE)
flip <- which(maqu < mpol)
E[, flip] <- 2 - E[, flip]                       # now 2 = aquilonia allele everywhere

## ---- recode to DIEM scheme (1 aqu, 2 het, 3 pol; 0 missing) --------------
code <- 3L - round(E)                            # 2->1 aqu, 1->2 het, 0->3 pol
code[is.na(code)] <- 0L
storage.mode(code) <- "integer"

## ---- positions: representative marker of each kept cluster --------------
rep_chr <- as.integer(sub("Chr", "", sub(":.*", "", gk$representative)))
rep_pos <- as.integer(sub(".*:", "", gk$representative))

## units x individuals, ordered by chr/pos (rows) and hybrid index (cols)
gt <- t(code)                                    # clusters x individuals
ord_u <- order(rep_chr, rep_pos)
gt <- gt[ord_u, ]; chr_num <- rep_chr[ord_u]

## hybrid index on the shared 1=aqu..3=pol code: 0 = aquilonia, 1 = polyctena
hi <- apply(gt, 2, function(x) { x <- x[x > 0]; if (!length(x)) NA_real_ else mean((x - 1) / 2) })
gt <- gt[, order(hi)]                            # inner = most aquilonia (matches SNP plot)

n_miss <- sum(gt == 0L)
ttl <- sprintf("LD-reduced (eMLG) @ DI = %d: %s clusters, %s missing (%.1f%%) dithered",
               DI_LABEL, format(nrow(gt), big.mark = ","),
               format(n_miss, big.mark = ","), 100 * n_miss / length(gt))

message("[diem-circos eMLG] rendering -> ", OUTPNG)
render_diem_circos(gt, chr_num, OUTPNG, title = ttl)
message("[diem-circos eMLG] done")
