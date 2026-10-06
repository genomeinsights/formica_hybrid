## =========================================================================
## sim_founder_fix -- shared functions for parental / founder LD checks
## (sourced by make_mosaic_founders.R, check_parent_ld.R, make_note_figures.R)
##
## LD statistic: within-species genotype correlation r^2 between the 20,807 DI25
## LD-unit representative SNPs (the unit set of module_allele_specific_sorting),
## for SNPs polymorphic in that species and data set. Classes:
##   cM bins (pairs with identical map position excluded), "unlinked" (all
##   cross-chromosome pairs, exact shortcut), and close pairs (< 5 kb) split by
##   whether both SNPs are in the same min_r2 = 0.2 simulation cluster
##   (group_info_new.rds), which is what the current founder set-up keys on.
## Haploid individuals (aquilonia founder males) are coded 0/2 so that their
## genotype correlation is the gametic r, comparable to the diploid composite r.
## With n individuals, unlinked pairs have E[r^2] ~ 1/(n - 1).
## =========================================================================
suppressMessages({ library(data.table); library(ggplot2) })
theme_set(theme_bw(base_size = 10) + theme(strip.background = element_blank(), panel.grid.minor = element_blank()))

CM_BREAKS <- c(0, 0.001, 0.01, 0.05, 0.2, 1, 5, Inf)
CM_LABELS <- c("<0.001", "0.001-0.01", "0.01-0.05", "0.05-0.2", "0.2-1", "1-5", ">5 cM")
CLASS_LEVELS <- c(CM_LABELS, "unlinked", "<5kb same sim cluster", "<5kb different sim cluster")

## ---- marker panel: all DI25 SNPs, their cM, and the evaluation units --------------
## include_neutral = TRUE adds near-neutral SNPs (DI <= -90 and pooled parental
## minor-allele frequency >= 0.15; ~14,100 SNPs) from the full genotype matrix, as a
## calibration anchor: neutral simulations must reproduce the empirical F_ST at these
## loci (~0.05) before a shortfall at ancestry-informative loci counts as evidence.
load_panel <- function(include_neutral = FALSE) {
  inp <- readRDS("module_di25/data/di25_inputs.rds")
  map <- as.data.table(inp$map)[, Pos := as.integer(Pos)][, class := "DI25"]
  parents <- inp$GTs_par
  if (include_neutral) {
    e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
    sdp <- e$sample_data_with_parents; m <- as.data.table(e$map_hyb_005)
    isp <- grepl("_parent$", sdp$Population[match(rownames(e$GTs_with_parents), sdp$Sample_ID)])
    i <- which(m$DiagnosticIndex <= -90 & !(m$marker %in% map$marker))
    Gp <- e$GTs_with_parents[isp, i, drop = FALSE]; rm(e)
    pf <- colMeans(Gp, na.rm = TRUE) / 2; keep <- pmin(pf, 1 - pf) >= 0.15
    stopifnot("parent IDs differ between the DI25 and full genotype sources" = setequal(rownames(Gp), rownames(parents)))
    Gp <- Gp[rownames(parents), keep, drop = FALSE]
    nm <- m[i[keep], .(Chr, Pos = as.integer(Pos), marker)][, class := "neutral"]
    map <- rbind(map[, .(Chr, Pos, marker, class)], nm)
    parents <- cbind(parents[, map$marker[map$class == "DI25"], drop = FALSE], Gp)
    o <- order(as.integer(sub("Chr", "", map$Chr)), map$Pos); map <- map[o]; parents <- parents[, o, drop = FALSE]
  }
  gm <- fread("data/Frufa_DTOL_PR.ref_genome.recmap"); gm[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
  map[, cM := NA_real_]
  for (ch in intersect(unique(map$Chr), unique(gm$Chr))) {
    g <- gm[Chr == ch & !is.na(cM)]; setorder(g, pos)
    if (nrow(g) < 2) next
    idx <- which(map$Chr == ch); inr <- map$Pos[idx] >= min(g$pos) & map$Pos[idx] <= max(g$pos)
    map$cM[idx] <- approxfun(g$pos, g$cM, rule = 2)(map$Pos[idx])
    map$cM[idx][!inr] <- NA                           # outside the map: no genetic distance
  }
  units <- readRDS("module_allele_specific_sorting/data/01_units.rds")[, .(Chr, Pos, unit_marker)]
  gi <- as.data.table(readRDS("group_info_new.rds"))
  gi[, mk := paste0("Chr", sub("chromosome_", "", chromosome), ":", as.integer(position))]
  units[, simcl := gi$group_id[match(unit_marker, gi$mk)]]
  units[, col := match(unit_marker, map$marker)]; stopifnot(!anyNA(units$col))
  units[, cM := map$cM[col]]
  list(map = map, units = units, parents = parents)
}

## pair structure for the evaluation units (computed once)
make_pairs <- function(units) {
  cols <- split(seq_len(nrow(units)), factor(units$Chr, levels = unique(units$Chr)))
  pr <- lapply(cols, function(ix) {
    m <- length(ix); ut <- which(upper.tri(diag(m)))
    i <- ix[row(diag(m))[ut]]; j <- ix[col(diag(m))[ut]]
    dcm <- abs(units$cM[j] - units$cM[i])
    cm <- cut(dcm, CM_BREAKS, labels = FALSE, right = FALSE); cm[which(dcm == 0)] <- NA
    close <- abs(units$Pos[j] - units$Pos[i]) < 5e3
    list(ut = ut, cm = cm, close = close, same = units$simcl[i] == units$simcl[j])
  })
  list(cols = cols, pairs = pr)
}

## LD profile for one species' genotype matrix (individuals x evaluation units, dosage 0/1/2)
ld_profile <- function(G, PS) {
  X <- sweep(G, 2, colMeans(G, na.rm = TRUE)); X[is.na(X)] <- 0
  nrm <- sqrt(colSums(X^2)); poly <- nrm > 0
  Z <- sweep(X, 2, ifelse(poly, nrm, NA), "/")
  w <- rbindlist(lapply(seq_along(PS$cols), function(k) {
    p <- PS$pairs[[k]]; r2 <- crossprod(Z[, PS$cols[[k]], drop = FALSE])[p$ut]^2
    rbind(data.table(cls = CM_LABELS[p$cm], r2 = r2)[!is.na(cls) & !is.na(r2), .(s = sum(r2), n = .N), by = cls],
          data.table(cls = ifelse(p$same, "<5kb same sim cluster", "<5kb different sim cluster")[p$close],
                     r2 = r2[p$close])[!is.na(r2) & !is.na(cls), .(s = sum(r2), n = .N), by = cls])
  }))[, .(s = sum(s), n = sum(n)), by = cls]
  Z0 <- Z; Z0[is.na(Z0)] <- 0
  M <- lapply(PS$cols, function(ix) tcrossprod(Z0[, ix, drop = FALSE])); np <- sapply(PS$cols, function(ix) sum(poly[ix]))
  half <- function(tot, parts) (tot - parts) / 2
  out <- rbind(w, data.table(cls = "unlinked", s = half(sum(Reduce(`+`, M)^2), sum(sapply(M, function(m) sum(m^2)))),
                             n = half(sum(np)^2, sum(np^2))))
  out[, .(cls = factor(cls, levels = CLASS_LEVELS), r2 = s / n, n_pairs = n)][order(cls)]
}
sample_info <- function(G) data.table(n_ind = nrow(G), n_poly = sum(apply(G, 2, function(v) length(unique(v[!is.na(v)])) > 1)),
                                      het = mean(G == 1, na.rm = TRUE), n_without_het = sum(rowSums(G == 1, na.rm = TRUE) == 0))

## ---- readers: every reader returns list(aquilonia = G, polyctena = G) on the evaluation units
read_empirical_parents <- function(P, units) {
  G <- P[, units$unit_marker]
  list(aquilonia = G[grepl("^Faqu", rownames(G)), ], polyctena = G[grepl("^Fpol", rownames(G)), ])
}
read_diem_parents <- function(bed, units) {
  h <- readLines(bed, n = 2); inds <- strsplit(strsplit(h[2], "\t")[[1]][10], "|", fixed = TRUE)[[1]]
  s <- fread(bed, skip = 2, header = FALSE, sep = "\t", select = c(1, 3, 10), colClasses = list(character = c(1, 10)), showProgress = FALSE)
  mk <- paste0("Chr", sub("ch", "", s$V1), ":", s$V3); keep <- match(units$unit_marker, mk)
  par <- grep("^(aq|pol)_", inds)
  S <- do.call(rbind, strsplit(sub("^S", "", s$V10[keep[!is.na(keep)]]), "", fixed = TRUE))[, par, drop = FALSE]
  G <- matrix(NA_integer_, length(par), nrow(units)); G[, !is.na(keep)] <- t(matrix(suppressWarnings(as.integer(S)), nrow(S)))
  sp <- grepl("^aq_", inds[par])
  list(aquilonia = G[sp, , drop = FALSE], polyctena = G[!sp, , drop = FALSE])
}
## founder VCF directory (founders_ch<id>.vcf; aq_* haploid "0"/"1", pol_* diploid "a|b").
## n_sample: evaluate a random subsample per species (default 15 = the empirical n).
read_founder_vcfs <- function(dir, units, n_sample = 15L, seed = 1L) {
  fs <- list.files(dir, pattern = "^founders_ch[0-9]+\\.vcf$", full.names = TRUE)
  v <- rbindlist(lapply(fs, function(f) fread(f, skip = "#CHROM", sep = "\t", colClasses = "character")))
  setnames(v, 3, "ID")
  keep <- match(units$unit_marker, v$ID)
  smp <- names(v)[-(1:9)]
  dos <- function(x) { x <- as.character(x); ifelse(grepl("\\|", x), as.integer(substr(x, 1, 1)) + as.integer(substr(x, 3, 3)), 2L * as.integer(x)) }
  set.seed(seed)
  pick <- function(rx) { s <- grep(rx, smp, value = TRUE); if (length(s) > n_sample) sample(s, n_sample) else s }
  getG <- function(ss) { G <- matrix(NA_integer_, length(ss), nrow(units))
                        G[, !is.na(keep)] <- t(sapply(ss, function(s) dos(v[[s]][keep[!is.na(keep)]]))); G }
  list(aquilonia = getG(pick("^aq_")), polyctena = getG(pick("^pol_")))
}
