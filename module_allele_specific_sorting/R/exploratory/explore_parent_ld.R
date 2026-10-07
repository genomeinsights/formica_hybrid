## =========================================================================
## module_allele_specific_sorting -- EXPLORATORY (not in the document, not committed):
## do LD patterns in the SIMULATED parents match those in the EMPIRICAL parents?
##
## Within-species LD (genotype correlation r^2) between the 20,807 unit
## representative SNPs, separately for F. aquilonia and F. polyctena parents:
##   empirical: the 15 + 15 reference parents (module_di25/data/di25_inputs.rds)
##   simulated: the 15 + 15 aq_/pol_ parents in each diem_boot replicate
## Only SNPs polymorphic within the species (in that data set) enter. Pairs binned
## by bp and by cM (zero-map-distance pairs excluded, as in 01); unlinked class via
## the exact shortcut sum_{a<b} <Z_a Z_a', Z_b Z_b'>. With n = 15, unlinked pairs have
## E[r^2] ~ 1/(n-1) ~ 0.07: compare each curve with its own unlinked level.
## Also reported per species: polymorphic units, heterozygosity, individuals
## without any heterozygote (haploid-like: SLiM's aquilonia founders are males).
## CAVEAT: existing replicates carry the founder-frequency bug.
## Output: data/explore_parent_ld.rds, Figures/explore_parent_ld.png
## Run from the formica_hybrid repo root:
##   Rscript module_allele_specific_sorting/R/exploratory/explore_parent_ld.R [N_REPS] [N_CORES]
## =========================================================================
source("module_allele_specific_sorting/R/00_utils.R")
args <- commandArgs(trailingOnly = TRUE)
N_REPS  <- if (length(args) >= 1) as.integer(args[1]) else 20L
N_CORES <- if (length(args) >= 2) as.integer(args[2]) else 8L
BEDFMT  <- "data/diem_outs_demo/diem_boot%d_output.bed"

u0 <- readRDS(file.path(OUT_DATA, "01_units.rds"))[, .(Chr, Pos, cM_pos, unit_marker)]
chrs <- unique(u0$Chr); cols <- split(seq_len(nrow(u0)), factor(u0$Chr, levels = chrs))
pairs <- lapply(cols, function(ix) {
  m <- length(ix); ut <- which(upper.tri(diag(m)))
  i <- ix[row(diag(m))[ut]]; j <- ix[col(diag(m))[ut]]
  dcm <- abs(u0$cM_pos[j] - u0$cM_pos[i])
  cm <- cut(dcm, CM_BREAKS, labels = FALSE, right = FALSE); cm[which(dcm == 0)] <- NA
  list(ut = ut, bp = cut(abs(u0$Pos[j] - u0$Pos[i]), BP_BREAKS, labels = FALSE, right = FALSE), cm = cm)
})

## LD summary for one species' genotype matrix (individuals x units)
ld_summary <- function(G) {
  X <- sweep(G, 2, colMeans(G, na.rm = TRUE)); X[is.na(X)] <- 0
  nrm <- sqrt(colSums(X^2)); poly <- nrm > 0
  Z <- sweep(X, 2, ifelse(poly, nrm, NA), "/")
  w <- rbindlist(lapply(seq_along(cols), function(k) {
    p <- pairs[[k]]; r2 <- crossprod(Z[, cols[[k]], drop = FALSE])[p$ut]^2
    ok <- !is.na(r2)
    rbind(data.table(scale = "bp", bin = p$bp[ok], r2 = r2[ok])[, .(s = sum(r2), n = .N), by = .(scale, bin)],
          data.table(scale = "cM", bin = p$cm[ok], r2 = r2[ok])[!is.na(bin), .(s = sum(r2), n = .N), by = .(scale, bin)])
  }))[, .(s = sum(s), n = sum(n)), by = .(scale, bin)]
  Z0 <- Z; Z0[is.na(Z0)] <- 0
  M <- lapply(cols, function(ix) tcrossprod(Z0[, ix, drop = FALSE])); np <- sapply(cols, function(ix) sum(poly[ix]))
  half <- function(tot, parts) (tot - parts) / 2
  unl <- data.table(bin = 99L, s = half(sum(Reduce(`+`, M)^2), sum(sapply(M, function(m) sum(m^2)))),
                    n = half(sum(np)^2, sum(np^2)))
  out <- rbind(w, unl[, scale := "bp"], copy(unl)[, scale := "cM"])
  out[, r2 := s / n]
  list(ld = out[, .(scale, bin, r2, n)],
       info = data.table(n_ind = nrow(G), n_poly = sum(poly), het = mean(G == 1, na.rm = TRUE),
                         n_no_het = sum(rowSums(G == 1, na.rm = TRUE) == 0)))
}
label_bins <- function(d) {
  d[, bin_lab := ifelse(bin == 99L, "unlinked", ifelse(scale == "bp", BP_LABELS[pmin(bin, length(BP_LABELS))],
                                                        CM_LABELS[pmin(bin, length(CM_LABELS))]))]
  d
}

## ---- empirical parents ---------------------------------------------------------------
inp <- readRDS("module_di25/data/di25_inputs.rds")
P <- inp$GTs_par[, u0$unit_marker]
emp_s <- lapply(c(aquilonia = "^Faqu", polyctena = "^Fpol"), function(rx) ld_summary(P[grepl(rx, rownames(P)), , drop = FALSE]))
emp <- rbindlist(lapply(emp_s, `[[`, "ld"), idcol = "species")
emp_info <- rbindlist(lapply(emp_s, `[[`, "info"), idcol = "species")

## ---- simulated parents -----------------------------------------------------------------
one_rep <- function(r) {
  f <- sprintf(BEDFMT, r)
  h <- readLines(f, n = 2); inds <- strsplit(strsplit(h[2], "\t")[[1]][10], "|", fixed = TRUE)[[1]]
  s <- fread(f, skip = 2, header = FALSE, sep = "\t", select = c(1, 3, 10), colClasses = list(character = c(1, 10)), showProgress = FALSE)
  mk <- paste0("Chr", sub("ch", "", s$V1), ":", s$V3); keep <- match(u0$unit_marker, mk)
  par <- grep("^(aq|pol)_", inds)
  S <- do.call(rbind, strsplit(sub("^S", "", s$V10[keep[!is.na(keep)]]), "", fixed = TRUE))[, par, drop = FALSE]
  G <- matrix(NA_integer_, length(par), nrow(u0)); G[, !is.na(keep)] <- t(matrix(suppressWarnings(as.integer(S)), nrow(S)))
  sp <- ifelse(grepl("^aq_", inds[par]), "aquilonia", "polyctena")
  res <- lapply(c("aquilonia", "polyctena"), function(x) ld_summary(G[sp == x, , drop = FALSE]))
  list(ld = rbindlist(lapply(1:2, function(k) res[[k]]$ld[, species := c("aquilonia", "polyctena")[k]]))[, rep := r],
       info = rbindlist(lapply(1:2, function(k) res[[k]]$info[, species := c("aquilonia", "polyctena")[k]]))[, rep := r])
}
sims <- parallel::mclapply(seq_len(N_REPS), one_rep, mc.cores = N_CORES)
sim_ld <- rbindlist(lapply(sims, `[[`, "ld")); sim_info <- rbindlist(lapply(sims, `[[`, "info"))

cat("\n[parents] sample summary (empirical; simulated = mean over replicates):\n")
print(rbind(emp_info[, data := "empirical"],
            sim_info[, .(n_ind = mean(n_ind), n_poly = mean(n_poly), het = mean(het), n_no_het = mean(n_no_het)), by = species][, data := "simulated"]),
      digits = 3)

emp <- label_bins(emp)[, data := "empirical"]
sim <- label_bins(sim_ld)[, .(r2_mean = mean(r2), r2_lo = quantile(r2, 0.025), r2_hi = quantile(r2, 0.975)), by = .(species, scale, bin, bin_lab)]
cmp <- merge(emp[, .(species, scale, bin, bin_lab, r2_emp = r2, n_emp = n)], sim, by = c("species", "scale", "bin", "bin_lab"))
setorder(cmp, species, scale, bin)
cat("\n[parents] mean within-species r^2 (empirical vs simulated mean [95% range]):\n")
print(cmp[, .(species, scale, bin_lab, r2_emp = round(r2_emp, 3), sim = sprintf("%.3f [%.3f, %.3f]", r2_mean, r2_lo, r2_hi))])
saveRDS(list(compare = cmp, emp_info = emp_info, sim_info = sim_info, sim_ld = sim_ld, n_reps = N_REPS),
        file.path(OUT_DATA, "explore_parent_ld.rds"))

d <- rbind(cmp[, .(species, scale, bin, bin_lab, data = "empirical parents", r2 = r2_emp, lo = NA_real_, hi = NA_real_)],
           cmp[, .(species, scale, bin, bin_lab, data = "simulated parents (buggy founders)", r2 = r2_mean, lo = r2_lo, hi = r2_hi)])
d[, bin_lab := factor(bin_lab, levels = unique(c(BP_LABELS, CM_LABELS, "unlinked")))]
p <- ggplot(d, aes(bin_lab, r2, colour = data, group = data)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = data), colour = NA, alpha = 0.2, show.legend = FALSE) +
  geom_line() + geom_point(size = 1.2) +
  facet_grid(species ~ scale, scales = "free_x") + scale_y_log10() +
  scale_colour_manual(values = c("black", "#d95f02"), aesthetics = c("colour", "fill"), name = NULL) +
  labs(x = "distance between unit SNPs", y = expression("mean within-species "*r^2*" (log)"),
       caption = sprintf("n = 15 per species: unlinked E[r^2] ~ 0.07. Simulated: mean and 95%% range over %d replicates.", N_REPS)) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
ggsave(file.path(OUT_FIG, "explore_parent_ld.png"), p, width = 11, height = 6.5, dpi = 200)
cat("[parents] done\n")
