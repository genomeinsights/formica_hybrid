## =========================================================================
## sim_founder_fix -- acceptance test: do simulated parents / founders reproduce the
## LD patterns of the EMPIRICAL parents?
##
## Usage (from the formica_hybrid repo root):
##   Rscript sim_founder_fix/check_parent_ld.R bed  <diem_boot_output.bed> [more .bed ...]
##   Rscript sim_founder_fix/check_parent_ld.R vcf  <founder_vcf_dir> [N_SETS=10]
## 'bed': the 15 + 15 aq_/pol_ parents of each DIEM output replicate (the existing
##        pipeline). 'vcf': random sets of 15 + 15 individuals drawn from a founder
##        VCF pool (e.g. make_mosaic_founders.R output); n = 15 matches the empirical
##        parents, so all data sets share the same finite-sample floor.
## Prints, per species, mean r^2 by class for the empirical parents and the tested
## source (mean and range over replicates / sets), sample summaries, and four checks:
##   1. no LD class is (near-)perfect: mean r^2 < 0.9 in every class
##   2. close pairs in DIFFERENT simulation clusters keep >= 75% of the empirical r^2
##   3. close pairs in the SAME simulation cluster are within +-0.10 of empirical
##   4. the 0.2-1 cM and 1-5 cM classes are within +-0.03 of empirical
## (thresholds are pragmatic; with only 15 donor parents, mosaic founders overshoot
##  short-range r^2 by ~0.05). Writes <OUT_PREFIX>.rds and <OUT_PREFIX>.png
##  (OUT_PREFIX defaults to sim_founder_fix/out/check_parent_ld_<mode>).
## =========================================================================
source("sim_founder_fix/parent_ld_lib.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: check_parent_ld.R bed <files...> | vcf <dir> [N_SETS]")
MODE <- args[1]
OUT_PREFIX <- Sys.getenv("OUT_PREFIX", sprintf("sim_founder_fix/out/check_parent_ld_%s", MODE))
dir.create(dirname(OUT_PREFIX), showWarnings = FALSE, recursive = TRUE)

panel <- load_panel(); units <- panel$units; PS <- make_pairs(units)
emp_G <- read_empirical_parents(panel$parents, units)
emp <- rbindlist(lapply(names(emp_G), function(sp) ld_profile(emp_G[[sp]], PS)[, species := sp]))
emp_info <- rbindlist(lapply(names(emp_G), function(sp) sample_info(emp_G[[sp]])[, species := sp]))

sources <- if (MODE == "bed") {
  lapply(args[-1], function(f) read_diem_parents(f, units))
} else if (MODE == "vcf") {
  n_sets <- if (length(args) >= 3) as.integer(args[3]) else 10L
  lapply(seq_len(n_sets), function(s) read_founder_vcfs(args[2], units, 15L, seed = s))
} else stop("MODE must be 'bed' or 'vcf'")

tst <- rbindlist(lapply(seq_along(sources), function(k) rbindlist(lapply(names(sources[[k]]), function(sp)
  ld_profile(sources[[k]][[sp]], PS)[, `:=`(species = sp, rep = k)]))))
tst_info <- rbindlist(lapply(seq_along(sources), function(k) rbindlist(lapply(names(sources[[k]]), function(sp)
  sample_info(sources[[k]][[sp]])[, `:=`(species = sp, rep = k)]))))[, lapply(.SD, mean), by = species, .SDcols = !"rep"]

cmp <- merge(emp[, .(species, cls, r2_emp = r2)],
             tst[, .(r2_test = mean(r2), lo = min(r2), hi = max(r2)), by = .(species, cls)], by = c("species", "cls"))
setorder(cmp, species, cls)
cat("\n[check] samples (test: mean over replicates/sets):\n")
print(rbind(emp_info[, source := "empirical"], tst_info[, source := MODE]), digits = 3)
cat("\n[check] mean within-species r^2:\n"); print(cmp, digits = 3)

g <- function(sp, cl, col) cmp[species == sp & cls == cl][[col]]
checks <- rbindlist(lapply(c("aquilonia", "polyctena"), function(sp) data.table(
  species = sp,
  check = c("1 no near-perfect LD class", "2 <5kb different cluster >= 75% of empirical",
            "3 <5kb same cluster within +-0.10", "4 0.2-1 cM within +-0.03", "4 1-5 cM within +-0.03"),
  pass = c(all(cmp[species == sp, r2_test] < 0.9),
           g(sp, "<5kb different sim cluster", "r2_test") >= 0.75 * g(sp, "<5kb different sim cluster", "r2_emp"),
           abs(g(sp, "<5kb same sim cluster", "r2_test") - g(sp, "<5kb same sim cluster", "r2_emp")) <= 0.10,
           abs(g(sp, "0.2-1", "r2_test") - g(sp, "0.2-1", "r2_emp")) <= 0.03,
           abs(g(sp, "1-5", "r2_test") - g(sp, "1-5", "r2_emp")) <= 0.03))))
cat("\n[check] acceptance:\n"); print(checks)
cat(sprintf("\n[check] OVERALL: %s\n", if (all(checks$pass)) "PASS" else "FAIL"))
saveRDS(list(compare = cmp, checks = checks, emp_info = emp_info, test_info = tst_info, per_rep = tst, mode = MODE, inputs = args[-1]),
        paste0(OUT_PREFIX, ".rds"))

d <- rbind(cmp[, .(species, cls, source = "empirical parents", r2 = r2_emp, lo = NA_real_, hi = NA_real_)],
           cmp[, .(species, cls, source = paste("tested:", MODE), r2 = r2_test, lo, hi)])
p <- ggplot(d, aes(cls, r2, colour = source, group = source)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) + geom_point(size = 1.6) +
  geom_line(data = d[!grepl("sim cluster", cls)]) +
  facet_wrap(~ species) + scale_y_log10() + scale_colour_manual(values = c("black", "#d95f02"), name = NULL) +
  labs(x = NULL, y = expression("mean within-species "*r^2*" (log), n = 15")) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom")
ggsave(paste0(OUT_PREFIX, ".png"), p, width = 10, height = 4.8, dpi = 200)
