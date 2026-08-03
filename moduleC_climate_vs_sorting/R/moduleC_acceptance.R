## =========================================================
## MODULE C (revised) -- pre-pilot acceptance checks (no BayPass needed)
## =========================================================
## Exercises the acceptance checklist that does not require a BayPass batch:
##   1  all scripts parse
##   2  annotation table = 32,840 unique eMLGs in canonical BayPass order
##   3  prop_fixed / DI / recomb / differentiated / directional correctly aligned
##   4  observed BF vectors reproduce the association RDS eBF1/eBF2
##   5  batch COVARIABLE/MRK ordering assertion accepts correct, rejects scrambled input
##   7  threshold labels resolve to exactly top 0.1% / 0.5% / 1% and all metrics
##   8  the primary sorting statistic uses differentiated eMLGs only
##   9  the continuous sorting-magnitude statistic uses prop_fixed, not uni_score
##   10 a deliberately altered reproduction check makes moduleC_analyse.R STOP
##   11 checkpoint fingerprint mismatch makes moduleC_null_regen.R STOP on resume
## Checks 6 (all pilot statistics finite) and the live batch ordering (5, on real
## output) are confirmed on the pilot's own log after moduleC_null_regen.R runs.
##
## Run from the formica_hybrid repo root. Fast (seconds). Prints PASS/FAIL per check.

suppressMessages({ library(data.table); library(digest) })
BASE <- "moduleC_climate_vs_sorting"
ok <- function(n, cond, msg = "") cat(sprintf("  [%s] %-2s %s %s\n",
        ifelse(cond, "PASS", "FAIL"), n, msg, ""))
FAILS <- 0L; check <- function(n, cond, msg) { if (!isTRUE(cond)) FAILS <<- FAILS + 1L; ok(n, isTRUE(cond), msg) }

source(file.path(BASE, "R", "moduleC_stat_functions.R"))
ann <- readRDS(file.path(BASE, "data", "moduleC_annotations.rds"))
grp <- readLines("aland_excluded_eMLG/eMLG_group_order.txt")
A   <- prepare_annotation_ranks(ann); NM <- A$N

## SAFETY: this harness writes throwaway null_stats/checkpoint files for the tamper
## tests (10, 11). Refuse to run if real ones exist, so it can never clobber a
## pilot/full run's outputs or a live checkpoint.
for (f in c("data/moduleC_null_stats.rds", "data/moduleC_null_ckpt.rds"))
  if (file.exists(file.path(BASE, f)))
    stop("refusing to run acceptance harness: real ", f, " present (would be clobbered by tamper tests)")

cat("\n=== Module C acceptance checks (pre-pilot) ===\n")

## 1 -- all scripts parse
scripts <- list.files(file.path(BASE, "R"), pattern = "\\.R$", full.names = TRUE)
parse_ok <- all(vapply(scripts, function(f)
  tryCatch({ parse(f); TRUE }, error = function(e) FALSE), logical(1)))
check(1, parse_ok, sprintf("all %d scripts parse", length(scripts)))

## 2 -- annotation table: 32,840 unique eMLGs in canonical BayPass order
check(2, nrow(ann) == 32840 && !any(duplicated(ann$group_id)) && identical(ann$group_id, grp),
      "annotation = 32,840 unique eMLGs, canonical order")

## 3 -- alignment / finiteness of annotations
pf_na <- is.na(ann$prop_fixed)
check(3, all(is.finite(ann$DI)) && all(is.finite(ann$recomb)) &&
         all(ann$directional %in% c(0L, 1L)) && all(!is.na(ann$differentiated)) &&
         identical(pf_na, is.na(ann$uni_score)) &&
         all(ann$prop_fixed[!pf_na] >= 0 & ann$prop_fixed[!pf_na] <= 1),
      "DI/recomb/differentiated/directional aligned; prop_fixed in [0,1] w/ established NAs")

## 4 -- observed BF reproduce association eBF1/eBF2
b1 <- fread("aland_excluded_eMLG/PC1_eMLG_withOmega_summary_betai_reg.out")$`BF(dB)`
b2 <- fread("aland_excluded_eMLG/PC2_eMLG_withOmega_summary_betai_reg.out")$`BF(dB)`
check(4, max(abs(b1 - ann$eBF1)) <= 1e-6 && max(abs(b2 - ann$eBF2)) <= 1e-6,
      sprintf("observed BF == eBF (max|d| = %.2e)", max(max(abs(b1 - ann$eBF1)), max(abs(b2 - ann$eBF2)))))

## 5 -- batch ordering assertion (overflow-tolerant): accepts covariate-major with a
##      '***'-overflowed COVARIABLE at cov >= 1000; rejects a scrambled MRK order.
##      Mirrors the logic in moduleC_null_regen.R on a small synthetic (NMs x Bs).
NMs <- 4L; Bs <- 1000L                                    # Bs>=1000 to exercise overflow
row_cov <- rep(seq_len(Bs), each = NMs); row_mrk <- rep(seq_len(NMs), times = Bs)
covchr <- as.character(row_cov); covchr[row_cov >= 1000] <- "***"   # simulate field overflow
ord_ok <- function(mrk, covc) {
  cn <- suppressWarnings(as.integer(covc)); ok <- !is.na(cn)
  all(mrk == rep(seq_len(NMs), times = Bs)) &&
    all(cn[ok] == row_cov[ok]) && all(row_cov[!ok] >= 1000)
}
accepts <- ord_ok(row_mrk, covchr)                        # correct + overflow -> accept
bad_mrk <- row_mrk; bad_mrk[1:2] <- bad_mrk[2:1]          # swap markers -> reject
rejects <- !ord_ok(bad_mrk, covchr)
check(5, accepts && rejects,
      "ordering check accepts covariate-major (with '***' overflow at cov>=1000), rejects scrambled MRK")

## 7 -- threshold labels resolve to exactly the three fractions + all metrics
frac_lookup <- c(top0001 = "top 0.1%", top0005 = "top 0.5%", top0010 = "top 1%")
metric_lookup <- c(DI_med = "median DI", rec_med = "median recomb",
                   pdiff = "prop. differentiated", pdir_diff = "prop. directional | diff.")
th_stats <- grep("^top", covariate_stat_names(), value = TRUE)
tags  <- sub("_.*", "", th_stats); mets <- sub("^top[0-9]+_", "", th_stats)
check(7, setequal(frac_lookup[unique(tags)], unname(frac_lookup)) &&
         all(!is.na(frac_lookup[tags])) && all(!is.na(metric_lookup[mets])) &&
         setequal(unname(metric_lookup[unique(mets)]), unname(metric_lookup)),
      "threshold tags -> {0.1%,0.5%,1%} x {DI,recomb,pdiff,pdir_diff}, none NA")

## 8 -- primary sorting statistic uses differentiated eMLGs only
##      (recompute the contrast manually on a random BF and match the function)
check(8, identical(PRIMARY_STATS[["sorting"]], "sort_gap_differentiated"), "PRIMARY sorting = sort_gap_differentiated")
set.seed(1); bf <- rnorm(NM)
s <- compute_covariate_stats(bf, A)
pct <- (rank(bf) - 0.5) / NM; dd <- ann$differentiated == TRUE; dir <- ann$directional == 1L
manual_diff <- mean(pct[dd & dir]) - mean(pct[dd & !dir])
manual_all  <- mean(pct[dir]) - mean(pct[!dir])
check("8b", abs(s[["sort_gap_differentiated"]] - manual_diff) < 1e-12 &&
            abs(s[["sort_gap_all"]] - manual_all) < 1e-12 &&
            manual_diff != manual_all,
      "sort_gap_differentiated matches differentiated-only contrast (!= all-eMLG)")

## 9 -- sorting magnitude uses prop_fixed (not uni_score)
ok9 <- !is.na(ann$prop_fixed)
manual_mag <- cor(rank(bf[ok9]), rank(ann$prop_fixed[ok9]))
manual_ori <- cor(rank(bf[!is.na(ann$uni_score)]), rank(ann$uni_score[!is.na(ann$uni_score)]))
check(9, abs(s[["rho_sort_magnitude"]] - manual_mag) < 1e-12 &&
         abs(s[["rho_sort_orientation"]] - manual_ori) < 1e-12,
      "rho_sort_magnitude = Spearman(BF, prop_fixed); orientation = Spearman(BF, uni_score)")

## 10 -- an altered reproduction check makes analyse.R STOP (no trusted output)
tmp <- tempfile(fileext = ".rds")
fake_null <- matrix(rnorm(200 * length(covariate_stat_names())), nrow = 200,
                    dimnames = list(NULL, covariate_stat_names()))
fake_obs  <- rbind(PC1 = colMeans(fake_null), PC2 = colMeans(fake_null))
saveRDS(list(observed = fake_obs, null_stats = fake_null, stat_names = colnames(fake_null),
             k_check = list(reproduced = FALSE, cor_k1 = 0.80, cor_k2 = 0.995,
                            rel1 = 0.11, rel2 = 0.00, max_abs_dk1 = 9000L, max_abs_dk2 = 40L,
                            k1r = integer(NM)),   # r below 0.99 & sum ratio off -> must STOP
             params = list(), fingerprint = list(), session = sessionInfo()),
        file.path(BASE, "data", "moduleC_null_stats.rds"))
rc <- system2("Rscript", file.path(BASE, "R", "moduleC_analyse.R"),
              stdout = FALSE, stderr = FALSE)
check(10, rc != 0, "analyse.R STOPS when k_check$reproduced == FALSE (nonzero exit)")
invisible(file.remove(file.path(BASE, "data", "moduleC_null_stats.rds")))   # remove the fake

## 11 -- checkpoint fingerprint mismatch makes regen STOP on resume (before BayPass)
bogus <- list(null_stats = matrix(0, 10000, length(covariate_stat_names()),
                                  dimnames = list(NULL, covariate_stat_names())),
              k1r = integer(NM), k2r = integer(NM), done = 1L,
              fingerprint = list(bogus = "does-not-match"), group_id = grp)
CKPT <- file.path(BASE, "data", "moduleC_null_ckpt.rds")
stopifnot("a real checkpoint exists -- refusing to overwrite" = !file.exists(CKPT))
saveRDS(bogus, CKPT)
rc11 <- system2("Rscript", file.path(BASE, "R", "moduleC_null_regen.R"),
                env = "MODC_NBATCH=1", stdout = FALSE, stderr = FALSE)
check(11, rc11 != 0, "null_regen.R STOPS on fingerprint mismatch (before running BayPass)")
if (file.exists(CKPT)) invisible(file.remove(CKPT))          # remove the bogus checkpoint

cat(sprintf("\n=== acceptance: %s (%d failures) ===\n", ifelse(FAILS == 0, "ALL PASS", "FAILURES"), FAILS))
quit(status = if (FAILS == 0) 0 else 1)
