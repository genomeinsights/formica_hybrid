## =========================================================
## Admixture-geometry exposure diagnostic (Module A)
## =========================================================
## Could a unit clear the near-fixation floor phi in a hybrid population through
## admixture geometry ALONE -- q_tilde = h*p_aqu + (1-h)*p_pol -- with no sorting,
## and if so is that exposure structured by the diagnostic index (DI)? A DI-structured
## artefact would compete with the reported DI-dependent sorting DIRECTION.
## DIAGNOSTIC ONLY: no null model, no drift correction, no significance test.
##
## Inputs (verified, names printed before use in the exploration step):
##   moduleA_sorting/data/moduleA_prep_snp.rds  -- prep$pop_means (22 pops x markers, mean dosage)
##   data/hybrids_and_parents_maf005.Rdata      -- map_hyb_005 (DI), sample_data (pop sizes)
## Orientation mirrors parallelism_stats(): higher oriented freq = aquilonia allele.
## No pre-computed hybrid index exists in the data, so h is DERIVED from the same
## oriented frequencies (orientation therefore consistent by construction).
## Run from the repo root.

suppressMessages(library(data.table))
PHI_GRID <- c(0.80, 0.85, 0.90, 0.95)
MIN_PARENT_MAF <- 0.15
AQU <- "aquilonia_parent"; POL <- "polyctena_parent"

prep <- readRDS("moduleA_sorting/data/moduleA_prep_snp.rds")
P    <- prep$pop_means / 2                                   # pops x markers, allele freq (unoriented)
mk   <- colnames(P)
e <- new.env(); load("data/hybrids_and_parents_maf005.Rdata", envir = e)
map <- as.data.table(e$map_hyb_005); smeta <- as.data.table(e$sample_data_with_parents); rm(e)
DI  <- setNames(map$DiagnosticIndex, map$marker)[mk]
nk  <- smeta[!grepl("_parent$", Population), .(n_k = .N), by = Population]     # diploids / hybrid pop
hyb <- setdiff(rownames(P), c(AQU, POL))

## ---- orientation + differentiated (retained) set ------------------------
f_aqu_par <- P[AQU, ]; f_pol_par <- P[POL, ]
sgn  <- sign(f_aqu_par - f_pol_par)                          # +1 => dosage allele is aquilonia
pbar <- (f_aqu_par + f_pol_par) / 2                          # pooled parental freq (parents 15 vs 15)
pmaf <- pmin(pbar, 1 - pbar)
keep <- which(is.finite(sgn) & sgn != 0 & pmaf >= MIN_PARENT_MAF)
cat(sprintf("markers total %d | differentiated (oriented, parent_maf>=%.2f) %d\n",
            length(mk), MIN_PARENT_MAF, length(keep)))
stopifnot("retained set violates parent_maf>=0.15" = all(pmaf[keep] >= MIN_PARENT_MAF))   # sanity 4

p_aqu <- ifelse(sgn[keep] > 0, f_aqu_par[keep], 1 - f_aqu_par[keep])   # oriented parental aqu-allele freq
p_pol <- ifelse(sgn[keep] > 0, f_pol_par[keep], 1 - f_pol_par[keep])
DIk   <- DI[keep]
qobs  <- t(P[hyb, keep, drop = FALSE])                       # markers x hybrid pops, oriented
flip  <- sgn[keep] < 0; qobs[flip, ] <- 1 - qobs[flip, ]
cat(sprintf("q_obs cells with data: %.4f | NA (all-missing pop x marker): %d\n",
            mean(!is.na(qobs)), sum(is.na(qobs))))

## ---- genome-wide hybrid index h per pop ---------------------------------
## least-squares admixture proportion: q_obs ~ p_pol + h*(p_aqu - p_pol).
dif    <- p_aqu - p_pol
h_reg  <- colSums(dif * (qobs - p_pol), na.rm = TRUE) / sum(dif^2)
h_mean <- colMeans(qobs, na.rm = TRUE)                       # cross-check (mean oriented aqu freq)
h      <- pmin(pmax(h_reg, 0), 1)                            # clamp to [0,1] for the geometry model
cat("\n[Step 1] hybrid index h (regression, oriented toward aquilonia), sorted:\n")
print(round(sort(setNames(h, hyb)), 3))
cat(sprintf("h range %.3f..%.3f | |h-0.5| range %.3f..%.3f | max|h_reg - h_mean| %.3f\n",
            min(h), max(h), min(abs(h - 0.5)), max(abs(h - 0.5)), max(abs(h_reg - h_mean))))

## ---- geometry expectation + observed near-fixation ----------------------
qtil <- outer(p_aqu, h) + outer(p_pol, 1 - h)               # markers x pops (phi-independent)
stopifnot("q_tilde outside [0,1]" = all(qtil >= 0 & qtil <= 1))              # sanity 3
DI_BREAKS <- quantile(DIk, 0:10 / 10, na.rm = TRUE)
dec <- cut(DIk, DI_BREAKS, include.lowest = TRUE, labels = FALSE)            # DI decile 1(low)..10(high)

exposure_by_phi <- function(phi) {
  nf_a <- qobs >= phi;      nf_p <- qobs <= 1 - phi          # OBSERVED near-fixed calls (the signal)
  ex_a <- qtil >= phi;      ex_p <- qtil <= 1 - phi          # geometry-exposed cells
  callA <- rowSums(nf_a, na.rm = TRUE); callP <- rowSums(nf_p, na.rm = TRUE)
  expA  <- rowSums(nf_a & ex_a, na.rm = TRUE); expP <- rowSums(nf_p & ex_p, na.rm = TRUE)
  list(phi = phi,
       any_exposed_unit = mean(rowSums(ex_a | ex_p) > 0),                    # context (Step 3.1)
       by_dec = data.table(dec, callA, callP, expA, expP)[
         , .(n_call = sum(callA + callP),
             frac_exp     = sum(expA + expP) / max(1, sum(callA + callP)),   # KEY column (Step 3.2)
             frac_exp_aqu = sum(expA) / max(1, sum(callA)),
             frac_exp_pol = sum(expP) / max(1, sum(callP))), by = dec][order(dec)],
       per_pop = data.table(pop = hyb, exp_cells = colSums((nf_a & ex_a) | (nf_p & ex_p), na.rm = TRUE)))
}
res <- lapply(PHI_GRID, exposure_by_phi); names(res) <- sprintf("phi%.2f", PHI_GRID)

## ---- Step 3 tables ------------------------------------------------------
tab_dec <- rbindlist(lapply(res, function(r) r$by_dec[, phi := r$phi][]))
cat("\n[Step 3] EXPOSURE by DI decile x phi (frac_exp = fraction of observed near-fixed calls",
    "explainable by admixture geometry):\n")
print(dcast(tab_dec, dec ~ phi, value.var = "frac_exp"))
cat("\noverall fraction of near-fixed calls exposed, and 'any-exposed-unit' context, per phi:\n")
print(data.table(phi = PHI_GRID,
                 frac_calls_exposed = sapply(res, function(r) sum(r$by_dec$frac_exp * r$by_dec$n_call) / sum(r$by_dec$n_call)),
                 any_exposed_unit   = sapply(res, function(r) r$any_exposed_unit)))
cat("\ndirectional breakdown at primary phi=0.85 (toward each parent, by DI decile):\n")
print(res[["phi0.85"]]$by_dec[, .(dec, n_call, frac_exp_aqu = round(frac_exp_aqu, 4), frac_exp_pol = round(frac_exp_pol, 4))])

## ---- Step 4 sanity: exposure concentrates in skewed-h pops; h=0.5 test ---
pp <- res[["phi0.85"]]$per_pop[, `:=`(h = h, absdev = abs(h - 0.5))][order(-exp_cells)]
cat("\n[Step 4.2] per-pop exposed cells vs |h-0.5| (Spearman rho =",
    round(cor(pp$exp_cells, pp$absdev, method = "spearman"), 3), "):\n"); print(head(pp, 6))
for (phi in PHI_GRID) {                                                     # sanity 1: force h=0.5
  qt0 <- outer(p_aqu, rep(0.5, length(hyb))) + outer(p_pol, rep(0.5, length(hyb)))
  cat(sprintf("[Step 4.1] h=0.5, phi=%.2f -> exposed cells = %d (expect 0 for phi>0.85; only p_bar boundary at 0.85)\n",
              phi, sum(qt0 >= phi | qt0 <= 1 - phi)))
}

## ---- Step 5 discreteness sub-check --------------------------------------
disc <- nk[data.table(Population = hyb), on = "Population"]
disc[, `:=`(gametes = 2 * n_k,
            min_count_085 = ceiling(0.85 * 2 * n_k),
            eff_bar_085 = ceiling(0.85 * 2 * n_k) / (2 * n_k))]
obs_rate <- data.table(pop = hyb,
                       nf_rate = colMeans((qobs >= 0.85) | (qobs <= 0.15), na.rm = TRUE))
disc <- obs_rate[disc, on = c(pop = "Population")]
cat("\n[Step 5] discreteness (phi=0.85): smallest achievable near-fixation count and observed rate:\n")
print(disc[order(n_k), .(pop, n_k, gametes, min_count_085, eff_bar_085 = round(eff_bar_085, 3),
                         nf_rate = round(nf_rate, 4))])

## do small (distorted-bar) pops concentrate in particular DI deciles?
small  <- hyb %in% disc[n_k <= 5, pop]
nf85   <- (qobs >= 0.85) | (qobs <= 0.15)
small_share <- data.table(dec, n_small = rowSums(nf85[, small, drop = FALSE], na.rm = TRUE),
                          n_all = rowSums(nf85, na.rm = TRUE))[
  , .(frac_calls_from_small_pops = round(sum(n_small) / max(1, sum(n_all)), 3)), by = dec][order(dec)]
cat("\n[Step 5] fraction of near-fixed calls (phi=0.85) contributed by small pops (n_k<=5), by DI decile:\n")
print(small_share)

saveRDS(list(tab_dec = tab_dec, per_pop = pp, discreteness = disc, small_share = small_share,
             h = setNames(h, hyb), h_mean = setNames(h_mean, hyb), n_diff = length(keep)),
        "moduleA_sorting/data/moduleA_admixture_exposure.rds")
cat("\nSaved: moduleA_sorting/data/moduleA_admixture_exposure.rds\n")
