## =========================================================
## module_manuscript_rho05 -- which populations drive the PC1-ancestry link?
## =========================================================
## Follows up moduleB_ancestry_vs_winter_climate.R's finding that PC1 (not
## PC2) has the more defensible ancestry correlation (aland_excluded r=-0.502,
## block-permutation p=0.033 across 17 independent lineage units). Identifies
## which population(s)/unit(s) are driving it via:
##   1. leave-one-population-out delta r (19 populations)
##   2. leave-one-INDEPENDENT-UNIT-out delta r (17 units -- LangholmenW+R and
##      Bunkkeri+Grundsund dropped together, matching the block-permutation's
##      own unit definition)
##   3. Cook's distance / hat leverage from lm(ancestry ~ PC1)
##
## Run from the repo root: Rscript module_manuscript_rho05/R/moduleB_ancestry_PC1_leverage.R
## =========================================================

suppressMessages({ library(data.table); library(ggplot2) })

obj <- readRDS("module_manuscript_rho05/data/moduleB_ancestry_vs_winter_climate.rds")
dt19 <- obj$dt19; unit <- obj$units
## rebuild the population->unit map exactly as before
pop_unit <- dt19$Population
pop_unit[dt19$Population %in% c("LangholmenW", "LangholmenR")] <- "unit_Lang"
pop_unit[dt19$Population %in% c("Bunkkeri", "Grundsund")]      <- "unit_BunGru"

obs_r <- cor(dt19$PC1, dt19$ancestry)
cat("observed cor(PC1, ancestry), aland_excluded n=19:", round(obs_r, 3), "\n\n")

## ---- 1. leave-one-population-out --------------------------------------
loo_pop <- rbindlist(lapply(dt19$Population, function(p) {
  d <- dt19[Population != p]
  data.table(Population = p, n = nrow(d), r_without = round(cor(d$PC1, d$ancestry), 3),
             delta = round(cor(d$PC1, d$ancestry) - obs_r, 3))
}))
setorder(loo_pop, delta)
cat("=== leave-one-POPULATION-out (sorted by delta = r_without - r_observed) ===\n")
cat("(most negative delta = removing it makes r MORE negative, i.e. that population was DAMPENING the signal;\n")
cat(" most positive delta = removing it makes r LESS negative, i.e. that population was DRIVING the signal)\n")
print(loo_pop)

## ---- 2. leave-one-UNIT-out (respects the 2 shared-origin pairs) --------
units <- unique(pop_unit)
loo_unit <- rbindlist(lapply(units, function(u) {
  d <- dt19[pop_unit != u]
  pops_in_unit <- dt19$Population[pop_unit == u]
  data.table(unit = u, populations = paste(pops_in_unit, collapse = "+"), n = nrow(d),
             r_without = round(cor(d$PC1, d$ancestry), 3),
             delta = round(cor(d$PC1, d$ancestry) - obs_r, 3))
}))
setorder(loo_unit, delta)
cat("\n=== leave-one-INDEPENDENT-UNIT-out (", length(units), " units) ===\n", sep = "")
print(loo_unit)

## ---- 3. Cook's distance / leverage from lm(ancestry ~ PC1) -------------
m <- lm(ancestry ~ PC1, data = dt19)
infl <- data.table(Population = dt19$Population, ancestry = dt19$ancestry, PC1 = dt19$PC1,
                   hat = round(hatvalues(m), 3), cooksD = round(cooks.distance(m), 3),
                   std_resid = round(rstandard(m), 2))
setorder(infl, -cooksD)
cat("\n=== regression influence diagnostics, lm(ancestry ~ PC1) ===\n")
cat("(hat > ", round(2 * 2 / nrow(dt19), 3), " = high leverage;  Cook's D > ", round(4 / nrow(dt19), 3), " = high influence, rule-of-thumb cutoffs)\n", sep = "")
print(infl)

saveRDS(list(loo_pop = loo_pop, loo_unit = loo_unit, infl = infl, obs_r = obs_r),
        "module_manuscript_rho05/data/moduleB_ancestry_PC1_leverage.rds")

## ---- figure: scatter, points sized by Cook's D, top driver(s) labelled --
dt19b <- copy(dt19)
dt19b[infl, on = "Population", `:=`(cooksD = i.cooksD)]
top_units <- loo_unit[order(-abs(delta))][1:2, unit]
dt19b[, unit := pop_unit]
dt19b[, is_top := unit %in% top_units]
p <- ggplot(dt19b, aes(PC1, ancestry)) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey70", linewidth = 0.5) +
  geom_point(aes(size = cooksD, colour = is_top)) +
  geom_text(aes(label = Population), size = 2.3, vjust = -0.8) +
  scale_colour_manual(values = c("FALSE" = "grey40", "TRUE" = "#D55E00"), guide = "none") +
  scale_size_continuous(name = "Cook's D") +
  labs(x = "climate PC1", y = "genome-wide aquilonia ancestry",
       title = sprintf("PC1 vs ancestry (aland_excluded, n=19): r=%.3f", obs_r),
       subtitle = sprintf("orange = the 1-2 independent units with the largest leave-one-out effect on r (%s)",
                          paste(top_units, collapse = ", "))) +
  theme_classic(base_size = 10)
ggsave("module_manuscript_rho05/Figures/moduleB_ancestry_PC1_leverage.png", p, width = 7.5, height = 5.5, dpi = 200)
cat("\nSaved figure + moduleB_ancestry_PC1_leverage.rds\n")
