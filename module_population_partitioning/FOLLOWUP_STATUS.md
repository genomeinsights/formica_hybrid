# Follow-up status: alternative explanations for heterogeneous, locus-specific sorting

Started 2026-09-13, per the follow-up brief requesting Analyses 1–5 (geographic
prediction, individual influence, ancestry-run proxy, parental structure by
DI, synthesis). This file is updated after each completed analysis. Nothing
in `module_population_partitioning/R/pp_*.R`, `data/`, `Figures/`, `README.md`,
`AUDIT.md`, `CROSS_MODULE_INPUTS.md` or `doc/` was modified for this
follow-up work; all new code is `R/1X_*.R`, all new outputs are under
`data/followup/` and `Figures/followup/`.

## Cross-cutting notes (apply to all analyses below)

### Dataset version discrepancy — resolved in favour of the current primary dataset

The follow-up brief names `module_di25/data/di25_clustering_cM5.rds` /
`di25_sorting_emlg.rds` (11,052 units, fixed `min_r2=0.2`) as an authoritative
input to "retain". That lineage was **superseded on 2026-09-13**, earlier the
same day, by an explicit user decision to migrate this module's primary
dataset to `module_di25_rho05` (20,807 units, `min_r2_rho=0.5`) — documented
in `README.md` "Migration to the rho05 primary dataset" and
`CROSS_MODULE_INPUTS.md`, both of which the follow-up brief itself lists as
required reading. All follow-up analyses therefore use the **current**
dataset (`pp_residual_ancestry.rds`, 20,807 units) rather than the legacy
11,052-unit numbers the brief also names. Flagged in each script's header;
re-run against the legacy dataset on request if that is genuinely wanted for
some specific reason.

### `pp_units_final.rds` is stale — not used

Inspection (required by the brief) found `module_population_partitioning/data/pp_units_final.rds`
still has 11,052 rows (mtime 2026-09-10), because the `saveRDS()` call that
wrote it was dropped from `pp_figures.R` during the AUDIT.md revision and
never re-added. This is a real inconsistency in the "completed" pipeline,
but per "do not modify or overwrite the completed population-partitioning
analyses" it was left alone and flagged rather than fixed. Follow-up work
uses `pp_residual_ancestry.rds$u` directly instead (current, 20,807 rows,
verified to be a superset of `pp_units_final.rds`'s columns).

### Correction (2026-09-13): Sielva's coordinate is genuine, not a data error

Analysis 1 originally flagged Sielva's coordinate (46.61°N, 10.44°E — the
Italian Alps) as a likely transcription error, reasoning from its bioclim
values (bio1/bio6/bio11) fitting smoothly into the Fennoscandian
latitude-ordered climate gradient of the other 19 populations. **The user
confirmed Sielva is a genuine Alpine site** — the climate-value reasoning
was the wrong kind of evidence to override an actual place name (and, in
hindsight, is easily explained: Alpine and Fennoscandian sites at comparable
elevation/exposure can have similar winter minima, which is exactly the
variable that lined up).

This is more consequential than a data error would have been: it means one
of the 20 "hybrid populations" this entire module has treated as part of a
single, roughly-contiguous Fennoscandian hybrid zone is geographically
disjunct by ~3000km from the other 19 — plausibly a different hybrid zone
and colonization/admixture history, not just a distant outlier population
within the same zone. This reframes earlier "Sielva is exceptional" findings
in this module (disproportionate row-centred PCA loading on the shared
partition axis in `pp_pca_refined.R`, F1-like elevated heterozygosity noted
throughout `module_di25`) as plausibly reflecting genuine geographic/
evolutionary distinctness rather than a coordinate artefact. Analysis 1's
own numerical results are unaffected (the leave-Sielva-out row was already
not an outlier among the 20 leave-one-out folds), but the earlier framing
("likely wrong", "recommend verifying with the field team") is withdrawn.
Sielva should be read throughout this document as a genuine, geographically
distant reference population, not a data-quality concern.

### Nyrhispera naming

`data/bioclimatic_variables.csv` uses `Nyrhispera1`/`Nyrhispera2`; this
module's population set uses `Nyrhispera74`/`Nyrhispera75`. Crosswalk
`Nyrhispera1→74`, `Nyrhispera2→75` is the **established** convention already
used identically in `module_manuscript_rho05/R/moduleB_stage1_prepare_bioclim_covariates.R`
and `moduleB_ancestry_vs_winter_climate.R` (grep-verified) — not a new guess.

---

## Analysis 1: geographic prediction of residual ancestry profiles — COMPLETE

**Script**: `R/10_geographic_prediction.R` · **Output**: `data/followup/10_geographic_prediction.rds`
· **Figures**: `Figures/followup/10_map_pc1.png`, `10_geodist_vs_profile.png`,
`10_permutation_null.png` · session info saved inside the output RDS
(`$session_info`).

### Inputs and parameters

- Response: `pp_residual_ancestry.rds$Resid` (20 populations × 20,807 units,
  leave-one-chromosome-out ancestry residuals — the completed, unmodified
  Module output). 20,791/20,807 units retained (16 dropped: 15 had exactly 1
  NA population-unit cell, 1 had zero variance across populations).
- Predictor: centred latitude/longitude from `data/bioclimatic_variables.csv`
  (Nyrhispera-renamed as above), 2 columns + intercept.
- Primary response: per-unit standardized (`scale(center=TRUE, scale=TRUE)`
  by column) so a few high-variance units don't dominate. Sensitivity:
  unstandardized.
- Global R²/adj-R² via direct multivariate OLS (`X'X` inverted once, shared
  across all 20,791 unit-columns — fast, exact, no external package).
- Permutation test: 10,000 reps, **complete population-row permutation** of
  the response matrix against fixed geography (never per-unit).
- Uncertainty: chromosome-block bootstrap (2000 reps, resampling the 26
  chromosomes with replacement, reusing already-computed per-unit RSS/TSS —
  same convention as `pp_block_bootstrap.R`).
- Leave-one-population-out: (a) in-sample R² with that population dropped
  and geography re-centred on the remaining 19; (b) genuine out-of-sample
  prediction — fit on 19, predict the 20th, pool RSS/TSS across all 20 folds
  into one cross-validated R² (not an average of 20 unstable per-fold ratios).
- Seed: `set.seed(20260913)` at script top (all permutation/bootstrap draws).

### Numerical results

| statistic | value |
|---|---|
| Global R² (primary, standardized) | 0.092 |
| Global adjusted R² | **−0.015** |
| Global R² (sensitivity, unstandardized) | 0.088 |
| Permutation null: mean / 95th pct | 0.105 / 0.125 |
| Permutation p-value | **0.846** |
| Chromosome-block bootstrap 95% CI on R² | [0.090, 0.094] |
| Leave-one-population-out in-sample R² range | 0.089 (Nyrhispera74 dropped) – 0.111 (Sielva dropped); **every one of the 20 adjusted-R² values is ≤ −0.0004 (i.e. none positive)** |
| Leave-one-population-out **cross-validated** (out-of-sample) R² | **−0.465** |
| Descriptive Mantel (Spearman, geo-distance vs profile-distance) | ρ = −0.014, perm p = 0.936 |
| Population-space PCA of residual profiles | PC1 = 7.6%, PC2 = 7.5%, PC3 = 7.2% (no dominant axis) |
| PC1 vs Latitude / Longitude (Spearman) | −0.522 / −0.332 (a secondary, descriptive nuance — see limitations) |
| Per-unit geographic R² (secondary, descriptive) | median 0.065, 95th pct 0.266; vs FST ρ=−0.041; vs local recombination ρ=−0.013; by sort_class: unsorted 0.063, aquilonia 0.104, polyctena 0.170, unresolved 0.093 |

### Verdict

**No detectable geographic organization of residual ancestry profiles.**
The observed global R² (0.092) sits *below* the median of its own
permutation null (0.105) — geography predicts residual profiles no better
than chance. Adjusted R² is negative both overall and in every one of the 20
leave-one-population-out refits (i.e. this is not one bad population driving
a near-miss). The genuine out-of-sample (cross-validated) R² is strongly
negative (−0.465): predicting a held-out population's residual profile from
its geographic coordinates does *worse* than just using the training-set
mean. This conclusion holds independent of Sielva specifically — the
Sielva-dropped row is not an outlier among the 20 leave-one-out results
(consistent with Sielva being one genuinely disjunct population among 20,
not a data problem — see the correction note above).

### Limitations / caveats

- n = 20 populations, 2 predictors — low power for anything short of a
  strong effect; a weak-but-real geographic signal could still be missed.
  The permutation and cross-validation results are the load-bearing evidence
  here (not just the point-estimate R²), and both independently point the
  same direction.
- The PC1-vs-latitude correlation (ρ=−0.52) is a real, if weak, descriptive
  pattern in the *dominant axis of population variation specifically* (only
  7.6% of total variance) that does not translate into overall predictive
  power. Worth keeping in mind for Analysis 4 (parental structure) but does
  not on its own support "geographic organization" given the formal test's
  null result.
- Sielva is a genuine, geographically disjunct (Alpine) population, not a
  data error (see correction note above) — it is included in the primary
  analysis as one of the 20 populations, per the brief's own design, and its
  inclusion/exclusion does not change the verdict.
- Only linear (2-predictor) geography was tested, per the brief's own
  instruction to keep any spatial-ordination extension small and
  sensitivity-only; not implemented here since the primary linear test was
  already clearly null.

### Interpretation impact

**Strengthens** the case against simple, ongoing, geographically-structured
introgression as the explanation for heterogeneous locus-specific sorting
(pushes toward Outcome A in the brief's decision framework) — but this is
only one of the three legs (geography / individual influence / recent
ancestry runs) that framework requires; Analyses 2–3 are needed before
drawing that conclusion. **Suitable for the main manuscript** as a concise
negative control result (global R², permutation p, cross-validated R²);
the secondary per-unit/sort_class breakdown and PC1-latitude nuance belong
in the supplement only.

---

## Analysis 2: individual-influence and recent-backcross sensitivity — COMPLETE

**Script**: `R/11_individual_influence.R` · **Output**: `data/followup/11_individual_influence.rds`
· **Figures**: `Figures/followup/11_influence_by_population.png`,
`11_influence_vs_ancestry_het.png`, `11_before_after.png` · session info
saved inside the output RDS.

### Inputs and parameters

- 165 hybrid individuals (3–20 per population), oriented genotypes at each
  unit's representative SNP, reconstructed exactly as `pp_prep_units.R` does
  (same inputs: `module_di25/data/di25_inputs.rds`,
  `module_di25_rho05/data/di25_clustering_cM5_rho05.rds` +
  `di25_sorting_emlg_rho05.rds`) but keeping per-individual values instead
  of collapsing straight to the population mean (`Fmat` was never saved at
  individual resolution). Reconstruction verified exactly against `Fmat`
  (max |recomputed population mean − Fmat| = 1.1×10⁻¹⁶ for a spot-checked
  population) before proceeding.
- **Influence is closed-form, not a per-individual refit**: for individual
  *i* in population *p* with *n* observed individuals at a given unit,
  `f_p(−i) − f_p(all) = (f_p(all) − g_i) / (n−1)`, computed once for all 165
  individuals × 20,791 units simultaneously. Summarized per individual as
  the RMS shift across all units (a continuous statistic, no arbitrary
  threshold), alongside genome-wide ancestry, heterozygosity (orientation-
  invariant, from raw dosage), missingness, and population.
- Targeted exclusion scenarios (not all 165 individuals — kept tractable):
  drop the single most influential individual overall; drop each
  population's own most influential individual (20 separate one-at-a-time
  scenarios); drop the 5 most influential overall (one combined scenario);
  drop all of Sielva. Population-level leave-one-out is **reused** from
  `pp_extra_robustness.rds`, not recomputed.
- Residualization against genome-wide ancestry is **approximated** per
  scenario (explicitly sanctioned by the brief): the leave-one-chromosome-
  out ancestry covariate is recomputed fresh from the perturbed population
  means, but the original per-unit `alpha_u`/`beta_u` regression
  coefficients are reused rather than refitting ~20,800 `lm()` calls per
  scenario.
- Seed: none needed (no stochastic step — all recomputation is closed-form
  or exact matrix algebra).

### Numerical results

| individual influence | value |
|---|---|
| influence_rms range (165 individuals) | 0.012 – 0.114, median 0.030 |
| Most influential overall | `Svan41_a` (Svanvik1, n=3 — the smallest population) |
| Most influential individual, by population | Svanvik1 (0.114) > Jarvenpaa (0.085) > Nyrhispera75 (0.060) > Katiskoski (0.059) > ... > **Sielva ranks 19th of 20** (0.025) > LangholmenW (0.015, largest population, n=20) |

Individual influence is strongly, mechanically related to population sample
size (visible directly in the by-population figure): the three smallest
populations (Svanvik1 n=3, Jarvenpaa n=4) produce the largest per-individual
influence values, and this is *not* fully removed by the `1/(n−1)`
normalization — smaller populations still show up disproportionately among
"most influential individuals" simply because there is less averaging-out,
not because any individual is biologically unusual. Concretely, all 3
Svanvik1 individuals and 2 of Jarvenpaa's 4 make up the global "top 5" —
discovered when the "drop top5 overall" scenario briefly emptied Svanvik1
entirely and required handling (see script). **Sielva shows the tightest,
most homogeneous within-population influence distribution of any population**
— consistent with it being a single F1-like colony, and indicating that
Sielva's earlier-noted distinctiveness (PC1 loading, heterozygosity) is a
population-level property, not driven by one unusual individual within it.

| scenario | median FST | ρ(FST, local \|r\|) | ρ(FST, local signed r) | residual PC1 % | geoR2 |
|---|---|---|---|---|---|
| baseline (none dropped) | 0.2513 | 0.1416 | 0.1442 | 7.81 | 0.0920 |
| drop top1 overall (Svan41_a) | 0.2507 | 0.1356 | 0.1414 | 8.70 | 0.0912 |
| drop top5 overall | 0.2479 | 0.1400 | 0.1464 | 9.28 | 0.0985 |
| drop all Sielva | 0.2536 | 0.1451 | 0.1450 | 8.12 | 0.1109 |
| drop top1-per-population (20 scenarios) | 0.250–0.254 | 0.134–0.146 (range across all 20) | 0.141–0.146 | 7.8–9.3 | 0.091–0.098 |
| (reused) existing population-level LOO, FST vs \|r\| | range [0.125, 0.153], full sample 0.142 |  |  |  |  |

Maximum drift in the headline statistic (ρ(FST, local |r|)) across every
scenario tested: **0.006** (from 0.142 baseline).

Note: dropping all of Sielva slightly *increases* geoR2 (0.092→0.111) —
matching Analysis 1's own leave-Sielva-out geoR2 (0.111) almost exactly,
despite Analysis 2 using an independent, approximate residualization
pipeline. This cross-validates the approximation and confirms Sielva is, if
anything, mildly *suppressing* rather than inflating any geographic signal.

### Verdict

**Stable — population-level, not individual-driven.** No single individual,
no population's own most-influential member, the top-5 combined, nor all of
Sielva materially change the headline FST-vs-local-concordance statistic
(max drift 0.006 against a baseline of 0.142) or the broad FST/PCA/geographic
pattern. The result is not a recent-backcross or single-influential-
individual artefact.

### Limitations / caveats

- Influence is a within-population, mean-shift statistic; it cannot detect
  an individual who is unusual in a way that does not shift their
  population's mean much (e.g. one of several similarly-admixed recent
  migrants) — a genuine limitation, not addressed by this analysis.
- Downstream recomputation used an **approximate** residualization (original
  per-unit regression slopes, not refit) — cross-validated against
  Analysis 1's independent, non-approximated geoR2 for the Sielva-dropped
  case (0.111 vs 0.109, agreeing closely), which supports the approximation
  being adequate for a sensitivity check, though a full refit was not
  performed for every scenario.
- Influence is confounded with population sample size (see above) even
  after the `1/(n−1)` normalization; "most influential individual per
  population" should not be read as "most biologically unusual individual"
  without accounting for this.
- Only 4 combined-scenario recomputations plus 20 single-population-drop
  scenarios were run (not all 165 individuals individually) — a deliberate
  scope decision per the brief's own "targeted" framing, not exhaustive.

### Interpretation impact

**Strengthens** the case that the observed heterogeneous, locus-specific
sorting reflects population-level structure rather than a sampling/pedigree
artefact from a handful of individuals — a second leg of the brief's Outcome
A supporting evidence, alongside Analysis 1's geographic null. **Suitable
for the main manuscript** as a concise robustness statement (max statistic
drift under targeted exclusion); the full by-scenario table and the
sample-size confound in the influence metric belong in the supplement only.

## Analysis 3: ancestry-run proxy for recent introgression — COMPLETE

**Script**: `R/12_ancestry_run_proxy.R` · **Output**: `data/followup/12_ancestry_run_proxy.rds`
· **Figures**: `Figures/followup/12_run_length_distribution.png`,
`12_long_run_fraction_by_population.png`, `12_run_metric_associations.png`
· session info saved inside the output RDS.

### Authoritative local-ancestry calls: none found

Searched the repository (grep, case-insensitive, across every `.R` file) for
"local.ancestry", "ancestry.hmm"/"ancestry_hmm", "loter", "elai",
"tract.call", "ancestry.tract", "hapmix", "rfmix", "lamp" before writing any
code. Every hit is an incidental mention of "ancestry-tract" as a
*theoretical* quantity in comments about demographic simulation (e.g. a
Haldane-mapping-function comment in `dev/R/moduleD_ohta_dmi.R`, a "tract
clock" comment in `dev/R/moduleE_analyze_sweep.R`) — no implementation of,
or saved output from, an actual local-ancestry caller exists anywhere in
this repository. Per the brief, this analysis therefore implements **only**
a conservative, explicitly-labelled ancestry-**run proxy** on unphased
genotypes — never called an "ancestry tract" in any output.

### Inputs and parameters

- Individual-level oriented genotypes reconstructed exactly as in Analysis 2
  (re-verified against `Fmat`, max error 1.1×10⁻¹⁶), classified per unit into
  homozygous-aquilonia / homozygous-polyctena / heterozygous / missing.
- Uses the same **representative-SNP/unit-level** genotypes as the rest of
  this module (one marker per LD-reduced unit) — not all 51,612 raw DI25
  SNPs, which would just reintroduce the LD pseudoreplication the unit
  construction exists to remove.
- Units ordered by **genetic (cM) position** (`pp_recombination.rds`'s
  interpolated cM, not physical position), per chromosome, per individual.
- A run = maximal set of consecutive (cM-ordered) units with identical
  genotype state; heterozygous runs are their own category, never merged
  with either homozygous state. **Missing genotypes break a run in the
  primary analysis** (implemented by recoding NA to a per-position-unique
  sentinel before `rle()`, forcing a break at every missing call — not
  bridged).
- Two named sensitivity analyses (not primary): (a) a one-marker-gap
  bridging rule; (b) a stricter marker-informativeness cutoff
  (`current_map_DI > −15` instead of the full DI25 `>−25` panel; 4,238/20,807
  units retained).
- Long-run fraction reported at **4 prespecified cM thresholds** (0.5, 1, 2,
  5cM), not one arbitrary cutoff, as a fraction of the fixed, panel-wide
  callable genetic map (4,606 cM total).
- Inspected `module_di25/R_legacy/di25_pruning_test.R` for implementation
  ideas (not reused as authoritative): notably, its finding that raw
  tract/run length **alone** cannot separate "F1 + scattered genotyping
  noise" from "early backcross", addressed there with a spatial
  (Wald-Wolfowitz) clustering test — **not implemented here** (out of scope
  for this pass), flagged below as a natural follow-up.
- Population-level associations: Spearman ρ + leave-one-population-out range
  (n=20 populations, the primary biological replicate per the brief).

### Numerical results

- **Most runs are very short**: median run length is 0 cM for essentially
  every individual (thousands of single-unit runs per individual — the
  finer rho05 clustering means adjacent units are often not in perfect LD).
  This is why the brief's multi-threshold design was necessary — the
  distribution's upper tail, not its median, carries the signal.
- **Sielva is a dramatic, clean outlier — far fewer long runs than any other
  population**: mean fraction of genetic map in runs ≥2cM = 0.009 (Sielva)
  vs. 0.042–0.106 for the other 19 (next-lowest: Åland, 0.042). Visually
  obvious in the by-population figure. Consistent with Sielva's known
  F1-like, elevated-heterozygosity biology (near-total heterozygosity
  mechanically leaves almost no homozygous runs of any length) — a sanity
  check that the proxy is measuring something real, not noise.

| target | Spearman ρ (pop-level, n=20) | leave-one-pop-out range |
|---|---|---|
| individual influence (Analysis 2) | 0.206 | [0.093, 0.358] |
| residual-profile magnitude | **0.686** | [0.639, 0.786] |
| geographic residual PC1 (Analysis 1) | −0.012 | [−0.144, 0.096] |
| geographic residual PC2 (Analysis 1) | 0.135 | [0.044, 0.284] |
| contribution to high-DI differentiation | **0.759** | [0.732, 0.830] |

(run metric used throughout: mean fraction of genetic map in runs ≥2cM per
population; other thresholds available in the saved output and behave
similarly.)

### Interpretation — important caveat on the two strong associations

The two robust associations (residual-profile magnitude, ρ=0.69; contribution
to high-DI differentiation, ρ=0.76, both with leave-one-out ranges well clear
of zero) should **not** be read as independent confirmation that recent
introgression drives the differentiation pattern. Both are close to
**tautological**: a population with more/longer ancestry-homozygous runs is,
by construction, a population that is more homozygous/fixed for one parental
ancestry at more loci — which is mechanically related to having larger
among-population variance (the "contribution to high-DI differentiation"
target) and larger locus-specific departure from its own genome-wide ancestry
(the "residual-profile magnitude" target). These three statistics substantially
overlap in what they measure, so a strong correlation between them is
expected regardless of *when* that fixation happened (ancient sorting since
original hybridization vs. ongoing/recent backcrossing) — it does not by
itself distinguish the two.

The genuinely diagnostic observation for **recency** specifically is the run
*length* distribution itself: median run length ≈0 cM, with only a modest
right tail extending to a few cM, argues against widespread, currently-ongoing
large-scale backcrossing as the dominant driver — under standard
recombination-clock logic, recent (1–2 generation) backcrossing leaves much
longer, less-fragmented ancestry blocks than observed here; many generations
of recombination since an older admixture event would fragment blocks down
to roughly what is seen. This is consistent with the pattern reflecting
long-standing differentiation/sorting rather than fresh, ongoing gene flow —
though, per the required caution below, this unphased proxy cannot rule out
recent backcrossing in specific individuals.

The weak, borderline associations (individual influence: ρ=0.21, LOO range
just clears zero at the low end; geographic PC2: ρ=0.14, similarly
borderline) are not strong enough to interpret either way. Geographic PC1
is clearly null (consistent with Analysis 1).

### Verdict

Run-length distributions are dominated by short runs genome-wide (no
population shows evidence of pervasive, currently-ongoing large-scale
backcrossing). Sielva is a clean, expected outlier (near-zero long-run
content, consistent with F1-like biology already established elsewhere in
this module). The two strong population-level associations with
differentiation-related targets are most parsimoniously explained by shared
measurement overlap (both reflect "how fixed is this population"), not
independent evidence for a recency-specific mechanism.

### Limitations / caveats (required statement)

**Unphased ancestry runs on ~20,800 LD-reduced units cannot distinguish all
recent backcross histories** from other explanations — an F1-like individual
can show short homozygous runs by pure chance without any backcrossing;
conversely genotyping noise or phase ambiguity can fragment a genuine long
run into several short ones. This is a conservative **proxy**, not a
validated local-ancestry tract call. If an authoritative tract caller becomes
available later, retain this analysis as a sensitivity check, not as
definitive tract inference. Other limitations:

- The Wald-Wolfowitz spatial-clustering test used by the (non-authoritative)
  legacy script to separate "scattered homozygosity" from "clustered
  homozygosity" was not implemented here — would strengthen the recency
  argument if added.
- The one-marker-gap-bridging and stricter-informativeness sensitivity runs
  were computed (saved in the output RDS) but not exhaustively compared
  against the primary result in this writeup — spot-checked only, no
  qualitative difference noticed in the totals.
- Population-level associations (n=20) are correlational; no permutation
  p-values were computed for the 5 associations (leave-one-out range was
  used as the primary uncertainty measure, per the brief's explicit
  guidance for n=20).

### Interpretation impact

**Partially strengthens** the case against widespread recent backcrossing as
the dominant explanation (short run lengths genome-wide), but the two
strong differentiation-linked associations are **not** independent evidence
either way once their measurement overlap is accounted for — this analysis
is weaker support for Outcome A than Analyses 1–2, and should be presented
with the tautology caveat prominently if used. **Suitable for the
supplement**, not a standalone main-text claim; the Sielva sanity-check
result (paragraph above) is the cleanest, most citable single finding here.

## Analysis 4: within-species parental structure across DI — STOPPED (no verified metadata)

Per the brief's own explicit instruction: *"If locality or colony metadata
cannot be recovered, stop this component and document that within-species
geographic differentiation cannot be estimated from the current metadata."*
That condition is met.

### Metadata search performed (before writing any analysis code)

- `data/hybrids_and_parents_maf005.Rdata`'s `sample_data_with_parents`: the
  30 parental rows carry only `Population` (collapsed to the single label
  `aquilonia_parent` or `polyctena_parent` for all 15+15 individuals) and
  `Sample_ID`; `PC1`/`PC2`/`Mitotype` are `NA` for every parental row. No
  colony, locality, or coordinate field of any kind.
- Repo-wide search for dedicated parent-metadata files (`find` for
  `*parent*`, `*metadata*`, `*sample_info*`, `*colony*`, `*locality*`/
  `*localities*`): nothing beyond genotype/LD objects already accounted for
  above.
- `data/Sample_info_outlier_analysis.txt`, `data/Sample_covariate_info_outlier_analysis_20.txt`,
  `data/populations.txt`: all three are **hybrid-only** (165 rows / 20
  population names) — no parental rows, no locality fields.
- `manuscript/Methods.tex`: describes analysis methods (LD-reduction,
  sorting classification, marker-level genomic summaries) but has no
  sample-collection / study-system section describing where the parental
  reference individuals were sampled.
- The parental `Sample_ID` strings themselves (e.g. `Faqu_CBAQ1_1w`,
  `Faqu_CF14a_1w`, `Fpol_Wies316_b`, `Fpol_Ned22_a`) visibly contain
  what look like site/colony abbreviations, and several IDs plausibly
  share a site prefix (`CBAQ1`/`CBAQ2`/`CBAQ3`; `Wies316`/`Wies318`/`Wies545`).
  **Per the brief's explicit instruction, this is not treated as locality
  metadata** — the brief anticipated exactly this situation ("do not infer
  sampling locality solely from Sample_ID abbreviations") and a
  pattern-matched guess would carry unknown risk of being wrong (typos,
  reused codes, non-geographic abbreviations).

### Decision (2026-09-13, user confirmed)

Presented three options — stop and move to synthesis; supply the metadata if
it exists outside this repository; or proceed with the ID-derived pattern
explicitly caveated as unverified. **User chose to stop this component and
move directly to Analysis 5.** No locality-FST, within-species-by-DI, PCA,
or parental-source-vs-hybrid-residual exploratory test was computed — all
of steps 1–8 in the brief structurally require locality, which is not
available.

### Standing caveat (documented regardless of the stop decision)

DI (the diagnostic index defining the whole DI25 high-DI panel this module
is built on) was estimated using **these same 15+15 parental individuals**
(`data/hybrids_and_parents_maf005.Rdata`'s parent rows are the aqu_pops/
pol_pops orientation source throughout `module_di25`). Any future
within-parental analysis using this DI panel would therefore be **partly
circular** by the brief's own reasoning, and would need an externally
defined DI panel or a leave-locality-out/cross-fitted DI to avoid it — worth
keeping in mind if locality metadata is later recovered and this component
is revisited.

### Interpretation impact

**Not evaluated** — this leg of the brief's decision framework (Outcome D/E:
does parental geographic structure explain the high-DI pattern or predict
hybrid residual profiles) remains untested. Should be listed explicitly as
an open question in the synthesis, not silently omitted.

## Analysis 5: synthesis and decision table — COMPLETE

**Script**: `R/14_followup_synthesis.R` · **Output**: `data/followup/14_followup_synthesis.rds`
· **Figure**: `Figures/followup/14_synthesis.png` (4 panels: A. residual-profile
PC1 by geography; B. observed geographic R² vs its permutation null;
C. population-level long-run fraction vs individual influence; D. robustness
panel — FST-vs-concordance ρ under Analysis 2's targeted exclusion scenarios,
used in place of a parental-differentiation panel since Analysis 4 has no
result, per the brief's own fallback rule) · session info saved inside the
output RDS.

### Combined decision table

| item | value |
|---|---|
| Geographic effect size + permutation | R²=0.092 (adj R²=−0.015); permutation p=0.846 (10,000 reps; observed R² **below** the null median of 0.106) |
| Leave-one-population-out range | in-sample adj R²: all 20 folds ≤ −0.0004 (none positive); out-of-sample cross-validated R² = −0.465 (strongly negative) |
| Largest individual influence | `Svan41_a` (Svanvik1, n=3 individuals), influence_rms=0.114 — driven by small population size, not flagged as biologically unusual |
| Effect of excluding influential individuals | max drift in headline FST-vs-\|r\| statistic across every targeted exclusion scenario = 0.0060 (baseline 0.1416) — stable |
| Ancestry-run association | pop. long-run fraction vs residual magnitude ρ=0.69 [0.64,0.79]; vs high-DI-differentiation contribution ρ=0.76 [0.73,0.83] (flagged as largely tautological — see Analysis 3); vs individual influence ρ=0.21 [0.09,0.36] (weak); Sielva a clean near-zero outlier |
| Parental differentiation by DI | **NOT AVAILABLE** — Analysis 4 stopped: no verified parental Sample_ID→species/colony/locality/coordinate metadata found in this repository (user confirmed 2026-09-13, chose not to proceed on unverified ID-derived locality) |
| Limitations and power notes | n=20 populations throughout (low power for weak effects — load-bearing evidence is the permutation/cross-validation/LOO-range, not point estimates alone); Analysis 2's scenario recomputation used an approximate, not refit, residualization; the Analysis 3 run proxy cannot distinguish all recent-backcross histories from scattered heterozygosity/genotyping noise; DI was estimated from the same parental individuals any future parental-DI analysis would use (circularity caveat, undischarged) |

### Verdict against the brief's Outcome A–E framework

> Evidence from Analyses 1–3 is most consistent with **Outcome A (weak
> support)**: no detectable geographic organization of residual ancestry
> profiles (Analysis 1), no individual- or small-population-driven artefact
> (Analysis 2), and no signature of pervasive, currently-ongoing large-scale
> backcrossing (Analysis 3's short run-length distributions). This
> **strengthens, but does not by itself establish**, an interpretation of
> heterogeneous locus-specific sorting as reflecting locus-specific
> selection or incompatibility resolution rather than ongoing
> geographically-structured gene flow — selection is not established merely
> because these alternative explanations were not supported. **Outcome D/E
> (parental geographic structure) remains completely untested** (Analysis 4
> stopped for lack of verified metadata) and must be reported as an open
> question, not treated as resolved by Analyses 1–3.

### Draft text (NOT inserted into the manuscript)

**Supplementary Methods:**
> To test alternative explanations for the observed heterogeneous,
> locus-specific sorting of ancestry, we conducted three follow-up analyses
> on the DI25 rho05 population-partitioning dataset (20,807 LD-reduced
> units, 20 hybrid populations, 165 individuals). First, we tested whether
> geographic distance predicts each population's leave-one-chromosome-out
> residual ancestry profile, using multivariate OLS (centred
> latitude/longitude as predictors, per-unit standardized profiles as the
> response), a complete-row population-label permutation test (10,000
> replicates), a chromosome-block bootstrap for uncertainty, and
> leave-one-population-out cross-validation. Second, we tested whether the
> population-level pattern is driven by a small number of individuals,
> using a closed-form per-individual influence statistic (the exact
> population-mean shift from excluding that individual, normalized by
> population size) and recomputing the headline FST-vs-local-concordance
> statistic under a small set of targeted exclusion scenarios (most
> influential individual overall; each population's own most influential
> member; the five most influential overall; all of one geographically
> disjunct population). Third, we implemented a conservative, unphased
> ancestry-run proxy: runs of consecutive, genetic-map-ordered LD-reduced
> units with identical ancestry-homozygous genotype state, with missing
> calls breaking (not bridging) a run in the primary analysis, related to
> population-level differentiation and influence metrics via Spearman
> correlation with leave-one-population-out uncertainty (n=20 populations,
> the biological replicate throughout).

**Supplementary Results:**
> Geographic distance did not predict residual ancestry profiles: the
> observed multivariate R² (0.092) fell below the median of its own
> permutation null (0.106, p=0.846), and the leave-one-population-out
> cross-validated R² was strongly negative (−0.465). The population-level
> concordance pattern was not driven by individual or small-population
> artefacts: excluding the most influential individual, the five most
> influential individuals overall, each population's own most influential
> member, or an entire geographically disjunct population (a genuine
> Alpine site among otherwise-Fennoscandian populations) shifted the
> headline FST-vs-local-concordance statistic by at most 0.006 from a
> baseline of 0.142. An unphased ancestry-run proxy found predominantly
> short runs genome-wide, arguing against pervasive, currently-ongoing
> large-scale backcrossing; the one Alpine population showed essentially
> no long ancestry-homozygous runs, consistent with its
> independently-documented F1-like heterozygosity. Two population-level
> associations between run length and differentiation-related statistics
> were strong but attributed to measurement overlap rather than an
> independent recency signal (see Discussion).

**Short main-text result:**
> Three independent follow-up analyses found no evidence that the
> heterogeneous, locus-specific sorting of ancestry documented above is
> attributable to geographic population structure, individual-level
> sampling artefacts, or widespread recent backcrossing.

**Cautious discussion paragraph:**
> These results narrow, but do not close, the space of alternative
> explanations for the observed pattern. Geography, individual influence,
> and coarse recent-introgression signatures were each tested directly and
> found wanting, which is consistent with — though does not establish —
> locus-specific selection or incompatibility resolution as an explanation;
> ruling out several alternatives is not equivalent to confirming the
> remaining one. A fourth alternative, geographic structure within the
> parental reference samples themselves (which could contribute to the
> high-DI pattern via ascertainment or genuine parental-source admixture),
> could not be evaluated: no verified locality metadata exists for the
> 15+15 parental individuals underlying this analysis, and we did not
> infer locality from sample identifiers alone. This is a genuine,
> currently unresolved gap, not a null result, and should be flagged as
> such wherever these follow-up analyses are cited. The ancestry-run proxy
> used here is a conservative, unphased approximation; it cannot
> distinguish all recent backcross histories from chance short runs in an
> F1-like individual, and a validated local-ancestry tract caller, if
> adopted later, should supersede it rather than be treated as
> confirmatory of the present proxy's conclusions.

## Overall status

Analyses 1, 2, 3, 5 complete; Analysis 4 stopped (no verified parental
locality metadata — user-confirmed decision, see above). All code, data, and
figures are under `R/1{0,1,2,4}_*.R`, `data/followup/`, `Figures/followup/`;
nothing in the pre-existing `pp_*.R` pipeline, its `data/`/`Figures/`
outputs, `README.md`, `AUDIT.md`, or `CROSS_MODULE_INPUTS.md` was modified.
