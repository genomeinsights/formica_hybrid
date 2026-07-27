# Module C — association with climate (the extrinsic arm)

**Question.** Do the regions that carry ancestry-diagnostic and/or directionally
sorted loci coincide with regions whose allele frequencies track **climate**? This
tests the *extrinsic* (climate-adaptation) hypothesis for predictable sorting.

**Headline result.** The climate signal is **static, not dynamic**:
ancestry-**diagnostic** loci are modestly but consistently **enriched** in
climate-associated regions on both climate axes, whereas directional **sorting**
shows **no** such association. The extrinsic link is with static differentiation, not
with dynamic directional sorting.

### Code & data

| Script (`R/`) | Role | Reads → Writes |
|---|---|---|
| `moduleC_analyse_baypass.R` | C1: outlier definition + Manhattans | clustering + BayPass `.out` → `manhattan_*.png` |
| `moduleC_diagnostic_index_enrichment.R` | C2: DI-enrichment of outlier clusters | clustering, `hybrids_only`, BayPass `.out` → `diagnostic_index_enrichment<tag>.csv`, forest/proportion figures |
| `moduleC_ancestry_confound.R` | ancestry↔climate confound check | `hybrids_and_parents` → `moduleC_ancestry_confound.rds`, figure |
| `moduleC_sorting_climate.R` | C3: sorting × outlier overlap (threshold/binary) | `hybrids_and_parents`, `hybrids_only`, clustering, `moduleA_snp.rds`, BayPass `.out` → `moduleC_C3_cl.rds` (**consensus checkpoint**), `moduleC_sorting_climate_<tag>.rds`, `moduleC_fig3_<tag>.{pdf,png}` |
| `moduleC_rate_based.R` | **primary analysis** (size-normalised) | `moduleC_C3_cl.rds`, clustering, `hybrids_and_parents`, BayPass `.out`, `diagnostic_index_enrichment<tag>.csv` → `moduleC_rate_based_<tag>.rds`, `moduleC_dose_response_<tag>.{pdf,png}` |

Run 2026-07-22/23. Threshold-tagged outputs carry the settings in the filename
(`_5_15` set is the one reported here). The BayPass inputs/runs are upstream on HPC:
`moduleC_prepare_{with_aland,aland_excluded,sielva_excluded}.R` →
`moduleC_write_baypass_inputs.R` → `<set>/run_baypass.sh` →
`<set>/*_summary_betai_reg.out`.

---

## Methods

### Genotype–environment association (BayPass)

Per-population allele frequencies were associated with **per-population climate
covariates** — two principal components (PC1, PC2) of the climate variables, constant
within a population — using **BayPass**, with evidence reported as **Bayes factors on
the deciban scale**, `BF(dB)`. Eight runs cross:

- **two population sets** (full; and **Åland-excluded**, the primary set — see the
  ancestry confound below),
- **two climate axes** (PC1, PC2),
- **two treatments of the among-population covariance matrix Ω** — either estimated
  beforehand from the **LD-pruned** marker set under the core model and supplied to
  the association run (`withOmega`, preferred, since LD among markers otherwise
  inflates apparent among-population covariance), or estimated internally from the
  full marker set (`noOmega`).

**Primary configuration: Åland-excluded × withOmega.** All eight runs are reported so
that dependence on either choice is visible.

**Positional matching, guarded.** The association output carries no marker identifier,
only a row index into the genotype file. Results are matched positionally, and every
import is guarded by assertions that the number of rows equals the number of markers
and that the index is strictly increasing — so a mismatch from regenerated inputs
raises an error rather than silently assigning evidence to the wrong markers.

### Ancestry as a potential confounder

Sorting *is* ancestry fixation, so a climate axis correlated with genome-wide ancestry
could co-localise with sorting for structural rather than ecological reasons. Ancestry
was estimated directly (per population, the mean *aquilonia*-allele frequency across
diagnostic markers, oriented from the parents) and correlated with each climate axis;
because a correlation across 20 populations can rest on very few, it was recomputed
with individual populations removed.

> **PC1** is only weakly ancestry-correlated (Spearman **ρ = −0.23**);
> **PC2** is moderately so (**ρ = +0.57**), driven predominantly by a single
> population, **Åland**.

Two controls follow and are applied throughout: (i) conditioning on Ω absorbs the
shared structure such an axis is correlated with; (ii) the **Åland-excluded** set is
the primary analysis, with the full set retained as a comparator.

### Association strength per cluster — a rate, not a count

Each cluster's climate association is summarised as its **significance rate**: the
proportion of its member markers reaching `BF(dB) ≥ 15`.

This is the central methodological point of Module C. A **count** of significant
members makes a cluster's *eligibility depend on its size*:

> a **9-marker** cluster (about the median) would need **56%** of its members
> individually significant to reach a count of five, whereas a **90-marker** cluster
> needs only **~6%**.

The same nominal criterion is therefore far stricter for small clusters. This
**cannot** be repaired by conditioning on cluster size, because the resulting
eligibility is a hard, strongly non-linear function of size — the adjustment would
have to extrapolate across a range where one outcome cannot occur at all. Expressing
association as a **rate** removes the dependence by construction; cluster size is
retained as a covariate for its independent effects.

Because a rate is estimated imprecisely when a cluster has few markers, results are
reported across **nested size strata** — all clusters, ≥ 20 markers, ≥ 50 markers —
rather than pooled, making the dependence of any estimate on measurement precision in
the predictor explicit.

### Testing at the cluster level

Clusters are the units of observation. For the **DI outcome**, the quantity modelled
is the proportion of a cluster's members that are ancestry-informative (**DI > −25**),
each cluster contributing one observation; a size-weighted model with an
**overdispersion correction** is reported alongside as a check. A marker-level model
is **not** the primary analysis — it treats LD-correlated markers within a cluster as
independent trials, which understates the standard errors substantially. The
sorting–climate relationship is tested the same way, with the cluster's sorting
classification as the outcome.

---

## Results

All findings are **descriptive** — climate associations are genotype–environment
*correlations*, and the sorting itself is licensed as selection only against the
forthcoming neutral null (Module E). *The identity of the climate axes — which
variables load on PC1 vs PC2 — is pending and is a separate author's task.*

**Ancestry-diagnostic loci concentrate in climate-associated regions.** Among clusters
large enough for the significance rate to be well estimated (**≥ 50 markers**), the
proportion of ancestry-diagnostic members increases with climate association:

> each **+10 percentage points of significance rate** corresponds to
> **+4.5 pp** of diagnostic loci for **PC1** (*p* = 8×10⁻⁴) and **+3.7 pp** for
> **PC2** (*p* = 0.013), against a baseline of roughly 10%.

A size-weighted, overdispersion-corrected model gives concordant estimates (**odds
ratios 1.50 and 1.47** per +10 pp). The effect is attenuated in the intermediate
stratum (**≥ 20 markers: +1.9 and +1.8 pp**) and is **not detectable when all clusters
are pooled**, exactly as expected when the predictor is measured with substantial error
in small clusters. **Both climate axes behave alike, at similar magnitude.** (The
informative ≥ 50-marker DI stratum comprises **878 clusters**.)

**Directional sorting shows no such association.** Applying the identical
size-normalised, cluster-level test to directional sorting yields **no relationship**
in the size-comparable stratum:

> **odds ratios 0.99 (PC1) and 1.03 (PC2), *p* > 0.95**; a member-level analysis is
> likewise null.

This is expected on mechanistic grounds: genotype–environment association is detected
through variation in allele frequency **among** populations, whereas unidirectional
sorting drives most populations toward the *same* near-fixed state and therefore
**minimises exactly that variance**. Climate-driven selection would in any case
predict populations diverging in *different* directions according to local conditions,
not the parallel fixation that unidirectional sorting describes. (The ≥ 50-marker
sorting stratum has only 451 clusters — "no evidence survives", not "proven absent".)

### Figures

- **Fig. 3 — `../Figures/moduleC_fig3_5_15.png`** (also `.pdf`; unsuffixed
  `moduleC_fig3.*` is an earlier version). Sorting × climate-outlier overlap under the
  cluster-level test.
- **Dose–response — `../Figures/moduleC_dose_response_5_15.png`** (also `.pdf`).
  Diagnostic-loci proportion (and the null sorting relationship) as a function of
  climate significance rate, across the nested size strata — the size-normalised,
  cluster-level enrichment that replaces the size-gated outlier count.
- **Ancestry confound — `../Figures/moduleC_ancestry_confound.png`** (also `.pdf`).
  Per-population ancestry vs each climate axis (PC1 ρ = −0.23; PC2 ρ = +0.57,
  Åland-driven) — motivates conditioning on Ω and taking the Åland-excluded set as
  primary.
- **DI-enrichment — `../Figures/diagnostic_index_enrichment_forest_5_15.png`,
  `../Figures/diagnostic_index_enrichment_proportions_5_15.png`.** Forest plot of the
  enrichment effect and the underlying proportions.
- **Manhattans — `../Figures/manhattan_{PC1,PC2}_aland_excluded_withOmega_eMLG_clustered_min5SigLoci_withBackground.png`**
  (primary config; the full 16-file grid spans both population sets × PC × Ω ×
  min5/min10). Cluster-level climate association along the genome, outliers over
  background.

---

## Interpretation and caveats

The climate signal in these data is a **static** one: regions whose allele frequencies
track climate are enriched for ancestry-diagnostic loci, on both climate axes, by a
modest but consistent margin. There is **no corresponding dynamic signal** —
directional sorting is not preferentially located in climate-associated regions.

Limits to keep in view: the informative stratum is **878 clusters**, the enrichment is
a few percentage points against a variable background, and these are **associations,
not demonstrations of selection**. Identifying the climate axes and evaluating the
sorting against the recombination-matched neutral null (Module E) remain the decisive
next steps; the intrinsic alternative is pursued separately through among-population LD
between unlinked diagnostic clusters (Module D).

> **A methodological note worth carrying forward:** never binarise on a size-gated
> count. The eligibility of a "≥ N significant members" rule scales with cluster size
> in a hard, non-linear way that no size covariate can undo. Use a rate
> (`n_sig / n_loci`), keep size as a covariate, and report nested size strata so the
> reader can see how the estimate depends on measurement precision in the predictor.

---

*This completes the descriptive arc (Stage 0 → A → B → C). The neutral null (Module E)
and the intrinsic-DMI screen (Module D) are separate workstreams — see the
[index](README.md#the-scientific-question).*
