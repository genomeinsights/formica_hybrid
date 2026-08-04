# Module A — the sorting-threshold (τ) sensitivity series

**Scripts:** `moduleA_sorting/R/moduleA_tau_sensitivity.R` (comparison),
`parallelism_stats.R::classify_sort()` (shared classifier), `moduleA_cluster_sorting.R`
(stamped eMLG annotations). **Figure:** `Figures/moduleA_tau_sensitivity.png`.

## The series

Ancestry sorting is reported over a **predefined three-point τ series**, with
everything else held fixed:

| τ | role | captures |
|---|---|---|
| **0.5** | relaxed diagnostic | moderately recurrent sorting, with more direction-unresolved calls |
| **0.6** | **primary (reported)** | the operating point |
| **0.8** | stringent | only highly repeated sorting |

Fixed across the series: **φ = 0.85** (`fix_th = 0.15`), **parental MAF ≥ 0.15**,
**cM05** clustering, **binomial** direction test at **α = 0.05**, **DI ungated**.
Every result is parameter-stamped `tau05` / `tau06` / `tau08` (e.g.
`moduleA_cluster_sorting_tau06.rds`). Because the per-unit fixation counts are
**τ-independent**, the series is a cheap post-hoc reclassification of stored counts —
no re-run of the expensive prep.

**Manuscript rule:** τ=0.6 defines the reported analysis; 0.5 and 0.8 show what is
gained (relaxed) or lost (stringent). We never select whichever τ gives the
strongest result. A result seen **only at τ=0.5** is described as threshold-sensitive
and exploratory; a result **surviving τ=0.8** is the strongest evidence (on a smaller
candidate set, so wider CIs and reduced enrichment power).

## Results (per-SNP and LD-reduced unit)

**Prevalence** — only prevalence should change:

| level | τ | aquilonia | polyctena | direction unresolved | unsorted | aqu:pol |
|---|---|--:|--:|--:|--:|--:|
| SNP | 0.5 | 105,979 | 126,375 | **12,851** | 410,332 | 0.84 |
| SNP | 0.6 | 70,431 | 83,015 | **1,171** | 500,920 | 0.85 |
| SNP | 0.8 | 15,307 | 18,578 | **8** | 621,644 | 0.82 |
| unit | 0.5 | 46,353 | 64,647 | 3,142 | 138,260 | 0.72 |
| unit | 0.6 | 32,697 | 45,068 | 238 | 174,399 | 0.73 |
| unit | 0.8 | 8,259 | 10,411 | 0 | 233,732 | 0.79 |
| eMLG | 0.5 | 2,304 | 1,939 | 16 | 22,947 | 1.19 |
| eMLG | 0.6 | 1,254 | 1,043 | 0 | 24,909 | 1.20 |
| eMLG | 0.8 | 152 | 119 | 0 | 26,935 | 1.28 |

(The aqu:pol ratio is stable *within* each resolution across τ; the eMLG consensus
leans aquilonia (~1.2) where individual SNPs lean polyctena (~0.83) — a resolution
effect, not a threshold effect. eMLG-level direction unresolved is 16 / 0 / 0.)

**Effect sizes** — direction model `P(aquilonia) ~ zDI + zr + zmaf + zcs` on the
directionally resolved units, and low-recombination depletion:

| τ | n unidir | zDI (95% CI) | z | zr | DI/recomb | sorted dec-1 | sorted rest |
|---|--:|---|--:|--:|--:|--:|--:|
| 0.5 | 111,000 | 1.39 [1.37, 1.41] | 127 | −0.09 | 16× | 0.259 | 0.467 |
| 0.6 | 77,765 | 1.57 [1.54, 1.59] | 112 | −0.09 | 17× | 0.184 | 0.318 |
| 0.8 | 18,670 | 1.81 [1.75, 1.87] | 57 | −0.10 | 18× | 0.052 | 0.075 |

(The magnitude slope `prop_fixed ~ recomb` is τ-independent: **+0.052, t = 117**.)

## The four questions

1. **Does direction / broad pattern persist?** Yes. The aquilonia:polyctena ratio
   is stable (~0.83 at SNP level across all τ), DI predicts direction with the same
   sign and increasing strength, and low-recombination regions are depleted for
   sorting at every τ.
2. **Do effect sizes stay comparable?** Yes. `zDI` = 1.39 / 1.57 / 1.81 and the
   DI-vs-recombination ratio 16× / 17× / 18× — the same effect, if anything cleaner
   at the stringent threshold.
3. **Does only prevalence change?** Yes — unsorted grows and the sorted set shrinks
   as τ tightens, as expected.
4. **Does anything appear only at the relaxed threshold?** Yes — **direction unresolved**
   (12,851 → 1,171 → 8 at SNP level; 16 → 0 → 0 at the eMLG level). It is therefore
   **threshold-sensitive / exploratory** (see `moduleA_direction_unresolved_peaks.md`: the
   τ=0.5 peaks are LD-cluster artifacts). Moreover, ~70% of the τ=0.5 direction-unresolved
   units are **power-limited directional leans** (median split 80:20 at n_fixed 10–11)
   for which the binomial test cannot resolve a direction at low n_fixed. This
   category is not evidence of genuine both-direction sorting. Everything that survives τ=0.8 — the
   DI-governed directional sorting — is the strongest evidence, on a ~19k-unit set
   (z=57, wider CIs) rather than the ~78k of the primary analysis.

**Bottom line:** τ=0.6 is the reported analysis; the flanks show that the directional
DI signal is robust to threshold choice, while direction-unresolved calls are
concentrated at the relaxed threshold and are exploratory. The manuscript reports
the τ=0.6 analysis.
