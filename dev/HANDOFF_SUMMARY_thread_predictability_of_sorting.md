# Handoff summary — "predictability of sorting" thread

Written for merging with the LD-clustering / eMLG / parallelism thread into
one pipeline. Distilled conclusions and validated code, not a conversation
replay. Companion doc: `dev/HANDOFF_SUMMARY_thread_LD_and_parallelism.md`
(shares the same eMLG-overwrite problem and the same `parallelism_stats.R`).

## Critical, time-sensitive issue (same overwrite as the other thread)

**`data/eMLG_5loci_0025.rds` was silently overwritten in place on
2026-07-21 01:05** (now `params$use_cM=TRUE, cM_threshold=0.5`, **474,014
groups / 32,840 has_eMLG clusters**, vs. 469,361 / 32,943 in the old
bp-backstop file with no `params` field). Group IDs are NOT stable across
this change — the same string (`"F10329"`, `"F59859"`, etc.) refers to a
different marker set now. Anything joining on `group_id` across the 01:05
boundary silently gives wrong answers.

**This thread's results are split across BOTH file versions**, so even
setting staleness aside, nothing here is currently internally consistent on
one clustering. Per-script verdict (by execution time, which I witnessed
live in-session, not just file mtime):

| script | last run | eMLG file used | verdict |
|---|---|---|---|
| `dev/R/parallelism_stats.R` | (library, not a run) | none | **stable** — pure statistic, no eMLG |
| `dev/R/parent_LD_diagnostic.R` | 2026-07-20 | none | **stable** — no eMLG |
| `dev/R/parent_LD_DI_sweep.R` | 2026-07-20 | none | **stable** — no eMLG |
| `dev/R/tier1_architecture.R` | 2026-07-20 | OLD | DI/recomb/π/d_xy **stable**; cluster-size numbers **STALE** |
| `dev/R/tier2_sorting_vs_recomb.R` | 2026-07-20 | OLD | unit-level results **STALE** (qualitative "flat" likely robust) |
| `dev/R/supplementary_analyses.R` | **never run whole**; doc numbers inherited from the two above (07-20, OLD) | OLD | eMLG-dependent numbers in the PDF/HTML **STALE** |
| `dev/R/sorting_baypass_crosscheck.R` | 2026-07-21 ~10:10 | **NEW** | **current** (but old-file-inconsistent with Tier1/2) |
| `dev/R/sorting_baypass_crosscheck_memberlevel.R` | 2026-07-21 ~10:20 | **NEW** | **current** |

Non-eMLG findings (all of §"Sorting statistics" below, plus DI↔recomb,
DI↔π, d_xy in Tier 1, plus both parent-LD scripts) are **not affected** —
they read only `GTs`, `map_hyb_005`, the recombination map, and BayPass BF
files. **For the merged pipeline: re-run everything eMLG-dependent on one
tagged clustering file** (endorse the other thread's version-tagging fix).

## Scope

What "predictability of sorting" investigates: given 20 replicate *F.
polyctena × F. aquilonia* hybrid populations plus allopatric parental
references, does parental ancestry sort **predictably** — the *same* parent
fixing across populations, as expected under directional selection / an
unequal-fitness Dobzhansky–Muller incompatibility (DMI) — or in **random
directions**, as expected under drift or an equal-fitness incompatibility?
This is the per-locus, 20-replicate generalisation of the 3-population
parallelism in Nouhaud et al. 2022 (*PLOS Biol*). The thread built the
per-locus sorting statistics, characterised the genomic architecture that
could confound them (differentiation vs recombination), and cross-checked
the sorting against the BayPass PC1/PC2 climate-association outliers.
Everything is **descriptive**; the neutral null that would license
selection/DMI claims is designed but not built (see open questions).

## Validated / reusable code (all uncommitted / exploratory, in `dev/`)

### `dev/R/parallelism_stats.R` — the core statistic (authored this thread)
Per locus, from a `prep` (= `ohta_fast_prepare(GTs, pops)`) plus named
parent populations:
- Orients dosage into an aquilonia-allele frequency from the parents.
- Per hybrid population: near-fixed toward aquilonia (freq ≥ 1−`fix_th`),
  polyctena (≤ `fix_th`), or unsorted (`fix_th`=0.1).
- **Decomposition** `prop_fixed = |uni_score| + bi_score`, with
  `uni_score=(n_aqu−n_pol)/n_obs`, `bi_score=2·min(n_aqu,n_pol)/n_obs`, both
  population fractions. One threshold `sort_th` (default 0.5) classifies
  `sort_class ∈ {aquilonia, polyctena, bidirectional, unsorted}`,
  comparably across SNPs, eMLGs and genome proportions.
- **Parallelism test**: `n_aqu ~ Binomial(n_fixed, null_prob)`, two-sided,
  exact (`.binom_two_sided_p`, vectorised). `null_prob=0.5` (symmetric).
  Reports `p_binom`/`q_binom`.
- **Differentiation gate**: `min_DI` (primary; pass `DI=setNames(map$DI,
  map$marker)`, aligned by name) and/or `min_parent_diff` (secondary,
  default 0). Prevents near-monomorphic loci from faking parallelism.
- **NB design point flagged by the other thread**: `sort_class` is assigned
  by the `sort_th` threshold, **not** gated by `p_binom`/`q_binom`. This is
  deliberate (descriptive class vs inferential test), but the two threads
  should agree a convention (see open questions).
- Requires parental genotypes → uses `data/hybrids_and_parents_maf005.Rdata`
  (built by `R/LD_decay_from_DIEM.R`, per the other thread).

### `dev/R/tier1_architecture.R` — genomic architecture of differentiation
Correlates DI with recombination (map-interpolated cM/Mb), within-species π,
absolute divergence d_xy, relative F_ST, and eMLG cluster size; recombination-
decile table; Cruickshank–Hahn (relative vs absolute divergence vs π). Fig
`Figures/tier1_architecture.png`. Cluster-size parts STALE (old eMLG).

### `dev/R/tier2_sorting_vs_recomb.R` — sorting vs recombination, unit level
Runs `parallelism_stats` on LD-pruning representatives (one tag SNP per
cluster = independent unit) vs a SNP-level sample; bins by recombination
decile; regresses `prop_fixed ~ recomb + DI (+ cluster_size collider)`. Fig
`Figures/tier2_sorting_vs_recomb.png`. Unit-level results STALE (old eMLG).

### `dev/R/supplementary_analyses.R` — consolidated reproduction script
Sections 0–6 sourcing `parallelism_stats.R`+`Ohta.R`, reproducing all
figures/numbers. **Never executed end-to-end**; if run now it uses the NEW
eMLG file and its eMLG-dependent numbers will differ from the PDF.

### `dev/R/sorting_baypass_crosscheck.R` / `_memberlevel.R` — sorting × BayPass
Overlap of directionally-sorted eMLG clusters with BayPass PC1/PC2 outlier
clusters (≥10 members BF(dB)≥20), size-adjusted. Representative-level
(`crosscheck.R`) and member-level (`_memberlevel.R`). **Ran on the NEW
file** → current, but inconsistent basis vs Tier1/2.

### `manuscript_notes/supplementary_methods.{html,pdf}`
6-page supplementary methods (self-contained, base64 figures, Chrome-print
PDF). Contains the STALE Tier1/2 numbers — regenerate after re-running.

### Also: `dev/R/parent_LD_diagnostic.R`, `dev/R/parent_LD_DI_sweep.R`
Within-species vs admixture composite-r² decay by recombination, and across
DI thresholds. No eMLG → stable. Feed the null-design decision (below).

## Key findings

**Sorting statistics & null (stable — no eMLG):**
- At `sort_th`=0.5 over 16,077 DI>−15 SNPs: **89.9% unsorted, 6.9%
  aquilonia, 2.7% polyctena, 0.4% bidirectional** — a directional,
  aquilonia-biased excess (matches the other thread's genome-wide ~1.8–2:1
  aqu:pol).
- **The per-locus "pooled" null is circular** — its null prob correlates
  0.96 with the observed outcome and yields 0 significant loci vs ~8,000
  under symmetric. Genome-wide mean aquilonia ancestry at diagnostic loci is
  ~0.50, so `null_prob=0.5` is the correct, non-circular choice. (This is
  the answer to the other thread's open question #3 about trying "pooled":
  don't — it's degenerate. A single global scalar is the valid
  admixture-skew correction, and here it ≈ 0.5.)
- **DI gate**: at DI>−15, 100% of loci have parental Δp ≥ 0.5 (medians 0.93/
  0.90/0.89/0.83 at DI>−15/−20/−25/−30) — DI and empirical parental
  differentiation agree across the working range.

**Tier 1 — architecture (DI↔recomb/π/d_xy stable; cluster-size STALE):**
- **DI is essentially orthogonal to recombination** (Spearman −0.03, flat
  across deciles). Its dominant correlate is **low within-species π** (−0.46,
  partly definitional), not elevated d_xy (+0.13). DI↔cluster_size +0.23
  (STALE number).
- A genuine differentiation-island signature (elevated F_ST *and* d_xy,
  reduced π, huge clusters) is confined to the **lowest recombination
  decile** and is largely orthogonal to DI.
- **Cluster size is a collider** of recombination and diversity → must not be
  conditioned on when relating DI/recomb/sorting. (Contrast: the DI-
  enrichment work below correctly *does* adjust for it, because there it is a
  confounder — see that section.)
- Consequence: **DI ⟂ recombination means they are separable predictors of
  sorting** — less confounding than initially feared.

**Tier 2 — sorting vs recombination (STALE basis; qualitative likely robust):**
- **At the independent-unit (eMLG) level, sorting is flat vs recombination**
  (`prop_fixed ~ recomb`: β≈0.000, p=0.97); the SNP-level relationship is the
  spatial-redundancy artifact that unit-counting removes. So single-unit
  sorting is **not** low-recombination-concentrated — consistent with the
  neutral expectation, and contrary to a naive window-based read.
- Small, separable DI effect (β=−0.014, p<1e-16), consistent across
  recombination tertiles. Adding cluster size (collider) spuriously induces a
  recombination effect — a clean demonstration of why to exclude it.
- Theory note established this thread: neutrally, a non-recombining autosomal
  block does **not** have lower Nₑ than surrounding regions (unlike mtDNA,
  whose low Nₑ is from ploidy/uniparental inheritance, not linkage). "Stronger
  drift in low recombination" is either (a) linked selection lowering local Nₑ
  (needs selection, not neutral) or (b) spatial book-keeping (larger blocks =
  fewer independent units = bigger per-event footprint). eMLG-level counting
  removes (b); a recombination-matched null is needed for (a).

**Sorting × BayPass cross-check (current — NEW file):**
- Directionally-sorted eMLG clusters overlap BayPass **PC2** (not PC1)
  climate-association outlier clusters. Representative-level (single tag SNP)
  gave size-adjusted **OR ≈ 16** — but this **overstated it**.
- **Member-level (fraction of member SNPs directionally sorted) is the
  honest number: OR ≈ 2.9 (all outliers) / 3.8 (PC2), p ≈ 0.004–0.012**,
  wide CIs, n=27 clusters. Outlier clusters are heterogeneous — several are
  fully sorted (e.g. small aquilonia clusters), but the largest mega-blocks
  (e.g. a 6,263-locus cluster, frac 0.04) barely sort. Predominantly
  aquilonia-leaning. The size confound at member level runs *against* the
  raw signal (bigger clusters co-sort more), so the adjusted OR is smaller
  than the raw rate difference (0.45 vs 0.09), not larger.
- Interpretation: a **modest, real** co-localisation of directional sorting
  with the PC2 climate axis — evidence that *extrinsic* (ecological)
  selection is associated with part of the sorting, distinct from the
  intrinsic-DMI hypothesis. Not a headline result; a lead.

## Elaboration: `R/analyse_baypass.R` + `R/diagnostic_index_enrichment.R`

These are pre-existing repo scripts (not mine) but I read them fully this
thread; fuller account than the commit-message level:

**What the BayPass layer is.** BayPass associates per-population allele
frequency with **environmental covariates** — `pc1_env`/`pc2_env` in
`R/write_baypass_inputs.R`, per-population climate PCs (constant within a
population). So BF(dB) is a genotype–environment association (climate
adaptation), **not** a genetic-structure axis derived from the same
genotypes — which is what makes the sorting cross-check non-circular. Runs
cover 2 population sets (`with_aland`, `aland_excluded`) × 2 PCs × 2 Omega
settings = 8 result files (`{PC}_DIEM_{omega}_summary_betai_reg.out`).
`withOmega` uses a pre-estimated population-covariance (structure) matrix;
`noOmega` lets BayPass estimate it internally — comparing the two is the
structure-control robustness check. Import is a **positional** join
(MRK==row index; the file has no marker ID), guarded by a `stopifnot` on row
count and 1:N order (has caught a real stale-data mismatch before).

**`analyse_baypass.R` = the outlier definition + visualisation.** Restricts
to eMLG-filtered markers (has_eMLG, n_loci≥5); an **outlier cluster** = an
LD cluster with **≥10 individually significant (BF≥20) members**. Rationale
(load-bearing): requiring broad within-cluster support filters out clusters
where only a few "lucky" markers cross threshold inside an otherwise
unremarkable LD block (a likely false positive); filtering on cluster *size*
alone (n_loci≥10) does not declutter, since size is independent of
association strength. Output = one background(grey)+foreground(colored-by-
cluster) Manhattan per PC×config, with **globally stable, position-ordered
cluster colors** so the same region is the same color across figures.

**`diagnostic_index_enrichment.R` = the DI-enrichment test on top.** Asks
whether those outlier clusters are enriched for high-DI (DI>−25) loci vs
background. The **confound it corrects**: the outlier definition *requires*
n_loci≥10 (mechanically size-biased), and cluster size independently
predicts DI genome-wide (OR=1.135 per doubling, Spearman 0.23 — the same
0.23 as Tier 1). A size-floor match is insufficient (outlier median ~147
loci vs size-matched background ~15); the primary test is a **marker-level
logistic regression** `hi ~ is_outlier + log(n_loci)` (the outlier
coefficient = enrichment surviving continuous size control), reported naive
vs adjusted in a forest plot. **Conclusion (from its last run, on OLD
clustering — needs re-verification): PC1 enrichment survives the size-
adjusted test in all 4 settings; PC2 splits (2 robust, 2 reverse below 1).**

**Cross-thread synthesis worth flagging.** The DI-*enrichment* (static
ancestry-informativeness ↔ climate outliers) is **PC1-driven**, whereas my
sorting-*overlap* (dynamic directional sorting ↔ climate outliers) is
**PC2-driven**. Different climate axes tag different phenomena. Both the
DI-enrichment and my cross-check treat cluster size — one as a confounder to
regress out (correct there), one as a collider to exclude (correct in Tier
2) — a nice illustration that the same variable needs opposite handling by
DAG position. Note both BayPass analyses read `eMLG_5loci_0025.rds`, so
their outputs are also on a pre-01:05 clustering and need re-running.

## Open questions / not-yet-done

1. **Re-run all eMLG-dependent code on one tagged clustering file** (Tier1
   cluster-size, Tier2, supplementary, and re-confirm the cross-checks) so
   the thread is internally consistent. Regenerate the supplementary PDF.
2. **The among-region (two-locus) DMI test is designed but NOT built** —
   screen *unlinked* eMLG pairs for elevated among-population LD (Ohta
   `D'2st`, via the repo's `dstat_unphased_*` / `run_C_eMLG_ohta`). This is
   the analysis that actually targets DMIs rather than generic sorting;
   Tier 2 showed single-unit sorting is *not* where the signal is.
3. **The neutral null is not built.** Design is settled (haplodiploid SLiM,
   empirical recombination map, founders seeded from **real haploid-male
   haplotypes** — not independent binomial draws, because within-species LD
   is 60–85% of admixture-scale LD in low-recombination regions and threshold-
   proof, per the parent-LD scripts; calibrate demography on LD-decay/π/F_ST,
   keeping sorting stats out of the calibration set). Until it exists, every
   "flat/consistent-with-neutral" statement is descriptive only.
4. **`sort_class` vs `p_binom` gating** — agree one convention across threads
   (descriptive threshold class vs binomial-significance-gated).
5. **Reconcile PC1 (DI-enrichment) vs PC2 (sorting-overlap)** divergence, and
   **identify what the PC1/PC2 climate axes actually are**.
6. `fix_th` sensitivity for the smallest hybrid populations (3–20 diploids);
   member-level cross-check size-effect sign flip; a `score_eMLG()`-style
   dilution check on the heterogeneous outlier clusters (ties to the other
   thread's open question #2).
