# Module D — the minimal pipeline (from-scratch spec)

*Distilled after full exploration (see `moduleD_wip_summary`, memory `moduled-results`).
This is the only-what-we-need version: one scan, four filters, one diagnostic, one
collapse, one candidate classification, one null. Everything before the null is
descriptive.*

## The one question
Is there **pair-specific two-locus incompatibility** (DMI / coadaptation) between *unlinked*
eMLG clusters, beyond (a) generic ancestry sorting (Modules A/B/C) and (b) **neutral
admixture LD given the real recombination map**?

## Inputs & locked definitions
- **Units** — eMLG consensus dosages, `data/eMLG_5loci_0025_cM05.rds` (164 hybrids × clusters).
- **Gate / covariates** — `data/moduleC_C3_cl.rds` (`differentiated`, `DI`, `sort_class`).
- **Recombination map** — `data/Frufa_DTOL_PR.ref_genome.recmap` (cumulative cM + cM/Mb).
- **Raw genotypes** — `data/hybrids_only_maf005.Rdata` (paralogy filter, heatmaps only).
- **Unlinked** = different chromosome **OR** same chromosome **> `LINK_CM = 10` cM**
  (genetic map, not physical distance — admixture LD decays on the map).
- **Trait / partner universe** = `differentiated` eMLGs.

## Step 1 — One arm: EMMAX structure-corrected two-locus scan
Each focal eMLG dosage as a quantitative trait, tested against every differentiated partner
in an LMM. Three conditioning choices, all load-bearing:
- **K = double-LOCO VanRaden GRM from *differentiated* eMLGs** — exclude the focal **and**
  partner chromosome. K must *contain* the confound it removes (the selected ancestry axis);
  a neutral/LD-pruned K under-corrects (λ inflation, ancestry leaks back).
- **Condition each focal on the top `N_PC = 10` genome PCs** (Covar) — strips secondary
  structure axes; random-control focal rate drops to ~0 hits with this on.
- **cis/trans sign = the GLS coefficient** (not mean-centred whitened variables).
- **Focal set** — a random genome-wide sample (honest rate / calibration) **plus** targeted
  hypotheses (e.g. bidirectional, high-marginal) for candidate enrichment. Always keep the
  random control. A full all-vs-all (~10⁹ tests) is unnecessary.

*We do **not** run the Ohta among-population arm as a discovery screen: the residual signal is
individual-level, so the among-population lens is dominated by the ancestry cline. `Ohta.R` is
**retained only** as the LD-measurement instrument for Module E's null.*

## Step 2 — Read-out filters (one honest procedure, in order)
1. **Unlinked only** (cM rule).
2. **Paralogy** — drop pairs with within-population |r| > `PARALOGY_R = 0.9` (assembly
   duplicates; a distance rule cannot catch cross-chromosome ones).
3. **Global BH-FDR q < `Q_FDR = 0.01`** across *all* unlinked tests (one whole-experiment
   correction, not per-focal Bonferroni).
4. **Leverage** — keep an edge only if its among-population |r| survives
   leave-one-population-out (`LEV_LOO = 0.3`). This dissolves single-deme / low-MAF
   pseudo-hubs on its own, so **no separate near-fixed MAF cutoff is needed**.

## Step 3 — Recombination annotation (NOT a filter)  ← revised
Low recombination is **not** evidence of neutrality. It is exactly the habitat where
coadapted/epistatic complexes accumulate and persist (recombination suppression is a
*hallmark* of coadaptation). For an unlinked pair the between-region recombination fraction
is 0.5 regardless of local rate, so low local recombination does not "hold" cross-chromosome
LD; it makes each region a **clean, low-noise ancestry block**, i.e. the highest-signal
marker of *whatever* latent factor produces residual covariation — neutral finer structure
**or** selection. After conditioning out global ancestry (Step 1 PCs), that residual is a
real puzzle, and selection is a live hypothesis.

Therefore:
- **Annotate** each locus/region by recombination percentile (`STRUCT_PCT = 0.10` → `structure`
  flag) as a **"high neutral-baseline"** marker — context for the null, **not** a discard.
- **Carry every locus, low-recombination included, into the candidate screen and the null.**
- Characterise the structure regions **per module** (Step 4), not pooled: pooling all
  structure loci mixes several distinct co-ancestry modules and has no coherent joint
  interpretation, so no genome-wide "structure PCA" is used.

## Step 4 — Collapse redundancy (honest counting / display)
- **Third-level merge** — `MERGE_METHOD = "average"`, within-chromosome ≤ `LINK_CM` cM,
  cut at |r| > `MERGE_R = 0.5`; >`LINK_CM` cM pairs given an un-mergeable distance. (Single
  linkage / connected-components **chains** uncorrelated ends — do not use.)
- **Cross-chromosome modules** — `MODULE_R = 0.4`, average-linkage on meta-node consensus
  across chromosomes. Shows the co-ancestry is several disjoint blocks, each a co-inherited
  haplotype block; per-module heatmaps for the coherent single-module view.

## Step 5 — Candidate classification
Candidates = leverage-stable, appreciable-MAF pairs, **regardless of recombination class**.
Classify by ancestry orientation:
- **same-ancestry cis = co-sorting** → hand to Module A/C (e.g. the F33454 aquilonia module).
- **ancestry-independent trans (both `unsorted`) = pairwise-DMI candidates** (F11431–F49480 …).
- **low-recombination modules** (the coherent multi-chromosome blocks) = **coadaptation
  candidates** to be tested against the recomb-matched null — *not* set aside.

## Step 6 — Module E: the recombination-matched neutral null (the only discriminant)
Neutral-relic co-ancestry vs selection-maintained LD is settled **only** by a null matched
to *local* recombination:
- Simulate neutral admixture with the **real recombination map** (SLiM on the true map, same
  admixture history, same founder ancestries, matched sample size / generations), then run the
  **identical** EMMAX read-out (Steps 1–2) on the simulated genotypes.
- A pair/module is a **candidate DMI/coadaptation** if its observed LD **exceeds** the
  neutral null **at its own local recombination rate** — so a low-recombination pair is judged
  against a legitimately *high* baseline and can still surface.
- Matching options (decide in E): (a) simulate directly on the real per-chromosome map so
  every region carries its own rate (preferred — no binning artefacts); (b) if per-region
  simulation is too costly, bin observed pairs by local cM/Mb and compare each bin to
  rate-matched simulated pairs.

### Corroboration available before E (both point at selection if present)
1. **Replicate-population parallelism** — 20 hybrid populations. Neutral relic LD drifts in
   *random directions* across replicates; epistatic selection maintains the *same* combination
   reproducibly. A module whose direction/ancestry-combination is parallel across many
   populations is a selection signature low recombination alone cannot produce.
2. **Specific exceedance over the low-recombination background** — is a given module *tighter*
   than the general low-recombination neutral expectation? Test the coherent multi-chromosome
   blocks first.

## Explicitly dropped (dead-ends)
- Ohta among-population arm as a **standalone screen/ranking** (structure-dominated; the
  D′2st ≠ R_st trap). `Ohta.R` kept only as E's instrument.
- Bidirectional SNP-level test (weak, diffuse; its permutation null overstates significance).
- F33028 climate-tie & inversion hypotheses (the recombination annotation + E subsume them).
- Per-focal Bonferroni; neutral/pruned GRM basis; a separate near-fixed MAF cutoff; any
  pre-null candidate-DMI network (hairball).
- Treating the `structure` label as a **filter** (the corrected point above).

## Code shape — the minimal file set
**Keep (core):** `emmax.R`, `moduleD_emmax.R`, `moduleD_paralogy.R`,
`moduleD_network_build.R` (FDR → paralogy → average-merge → leverage → recomb annotation →
modules; **remove the near-fixed/MIN_MAF check**),
`moduleD_module_heatmaps.R`, `moduleD_network_circos.R`.
**Keep for E only:** `Ohta.R`, `moduleD_ohta_dmi.R` (null instrument — not a screen).
**Drop after the WIP restructure:** `moduleD_bidirectional.R`, `moduleD_climate_tie.R`,
`moduleD_recomb_env.R` (logic folded into the build), `moduleD_nonhub_trans.R`,
`moduleD_structure_heatmap.R` (superseded by per-module), and — decide — `moduleD_emmax_independence.R`
(QC), `moduleD_trans_heatmaps.R` (retain only if the WIP keeps those two heatmaps).
**Never touch:** `*_PK.R` (PK's copies), Module E's shared-tree files, `dev/run_cM1.R`.
