# DI25 among-region association — a pairwise two-locus screen on the high-DI markers

The **intrinsic-side** test of the headline question, re-run on the **high-DI diagnostic
markers only** (DI > −25) and their from-scratch DI25 LD clustering. This mirrors
[`moduleD_among_region_association/`](../../moduleD_among_region_association/) — the
pair-level question of whether **unlinked** LD-reduced units are associated *beyond* what
shared ancestry and relatedness predict (a two-locus Dobzhansky–Muller incompatibility) —
but every unit here summarises **high-DI variation exclusively**, built from
[`module_di25/`](../)'s DI25 clustering + sorting rather than the canonical
all-marker eMLGs and the moduleC gate.

**Headline result (descriptive, pre-simulation): largely negative.** After the structure-
corrected pairwise scan and the honest read-out, the DI25 signal reduces to **two leads** —
one trans/repulsion pair (Chr13–Chr15, opposite ancestry directions) and one low-recomb
cis/coupling pair (Chr2) — on a **diffuse residual-co-ancestry background** (genomic
inflation λ = 1.34 among the sorted candidates). Both leads are strongly among-population, so
the DMI-vs-neutral call is **deferred to Module E's recombination-matched neutral null**.

*Unlinked* = different chromosome, or same chromosome > `LINK_CM = 10` cM on the genetic map.
Scripts run **from the repository root** (`~/gitlab/formica_hybrid`); this module's inputs are
under `module_di25/` and shared `data/`, outputs under
`module_di25/among_region_association/{data,Figures}/`.

---

## Design (how it differs from the canonical Module D)

- **Units** — `module_di25/data/di25_clustering_cM5.rds` folded into one aqu-oriented
  165 × 11,052 dosage matrix (eMLG consensus for >2-marker units, representative SNP for 1–2).
  The 165th hybrid `Hei159_h` (absent from `sample_data`) is recovered via its colony prefix.
- **Gate** — `module_di25/data/di25_sorting_emlg.rds`, a **superset** of `moduleC_C3_cl`
  (`differentiated`, `DI`, `sort_class`, `prop_fixed`, `uni_score`, `bi_score`). `differentiated`
  is near-vacuous here (11,042/11,052) since every DI25 marker is diagnostic.
- **Pairwise CANDIDATE NETWORK, not focal-vs-all.** Rather than test a few foci against all
  ~11k units, a marginal candidate set is tested **only against each other**: **`bi_score ≥ 0.2`
  = 2,083 units** — BIDIRECTIONAL loci (strong sorting in *inconsistent* directions across the
  replicate populations = the equal-fitness/DMI-plausible signature). Because a bidirectional
  locus must segregate to sort both ways, this **also enforces appreciable MAF** (min within-
  hybrid MAF 0.16), which is essential: the earlier `prop_fixed ∈ [0.5,0.8)` gate anti-correlated
  with MAF (cor −0.94 with within-pop MAF) and **selected near-fixed low-MAF loci** (median MAF
  0.15) whose within-pop LD was inflated by a thin heterozygote layer. All **~2.17M unlinked
  pairs** among the bidirectional units are tested; the double-LOCO GRM and top-10 genome PCs are
  still built from **all** differentiated units. *Trade-off:* a marginal-sorting gate is blind to
  ancestry-independent trans-DMIs where both loci look marginally unsorted; a MAF-only filter is
  not selective on DI25 (median MAF 0.38 → ~50M all-vs-all), so bidirectionality is the gate that
  is both DMI-motivated and MAF-enforcing.
- **Ohta functions come from `LDscnR`** (`devtools::load_all`); EMMAX is not in LDscnR, so a
  local `R/emmax.R` copy is used.

## Pipeline (`R/`)

Run e.g. `Rscript module_di25/among_region_association/R/di25D_emmax.R`.

| script | role |
|----|----|
| **`di25D_assemble_units.R`** | STEP 0. Builds the input contract `data/di25D_units.rds` — the aqu-oriented dosage matrix (orientation verified: r = 0.99 vs the gate's `f_aqu_pooled`), the unit-aligned gate, per-unit chr / cM / local cM-Mb, and paralogy inputs. |
| **`di25D_ohta_dmi.R`** | Among-population Ohta decomposition on unlinked pairs → `data/di25D_ohta.rds` + `Figures/di25D_fig4.{pdf,png}`. **Descriptive / Module-E instrument, not a screen.** Reproduces Module D's structure signature on DI25: `D2st ≫ D2is`, `D'2st ≈ 0` = drift. NB **R_st is the D'2st axis** (cor(R_st², Dp2st) = 0.74), **not** D2st (0.02); the epistasis-relevant term is the *systematic* Dp2st, which is ≈ 0 pre-null. Stores both Ohta contrasts (`Ohta_D`, `Ohta_Dprime`). |
| **`di25D_emmax.R`** | The pairwise-candidate structure-corrected scan → `data/di25D_emmax.rds` + `Figures/di25D_emmax_network_summary.png`. Double-LOCO VanRaden GRM + 10 genome PCs (from all differentiated units); pooled BH-FDR `q < Q_FDR = 0.05` over all ~2.17M unlinked pairs → `r2crit = 0.159`, **7 edges** (→ 5 survive leverage); the stringent `q < 0.01` (`r2crit = 0.245`) isolates the single strongest pair (F9614–F9879); cis/trans by the GLS coefficient sign; paralogy-flagged. Also emits three **cheap uncorrected estimates** per pair (`r2_raw` = plain r²; `R_st` = among-population correlation; `within_pop_r` = median within-population correlation) so the correction is transparent. Reports **λ = 1.28** and does **not** genomic-control the FDR (it is a pre-filter; the discriminant is Module E). |
| **`di25D_estimate_pairs.R`** | Scatter-matrix (GGally) of `Rsq`, `r2_raw`, `R_st`, `within_pop_r` → `Figures/di25D_estimate_pairs.png`. Subset (random + all tails + FDR edges highlighted). Shows raw r² ≈ among-pop structure, while EMMAX `Rsq` is orthogonal to `R_st` (0.13) and tracks `within_pop_r` (0.36). |
| **`di25D_emmax_circos.R`** | Candidate circos → `Figures/di25D_emmax_candidate_circos.png`. 2,083 nodes by sorting direction, all chords at r² > 0.10 faint (the diffuse residual co-ancestry); leverage-**surviving** leads bright (by coupling), FDR edges dropped by leverage grey dashed. Uses `circos.link` (point curves), not `circos.genomicLink` (ribbons). |
| **`di25D_lead_snp_heatmaps.R`** | All member-SNP genotypes (raw, not consensus) of the top leads → `Figures/di25D_lead_snp_heatmaps.png`; exposes the true MAF/het structure the consensus hides. |
| **`di25D_rst_heatmaps.R`** | Genotype heatmaps of the top **among-population** (R_st) different-chromosome pairs → `Figures/di25D_rst_heatmaps.png`: parallel/anti-parallel co-sorting across populations (|R_st| 0.86–0.91) but within_pop_r ≈ 0 / EMMAX ≈ 0 — the structure-confounded D′2st candidates (drift vs parallel-DMI, needs Module E). The among-pop counterpart to the within-pop lead heatmaps. |
| **`di25D_detectability_schematic.R`** | Conceptual figure of the DMI **detectability window** → `Figures/di25D_detectability_schematic.png`: within-pop LD (EMMAX) peaks mid-sorting then vanishes at fixation; among-pop LD rises into the shared-founding-confounded blind spot. |
| **`di25D_detectability_data.R`** | Data-grounded companion → `Figures/di25D_detectability_data.png`: (a) within-pop segregation collapses as loci sort (blind spot ~1% here); (b) same-chr LD decay — among-pop flattens by 10 cM, within-pop extends further, leads are long-range within-chr outliers. |
| **`di25D_network_build.R`** | The honest read-out → `data/di25D_network.rds`. From the clean FDR edges: within-chromosome average-linkage **merge** → **WITHIN-population leverage** (leave-one-population-out `min |median within-pop r| ≥ LEV_LOO = 0.3`) → low-recomb **annotation** (`rec_pct < 0.10`; carried to E, never a filter) → cross-chromosome **modules**. *The leverage is within-population, NOT among-population (R_st): the decomposition showed EMMAX `Rsq` tracks the within-pop signal, so an among-pop leverage would select for structure and discard the DMI candidates.* Emits the candidate DMI **lead table** (among-pop R_st + its LOO still reported). |
| **`di25D_lead_heatmaps.R`** | Genotype view of the surviving leads → `Figures/di25D_lead_heatmaps.png`. Per lead: individuals × the two loci (grouped BY POPULATION, labelled) + the 3×3 joint-genotype contingency (obs : exp). Descriptive only. |
| **`di25D_paralogy.R`** | Shared cross-chromosome duplication filter (`flag_paralogy`; median within-population \|r\| > 0.9). Copied from Module D; cluster-generic. |
| **`emmax.R`** | EMMAX LMM engine (Li, Kemppainen, Rastas & Merilä), copied verbatim (not exported by LDscnR). |

## Results

- **Ohta arm:** the high-DI data reproduces Module D's `D2st ≫ D2is` with `D'2st ≈ 0` — drift/
  structure, not systematic epistatic LD, pre-null.
- **Pairwise EMMAX:** **λ = 1.28** on the bidirectional candidate set. This is **not** a broken
  K+PC model — it is the candidate *selection*: strongly-sorted loci share **diffuse residual
  co-ancestry** (the whole null lifted, no fat tail). At the working `q < 0.05` over ~2.17M pairs
  (`r2crit = 0.159`) **7 pairs** clear the pooled FDR (→ 5 survive leverage); the stringent
  `q < 0.01` (`r2crit = 0.245`) leaves just the single strongest pair.
- **Estimate decomposition** (`di25D_estimate_pairs.R`): raw r² ≈ among-pop structure
  (`cor(|r|raw, |R_st|) ≈ 0.84`), EMMAX `Rsq` orthogonal to `R_st` (~0.12) and tracking
  `within_pop_r` (~0.35). This motivated the **within-population** leverage (an among-pop R_st
  leverage would select for structure and discard the DMI candidates).
- **The low-MAF lesson:** an earlier `prop_fixed ∈ [0.5,0.8)` candidate gate produced ~4 leads
  that were all **low-MAF, near-fixed one direction** (MAF 0.05–0.13), with within-pop LD carried
  by a thin heterozygote layer (`cor(min-MAF, |within_pop_r|) = −0.33`; 74% of edges had min-MAF
  < 0.15) — a **low-MAF inflation artifact**, which the within-pop leverage does *not* catch.
  Switching the gate to `bi_score ≥ 0.2` enforces MAF ≥ 0.16 and removes the artifact.
- **Read-out (bidirectional set):** at `q < 0.05` the within-population leverage keeps **5 of 7**
  FDR edges → **5 leads**, all appreciable-MAF (0.27–0.48) and robust across **11–13 populations**,
  and — see below — **all same-chromosome**. The standout, and the sole `q < 0.01` survivor, is
  **F9614 – F9879 (both Chr7, > 10 cM apart), coupling:** within-pop r = **0.68** across 11
  populations (LOO 0.59), MAF 0.27 / 0.45, R_st 0.50 (both within- and among-pop signal), both
  marginally `unsorted` but opposite ancestry directions — the equal-fitness bidirectional signature;
  its `Rsq = 0.245` sits at the stringent r²crit (borderline). The DMI-vs-neutral call is deferred to
  **Module E**.

## Detectability, linkage scale, and power

- **All EMMAX leads are same-chromosome; no robust trans-DMI.** Across FDR thresholds, every
  within-population lead is a **same-chromosome** pair > 10 cM apart (up to 92 cM); the single
  different-chromosome (true trans) FDR edge fails the leverage filter. `di25D_detectability_data.R`
  panel b shows why: **among-population** LD (R_st²) decays to the inter-chromosomal baseline by
  ~10 cM (the admixture-LD scale that justifies the "unlinked" rule, as the Ohta arm validated),
  but **within-population** LD keeps declining well past 10 cM, and the leads sit far above the
  average as long-range within-chromosome outliers (low-recomb / structural). So the 10 cM rule is
  valid for the among-pop axis but does **not** make a same-chr pair trans on the within-pop axis.
  Reading: the surviving leads are **long-range within-chromosome LD** (structural / intra-
  chromosomal coadaptation), not trans-DMIs; the genuine trans (different-chromosome) signal is empty.
  *(A clean trans-DMI screen would restrict `unlinked` to different-chromosome pairs only.)*

- **Detectability window (power).** An equal-fitness DMI builds **within-population** LD while both
  loci still segregate → EMMAX detects it; once both loci fix within every deme the within-pop
  variation (and thus the detectable LD) vanishes, leaving only the **among-population** co-sorting
  direction, which is confounded by the shared founding LD common to all 20 replicates. So we are
  powered for **in-progress** DMIs and blind to **completed** ones (`di25D_detectability_schematic.R`).
  Crucially, `di25D_detectability_data.R` panel a shows the blind spot is **currently small** — only
  ~1% of loci are near-fixed (`prop_fixed > 0.8`); sorting is mostly incomplete, so most loci (and all
  candidates/leads) sit in the detectable window. The near-empty result therefore reads as **"no
  detectable in-progress pairwise DMIs"** and is **not** primarily a blind-spot artifact — though
  completed DMIs remain the principled blind spot.

## Inputs

- `module_di25/data/di25_clustering_cM5.rds`, `di25_sorting_emlg.rds`, `di25_inputs.rds`
  (the DI25 clustering, sorting gate, and 012 genotypes).
- `data/hybrids_and_parents_maf005.Rdata` (`sample_data_with_parents` — population labels),
  `data/Frufa_DTOL_PR.ref_genome.recmap` (genetic map).
- `LDscnR` via `devtools::load_all("~/gitlab/LDscnR/")` (Ohta + `consensus_dosage`).

## Status

Descriptive stage **complete**. The two leads are queued for **Module E**'s recombination-
matched neutral null, which makes the DMI/coadaptation call.
