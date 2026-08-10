# module_di25 (high-DI analyses)

Figures for the plenary talk, focused on **Module A** (ancestry sorting) with the
high-DI diagnostic markers (DI > −25).

Subfolders mirror the Module A layout: `R/` `data/` `doc/` `Figures/`.
Scripts run **from the repo root** (`~/gitlab/formica_hybrid`).

**Manuscript figures:** `Figures.R` (at the folder root, not `R/`) is the home for
publication-ready figures — **no descriptive titles/subtitles**, only axes, legend,
panel tags and minimal condition labels; each `fig_*()` writes a PDF (cairo, vector)
and PNG (300 dpi) via `save_fig()`. Conventions: interpretable τ as "≥ N% of
populations near-fixed" (plotmath ≥), italic species names, wide line-type legend keys
(`theme_ms`). First figure `fig_recomb_directional` = directional sorting vs
recombination at two thresholds (tagged **a**, two sub-panels) + directional
block-bootstrap recombination coefficient per τ (**b**). The exploratory scripts below
(in `R/`) keep their descriptive titles.

**Deck palette (species):** *F. aquilonia* = teal `#21918C`; *F. polyctena* = darker
yellow `#D3C93B`; heterozygous / unresolved = dark viridis `#440154`; unsorted /
not-fixed / missing = grey. The continuous per-population ramp runs teal→yellow
(through green at "balanced"). All figure scripts use these; to re-skin, change the
hex codes in the scripts (`diem_circos_core.R` `DIEM_COL` drives every DIEM circos).

## Pipeline (`R/`)

| script | role |
|----|----|
| `diem_circos_core.R` | Shared renderer `render_diem_circos()` — draws a DIEM genotype circos by direct polar rasterisation. Sourced by both plot scripts so they stay pixel-comparable. |
| `diem_circos.R` | **Per-SNP** DIEM circos. Reads `data/species_diagnostic_markers_DI25_20pops.tsv.gz`; writes `Figures/diem_circos.png`. |
| `diem_circos_eMLG.R` | **Per-eMLG (LD-reduced)** DIEM circos: diagnostic SNPs collapsed into their LD-cluster consensus. Writes `Figures/diem_circos_eMLG.png`. |
| `di25_ld_clustering.R` | **From-scratch** LD-decay → complexity reduction on the DI25 markers ONLY (so clusters summarise high-DI variation exclusively). Hybrids-only LD; cM-cap sensitivity sweep. Writes `data/di25_clustering_cM*.rds`. |
| `diem_circos_di25_eMLG.R` | DIEM circos on the from-scratch DI25 clustering (5 cM), all units: eMLG consensus for >2-marker clusters, representative SNP for 1–2. Writes `Figures/diem_circos_di25_eMLG_cM5.png`. |
| `diem_circos_compare.R` | Side-by-side per-SNP vs LD-reduced (5 cM) DIEM circos, same markers and a shared individual order (ring _k_ = same ant in both). Writes `Figures/diem_circos_compare_snp_vs_ldreduced.png`. Uses `render_diem_circos(new_device = FALSE)` to draw both panels into one canvas. |
| `di25_sorting.R` | Ancestry sorting à la Module A on the DI25 data, per SNP and per eMLG (5 cM). φ=0.85 fixed, τ swept {0.5,0.6,0.7,0.8}, `sort_rule="binom"`, parent-MAF gate 0.15, DI ungated. Writes `data/di25_sorting_{snp,emlg,sweep}.rds`. |
| `diem_circos_sorting_sweep.R` | The τ sweep drawn around the genome in the DIEM circos layout: 4 concentric rings = τ (0.5→0.8 inner→outer), each unit coloured by sort class (aqu purple / pol yellow / unresolved teal / unsorted grey). `SWEEP_LEVEL=SNP|eMLG` env var. Writes `Figures/diem_circos_sorting_sweep_{snp,emlg}.png`. Uses `render_circos_raster()` (general categorical circos added to the core). |
| `diem_circos_di25_sortclass.R` | LD-reduced DIEM circos (genotype rings) with the sort class at τ=0.6 as an **outer annotation ring**, aligned unit-for-unit. Writes `Figures/diem_circos_di25_eMLG_sortclass_cM5.png`. Composites `render_diem_circos(..., r_out=0.80, draw_chr_labels=FALSE)` with `render_circos_raster(..., add=TRUE, bg_col=NA)` (transparent overlay; core gained `add`/`draw_chr_labels`/transparent-bg support). |
| `diem_circos_population.R` | Per-**population** DIEM circos: the 195 individual rings collapse to 20 hybrid-population rings, each cell = the population's mean *F. aquilonia*-allele frequency at the unit (continuous viridis via a 100-step palette). Both levels. Rings inner→outer by overall ancestry. Writes `Figures/diem_circos_population_{snp,emlg}.png`. |

### `di25_architecture.R` (is sorting stronger in low recombination?)

Module A's architecture analysis (Table 3 / Fig S4b / Table 4) restricted to the
high-DI loci. Interpolates cM/Mb per chromosome (`approx(..., rule=2)`), deciles
recombination, compares fraction sorted at SNP vs eMLG, and regresses
`prop_fixed ~ scale(log10(recomb+0.1)) + scale(DI)` (τ=0.6). **Finding:** among
high-DI loci sorting is stronger at *lower* recombination — the **reverse** of the
full-genome result (magnitude recomb coef −0.032 SNP / −0.015 eMLG vs full-data
+0.05). Caveat baked into the figure: the lowest-recomb deciles hold thousands of
SNPs but only 8–41 independent eMLG units (one giant LD block), so point size ∝ n;
the robust trend is deciles 4–10 (0.12→0.05). Direction still governed by DI
(z=6.1), not recombination (z=0.5, ns). **Statistical test** (eMLG units are
~independent, so one logistic obs/unit is valid): `P(sorted) ~ log10(recomb)`
= −1.24 (p≈1e-18) all units, −1.35 (+DI), and −1.08 (z=−6.2, **p=7e-10**) in the
reliable range (deciles ≥4, n≥100) — so it is NOT a low-n artifact. The figure
puts n-per-decile as scaled background bars (bottom 3 deciles = 10,244 SNPs but
only 57 independent units), Wilson CIs on the points, and shades the low-n
deciles. Writes `data/di25_architecture{,_deciles}.rds`, `Figures/di25_sorting_vs_recomb.png`.

### `di25_direction_flip.R` (the SNP→eMLG direction flip, with a control)

Quantifies the flip (per-SNP the strongly sorted loci lean *F. polyctena*, ~34% aqu
at τ=0.8; as independent units the bias inverts to a stable ~80% *F. aquilonia*). Additive —
reads `di25_sorting.R` outputs, modifies nothing. Panel a plots
`pct_aqu_of_resolved` vs τ for **per-SNP**, **per-eMLG**, a **LD-thinned control**
(one raw representative SNP per cluster, no consensus), and per-eMLG **consensus-only**;
Wilson CIs on every series except per-SNP (which is pseudoreplicated → no valid CI);
`n` shown as a table row beneath the axis; the thinned≈eMLG coincidence is annotated
in-panel. Panel b: cumulative share of each direction's sorted SNPs by cluster rank.
The sorted-set composition strip is a **separate backup figure**
(`di25_direction_flip_composition.png`). **Result:** the thinned control tracks
eMLG, not per-SNP → the flip is driven by **independence** (LD pseudoreplication
removed), **not** by the eMLG consensus construction — confirming the claim.
Mechanism: top-10 clusters supply **70%** of polyctena- vs 49% of aquilonia-directed
sorted SNPs (Chr 5/25/26 blocks). Nuance: consensus-only units (the big blocks)
stay polyctena-leaning (~49% aqu at τ=0.8); the unit-level aquilonia bias is carried
by the many small units. Writes `data/di25_direction_flip.rds`,
`Figures/di25_direction_flip.png`.

### `di25_recomb_tau_sweep.R` (recombination effect swept over τ, honest inference)

Sweeps the sorting-vs-recombination test over τ∈{0.5,0.6,0.7,0.8}. **Inference
note baked in:** LD-reduced units are *not* independent, only much more so than
SNPs — so the per-unit logistic SE is anticonservative. CIs on the recombination
coefficient come from a **genomic block bootstrap resampling whole chromosomes**
(preserves within-chromosome residual LD). Result: the effect (sorting weakly
higher at low recombination) is negative at **every** τ and the block-bootstrap
95% CI excludes 0 at every τ (coef −0.82→−1.68 as τ rises), but stays **small**
(McFadden R² 0.005–0.014; P(sorted) e.g. 16%→5% at τ=0.6). Writes
`data/di25_recomb_tau_sweep.rds`, `Figures/di25_recomb_tau_sweep.png`.

### `di25_region_dmi.R` (LD / missing genotypes between the 3 polyctena blocks)

Tests the epistasis/DMI idea for the Chr 5/25/26 polyctena blocks: do their rare
*aquilonia* variants co-occur across individuals, and are three-region combinations
missing? Per-hybrid oriented *aquilonia* dosage per block (top polyctena cluster:
F8626/F6986/F7174). **Answer: cannot be tested — the signal is population structure.**
The all-aquilonia class is enriched (5 vs 0.03) and mixed classes absent, *but* every
carrier belongs to one of 3 populations (Åland→Chr26, Nyrhisperä75→Chr5, Sielva→all
three); Chr5–Chr25 r = 0.83 → −0.03 once Sielva+Nyrhisperä75 are dropped. Mixed-class
expected counts are <1, so their absence is not evidence. Writes
`data/di25_region_dmi.rds`, `Figures/di25_region_dmi.png`.

### `di25_recomb_directional.R` (sorting vs recombination, by direction)

The recombination-effect figure with the sorted fraction split by direction:
**colour = species** (aquilonia teal / polyctena yellow); **line type = level**
(solid eMLG / dashed SNP). Grey bars = n independent eMLG units per decile; Wilson
CIs on eMLG; low-n deciles shaded. **`TAU` env-overridable** (tau-stamped outputs) —
render one per τ with
`for t in 0.5 0.6 0.7 0.8; do TAU=$t Rscript module_di25/R/di25_recomb_directional.R; done`.
Shows polyctena sorting concentrated in low recombination (the big blocks; shaded,
huge CIs) vs aquilonia spread across, and the low-recomb polyctena dominance growing
with τ. Writes `data/di25_recomb_directional_<stamp>.rds`, `Figures/di25_recomb_directional_<stamp>.png`.

### `di25_ld_effect_recomb_DI.R` (the LD-clustering effect, one figure)

Amount of sorting toward each species vs recombination decile, split by DI bin,
SNP-based vs cluster-based. `facet_grid(DI bin ~ level)`; two species lines per
panel (aqu purple / pol yellow); point size ~ n loci; lines only where n≥15.
Shows the SNP-based low-recombination inflation (peak ~1.5 cM/Mb) that the eMLG
level removes — and that the low-recomb deciles hold almost no independent units.
Direction gradient down the DI rows. Writes `data/di25_ld_effect_recomb_DI.rds`,
`Figures/di25_ld_effect_recomb_DI.png`.

### `di25_sielva_pruning.R` (early allele pruning in Sielva)

Tests whether Sielva (F1-like colony) has already lost heterozygosity at sorted
loci. Bins units by sorting strength (fraction of the other 19 pops near-fixed)
and shows Sielva's genotype composition. **Result:** het falls ~83%→~52% as loci
sort, and the homozygous fraction is almost entirely TOWARD the fixed allele
(121:1 at strongly-sorted loci) — directional early pruning. Survives recombination
control (logistic sortedness −2.26, p≈1e-72; recomb −0.21). Exception: the 3 huge
blocks (95–100% bin) are still held as intact F1 het (pruning needs recombination).
Caveats: Sielva = 1 colony; selection vs shared-local-ancestry not separable.
Writes `data/di25_sielva_pruning.rds`, `Figures/di25_sielva_pruning.png`.

### `di25_region_segregation.R` (do the blocks segregate within populations?)

Tests two co-adaptation predictions on the Chr 5/25/26 blocks: (1) no within-population
segregation (recombinants inviable) → within-population dosage SD ≈ 0 for every
population (only Nyrhisperä75/Chr5 and Åland/Chr26 mildly segregate, single-block);
(2) Sielva is the F1-like exception (81% het genome-wide vs ~0.36; retains the
heterozygous state at all three blocks). Both hold descriptively but can't
distinguish epistasis from directional sorting (blocks near-fixed → no segregating
variation to test), and Sielva = one colony (`Siel338_*`, n≈1 genotype). Writes
`data/di25_region_segregation.rds`, `Figures/di25_region_segregation.png`.

### `di25_population_fixation.R` (per-population near-fixation)

The population analog of the aggregate sorting sweep: 20 population rings, each
cell = that population's OWN near-fixation direction at the unit (φ=0.85; purple
aqu-fixed / yellow pol-fixed / grey not-fixed), rather than one cross-population
sort class. Exposes repeatability directly — the three *polyctena* blocks
(Chr 5/25/26) are solid yellow across all 20 rings (fixed in every population),
whereas *aquilonia* near-fixation is patchier large blocks + many small units.
Both levels. Writes `Figures/di25_popfix_{snp,emlg}.png`.

### `diem_circos_population.R` (per-population summary)

Collapses individuals to the 20 hybrid populations (Åland most *aquilonia* →
Heinamäki most *polyctena*). Continuous population-mean *aquilonia* frequency is
rendered by binning to a 100-step viridis palette (`freq_to_code()`), so
`render_circos_raster()` handles it like any categorical track. The per-SNP panel
shows strongly-sorted regions as broad consistent bands across all 20 rings
(the large Chr 25/26/7/16/17 LD blocks); the per-eMLG panel collapses those
blocks. `orient_aqu()` aligns each unit so 2 = *aquilonia* via the parent rows.

### `diem_circos_sorting_sweep.R` (τ-sweep circos)

No DIEM genotype data — the rings are the four sorting thresholds, so sorting
visibly **thins out** as τ increases outward, and the direction colours show
where each ancestry prevails. The per-SNP panel's large **yellow (polyctena)
blocks** on Chr 5/25/26 persist to the outer τ=0.8 ring; in the per-eMLG panel
those blocks collapse to single thin units and the residual signal is
predominantly **purple (aquilonia)** — the visual form of the SNP-vs-eMLG
direction flip quantified in `di25_sorting.R`.

### `di25_sorting.R` (sorting, per SNP + per eMLG)

Runs `parallelism_stats()` once per level (τ-independent counts), then
`classify_sort()` at each τ. Individuals: 164 hybrids (20 pops) + 30 parents from
`sample_data` (drops `Hei159_h`, the one tsv hybrid with no population label);
parents supply aquilonia/polyctena orientation. Per-eMLG units use the consensus
for >2-marker clusters and the representative SNP otherwise. Key results at τ=0.6:
per-SNP 12.6% sorted vs per-eMLG 7.2% (LD pseudoreplication removed); and the
direction of the most strongly sorted loci flips from polyctena-leaning at the SNP
level (34% aqu at τ=0.8) to a stable ~80% aquilonia bias at the eMLG level — the
per-SNP polyctena signal is a few large low-recombination blocks counted many times.

### `diem_circos_di25_eMLG.R` (the high-DI-only LD-reduced circos)

Directly comparable to the per-SNP `diem_circos.R` (same DI25 markers, same 195
individuals, same palette), but with the 51,612 SNPs collapsed to the 11,052
high-DI-only LD units. Each unit is its **eMLG consensus** (>2 markers) or its
**representative SNP** (1–2). Consensus is computed over the combined 165
hybrids + 30 parents in one pass (single polarization per unit), then each unit
is oriented to *F. aquilonia* via the parent rows and recoded to the shared
1/2/3 scheme. `CM_STAMP` at the top selects which sweep clustering to plot.

### `di25_ld_clustering.R` (high-DI-only clustering)

The canonical module0 eMLGs are built on all markers, so a cluster anchored at a
high-DI marker also averages in linked lower-DI SNPs. This rebuilds the LD
reduction using **only** the DI25 markers.

* Reads the DI25 tsv, translates to 012 (`1→0, 2→1, 3→2, 0→NA`).
* **Hybrids only** for all LD estimation (parents are near-fixed at diagnostic
  markers and would manufacture artificial LD); parent 012 matrix saved aside.
* `compute_LD_decay` (slide 100, `corr`, keep_el) → `compute_ld_w` (rho 0.95) →
  Stage 1 `ld_complexity_reduction` (rho 0.5) — computed once, cached.
* Stage 2 `ld_prune_and_eMLG`, **no loci limit** entering the merge
  (`ld_w_threshold=0`, `min_n_loci_flag=1`), eMLG stored for clusters with
  **>2 markers** (`min_n_loci_eMLG=3`; 1–2-marker clusters → representative SNP),
  gates `score≥0.80`, `min_r2=0.2`.
* **cM-cap sensitivity sweep** `{0.5, 1, 2, 5, 10}` (fewer SNPs → a tight cap
  over-splits low-density regions). 51,612 markers → 11–13k clusters (75–79%
  reduction); ~4k eMLG units; singletons fall 6,557→4,995 as the cap relaxes,
  saturating by ~5 cM. Outputs one `data/di25_clustering_<stamp>.rds` per cap
  plus `di25_clustering_sweep_summary.rds`.

### `diem_circos_eMLG.R` (LD-reduced)

Maps the DI25 diagnostic markers onto LD-clusters and plots the cluster
**eMLG consensus** instead of individual SNPs:

* Hybrid consensus from `moduleA_sorting/data/eMLG_sorted_cM05.rds` (`$eMLG`,
  164 hybrids × 96,461 clusters) and parent consensus from
  `eMLG_sorted_cM05_parents.rds` (30 parents × the same clusters) → rbind to 194
  individuals.
* A cluster is kept if any member marker is a DI25 diagnostic marker
  (`Chr<n>:<pos>`) → **5,178 clusters**, covering 55.7 % of the DI25 markers
  (the sorted-touched set only carries clusters hit by a *sorted* marker; the
  full canonical partition would be needed for 100 % coverage).
* **Orientation:** the eMLG consensus is oriented to an arbitrary within-cluster
  reference allele, so each cluster is re-oriented to *F. aquilonia* using the
  parent side (mean Faqu dosage ≥ mean Fpol ⇒ keep, else `2 − x`), then rounded
  and recoded to the shared 1/2/3 scheme via `3 − round(a)`.
* Clusters are positioned by their representative marker; individuals ordered by
  hybrid index (inner = most *F. aquilonia*, matching the SNP plot). Parents
  anchor the extremes (Faqu ≈ 0.17, Fpol ≈ 0.73); hybrids sit tightly mid-range.

### `diem_circos.R` (per-SNP)

Input genotypes are coded relative to the two parental species:

| code | meaning | colour (viridis) |
|----|----|----|
| 0 | missing (the "dithered" cells) | dithered → random state; grey if `DITHER=FALSE` |
| 1 | homozygous *F. aquilonia* | purple `#440154` |
| 2 | heterozygous | teal `#21918C` |
| 3 | homozygous *F. polyctena* | yellow `#FDE725` |

(Confirmed against the 15 Faqu parent columns ≈ all 1 and the Fpol parents ≈ all 3.)

**Layout.** One sector per chromosome, angular span ∝ marker count (chr 23 is
absent from the data). Each individual is a concentric ring; individuals are
ordered by genome-wide hybrid index, so the innermost rings are the most
*F. aquilonia* (solid purple core) and the outermost the most *F. polyctena*.

**Rendering.** The genotype rings are drawn by direct polar rasterisation (each
output pixel mapped to chromosome/marker/individual, then coloured), which
renders the full 51,612-marker × 195-individual matrix in ~2 s. This replaces
`circos.heatmap()`, which draws every cell as a vector polygon (~10 M polygons)
and does not scale here. Additional annular tracks (DI, sorting call, …) can be
layered on the same `[-1,1]` canvas.

Key parameters at the top of the script: `DITHER`, palette `COL`, canvas `NPX`,
ring radii `R_IN`/`R_OUT`, `OPEN_DEG`/`CHR_GAP_DEG`.

## Inputs

`data/species_diagnostic_markers_DI25_20pops.tsv.gz` — 195 individual genotype
columns + `chromosome` + `position`; species-diagnostic markers at DI > −25
across the 20 hybrid populations (plus the parental reference individuals).
