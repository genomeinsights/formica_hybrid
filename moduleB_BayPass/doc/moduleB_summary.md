# Module B — Climate association: summary for the manuscript

**Scope.** Module B (`moduleB_BayPass/`; the cleaned re-organisation of the old `moduleC_*`) tests
whether population allele frequencies are associated with the climate covariates
**PC1** and **PC2**, at two resolutions: (i) individual SNPs (standard BayPass
outlier scan) and (ii) LD-complexity-reduced clusters (eMLGs; one consensus
genotype per LD cluster, *sensu* Li et al. 2018). It then calibrates the
cluster-level signal against a population-structure null to obtain a valid
genome-wide FDR.

**Primary configuration throughout:** Åland-excluded (19 populations, 154
individuals), fixed LD-pruned Ω covariance (`withOmega`); marker-level
significance = **BF(dB) ≥ 15**. Full (20-population) and internally-estimated-Ω
runs are sensitivity comparisons. Data: 1,114,340 MAF-filtered SNPs → 32,840 LD
clusters with an eMLG (≥5 loci; `ld_w=0.025`, `cM=0.5`); genome-wide median
cluster size = 9 markers.

---

## 1. Main-text conclusions (short — headline only)

1. **Standard single-SNP climate-outlier scans over-call large, differentiated
   regions.** Under the size-dependent member-count criterion (a cluster is a
   "candidate" when ≥5 member SNPs reach BF ≥ 15), candidate clusters are 7–350×
   larger than the genome-wide median and sit near the background significant-marker
   rate — i.e. differentiated blocks, not focal climate peaks.
2. **LD-complexity reduction increases power by collapsing ~1.1 M SNP tests to
   32,840 cluster tests (~34× lower multiple-testing burden).** Testing the eMLG
   (cluster consensus) recovers coherent cluster-level climate signals that the
   SNP scan cannot: **28 "eMLG-only" clusters** reach eMLG BF ≥ 15 with *no*
   individually significant member SNP.
3. **After a population-structure null and a valid genome-wide FDR, the climate
   signal is largely a structure/differentiation by-product.** Against 10,000
   Ω-structured null covariates, **none of the SNP-based candidates survive**
   (the flagship Chr24 block F63637, BF = 38.4, is beaten by 12/10,000 structure
   nulls); only **5 PC2 "eMLG-only" clusters** resist structure (set FDR ≈ 0.66 —
   suggestive, not significant), and **PC1 shows nothing**.
4. **Net:** no locus reaches a defensible FDR for climate association after
   structure control. This is a clean negative that both (a) corroborates the
   weak climate link seen elsewhere and (b) is itself a demonstration of the value
   of LD-complexity reduction — single-SNP outliers behave as expected false
   positives, and the only residual signal is the reduced-test, aggregated eMLG
   signal.

---

## 2. Methods (supplementary by default)

### 2.1 SNP-level BayPass climate association
Association of per-population allele frequencies with PC1/PC2 was evaluated with
the standard covariate (IS) model in BayPass v3 (`-nocovscaling -nval 500 -burnin
5000 -thin 25`), evidence reported as Bayes factors on the deciban scale, BF(dB).
Eight configurations crossed {full, Åland-excluded} × {PC1, PC2} × {Ω estimated
internally from the full set (`noOmega`), Ω pre-estimated from the LD-pruned set
and supplied fixed (`withOmega`)}. Output rows were matched positionally to the
input markers (row-count and strictly-increasing-index asserts). Åland-excluded ×
withOmega is primary. Background significant-marker rate (BF ≥ 15): **PC1 0.17 %,
PC2 0.15 %**.
Scripts: `moduleB_BayPass/R/moduleB_prepare_{with_aland,aland_excluded}.R`,
`moduleB_write_baypass_inputs.R`, `{with_aland,aland_excluded}/run_baypass.sh`,
`moduleB_analyse_baypass.R`.

### 2.2 PC provenance and the ancestry confound
PC1/PC2 are per-population **climate** covariates (within-population SD = 0, i.e.
not a per-individual genetic PCA). Rank correlation with per-population *aquilonia*
ancestry (≈49k diagnostic markers): **PC1 = −0.229**, **PC2 = +0.568**. The PC2–
ancestry link is driven by **Åland** (+0.568 → +0.384 when dropped), not Sielva
(+0.568 → +0.642 when dropped). Controls: fixed Ω (conditions out shared structure,
present in both runs) and the Åland-excluded primary analysis. Associations are
therefore interpreted as **conditional** associations, not direct evidence of local
adaptation. Script: `moduleB_BayPass/R/moduleB_ancestry_confound.R`.

### 2.3 LD-complexity reduction and the eMLG association (Li et al. 2018)
Each ≥5-locus LD cluster was reduced to one **eMLG** (per-individual consensus
dosage, rounded to {0,1,2}), aggregated to per-population diploid allele counts,
and run through the *same* BayPass covariate model with the *same* fixed Ω,
covariates and pool sizes (19-population order asserted identical). This tests
whether a cluster tracks climate **as a unit** rather than via a handful of
individually significant members. Because Ω is a property of the populations and
BayPass co-estimates a genome-wide allele-frequency prior from *all* markers in a
run, this must be run on the full 32,840-eMLG set (verified: subsetting to
candidates shifts BF by up to tens of dB). Scripts:
`moduleB_BayPass/R/moduleB_write_eMLG_baypass_inputs.R`, `run_baypass_eMLG.sh`,
`moduleB_eMLG_manhattan.R`. (Note: BayPass `eBPis` is degenerate here — essentially
all eMLGs p < 0.05 — so BF(dB), calibrated by the sim-null below, is used, not eBPis.)

### 2.4 Ω-structured simulation null and genome-wide sim-FDR
Population structure inflates BF for differentiated clusters (structure-only null
covariates reach eMLG BF up to 44). To calibrate, **10,000 null "climate"
covariates** were drawn from the fixed Ω (structure autocorrelation, no genotype
link) and run through the eMLG covariate model in 10 batches of 1,000 on the full
set (validated equivalent to a single 10,000-covariate run and to separate runs;
marker-level BF depends on the marker set but not the covariate count). For each
eMLG we count `k` = number of the 10,000 structure nulls with BF ≥ the observed
climate BF (`k1` for PC1, `k2` for PC2). A **floor-survivor** (`k = 0`) is a cluster
whose real climate BF beats *all* 10,000 structure nulls. This is a valid
permutation-style genome-wide FDR: under the pure-structure null the observed BF is
exchangeable with its nulls, so the expected number of floor-survivors is
N/(NSIM+1) = 32,840/10,001 ≈ **3.28** per axis, and the FDR of the observed
floor-survivor set = 3.28 / (number observed). Any cluster above the floor (k ≥ 1)
cannot survive FDR and is discarded. Script: `moduleB_BayPass/R/moduleB_eMLG_null.R` →
`moduleB_BayPass/data/moduleB_eMLG_null.rds`. (Scoping FDR to a pre-screened subset is invalid —
selection on the same statistic — so the genome-wide denominator is used.)

---

## 3. Figures (primary analyses)

![**SNP-level Manhattan, PC1** (primary config). Grey = eMLG-filtered background; coloured = clusters with ≥5 member SNPs at BF(dB) ≥ 15 (the member-count candidate set), coloured by LD cluster; dashed line BF = 15.](../Figures/manhattan_PC1_aland_excluded_withOmega_eMLG_clustered_min5SigLoci_withBackground.png){width=100%}

![**SNP-level Manhattan, PC2** (primary config). Colouring as above.](../Figures/manhattan_PC2_aland_excluded_withOmega_eMLG_clustered_min5SigLoci_withBackground.png){width=100%}

![**eMLG-level Manhattan (complexity-reduced)**, primary config. One point per eMLG at its representative position. Red = eMLG BF ≥ 15 with 0 significant member SNPs (eMLG-only); blue = 1–4 sig SNPs; black = ≥5 (member-count candidates); grey = ns; dashed line BF = 15. **Larger triangles = the sim-FDR floor-survivors** (real eMLG BF exceeds all 10,000 Ω-structured population-structure nulls): F7388 on PC1; F14310, F49926, F53948, F54421, F58676 on PC2.](../Figures/moduleB_eMLG_manhattan.png){width=100%}

![**Final Manhattan — FDR-passing eMLG clusters over all SNPs.** All 1,114,340 SNPs in grey (BF per axis); the member SNPs of the sim-FDR floor-surviving eMLG clusters are coloured, one colour per cluster (PC1: F7388; PC2: F14310, F49926, F53948, F54421, F58676). Most member SNPs lie *below* BF = 15 — the cluster-level signal is not driven by individually significant SNPs. Dashed line BF = 15.](../Figures/moduleB_fdr_snp_manhattan.png){width=100%}

![**PC–ancestry confound.** Per-population climate PC vs genome-wide *aquilonia* ancestry; Åland/Sielva highlighted.](../Figures/moduleC_ancestry_confound.png){width=68%}

Compiled candidate tables (standalone PDFs): `tables/climate_candidate_clusters.pdf`,
`tables/climate_candidate_full.pdf`.

---

## 4. Tables (complete)

### Table S1 — Putative candidate counts across the 8 BayPass configurations
Member-count criterion only (a cluster is a candidate when ≥5, or ≥10, member SNPs
reach BF ≥ 15). Genome-wide median cluster size = 9 markers.

| Population set | Ω (covariance) | Axis | Candidates (≥5) | Candidates (≥10) | Median candidate size |
|---|---|---|---:|---:|---:|
| Åland-excluded (19) | Internal (full data) | PC1 | 14 | 7 | 68 |
| Åland-excluded (19) | Internal (full data) | PC2 | 16 | 7 | 116 |
| **Åland-excluded (19)** | **LD-pruned (fixed) [PRIMARY]** | **PC1** | **16** | **10** | **68** |
| **Åland-excluded (19)** | **LD-pruned (fixed) [PRIMARY]** | **PC2** | **21** | **8** | **118** |
| Full (20) | Internal (full data) | PC1 | 9 | 5 | 272 |
| Full (20) | Internal (full data) | PC2 | 43 | 24 | 94 |
| Full (20) | LD-pruned (fixed) | PC1 | 7 | 5 | 486 |
| Full (20) | LD-pruned (fixed) | PC2 | 53 | 36 | 89 |

### Table S2 — Member-count candidates, primary config (35 clusters)
Entry = ≥5 member SNPs at BF ≥ 15 on PC1 or PC2. Columns: member SNP count/% on the
entry axis; eMLG BF(dB) and β (cluster-consensus association); XtX (among-population
differentiation); DI% (proportion of members with DiagnosticIndex > −25); `simnull_k`
= # of 10,000 structure nulls ≥ observed eMLG BF on the entry axis; `gate15` = eMLG
BF ≥ 15; `simFDR` = floor-survivor (k = 0). Ordered by eMLG BF.

| cluster | Chr | pos (Mb) | size | axis | SNP nsig | SNP % | eMLG BF | eMLG β | XtX | DI% | simnull_k | gate15 | simFDR |
|---|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|:--:|:--:|
| F63637 | Chr24 | 1.3–1.8 | 272 | PC1 | 205 | 75.4 | 38.4 | 0.094 | 37 | 59 | 12 | ✓ | — |
| F63772 | Chr24 | 1.8–1.8 | 83 | PC2 | 60 | 72.3 | 36.5 | 0.122 | 26 | 0 | 1 | ✓ | — |
| F17762 | Chr5 | 3.0–3.3 | 11 | PC2 | 7 | 63.6 | 30.7 | 0.123 | 26 | 27 | 1 | ✓ | — |
| F41962 | Chr13 | 5.8–6.2 | 50 | PC2 | 7 | 14.0 | 21.8 | 0.095 | 26 | 0 | 8 | ✓ | — |
| F64355 | Chr24 | 1.9–5.4 | 35 | PC1 | 7 | 20.0 | 21.2 | 0.051 | 23 | 6 | 9 | ✓ | — |
| F6554 | Chr2 | 10.0–10.1 | 24 | PC1 | 12 | 50.0 | 19.7 | 0.067 | 20 | 0 | 9 | ✓ | — |
| F63812 | Chr24 | 2.0–2.3 | 23 | PC1 | 10 | 43.5 | 18.9 | 0.057 | 23 | 22 | 18 | ✓ | — |
| F27646 | Chr8 | 6.2–6.8 | 486 | PC1 | 108 | 22.2 | 14.9 | 0.058 | 23 | 46 | 119 | — | — |
| F67035 | Chr26 | 1.6–4.7 | 397 | PC2 | 19 | 4.8 | 14.4 | 0.076 | 23 | 73 | 6 | — | — |
| F49451 | Chr16 | 0.4–5.9 | 522 | PC1 | 69 | 13.2 | 13.7 | 0.055 | 21 | 0 | 16 | — | — |
| F17747 | Chr5 | 2.8–4.3 | 690 | both | 24 | 3.5 | 11.5 | 0.088 | 23 | 14 | 296 | — | — |
| F11460 | Chr3 | 6.2–8.5 | 1143 | PC2 | 25 | 2.2 | 10.5 | 0.079 | 21 | 1 | 45 | — | — |
| F66997 | Chr26 | 1.4–4.7 | 3042 | PC2 | 37 | 1.2 | 10.1 | 0.082 | 23 | 17 | 231 | — | — |
| F48626 | Chr15 | 8.2–8.2 | 19 | PC1 | 5 | 26.3 | 8.2 | 0.050 | 25 | 0 | 401 | — | — |
| F14815 | Chr4 | 6.4–11.0 | 99 | PC2 | 19 | 19.2 | 7.6 | 0.078 | 30 | 28 | 557 | — | — |
| F27668 | Chr8 | 1.0–6.9 | 400 | PC2 | 12 | 3.0 | 7.2 | 0.061 | 24 | 0 | 568 | — | — |
| F64102 | Chr24 | 3.1–4.4 | 28 | both | 7 | 25.0 | 5.7 | 0.074 | 26 | 4 | 613 | — | — |
| F67238 | Chr26 | 6.6–7.5 | 42 | PC2 | 7 | 16.7 | 5.0 | 0.062 | 21 | 12 | 453 | — | — |
| F51523 | Chr17 | 0.0–4.8 | 3183 | PC2 | 24 | 0.8 | 4.3 | 0.058 | 19 | 0 | 561 | — | — |
| F57265 | Chr20 | 0.4–4.4 | 87 | PC2 | 5 | 5.7 | 3.7 | 0.063 | 21 | 0 | 985 | — | — |
| F40144 | Chr12 | 6.3–10.4 | 22 | PC1 | 5 | 22.7 | 1.6 | 0.039 | 20 | 5 | 746 | — | — |
| F13363 | Chr3 | 13.4–15.7 | 131 | PC2 | 6 | 4.6 | 0.7 | 0.048 | 23 | 1 | 2012 | — | — |
| F67071 | Chr26 | 1.5–5.9 | 379 | PC1 | 13 | 3.4 | −2.5 | 0.030 | 23 | 20 | 3400 | — | — |
| F51449 | Chr17 | 5.1–5.7 | 893 | PC1 | 6 | 0.7 | −4.0 | 0.028 | 20 | 5 | 4243 | — | — |
| F64194 | Chr24 | 1.9–5.1 | 182 | PC1 | 14 | 7.7 | −6.4 | 0.024 | 26 | 8 | 6390 | — | — |
| F68272 | Chr27 | 0.1–3.2 | 496 | PC2 | 8 | 1.6 | −10.9 | 0.018 | 19 | 12 | 9991 | — | — |
| F24440 | Chr6 | 10.1–10.2 | 82 | PC2 | 5 | 6.1 | −11.2 | 0.018 | 23 | 33 | 9989 | — | — |
| F67028 | Chr26 | 1.8–6.7 | 1236 | PC2 | 5 | 0.4 | −11.9 | 0.016 | 22 | 41 | 9997 | — | — |
| F67081 | Chr26 | 6.3–7.0 | 449 | PC2 | 5 | 1.1 | −12.6 | 0.014 | 22 | 51 | 10000 | — | — |
| F65755 | Chr25 | 12.7–12.7 | 20 | PC2 | 6 | 30.0 | −13.4 | 0.012 | 21 | 5 | 10000 | — | — |
| F67249 | Chr26 | 7.4–7.5 | 71 | PC2 | 26 | 36.6 | −14.1 | 0.011 | 22 | 62 | 10000 | — | — |
| F54612 | Chr18 | 3.0–3.1 | 118 | PC2 | 8 | 6.8 | −14.5 | 0.009 | 22 | 42 | 10000 | — | — |
| F7206 | Chr2 | 9.7–10.1 | 64 | PC1 | 28 | 43.8 | −15.6 | 0.007 | 26 | 2 | 10000 | — | — |
| F7494 | Chr2 | 9.9–10.2 | 71 | PC1 | 22 | 31.0 | −16.4 | 0.006 | 23 | 0 | 10000 | — | — |
| F6596 | Chr2 | 10.0–10.1 | 14 | PC1 | 5 | 35.7 | −17.2 | 0.005 | 23 | 0 | 10000 | — | — |

*Note the demotion pattern:* several high member-% clusters (F6596, F7494, F7206,
F67249, F65755) have strongly **negative** eMLG BF — a handful of member SNPs cross
the threshold but the cluster consensus anti-associates. **None of the 35 survive
the sim-FDR** (`simFDR` all —).

### Table S3 — eMLG outliers (BF ≥ 15), primary config (51 clusters)
The complexity-reduced view. `support`: eMLG-only = 0 significant member SNPs;
1–4 SNPs; ≥5 SNPs = also a member-count candidate. `simnull_k` = # of 10,000
structure nulls ≥ observed; **floor** = k = 0 (survives the structure null).

| cluster | Chr | pos (bp) | size | axis | SNP nsig | SNP % | eMLG BF | support | simnull_k | floor |
|---|---|---:|---:|---|---:|---:|---:|---|---:|:--:|
| F63637 | Chr24 | 1673178 | 272 | PC1 | 205 | 75.4 | 38.4 | ≥5 SNPs(cand) | 12 | — |
| F7388 | Chr2 | 10138607 | 20 | PC1 | 1 | 5.0 | 31.0 | 1–4 SNPs | 0 | **✓** |
| F43576 | Chr13 | 8123831 | 6 | PC1 | 2 | 33.3 | 24.6 | 1–4 SNPs | 5 | — |
| F7404 | Chr2 | 12899649 | 28 | PC1 | 2 | 7.1 | 21.6 | 1–4 SNPs | 12 | — |
| F64355 | Chr24 | 4784788 | 35 | PC1 | 7 | 20.0 | 21.2 | ≥5 SNPs(cand) | 9 | — |
| F6554 | Chr2 | 10056976 | 24 | PC1 | 12 | 50.0 | 19.7 | ≥5 SNPs(cand) | 9 | — |
| F63812 | Chr24 | 2032351 | 23 | PC1 | 10 | 43.5 | 18.9 | ≥5 SNPs(cand) | 18 | — |
| F31722 | Chr9 | 7966699 | 5 | PC1 | 2 | 40.0 | 18.2 | 1–4 SNPs | 11 | — |
| F45693 | Chr14 | 5672151 | 10 | PC1 | 0 | 0.0 | 16.6 | eMLG-only | 34 | — |
| F41301 | Chr13 | 4478466 | 12 | PC1 | 0 | 0.0 | 16.3 | eMLG-only | 48 | — |
| F7080 | Chr2 | 9189106 | 41 | PC1 | 0 | 0.0 | 16.1 | eMLG-only | 97 | — |
| F23332 | Chr6 | 7613474 | 8 | PC1 | 0 | 0.0 | 16.0 | eMLG-only | 75 | — |
| F61797 | Chr22 | 2726658 | 42 | PC1 | 0 | 0.0 | 15.4 | eMLG-only | 64 | — |
| F20510 | Chr5 | 8946084 | 13 | PC1 | 0 | 0.0 | 15.3 | eMLG-only | 110 | — |
| F63772 | Chr24 | 1777883 | 83 | PC2 | 60 | 72.3 | 36.5 | ≥5 SNPs(cand) | 1 | — |
| F54421 | Chr18 | 4661936 | 7 | PC2 | 0 | 0.0 | 35.9 | eMLG-only | 0 | **✓** |
| F53948 | Chr18 | 6499523 | 11 | PC2 | 0 | 0.0 | 34.2 | eMLG-only | 0 | **✓** |
| F14310 | Chr4 | 9075597 | 24 | PC2 | 0 | 0.0 | 31.7 | eMLG-only | 0 | **✓** |
| F17762 | Chr5 | 3040193 | 11 | PC2 | 7 | 63.6 | 30.7 | ≥5 SNPs(cand) | 1 | — |
| F21890 | Chr6 | 3096997 | 5 | PC2 | 2 | 40.0 | 27.3 | 1–4 SNPs | 6 | — |
| F58676 | Chr20 | 8129042 | 5 | PC2 | 0 | 0.0 | 23.9 | eMLG-only | 0 | **✓** |
| F12448 | Chr3 | 10504842 | 6 | PC2 | 2 | 33.3 | 23.3 | 1–4 SNPs | 8 | — |
| F49926 | Chr16 | 7271131 | 17 | PC2 | 0 | 0.0 | 22.6 | eMLG-only | 0 | **✓** |
| F35840 | Chr11 | 6084388 | 14 | PC2 | 0 | 0.0 | 22.0 | eMLG-only | 1 | — |
| F41962 | Chr13 | 5926761 | 50 | PC2 | 7 | 14.0 | 21.8 | ≥5 SNPs(cand) | 8 | — |
| F32802 | Chr9 | 9005754 | 9 | PC2 | 4 | 44.4 | 21.5 | 1–4 SNPs | 14 | — |
| F3838 | Chr1 | 14393778 | 12 | PC2 | 0 | 0.0 | 20.7 | eMLG-only | 3 | — |
| F6754 | Chr2 | 11986955 | 10 | PC2 | 0 | 0.0 | 20.1 | eMLG-only | 1 | — |
| F60466 | Chr21 | 9159798 | 5 | PC2 | 1 | 20.0 | 20.0 | 1–4 SNPs | 13 | — |
| F13395 | Chr3 | 14731412 | 9 | PC2 | 4 | 44.4 | 19.7 | 1–4 SNPs | 3 | — |
| F12565 | Chr3 | 10321642 | 6 | PC2 | 1 | 16.7 | 19.5 | 1–4 SNPs | 7 | — |
| F16915 | Chr4 | 13266783 | 15 | PC2 | 0 | 0.0 | 19.3 | eMLG-only | 1 | — |
| F41903 | Chr13 | 6175093 | 7 | PC2 | 1 | 14.3 | 19.0 | 1–4 SNPs | 34 | — |
| F31253 | Chr9 | 10491866 | 11 | PC2 | 0 | 0.0 | 18.6 | eMLG-only | 6 | — |
| F53378 | Chr18 | 6514541 | 5 | PC2 | 0 | 0.0 | 18.3 | eMLG-only | 2 | — |
| F40827 | Chr13 | 2402018 | 57 | PC2 | 1 | 1.8 | 18.1 | 1–4 SNPs | 14 | — |
| F35768 | Chr11 | 6010746 | 6 | PC2 | 0 | 0.0 | 18.1 | eMLG-only | 13 | — |
| F36504 | Chr11 | 8883438 | 9 | PC2 | 0 | 0.0 | 17.8 | eMLG-only | 2 | — |
| F30134 | Chr9 | 3458646 | 10 | PC2 | 1 | 10.0 | 17.5 | 1–4 SNPs | 2 | — |
| F13266 | Chr3 | 14003903 | 16 | PC2 | 0 | 0.0 | 17.4 | eMLG-only | 22 | — |
| F11525 | Chr3 | 5371732 | 5 | PC2 | 0 | 0.0 | 17.0 | eMLG-only | 2 | — |
| F61514 | Chr21 | 9665372 | 16 | PC2 | 0 | 0.0 | 17.0 | eMLG-only | 36 | — |
| F57853 | Chr20 | 6463479 | 6 | PC2 | 1 | 16.7 | 16.5 | 1–4 SNPs | 61 | — |
| F48934 | Chr15 | 8997189 | 6 | PC2 | 0 | 0.0 | 16.0 | eMLG-only | 27 | — |
| F34281 | Chr10 | 10908690 | 8 | PC2 | 0 | 0.0 | 16.0 | eMLG-only | 33 | — |
| F34785 | Chr10 | 11485945 | 8 | PC2 | 0 | 0.0 | 15.6 | eMLG-only | 43 | — |
| F26091 | Chr7 | 12741603 | 5 | PC2 | 0 | 0.0 | 15.6 | eMLG-only | 13 | — |
| F32290 | Chr9 | 6756230 | 8 | PC2 | 1 | 12.5 | 15.1 | 1–4 SNPs | 11 | — |
| F39601 | Chr12 | 8679100 | 12 | PC2 | 0 | 0.0 | 15.1 | eMLG-only | 16 | — |
| F48676 | Chr15 | 8634509 | 5 | PC2 | 0 | 0.0 | 15.1 | eMLG-only | 22 | — |
| F52036 | Chr17 | 6835729 | 8 | PC2 | 2 | 25.0 | 15.0 | 1–4 SNPs | 54 | — |

Composition of the 51 outliers: PC1 = 14 (4 member-count candidates, 4 with 1–4
sig SNPs, 6 eMLG-only); PC2 = 37 (3 candidates, 12 with 1–4 sig SNPs, 22 eMLG-only).

### Table S4 — Ω-structured sim-FDR (10,000 nulls), primary config
Floor-survivor = eMLG whose real climate BF exceeds all 10,000 population-structure
nulls (k = 0). Expected floor count under the pure-structure null = 3.28/axis;
set FDR = 3.28 / (observed).

| axis | floor-survivors | expected (null) | set FDR | survivors |
|---|---:|---:|---:|---|
| PC1 | 1 | 3.28 | ≈ 3.28 (→ 1) | F7388 |
| PC2 | 5 | 3.28 | **0.66** | F54421, F53948, F14310, F58676, F49926 |

All five PC2 survivors are **eMLG-only** clusters (no individually significant member
SNP), max eMLG BF 35.9 (F54421). None of the 7 SNP-based member-count candidates
survive (Table S2, `simnull_k` column: F63637 k = 12, F63812 k = 18, …).

---

## 5. Supplementary text / interpretation notes

- **Size dependence of the member-count criterion.** Median candidate cluster size
  (68–486 markers) vastly exceeds the genome-wide median (9), because an *absolute*
  count of significant members scales with cluster size. This is why the member-count
  set is treated as descriptive/putative and the quantitative inference uses the
  eMLG.
- **Why the sim-null, and why it must be genome-wide.** The eMLG set is enriched for
  strong-LD / differentiated blocks, exactly where climate PCs (partly confounded
  with ancestry, §2.2) leak through Ω. Structure-only null covariates reach eMLG BF
  up to 44, so BF ≥ 15 alone is not self-sufficient. FDR must use the full 32,840-eMLG
  denominator; scoping it to the pre-screened outliers is selection on the test
  statistic and is invalid (equivalently, it is the same error as subsetting the
  BayPass input, which also changes the co-estimated allele-frequency prior and the
  BF).
- **Resolution.** At 1,000 sims the empirical-p floor (1/1001) equals the null
  expectation (32.8 ≈ observed 32) → FDR ≈ 1 (nothing). 10,000 sims lower the null
  expectation to 3.28, so `set FDR = 3.28/m` and q < 1 becomes attainable (PC2 = 0.66).
  Reaching q < 0.05 would require ≥ 66 floor-survivors, which the data do not deliver.
- **Interpretive framing (positive case for LD-complexity reduction).** The
  single-SNP/member-count candidates failing while the reduced-test aggregated eMLG
  signal is the only residual is *expected*: outlier scans are prone to
  differentiation-driven false positives. The eMLG reduction is designed for higher
  power (≈34× fewer tests) and recovers coherent cluster-level signals invisible to
  the SNP scan (the eMLG-only class). The residual PC2 eMLG-only set (FDR ≈ 0.66) is
  suggestive only.

---

## 6. Reproducibility (scripts & data objects)

Self-contained module folder `moduleB_BayPass/` with `R/` (scripts), `data/`
(output objects), `Figures/` (plots) and `doc/` (this summary + `doc/tables/`
candidate tables). Scripts are run **from the repo root** in the order below;
they read shared/raw inputs from the repo root — the LD clustering
`module0_ld_pruning/data/eMLG_5loci_0025_cM05.rds` (a Module 0 output), the
genotype `data/hybrids_*_maf005.Rdata` and `data/pruned_markers.rds`, and
the BayPass run directories `with_aland/`, `aland_excluded/`, `aland_excluded_eMLG/`
(which hold the per-run inputs, `run_baypass*.sh`, and the large `*_summary_*.out`)
— and write all module products under `moduleB_BayPass/`.

| # | Step | Script (in `moduleB_BayPass/R/`) | Key output |
|---|---|---|---|
| 1 | BayPass inputs | `moduleB_prepare_{with_aland,aland_excluded}.R` → `moduleB_write_baypass_inputs.R` | `{with_aland,aland_excluded}/` BayPass input files |
| — | BayPass runs (external) | `{with_aland,aland_excluded}/run_baypass.sh` | `…/PC{1,2}_DIEM_*_summary_*.out` |
| 2 | SNP Manhattans | `moduleB_analyse_baypass.R` | `Figures/manhattan_*_min5SigLoci_*.png` |
| 3 | PC/ancestry confound | `moduleB_ancestry_confound.R` | `data/moduleC_ancestry_confound.{rds,png}` |
| 4 | Candidate definition + tables | `moduleB_climate_candidates.R` | `data/moduleB_climate_candidates.rds`; `doc/tables/climate_candidate_*.{tex,pdf}` |
| 5 | eMLG BayPass inputs | `moduleB_write_eMLG_baypass_inputs.R` | `aland_excluded_eMLG/u_eMLG.geno` (+ staged Ω/covariates) |
| — | eMLG BayPass runs (external) | `run_baypass_eMLG.sh` | `aland_excluded_eMLG/PC*_eMLG_*_betai_reg.out` |
| 6 | eMLG Manhattan + outliers | `moduleB_eMLG_manhattan.R` | `data/moduleB_eMLG_association.rds`; `data/moduleB_eMLG_outliers.csv`; `Figures/moduleB_eMLG_manhattan.png` |
| 7 | Ω-structured sim-FDR (10k) | `moduleB_eMLG_null.R` | `data/moduleB_eMLG_null.rds` |
| 8 | FDR-overlay figures | re-run `moduleB_eMLG_manhattan.R` (reads step 7) | `Figures/moduleB_eMLG_manhattan.png` (triangles) + `Figures/moduleB_fdr_snp_manhattan.png` |

Step 8: `moduleB_eMLG_manhattan.R` picks up `data/moduleB_eMLG_null.rds` if present,
so re-running it after step 7 adds the FDR floor-survivor overlays (Figs 3–4). This
summary PDF is rebuilt from `doc/` with `doc/render_moduleB_summary.sh`.

**Reference for the LD-cluster / eMLG association approach:** Li Z, Kemppainen P,
Rastas P, Merilä J. 2018. Linkage disequilibrium clustering-based approach for
association mapping with tightly linked genomewide data. *Mol Ecol Resour* 18:809–824.
