# Module A — are the per-SNP "bidirectional" manhattan peaks real?

**Script:** `moduleA_sorting/R/moduleA_bidirectional_peaks.R`
**Figure:** `moduleA_sorting/Figures/moduleA_bidirectional_peaks.png`
**Table:** `moduleA_sorting/data/moduleA_bidirectional_peaks.csv`

## Question

The binom sorting manhattan (`moduleA_fig_sorting_manhattan.R`) shows local peaks
in the **per-SNP bidirectional fraction at τ = 0.5**, concentrated on Chr1, Chr8,
Chr17, Chr22, Chr24, Chr26. Several coincide with high‑`ld_w` regions, so they
could be single large LD clusters counted once per SNP rather than genomic regions
of genuinely heterogeneous ancestry outcomes. We test four things (bidirectional
defined at τ = 0.5: `prop_fixed ≥ 0.5 & p_binom ≥ 0.05 & n_fixed ≥ 6`).

## Findings

### 1. Per SNP vs per LD‑reduced unit — the peaks collapse

Across all peaks, **1,946 per‑SNP bidirectional calls reduce to 192 distinct LD
clusters, and only 61 bidirectional cluster representatives.** The prominent peaks
are each essentially *one* high‑`ld_w` cluster whose representative is not even
bidirectional:

| peak | n_bi SNP | LD clusters | bi reps | LD‑inflation | mean `ld_w` |
|---|--:|--:|--:|--:|--:|
| Chr17:1.7–2.4 | 270 | **1** | 0 | 270× | 0.32 |
| Chr24:1.3–1.7 | 209 | 2 | 1 | 105× | 0.74 |
| Chr26:4.9–5.8 | 162 | 2 | 0 | 81× | 0.83 |
| Chr1:5.2–5.8 | 111 | **1** | 0 | 111× | 0.62 |
| Chr8:6.1–6.4 | 109 | 6 | 1 | 18× | 0.63 |
| Chr22:2.0–2.2 | 93 | **1** | 0 | 93× | 0.83 |
| Chr12:6.2–6.3 | 38 | 16 | 4 | 2.4× | 0.06 |
| Chr12:6.4–6.5 | 38 | 16 | 7 | 2.4× | 0.07 |

Only the **low‑`ld_w` peaks (Chr3, Chr4, Chr12)** span several independent units;
the visually striking high‑`ld_w` peaks are single blocks.

### 2. `n_fixed` inside the peaks

**1,641 of 1,946 (84%) peak‑bidirectional SNPs have `n_fixed` = 10 or 11**
(`prop_fixed` 0.50–0.55). The primary τ = 0.6 gate needs `n_fixed ≥ 12`, so almost
all of them drop out between τ = 0.5 and τ = 0.6 — this is exactly why the peaks
are blue (τ = 0.5) only in the manhattan.

### 3. Peak intervals

~70 enriched 100‑kb windows, dominated by **Chr1:2.4–5.8, Chr8:6.1–6.7,
Chr17:1.7–2.4 & 5.1, Chr22:2.0–2.2, Chr24:1.3–2.1, Chr26:4.9–5.8** (single/few‑cluster,
high‑`ld_w`), plus scattered low‑`ld_w` windows on Chr3/4/12. Full list in the CSV.

### 4. Population groupings

- **Within a peak the split is coherent** (mean coherence ≈ **1.00**): because each
  peak is one LD block, every population sorts the *same* way across all its SNPs
  (clean `52/0` / `0/47` counts). A "peak" is really one locus where roughly half
  the populations fixed aquilonia and half fixed polyctena.
- **Different peaks partition the populations differently.** Adjacent same‑block
  peaks share a split, but across regions the population net‑direction vectors are
  weakly or negatively correlated (e.g. Chr22 vs Chr24 ≈ −0.19; Chr24 vs Chr6 ≈
  −0.47). There is **no global population structure** driving the peaks — see the
  figure: each column is a clean but region‑specific partition.

## Conclusion

The peaks **collapse under LD reduction**: at the eMLG/consensus level bidirectional
is already 0 (`moduleA_cluster_sorting.rds`), and the per‑SNP peaks are LD‑inflated
single blocks of marginally sorted loci (`n_fixed` 10–11) with region‑specific
population splits. They are therefore a **demonstration of why the primary τ = 0.6
threshold and LD‑aware (per‑unit) summaries are necessary**, not a set of genuinely
heterogeneous‑ancestry regions. The residue of genuine both‑direction units (~61
genome‑wide, mostly the low‑`ld_w` Chr3/4/12 windows) is marginal and would require
Module E's neutral null before it could be interpreted.

**Figure — `moduleA_bidirectional_peaks.png`:** hybrid population (rows) × the ten
largest bidirectional peaks (columns); tile colour = mean per‑locus fixation
direction (blue = aquilonia, orange = polyctena), label = #loci aquilonia‑fixed /
polyctena‑fixed. Each column is a clean partition; the partitions differ between
peaks.
