# Founder AIM frequencies: indexing bug in `AIMs_for_SLiM.R`

Repository: github.com/bcportinha/Replicate-hybrid-evolution, `SLiM/R_scripts/AIMs_for_SLiM.R`
(commit c0322fc). Found while comparing the `diem_boot*` outputs with the empirical data.

## The bug

Inside the per-chromosome loop, the cluster means are computed as

```r
cluster_i <- which(tmp$marker %in% cluster_tmp$marker)      # row numbers WITHIN this chromosome (tmp)
tmp[cluster_i, "aq_fix"]  <- mean(aims_raw$faqu_fix[cluster_i], na.rm=T)   # ...used on the GENOME-WIDE table
tmp[cluster_i, "pol_fix"] <- mean(aims_raw$fpol_fix[cluster_i], na.rm=T)
```

`cluster_i` indexes rows of `tmp` (one chromosome), but it is applied to `aims_raw` (all
chromosomes). Chromosome 1 comes first in `aims_raw`, so it is correct. **Every other
chromosome gets chromosome 1's fixation levels at the same row numbers.** On chromosomes with
more markers than chromosome 1 (7, 17, 25), rows past chromosome 1's end read into whichever
chromosome follows it in the input table.

The SLiM model reads the AIMs files correctly. The LD-cluster structure and all positions are also correct
(using the `min_r2 = 0.2` clustering is intentional, see below). Only the per-cluster founder
probabilities (`pm10_aq`, `pm10_pol`) are wrong.

## Evidence

The input table used for the repo's AIMs files (51,612 markers = `group_info_new.rds`) is not
in our repositories, so it was rebuilt from the 30 reference parents of the DI25 panel. It uses the same
definition as `species_diagnostic_marker_fixation_levels.tsv`: `faqu_fix` = fraction of *aquilonia*
parents homozygous for the aquilonia allele, and `fpol_fix` likewise. This matches the legacy table on 98.6–98.9% of shared markers.
The rebuilt table is `species_diagnostic_marker_fixation_levels_DI25_rebuilt.tsv`.

- Running the **original** script on the rebuilt table reproduces the repo's `SLiM/AIMs/*.txt`
  **exactly** (same cluster ids, identical values) on 23 of 26 chromosomes. The 3 exceptions are
  7, 17 and 25, where the index runs past chromosome 1 (see above).
- In the repo's AIMs files, only 0–0.7% of founder probabilities on chromosomes 2–27 equal the
  intended cluster mean. Their correlation with the marker's own value is ≈ 0 (−0.07 to 0.06).
  Their correlation with chromosome 1's value at the same row number is 0.33–0.72.
- The simulations follow the buggy input. In `diem_boot1_output.bed`, the simulated parental
  |p_aq − p_pol| tracks the intended founder difference on chromosome 1 (Spearman 0.74) but not
  elsewhere (mostly near 0, range −0.48 to 0.60; the 0.60 is chromosome 17, one of the chromosomes where the index runs past chromosome 1).

## The fix (`AIMs_for_SLiM_fixed.R`)

Same inputs, same output format (`chr_id pos_0based cluster pm10_aq pm10_pol`), same cluster
numbering, so the SLiM model needs no change. Changes:

1. Fixation levels are taken from the chromosome's own rows, looked up by marker (`chr:pos`),
   never by a row index into the genome-wide table.
2. AIMs past the chromosome end are now actually trimmed. The original computed `aims_trimmed`
   but wrote the untrimmed table. This had no effect on the current data, because no AIM lies
   past a chromosome end.
3. `options(scipen=999)`, as in `convert_recmap.R`, so positions are never written as `6e+06`.
4. The script fails loudly if a cluster would get `NaN` (all members `NA`).

## The check (`check_AIMs.R`)

```bash
# Part 1: AIMs files vs the inputs they were built from (run BEFORE SLiM)
Rscript check_AIMs.R <AIMs_dir> species_diagnostic_marker_fixation_levels.tsv group_info_new.rds
# Part 2 (optional): one bootstrap DIEM output AFTER the rerun
Rscript check_AIMs.R <AIMs_dir> species_diagnostic_marker_fixation_levels.tsv group_info_new.rds diem_boot1_output.bed
```

Part 1 must print `PART 1 (AIMs files): PASS`. That requires:
- every input marker appears exactly once;
- the cluster partition equals `group_info_new.rds`;
- `pct_values_correct` is 100 on every chromosome;
- `cor_own_marker` is high (< 1 only because markers in a cluster share the cluster mean) and
  `cor_chr1_same_row` ≈ 0 except on chromosome 1.

Part 2 compares the simulation with the *intended* values from the input table, not with the AIMs
files, so a buggy AIMs file cannot pass it. Spearman should be clearly positive on every chromosome.

Results already obtained here:

| Run | Part 1 | Part 2 |
|---|---|---|
| Fixed script, rebuilt 51,612-marker table | **PASS**: 100% on all 26 chromosomes, `cor_chr1_same_row` −0.06 to 0.05 | — |
| Repo `SLiM/AIMs/*.txt` + current `diem_boot1` | **FAIL**: 0–0.7% correct off chromosome 1 | ≈ 0 except chromosome 1 (0.74) |

## Secondary point (not a bug; worth deciding before the rerun)

`faqu_fix` / `fpol_fix` are the fractions of parents **homozygous** for the
species allele, and they are used directly as the probability that a founder
haplosome carries it. The allele frequency is higher (homozygotes plus half the
heterozygotes). On the DI25 panel, median parental Δ is 0.80 from fixation levels
vs 0.90 from allele frequencies. Most of the lower differentiation of the simulated
parents (median |Δp| 0.77) comes from this, plus cluster averaging. If founders
are meant to carry the empirical parental allele frequencies, use
`(2·n_hom + n_het) / (2·n_valid)` per species. That needs the heterozygote counts,
which the fixation table doesn't contain but the parental genotypes do.

## What needs rerunning

Regenerate the AIMs files with the fixed script, run the check, then rerun SLiM and the
VCF→DIEM bootstrap pipeline. Everything downstream of `diem_outs_demo` should be redone on the new
runs: the neutral sorting null, the ld_w landscape comparison, and FST vs DI.

## Note on the clustering used (intentional, not a bug)

The simulations use the `min_r2 = 0.2` clustering (`group_info_new.rds`), not the rho05 one used
in the empirical analyses. This is deliberate. The founder set-up treats clusters as independent,
although physically adjacent clusters are in reality still correlated. Looser clusters (more
weakly correlated SNPs included) offset that, so founders start with realistic LD without a full
burn-in. The LD patterns are then regenerated by the simulation itself.
