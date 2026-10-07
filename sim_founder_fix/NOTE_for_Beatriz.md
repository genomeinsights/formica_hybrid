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

## Second issue: founder LD does not match the empirical parents (see `note_parent_ld/NOTE_parent_ld.pdf`)

Independently of the indexing bug, the LD among the simulated parents (founders) is
*binary*: SNPs in the same `min_r2 = 0.2` cluster are in perfect LD (r² = 1.000), SNPs in
different clusters are independent (r² ≈ 0.07, the floor for n = 15). Empirical parental LD
is graded (0.27–0.31 and 0.16–0.20 for the same close pairs), so simulated short-range LD is
2–3× too high and extends to 0.2–1 cM. Simulated *F. aquilonia* parents are haploid males;
simulated *F. polyctena* parents are twice as heterozygous as the empirical ones.

Proposed alternative: founders as genotype mosaics of the 15 empirical parents per species
(switch donor at ~1 per cM). They reproduce the empirical LD decay, graded close-pair LD,
allele frequencies and heterozygosity, and pass an acceptance test that the current
founders fail. The first version split heterozygous sites at random between haplosomes; the
recommended version copies **phased** parental haplotypes (Beagle 5), which also keeps the
within-species LD among polymorphisms shared by the two species (see the handoff section
below). Files:

| File | Purpose |
|---|---|
| `make_mosaic_founders.R` | writes `founders_ch<id>.vcf` (haploid `aq_hapNN`, phased diploid `pol_femNN`; the format of Petri's earlier real-founder SLiM script, read with `readHaplosomesFromVCF`) |
| `check_parent_ld.R` | acceptance test: `bed` mode for DIEM outputs, `vcf` mode for founder pools |
| `phase_parents.R` | Beagle phasing of the 30 parents (already run; output committed as `data/parents_phased.rds`) |
| `data/parents_phased.rds` | phased parental haplotypes + marker map (Chr, Pos, cM, DI25/near-neutral class); the only input the founder generator needs |
| `data/neutral_markers.txt` | IDs (`ChrN:pos`) of the 14,093 near-neutral SNPs |
| `slim/` | SLiM model with VCF founders and the job scripts (handoff section below) |
| `parent_ld_lib.R` | shared LD functions |
| `make_note_figures.R` | figures for the note |
| `note_parent_ld/` | the note (LaTeX source, PDF, figures) |

Recommended addition: include ~14,100 near-neutral SNPs (DI <= -90, pooled parental MAF >= 0.15) as a
calibration anchor — the null must reproduce the empirical background F_ST (~0.05) at these loci before a
shortfall at ancestry-informative loci counts as evidence. `make_mosaic_founders.R ... DI25+neutral` adds them
(see section 5 of the note).

Generated founder pools and check outputs go to `out/` (git-ignored).

## Handoff: rerunning the neutral simulations with phased mosaic founders

Everything needed is in this folder (`sim_founder_fix/`); the founder generator in phased mode
reads only `data/parents_phased.rds` (no other project data). Tested end to end in a clean
directory containing only `sim_founder_fix/` and the SLiM inputs from your repository.

### What changes compared with your current set-up

- **No AIM files.** Founders are not drawn per cluster from `pm10_aq`/`pm10_pol`; each founder
  copies stretches of real parental haplotypes. The indexing bug, the binary founder LD and the
  `min_r2 = 0.2` clustering all drop out.
- **Markers.** Every founder carries 65,705 SNPs at their empirical positions: the 51,612
  ancestry-informative SNPs (DI > -25) plus 14,093 near-neutral SNPs (DI <= -90, pooled parental
  MAF >= 0.15; listed in `data/neutral_markers.txt`). Intermediate DI (-90 < DI <= -25) is not
  included. The near-neutral SNPs are a calibration anchor: neutral simulations should reproduce
  their empirical F_ST before a shortfall at ancestry-informative SNPs counts as evidence.
- **SLiM model** `slim/SpecIAnt_rufa_mosaic_founders.slim` = your
  `SpecIAnt_rufa_genome_LDclusters_demo_server.slim` with only these changes (listed in its
  header): founders read with `readHaplosomesFromVCF` (initialN/2 haploid *F. aquilonia* males,
  initialN/2 diploid *F. polyctena* queens); paths and IDs as command-line constants instead of
  `sim_table.txt`; optional `CHROMS` (chromosome subset) and `RECSCALE` (recombination-map
  multiplier, default 1). Neutral as before (no QTN targets; SIGMA_W = 1000, DMI_W = 1).

### Running it

Directory layout (`BASE`): `sim_founder_fix/` plus your `slim_inputs/recombination_maps/`
(`ch_<id>.recmap`) and `slim_inputs/climate/climate_rep1.txt` (paths can be overridden).
Requirements: SLiM 5, R with `data.table` and `ggplot2`.

```bash
# one founding group x replicate (lan = lanR + lanW, bungrund = bun + grund share founders)
K=12500 INITN=100 bash sim_founder_fix/slim/run_group_job.sh aland 1 $BASE
# all 18 groups for replicates 1..N, e.g. as a SLURM array over "group run" lines
for r in $(seq 1 200); do for g in aland katis lan svan1 svan2 tvar bungrund pik nyr1 nyr2 \
  heina pari hiiv vuos kumm karsi jarven sielva; do echo "$g $r"; done; done > jobs.txt
```

Per job: the group's founder pool is generated (seed = 100000 + 1000 x group index + run; also
used as FOUNDER_SEED), SLiM runs once per member population and stops at the sampled cycle
(125; Sielva 11), and only `results/females_ckl<cycle>_<pop>_<RR>.vcf.gz` (100 sampled new
queens) and `results/<pop>_<RR>.anc` are kept (`KEEP=both` also keeps the males). Everything
else, including the founder pool, is deleted. Re-running a job skips finished populations.

Output VCFs: `#CHROM` = `ch<id>`, positions = empirical positions (`ChrN:pos` = marker ID),
founder SNPs have `MT=10`; allele 1 = the coded allele of the empirical genotype matrices. A
panel SNP missing from a VCF has lost allele 1 in that sample. For the DIEM pipeline: markers are
now individual SNPs, not AIM-cluster positions; drop the near-neutral SNPs
(`data/neutral_markers.txt`) for DI25 analyses, or keep them as the anchor.

### Settings, cost and checks

| | your current default | recommended |
|---|---|---|
| K (carrying capacity) | 6,250 | 12,500 |
| initialN (founders) | 100 | 100 |
| sampled cycle | 125 (Sielva 11) | 125 (Sielva 11) |

Recommendation from a calibration grid (chromosomes 1-6; K 6,250/12,500 x 100/1,000 founders x
cycles 60-1,000; 6 independent populations per setting): K 12,500 with 100 founders at cycle 125
matches the empirical near-neutral F_ST best (0.033 vs 0.031); your default overshoots it
(0.047). Under every setting the neutral model leaves ancestry-informative F_ST and within-population
fixation far below the empirical values (max 0.086 vs 0.26; max 0.9% vs 17.6% of unit x population
combinations monomorphic), which is the main result of `module_allele_specific_sorting` (script 07).

Cost (full genome, measured on a 14-core Mac mini): K 6,250 takes 13-22 min and about 19 GB peak
per population; memory scales roughly with K, so expect about 35-40 GB at K 12,500 (not yet
measured full genome). One replicate = 20 populations (Sielva runs only 11 cycles), i.e. about
5 core-hours at K 6,250; 200 replicates give permutation p-values down to 0.005, 1,000 down to 0.001.

Checks: each pool's `provenance.txt` should report 65,705 SNPs on 26 chromosomes and the SLiM
log `MOSAIC: loaded 65705 founder marker sites`. After the first few replicates, send us the
`results/` folder: `module_allele_specific_sorting/R/exploratory/explore_mosaic_sims.R` computes
near-neutral and DI25 F_ST, % sorted and the LD profiles against the empirical data.

### Known limitation

The neutral model overestimates short-range within-population LD for *all* loci, including the
near-neutral ones (about 3x at 0.05-0.2 cM even with the recombination map scaled tenfold); the
founders themselves match the parents, so the excess arises during the hybrid phase. Analyses
that compare LD itself (e.g. the ld_w landscape) against these simulations need this caveat; the
sorting null and F_ST versus DI are not affected.
