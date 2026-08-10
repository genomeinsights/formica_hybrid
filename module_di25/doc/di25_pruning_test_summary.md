# DI25 "early pruning" test — colony-level result (revised)

**Script:** `module_di25/R/di25_pruning_test.R` (additive; reads `di25_*` outputs)
**Run:** full N_BOOT = N_PERM = 1000, colony (nest) as the sampling unit
**Date:** 2026-08-08

## One-sentence verdict

Heterozygosity is directionally lost toward the fixed allele (135:0), and this
scales across the 19 modelled populations — with Sielva the endpoint — but the loss
is **concentrated in low-recombination regions across 18 of 19 populations**, which
is the signature of **linked purging**, the *opposite* of active recombination-
exposure pruning. So: the phenomenon is real; the "active pruning" mechanism is not.

- **Directional loss toward the fixed allele:** REAL (135:0, unit level).
- **19-population continuum:** real but **marginal** (lm p = 0.021, n = 19).
- **Mechanism:** het loss steeper at *low* recombination in 18/19 pops
  (sign-test p = 7.6e-5) = **linked purging**, not recomb-exposure pruning.

## What changed from the previous version (5 corrections)

1. **"Three large blocks" is really one.** Chr26 = 5.47 Mb (65 markers/Mb) is the
   only large block; Chr5 = 0.33 Mb (796 markers/Mb) and Chr25 = 0.23 Mb (1725
   markers/Mb) are sub-megabase but marker-dense → wide circos arcs, narrow regions.
   Say "~6 Mb total" out loud, and "one large Chr26 block plus two small dense
   regions," not "three large blocks." (A separate dense consensus-polyctena region
   sits on Chr5 at ~2.5–3.2 Mb, across a real 1.3 Mb gap from the anchored block.)
2. **Regression CI was measuring the wrong thing.** The 19-slope regression is over
   POPULATIONS, so it is now resampled at the population level (bootstrap +
   permutation), with the lm p as the honest headline. The chromosome block
   bootstrap is kept only for Sielva's within-population M1, where it belongs.
   Result: lm p = 0.021, pop-bootstrap CI [0.64, 9.37], pop-perm p = 0.018 —
   marginal, not "perm p ≈ 0."
3. **The mechanism is a positive result, not a null.** 18 of 19 populations show a
   POSITIVE interaction (het loss steeper at low recombination) — sign-test
   p = 7.6e-5. That is linked purging, consistent with Nouhaud's low-recombination
   purging and Result 3 (sorting stronger at low recombination). Sielva's lone
   negative rests on 1 colony and is best read as noise.
4. **F1 verdict retracted.** Tract length cannot separate "F1 + scattered error"
   from "early backcross." A Wald–Wolfowitz runs test on Sielva's homozygous units
   shows significant but modest spatial clustering (Z = −4.69, p = 3e-6): an
   F1-like background carrying some real ancestry tracts — not a clean F1, not
   clearly a backcross. "Leans early-backcross" is dropped.
5. **19 populations, and Sielva is n = 1.** Svanvik1 (n = 3) is excluded, so 19 are
   modelled, not 20. Nine of the 20 are single-nest; Sielva is one, with binomial
   weight 1. Its position on the regression line is genuine but rests on a single
   genome, not on within-population replication.

## Falsification criteria (Section 9)

```
  (1) effect gone at unit level ...................... no
  (2) effect gone with blocks excluded ............... no
  (3) Sielva off-line outlier, not endpoint .......... no
  (4) recomb-exposure PRUNING mechanism ............. UNSUPPORTED (interaction is +ve: linked purging)
  (5) missingness reproduces the gradient alone ...... no
  (6) 5 Sielva workers are one genotype .............. no
```

## Console readouts

### Setup

```
Individuals: 164 hybrids in 20 populations + 30 parents
Units: 11052 (4010 eMLG consensus + 7042 rep-SNP); 10973 oriented (flipped)
Modelled populations: 19 of 20 (MIN_N_IND = 4). Excluded (too few individuals): Svanvik1 n=3.
Sampling unit = COLONY (nest): 9 of 20 pops are single-nest (Sielva = 1 nest, 5 workers).
  Nests per pop: Karsikas:1 Heinamaki:2 Nyrhispera75:1 Nyrhispera74:1 Kummunmaki:3
  Jarvenpaa:1 Parikkala:1 Katiskoski:2 Hiivola:2 Vuosaari:1 Pikkala:2 Aland:3
  Grundsund:2 LangholmenW:7 LangholmenR:3 Tvarminne:1 Bunkkeri:2 Svanvik1:1 Svanvik2:2 Sielva:1
```

### [2a] Sielva 5-worker IBS matrix (is it one genotype?)

```
          Siel338_a Siel338_b Siel338_c Siel338_d Siel338_k
Siel338_a     1.000     0.951     0.950     0.952     0.956
Siel338_b     0.951     1.000     0.955     0.962     0.957
Siel338_c     0.950     0.955     1.000     0.956     0.953
Siel338_d     0.952     0.962     0.956     1.000     0.961
Siel338_k     0.956     0.957     0.953     0.961     1.000
     Sielva mean within-pop IBS = 0.955 ; rank among 20 pops (1=most clonal): 1
```

### [2b] Missingness confound

```
     Spearman(unit missingness, mean sortedness) = 0.239
     Sielva missingness-ONLY het slope: -1.525 (M1 full missingness coef: -0.321)
```

### [2c] F1? — tract lengths NOT diagnostic; spatial runs test is

```
   n_tracts med_units med_cM max_cM frac_ge5cM
1:     9398         3   0.64  18.04       0.08
     Wald-Wolfowitz runs on Sielva-consensus homozygous units: 22% hom; runs Z=-4.69, p=2.8e-06
     -> homozygous units are SIGNIFICANTLY clustered (real tract structure), but modestly --
        F1-like background carrying some tracts; NOT a clean F1, NOT clearly a backcross.
```

### [3/5] 19-slope regression (observations = POPULATIONS)

```
               Estimate Std. Error t value Pr(>|t|)
(Intercept)      -3.294      0.449  -7.339    0.000
own_sortedness    3.738      1.476   2.532    0.021
     lm p = 0.021 (n=19, the HONEST test); pop-bootstrap 95% CI [0.642, 9.366]; pop-perm p = 0.018
     Sielva studentized residual: -0.00 (on the line); Sielva = 1 colony (weight 1).
```

### [4] M1 recombination interaction across the 19 populations (the mechanism test)

```
     18 POSITIVE / 1 negative; sign-test p = 7.6e-05; mean interaction (excl Sielva) = +0.876
     POSITIVE = het loss steeper at LOW recombination (linked purging); recomb-EXPOSURE pruning predicts NEGATIVE.
     Sielva is the lone negative (-1.015; block-boot CI [-2.09, -0.04]; perm p=0.066) -- 1 colony, noise.
     Sielva het-loss slope -2.685 (block-boot CI [-2.97, -2.44]; perm p=0.000) -- direction, not mechanism.
```

### [4] M2 direction of Sielva homozygotes (sortedness > 0.75, 1 nest = 1 obs/unit)

```
      pop n_hom n_toward n_against frac_toward      p_binom
1: Sielva   135      135         0           1 4.591775e-41
```

> Caveat: the 0.5 null holds only if the 20 populations did not draw from one shared,
> already-skewed founding pool. The absent pruning interaction makes shared skew *more*
> plausible, so this is worth testing (Module E shared-source), not assuming.

### [7] Blocks — one large + two small dense regions (say Mb, not marker count)

```
     chr n_units n_markers span_Mb markers_per_Mb
1:    26       4       356    5.47             65
2:     5       4       262    0.33            796
3:    25       2       395    0.23           1725
     ~6 Mb total. Sielva het: 80% IN the blocks vs 51% in NON-block units at matched sortedness (>0.75).
     DESCRIPTIVE only: 10 block units split across 2 bins (~5/bin) -- describe, do not test.
     M1 with blocks EXCLUDED -- Sielva sortedness slope -2.698 (was -2.685); interaction -0.832 (was -1.015)
     19-slope regression (blocks excluded): slope +3.768 (lm p=0.021)
```

## Outputs

- Figures: `module_di25/Figures/di25_pruning_{slopes,interaction,sielva}.png`
  (interaction is now a 19-population forest plot; slopes carries the population-level CI)
- Caches: `module_di25/data/di25_pruning_{models,slopes,boot,perm}.rds`,
  `di25_three_blocks.rds` (now includes `markers_per_Mb`)
