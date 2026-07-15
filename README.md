# formica_hybrid

Analyses of climate adaptation in hybrid populations of F. polyctena and F. aquilonia

## Pipeline

## 1. Parse and save DIEM data

./R/Parse_DIEM_data.R

## 2. Estimate LD-decay

./R/LD_decay.R

## 3. BayPass analyses

./R/baypass.R

Uses LDscnR::ld_complexity_reduction() to prune markers for GRM construction.

## 4. LDscnR of PC1 and PC2

./R/LDscnR_PC12.R

## 5. Joint LD-weight / XtX outlier filtering

./R/LD_XtX_outlier_filtering.R

Combines BayPass (step 3) and LD-decay (step 2) results to find outlier SNPs that are stable across LD-weight and XtX significance thresholds.

## 6. eMLG complexity reduction

./R/eMLG_complexity_reduction.R

Collapses SNPs into effective marker-locus groups (eMLGs) using LD clustering.

## 7. Haplotype blocks and LD-circular plots

./R/hap_blocks.R

## Shared functions

./R/Functions.R holds plotting/analysis helpers (Manhattan+DIEM overlay plots, haplotype block detection) used across the pipeline scripts above.
