# module_external_validation

External validation of the hybrid-zone BayPass outlier regions in an independent set of **21 hybrid
workers** (one per site; Finland 2, Sweden 16, Germany 3), using EMMAX on the same climate covariates.
Run every script from the repo root (`~/gitlab/formica_hybrid`). Data inputs live in `data/`:
`NooraNewHybridSamples.xlsx` (site metadata) and the 21-sample VCF
`HybridSamples_SNPQ30...FiSeDe.vcf.gz` (same assembly/coordinates as the hybrid data; INFO NS/AC/AF
are stale from a 384-sample call, MAF is recomputed). Large intermediates go to
`module_external_validation/data/` (git-ignored by the repo's `/data/` rule).

| Script | What it does | Key outputs |
|---|---|---|
| `R/00_audit_external_vcf.R` | VCF vs original SNP set; outlier units whose best SNP / LD-cluster members exist in the VCF | `results/audit_summary.txt`, `outlier_units_vs_vcf.tsv` |
| `R/01_project_climate_pca.R` | Tries to reproduce original bio1-19 from WorldClim and project the new sites (diagnostic: **original climate is not plain WorldClim 2.1**, ~1.65 C warmer in bio1) | `results/climate_pca_projection_checks.txt` |
| `R/01b_climate_pca_worldclim_consistent.R` | PCA rebuilt from WorldClim 30 s alone on the 20 original sites (fidelity r = 0.997 / 0.977 to stored PC1/PC2); new samples projected | `data/new_sample_climate_pca_wc.tsv` |
| `R/02_ld_decay_external.R [slide]` | LD decay for MAF>0.1 SNPs (230,392 SNPs after dropping duplicated positions); `slide=100` is the working fit (~128 kb window) | `data/ld_decay_ext_maf10_slide100.rds` |
| `R/03_stage1_clustering_external.R` | Stage-1 clustering (rho = 0.5): 155,665 clusters -> LD-pruned set for the GRM | `data/ext_pruned_markers_maf10.txt` |
| `R/04_compare_decay_rates.R` | per-chromosome decay: external vs original (raw a r = 0.64; a_pred only reflects chromosome size) | `results/decay_rate_comparison.*` |
| `R/05_compare_ldw095.R` | ld_w_095 external vs original (per SNP Spearman 0.37; 100 kb bins 0.80) | `results/ldw095_comparison.*` |
| `R/06_subsample_refit_original.R` | 21-sample subsamples of the original: background LD reproduced (0.214); decay rate `a` depends strongly on marker density/window | `results/subsample_refit_original.*` |
| `R/07_build_tested_snps.R` | one SNP per outlier region (best marker of top unit; else other unit's best; else LD proxy |r|>=0.8) | `results/tested_snps.tsv` |
| `R/09_orientation_ancestry.R` | ancestry-axis inference of whether VCF ALT == original coded allele (orientation-free PCA + parental AFs) | `data/orientation_ancestry.rds`, `results/orientation_ancestry_summary.txt` |
| `R/10_orientation_polarity_rule.R` | rule **ALT == coded allele iff DIEM Polarity == 0**; validation (ancestry inference, LD-phase, AF); oriented effects | `results/orientation_validation.txt`, `emmax_tested_snps_oriented.tsv` |
| `R/11_concordance_test.R` | sign concordance of external effects with original BayPass direction (Beta_is convention checked empirically: refers to the NON-coded allele in 101/101 crossing units) | `results/concordance_summary.txt`, `concordance_tested_snps.tsv`, `concordance.png` |
| `R/12_outlier_snps_to_external_units.R` | ALL outlier SNPs (members of raw-crossing original units) -> external Stage-1 units -> best SNP per unit, direction predicted from the original data, region-level permutation | `results/extunits_summary.txt`, `extunits_tested.tsv` |
| `R/utils_emmax.R` | shared EMMAX helpers (beta/SE, rank-normal, imputation) | - |
| `R/08_grm_emmax.R` | GRM from pruned markers; genome-wide EMMAX (lambda_GC ~1) + tested-SNP table, PC1/PC2 raw and rank-normal | `results/emmax_tested_snps.tsv`, `grm_emmax_summary.txt`, `emmax_qq.png` |

## Status / open issues
* **Allele orientation resolved (09/10).** Original coded allele = F. aquilonia-associated allele (100% of
  DI > -25 markers); VCF ALT == coded allele iff DIEM `Polarity == 0`. Validated by the ancestry axis of the 21
  hybrids (99.8% agreement at |cor| > 0.6, Polarity not used), by LD phase of 21,394 neighbouring strong-LD pairs
  (100% agreement with the rule, 94.5% without orientation; 100% also for 1,180 different-polarity pairs), and
  AF consistency. `emmax_tested_snps_oriented.tsv` has `beta_coded` (per coded allele). Original `Beta_is` refers to the NON-coded (dosage-0)
  allele (101/101 BF>=15 units); the `flipped` flag only affects filling of missing calls.
* **Concordance result (11).** PC1 8/16 SNPs concordant (binomial p = 0.60), PC2 6/9 (p = 0.25); sign-flip
  permutation of sum(e*z) p = 0.66 / 0.21; rank-normal PCs the same. No evidence of replication, but power is very
  low (>=11/16 or >=8/9 concordant would be needed for p < 0.05).
* **All-outlier-SNP extension (12).** 963 PC1 / 472 PC2 outlier SNPs -> 106 / 70 present in the external set -> 39 / 20
  external Stage-1 units (38 / 15 with a predicted direction, in 18 / 9 original regions -- vs 16 / 9 regions with the
  best-SNP design). Concordance: PC1 18/38 units, 9/18 regions (perm p = 0.44); PC2 10/15 units, 6/9 regions (p = 0.25).
  Same conclusion; the limit is how many outlier SNPs exist (MAF>0.1, in the VCF) in the external data, not the unit definition.
* **Climate source mismatch.** The new samples use WorldClim 2.1 (30 s) PCs; 11/21 lie outside the original
  PC range (PC2 up to ~19 for SW-Sweden samples, original max 2.7) -> rank-based PCs reported alongside raw.
* **Testable regions** (n = 21, MAF>0.1, <=20% missing): PC1 16/43, PC2 9/22, mitotype 4/13 (no mitotype
  in the nuclear VCF), coastal-inland 5/22, heat tolerance 10/27 (no phenotype). Only PC1/PC2 analysed.
* Positions occurring in >1 VCF record (split multiallelics) were dropped (52,110 of the MAF>0.1 records).
