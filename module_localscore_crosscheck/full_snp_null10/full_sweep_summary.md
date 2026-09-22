# Full 10-replicate x 2-mode x 2-null-type sweep -- completion summary

Status: 40/40 BayPass runs completed successfully on mini2, raw output on
`/Volumes/T9/module_localscore_crosscheck_full_snp_null10/raw/` (not synced
to this repo -- see `.gitignore`). Compact per-SNP summary files (
`_summary_betai_reg.out` for continuous, `_summary_contrast.out` for
mitoC2) copied back to `full_snp_null10/raw/` here (gitignored, kept on
disk locally); the aggregated table below (`sweep_comparison_stats.rds`) is
the only sweep output committed to git.

## Timing

- Launched: 2026-09-17 09:33:53
- Completed: 2026-09-19 17:02:06
- Total wall-clock: ~55.5h (estimate before launch was ~57h)
- Per-run range: ~1h20m-1h30m, consistent across all 40 runs, no crashes or
  restarts.

## Validation

- All 40 output files present (10 replicates x 2 modes x 2 null types).
- Every file has exactly 1,114,423 rows (= the full SNP count in
  `u_DIEM.geno`), confirmed for all 40 files individually, not just
  spot-checked.

## Continuous mode (BF(dB), threshold BF>=15 and BF>=20 as in precedent usage)

| null type | mean N(BF>=15) across 10 draws | mean N(BF>=20) |
|---|---|---|
| structured | 2004.8 | 395.4 |
| unstructured | 2140.5 | 467.5 |

Per-draw range: structured N(BF>=15) 1593-2723; unstructured N(BF>=15)
1880-2487. The pilot's single draws (structured 2047, unstructured 2384)
sit within this replicate-to-replicate range -- consistent with the full
10-replicate result, not an outlier.

## mitoC2 mode (log10(1/pval), threshold >=3 as in precedent usage)

| null type | mean N(log10(1/pval)>=3) across 10 draws |
|---|---|
| structured | 986.4 |
| unstructured | 877.5 |

Per-draw range: structured 736-1360; unstructured 757-991.

## What this does NOT yet establish

This is validation of the raw sweep output only (dimensions, completeness,
descriptive exceedance counts at two conventional thresholds). It is not
yet the local-score analysis itself, and not yet a formal structured-vs-
unstructured comparison (e.g. paired test across the 10 replicates,
distributional comparison, or false-positive-rate calibration). That
belongs in the next stage: running `compute.local.scores()` independently
on each of these 40 outputs and building the confirmatory document around
the resulting window-level comparison, per the module's original
instructions.
