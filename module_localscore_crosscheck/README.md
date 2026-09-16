# module_localscore_crosscheck

**Status: clean, confirmatory implementation.** This module independently regenerates and
validates the BayPass local-score outlier-region cross-check (Fariello et al. 2017; Bonhomme
et al. 2019) applied to the Stage-1-direct climate/mitotype scan, including an explicit test
of whether local-score windows are biased toward large, low-recombination Stage-2 LD
clusters.

## Relationship to `module_localscore_crosscheck_exploration`

The original exploratory work on this question lives in
[`module_localscore_crosscheck_exploration/`](../module_localscore_crosscheck_exploration/)
(archived by directory rename, full history preserved). That module is **read-only reference
material**:

- it is not modified, patched, rerun, or regenerated from this point forward;
- its code may be inspected to understand file formats, conventions, or locate source data;
- **no result, table, RDS object, or figure from the exploratory module is copied into this
  one.** Every number reported here is independently regenerated from the authoritative
  inputs in `module_manuscript_rho05/` and validated with its own checks, not inherited from
  the exploratory scripts' assumptions.

The exploratory module found (among other things) a real but mechanistically unresolved
resolution-dependent anomaly in the PC1/bio\_winter Stage-1-cluster BayPass comparison, and a
suggestive but not rigorously quantified link between local-score window detection and large,
low-recombination Stage-2 LD clusters. This module's job is to either confirm or revise those
observations with a clean, pre-specified, adequately-powered design (10 structured + 10
paired-unstructured full-SNP null replicates, rather than 1 informally-selected draw; explicit
window-level LD/redundancy annotation; genome-wide null-susceptibility analysis rather than
detected-windows-only).

## Structure

- `R/` -- analysis scripts (numbered/staged; each documents its inputs/outputs in its own header)
- `config/` -- run parameters, covariate/draw selection, checksummed input manifests
- `manifests/` -- provenance records (source file checksums, parameter logs, validation results)
- `data/` -- compact derived summaries (TSVs, small RDS); never raw multi-GB BayPass output
- `Figures/` -- final figures (repository-policy gitignored like every other module's `Figures/`)
- `doc/` -- the confirmatory report (LaTeX + compiled PDF)
- `full_snp_null10/raw/` -- raw BayPass output for the 10-replicate full-SNP null comparison
  (multi-GB; gitignored; regenerable from `R/` + `config/`)

## Ground rules

- Every reported result is regenerated from the authoritative inputs in
  `module_manuscript_rho05/baypass_stage1/` (observed BayPass scans, preserved null
  covariate/contrast files) and from the canonical Stage-2 rho05 clustering object -- not from
  any file under `module_localscore_crosscheck_exploration/`.
- Stage-1-resolution scripts use `best_marker` positions from the outset (not inherited/patched
  from the exploratory Stage-1 scripts).
- Explicit validation checks (marker order, population order, join integrity, dimension
  checks) are run and recorded before any result is interpreted -- see `manifests/`.
