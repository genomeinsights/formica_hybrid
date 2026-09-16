# Stage-1 unit position table -- manifest
Source: `module_manuscript_rho05/data/moduleB_stage1_units_bestsnp.rds`
N Stage-1 units: 18361
Units where representative == best_marker: 10313 (56.2%)

Validation performed:
- every best_marker exists in map_hyb_005 (the authoritative marker map)
- map-derived Chr/Pos agrees exactly with best_marker's own Chr:Pos string
- Stage-1 unit group_id set matches the BayPass run's S1units_group_order.txt exactly
- no duplicate group_id, no duplicate (Chr,Pos), no missing Pos

Output: module_localscore_crosscheck/data/stage1_unit_positions.tsv
(group_id, Chr, Pos, best_marker, representative, best_r, rep_is_best, n_loci)
