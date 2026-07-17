# Build BayPass input files for the Sielva-excluded comparison analyses --
# Sielva is a very young population with unusually high heterozygosity
# genome-wide, a similar potential confound to Aland (see
# R/prepare_aland_excluded.R). Writes TWO configurations, both repeating
# the existing with_aland/ vs. aland_excluded/ comparison but with Sielva
# also excluded from each side:
#   - ./sielva_excluded/       -- Sielva out, Aland IN  (parallel to with_aland/)
#   - ./aland_sielva_excluded/ -- Sielva out, Aland OUT (parallel to aland_excluded/)
#
# See write_baypass_inputs.R for the shared logic (including the
# poolsize-ordering fix). Uses the SAME pruned marker set
# (./data/pruned_markers.rds) as with_aland/ and aland_excluded/ so all
# four population-exclusion configurations stay directly comparable.

source("./R/write_baypass_inputs.R")

load("./data/hybrids_only_maf005.Rdata")
pruned_markers <- readRDS("./data/pruned_markers.rds")

write_baypass_inputs(
  GTs = GTs_hybrids_005, map = map_hyb_005, sample_data = sample_data,
  pruned_markers = pruned_markers, out_folder = "./sielva_excluded/",
  exclude_population = "Sielva"
)

write_baypass_inputs(
  GTs = GTs_hybrids_005, map = map_hyb_005, sample_data = sample_data,
  pruned_markers = pruned_markers, out_folder = "./aland_sielva_excluded/",
  exclude_population = c("Aland", "Sielva")
)
