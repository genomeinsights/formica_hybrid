# Build BayPass input files for the full dataset (all populations,
# including Aland), regenerated with the poolsize-ordering fix so it's a
# valid comparison against ./aland_excluded/ -- the existing out_baypass_2
# results were generated before that fix and are not directly comparable.
# See write_baypass_inputs.R for the shared logic.

source("./R/write_baypass_inputs.R")

load("./data/hybrids_only_maf005.Rdata")
pruned_markers <- readRDS("./data/pruned_markers.rds")

write_baypass_inputs(
  GTs = GTs_hybrids_005, map = map_hyb_005, sample_data = sample_data,
  pruned_markers = pruned_markers, out_folder = "./with_aland/",
  exclude_population = NULL
)
