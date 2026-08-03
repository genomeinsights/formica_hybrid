# Build BayPass input files for the Aland-excluded comparison analysis.
# See moduleC_write_baypass_inputs.R for the shared logic (including the
# poolsize-ordering fix) and R/moduleC_prepare_with_aland.R for the full-dataset
# counterpart used for comparison.

source("./moduleB_climate_GEA/R/moduleB_write_baypass_inputs.R")

load("./data/hybrids_only_maf005.Rdata")
pruned_markers <- readRDS("./data/pruned_markers.rds")

write_baypass_inputs(
  GTs = GTs_hybrids_005, map = map_hyb_005, sample_data = sample_data,
  pruned_markers = pruned_markers, out_folder = "./aland_excluded/",
  exclude_population = "Aland"
)
