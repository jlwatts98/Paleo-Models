##### Calculating percent change #####
# source
source("Header.R")

# read in data
dist_change = readr::read_rds("objects/dist_change.rds")

# change function over all species
perc_change = lapply(dist_change, perc_change)

# make into dataframe
library(data.table)
perc_change_df = rbindlist(perc_change) |>
  as.data.frame() |>
  rename(species = V1,
         pp_stable = V2,
         pp_loss = V3,
         pp_gain = V4,
         pf_stable = V5,
         pf_loss = V6,
         pf_gain = V7)

# write rds
readr::write_rds(perc_change_df, "objects/perc_change.rds")
