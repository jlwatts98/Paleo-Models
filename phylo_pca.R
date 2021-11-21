##### Phylogenetic PCA #####
# https://cran.r-project.org/web/packages/geomorph/vignettes/geomorph.PCA.html
# source header
source("Header.R")

# load in cleaned presence records
cleaned_pres = readr::read_rds("objects/cleaned_pres.rds")

# load in predictors
predictors = readr::read_rds("objects/predictors.rds")

# extract function
env_vals = lapply(cleaned_pres, get_env_vals, predictors)

env_vals_df = rbindlist(env_vals) |>
  as.data.frame()

# write rds
write_rds(env_vals_df, "objects/env_vals.rds")

##### Phylogenetic tree #####
library(geomorph)


