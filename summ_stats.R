##### Summary Statistics #####
# source
source("Header.R")

# read in data
dist_change = readr::read_rds("objects/dist_change.rds")
predictors = readr::read_rds("objects/predictors.rds")

# define dem
dem = predictors[[21]]

# apply custom function summ_stats
summ_stats = summ_stats(dist_change[[1]], dem)


