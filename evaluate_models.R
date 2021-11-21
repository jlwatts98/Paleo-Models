##### Evaluate Models #####
source("Header.R")
# read rds object
models = readr::read_rds("objects/models.rds")
predictors = readr::read_rds("objects/predictors.rds")

# evaluate
evaluations = lapply(models, evaluate_model, predictors)

# make into dataframe
library(data.table)
evaluations_df = rbindlist(evaluations) |>
  as.data.frame() |>
  rename(species = V1,
         AUC = V2,
         TSS = V3,
         AICc = V4)
evaluations_df$AUC = round(evaluations_df$AUC, digits = 2)
evaluations_df$TSS = round(evaluations_df$TSS, digits = 2)
evaluations_df$AICc = round(evaluations_df$AICc, digits = 0)

# write rds object
readr::write_rds(evaluations_df, "objects/evaluations.rds")
