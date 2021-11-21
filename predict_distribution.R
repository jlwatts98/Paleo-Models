##### Predict Models #####
source("Header.r")
# read in model data and predictors
models = readr::read_rds("objects/models.rds")
predictors = readr::read_rds("objects/predictors.rds")
fut_predictors = readr::read_rds("objects/fut_predictors.rds")
past_predictors = readr::read_rds("objects/past_predictors.rds")

# retrain models
retrained = lapply(models, retrain_model)

# write rds object
readr::write_rds(retrained, "objects/retrained.rds")

# read rds object
retrained = readr::read_rds("objects/retrained.rds")

# predict
predictions = lapply(retrained, predict_model, predictors)
fut_predictions = lapply(retrained, predict_model, fut_predictors)
past_predictions = lapply(retrained, predict_model, past_predictors)

# write rds object
readr::write_rds(predictions, "objects/predictions.rds")
readr::write_rds(fut_predictions, "objects/fut_predictions.rds")
readr::write_rds(past_predictions, "objects/past_predictions.rds")
