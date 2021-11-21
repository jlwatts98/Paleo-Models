##### Build Models #####
source("Header.R")

# load in SWD objects and raster data
SWDs = readr::read_rds("objects/SWDs.rds")
predictors = readr::read_rds("objects/predictors.rds")

# run the models
model1 = model_tune_maxent(SWDs[[1]], 15, predictors)
model2 = model_tune_maxent(SWDs[[2]], 15, predictors)
model3 = model_tune_maxent(SWDs[[3]], 15, predictors)
model4 = model_tune_maxent(SWDs[[4]], 15, predictors)
model5 = model_tune_maxent(SWDs[[5]], 15, predictors)
model6 = model_tune_maxent(SWDs[[6]], 15, predictors)
model7 = model_tune_maxent(SWDs[[7]], 15, predictors)
model8 = model_tune_maxent(SWDs[[8]], 15, predictors)
model9 = model_tune_maxent(SWDs[[9]], 15, predictors)

models = list(model1, model2, model3, model4,
              model5, model6, model7, model8,
              model9)

# write an rds object
readr::write_rds(models, "objects/models.rds")
