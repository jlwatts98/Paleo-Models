##### Functions #####
##### Remember to thinData via SDMtune #####

# function that builds an SDM from occurence points



triad_SWD_test2 = triad_SWD(par1 = "Asplenium montanum",
                           par2 = "Asplenium platyneuron",
                           hyb = "Asplenium bradleyi",
                           genus = "Asplenium",
                           specific_ep1 = "montanum*",
                           specific_ep2 = "platyneuron*",
                           specific_ep_hyb = "bradleyi",
                           year = 1900,
                           samp_dist = 200000,
                           samp_num = 3000,
                           predictors = clim_crop_test,
                           thin = TRUE)

triad_SDM_test2 = triad_SDM(triad_SWD_test2, 1)
triad_SDM_test2


evaluations_test2 = evaluate_triad_models(triad_SDM_test2, clim_crop)
evaluations_test2

triad_retrain1 = retrain_triad(triad_SDM_out = triad_SDM_test2)
triad_retrain1

triad_predict1 = predict_triad(triad_retrain1, clim_crop_test)
triad_predict1

hyb_test = bin_hybrid_zone(predict_triad_out = triad_predict)
raster::plot(hyb_test)

# asplenium montanum

SWD_test_thin = create_SWD("Asplenium montanum",
                     "Asplenium",
                     "montanum*",
                     1900,
                     200000,
                     3000,
                     clim_crop,
                     thin = TRUE)

model_test9 = model_tune(SWD_test_thin, 1)
model_test9

evaluation6 = evaluate_model(model_test9, clim_crop)
evaluation6

retrained = retrain_model(model_test9)

preds = predict_model(retrained, clim_crop)

raster::plot(preds[[2]])
raster::plot(preds[[3]])

# asplenium platyneuron

SWD_test1 = create_SWD("Asplenium platyneuron",
                      "Asplenium",
                      "platyneuron*",
                      1900,
                      200000,
                      3000,
                      clim_crop,
                      thin = TRUE)

model_test7 = model_tune(SWD_test1, 3)
model_test7

evaluation4 = evaluate_model(model_test7, clim_crop)
evaluation4

retrained1 = retrain_model(model_test7)

preds1 = predict_model(retrained1, clim_crop)

raster::plot(preds1[[2]])
raster::plot(preds1[[3]])

# load from gbif.

load = function(
  genus,
  specific_ep
){
  
  #checks
  checkmate::assert_character(genus)
  checkmate::assert_character(specific_ep)
  
  # load the data from gbif
  presence = gbif(genus = genus, species = specific_ep)
  
  #return
  return(presence)
  
}

# helps in omitting NAs

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

# cleans presence points

# clean
clean_pres = function(
  presence,
  year,
  predictors,
  thin
){
  
  # checks
  checkmate::assert_data_frame(presence)
  checkmate::assert_numeric(year)
  checkmate::assert_flag(thin)
  
  # clean presence data
  presence_complete = completeFun(presence, c("lon", "lat"))
  
  # remove duplicates
  presence_dup = presence_complete[!duplicated(presence_complete[c("lon","lat")]),]
  
  # remove records before year
  presence_year = presence_dup %>% filter(year > 1900)
  
  # coordinate cleaner package
  flags = clean_coordinates(presence_year, lon = "lon", lat = "lat",
                            species = "species", 
                            tests = c("capitals", "centroids", "equal", 
                                      "gbif", "institutions", "outliers", 
                                      "seas", "zeros"))
  
  #Exclude problematic records
  presence_final = presence_year[flags$.summary,]
  
  # get rid of metadata
  presence_xy = presence_final[,c("lon", "lat")]
  
  # if else statement for thin
  
  if (thin == TRUE) {
    # thin data
    presence_thin = SDMtune::thinData(presence_xy, predictors, "lon", "lat")
    
    # return
    return(presence_thin)
    
  } else {
    return(presence_xy)
    
  }
}

# gets background points for a certain species

get_backgr = function(
  presence,
  samp_dist,
  samp_num,
  predictors
){
  
  #checks
  checkmate::assert_data_frame(presence)
  checkmate::assert_numeric(samp_dist)
  checkmate::assert_numeric(samp_num)
  
  # selecting background points
  # make into spatial
  coordinates(presence) = ~ lon + lat
  projection(presence) = CRS('+proj=longlat +datum=WGS84')
  
  # circles with a radius of d
  x = circles(presence, d = samp_dist, lonlat = TRUE)
  pol = polygons(x)
  
  # sample randomly from all circles
  samp1 = spsample(pol, samp_num, type = 'random', iter = 25)
  ## Warning in proj4string(obj): CRS object has comment, which is lost in output
  # get unique cells
  cells = cellFromXY(predictors[[1]], samp1)
  
  cells_unique = unique(cells)
  
  xy = xyFromCell(predictors[[1]], cells_unique)
  
  spxy = SpatialPoints(xy, proj4string = CRS('+proj=longlat +datum=WGS84'))
  o = over(spxy, geometry(x))
  xyInside = xy[!is.na(o), ]
  
  # make into dataframe
  background = as.data.frame(xyInside)
  
  # return the result
  return(background)
  
}

# a function that creates an SWD object from a species and other parameters

create_SWD = function(
  species,
  genus,
  specific_ep,
  year,
  samp_dist,
  samp_num,
  predictors,
  thin
) {
  
  #checks
  checkmate::assert_character(species)
  checkmate::assert_character(genus)
  checkmate::assert_character(specific_ep)
  checkmate::assert_numeric(year)
  checkmate::assert_numeric(samp_dist)
  checkmate::assert_numeric(samp_num)
  checkmate::assert_flag(thin)
  
  # load data
  presence = load(genus = genus, specific_ep = specific_ep)
  
  # clean data
  presence_clean = clean_pres(presence, year, predictors, thin)
  
  # generate psuedo absences
  background = get_backgr(presence = presence_clean, samp_dist, samp_num, predictors)
  
  #build SWD
  SWD = prepareSWD(species = species, 
                              p = presence_clean, 
                              a = background, 
                              env = predictors)
  
  # return
  return(SWD)
  
}

# function that takes an SWD object and returns a 
# tuned MaxEnt model by selecting variables

model_tune = function(
  SWD_object,
  threshold,
  predictors
){
  
  # checks
  checkmate::assert_numeric(threshold)
  
  # split model into train and test
  list = trainValTest(SWD_object, test = 0.2, 
                                   only_presence = TRUE, seed = 25)
  
  # define train and test
  train = list[[1]]
  test = list[[2]]
  
  # run default model
  default_model = SDMtune::train(method = "Maxent", data = train)
  
  # tune model by removing correlated variables using jacknife removal
  reduced_variables_model <- reduceVar(default_model, th = threshold, metric = "auc", 
                                       test = test, permut = 1, use_jk = TRUE)
  
  # extract new SWD_object from model
  new_SWD = reduced_variables_model@data
  
  # recover testing data via a merge
  merged_SWD = SDMtune::mergeSWD(new_SWD, test, only_presence = TRUE)
  
  # use merged SWD to tune hyperparameters
  # split into training, validation, and testing sets
  reduced_list = trainValTest(merged_SWD, val = 0.2, test = 0.2, 
                              only_presence = TRUE, seed = 61516)
  
  # define train, validation, and test
  reduced_train = reduced_list[[1]]
  reduced_validation = reduced_list[[2]]
  reduced_test = reduced_list[[3]]
  
  # run a model
  maxent_model <- SDMtune::train("Maxent", data = reduced_train)
  
  # define hyperparameters by which to tune
  h <- list(fc = c("l", "lq", "lh", "lqp", "lqph", "lqpht"), 
            iter = seq(300, 1100, 200), 
            reg = seq(0.2, 2, 0.2))
  
  # optimize
  optimum_params = optimizeModel(maxent_model, hypers = h, metric = "auc", 
                                test = val, pop = 15, gen = 2, seed = 798)
  
  # build new model, collapsing validation back into training
  index <- which.max(optimum_params@results$test_AUC)  # Index of the best model in the experiment
  train_val <- mergeSWD(reduced_train, reduced_validation, only_presence = TRUE)
  
  final_model <- SDMtune::train("Maxent", data = train_val, 
                       fc = optimum_params@results[index, 1], 
                       reg = optimum_params@results[index, 2],
                       iter = optimum_params@results[index, 3])
  
  maxent_auc = max(optimum_params@results$test_AUC)
  
  ### run a maxnet model
  
  # run default model
  maxnet_model = SDMtune::train(method = "Maxnet", data = train)
  
  # tune model by removing correlated variables using jacknife removal
  reduced_variables_model <- reduceVar(maxnet_model, th = threshold, metric = "auc", 
                                       test = test, permut = 1, use_jk = TRUE)
  
  # extract new SWD_object from model
  new_SWD_maxnet = reduced_variables_model@data
  
  # recover testing data via a merge
  merged_SWD_maxnet = SDMtune::mergeSWD(new_SWD_maxnet, test, only_presence = TRUE)
  
  # use merged SWD to tune hyperparameters
  # split into training, validation, and testing sets
  reduced_list_maxnet = trainValTest(merged_SWD_maxnet, val = 0.2, test = 0.2, 
                              only_presence = TRUE, seed = 61516)
  
  # define train, validation, and test
  reduced_train_maxnet = reduced_list_maxnet[[1]]
  reduced_validation_maxnet = reduced_list_maxnet[[2]]
  reduced_test_maxnet = reduced_list_maxnet[[3]]
  
  # run a model
  maxnet_model_reduced <- SDMtune::train("Maxnet", data = reduced_train)
  
  # define hyperparameters by which to tune
  h <- list(fc = c("l", "lq", "lh", "lqp", "lqph", "lqpht"), 
            reg = seq(0.2, 2, 0.2))
  
  # optimize
  optimum_params_maxnet = optimizeModel(maxnet_model_reduced, hypers = h, metric = "auc", 
                                 test = val, pop = 15, gen = 2, seed = 798)
  
  # build new model, collapsing validation back into training
  index_maxnet <- which.max(optimum_params_maxnet@results$test_AUC)  # Index of the best model in the experiment
  train_val_maxnet <- mergeSWD(reduced_train_maxnet, reduced_validation_maxnet, 
                        only_presence = TRUE)
  
  final_model_maxnet <- SDMtune::train("Maxnet", data = train_val, 
                                fc = optimum_params_maxnet@results[index_maxnet, 1], 
                                reg = optimum_params_maxnet@results[index_maxnet, 2])
  
  maxnet_auc = max(optimum_params_maxnet@results$test_AUC)
  
  # test which model is better and return that one
  
  if (maxent_auc >= maxnet_auc) {
    # return tuned maxent model
    final_list = list(train_val, reduced_test, final_model)
    return(final_list)
  } else {
    # return maxnet model
    final_list = list(train_val_maxnet, reduced_test_maxnet, final_model_maxnet)
    return(final_list)
  }
  
}

# a function that get the environment for any model

get_env = function(
  model,
  predictors
){

  # get layers
  layers = names(model@data@data)
  
  # subset rasterStack
  env = raster::subset(predictors, subset = layers)
  
  # return
  return(env)
  
}

# a function that evaluates any model

evaluate_model = function(
  model_tune_out,
  predictors
){
  
  # checks
  checkmate::assert_list(model_tune_out)
  
  # define
  test = model_tune_out[[2]]
  model = model_tune_out[[3]]
  species = model@data@species
  
  # calculate aicc
  aicc = SDMtune::aicc(model, env = get_env(model, predictors))
  
  # calculate auc
  auc = SDMtune::auc(model, test = test)
  
  # calculate tss
  tss = SDMtune::tss(model, test = test)
  
  # return species and model evaluations
  # make a list
  list = list(species, auc, tss, aicc)
  
  return(list)
  
}

# retrain model with full dataset

retrain_model = function(
  model_tune_out
){
  
  # checks
  checkmate::assert_list(model_tune_out)
  
  # define
  train = model_tune_out[[1]]
  test = model_tune_out[[2]]
  model = model_tune_out[[3]]
  
  # merge
  all_occ = mergeSWD(train, test)
  
  # model class
  class = class(model@model)
  
  # if else statement to retrain model
  if (class == "Maxnet") {
    # retrain
    final_model <- SDMtune::train("Maxnet", data = all_occ, 
                                  fc = model@model@fc, 
                                  reg = model@model@reg)
    
    # return all_occ and final_model
    # make a list
    list = list(final_model, all_occ)
    
    return(list)
  } else {
    
    # retrain
    final_model <- SDMtune::train("Maxent", data = all_occ, 
                                  fc = model@model@fc, 
                                  reg = model@model@reg,
                                  iter = model@model@iter)
    
    # return all_occ and final_model
    # make a list
    list = list(final_model, all_occ)
    
    return(list)
  }
  
}

# a function that projects the model onto a certain climate

predict_model = function(
  retrained_out,
  predictors
){
  
  # checks
  checkmate::assert_list(retrained_out)
  
  # define
  model = retrained_out[[1]]
  env = get_env(model, predictors)
  species = model@data@species
  
  # predict
  # prediction
  pred = SDMtune::predict(model, data = env, type = "cloglog")
  
  ths = thresholds(model, type = "cloglog")
  ths_val = ths[2,2]
  
  reclass_df = c(0, ths_val, 0,
                  ths_val, 1, 1)
  reclass_m = matrix(reclass_df,
                      ncol = 3,
                      byrow = TRUE)
  
  # reclassify the raster using the reclass object - reclass_m
  presence_absence = reclassify(pred, reclass_m)
  
  # return binary and prediction
  # make a list
  list = list(species, pred, presence_absence)
  
  return(list)
  
}

##### Polyploid Triads #####

# a function that builds three models, one for each parent, one for hybrids

triad_SWD = function(
  par1,
  par2,
  hyb,
  genus,
  specific_ep1,
  specific_ep2,
  specific_ep_hyb,
  year,
  samp_dist,
  samp_num,
  predictors,
  thin
){
  
  # checks
  checkmate::assert_character(par1)
  checkmate::assert_character(par2)
  checkmate::assert_character(hyb)
  checkmate::assert_character(genus)
  checkmate::assert_character(specific_ep1)
  checkmate::assert_character(specific_ep2)
  checkmate::assert_character(specific_ep_hyb)
  checkmate::assert_numeric(year)
  checkmate::assert_numeric(samp_dist)
  checkmate::assert_numeric(samp_num)
  checkmate::assert_flag(thin)
  
  # build SWDs
  par1 = create_SWD(species = par1, genus = genus, specific_ep = specific_ep1, year = year,
                    samp_dist = samp_dist, samp_num = samp_num, predictors = predictors,
                    thin = thin)
  
  par2 = create_SWD(species = par2, genus = genus, specific_ep = specific_ep2, year = year,
                    samp_dist = samp_dist, samp_num = samp_num, predictors = predictors,
                    thin = thin)
  
  hyb = create_SWD(species = hyb, genus = genus, specific_ep = specific_ep_hyb, year = year,
                   samp_dist = samp_dist, samp_num = samp_num, predictors = predictors,
                   thin = thin)
  
  # create list and return
  list = list(par1, par2, hyb)
  return(list)
  
}

# a function that builds tuned SDMs for each species in the triad

triad_SDM = function(
  triad_SWD,
  threshold
){
  
  # checks
  checkmate::assert_list(triad_SWD)
  checkmate::assert_numeric(threshold)
  
  # do SDM_tune function over a list
  
  SDMs = lapply(triad_SWD, model_tune, threshold = threshold)
  
  # return result
  return(SDMs)
  
}

# evaluate models

evaluate_triad_models = function(
  triad_SDM_out,
  predictors
){
  
  # checks
  checkmate::assert_list(triad_SDM_out)
  
  # do evaluate model function over a list
  
  evals = lapply(triad_SDM_out, evaluate_model, predictors = predictors)
  
  # return
  return(evals)
  
}

# retrain models

retrain_triad = function(
  triad_SDM_out
){
  
  # checks
  checkmate::assert_list(triad_SDM_out)
  
  # run retrain function over a list
  retrained_triads = lapply(triad_SDM_out, retrain_model)
  
  # return
  return(retrained_triads)
  
}

# a function that predicts models for triads

predict_triad = function(
  retrain_triad_out,
  predictors
){
  
  # checks
  checkmate::assert_list(retrain_triad_out)
  
  # run predict function over a list
  predictions = lapply(retrain_triad_out, predict_model, predictors = predictors)
  
  # return predictions
  return(predictions)
  
}

# Identify potential hybridization zones

bin_hybrid_zone = function(
  predict_triad_out
){
  
  # checks
  checkmate::assert_list(predict_triad_out)
  
  # define
  parent1_bin = predict_triad_out[[1]][[3]]
  parent2_bin = predict_triad_out[[2]][[3]]
  
  # raster multiplication
  hybrid_zone = parent1_bin * parent2_bin
  
  return(hybrid_zone)
  
}

hyb = hybrid_zone(preds, preds1)

raster::plot(preds[[3]])
raster::plot(preds[[2]])
raster::plot(preds1[[2]])
raster::plot(preds1[[3]])
raster::plot(hyb)

pred = SDMtune::predict(retrained, data = clim_crop, type = "cloglog")
h <- list(fc = c("l", "lq", "lh", "lqp", "lqph", "lqpht"), 
          iter = seq(300, 1100, 200), 
          reg = seq(0.2, 2, 0.2))
maxent_model = SDMtune::train(method = "Maxent", data = new_SWD)
optimum_params = optimizeModel(maxent_model, hypers = h, metric = "auc", 
                               test = val, pop = 15, gen = 2, seed = 798)

# build new model, collapsing validation back into training
# Index of the best model in the experiment
index <- which.max(optimum_params@results$test_AUC)  
train_val <- mergeSWD(reduced_train, reduced_validation, only_presence = TRUE)

final_model <- SDMtune::train("Maxent", data = train_val, 
                              fc = optimum_params@results[index, 1], 
                              reg = optimum_params@results[index, 2],
                              iter = optimum_params@results[index, 3])

maxent_auc = max(optimum_params@results$test_AUC)
index

class = class(model_test9[[3]]@model)
class
##### Predictions Through Time + Niche Comparisons #####
##### Climate Stability #####

##### Niche Distance from Parental Hybrid Zones #####

# schoener's D statistic
