##### Test Push File #####



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
class(a_montanum)
class(a_montanum$lat)
flags = clean_coordinates(test, lon = "lon", lat = "lat",
                          species = "species", 
                          tests = c("capitals", "centroids", "equal", 
                                    "gbif", "institutions", "outliers", 
                                    "seas", "zeros"))
