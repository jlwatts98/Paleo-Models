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