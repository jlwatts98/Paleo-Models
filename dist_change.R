##### Overlap Current and Future Rasters #####
# source header
source("Header.R")

# load in model predictions
predictions = readr::read_rds("objects/predictions.rds")
fut_predictions = readr::read_rds("objects/fut_predictions.rds")
past_predictions = readr::read_rds("objects/past_predictions.rds")

# raster addition
sp1 = predictions[[1]][[3]]
sp2 = predictions[[2]][[3]]
sp3 = predictions[[3]][[3]]
sp4 = predictions[[4]][[3]]
sp5 = predictions[[5]][[3]]
sp6 = predictions[[6]][[3]]
sp7 = predictions[[7]][[3]]
sp8 = predictions[[8]][[3]]
sp9 = predictions[[9]][[3]]

fut_sp1 = fut_predictions[[1]][[3]] * 2
fut_sp2 = fut_predictions[[2]][[3]] * 2
fut_sp3 = fut_predictions[[3]][[3]] * 2
fut_sp4 = fut_predictions[[4]][[3]] * 2
fut_sp5 = fut_predictions[[5]][[3]] * 2
fut_sp6 = fut_predictions[[6]][[3]] * 2
fut_sp7 = fut_predictions[[7]][[3]] * 2
fut_sp8 = fut_predictions[[8]][[3]] * 2
fut_sp9 = fut_predictions[[9]][[3]] * 2

past_sp1 = past_predictions[[1]][[3]] * 2
past_sp2 = past_predictions[[2]][[3]] * 2
past_sp3 = past_predictions[[3]][[3]] * 2
past_sp4 = past_predictions[[4]][[3]] * 2
past_sp5 = past_predictions[[5]][[3]] * 2
past_sp6 = past_predictions[[6]][[3]] * 2
past_sp7 = past_predictions[[7]][[3]] * 2
past_sp8 = past_predictions[[8]][[3]] * 2
past_sp9 = past_predictions[[9]][[3]] * 2

fut_overlap1 = sp1 + fut_sp1
fut_overlap2 = sp2 + fut_sp2
fut_overlap3 = sp3 + fut_sp3
fut_overlap4 = sp4 + fut_sp4
fut_overlap5 = sp5 + fut_sp5
fut_overlap6 = sp6 + fut_sp6
fut_overlap7 = sp7 + fut_sp7
fut_overlap8 = sp8 + fut_sp8
fut_overlap9 = sp9 + fut_sp9

# add past overlap too
past_overlap1 = sp1 + past_sp1
past_overlap2 = sp2 + past_sp2
past_overlap3 = sp3 + past_sp3
past_overlap4 = sp4 + past_sp4
past_overlap5 = sp5 + past_sp5
past_overlap6 = sp6 + past_sp6
past_overlap7 = sp7 + past_sp7
past_overlap8 = sp8 + past_sp8
past_overlap9 = sp9 + past_sp9

# write objects per species
# load in evaluations
evaluations = readr::read_rds("objects/evaluations.rds")
species = evaluations$species
species[1]
# make lists
list1 = list(species[1], fut_overlap1, past_overlap1)
list2 = list(species[2], fut_overlap2, past_overlap2)
list3 = list(species[3], fut_overlap3, past_overlap3)
list4 = list(species[4], fut_overlap4, past_overlap4)
list5 = list(species[5], fut_overlap5, past_overlap5)
list6 = list(species[6], fut_overlap6, past_overlap6)
list7 = list(species[7], fut_overlap7, past_overlap7)
list8 = list(species[8], fut_overlap8, past_overlap8)
list9 = list(species[9], fut_overlap9, past_overlap9)

dist_change = list(list1, list2, list3,
                   list4, list5, list6,
                   list7, list8, list9)
# write rds
readr::write_rds(dist_change, "objects/dist_change.rds")
