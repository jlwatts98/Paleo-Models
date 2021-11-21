##### Idea #####

# calculate hybridization viability in different times by using raster mathematics
# calculate climate stability for each time zone
# calculate probability of occurrence for each species for each time zone
# calculate climatic distance from hybrid for each pixel of each time zone
# multiply all these rasters to get a hybridization index
# figure out hybridization through time for each triad
##### Building Present Models #####
# use presence data from Suissa, Sundue, and Testo
# https://onlinelibrary-wiley-com.eres.library.manoa.hawaii.edu/doi/full/10.1111/jbi.14076
# follow methods from Vila-Vicosa et al., 2020
# https://www.nature.com/articles/s41598-020-78576-9#Sec8
# other important climate resources
# http://www.paleoclim.org
# https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.03031
# https://journals.ku.edu/jbi/article/view/9786

# load fern occurences

ferns = read.csv("fern_occurences.csv")
str(ferns)
ferns$Accepted_binomial = as.factor(ferns$Accepted_binomial)

# asplenium triad
a_montanum = gbif(genus = "Asplenium", species = "montanum*")
a_platyneuron = gbif(genus = "Asplenium", species = "platyneuron*")
a_bradleyi = gbif(genus = "Asplenium", species = "bradleyi*")

# dryopteris triad
d_goldiana = gbif(genus = "Dryopteris", species = "goldiana*")
d_ludoviciana = gbif(genus = "Dryopteris", species = "ludoviciana*")
d_celsa = gbif(genus = "Dryopteris", species = "celsa*")

# load explanatory variables and country boundaries
usa = raster::getData("GADM", country = "USA", level = 0)
can = raster::getData("GADM", country= "CAN", level=0)
mex = raster::getData("GADM", country= "MEX", level=0)
row.names(usa) <- paste("usa", row.names(usa), sep="_")
row.names(can) <- paste("can", row.names(can), sep="_")
row.names(mex) <- paste("mex", row.names(mex), sep="_")
countries <- rbind(usa, can, mex)


clim_layers = raster::getData("worldclim", var="bio", res = 2.5)
#plot(clim_layers[[1]])

clim_crop = crop(clim_layers, extent(countries))
e = extent(-180, -50, 10, 90)
clim_crop = mask(clim_crop, countries)
clim_crop = crop(clim_crop, e)
clim_crop = raster::stack(clim_crop)

plot(clim_crop[[1]])
names(clim_crop)

##### Soil + Elevation layers #####
# dem downloaded from here
#https://www.worldclim.org/data/worldclim21.html
library(raster)

dem30 = raster("elevation_30.tif")
dem2.5 = raster("elevation_2.5.tif")

# topographic ruggedness index
library(spatialEco)
tri2.5 = tri(dem2.5)
plot(tri2.5)

# topographic wetness index
library(dynatopmodel)
twi2.5 = upslope.area(dem2.5, atb = TRUE)
plot(twi2.5[[2]])
layerStats(twi2.5, stat = "cov")

# stack layers
dem_layers = addLayer(twi2.5, dem2.5)
dem_layers = addLayer(dem_layers, tri2.5)

# crop everything
dem_crop = crop(dem_layers, extent(countries))
e = extent(-180, -50, 10, 90)
dem_crop = mask(dem_crop, countries)
dem_crop = crop(dem_crop, e)
dem_crop = raster::stack(dem_crop)

plot(dem_crop[[3]])
names(dem_crop) = c("upslope", "twi", "dem", "tri")
names(dem_crop)

extent(dem_crop)
clim_crop = crop(clim_crop, e)
extent(clim_crop)

dem_crop_test = extend(dem_crop, clim_crop)
dem_crop = stack(dem_crop_test)
res(dem_crop)
res(clim_crop)
dem_crop_test = resample(dem_crop, clim_crop)
dem_crop = dem_crop_test

# add topographic data to clim_crop

clim_crop_test = stack(clim_crop, dem_crop)
names(clim_crop_test)
plot(clim_crop_test[[20]])


#soilph = raster("sdat_1242_6_20210709_204019824.tif")
#plot(soilph)
#weights = matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), nrow = 5)
#soilph_focal = raster::focal(soilph, w = weights, fun = mean, NAonly = T, na.rm = T)
#plot(soilph_focal)

#first import all files in a single folder as a list 
rastlist <- list.files(path = "/Users/muirlab/Documents/Jacob/Fern Hybrid Paleo SDM/NACP_MSTMIP_UNIFIED_NA_SOIL_MA_1242/data", pattern='.TIF$', all.files=TRUE, full.names=FALSE)

library(raster)
allrasters <- stack(rastlist)
# aggregate to a larger resolution
bio1 = clim_crop[[1]]
bio1_ag = aggregate(bio1, fact = 2)
e = extent(-180, -50, 10, 90)
plot(bio1_ag, ext = e)
plot(bio1, ext = e)

clim_resamp = resample(clim_crop, bio1_ag)
clim_resamp = crop(clim_resamp, e)
clim_resamp = stack(clim_resamp)


# plot with rasterVis
gplot(clim_resamp$bio1) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  scale_fill_gradientn(colours = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"),
                       na.value = "transparent",
                       name = "°C x 10") +
  labs(title = "Annual Mean Temperature",
       x = "longitude",
       y = "latitude") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())

# clean presence data
a_montanum = a_montanum %>% filter(!is.na(lon) & !is.na(lat))

# remove duplicates
a_montanum[!duplicated(a_montanum[c("lon","lat")]),]

# coordinate cleaner package
flags = clean_coordinates(a_montanum, lon = "lon", lat = "lat",
                         species = "species", 
                         tests = c("capitals", "centroids", "equal", 
                                  "gbif", "institutions", "outliers", 
                                  "seas", "zeros"))
summary(flags)

#Exclude problematic records
a_montanum = a_montanum[flags$.summary,]

# basis of record
table(a_montanum$basisOfRecord)

# year
table(a_montanum$year)

# remove records before 1900
a_montanum = a_montanum %>% filter(year > 1900)

str(a_montanum)

# thin
a_montanum_thin = thin(a_montanum, lat.col = "lat", long.col = "lon",
                       spec.col = "species", thin.par = 10,
                       reps = 1, write.files = FALSE)

# plot
wm = borders("world", colour="gray50", fill="gray50")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = a_montanum, aes(x = lon, y = lat),
             colour = "darkred", size = 0.5)+
  theme_bw() +
  scale_x_continuous(limits = c(-180, -50)) +
  scale_y_continuous(limits = c(10, 90))


a_montanum_sp = SpatialPointsDataFrame(a_montanum[,c("lon", "lat")], a_montanum[,1:147])


# https://rspatial.org/raster/sdm/3_sdm_absence-background.html
# selecting background points
# make into spatial
coordinates(a_montanum) = ~lon+lat
projection(a_montanum) = CRS('+proj=longlat +datum=WGS84')

# circles with a radius of 50 km
x <- circles(a_montanum, d=200000, lonlat=TRUE)
pol <- polygons(x)

# sample randomly from all circles
samp1 <- spsample(pol, 2500, type='random', iter=25)
## Warning in proj4string(obj): CRS object has comment, which is lost in output
# get unique cells
cells <- cellFromXY(clim_resamp[[1]], samp1)
length(cells)
## [1] 250
cells <- unique(cells)
length(cells)
## [1] 161
xy <- xyFromCell(clim_resamp[[1]], cells)
# plot to inspect the results
plot(pol, axes=TRUE)
points(xy, cex=0.75, pch=20, col='blue')
spxy <- SpatialPoints(xy, proj4string=CRS('+proj=longlat +datum=WGS84'))
o <- over(spxy, geometry(x))
xyInside <- xy[!is.na(o), ]
# make into dataframe
background = as.data.frame(xyInside)

# change a_montanum back to dataframe
a_montanum = as.data.frame(a_montanum)

a_montanum_xy = a_montanum[,c("lon","lat")]

# plot background
ggplot(map_data("world"), aes(long, lat)) +
  geom_polygon(aes(group = group), fill = "grey95", color = "gray40", size = 0.2) +
  geom_jitter(data = background, aes(x = x, y = y), color = "red",
              alpha = 0.4, size = 1) +
  labs(x = "longitude", y = "latitude") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_fixed() +
  scale_x_continuous(limits = c(-180, -50)) +
  scale_y_continuous(limits = c(10, 90))

# prepare SWD object
a_montanum_SWD = prepareSWD(species = "Asplenium montanum", 
                            p = a_montanum_xy, 
                            a = background, 
                            env = clim_resamp)

# visualize
a_montanum_SWD

head(a_montanum_SWD@data)
head(a_montanum_SWD@coords)

# get elevation raster and calculate topographic ruggedness and topographic wetness
# https://www.rdocumentation.org/packages/spatialEco/versions/1.3-6/topics/tri
dem = raster::getData("SRTM", lon = -150, lat = 60)
plot(dem)
dem1 = raster::getData("SRTM", lon = -100, lat = 50)
plot(dem1)
dem2 = raster::getData("SRTM", lon = -95, lat = 50)
plot(dem2)


# prepare data for modeling
## biomod2 package
# help
# https://rstudio-pubs-static.s3.amazonaws.com/47117_41f8cf66a4f24f5aa52d810ad307141f.html
# more help
# https://griffithdan.github.io/pages/outreach/SDM-Workshop-OSU-FALL2017.pdf
# video
# https://www.youtube.com/watch?v=QrwqhJgRbnY



# split species occurrence data into building and evaluation datasets
# drop NA
a_montanum = a_montanum %>% filter(!is.na(lon) & !is.na(lat))
## 75% of the sample size
smp_size = floor(0.75 * nrow(a_montanum))

## set the seed to make your partition reproducible
set.seed(1)
train_ind = sample(seq_len(nrow(a_montanum)), size = smp_size)

a_montanum_train = a_montanum[train_ind, ]
coords = a_montanum_train[, c("lon", "lat")]
a_montanum_train = SpatialPoints(coords = coords)

a_montanum_test = a_montanum[-train_ind, ]
coords = a_montanum_test[, c("lon", "lat")]
a_montanum_test = SpatialPoints(coords = coords)


a_montanum_format = BIOMOD_FormatingData(resp.var = a_montanum_train, 
                     expl.var = clim_crop, 
                     #eval.resp.var = a_montanum_test,
                     resp.name = "Asplenium montanum",
                     PA.nb.rep = 1,
                     PA.nb.absences = 6000,
                     PA.strategy = "random",
                     na.rm = TRUE)

a_montanum_options = BIOMOD_ModelingOptions()

a_montanum_model = BIOMOD_Modeling(data = a_montanum_format, 
                                   models = c("GLM", "GAM","ANN","RF"),
                                   models.options = a_montanum_options
                                   )
a_montanum_model

# evaluate the model
a_montanum_eval = get_evaluations(a_montanum_model)
a_montanum_eval


# project the model onto present climate

a_montanum_projection = BIOMOD_Projection(modeling.output = a_montanum_model,
                                  new.env = stack(clim_layers), # modern environment
                                  proj.name = 'current',
                                  selected.models = 'all',
                                  binary.meth = 'TSS',
                                  compress = 'xz',
                                  clamping.mask = F,
                                  output.format = '.grd')
plot(a_montanum_projection)


# test
## Not run: 
# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)

# the name of studied species
myRespName <- c("ConnochaetesGnou", "GuloGulo", "PantheraOnca", 
                "PteropusGiganteus", "TenrecEcaudatus", "VulpesVulpes")

# the presence/absences data for our species 
myResp <- DataSpecies[,myRespName]

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]

multiple.plot(Data = myResp,
              coor = myRespXY )

## End(Not run)

##### SDMtune #####
#https://consbiol-unibern.github.io/SDMtune/articles/articles/make_predictions.html
default_model = SDMtune::train(method = "Maxent", data = SWD_test)
default_model
slotNames(default_model)
slotNames(default_model@data)
slotNames(default_model@model)

# change the settings
default_model <- SDMtune::train(method = "Maxent", data = SWD_test, fc = "lqph", reg = 1, iter = 500)

model <- SDMtune::train(method = "Maxent", data = SWD_test, fc = "lh", reg = 0.5, iter = 700)


# maxnet model
maxnet_model <- SDMtune::train("Maxnet", data = SWD_test)

# prediction
pred <- SDMtune::predict(default_model, data = SWD_test, type = "cloglog")

map <- predict(default_model, data = clim_resamp, type = "cloglog")

# plot
plotPred(map)

plotPred(map, lt = "Habitat\nsuitability",
         colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))

thresholds(default_model, type = "cloglog")

plotPA(map, th = 0.314)

# evaluate the model

SDMtune::auc(default_model)

SDMtune::plotROC(default_model)

SDMtune::tss(default_model)

SDMtune::aicc(default_model, env = clim_resamp)

# testing and training

c(train, test) %<-% trainValTest(SWD_test, test = 0.2, only_presence = TRUE, seed = 25)
maxnet_model <- SDMtune::train("Maxnet", data = train)

# evaluate
cat("Training auc: ", SDMtune::auc(maxnet_model))
cat("Testing auc: ", SDMtune::auc(maxnet_model, test = test))

# test 
output <- data.frame(matrix(NA, nrow = 10, ncol = 3)) # Create an empty data.frame
colnames(output) <- c("seed", "trainAUC", "testAUC")
set.seed(25)
seeds <- sample.int(1000, 10) # Create 10 different random seeds
for (i in 1:length(seeds)) { # Loop through the seeds
  c(train, test) %<-% trainValTest(SWD_test, test = 0.2, seed = seeds[i], only_presence = TRUE) # Make the train/test split
  m <- SDMtune::train("Maxnet", data = train) # train the model
  # Populate the output data.frame
  output[i, 1] <- seeds[i]
  output[i, 2] <- SDMtune::auc(m)
  output[i, 3] <- SDMtune::auc(m, test = test)
}

output

# cross validation

folds <- SDMtune::randomFolds(SWD_test, k = 4, only_presence = TRUE, seed = 25)

# 4 fold cross validation using the Maxnet model
cv_model <- SDMtune::train("Maxnet", data = SWD_test, folds = folds)
cv_model

# AUC
cat("Training AUC: ", SDMtune::auc(cv_model))
#> Training AUC:  0.8734151
cat("Testing AUC: ", SDMtune::auc(cv_model, test = TRUE))
#> Testing AUC:  0.8553538

# getting the variable importance from Maxent
default_model@model@results

# change to a readable format
vi <- maxentVarImp(default_model)
vi

plotVarImp(vi[, 1:2])
# The function accepts a data.frame with 2 columns: one with the variable name
# and one with the values, so it is enough to select the first and the third
# columns from the vi data.frame
plotVarImp(vi[, c(1,3)])


c(train, test) %<-% trainValTest(SWD_test, test = 0.2, only_presence = TRUE, seed = 25)
maxnet_model <- SDMtune::train("Maxnet", data = train)

vi_maxnet <- SDMtune::varImp(maxnet_model, permut = 5)
vi_maxnet

plotVarImp(vi_maxnet)

# Compute the permutation importance
vi_maxent <- SDMtune::varImp(default_model, permut = 10)
# Print it
vi_maxent
# Compare with Maxent output
maxentVarImp(default_model)

jk <- doJk(maxnet_model, metric = "auc", test = test)
jk

plotJk(jk, type = "train", ref = SDMtune::auc(maxnet_model))

plotJk(jk, type = "test", ref = SDMtune::auc(maxnet_model, test = test))

plotResponse(maxnet_model, var = "bio1", type = "cloglog", 
             only_presence = TRUE, marginal = FALSE, rug = TRUE)

plotResponse(maxnet_model, var = "bio17", type = "logistic", 
             only_presence = TRUE, marginal = TRUE, fun = mean, color = "blue")

plotResponse(cv_model, var = "bio1", type = "cloglog", 
             only_presence = TRUE, marginal = TRUE, fun = mean, rug = TRUE)

modelReport(maxnet_model, type = "cloglog", folder = "virtual-sp", test = test, 
            response_curves = TRUE, only_presence = TRUE, jk = TRUE, env = clim_resamp)


## data driven variable selection
plotCor(SWD_test, method = "spearman", cor_th = 0.7)
corVar(SWD_test, method = "spearman", cor_th = 0.7)

selected_variables_model <- varSel(maxnet_model, metric = "auc", 
                                   test = test, bg4cor = SWD_test, 
                                   method = "spearman", 
                                   cor_th = 0.7, permut = 1)
selected_variables_model


# You need to pass the env argument for the AICc and the use_pc argument to use the percent contribution
selected_variables_model <- varSel(default_model, metric = "aicc", bg4cor = SWD_test, 
                                   method = "spearman", cor_th = 0.7, env = clim_resamp, 
                                   use_pc = TRUE)

cat("Testing AUC before: ", SDMtune::auc(maxnet_model, test = test))

reduced_variables_model <- reduceVar(maxnet_model, th = 6, metric = "auc", test = test, 
                                     permut = 1)

cat("Testing AUC after: ", SDMtune::auc(reduced_variables_model, test = test))


cat("Testing AUC before: ", SDMtune::auc(maxnet_model, test = test))

reduced_variables_model <- reduceVar(maxnet_model, th = 15, metric = "auc", 
                                     test = test, permut = 1, use_jk = TRUE)

cat("Testing AUC after: ", SDMtune::auc(reduced_variables_model, test = test))

# You need to pass TRUE to the test argument
selected_variables_model <- SDMtune::reduceVar(cv_model, th = 6, metric = "tss", 
                                      test = TRUE, permut = 1)

# tuning hyperparameters

library(zeallot)  # For unpacking assignment
c(train, val, test) %<-% trainValTest(SWD_test, val = 0.2, test = 0.2, 
                                      only_presence = TRUE, seed = 61516)
cat("# Training  : ", nrow(train@data))
#> # Training  :  5240
cat("# Validation: ", nrow(val@data))
#> # Validation:  5080
cat("# Testing   : ", nrow(test@data))

model <- SDMtune::train("Maxnet", data = train)

# Define the values for bg
h <- list(reg = seq(0.2, 1, 0.1))
# Call the gridSearch function
exp_1 <- gridSearch(model, hypers = h, metric = "auc", test = val)

exp_1

SDMtune::plot(exp_1, title = "Experiment 1")

plot(exp_1, title = "Experiment 1", interactive = TRUE)

exp_1@results

exp_1@results[order(-exp_1@results$test_AUC), ]

# Define the values for reg
h <- list(reg = 1:4)
# Call the gridSearch function
exp_2 <- gridSearch(model, hypers = h, metric = "tss", test = val)

# Define the values for fc
h <- list(fc = c("l", "lq", "lh", "lqp", "lqph", "lqpht"))
# Call the gridSearch function
exp_3 <- gridSearch(model, hypers = h, metric = "auc", test = val)


maxent_model <- SDMtune::train("Maxent", data = SWD_test)
# Define the values for fc
h <- list("iter" = seq(300, 1100, 200))
# Call the gridSearch function
exp_4 <- gridSearch(maxent_model, hypers = h, metric = "auc", test = val)

SDMtune::get_tunable_args(maxent_model)

h <- list(reg = seq(0.2, 2, 0.2), fc = c("l", "lq", "lh", "lqp", "lqph", "lqpht"))
exp_5 <- gridSearch(model, hypers = h, metric = "auc", test = val)

# random search
h <- list(reg = seq(0.2, 5, 0.2), fc = c("l", "lq", "lh", "lp", "lqp", "lqph"))
exp_6 <- randomSearch(model, hypers = h, metric = "auc", test = val, pop = 10, seed = 65466)

exp_7 <- optimizeModel(model, hypers = h, metric = "auc", test = val, pop = 15, 
                       gen = 2, keep_best = 0.4, keep_random = 0.2, mutation_chance = 0.4, seed = 798)

# merge training and validation datasets
merged_data <- mergeSWD(train, val)

final_model <- SDMtune::train("Maxnet", data = merged_data, fc = exp_6@results[1, 1], 
                     reg = exp_6@results[1, 2])

SDMtune::auc(final_model, test = test)

map <- predict(final_model, data = clim_resamp, type = "cloglog")
plotPred(map, lt = "Habitat\nsuitability",
         colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))


##### SDMtune example code #####
library(SDMtune)

# Set general seed for all experiments

seed = 186546

## Load and prepare data--------------------------------------------------------

files = list.files(path = file.path(system.file(package = "dismo"), "ex"), pattern = "grd", full.names = TRUE)

predictors = raster::stack(files)

p_coords = virtualSp$presence

a_coords = virtualSp$absence

data = prepareSWD(species = "Virtual species", p = p_coords, a = a_coords, env = predictors[[1:8]])

folds = randomFolds(data, k = 10, seed = seed)

## ANN experiment 1200 hyperparameters------------------------------------------

# Train starting model ==> size inner layer = number of variables

set.seed(seed)

model_ann = SDMtune::train("ANN", data = data, size = 8, folds = folds)

SDMtune::auc(model_ann)

SDMtune::auc(model_ann, test = TRUE)

# 1200 hyperparameters' combinations

h_ann = list(size = 2:81, decay = c(0.01, 0.05, 0.1, 0.3, 0.5),
              
              maxit = c(100, 500, 1000))

nrow(expand.grid(h_ann)) == 1200 # Make sure there are 1200 combinations

# Genetic Algorithm

om_ann = optimizeModel(model_ann, hypers = h_ann, metric = "auc", seed = seed)

om_ann@results[1:5,]

# Grid Search

set.seed(seed)

gs_ann = gridSearch(model_ann, hypers = h_ann, metric = "auc", save_models = FALSE)

head(gs_ann@results[order(-gs_ann@results$test_AUC),])

## BRT experiment 1200 hyperparameters------------------------------------------

# Train starting model

set.seed(seed)

model_brt = train("BRT", data = data, folds = folds)

auc(model_brt)

auc(model_brt, test = TRUE)

# 1200 hyperparameters' combinations

h_brt = list(n.trees = seq(40, 1020, 20), interaction.depth = 1:4, shrinkage = seq(0.05, 0.1, 0.01))

nrow(expand.grid(h_brt)) == 1200 # Make sure there are 1200 combinations

# Genetic Algorithm

om_brt = optimizeModel(model_brt, hypers = h_brt, metric = "auc", seed = seed)

om_brt@results[1:5,]

# Grid Search

gs_brt = gridSearch(model_brt, hypers = h_brt, metric = "auc", save_models = FALSE)

head(gs_brt@results[order(-gs_brt@results$test_AUC),])

## RF experiment 1200 hyperparameters------------------------------------------

# Train starting model

set.seed(seed)

model_rf = train("RF", data = data, folds = folds)

auc(model_rf)

auc(model_rf, test = TRUE)

# 1200 hyperparameters' combinations

h_rf = list(ntree = seq(420, 1000, 20), mtry = 3:6, nodesize = 1:10)

nrow(expand.grid(h_rf)) == 1200 # Make sure there are 1200 combinations

# Genetic Algorithm

om_rf = optimizeModel(model_rf, hypers = h_rf, metric = "auc", seed = seed)

om_rf@results[1:5,]

# Grid Search

gs_rf = gridSearch(model_rf, hypers = h_rf, metric = "auc", save_models = FALSE)

head(gs_rf@results[order(-gs_rf@results$test_AUC),])

## Maxnet experiment 1200 hyperparameters---------------------------------------

# Train starting model

bg_coords = virtualSp$background

data = prepareSWD(species = "Virtual species", p = p_coords, a = bg_coords, env = predictors[[1:8]])

folds = randomFolds(data, k = 10, only_presence = TRUE, seed = seed)

model_mx = train("Maxnet", data = data, folds = folds)

auc(model_mx)

auc(model_mx, test = TRUE)

# 1200 hyperparameters' combinations

h_mx = list(reg = seq(0.1, 4.88, 0.02), fc = c("l", "lh", "lqp", "lqph", "lqpht"))

nrow(expand.grid(h_mx)) == 1200 # Make sure there are 1200 combinations

# Genetic Algorithm

om_mx = optimizeModel(model_mx, hypers = h_mx, metric = "auc", seed = seed)

om_mx@results[1:5,]

# Grid Search

gs_mx = gridSearch(model_mx, hypers = h_mx, metric = "auc", save_models = FALSE)

head(gs_mx@results[order(-gs_mx@results$test_AUC),])



# split model into train and test
list = trainValTest(SWD_test, test = 0.2, 
                    only_presence = TRUE, seed = 25)

# define train and test
train = list[[1]]
test = list[[2]]

# run default model
default_model = SDMtune::train(method = "Maxent", data = train)

# tune model by removing correlated variables using jacknife removal
reduced_variables_model <- reduceVar(default_model, th = 5, metric = "auc", 
                                     test = test, permut = 1, use_jk = FALSE)

# extract new SWD_object from model
new_SWD = reduced_variables_model@data

# recover testing data via a merge
merged_SWD = SDMtune::mergeSWD(new_SWD, test, only_presence = TRUE)
# split into training, validation, and testing sets
reduced_list = trainValTest(merged_SWD, val = 0.2, test = 0.2, 
                            only_presence = TRUE, seed = 61516)

# define train, validation, and test
reduced_train = reduced_list[[1]]
reduced_validation = reduced_list[[2]]
reduced_test = reduced_list[[3]]

# run a model
maxent_model <- train("Maxent", data = reduced_train)

presence_thin = SDMtune::thinData(a_montanum, clim_crop[[1]], "lon", "lat")

plot(clim_crop[[1]])

