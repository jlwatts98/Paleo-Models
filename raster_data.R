##### Downloading and cleaning raster data #####

source("Header.R")
# define NZ boundaries
nz = raster::getData("GADM", country = "NZL", level = 0)

# download climate data
clim_layers = raster::getData("worldclim", var="bio", res = 2.5)

# download future climate data
# http://worldclim.com/cmip5_2.5m
# CNRM-CM5 2070, rcp 85
fut_bio1 = raster("cn85bi70/cn85bi701.tif")
fut_bio2 = raster("cn85bi70/cn85bi702.tif")
fut_bio3 = raster("cn85bi70/cn85bi703.tif")
fut_bio4 = raster("cn85bi70/cn85bi704.tif")
fut_bio5 = raster("cn85bi70/cn85bi705.tif")
fut_bio6 = raster("cn85bi70/cn85bi706.tif")
fut_bio7 = raster("cn85bi70/cn85bi707.tif")
fut_bio8 = raster("cn85bi70/cn85bi708.tif")
fut_bio9 = raster("cn85bi70/cn85bi709.tif")
fut_bio10 = raster("cn85bi70/cn85bi7010.tif")
fut_bio11 = raster("cn85bi70/cn85bi7011.tif")
fut_bio12 = raster("cn85bi70/cn85bi7012.tif")
fut_bio13 = raster("cn85bi70/cn85bi7013.tif")
fut_bio14 = raster("cn85bi70/cn85bi7014.tif")
fut_bio15 = raster("cn85bi70/cn85bi7015.tif")
fut_bio16 = raster("cn85bi70/cn85bi7016.tif")
fut_bio17 = raster("cn85bi70/cn85bi7017.tif")
fut_bio18 = raster("cn85bi70/cn85bi7018.tif")
fut_bio19 = raster("cn85bi70/cn85bi7019.tif")

# download past data
# last glacial maximum 22,000 years ago
# https://www.worldclim.org/data/v1.4/paleo1.4.html#google_vignette
# CCSM4

past_bio1 = raster("cclgmbi_2-5m/cclgmbi1.tif")
past_bio2 = raster("cclgmbi_2-5m/cclgmbi2.tif")
past_bio3 = raster("cclgmbi_2-5m/cclgmbi3.tif")
past_bio4 = raster("cclgmbi_2-5m/cclgmbi4.tif")
past_bio5 = raster("cclgmbi_2-5m/cclgmbi5.tif")
past_bio6 = raster("cclgmbi_2-5m/cclgmbi6.tif")
past_bio7 = raster("cclgmbi_2-5m/cclgmbi7.tif")
past_bio8 = raster("cclgmbi_2-5m/cclgmbi8.tif")
past_bio9 = raster("cclgmbi_2-5m/cclgmbi9.tif")
past_bio10 = raster("cclgmbi_2-5m/cclgmbi10.tif")
past_bio11 = raster("cclgmbi_2-5m/cclgmbi11.tif")
past_bio12 = raster("cclgmbi_2-5m/cclgmbi12.tif")
past_bio13 = raster("cclgmbi_2-5m/cclgmbi13.tif")
past_bio14 = raster("cclgmbi_2-5m/cclgmbi14.tif")
past_bio15 = raster("cclgmbi_2-5m/cclgmbi15.tif")
past_bio16 = raster("cclgmbi_2-5m/cclgmbi16.tif")
past_bio17 = raster("cclgmbi_2-5m/cclgmbi17.tif")
past_bio18 = raster("cclgmbi_2-5m/cclgmbi18.tif")
past_bio19 = raster("cclgmbi_2-5m/cclgmbi19.tif")

# stack
fut_clim_layers = raster::stack(fut_bio1, fut_bio2,
                                fut_bio3, fut_bio4,
                                fut_bio5, fut_bio6,
                                fut_bio7, fut_bio8,
                                fut_bio9, fut_bio10,
                                fut_bio11, fut_bio12,
                                fut_bio13, fut_bio14,
                                fut_bio15, fut_bio16,
                                fut_bio17, fut_bio18,
                                fut_bio19)

# stack
past_clim_layers = raster::stack(past_bio1, past_bio2,
                                past_bio3, past_bio4,
                                past_bio5, past_bio6,
                                past_bio7, past_bio8,
                                past_bio9, past_bio10,
                                past_bio11, past_bio12,
                                past_bio13, past_bio14,
                                past_bio15, past_bio16,
                                past_bio17, past_bio18,
                                past_bio19)

# import dem from bioclim
# https://www.worldclim.org/data/worldclim21.html
dem = raster("dem25.tif")

# crop for nz
# -48, 165, -33, 180
e = extent(165, 180, -48, -33)
clim_crop = crop(clim_layers, extent(nz))
clim_crop = mask(clim_crop, nz)
clim_crop = crop(clim_crop, e)
clim_crop = raster::stack(clim_crop)

past_clim_crop = crop(past_clim_layers, extent(nz))
past_clim_crop = mask(past_clim_crop, nz)
past_clim_crop = crop(past_clim_crop, e)
past_clim_crop = raster::stack(past_clim_crop)

fut_clim_crop = crop(fut_clim_layers, extent(nz))
fut_clim_crop = mask(fut_clim_crop, nz)
fut_clim_crop = crop(fut_clim_crop, e)
fut_clim_crop = raster::stack(fut_clim_crop)

dem_crop = crop(dem, extent(nz))
dem_crop = mask(dem_crop, nz)
dem_crop = crop(dem_crop, e)

# topographic ruggedness index
library(spatialEco)
tri = tri(dem_crop)
plot(tri)

# stack everything
dem_layers = addLayer(tri, dem_crop)
names(dem_layers) = c("tri", "dem")
names(fut_clim_crop) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7",
                         "bio8", "bio9", "bio10", "bio11", "bio12", "bio13",
                         "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
names(past_clim_crop) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7",
                         "bio8", "bio9", "bio10", "bio11", "bio12", "bio13",
                         "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
predictors = stack(clim_crop, dem_layers)
fut_predictors = stack(fut_clim_crop, dem_layers)
past_predictors = stack(past_clim_crop, dem_layers)

# save as an RDS object
readr::write_rds(predictors, "objects/predictors.rds")
readr::write_rds(fut_predictors, "objects/fut_predictors.rds")
readr::write_rds(past_predictors, "objects/past_predictors.rds")
