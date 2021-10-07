##### Downloading Soil Data #####
# https://rpubs.com/ials2un/soilgrids_webdav
library(devtools)
install.packages('XML')
install.packages('rgdal')
install.packages('gdalUtils')
install.packages('sf')
install.packages('dplyr')
install_github("envirometrix/landmap")
install.packages('leaflet')
install.packages('mapview')

devtools::install_github("https://github.com/be-marc/leaflet.opacity", dependencies=TRUE)

library(XML)
library(rgdal)
library(gdalUtils)
library(raster)
library(sf)
library(dplyr)
library(RColorBrewer)
library(leaflet.opacity)
library(leaflet)
library(leaflet.opacity)
library(mapview)

url = "https://files.isric.org/soilgrids/latest/data/" # Path to the webDAV data.

voi = "nitrogen" # variable of interest
depth = "15-30cm"
quantile = "mean"

variable = paste(url, voi, sep="")

layer = paste(variable,depth,quantile, sep="_") # layer of interest 

vrt_layer = paste(layer, '.vrt', sep="")

nitro = raster("https://files.isric.org/soilgrids/latest/data/nitrogen/nitrogen_15-30cm_mean.vrt")
nitro1 = raster(vrt_layer)
raster::plot(nitro)
nitro
options(rasterMaxMemory = 1e11)
nitro_resamp = resample(nitro, clim_crop[[1]])
nitro_agg = aggregate(nitro, fact=5)

# this raster seems to be too big to work with this computer

library(sf)
sp_df <- readOGR("/Users/muirlab/Documents/Jacob/Fern Hybrid Paleo SDM/DSMW/DSMW.shp") # reads the shapefile into a large polygon
# data frame

extent(sp_df) # shows the extent of sp_df

summary(sp_df) # simple summary

# Object of class SpatialPolygonsDataFrame
# Coordinates:
#        min      max
# x 32539350 32544695
# y  5732607  5735439
# Is projected: TRUE 
# proj4string :
# [+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=32500000 +y_0=0 +ellps=GRS80+units=m+no_defs]
# Data attributes:
# X_CENTR           Y_CENTR            ATTRIBUTE 
# Min.   :3539714   Min.   :5734725   9.4-.2.3 :21  
# 1st Qu.:3540436   1st Qu.:5735323   23.4-.2.3:19
# Median :3541226   Median :5735830   9.4.2.3  :19  
# Mean   :3541263   Mean   :5735824   23.4.3.5 :18  
# 3rd Qu.:3542031   3rd Qu.:5736358   23.4.2.3 :16  
# Max.   :3543461   Max.   :5737024   19.4.3.5 :12  
#                                     (Other)  :89  

r <- raster(extent(sp_df)) # creates a raster with the extent of sp_df
projection(r) <- proj4string(sp_df) # uses the projection of the shapefile
# for the raster    
res(r)=10000 # sets the resolution of the raster to 10 m

r10 <- rasterize(sp_df, field="SNUM", r) # converts sp_df to a raster
# assigning the last attribute
# variable in a cell (default)

writeRaster(r10, filename="stok_subs.tif", format="GTiff") # writes the raster
# to a file

plot(r10) # plots raster
