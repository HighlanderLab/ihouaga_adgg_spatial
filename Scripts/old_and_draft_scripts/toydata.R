rm(list=ls())

data = read.csv("toydata.csv")
head(data)
library(sf)
library(rgdal)

#library(raster)

mapwa = readOGR("2012_Wards_Shapefiles/2012 Wards Shapefiles/TZwards.shp")
# plot(mapwa)

coordinates(data) = ~long+lat

proj4string(mapwa) = CRS("+proj=longlat +datum=WGS84")
proj4string(data) = CRS("+proj=longlat +datum=WGS84")

# same results if that's how projection is defined
#proj4string(data) = CRS("+proj=longlat +a=6378249.145 +rf=293.465 +no_defs +type=crs")
#projection(mapwa)=CRS("+proj=longlat +a=6378249.145 +rf=293.465 +no_defs +type=crs")

return = over(geometry(data), mapwa, returnList = F)

return
getwd()
