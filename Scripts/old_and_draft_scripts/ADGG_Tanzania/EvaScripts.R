# Read phenotype data
setwd("/Users/marialse/Box Sync/RoslinProjects/HerdDataProject/RealHerdData/BSW_data/")


phenoData = read.csv(file = "merit_name.csv")
names(phenoData)
dim(phenoData)
# The phenotype is last004. Use vpliv1, vpliv2 and vpliv3 as fixed effects, idlok as iid random effect, idziva as random effect with relationship matrix from pedigree.

phenoData$last004scaled = scale(phenoData$last004)





# Build relationship matrix from pedigree
ped = read.csv("sorod.csv")
head(ped)


# Relationship matrix from pedigree should be built from idziva, idocea and idmata
ped2 = data.frame(label = ped$idziva, dam = ped$idmata, sire = ped$idocea)
library("pedigree")

ped3 = editPed(sire = ped2$sire, dam = ped2$dam, label = ped2$label)

# Pedigree object for pedigreemm functions
library(package="pedigreemm")
pedPMM = pedigree(label = ped3$label,dam=ped3$dam,sire=ped3$sire)



# Precision matrix (A inverse)
Tinv    = as(pedPMM, "sparseMatrix") ## T^{-1} in A^{-1} = (T^{-1})' D^{-1} T^{-1}
DmatTemp = pedigreemm::Dmat(pedPMM)
D       = Diagonal(x=DmatTemp)   ## D      in A = TDT'
Dinv    = solve(D)                   ## ...
Ainv = t(Tinv) %*% Dinv %*% Tinv  ## ...

# Labelling of the A matrix is according to pedPMM
labelAinv = as.numeric(pedPMM@label)
mappingAinv = data.frame(idziva = labelAinv, rowNumberAinv = 1:length(labelAinv))

phenoData = merge(phenoData, mappingAinv)

# # tiny test example 
# pedX = head(ped2)
# pedXedit = editPed(sire = pedX$sire, dam = pedX$dam, label = pedX$label)
# pedY = pedigree(label = pedXedit$label,dam=pedXedit$dam,sire=pedXedit$sire)
# Tinv    = as(pedY, "sparseMatrix") ## T^{-1} in A^{-1} = (T^{-1})' D^{-1} T^{-1}
# DmatTemp = pedigreemm::Dmat(pedY)
# D       = Diagonal(x=DmatTemp)   ## D      in A = TDT'
# Dinv    = solve(D)                   ## ...
# Ainv = t(Tinv) %*% Dinv %*% Tinv  ## ...





# Read spatial data
#install.packages("rgdal")
library(package = "rgdal")
shape = readOGR(dsn = ".", layer = "LOK_KO_UE_OB")
spatialData = data.frame(idlok = shape$IDLOK, xpos = shape$X/1000, ypos = shape$Y/1000, UE = shape$UE_MID, OB = shape$OB_MID, KO = shape$SIFKO)
plot(spatialData$xpos, spatialData$ypos)


# Add locations to phenoData
library(tidyverse)
phenoData = merge(phenoData, spatialData)


head(phenoData)





plot(phenoData$xpos, phenoData$ypos)
slovenia = readOGR("Slovenija/Boundary_Lines.shp")
slovenia@proj4string


library(rworldmap)
# newMap <- getMap(resolution="low")
# prec_stations_geographical<-SpatialPoints(newMap@polygons[[205]]@Polygons[[1]]@coords, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# prec_utmLow<-spTransform(prec_stations_geographical,CRS("+proj=tmerc +lat_0=0 +lon_0=15 +k=0.9999 +x_0=500000 +y_0=-5000000 +ellps=GRS80 +units=m +no_defs"))
# sloveniaBoundary = cbind(prec_utm@coords[,1]/1000, prec_utm@coords[,2]/1000)

# 
newMap <- getMap(resolution="high")
prec_stations_geographical<-SpatialPoints(newMap@polygons[[212]]@Polygons[[1]]@coords, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
prec_utm<-spTransform(prec_stations_geographical,CRS("+proj=tmerc +lat_0=-0.01 +lon_0=14.97 +k=0.9999 +x_0=500000 +y_0=-5000000 +ellps=GRS80 +units=m +no_defs"))
sloveniaBoundary = cbind(prec_utm@coords[,1]/1000, prec_utm@coords[,2]/1000)
par(mfrow = c(1,1))
plot(x = prec_utm@coords[,1]/1000, prec_utm@coords[,2]/1000,type="o")
#lines(x = prec_utmLow@coords[,1]/1000, prec_utmLow@coords[,2]/1000,type="o", col = 2)
points(shape$X/1000, shape$Y/1000, col = 3)
# 


setwd("../ProcessedData/")
save(phenoData, Ainv,sloveniaBoundary,ped2, file = "ProcessedPhenoData.RData" )


# 
# "+proj=utm +zone=33 +datum=WGS84"
# 
# 
# plot(prec_stations_geographical)
# points(phenoData$xpos/10, phenoData$ypos/10, col = "green")
# 
# 
# 
# 
# newMap <- getMap(resolution="low")
# prec_stations_geographical<-SpatialPoints(newMap@polygons[[17]]@Polygons[[1]]@coords, proj4string=CRS("+proj=longlat +ellps=intl +units=m +towgs84=-87,-98,-121 +nodefs"))
# prec_utm<-spTransform(prec_stations_geographical,CRS("+proj=lcc +lat_1=49 +lat_2=46 +lat_0=47.5 +lon_0=13.33333333333333 +x_0=400000 +y_0=400000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs"))
 
# 


# Add map of Slovenia to the herd locations 
# library(rworldmap)
# 
# newMap <- getMap(resolution="low")
# prec_stations_geographical<-SpatialPoints(newMap@polygons[[205]]@Polygons[[1]]@coords, proj4string=CRS("+proj=longlat +ellps=intl +units=m +towgs84=-87,-98,-121 +nodefs"))
# prec_utm<-spTransform(prec_stations_geographical,CRS("+proj=lcc +lat_1=49 +lat_2=46 +lat_0=47.5 +lon_0=13.33333333333333 +x_0=400000 +y_0=400000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs"))
# 
# 
# prec_utm<-spTransform(prec_stations_geographical,CRS(paste0("+proj=utm +zone=",33," +datum=WGS84")))
# lines(prec_utm@coords)
# 
# 

