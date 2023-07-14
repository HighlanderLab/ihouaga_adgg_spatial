#===========================================================================================================================================================
# ADGG data analysis
# Construct spatial graph
# Author: Isidore Houaga                                  #
# Version 1.0.0                                                        #
# Date: 25/02/2023                                                    #
#===========================================================================================================================================================
# Setting working directory
setwd("C:/Users/Lenovo/OneDrive/Documents/adgg")
getwd()
rm(list=ls())
data<- read.table(file = "data2.dat", header = TRUE)
head(data)
colnames(data) <- c("cow","milk","ward","herd","cyrsn","tyrmn","dgrp","lac","lacgr","age","long","lat", "intercept", "leg1","leg2", "mean")
head(data)
#write.csv(toydata,"/Users/Lenovo/OneDrive/Documents/adgg/toydata.csv" )

library(sf)
library(rgdal)
library(raster)
library(INLA)
library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(spData)# Required for spdep

mapwa = readOGR("/Users/Lenovo/OneDrive/Documents/adgg/TZwards.shp")
 plot(mapwa)
data_wa=data
coordinates(data_wa) = ~long+lat # Cordinates are in columns long and lat


projection(mapwa) = CRS("+proj=longlat +datum=WGS84")
projection(data_wa) = CRS("+proj=longlat +datum=WGS84")# Map of my data using same projection(language)

# same results if that's how projection is defined
#proj4string(data) = CRS("+proj=longlat +a=6378249.145 +rf=293.465 +no_defs +type=crs")
#projection(mapwa)=CRS("+proj=longlat +a=6378249.145 +rf=293.465 +no_defs +type=crs")

return = over(geometry(data_wa), mapwa, returnList = F) # returnList = F gives dataframe instead list
head(return)
return
head(return)
#library(reshape2) # Reshape the dataframe
data$ward_code <- return$Ward_Code
head(data)

################################################################################
#Neigbouring relationship
################################################################################
# https://github.com/r-spatial/sf/issues/1762
mapwa2 <- st_read("/Users/Lenovo/OneDrive/Documents/adgg/TZwards.shp", quiet = TRUE)
class(mapwa2) # class sf
head(mapwa2)
plot(mapwa)
ggplot(mapwa2) + geom_sf()+ geom_sf_text(data = mapwa2, aes(label = Ward_Code), size = 3, colour = "black")

# construct graph of neighbours (Library INLA)
library(spdep) # To identify if two regions are neighbours
#remotes::install_github("r-spatial/s2")
library(s2)
sf_use_s2(FALSE)

nb.mapwa2 <- poly2nb(mapwa2)# Look at map and tell neighbours
nb2INLA("map.graphwa2",nb.mapwa2)
gwa <- inla.read.graph(filename = "map.graphwa2") # Lines 58 and 59 run together
# then: construct (inverse) of neighbour matrix

Rwa = -inla.graph2matrix(gwa) # Convert graph to matrix
dim(Rwa) #3644*3644 #
#abs(det(Rwa)) # Determinant is zero
#diag(Rwa) = gwa$nnbs # diagonal = nr of wards?????
summary(gwa)
#diag(Rwa) = diag(Rwa)+0.001 # Add 0.001 to the diagonal
#Rwainv<- solve(Rwa) # Calculate the inverse
#Rwainv
rm(Rwa)
diag(Rwa)
head(data_wa)
#Scaling Rwa matrix
scaledRwa = inla.scale.model(Rwa, constr = list(A = matrix(1, 1, 3644), e=0))

# triangular form
Qlowwa <- tril(scaledRwa)
triowa <- summary(Qlowwa) # export that - neighbours_scaledwa2.txt file -> blupf90
summary(triowa$i)

# triangular form of non-scaled
Qlowwa2 <- tril(Rwa)
triowa2 <- summary(Qlowwa2) # export that - neighbours_scaledwa2.txt file -> blupf90
summary(triowa2$i)
write.table(triowa, "neighbours_scaledwa2.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
head(data_wa)
getwd()
#Export non sclaed neighbours matrix
write.table(triowa2, "neighbours.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
head(data_wa)

#Export data with ward_code for blupf90
write.table(data, "data3.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")







