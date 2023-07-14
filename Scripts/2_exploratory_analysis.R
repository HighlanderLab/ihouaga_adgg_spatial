# ---- Header ------------------------------------------------------------------

# Exloratory analysis of this data - find key patterns and errors - these error
#   correcting codes have then been ported to the quality control section of the
#   1_first_look_at_the_data_and_quality_control.R script
# Author: Isidore
# Version 1.0.0
# Date: 2023-07-13

# ---- Setup -------------------------------------------------------------------

# Working directory
# ... Isidore's laptop
baseDir <- ""
# ... Isidore's Eddie workspace
baseDir <- ""

# ... path on ILRI hpc
baseDir <- ""

# ... Gregor's office computer
baseDir <- "/Users/ggorjanc/Storages/GitBox/HighlanderLab/ihouaga_adgg_spatial"

# Change working directory
setwd(dir = baseDir)
dir()

# ---- Installing and load packages --------------------------------------------

if (FALSE) {
  requiredPackages <- c(
    "tidyverse", # for data manipulation
  )
  install.packages(pkgs = requiredPackages)
}
library(tidyverse)

(.packages()) # Check loaded packages

# ---- Import data -------------------------------------------------------------

# Clear the environment
rm(list = ls())

# Read in the data
data1 <- read.csv(file = "data/cleaned_data/milk_yield_pheno_cleaned.csv")
str(data1)

# Factors
data1 <- data1 %>%
  mutate_at(.vars = c("cow", "ward", "herd", "cyrsn", "tyrmn", "dgrp", "lac",
                      "lacgr"),
            .funs = as.factor)
str(data1)

# ---- Descriptive stats -------------------------------------------------------

summary(data1$milk)
hist(data1$milk, breaks = 0:41)
# tail perhaps a bit too long, but likely OK!

# ---- ... by lactation --------------------------------------------------------

summary(data1$lac)
tmp <- data1 %>%
  group_by(lac) %>%
  dplyr::summarise(
    n = n(),
    mean = mean(milk, na.rm = TRUE),
    sd =  sd(milk, na.rm = TRUE)
  )
tmp

tmp %>%
  ggplot(aes(y = n, x = lac, colour = lac)) +
  geom_point(aes(label = lac)) +
  xlab("Lactation") +
  ylab("Number") +
  theme_bw()

tmp %>%
  ggplot(aes(y = mean, x = lac, colour = lac)) +
  geom_point(aes(label = lac)) +
  xlab("Lactation") +
  ylab("Mean") +
  theme_bw()

tmp %>%
  ggplot(aes(y = sd, x = lac, colour = lac)) +
  geom_point(aes(label = lac)) +
  xlab("Lactation") +
  ylab("Standard deviation") +
  theme_bw()

# ---- ... by TODO --------------------------------------------------------


table(data1$lac)
table(data2$lacgr)

length(unique(data2$herd))
#Number of rows
nrow(data2) # 19494 records
#number of variables
ncol(data2)#  16
#Number of cows
length(unique(data2$cow)) # 1906
# Milk Yield
summary(data2$milk) #min=1, median=8.00 max= 40
mean(data2$milk) #8.33
var(data2$milk) # 18.66
sd(data2$milk) # 4.32
hist(na.omit(data2$milk))

################################################################################
#Ward
################################################################################

#Number of levels of wards
length(unique(data2$ward)) # 156 wards
table(data2$ward)
table(table(data2$ward)) # number of records per ward ranged from 3 to 979
plot(data2$ward,data2$milk)

# Load library dplyr
library(dplyr)

#calculate mean and sd of milk yield  by ward
data2_ward_mean_std <- data2 %>%
  group_by(ward) %>%
  summarise_at(vars(milk), list(mean=mean, sd=sd)) %>%
  as.data.frame()

#view results
data2_ward_mean_std
#plot mean and standard deviation of milk yield by ward
ggplot(data2_ward_mean_std , aes(x=ward, y=mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3) +
  geom_point(size=2)

################################################################################
#Calving year-season
################################################################################

#Number of levels of CYS
length(unique(data2$cyrsn)) # 29 levels
table(data2$cyrsn)
table(table(data2$cyrsn)) # Ranges from 1 to 2663

# Plot mean and sd of milk yield  per calving year season level

data2_cyrsn_mean_std <- data2 %>%
  group_by(cyrsn) %>%
  summarise_at(vars(milk), list(mean=mean, sd=sd)) %>%
  as.data.frame()

#view results
data2_cyrsn_mean_std
#plot mean and standard deviation of milk yield by ward
ggplot(data2_cyrsn_mean_std , aes(x=cyrsn, y=mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3) +
  geom_point(size=1)

################################################################################
#Lactation group
################################################################################

length(unique(data2$lacgr)) # 3 levels
table(data2$lacgr)
table(table(data2$lacgr))

# Plot mean and sd of milk yield  per lactation group

data2_lacgr_mean_std <- data2 %>%
  group_by(lacgr) %>%
  summarise_at(vars(milk), list(mean=mean, sd=sd)) %>%
  as.data.frame()

#view results
data2_lacgr_mean_std
#plot mean and standard deviation of milk yield by ward
ggplot(data2_lacgr_mean_std , aes(x=lacgr, y=mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3) +
  geom_point(size=2)

################################################################################
# Herd
################################################################################
data2$herd =as.factor(data2$herd)
length(unique(data2$herd)) # 1394
table(data2$herd)
table(table(data2$herd))
# Plot mean and sd of milk yield  per herd

data2_herd_mean_std <- data2 %>%
  group_by(herd) %>%
  summarise_at(vars(milk), list(mean=mean, sd=sd)) %>%
  as.data.frame()

#view results
data2_herd_mean_std
#plot mean and standard deviation of milk yield by ward
ggplot(data2_herd_mean_std , aes(x=herd, y=mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3) +
  geom_point(size=2)

################################################################################
# Explore coordinates
################################################################################

# 1. How many records we have?

length(data1$cow)
# 19494 records

# 2. How many unique animals with records we have?
length(unique(data2$cow))
#1906
# 3. How many herds we have?
length(unique(data2$herd))
# 1394 herds
# 4. How many records per herd we have?
table(table(data2$herd))
# Number of record per herd ranged from 3 to 71

nrecords_byherd<- data2 %>%
  group_by(herd) %>%
  summarise(recordherd=n())
summary(nrecords_byherd$recordherd)

summary(nrecords_byherd$recordherd)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 3.00    7.00   11.00   13.98   18.00   71.00

# 5. How many animals per herd we have?
# Number of records per cow per herd
ncowherd<- data2 %>%
  group_by(herd,cow)%>%
  summarise(ncow_byherd= n())
summary(ncowherd$ncow_byherd)
# number of cows per herd
ncowherd_2<- ncowherd %>%
  select(herd) %>%
  summarise(ncowherd=n())
summary(ncowherd_2$ncowherd)

#Min. 1st Qu.  Median  Mean    3rd Qu.    Max.
# 3.00    6.00    9.00   10.23   14.00      37.00

# Descriptive stats of age and days in milk
summary(data2$lac)
summary(data2$age) # age in months
#Min.   1st Qu.  Median Mean    3rd Qu.    Max.
#18.05   38.79   47.86   50.29   58.88     162.10
# Histogram of age
hist(na.omit(data2$age))
head(data2, n=1)

################################################################################
# Descriptive statistics of coordinates
################################################################################
setwd("C:/Users/Lenovo/OneDrive/Documents/adgg")
getwd()
dir()
# Import data
data3<- read.table(file = "data3.dat", header = FALSE)
colnames(data3) <- c("cow","milk","ward","herd","cyrsn","tyrmn","dgrp","lac","lacgr","age","long","lat", "intercept", "leg1","leg2", "mean", "ward_code")

head(data3)
nrow(data3) # 19494
# number of herds
length(unique(data3$herd)) # 1394
dim(data3) #19494    17
dim(unique(data3[, c("long","lat")])) # 1397 2
# 1394 herds with 1397 unique longitude and latitude.The difference is 3.  This means that some herds
#are associated with 2 or 3 couples of longitude and latitude.

#19494 milk yield recorded at 1397 different locations.

#summary of coordinates
summary(data3$long)
#Min.   1st Qu.   Median    Mean    3rd Qu.    Max.
# 0.0001 35.2967  36.7329   36.2023 37.5862    39.1247
summary(data3$lat)
# Min.      1st Qu.    Median     Mean      3rd Qu.  Max.
# -9.45992 -8.30842   -5.03902   -5.52684 -3.35278  0.09999
hist(data3$long)
hist(data3$lat)
library(geoR)
library(sp)
library(rgdal)
library(sf)
library(ggplot2)
library(dplyr)
gps1<- unique(data3[, c("long","lat")])
gps1
length(gps1$long)
length(gps1$lat)
gps1_lat_long <- gps1 %>% select(lat, everything())
data4 <- data3
# Create new column longxlat
data4$longxlat  <- paste(data4$long, data4$lat, sep= "-")
head(data4)

table(table(data4$longxlat))

#  Remove wrong GPS
data5 <- subset(data4, !(long %in% c(0.0001, 0.9999)))
# New number of records
length(data5$cow) # 19310
#lost records after removing wrong GPS
19494-19310 # 184
# New number of cows
length(unique(data5$cow))# 1888
#lost number of cows after removing wrong GPS
1906-1888# 18 cows lost after removing wrong GPS
#New number of herds
length(unique(data5$herd)) # 1386
# lost herds after removing wrong GPS
1394-1386 # 8 herds
# New number of unique couple of GPS
gps2<- unique(data5[, c("long","lat")])
length(gps2$long)
length(gps2$lat)

# Identifying herds with wrong GPS
herd_Wrong_GPS <- subset(data4, (long %in% c(0.0001, 0.9999)))
summary(herd_Wrong_GPS)
length(unique(herd_Wrong_GPS$cow)) # cows
length(herd_Wrong_GPS$cow) # 184 records
summary(herd_Wrong_GPS$herd)
length(unique(herd_Wrong_GPS$herd))
herd_Wrong_GPS$herd <- as.factor(herd_Wrong_GPS$herd)
summary(herd_Wrong_GPS$herd)

#The wrong GPS
Wrong_GPS_herd<-  unique(herd_Wrong_GPS[, c("herd","long","lat")])

# herds with wrong records: 34 (Right too)   99  108  112  113  278  284  393  645 (Right too)  763(Right too)  872(Right too)  928  952(Right too) 1115 (Right too)
datax<- data5

datax$herd <- as.factor(datax$herd)
summary(datax$herd)
length(unique(data5$herd))

data_herd_mist<- subset(data5, (herd %in% c(34,645,763,872,952,1115)))

k<- sort(unique(data5$herd))
summary(k)
k<- unlist(k)
k
tail(k, n=400)

# herds with duplicated GPS (Wrong vs right)

duplicated_Wrong_RGPS_herd<-  unique(data_herd_mist[, c("herd","long","lat")])

duplicated_Wrong_RGPS_herd
length(unique(data5$ward)) #156 wards in
# 1386 herds from 1385 locations. 2 herds from same location. Very common in smallholder systems where herds from same location are treated as different because of different ownership
dim(data5)
length(unique(data5$cow))
length(unique(data5$herd))
length(unique(data5$ward))
#
data5<- read.table(file = "data5.dat", header = FALSE)
colnames(data5) <- c("cow","milk","ward","herd","cyrsn","tyrmn","dgrp","lac","lacgr","age","long","lat", "intercept", "leg1","leg2", "mean", "ward_code", "Region_code")
head(data5)
nrow(data5)
model1 <- d

############################################################################
#Mapping  new data (data5) without wrong GPS in R

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
data_wa=data5
head(data_wa)
coordinates(data_wa) = ~long+lat # Cordinates are in columns long and lat
getwd()

projection(mapwa) = CRS("+proj=longlat +datum=WGS84")
projection(data_wa) = CRS("+proj=longlat +datum=WGS84")# Map of my data using same projection(language)

# same results if that's how projection is defined
#proj4string(data) = CRS("+proj=longlat +a=6378249.145 +rf=293.465 +no_defs +type=crs")
#projection(mapwa)=CRS("+proj=longlat +a=6378249.145 +rf=293.465 +no_defs +type=crs")

return1 = over(geometry(data_wa), mapwa, returnList = F) # returnList = F gives dataframe instead list
head(return1)
return1
# binding war_code to data5
#Remove old war_code and longxlat
data5 <- data5[,c(-17,-18)]
head(data5)
dim(data5)
# Adding ward_code and Region_code
data5$ward_code <- return1$Ward_Code
data5$Region_Cod <- return1$Region_Cod
head(data5)
#Ploting data on map
mapwa2 <- st_read("/Users/Lenovo/OneDrive/Documents/adgg/TZwards.shp", quiet = TRUE)
class(mapwa2) # class sf
head(mapwa2)
plot(mapwa2)

#ggplot(mapwa2) + geom_sf()+ geom_sf_text(data = mapwa2, aes(label = Ward_Code), size = 3, colour = "black")

### Get long and lat from your data.frame. Make sure that the order is in lon/lat.
head(data5)
xy <- data5[,c(11,12)]

#Plotting the coordinates on my data to map
ggplot(data = mapwa2) +
  geom_sf() +
  geom_point(data = xy, aes(x = long, y = lat))


#Let's plot data containing wrong GPS.
x2y2 <- data3[,c(11,12)]


ggplot(data = mapwa2) +
  geom_sf() +
  geom_point(data = x2y2, aes(x = long, y = lat))
################################################################################
# Ward spatial
################################################################################
# Ward Neigbouring relationship
################################################################################

# construct graph of neighbours (Library INLA)
library(spdep) # To identify if two regions are neighbours
#remotes::install_github("r-spatial/s2")
library(s2)
sf_use_s2(FALSE)

nb.mapwa2 <- poly2nb(mapwa2)# Look at map and tell neighbours
nb2INLA("map.graphwa2",nb.mapwa2)
gwa <- inla.read.graph(filename = "map.graphwa2") # Lines 334 and 335 run together
# then: construct (inverse) of neighbour matrix

Rwa = -inla.graph2matrix(gwa) # Convert graph to matrix
dim(Rwa) #3644*3644 #
#abs(det(Rwa)) # Determinant is zero
diag(Rwa) = gwa$nnbs # diagonal = nr of wards?????
summary(gwa)
diag(Rwa)
#Scaling Rwa matrix
scaledRwa = inla.scale.model(Rwa, constr = list(A = matrix(1, 1, 3644), e=0))# scaling failed.

# triangular form
Qlowwa <- tril(scaledRwa)
triowa <- summary(Qlowwa)
summary(triowa$i)
# export that - neighbours_scaledwa.txt file -> blupf90
write.table(triowa, "neighbours_scaledwa.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
head(data_wa)

# triangular form of non-scaled
Qlowwa2 <- tril(Rwa)
triowa2 <- summary(Qlowwa2) # export that - neighbours_ns.txt file -> blupf90
summary(triowa2$i)
write.table(triowa2, "neighbours_ns.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
head(data_wa)

#Export data3 with ward_code for blupf90
write.table(data, "data3.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")

#Export data5 (no wrong GPS) with ward_code for blupf90
write.table(data5, "data5.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
colnames(data5)
#data6= data3 cows ordered  from 1
data6<- data3[order(data3$cow),]
#Export data6 for blupf90
write.table(data6, "data6.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
# Save data5 as .csv
write.csv(data5,"/Users/Lenovo/OneDrive/Documents/adgg/data5.csv")
getwd()

nrow(data5)
length(unique(data5$cow))
length(unique(data5$ward))
getwd()


