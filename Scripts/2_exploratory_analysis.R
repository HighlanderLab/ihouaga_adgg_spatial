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
data1 <- read.csv(file = "data/cleaned_data/milk_yield_pheno_raw.csv")
str(data1)

# Factors
data1 <- data1 %>%
  mutate_at(.vars = c("cow", "ward", "herd", "cyrsn", "tyrmn", "dgrp", "lac", "ward_code"),
            .funs = as.factor)
str(data1)

# ---- Data summaries ----------------------------------------------------------

summary(data1)

# number of records
nrow(data1) # 19538

# number of cows
length(unique(data1$cow)) # 1911

# number of herds
length(unique(data1$herd)) # 1396, so out of 1400 herd codes some are not present

# number of wards
length(unique(data1$ward)) # 156, so out of 157 herd codes some are not present

# number of records per herds per ward
(tmp <- table(data1$herd, data1$ward)) # 3 to 71
(tmp2 <- table(tmp[tmp > 0]))
hist(tmp2)

# number of herds per ward
tmp <- data1 %>%
  group_by(ward, herd) %>%
  dplyr::summarise(n = n())
tmp2 <- tmp %>%
  select(ward) %>%
  dplyr::summarise(nHerdInWard = n())
summary(tmp2$nHerdInWard) # 1 to 42

# Number of records per cow per herd
(tmp <- table(data1$herd, data1$cow))
(tmp2 <- table(tmp[tmp > 0]))
hist(tmp2) # 3 to 37

# number of records per herd
(tmp <- table(data1$herd))
(tmp2 <- table(tmp[tmp > 0]))
hist(tmp2) # 3 to 71

nrecords_byherd<- data1 %>%
  group_by(herd) %>%
  dplyr::summarise(recordherd=n())
summary(nrecords_byherd$recordherd)

#number of records per ward
nrecords_byward <- data1 %>%
  group_by(ward) %>%
  dplyr::summarise(nrecord=n())
summary(nrecords_byward$nrecord)

# Number of records per cow per ward
ncowward<- data1 %>%
  group_by(ward,cow)%>%
  dplyr::summarise(ncow_byward= n())
summary(ncowward$ncow_byward)
# number of cows per ward
ncowward_2<- ncowward %>%
  select(ward) %>%
  dplyr::summarise(ncowward=n())
summary(ncowward_2$ncowward)

# Number of parities
summary(data1$lac)

# Number of cows per parity
(tmp <- data1 %>%
    group_by(lac, cow) %>%
    dplyr::summarise(n = n()) %>%
    select(lac) %>%
    dplyr::summarise(n = n()))
#  lac     no. cows
# 1 1      1031
# 2 2      1217
# 3 3       735
# 4 4       279
# 5 5        69
# 6 6        16
# 7 7         6
# 8 8         2
# 9 9         1

(tmp <- data1 %>%
    group_by(lacest, cow) %>%
    dplyr::summarise(n = n()) %>%
    select(lacest) %>%
    dplyr::summarise(n = n()))
# 1 1       1031
# 2 2       1465

# Summary age
summary(data1$age)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 18.05   38.79   47.83   50.28   58.88  162.10

#number of levels of calving seasons
nlevels(data1$ksea) # 6

#number of levels of calving year season
nlevels(data1$cyrsn) # 29

#levels of test-year-month
nlevels(data1$tyrmn) # 61

# number of levels of breed proportion (exotic)
nlevels(data1$dgrp) # 4

# ---- QC coordinates ----------------------------------------------------------

# 184 missing wards - it seems from the wrong coordinates
sel <- is.na(data1$ward_code)
table(factor(data1[sel, "herd"]))
#   34   99  108  112  113  278  284  393  645  763  872  928  952 1115
#    6   12   10   15    9    7   14   41   10   13   10    9   15   13

summary(data1$lat)
#.    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -9.45992 -8.30842 -5.03902 -5.52181 -3.35278  0.09999
summary(data1$lat[!sel]) # animals with known ward
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -9.460  -8.305  -5.038  -5.512  -3.353  -3.070
summary(data1$lat[sel]) # animals with unknown ward
#.    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -9.41693 -9.30934 -9.30494 -6.54163  0.09999  0.09999
# they go to ~0, while with known ward they go to ~-3

summary(data1$long)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0001 35.2970 36.7329 36.2039 37.5846 39.1247
summary(data1$long[!sel]) # animals with known ward
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 32.86   35.32   36.73   36.55   37.59   39.12
summary(data1$long[sel]) # animals with unknown ward
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0001  0.0001  0.0001  0.2935  0.9999  0.9999

tmpHerds <- unique(data1[sel, "herd"])
selTmpHerds <- data1$herd %in% tmpHerds
sum(selTmpHerds) # 231 records
tmp <- table(factor(data1[selTmpHerds, "herd"]),
             paste0(data1[selTmpHerds, "lat"], "/", data1[selTmpHerds, "long"]))
selTmp <- rowSums(tmp > 0) > 1
tmp[selTmp, ]
# for herd 34 we have two locations
# -9.41693/1e-04 & -9.41693/34.7935
data1[data1$herd == "34", "long"] <- 34.7935
# for herd 108 we have two locations
# -9.30934/1e-04 & -9.30934/34.7776
data1[data1$herd == "108", "long"] <- 34.7776
# for herd 645 we have two locations
# -9.30934/1e-04 & -9.30934/34.7818
data1[data1$herd == "645", "long"] <- 34.7818
# for herd 763 we have two locations
# -9.41181/1e-04 & -9.41181/34.8018
data1[data1$herd == "763", "long"] <- 34.8018
# for herd 872 we have two locations
# -9.30876/1e-04 & -9.30876/34.7789
data1[data1$herd == "872", "long"] <- 34.7789
# for herd 952 we have two locations
# -9.30494/1e-04 & -9.30494/34.7771
data1[data1$herd == "952", "long"] <- 34.7771

# 184 missing wards - it seems from the wrong coordinates
sel <- is.na(data1$ward_code) & !(data1$herd %in% c("34", "108", "645", "763", "872", "952"))
table(factor(data1[sel, "herd"]))
# 99  112  113  278  284  393  928 1115
# 12   15    9    7   14   41    9   13
sum(sel)
# 120 still to be fixed - in 8 herds

summary(data1$lat)
#.    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -9.45992 -8.30842 -5.03902 -5.52181 -3.35278  0.09999
summary(data1$lat[!sel]) # animals with known ward
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -9.460  -8.307  -5.039  -5.525  -3.353  -3.070
summary(data1$lat[sel]) # animals with unknown ward
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -9.41240 -9.30673 -9.07640 -5.04964  0.09999  0.09999
# they go to ~0, while with known ward they go to ~-3

summary(data1$long)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0001 35.2970 36.7329 36.2039 37.5846 39.1247
summary(data1$long[!sel]) # animals with known ward
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 32.86   35.32   36.73   36.55   37.59   39.12
summary(data1$long[sel]) # animals with unknown ward
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0001  0.0001  0.0001  0.2935  0.9999  0.9999
# again all around zero!

# We don't know what to do with them so we will remove them, BUT these farms
# should have coordinates fixed in the future
(herdsToRemove <- unique(data1[sel, "herd"]))
herdsToRemove <- c(99, 112, 113, 278, 284, 393, 928, 1115)
sel <- data1$herd %in% herdsToRemove
sum(sel) # 120
data1 <- data1[!sel,]
nrow(data1) # 19418
length(unique(data1$cow)) # 1899

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
table(data1$lacgr)

length(unique(data1$herd))
#Number of rows
nrow(data1) # 19494 records
#number of variables
ncol(data1)#  16
#Number of cows
length(unique(data1$cow)) # 1906
# Milk Yield
summary(data1$milk) #min=1, median=8.00 max= 40
mean(data1$milk) #8.33
var(data1$milk) # 18.66
sd(data1$milk) # 4.32
hist(na.omit(data1$milk))

# ---- Ward --------------------------------------------------------------------

#Number of levels of wards
length(unique(data1$ward)) # 156 wards
table(data1$ward)
table(table(data1$ward)) # number of records per ward ranged from 3 to 979
plot(data1$ward,data1$milk)

# Load library dplyr
library(dplyr)

#calculate mean and sd of milk yield  by ward
data1_ward_mean_std <- data1 %>%
  group_by(ward) %>%
  summarise_at(vars(milk), list(mean=mean, sd=sd)) %>%
  as.data.frame()

#view results
data1_ward_mean_std
#plot mean and standard deviation of milk yield by ward
ggplot(data1_ward_mean_std , aes(x=ward, y=mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3) +
  geom_point(size=2)

# ---- Calving year-season -----------------------------------------------------

#Number of levels of CYS
length(unique(data1$cyrsn)) # 29 levels
table(data1$cyrsn)
table(table(data1$cyrsn)) # Ranges from 1 to 2663

# Plot mean and sd of milk yield  per calving year season level

data1_cyrsn_mean_std <- data1 %>%
  group_by(cyrsn) %>%
  summarise_at(vars(milk), list(mean=mean, sd=sd)) %>%
  as.data.frame()

#view results
data1_cyrsn_mean_std
#plot mean and standard deviation of milk yield by ward
ggplot(data1_cyrsn_mean_std , aes(x=cyrsn, y=mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3) +
  geom_point(size=1)

# ---- Lactation group ---------------------------------------------------------

length(unique(data1$lacgr)) # 3 levels
table(data1$lacgr)
table(table(data1$lacgr))

# Plot mean and sd of milk yield  per lactation group

data1_lacgr_mean_std <- data1 %>%
  group_by(lacgr) %>%
  summarise_at(vars(milk), list(mean=mean, sd=sd)) %>%
  as.data.frame()

#view results
data1_lacgr_mean_std
#plot mean and standard deviation of milk yield by ward
ggplot(data1_lacgr_mean_std , aes(x=lacgr, y=mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3) +
  geom_point(size=2)

# ---- Herd --------------------------------------------------------------------

data1$herd = as.factor(data1$herd)
length(unique(data1$herd)) # 1394
table(data1$herd)
table(table(data1$herd))
# Plot mean and sd of milk yield  per herd

data1_herd_mean_std <- data1 %>%
  group_by(herd) %>%
  summarise_at(vars(milk), list(mean=mean, sd=sd)) %>%
  as.data.frame()

#view results
data1_herd_mean_std
#plot mean and standard deviation of milk yield by ward
ggplot(data1_herd_mean_std , aes(x=herd, y=mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3) +
  geom_point(size=2)

# ---- Explore coordinates -----------------------------------------------------

# 1. How many records we have?

length(data1$cow)
# 19538 records
# 19494 records

# 2. How many unique animals with records we have?
length(unique(data1$cow))
# 1911
# 1906

# 3. How many herds we have?
length(unique(data1$herd))
# 1394 herds
# 4. How many records per herd we have?
table(table(data1$herd))
# Number of record per herd ranged from 3 to 71

nrecords_byherd<- data1 %>%
  group_by(herd) %>%
  summarise(recordherd=n())
summary(nrecords_byherd$recordherd)

summary(nrecords_byherd$recordherd)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 3.00    7.00   11.00   13.98   18.00   71.00

# 5. How many animals per herd we have?
# Number of records per cow per herd
ncowherd<- data1 %>%
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
summary(data1$lac)
summary(data1$age) # age in months
#Min.   1st Qu.  Median Mean    3rd Qu.    Max.
#18.05   38.79   47.86   50.29   58.88     162.10
# Histogram of age
hist(na.omit(data1$age))
head(data1, n=1)

# ----- Descriptive statistics of coordinates ----------------------------------

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

# ---- Mapping  new data (data5) without wrong GPS in R ------------------------

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

# ---- Ward spatial ------------------------------------------------------------

# Ward Neigbouring relationship

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


