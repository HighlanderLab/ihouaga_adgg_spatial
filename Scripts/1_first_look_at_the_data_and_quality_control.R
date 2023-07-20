# ---- Header ------------------------------------------------------------------

# First look at the data, quality control, preparation, and export for later
#   analyses with other tools
# Trait: Test day milk yield
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
baseDir <- "/Users/ggorjanc/Storages/GitBox/HighlanderLab/ihouaga_adgg_spatial_upstream/"

# Change working directory
setwd(dir = baseDir)
dir()

# ---- Installing and load packages --------------------------------------------

if (FALSE) {
  requiredPackages <- c(
    "tidyverse", # for data manipulation
    "rgdal", # for shapefile work
    "sf", # for shapefile work
    "raster", # for shapefile work
    "spdep", # for shapefile work - to identify neighbouring regions
    "orthopolynom" # for legendre.polynomials()
  )
  install.packages(pkgs = requiredPackages)
  install.packages(pkgs = "INLA",
                   repos=c(getOption("repos"),
                           INLA = "https://inla.r-inla-download.org/R/stable"),
                   dep = TRUE)
  inla.upgrade(testing = TRUE)
}
library(tidyverse)
library(rgdal)
library(sf)
library(raster)
library(spdep)
library(INLA)
library(orthopolynom)

(.packages()) # Check loaded packages

# ---- Import pheno data -------------------------------------------------------

# Clear the environment
rm(list = ls())

data1 <- read.table(file = "data/original_data/milk4d-gps-pre.txt", header = FALSE)
nrow(data1) # 19538
ncol(data1) # 22
head(data1)
# V1 V2 V3  V4   V5 V6 V7 V8  V9   V10 V11 V12  V13 V14     V15   V16 V17 V18     V19 V20    V21      V22
# 1 296  4  1 394 4001  1  2  4 143 14524 139 140 4001   4 0.92618 82.45  10 160 291.190   3 36.745 -3.37952
# 2 296  4  1 394 3910  1  2  4  55 14435 139 139 3910   4 0.92618 82.45  10 160 291.190   3 36.745 -3.37952
# 3 296  4  1 394 3911  1  2  4  90 14470 139 139 3911   4 0.92618 82.45  10 159 286.913   4 36.745 -3.37952
# 4 296  4  1 394 4002  1  2  4 183 14564 139 140 4002   4 0.92618 82.45  10 160 291.190   3 36.745 -3.37952
# 5 296  4  1 373 3805  1  2  2 369 13890 137 138 3805   3 0.92618 54.50   8 156 274.082   3 36.745 -3.37952
# 6 296  4  1 373 3802  1  2  2 282 13799 137 138 3802   3 0.92618 54.50   8 158 282.636   3 36.745 -3.37952

summary(data1$V22)
table(data1$V20)
length(table(data1$V15))

colnames(data1) <- c("cow", # 1 cow ID in numeric form 1:1911
                     "ward", # 2 ward ID in numeric form 1:157
                     "herd", # 3 herd ID in numeric form 1:1400
                     "cyrsn", # 4 calving year season
                     "tyrmn", # 5 test year month
                     "dgrp", # 6 breed proportion class 1:4
                     "lacest", # 7 lactation 1 and 2+
                     "lac", # 8 lactation not combined 1:9
                     "dim", # 9 days in milk 4:500
                     "htd", # 10 herd-test-day (16352 htd)
                     "hcyr", # 11 herd calving-year (2729 hcyr)
                     "htyr", # 12 herd-test-year (3835 htyr)
                     "pym", # 13 milk-year-month (61 pym - strange coding!?)
                     "ksea", # 14 calving season 1:6
                     "bprop", # 15 proportion of DNA from exotic breeds [0, 1]
                     "age", # 16 age at test in months [18, 162]
                     "milk", # 17 test-day milk yield [1, 40]
                     "hgirth", # 18 hearth-girth [0, 9826.0] - errors or missing data coded as 0!?
                     "bodywt", # 19 body weight [-29, 41632.67] - errors too!?
                     "bcs", # 20 body condition score [1, 5]
                     "long", # 21 longitude (east-west - horizontal dimension) [0, 39]
                     "lat") # 22 latitude (north-south - vertical dimension) [-9, 0]

# ---- Fix coordinates -------------------------------------------------------

# Based on preliminary exploratory analysis we found these mistakes
# (see 2_exploratory_analysis.R). IMPORTANT, to be able to find these mistakes,
# you have to comment this part of the code, export the RAW (uncleaned) data
# and use 2_exploratory_analysis.R.

# if (FALSE) {
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

# For these herds we only had one coordinates so we could not fix them and we
# have hence, sadly, removed them from the analysis
herdsToRemove <- c(99, 112, 113, 278, 284, 393, 928, 1115)
sel <- data1$herd %in% herdsToRemove
data1 <- data1[!sel, ]
# }

# ---- Import shapefile geo data -----------------------------------------------

map <- readOGR(dsn = "data/shapefiles/wards_2011/TZwards.shp")
# plot(map) # this takes quite a bit of time so we comment it out by default
# coordinates are in columns long and lat
data1SpatialPointsDataFrame <- data1
coordinates(data1SpatialPointsDataFrame) <- ~long + lat

# Define projection of coordinates
projection(map) <- CRS(projargs = "+proj=longlat +datum=WGS84")
projection(data1SpatialPointsDataFrame) <- CRS(projargs = "+proj=longlat +datum=WGS84")

# same results if that's how projection is defined
#proj4string(data) <- CRS("+proj=longlat +a=6378249.145 +rf=293.465 +no_defs +type=crs")
#projection(mapwa) <- CRS("+proj=longlat +a=6378249.145 +rf=293.465 +no_defs +type=crs")

# Overlay data points with the map to allocate data points to map features
tmp <- over(x = geometry(data1SpatialPointsDataFrame), y = map, returnList = FALSE)
head(tmp)
data1$ward_code <- tmp$Ward_Code
# Note that we could bring in also region and district code, but to do Besag type
# modelling, we would also need to get shapefiles of these spatial organisation
# units. There are 6 regions, so that will not be very fine-grained modelling.
head(data1)

# ---- Create Besag neighbourhood matrix/file ----------------------------------

# Read in the shapefile
map <- st_read(dsn = "data/shapefiles/wards_2011/TZwards.shp")
class(map) # class sf and data.frame
head(map)
if (FALSE) {
  plot(map) # this takes quite a bit of time so we comment it out by default
  ggplot(map) +
    geom_sf() +
    geom_sf_text(data = map, aes(label = Ward_Code), size = 3, colour = "black") +
    geom_point(data = data1, aes(x = long, y = lat))
  plot(data1$long, data1$lat) # looks about alright!
}

# Construct graph of neighbouring regions for later INLA modelling

# Following the advice from https://github.com/r-spatial/sf/issues/1762
sf_use_s2(use_s2 = FALSE)
nb.map <- poly2nb(map) # Construct neighbours list from polygon list
nb2INLA(file = "data/cleaned_data/ward_neighbours.txt", nb = nb.map)
nb.map <- inla.read.graph(filename = "data/cleaned_data/ward_neighbours.txt")

nb.matrix <- -inla.graph2matrix(nb.map) # Convert graph to matrix
dim(nb.matrix) # 3644 by 3644
nb.matrix[1:5, 1:5]
diag(nb.matrix) <- 0
diag(nb.matrix) <- -rowSums(nb.matrix)
nb.matrix[1:5, 1:5]
n <- dim(nb.matrix)[1]
nb.matrixScaled <- inla.scale.model(nb.matrix,
                                    constr = list(A = matrix(1, 1, n), e = 0))
# uncomment the below to see what scaling does, but best to keep it commented out
# because it takes a lot of time to run, ginv() part is slow:(
# print(diag(MASS::ginv(as.matrix(nb.matrixScaled))))

# Saving matrix later for INLA
save(nb.matrix, nb.matrixScaled,
     file = "data/cleaned_data/ward_neighbours_precision_matrix.RData")

# Saving triplets form for blupf90
nb.matrixTril <- summary(tril(nb.matrix))
nb.matrixScaledTril <- summary(tril(nb.matrixScaled))
write.table(x = nb.matrixTril,
            file = "data/cleaned_data/ward_neighbours_precision_matrix_tril.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
write.table(x = nb.matrixScaledTril,
            file = "data/cleaned_data/ward_neighbours_precision_matrixScaled_tril.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")

# ---- Regions based on clustering of the coordinates --------------------------

# We did the following based on looking at this plot
plot(data1$long, data1$lat) # looks about alright!

data1$region <- 1
data1$region[data1$lat > -4] <- 2
data1$region[data1$lat < -7 & data1$long > 34] <- 3
data1$region[data1$lat < -7 & data1$long < 34] <- 4
table(data1$region)
#    1    2    3    4
# 4547 8525 3674 2672

# ---- Factor and numeric variables --------------------------------------------

str(data1)

# Factors
data1 <- data1 %>%
  mutate_at(.vars = c("cow", "ward", "herd", "cyrsn", "tyrmn", "dgrp", "lacest",
                      "lac", "htd", "hcyr", "htyr", "pym", "ksea", "ward_code",
                      "region"),
            .funs = as.factor)
str(data1)

# Numeric
data1 <- data1 %>%
  mutate_at(.vars = c("dim", "bprop", "age", "milk", "hgirth", "bodywt", "bcs",
                      "long", "lat"),
            .funs = as.numeric)
str(data1)

# ---- EXPORT raw data for exploratory analysis --------------------------------

# Uncomment this only if you want to get the raw data, but then also make sure
# that any data edits above have been done too!
if (FALSE) {
  selectColumns <- c(
    "cow",
    "milk",
    "ward",
    "herd",
    "cyrsn",
    "tyrmn",
    "dgrp",
    "lac",
    "age",
    "long",
    "lat",
    "ward_code",
    "region"
  )
  write.csv(x = data1[, selectColumns],
            file = "data/cleaned_data/milk_yield_pheno_raw.csv",
            row.names = FALSE)
}

# ---- Remove duplicates -------------------------------------------------------

nrow(data1) # 19418

# Remove potential duplicated rows
data1 <- dplyr::distinct(data1)
nrow(data1) # 19418, no duplicates!

# Check for missing data (NA)
sum(is.na(data1)) # 0

# NOTE: we have seen some 0 for hgirth, which looks like a missing value!

# ---- EXPORT cleaned data for blupf90 & INLA as in Mrode et al. (2021) --------

dir.create(path = "data/cleaned_data/")

# Create column for grouped parity (1,2 and 3+)
data1$lac <- as.numeric(data1$lac)
data1 <- data1 %>%
  mutate(lacgr = lac) %>%
  mutate(lacgr = ifelse(lacgr >= 3, 3, lacgr))
table(data1$lac, data1$lacgr)
data1$lacgr <- as.factor(data1$lacgr)

# Create columns for intercept and legendre polynomial
m1 <- legendre.polynomials(n = 2, normalized = TRUE)
# u + ((v - u) * ((x - min(x))/(max(x) - min(x)))
# if u=-1 and v=1 than -1 + 2*((x - min(x))/(max(x) - min(x)))
m1 <- polynomial.values(polynomials = m1,
                        x = scaleX(data1$dim, u=-1, v=1))
dataLegendre <- as.matrix(as.data.frame(m1))
colnames(dataLegendre) <- c("leg0", "leg1", "leg2")
data1 <- cbind(data1, dataLegendre)
head(data1)

# Exporting data for blupf90, INLA, etc.
selectColumns <- c(
  "cow", # 1
  "milk", # 2
  "ward", # 3
  "herd", # 4
  "cyrsn", # 5
  "tyrmn", # 6
  "dgrp", # 7
  "lac", # 8
  "lacgr", # 9
  "age", # 10
  "long", # 11
  "lat", # 12
  "leg0", # 13
  "leg1", # 14
  "leg2", # 15
  "ward_code", # 16
  "region" # 17
)
write.table(x = data1[, selectColumns],
            file = "data/cleaned_data/milk_yield_pheno_cleaned_for_blupf90_no_header.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")

write.csv(x = data1[, selectColumns],
          file = "data/cleaned_data/milk_yield_pheno_cleaned.csv",
          row.names = FALSE)

# ---- Toydata -----------------------------------------------------------------

# write.csv(dataHST, "/Users/ihouaga/Documents/Project1/dataHSTRaw.csv") ## save .csv copy of raw data
# My Toy data to test scripts
# cow<- c(1,1,1,1,1,2,2,2,3,3,4,4,5,5,5)
# herd<- c(1,1,1,1,1,3,3,3,2,2,3,3, 2,2,2)
# ward <- c(1,1,1,1,1,2,2,2,2,2,2,2,2,2,2)
# lac <- c(1,2,3,4,5,1,2,3,1,2,1,2,1,2,3)
# milk<- c(5,6,7,8,9,4,5,6,3,4,5,6,5,6,7)
# toy<- data.frame(cow,herd,ward,lac,milk)
# toy
# toy$milk <- as.numeric(toy$milk)
# toy$cow <- as.factor(toy$cow)
# toy$herd <- as.factor(toy$herd)
# toy$ward <- as.factor(toy$ward)
# toy$lac <- as.factor(toy$lac)
