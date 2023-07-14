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
baseDir <- "/Users/ggorjanc/Storages/GitBox/HighlanderLab/ihouaga_adgg_spatial"

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
    "orthopolynom" # for legendre.polynomials()
  )
  install.packages(pkgs = requiredPackages)
}
library(tidyverse)
library(rgdal)
library(sf)
library(raster)
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

# ---- Import shapefile geo data -----------------------------------------------

map <- readOGR(dsn = "data/shapefiles/wards_2011/TZwards.shp")
plot(map)
# coordinates are in columns long and lat
coordinates(data1) <- ~long + lat

# Define projection of coordinates
projection(map) <- CRS("+proj=longlat +datum=WGS84")
projection(data1) <- CRS("+proj=longlat +datum=WGS84")

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

# ---- Factor and numeric variables --------------------------------------------

str(data1)

# Factors
data1 <- data1 %>%
  mutate_at(.vars = c("cow", "ward", "herd", "cyrsn", "tyrmn", "dgrp", "lacest",
                      "lac", "htd", "hcyr", "htyr", "pym", "ksea"),
            .funs = as.factor)
str(data1)

# Numeric
data1 <- data1 %>%
  mutate_at(.vars = c("dim", "bprop", "age", "milk", "hgirth", "bodywt", "bcs",
                      "long", "lat"),
            .funs = as.numeric)
str(data1)

# ---- Check for duplicates and NA ---------------------------------------------

nrow(data1) # 19538

# Remove potential duplicated rows
data1 <- dplyr::distinct(data1)

nrow(data1) # 19538, no duplicates it seems

# Check for missing data (NA)
sum(is.na(data1)) # 0 nothing missing as a NA

# NOTE: we have seen some 0 for hgirth, which looks like a missing value!

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

# ---- Milk stats --------------------------------------------------------------

# Milk test-day milk yield across 9 lactations
min(data1$milk)
max(data1$milk)
mean(data1$milk, na.rm = TRUE)
sd(data1$milk, na.rm = TRUE)
hist(na.omit(data1$milk))

# Mean and Variance - grouped by parity (lactation)
SMilk <- data1 %>%
  group_by(lac) %>%
  dplyr::summarise(
    MeanMilk = mean(milk, na.rm = TRUE),
    MedianMilk = median(milk, na.rm = TRUE),
    VarMilk = var(milk, na.rm = TRUE),
    SdMilk =  sqrt(VarMilk),
    nMilk = n()
  )
SMilk

SMilk %>%
  ggplot(aes(y = MeanMilk, x = lac,  colour = lac)) +
  geom_point(aes(label = lac)) +
  xlab("Parity") +
  ylab("E(Milk)") +
  theme_bw()

SMilk %>%
  ggplot(aes(y = SdMilk, x = lac,  colour = lac)) +
  geom_point(aes(label = lac)) +
  xlab("Parity") +
  ylab("Sd(Milk)") +
  theme_bw()

# We can see the milk production increases up to parity xxxx 3-4 and then
# starts to decrease. This is related to the combination of effect of
# parity and number of observations???? It's fine - this is exepected from such
# data - senescence and selection.

data1 %>%
  mutate_at("lac", factor) %>%
  mutate(MYc = milk-mean(milk)) %>%
  ggplot(aes(y = MYc, x = lac)) +
  geom_boxplot() +
  xlab("Levels of Parity") +
  ylab("Test-day Milk Yield") +
  theme_bw()

# ---- Export cleaned data for blupf90 & INLA as in Mrode et al. (2021) --------

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

# Exporting data for blupf90
selectColumns <- c(
  "cow", # 1.
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
  "leg2" # 15
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
