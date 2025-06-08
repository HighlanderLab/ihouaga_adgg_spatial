# ---- Header ------------------------------------------------------------------
# Spatial modelling
# Trait: Test day milk yield
# Author: Isidore Houaga, Ivan Pocrnic and Gregor Gorjanc.
# Version 1.1.0
# Date: 2025-06-08
# ---- Setup--------------------------------------------------------------------

# Working directory
baseDir <- "/Users/ihouaga2/ihouaga_adgg_spatial"
# Change working directory
setwd(dir = baseDir)
getwd()

# ---- Installing and loading packages--------------------------------------------

if (FALSE) {
  requiredPackages <- c(
    "tidyverse", # for data manipulation
    "fmesher",
    "inlabru",
    "verification", # for Continuous Ranked Probability Score
    "irlba" # for fast PCA
  )
  install.packages(pkgs = requiredPackages)
  install.packages(pkgs = "INLA",
                   repos=c(getOption("repos"),
                           INLA = "https://inla.r-inla-download.org/R/stable"),
                   dep = TRUE)
  inla.upgrade(testing = TRUE)
}
library(tidyverse)
library(INLA)
library(inlabru)
library(sf)
library(sp)
library(fmesher)
library(gridExtra) # Visualize random field grid
library(verification) # Visualize random field grid
library(lattice)
library(raster) # GetData to plot country map
library(rgeos)
library(geodata)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(viridis) # Plot mean and sd spatial effect
library(fields)# Plot mean and sd spatial effect
library(ggpubr) # Plot mean and sd spatial effect
library(irlba)
library(Hmisc) # Correlation matrix with P-values
#library("MASS") # Pca
#library("factoextra") # Pca
library(easyGgplot2)
library(psych) # scatter-plot matrix

(.packages()) # Check loaded packages

# ---- Import data -------------------------------------------------------------
# Clear the environment
rm(list = ls())
#-------------------Accuracy of Prediction: Cross validation--------------------

#----------Cross-validation_Breed_Proportion_Paper_Predicted Pheno Vs Obs_Pheno------------------
#-----------------Cross-validation_Breed_Proportion fitGP-----------------------------------------
# Masking phenotypes of cows in dgrp (Breed proportion class)
#Creating data with missing penotypes
#dgrp1 masked
data2_1_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "1", NA, milkZ))
#dgrp2 masked
data2_2_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "2", NA, milkZ))
#dgrp3 masked
data2_3_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "3", NA, milkZ))
#dgrp4 masked
data2_4_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "4", NA, milkZ))

#Creating data for each dgrp
data2_1 <- subset(data2, dgrp=="1") #731 cows
data2_2 <- subset(data2, dgrp=="2") #770 cows
data2_3<- subset(data2, dgrp=="3") #309
data2_4<- subset(data2, dgrp=="4") #84
length(unique(data2_1$cow))

# Fitting models on data having missing data
fitGP_1_NA <- bru(modelGP,
                 like(family = "Gaussian",
                      modelGPFormula,
                      data = data2_1_NA))
save(fitGP_1_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGP/fitGP_1_NA.RData")


fitGP_2_NA <- bru(modelGP,
                 like(family = "Gaussian",
                      modelGPFormula,
                      data = data2_2_NA))
save(fitGP_2_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGP/fitGP_2_NA.RData")

fitGP_3_NA <- bru(modelGP,
                 like(family = "Gaussian",
                      modelGPFormula,
                      data = data2_3_NA))
save(fitGP_3_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGP/fitGP_3_NA.RData")

fitGP_4_NA <- bru(modelGP,
                 like(family = "Gaussian",
                      modelGPFormula,
                      data = data2_4_NA))
# Predicting with-held phenotypes
predictGP_1 <- predict(fitGP_1_NA, newdata=data2_1,
                      formula= ~ fixed_effects + perm  + animal)

predictGP_2 <- predict(fitGP_2_NA, newdata=data2_2,
                      formula= ~ fixed_effects + perm  + animal)

predictGP_3 <- predict(fitGP_3_NA, newdata=data2_3,
                      formula= ~ fixed_effects + perm  + animal)

predictGP_4 <- predict(fitGP_4_NA, newdata=data2_4,
                      formula= ~ fixed_effects + perm  + animal)

# Accuracy of fitGP (Correlating predicted phenotypes with observed (true) phenotypes
Accu1 <- round(cor(predictGP_1$milkZ,predictGP_1$mean),2)
Accu1 #0.29

Accu2 <- round(cor(predictGP_2$milkZ,predictGP_2$mean),2)
Accu2 #0.35

Accu3 <- round(cor(predictGP_3$milkZ,predictGP_3$mean),2)
Accu3 # 0.39

Accu4 <- round(cor(predictGP_4$milkZ,predictGP_4$mean),2)
Accu4 # 0.40
# Cross-validation accuracy of fitGP
Accu_GP <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GP # 0.36

#-----------------Cross-validation fitGPH-----------------------------------------
# Fitting models on data having missing data
fitGPH_1_NA <- bru(modelGPH,
                  like(family = "Gaussian",
                       modelGPHFormula,
                       data = data2_1_NA))
save(fitGPH_1_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGPH/fitGPH_1_NA.RData")

fitGPH_2_NA <- bru(modelGPH,
                  like(family = "Gaussian",
                       modelGPHFormula,
                       data = data2_2_NA))
save(fitGPH_2_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGPH/fitGPH_2_NA.RData")


fitGPH_3_NA <- bru(modelGPH,
                  like(family = "Gaussian",
                       modelGPHFormula,
                       data = data2_3_NA))
save(fitGPH_3_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGPH/fitGPH_3_NA.RData")


fitGPH_4_NA <- bru(modelGPH,
                  like(family = "Gaussian",
                       modelGPHFormula,
                       data = data2_4_NA))
save(fitGPH_4_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGPH/fitGPH_4_NA.RData")


# Predicting with-held phenotypes
predictGPH_1 <- predict(fitGPH_1_NA, newdata=data2_1,
                       formula= ~ fixed_effects + perm + animal + herd)

predictGPH_2 <- predict(fitGPH_2_NA, newdata=data2_2,
                       formula= ~ fixed_effects + perm + animal + herd)

predictGPH_3 <- predict(fitGPH_3_NA, newdata=data2_3,
                       formula= ~ fixed_effects + perm +  animal + herd)

predictGPH_4 <- predict(fitGPH_4_NA, newdata=data2_4,
                       formula= ~ fixed_effects + perm + animal + herd)


# Accuracy of fitGPH

Accu1 <- round(cor(predictGPH_1$milkZ,predictGPH_1$mean),2)
Accu1 #0.36

Accu2 <- round(cor(predictGPH_2$milkZ,predictGPH_2$mean),2)
Accu2 # 0.45

Accu3 <- round(cor(predictGPH_3$milkZ,predictGPH_3$mean),2)
Accu3 # 0.46

Accu4 <- round(cor(predictGPH_4$milkZ,predictGPH_4$mean),2)
Accu4 # #0.38
# Cross-validation accuracy of fitGPH
Accu_GPH <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GPH #0.41


#-----------------Cross-validation_breed_proportion fitGPS-----------------------------------------
# Fitting models on data having missing data
fitGPS_1_NA <- bru(modelGPS,
                  like(family = "Gaussian",
                       modelGPSFormula,
                       data = data2_1_NA))
save(fitGPS_1_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGPS/fitGPS_1_NA.RData")

fitGPS_2_NA <- bru(modelGPS,
                  like(family = "Gaussian",
                       modelGPSFormula,
                       data = data2_2_NA))
save(fitGPS_2_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGPS/fitGPS_2_NA.RData")

fitGPS_3_NA <- bru(modelGPS,
                  like(family = "Gaussian",
                       modelGPSFormula,
                       data = data2_3_NA))
save(fitGPS_3_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGPS/fitGPS_3_NA.RData")

fitGPS_4_NA <- bru(modelGPS,
                  like(family = "Gaussian",
                       modelGPSFormula,
                       data = data2_4_NA))
save(fitGPS_4_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGPS/fitGPS_4_NA.RData")

# Predicting with-held phenotypes
predictGPS_1 <- predict(fitGPS_1_NA, newdata=data2_1,
                       formula= ~ fixed_effects + perm + animal + field)

predictGPS_2 <- predict(fitGPS_2_NA, newdata=data2_2,
                       formula= ~ fixed_effects + perm + animal + field)

predictGPS_3 <- predict(fitGPS_3_NA, newdata=data2_3,
                       formula= ~ fixed_effects + perm + animal + field)

predictGPS_4 <- predict(fitGPS_4_NA, newdata=data2_4,
                       formula= ~ fixed_effects + perm + animal + field)

# Accuracy of fitGPS

Accu1 <- round(cor(predictGPS_1$milkZ,predictGPS_1$mean),2)
Accu1 # 0.56

Accu2 <- round(cor(predictGPS_2$milkZ,predictGPS_2$mean),2)
Accu2 # 0.65

Accu3 <- round(cor(predictGPS_3$milkZ,predictGPS_3$mean),2)
Accu3 # 0.64

Accu4 <- round(cor(predictGPS_4$milkZ,predictGPS_4$mean),2)
Accu4 # 0.57

# Cross-validation accuracy of fitGPS
Accu_GPS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GPS #  0.60

#-----------------Cross-validation_breed_proportion fitGPHS-----------------------------------------
# Fitting models on data having missing data
fitGPHS_1_NA <- bru(modelGPHS,
                   like(family = "Gaussian",
                        modelGPHSFormula,
                        data = data2_1_NA))
save(fitGPHS_1_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGPHS/fitGPHS_1_NA.RData")

fitGPHS_2_NA <- bru(modelGPHS,
                   like(family = "Gaussian",
                        modelGPHSFormula,
                        data = data2_2_NA))
save(fitGPHS_2_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGPHS/fitGPHS_2_NA.RData")

fitGPHS_3_NA <- bru(modelGPHS,
                   like(family = "Gaussian",
                        modelGPHSFormula,
                        data = data2_3_NA))
save(fitGPHS_3_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGPHS/fitGPHS_3_NA.RData")

fitGPHS_4_NA <- bru(modelGPHS,
                   like(family = "Gaussian",
                        modelGPHSFormula,
                        data = data2_4_NA))
save(fitGPHS_4_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGPHS/fitGPHS_4_NA.RData")

# Predicting with-held phenotypes
predictGPHS_1 <- predict(fitGPHS_1_NA, newdata=data2_1,
                        formula= ~ fixed_effects + perm + animal + herd +field)

predictGPHS_2 <- predict(fitGPHS_2_NA, newdata=data2_2,
                        formula= ~ fixed_effects + perm + animal + herd + field)

predictGPHS_3 <- predict(fitGPHS_3_NA, newdata=data2_3,
                        formula= ~ fixed_effects + perm + animal + herd + field)

predictGPHS_4 <- predict(fitGPHS_4_NA, newdata=data2_4,
                        formula= ~ fixed_effects + perm + animal + herd + field)

# Accuracy of fitGPHS

Accu1 <- round(cor(predictGPHS_1$milkZ,predictGPHS_1$mean),2)
Accu1 #0.58

Accu2 <- round(cor(predictGPHS_2$milkZ,predictGPHS_2$mean),2)
Accu2 # 0.67

Accu3 <- round(cor(predictGPHS_3$milkZ,predictGPHS_3$mean),2)
Accu3 # 0.66

Accu4 <- round(cor(predictGPHS_4$milkZ,predictGPHS_4$mean),2)
Accu4 # 0.57
# Cross-validation accuracy of fitGPH
Accu_GPHS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GPHS # 0.62

#-----------------Cross_validation_Region_Paper-----------------------------------------
# Masking phenotypes of cows in specific region
#Creating data with missing penotypes
#region 1 masked
data2_1_NA<- data2 %>% mutate(milkZ = ifelse(region == "1", NA, milkZ))
#region 2 masked
data2_2_NA<- data2 %>% mutate(milkZ = ifelse(region == "2", NA, milkZ))
#region3 masked
data2_3_NA<- data2 %>% mutate(milkZ = ifelse(region == "3", NA, milkZ))
#region4 masked
data2_4_NA<- data2 %>% mutate(milkZ = ifelse(region == "4", NA, milkZ))
summary(data2_1_NA$milkZ)
#Creating data for each region
data2_1 <- subset(data2, region=="1") # 437 cows and 4547 record
data2_2 <- subset(data2, region=="2") # 775 cows and 8482 records
data2_3<- subset(data2, region=="3") # 414 Cows and 3674 records
data2_4<- subset(data2, region=="4") # 268 cows and 2672 records
dim(data2_4)
length(unique(data2_4$cow))

#-----------------Cross-validation_region fitGP-----------------------------------------
#----------------Models with Permanent env effect----------------------------------------
# Fitting models on data having missing data
fitGP_1_NA <- bru(modelGP,
                 like(family = "Gaussian",
                      modelGPFormula,
                      data = data2_1_NA))
save(fitGP_1_NA,file = "data/cleaned_data/cross_validation_region/fitGP/fitGP_1_NA.RData")

fitGP_2_NA <- bru(modelGP,
                 like(family = "Gaussian",
                      modelGPFormula,
                      data = data2_2_NA))
save(fitGP_2_NA,file = "data/cleaned_data/cross_validation_region/fitGP/fitGP_2_NA.RData")

fitGP_3_NA <- bru(modelGP,
                 like(family = "Gaussian",
                      modelGPFormula,
                      data = data2_3_NA))
save(fitGP_3_NA,file = "data/cleaned_data/cross_validation_region/fitGP/fitGP_3_NA.RData")

fitGP_4_NA <- bru(modelGP,
                 like(family = "Gaussian",
                      modelGPFormula,
                      data = data2_4_NA))
save(fitGP_4_NA,file = "data/cleaned_data/cross_validation_region/fitGP/fitGP_4_NA.RData")

# Predicting with-held phenotypes
predictGP_1 <- predict(fitGP_1_NA, newdata=data2_1,
                      formula= ~ fixed_effects + perm  + animal)

predictGP_2 <- predict(fitGP_2_NA, newdata=data2_2,
                      formula= ~ fixed_effects + perm + animal)

predictGP_3 <- predict(fitGP_3_NA, newdata=data2_3,
                      formula= ~ fixed_effects + perm + animal)

predictGP_4 <- predict(fitGP_4_NA, newdata=data2_4,
                      formula= ~ fixed_effects + perm + animal)

# Accuracy-region of fitGP

Accu1 <- round(cor(predictGP_1$milkZ,predictGP_1$mean),2)
Accu1 #0.08

Accu2 <- round(cor(predictGP_2$milkZ,predictGP_2$mean),2)
Accu2 #-0.01

Accu3 <- round(cor(predictGP_3$milkZ,predictGP_3$mean),2)
Accu3 #0.27

Accu4 <- round(cor(predictGP_4$milkZ,predictGP_4$mean),2)
Accu4 #0.15
# Cross-validation accuracy of fitGP
Accu_GP <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GP #0.12


#-----------------Cross-validation_region fitGPH---------------------------------
# Fitting models on data having missing data
fitGPH_1_NA <- bru(modelGPH,
                  like(family = "Gaussian",
                       modelGPHFormula,
                       data = data2_1_NA))
save(fitGPH_1_NA,file = "data/cleaned_data/cross_validation_region/fitGPH/fitGPH_1_NA.RData")
summary(data2_1_NA$dgrp)
fitGPH_2_NA <- bru(modelGPH,
                  like(family = "Gaussian",
                       modelGPHFormula,
                       data = data2_2_NA))
save(fitGPH_2_NA,file = "data/cleaned_data/cross_validation_region/fitGPH/fitGPH_2_NA.RData")


fitGPH_3_NA <- bru(modelGPH,
                  like(family = "Gaussian",
                       modelGPHFormula,
                       data = data2_3_NA))
save(fitGPH_3_NA,file = "data/cleaned_data/cross_validation_region/fitGPH/fitGPH_3_NA.RData")


fitGPH_4_NA <- bru(modelGPH,
                  like(family = "Gaussian",
                       modelGPHFormula,
                       data = data2_4_NA))
save(fitGPH_4_NA,file = "data/cleaned_data/cross_validation_region/fitGPH/fitGPH_4_NA.RData")


# Predicting with-held phenotypes
predictGPH_1 <- predict(fitGPH_1_NA, newdata=data2_1,
                       formula= ~ fixed_effects + perm + animal + herd)

predictGPH_2 <- predict(fitGPH_2_NA, newdata=data2_2,
                       formula= ~ fixed_effects + perm  + animal + herd)

predictGPH_3 <- predict(fitGPH_3_NA, newdata=data2_3,
                       formula= ~ fixed_effects + perm + animal + herd)

predictGPH_4 <- predict(fitGPH_4_NA, newdata=data2_4,
                       formula= ~ fixed_effects + perm  + animal + herd)

# Accuracy-region of fitGPH

Accu1 <- round(cor(predictGPH_1$milkZ,predictGPH_1$mean),2)
Accu1 #0.03

Accu2 <- round(cor(predictGPH_2$milkZ,predictGPH_2$mean),2)
Accu2 #0.06

Accu3 <- round(cor(predictGPH_3$milkZ,predictGPH_3$mean),2)
Accu3 #0.27

Accu4 <- round(cor(predictGPH_4$milkZ,predictGPH_4$mean),2)
Accu4 #0.19
# Cross-validation accuracy of fitGPH
Accu_GPH <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GPH #0.14

#-----------------Cross-validation_region fitGPS-----------------------------------------
# Fitting models on data having missing data
fitGPS_1_NA <- bru(modelGPS,
                  like(family = "Gaussian",
                       modelGPSFormula,
                       data = data2_1_NA))
save(fitGPS_1_NA,file = "data/cleaned_data/cross_validation_region/fitGPS/fitGPS_1_NA.RData")

fitGPS_2_NA <- bru(modelGPS,
                  like(family = "Gaussian",
                       modelGPSFormula,
                       data = data2_2_NA))
save(fitGPS_2_NA,file = "data/cleaned_data/cross_validation_region/fitGPS/fitGPS_2_NA.RData")

fitGPS_3_NA <- bru(modelGPS,
                  like(family = "Gaussian",
                       modelGPSFormula,
                       data = data2_3_NA))
save(fitGPS_3_NA,file = "data/cleaned_data/cross_validation_region/fitGPS/fitGPS_3_NA.RData")

fitGPS_4_NA <- bru(modelGPS,
                  like(family = "Gaussian",
                       modelGPSFormula,
                       data = data2_4_NA))
save(fitGPS_4_NA,file = "data/cleaned_data/cross_validation_region/fitGPS/fitGPS_4_NA.RData")

# Predicting with-held phenotypes
predictGPS_1 <- predict(fitGPS_1_NA, newdata=data2_1,
                       formula= ~ fixed_effects + perm + animal + field)

predictGPS_2 <- predict(fitGPS_2_NA, newdata=data2_2,
                       formula= ~ fixed_effects + perm + animal + field)

predictGPS_3 <- predict(fitGPS_3_NA, newdata=data2_3,
                       formula= ~ fixed_effects + perm + animal + field)

predictGPS_4 <- predict(fitGPS_4_NA, newdata=data2_4,
                       formula= ~ fixed_effects + perm + animal + field)


# Accuracy_region of fitGPS

Accu1 <- round(cor(predictGPS_1$milkZ,predictGPS_1$mean),2)
Accu1 # 0.10 

Accu2 <- round(cor(predictGPS_2$milkZ,predictGPS_2$mean),2)
Accu2 # -0.08 

Accu3 <- round(cor(predictGPS_3$milkZ,predictGPS_3$mean),2)
Accu3 # 0.27 

Accu4 <- round(cor(predictGPS_4$milkZ,predictGPS_4$mean),2)
Accu4 # 0.17 
# Cross-validation accuracy of fitGPS
Accu_GPS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GPS # 0.12 
#-----------------Cross-validation_region fitGPHS-----------------------------------------
# Fitting models on data having missing data
fitGPHS_1_NA <- bru(modelGPHS,
                   like(family = "Gaussian",
                        modelGPHSFormula,
                        data = data2_1_NA))
save(fitGPHS_1_NA,file = "data/cleaned_data/cross_validation_region/fitGPHS/fitGPHS_1_NA.RData")

fitGPHS_2_NA <- bru(modelGPHS,
                   like(family = "Gaussian",
                        modelGPHSFormula,
                        data = data2_2_NA))
save(fitGPHS_2_NA,file = "data/cleaned_data/cross_validation_region/fitGPHS/fitGPHS_2_NA.RData")

fitGPHS_3_NA <- bru(modelGPHS,
                   like(family = "Gaussian",
                        modelGPHSFormula,
                        data = data2_3_NA))save(fitGPHS_3_NA,file = "data/cleaned_data/cross_validation_region/fitGPHS/fitGPHS_3_NA.RData")

fitGPHS_4_NA <- bru(modelGPHS,
                   like(family = "Gaussian",
                        modelGPHSFormula,
                        data = data2_4_NA))
save(fitGPHS_4_NA,file = "data/cleaned_data/cross_validation_region/fitGPHS/fitGPHS_4_NA.RData")

# Predicting with-held phenotypes
predictGPHS_1 <- predict(fitGPHS_1_NA, newdata=data2_1,
                        formula= ~ fixed_effects + perm + animal + herd +field)

predictGPHS_2 <- predict(fitGPHS_2_NA, newdata=data2_2,
                        formula= ~ fixed_effects + perm + animal + herd + field)

predictGPHS_3 <- predict(fitGPHS_3_NA, newdata=data2_3,
                        formula= ~ fixed_effects + perm  + animal + herd + field)

predictGPHS_4 <- predict(fitGPHS_4_NA, newdata=data2_4,
                        formula= ~ fixed_effects + perm + animal + herd + field)

# Accuracy of fitGPHS

Accu1 <- round(cor(predictGPHS_1$milkZ,predictGPHS_1$mean),2)
Accu1 # 0.04 

Accu2 <- round(cor(predictGPHS_2$milkZ,predictGPHS_2$mean),2)
Accu2 #-0.03 

Accu3 <- round(cor(predictGPHS_3$milkZ,predictGPHS_3$mean),2)
Accu3 # 0.16 

Accu4 <- round(cor(predictGPHS_4$milkZ,predictGPHS_4$mean),2)
Accu4 # 0.13 
# Cross-validation accuracy of fitGPHS
Accu_GPHS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GPHS # 0.08 

#-----------------Foward-validation_Model with permanent environmental effect---------------------------------------------
#--Predicting phenotypes of young cows born in 2016(141 cows) and 2017(5 cows)--

#importing cow IDS and birth dates

cowbdate<- read.table(file = "data/original_data/isibdate.txt",header = FALSE)
# isibdate.txt (see file description in README.txt for data)
colnames(cowbdate) <- c("cow","birthdaymonthyear")
summary(cowbdate$birthdaymonthyear)
# extract characters from 5th index to 8th index
cowbdate['birthyear'] <-str_sub(cowbdate$birthdaymonthyear, 5, 8)  # prints "last 4 characters ie year"
cowbdate <- cowbdate[c('cow','birthyear')]
str(cowbdate)
table(cowbdate$birthyear)
# 007  008  010  011  012  013  014  015  016  017 2004 2005 2006 2007 2008 2009
#2    1    4    1   90    7  319   33   58    2    1    1    1    2    5    2
#2010 2011 2012 2013 2014 2015 2016 2017
#10   22  219   93  860   97   83    3

# Standardise year format
# Replace Values Based on Condition
cowbdate$birthyear[cowbdate$birthyear == "014"] <- "2014"
cowbdate$birthyear[cowbdate$birthyear == "007"] <- "2007"
cowbdate$birthyear[cowbdate$birthyear == "008"] <- "2008"
cowbdate$birthyear[cowbdate$birthyear == "010"] <- "2010"
cowbdate$birthyear[cowbdate$birthyear == "011"] <- "2011"
cowbdate$birthyear[cowbdate$birthyear == "012"] <- "2012"
cowbdate$birthyear[cowbdate$birthyear == "013"] <- "2013"
cowbdate$birthyear[cowbdate$birthyear == "015"] <- "2015"
cowbdate$birthyear[cowbdate$birthyear == "016"] <- "2016"
cowbdate$birthyear[cowbdate$birthyear == "017"] <- "2017"
table(cowbdate$birthyear)
#2016 2017
#141  5
length(unique(cowbdate$cow))

# Create birthyear column in data2 by merging with cowbdate

data2 <-  merge(data2, cowbdate, by= "cow", all.x = TRUE)
summary(data2$birthyear)
#2004  2005  2006  2007  2008  2009  2010  2011  2012  2013  2014  2015  2016  2017
#17    12     7    38    41     7   122   232  3044  1066 12017  1298  1444    30
data2$birthyear <- as.factor(data2$birthyear)
summary(data2$birthyear)
#data2_2016_2017 <- subset(data1_birthdate, birthyear==2016 | birthyear==2017 )


#Let's make "NA" milkz of cows born in 2016 and 2017
data2_NA<- data2%>% mutate(milkZ = ifelse(birthyear=="2016" | birthyear=="2017", NA, milkZ))
sum(is.na(data2_NA$milkZ))   #1424 Expected
length(unique(data2_NA$cow)) #1894 Expected

# Let's prepare data for young cows
data2_young <- subset(data2, birthyear=="2016" | birthyear=="2017") # 146 cows 1424 records
length(unique(data2_young$cow)) #146
#1  2  3  4
#60 56 24  6

data2_young_sub <- data2_young[, c("cow","region","dgrp")]
data2_young_sub <- distinct(data2_young_sub)
table(data2_young_sub$dgrp)
table(data2_young_sub$region)


# Fitting models on data having missing data
fitGP_NA <- bru(modelGP,
                like(family = "Gaussian",
                     modelGPFormula,
                     data = data2_NA))
save(fitGP_NA,file = "data/cleaned_data/foward_validation/fitGP/fitGP_NA.RData")

fitGPH_NA <- bru(modelGPH,
                 like(family = "Gaussian",
                      modelGPHFormula,
                      data = data2_NA))
save(fitGPH_NA,file = "data/cleaned_data/foward_validation/fitGPH/fitGPH_NA.RData")



fitGPS_NA <- bru(modelGPS,
                 like(family = "Gaussian",
                      modelGPSFormula,
                      data = data2_NA))
save(fitGPS_NA,file = "data/cleaned_data/foward_validation/fitGPS/fitGPS_NA.RData")



fitGPHS_NA <- bru(modelGPHS,
                  like(family = "Gaussian",
                       modelGPHSFormula,
                       data = data2_NA))
save(fitGPHS_NA,file = "data/cleaned_data/foward_validation/fitGPHS/fitGPHS_NA.RData")

# Predicting with-held phenotypes

predictGP <- predict(fitGP_NA, newdata=data2_young,
                     formula= ~ fixed_effects + perm + animal)


predictGPH <- predict(fitGH_NA, newdata=data2_young,
                      formula= ~ fixed_effects + perm + animal + herd)


predictGPS <- predict(fitGPS_NA, newdata=data2_young,
                      formula= ~ fixed_effects + perm + animal + field)


predictGPHS <- predict(fitGPHS_NA, newdata=data2_young,
                       formula= ~ fixed_effects + perm + animal + herd + field)

# Accuracy forward validation GP

AccuGP <- round(cor(predictGP$milkZ,predictGP$mean),2)
AccuGP 

# Accuracy forward validation GPH
AccuGPH <- round(cor(predictGPH$milkZ,predictGPH$mean),2)
AccuGPH

# Accuracy forward validation GPS
AccuGPS <- round(cor(predictGPS$milkZ,predictGPS$mean),2)
AccuGPS

# Accuracy forward validation GPHS
AccuGPHS <- round(cor(predictGPHS$milkZ,predictGPHS$mean),2)
AccuGPHS

#----Summarise forward validation accuracy by dgrp and region for GP models-------------------
#Fit GPS dgrp
predictGPS <- predictGPS[, c("cow","region","dgrp", "milkZ", "mean")]
predictGPS_dgrp1 <- subset(predictGPS, dgrp=="1")
length(unique(predictGPS_dgrp1$cow)) # 60 cows
predictGPS_dgrp2 <- subset(predictGPS, dgrp=="2")
length(unique(predictGPS_dgrp2$cow)) #  56 cows

predictGPS_dgrp3 <- subset(predictGPS, dgrp=="3")
length(unique(predictGPS_dgrp3$cow)) # 24 cows

predictGPS_dgrp4 <- subset(predictGPS, dgrp=="4")
length(unique(predictGPS_dgrp4$cow)) #  6 cows

Accu_dgrp1 <-  round(cor(predictGPS_dgrp1$milkZ,predictGPS_dgrp1$mean),2)
Accu_dgrp1 # 0.36

Accu_dgrp2 <-  round(cor(predictGPS_dgrp2$milkZ,predictGPS_dgrp2$mean),2)
Accu_dgrp2 #0.36

Accu_dgrp3 <-  round(cor(predictGPS_dgrp3$milkZ,predictGPS_dgrp3$mean),2)
Accu_dgrp3 #0.60

Accu_dgrp4 <-  round(cor(predictGPS_dgrp4$milkZ,predictGPS_dgrp4$mean),2)
Accu_dgrp4 #0.70

Accu_dgrp <- round((Accu_dgrp1 + Accu_dgrp2 + Accu_dgrp3 + Accu_dgrp4)/4,2)
Accu_dgrp #0.50

#Fit GPHS dgrp
predictGPHS <- predictGPHS[, c("cow","region","dgrp", "milkZ", "mean")]
predictGPHS_dgrp1 <- subset(predictGPHS, dgrp=="1")
length(unique(predictGPHS_dgrp1$cow)) # 60 cows

predictGPHS_dgrp2 <- subset(predictGPHS, dgrp=="2")
length(unique(predictGPHS_dgrp2$cow)) #  56 cows

predictGPHS_dgrp3 <- subset(predictGPHS, dgrp=="3")
length(unique(predictGPHS_dgrp3$cow)) # 24 cows

predictGPHS_dgrp4 <- subset(predictGPHS, dgrp=="4")
length(unique(predictGPHS_dgrp4$cow)) #  6 cows

Accu_dgrp1 <-  round(cor(predictGPHS_dgrp1$milkZ,predictGPHS_dgrp1$mean),2)
Accu_dgrp1 #0.48

Accu_dgrp2 <-  round(cor(predictGPHS_dgrp2$milkZ,predictGPHS_dgrp2$mean),2)
Accu_dgrp2 #o.7

Accu_dgrp3 <-  round(cor(predictGPHS_dgrp3$milkZ,predictGPHS_dgrp3$mean),2)
Accu_dgrp3 #0.75

Accu_dgrp4 <-  round(cor(predictGPHS_dgrp4$milkZ,predictGPHS_dgrp4$mean),2)
Accu_dgrp4 #0.76

Accu_dgrp <- round((Accu_dgrp1 + Accu_dgrp2 + Accu_dgrp3 + Accu_dgrp4)/4,2)
Accu_dgrp #0.67


#FitGPS region

predictGPS <- predictGPS[, c("cow","region","region", "milkZ", "mean")]
predictGPS_region1 <- subset(predictGPS, region=="1")
length(unique(predictGPS_region1$cow)) # 30 cows

predictGPS_region2 <- subset(predictGPS, region=="2")
length(unique(predictGPS_region2$cow)) #  60 cows

predictGPS_region3 <- subset(predictGPS, region=="3")
length(unique(predictGPS_region3$cow)) # 20 cows

predictGPS_region4 <- subset(predictGPS, region=="4")
length(unique(predictGPS_region4$cow)) #  36 cows

Accu_region1 <-  round(cor(predictGPS_region1$milkZ,predictGPS_region1$mean),2)
Accu_region1 #0.61

Accu_region2 <-  round(cor(predictGPS_region2$milkZ,predictGPS_region2$mean),2)
Accu_region2 #0.56

Accu_region3 <-  round(cor(predictGPS_region3$milkZ,predictGPS_region3$mean),2)
Accu_region3 #0.41

Accu_region4 <-  round(cor(predictGPS_region4$milkZ,predictGPS_region4$mean),2)
Accu_region4 #0.38

Accu_region <- round((Accu_region1 + Accu_region2 + Accu_region3 + Accu_region4)/4,2)
Accu_region #0.49

#Fit GPHS region

predictGPHS <- predictGPHS[, c("cow","region","region", "milkZ", "mean")]
predictGPHS_region1 <- subset(predictGPHS, region=="1")
length(unique(predictGPHS_region1$cow)) # 30 cows

predictGPHS_region2 <- subset(predictGPHS, region=="2")
length(unique(predictGPHS_region2$cow)) #  60 cows

predictGPHS_region3 <- subset(predictGPHS, region=="3")
length(unique(predictGPHS_region3$cow)) # 20 cows

predictGPHS_region4 <- subset(predictGPHS, region=="4")
length(unique(predictGPHS_region4$cow)) #  36 cows

Accu_region1 <-  round(cor(predictGPHS_region1$milkZ,predictGPHS_region1$mean),2)
Accu_region1 #0.66

Accu_region2 <-  round(cor(predictGPHS_region2$milkZ,predictGPHS_region2$mean),2)
Accu_region2 #0.63

Accu_region3 <-  round(cor(predictGPHS_region3$milkZ,predictGPHS_region3$mean),2)
Accu_region3 #0.44

Accu_region4 <-  round(cor(predictGPHS_region4$milkZ,predictGPHS_region4$mean),2)
Accu_region4 #0.47

Accu_region <- round((Accu_region1 + Accu_region2 + Accu_region3 + Accu_region4)/4,2)
Accu_region #0.55


#----------------Models without Permanent env effect (G models)-----------------

#------------- Cross-validation--------------------------------------------------
# Masking phenotypes of cows in dgrp (Breed proportion class)
#Creating data with missing penotypes
#dgrp1 masked
data2_1_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "1", NA, milkZ))
#dgrp2 masked
data2_2_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "2", NA, milkZ))
#dgrp3 masked
data2_3_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "3", NA, milkZ))
#dgrp4 masked
data2_4_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "4", NA, milkZ))

#Creating data for each dgrp
data2_1 <- subset(data2, dgrp=="1") #731 cows
data2_2 <- subset(data2, dgrp=="2") #770 cows
data2_3<- subset(data2, dgrp=="3") # 309
data2_4<- subset(data2, dgrp=="4") # 84

# Fitting models on data having missing data
fitG_1_NA <- bru(modelG,
                  like(family = "Gaussian",
                       modelGFormula,
                       data = data2_1_NA))
save(fitG_1_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitG/fitG_1_NA.RData")


fitG_2_NA <- bru(modelG,
                  like(family = "Gaussian",
                       modelGFormula,
                       data = data2_2_NA))
save(fitG_2_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitG/fitG_2_NA.RData")

fitG_3_NA <- bru(modelG,
                  like(family = "Gaussian",
                       modelGFormula,
                       data = data2_3_NA))
save(fitG_3_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitG/fitG_3_NA.RData")

fitG_4_NA <- bru(modelG,
                  like(family = "Gaussian",
                       modelGFormula,
                       data = data2_4_NA))
save(fitG_4_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitG/fitG_4_NA.RData")


# Predicting with-held phenotypes
predictG_1 <- predict(fitG_1_NA, newdata=data2_1,
                       formula= ~ fixed_effects  + animal)

predictG_2 <- predict(fitG_2_NA, newdata=data2_2,
                       formula= ~ fixed_effects  + animal)

predictG_3 <- predict(fitG_3_NA, newdata=data2_3,
                       formula= ~ fixed_effects + animal)

predictG_4 <- predict(fitG_4_NA, newdata=data2_4,
                       formula= ~ fixed_effects  + animal)


# Accuracy of fitG

Accu1 <- round(cor(predictG_1$milkZ,predictG_1$mean),2)
Accu1 #0. 0.27

Accu2 <- round(cor(predictG_2$milkZ,predictG_2$mean),2)
Accu2 #0.31

Accu3 <- round(cor(predictG_3$milkZ,predictG_3$mean),2)
Accu3 #0.34

Accu4 <- round(cor(predictG_4$milkZ,predictG_4$mean),2)
Accu4 # 0.40
# Cross-validation accuracy of fitG
Accu_G <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_G #0.33

#-----------------Cross-validation fitGH-----------------------------------------
# Fitting models on data having missing data
fitGH_1_NA <- bru(modelGH,
                   like(family = "Gaussian",
                        modelGHFormula,
                        data = data2_1_NA))
save(fitGH_1_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGH/fitGH_1_NA.RData")

fitGH_2_NA <- bru(modelGH,
                   like(family = "Gaussian",
                        modelGHFormula,
                        data = data2_2_NA))
save(fitGH_2_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGH/fitGH_2_NA.RData")


fitGH_3_NA <- bru(modelGH,
                   like(family = "Gaussian",
                        modelGHFormula,
                        data = data2_3_NA))
save(fitGH_3_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGH/fitGH_3_NA.RData")


fitGH_4_NA <- bru(modelGH,
                   like(family = "Gaussian",
                        modelGHFormula,
                        data = data2_4_NA))
save(fitGH_4_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGH/fitGH_4_NA.RData")


# Predicting with-held phenotypes
predictGH_1 <- predict(fitGH_1_NA, newdata=data2_1,
                        formula= ~ fixed_effects + animal + herd)

predictGH_2 <- predict(fitGH_2_NA, newdata=data2_2,
                        formula= ~ fixed_effects + animal + herd)

predictGH_3 <- predict(fitGH_3_NA, newdata=data2_3,
                        formula= ~ fixed_effects + animal + herd)

predictGH_4 <- predict(fitGH_4_NA, newdata=data2_4,
                        formula= ~ fixed_effects + animal + herd)


# Accuracy of fitGH

Accu1 <- round(cor(predictGH_1$milkZ,predictGH_1$mean),2)
Accu1 #0.36

Accu2 <- round(cor(predictGH_2$milkZ,predictGH_2$mean),2)
Accu2 # 0.45

Accu3 <- round(cor(predictGH_3$milkZ,predictGH_3$mean),2)
Accu3 #0.45

Accu4 <- round(cor(predictGH_4$milkZ,predictGH_4$mean),2)
Accu4 #0.35
# Cross-validation accuracy of fitGH
Accu_GH <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GH #0.40
#-----------------Cross-validation_breed_proportion fitGS-----------------------------------------
# Fitting models on data having missing data
fitGS_1_NA <- bru(modelGS,
                   like(family = "Gaussian",
                        modelGSFormula,
                        data = data2_1_NA))
save(fitGS_1_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGS/fitGS_1_NA.RData")

fitGS_2_NA <- bru(modelGS,
                   like(family = "Gaussian",
                        modelGSFormula,
                        data = data2_2_NA))
save(fitGS_2_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGS/fitGS_2_NA.RData")

fitGS_3_NA <- bru(modelGS,
                   like(family = "Gaussian",
                        modelGSFormula,
                        data = data2_3_NA))
save(fitGS_3_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGS/fitGS_3_NA.RData")

fitGS_4_NA <- bru(modelGS,
                   like(family = "Gaussian",
                        modelGSFormula,
                        data = data2_4_NA))
save(fitGS_4_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGS/fitGS_4_NA.RData")

# Predicting with-held phenotypes
predictGS_1 <- predict(fitGS_1_NA, newdata=data2_1,
                        formula= ~ fixed_effects + animal + field)

predictGS_2 <- predict(fitGS_2_NA, newdata=data2_2,
                        formula= ~ fixed_effects + animal + field)

predictGS_3 <- predict(fitGS_3_NA, newdata=data2_3,
                        formula= ~ fixed_effects + animal + field)

predictGS_4 <- predict(fitGS_4_NA, newdata=data2_4,
                        formula= ~ fixed_effects + animal + field)

#-----------------Cross-validation_breed_proportion fitGS-----------------------------------------

Accu1 <- round(cor(predictGS_1$milkZ,predictGS_1$mean),2)
Accu1 #  0.53
Accu2 <- round(cor(predictGS_2$milkZ,predictGS_2$mean),2)
Accu2 # 0.62
Accu3 <- round(cor(predictGS_3$milkZ,predictGS_3$mean),2)
Accu3 # 0.60
Accu4 <- round(cor(predictGS_4$milkZ,predictGS_4$mean),2)
Accu4 # 0.48
# Cross-validation accuracy of fitGS
Accu_GS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GS # 0.56
#-----------------Cross-validation_breed_proportion fitGHS-----------------------------------------
# Fitting models on data having missing data
fitGHS_1_NA <- bru(modelGHS,
                    like(family = "Gaussian",
                         modelGHSFormula,
                         data = data2_1_NA))
save(fitGHS_1_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGHS/fitGHS_1_NA.RData")

fitGHS_2_NA <- bru(modelGHS,
                    like(family = "Gaussian",
                         modelGHSFormula,
                         data = data2_2_NA))
save(fitGHS_2_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGHS/fitGHS_2_NA.RData")

fitGHS_3_NA <- bru(modelGHS,
                    like(family = "Gaussian",
                         modelGHSFormula,
                         data = data2_3_NA))
save(fitGHS_3_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGHS/fitGHS_3_NA.RData")

fitGHS_4_NA <- bru(modelGHS,
                    like(family = "Gaussian",
                         modelGHSFormula,
                         data = data2_4_NA))
save(fitGHS_4_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGHS/fitGHS_4_NA.RData")


# Predicting with-held phenotypes
predictGHS_1 <- predict(fitGHS_1_NA, newdata=data2_1,
                         formula= ~ fixed_effects + animal + herd +field)

predictGHS_2 <- predict(fitGHS_2_NA, newdata=data2_2,
                         formula= ~ fixed_effects + animal + herd + field)

predictGHS_3 <- predict(fitGHS_3_NA, newdata=data2_3,
                         formula= ~ fixed_effects + animal + herd + field)

predictGHS_4 <- predict(fitGHS_4_NA, newdata=data2_4,
                         formula= ~ fixed_effects + animal + herd + field)
# Accuracy of fitGHS
Accu1 <- round(cor(predictGHS_1$milkZ,predictGHS_1$mean),2)
Accu1 # 0.57
Accu2 <- round(cor(predictGHS_2$milkZ,predictGHS_2$mean),2)
Accu2 # 0.65
Accu3 <- round(cor(predictGHS_3$milkZ,predictGHS_3$mean),2)
Accu3 # 0.64
Accu4 <- round(cor(predictGHS_4$milkZ,predictGHS_4$mean),2)
Accu4 # 0.53
# Cross-validation accuracy of fitGH
Accu_GHS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GHS #  0.60
#-----------------Cross_validation_Region_Paper_G Models-----------------------------------------
# Masking phenotypes of cows in region
#Creating data with missing penotypes
#region 1 masked
data2_1_NA<- data2 %>% mutate(milkZ = ifelse(region == "1", NA, milkZ))
#region 2 masked
data2_2_NA<- data2 %>% mutate(milkZ = ifelse(region == "2", NA, milkZ))
#region3 masked
data2_3_NA<- data2 %>% mutate(milkZ = ifelse(region == "3", NA, milkZ))
#region4 masked
data2_4_NA<- data2 %>% mutate(milkZ = ifelse(region == "4", NA, milkZ))
summary(data2_1_NA$milkZ)
#Creating data for each region
data2_1 <- subset(data2, region=="1") # 437 cows and 4547 record
data2_2 <- subset(data2, region=="2") # 775 cows and 8482 records
data2_3<- subset(data2, region=="3") # 414 Cows and 3674 records
data2_4<- subset(data2, region=="4") # 268 cows and 2672 records
dim(data2_4)
length(unique(data2_4$cow))

#-----------------Cross-validation_region fitG-----------------------------------------
#----------------Models without Permanent env effect----------------------------------------
# Fitting models on data having missing data
fitG_1_NA <- bru(modelG,
                  like(family = "Gaussian",
                       modelGFormula,
                       data = data2_1_NA))
save(fitG_1_NA,file = "data/cleaned_data/cross_validation_region/fitG/fitG_1_NA.RData")

fitG_2_NA <- bru(modelG,
                  like(family = "Gaussian",
                       modelGFormula,
                       data = data2_2_NA))
save(fitG_2_NA,file = "data/cleaned_data/cross_validation_region/fitG/fitG_2_NA.RData")

fitG_3_NA <- bru(modelG,
                  like(family = "Gaussian",
                       modelGFormula,
                       data = data2_3_NA))
save(fitG_3_NA,file = "data/cleaned_data/cross_validation_region/fitG/fitG_3_NA.RData")

fitG_4_NA <- bru(modelG,
                  like(family = "Gaussian",
                       modelGFormula,
                       data = data2_4_NA))
save(fitG_4_NA,file = "data/cleaned_data/cross_validation_region/fitG/fitG_4_NA.RData")

# Predicting with-held phenotypes
predictG_1 <- predict(fitG_1_NA, newdata=data2_1,
                       formula= ~ fixed_effects  + animal)

predictG_2 <- predict(fitG_2_NA, newdata=data2_2,
                       formula= ~ fixed_effects + animal)

predictG_3 <- predict(fitG_3_NA, newdata=data2_3,
                       formula= ~ fixed_effects + animal)

predictG_4 <- predict(fitG_4_NA, newdata=data2_4,
                       formula= ~ fixed_effects + animal)

# Accuracy-region of fitG

Accu1 <- round(cor(predictG_1$milkZ,predictG_1$mean),2)
Accu1 #0.06

Accu2 <- round(cor(predictG_2$milkZ,predictG_2$mean),2)
Accu2 #0.07

Accu3 <- round(cor(predictG_3$milkZ,predictG_3$mean),2)
Accu3 #0.27

Accu4 <- round(cor(predictG_4$milkZ,predictG_4$mean),2)
Accu4 #0.15
# Cross-validation accuracy of fitG
Accu_G <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_G #0.11

#-----------------Cross-validation_region fitGH---------------------------------
# Fitting models on data having missing data
fitGH_1_NA <- bru(modelGH,
                   like(family = "Gaussian",
                        modelGHFormula,
                        data = data2_1_NA))
save(fitGH_1_NA,file = "data/cleaned_data/cross_validation_region/fitGH/fitGH_1_NA.RData")
summary(data2_1_NA$dgrp)
fitGH_2_NA <- bru(modelGH,
                   like(family = "Gaussian",
                        modelGHFormula,
                        data = data2_2_NA))
save(fitGH_2_NA,file = "data/cleaned_data/cross_validation_region/fitGH/fitGH_2_NA.RData")


fitGH_3_NA <- bru(modelGH,
                   like(family = "Gaussian",
                        modelGHFormula,
                        data = data2_3_NA))
save(fitGH_3_NA,file = "data/cleaned_data/cross_validation_region/fitGH/fitGH_3_NA.RData")


fitGH_4_NA <- bru(modelGH,
                   like(family = "Gaussian",
                        modelGHFormula,
                        data = data2_4_NA))
save(fitGH_4_NA,file = "data/cleaned_data/cross_validation_region/fitGH/fitGH_4_NA.RData")


# Predicting with-held phenotypes
predictGH_1 <- predict(fitGH_1_NA, newdata=data2_1,
                        formula= ~ fixed_effects + animal + herd)

predictGH_2 <- predict(fitGH_2_NA, newdata=data2_2,
                        formula= ~ fixed_effects  + animal + herd)

predictGH_3 <- predict(fitGH_3_NA, newdata=data2_3,
                        formula= ~ fixed_effects + animal + herd)

predictGH_4 <- predict(fitGH_4_NA, newdata=data2_4,
                        formula= ~ fixed_effects  + animal + herd)

# Accuracy-region of fitGH

Accu1 <- round(cor(predictGH_1$milkZ,predictGH_1$mean),2)
Accu1 # 0.10
Accu2 <- round(cor(predictGH_2$milkZ,predictGH_2$mean),2)
Accu2 #0.10
Accu3 <- round(cor(predictGH_3$milkZ,predictGH_3$mean),2)
Accu3 #0.26
Accu4 <- round(cor(predictGH_4$milkZ,predictGH_4$mean),2)
Accu4 #0.17
# Cross-validation accuracy of fitGH
Accu_GH <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GH #0.16

#-----------------Cross-validation_region fitGS---------------------------------
# Fitting models on data having missing data
fitGS_1_NA <- bru(modelGS,
                   like(family = "Gaussian",
                        modelGSFormula,
                        data = data2_1_NA))
save(fitGS_1_NA,file = "data/cleaned_data/cross_validation_region/fitGS/fitGS_1_NA.RData")

fitGS_2_NA <- bru(modelGS,
                   like(family = "Gaussian",
                        modelGSFormula,
                        data = data2_2_NA))
save(fitGS_2_NA,file = "data/cleaned_data/cross_validation_region/fitGS/fitGS_2_NA.RData")

fitGS_3_NA <- bru(modelGS,
                   like(family = "Gaussian",
                        modelGSFormula,
                        data = data2_3_NA))
save(fitGS_3_NA,file = "data/cleaned_data/cross_validation_region/fitGS/fitGS_3_NA.RData")

fitGS_4_NA <- bru(modelGS,
                   like(family = "Gaussian",
                        modelGSFormula,
                        data = data2_4_NA))
save(fitGS_4_NA,file = "data/cleaned_data/cross_validation_region/fitGS/fitGS_4_NA.RData")

# Predicting with-held phenotypes
predictGS_1 <- predict(fitGS_1_NA, newdata=data2_1,
                        formula= ~ fixed_effects  + animal + field)

predictGS_2 <- predict(fitGS_2_NA, newdata=data2_2,
                        formula= ~ fixed_effects + animal + field)

predictGS_3 <- predict(fitGS_3_NA, newdata=data2_3,
                        formula= ~ fixed_effects + animal + field)

predictGS_4 <- predict(fitGS_4_NA, newdata=data2_4,
                        formula= ~ fixed_effects + animal + field)

# Accuracy_region of fitGS
Accu1 <- round(cor(predictGS_1$milkZ,predictGS_1$mean),2)
Accu1 # 0.02
Accu2 <- round(cor(predictGS_2$milkZ,predictGS_2$mean),2)
Accu2 # 0.01
Accu3 <- round(cor(predictGS_3$milkZ,predictGS_3$mean),2)
Accu3 # 0.24
Accu4 <- round(cor(predictGS_4$milkZ,predictGS_4$mean),2)
Accu4 # 0.15
# Cross-validation region accuracy of fitGS
Accu_GS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GS # 0.13
#-----------------Cross-validation_region fitGHS-----------------------------------------
# Fitting models on data having missing data
fitGHS_1_NA <- bru(modelGHS,
                    like(family = "Gaussian",
                         modelGHSFormula,
                         data = data2_1_NA))
save(fitGHS_1_NA,file = "data/cleaned_data/cross_validation_region/fitGHS/fitGHS_1_NA.RData")

fitGHS_2_NA <- bru(modelGHS,
                    like(family = "Gaussian",
                         modelGHSFormula,
                         data = data2_2_NA))
save(fitGHS_2_NA,file = "data/cleaned_data/cross_validation_region/fitGHS/fitGHS_2_NA.RData")

fitGHS_3_NA <- bru(modelGHS,
                    like(family = "Gaussian",
                         modelGHSFormula,
                         data = data2_3_NA))
save(fitGHS_3_NA,file = "data/cleaned_data/cross_validation_region/fitGHS/fitGHS_3_NA.RData")

fitGHS_4_NA <- bru(modelGHS,
                    like(family = "Gaussian",
                         modelGHSFormula,
                         data = data2_4_NA))
save(fitGHS_4_NA,file = "data/cleaned_data/cross_validation_region/fitGHS/fitGHS_4_NA.RData")

# Predicting with-held phenotypes
predictGHS_1 <- predict(fitGHS_1_NA, newdata=data2_1,
                         formula= ~ fixed_effects + animal + herd +field)

predictGHS_2 <- predict(fitGHS_2_NA, newdata=data2_2,
                         formula= ~ fixed_effects + animal + herd + field)

predictGHS_3 <- predict(fitGHS_3_NA, newdata=data2_3,
                         formula= ~ fixed_effects  + animal + herd + field)

predictGHS_4 <- predict(fitGHS_4_NA, newdata=data2_4,
                        formula= ~ fixed_effects + animal + herd + field)
#rm (predictGHS_1, predictGHS_2,predictGHS_3,predictGHS_4,predictGS_1, predictGS_2,predictGS_3,predictGS_4)

# Accuracy regions of fitGHS

Accu1 <- round(cor(predictGHS_1$milkZ,predictGHS_1$mean),2)
Accu1 # 0.11
Accu2 <- round(cor(predictGHS_2$milkZ,predictGHS_2$mean),2)
Accu2 # -0.02

Accu3 <- round(cor(predictGHS_3$milkZ,predictGHS_3$mean),2)
Accu3 # 0.18

Accu4 <- round(cor(predictGHS_4$milkZ,predictGHS_4$mean),2)
Accu4 # 0.14
# Cross-validation accuracy of fitGHS
Accu_GHS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GHS # 0.10
#-------- Foward validation G models_Paper----------------------------------------------------
#--Predicting phenotypes of young cows born in 2016(141 cows) and 2017(5 cows)--
#importing cow IDS and birth dates
cowbdate<- read.table(file = "data/original_data/isibdate.txt",header = FALSE)
# isibdate.txt (see file description in README.txt for data)
colnames(cowbdate) <- c("cow","birthdaymonthyear")
summary(cowbdate$birthdaymonthyear)
# extract characters from 5th index to 8th index
cowbdate['birthyear'] <-str_sub(cowbdate$birthdaymonthyear, 5, 8)  # prints "last 4 characters ie year"
cowbdate <- cowbdate[c('cow','birthyear')]
str(cowbdate)
table(cowbdate$birthyear)
# 007  008  010  011  012  013  014  015  016  017 2004 2005 2006 2007 2008 2009
#2    1    4    1   90    7  319   33   58    2    1    1    1    2    5    2
#2010 2011 2012 2013 2014 2015 2016 2017
#10   22  219   93  860   97   83    3

# Standardise year format
# Replace Values Based on Condition
cowbdate$birthyear[cowbdate$birthyear == "014"] <- "2014"
cowbdate$birthyear[cowbdate$birthyear == "007"] <- "2007"
cowbdate$birthyear[cowbdate$birthyear == "008"] <- "2008"
cowbdate$birthyear[cowbdate$birthyear == "010"] <- "2010"
cowbdate$birthyear[cowbdate$birthyear == "011"] <- "2011"
cowbdate$birthyear[cowbdate$birthyear == "012"] <- "2012"
cowbdate$birthyear[cowbdate$birthyear == "013"] <- "2013"
cowbdate$birthyear[cowbdate$birthyear == "015"] <- "2015"
cowbdate$birthyear[cowbdate$birthyear == "016"] <- "2016"
cowbdate$birthyear[cowbdate$birthyear == "017"] <- "2017"
table(cowbdate$birthyear)
#2016 2017
#141  5
length(unique(cowbdate$cow))

# Create birthyear column in data2 by merging with cowbdate

data2 <-  merge(data2, cowbdate, by= "cow", all.x = TRUE)
summary(data2$birthyear)
#2004  2005  2006  2007  2008  2009  2010  2011  2012  2013  2014  2015  2016  2017
#17    12     7    38    41     7   122   232  3044  1066 12017  1298  1444    30
data2$birthyear <- as.factor(data2$birthyear)
summary(data2$birthyear)
#data2_2016_2017 <- subset(data1_birthdate, birthyear==2016 | birthyear==2017 )

#Let's make NA milkz of cows born in 2016 and 2017
data2_NA<- data2%>% mutate(milkZ = ifelse(birthyear=="2016" | birthyear=="2017", NA, milkZ))
sum(is.na(data2_NA$milkZ))   #1424 Expected
length(unique(data2_NA$cow)) #1894 Expected

# Let's prepare data for young cows
data2_young <- subset(data2, birthyear=="2016" | birthyear=="2017") # 146 cows 1424 records
length(unique(data2_young$cow)) #146
#1  2  3  4
#60 56 24  6

data2_young_sub <- data2_young[, c("cow","region","dgrp")]
data2_young_sub <- distinct(data2_young_sub)
table(data2_young_sub$dgrp)
table(data2_young_sub$region)

# Fitting models on data having missing data
fitG_NA <- bru(modelG,
                like(family = "Gaussian",
                     modelGFormula,
                     data = data2_NA))
save(fitG_NA,file = "data/cleaned_data/foward_validation/fitG/fitG_NA.RData")

fitGH_NA <- bru(modelGH,
                 like(family = "Gaussian",
                      modelGHFormula,
                      data = data2_NA))
save(fitGH_NA,file = "data/cleaned_data/foward_validation/fitGH/fitGH_NA.RData")

fitGS_NA <- bru(modelGS,
                 like(family = "Gaussian",
                      modelGSFormula,
                      data = data2_NA))
save(fitGS_NA,file = "data/cleaned_data/foward_validation/fitGS/fitGS_NA.RData")

fitGHS_NA <- bru(modelGHS,
                  like(family = "Gaussian",
                       modelGHSFormula,
                       data = data2_NA))
save(fitGHS_NA,file = "data/cleaned_data/foward_validation/fitGHS/fitGHS_NA.RData")

# Predicting with-held phenotypes
predictG <- predict(fitG_NA, newdata=data2_young,
                     formula= ~ fixed_effects + animal)
predictGH <- predict(fitGH_NA, newdata=data2_young,
                      formula= ~ fixed_effects + animal + herd)
predictGS <- predict(fitGS_NA, newdata=data2_young,
                      formula= ~ fixed_effects + animal + field)
predictGHS <- predict(fitGHS_NA, newdata=data2_young,
                       formula= ~ fixed_effects + animal + herd +field)
# Accuracy forward validation G
AccuG <- round(cor(predictG$milkZ,predictG$mean),2)
AccuG 
# Accuracy forward validation GH
AccuGH <- round(cor(predictGH$milkZ,predictGH$mean),2)
AccuGH
# Accuracy forward validation GS
AccuGS <- round(cor(predictGS$milkZ,predictGS$mean),2)
AccuGS
# Accuracy forward validation GHS
AccuGHS <- round(cor(predictGHS$milkZ,predictGHS$mean),2)
AccuGHS
#----Summarise forward validation accuracy by dgrp and region for G models-------------------
# dgrp
#Fit G dgrp
predictG <- predictG[, c("cow","region","dgrp", "milkZ", "mean")]
predictG_dgrp1 <- subset(predictG, dgrp=="1")
length(unique(predictG_dgrp1$cow)) # 60 cows

predictG_dgrp2 <- subset(predictG, dgrp=="2")
length(unique(predictG_dgrp2$cow)) #  56 cows

predictG_dgrp3 <- subset(predictG, dgrp=="3")
length(unique(predictG_dgrp3$cow)) # 24 cows

predictG_dgrp4 <- subset(predictG, dgrp=="4")
length(unique(predictG_dgrp4$cow)) #  6 cows

Accu_dgrp1 <-  round(cor(predictG_dgrp1$milkZ,predictG_dgrp1$mean),2)
Accu_dgrp1 # 0.29

Accu_dgrp2 <-  round(cor(predictG_dgrp2$milkZ,predictG_dgrp2$mean),2)
Accu_dgrp2 #0.28

Accu_dgrp3 <-  round(cor(predictG_dgrp3$milkZ,predictG_dgrp3$mean),2)
Accu_dgrp3 # 0.54

Accu_dgrp4 <-  round(cor(predictG_dgrp4$milkZ,predictG_dgrp4$mean),2)
Accu_dgrp4 # 0.83

Accu_dgrp <- round((Accu_dgrp1 + Accu_dgrp2 + Accu_dgrp3 + Accu_dgrp4)/4,2)
Accu_dgrp # 0.49

#Fit GH dgrp
predictGH <- predictGH[, c("cow","region","dgrp", "milkZ", "mean")]
predictGH_dgrp1 <- subset(predictGH, dgrp=="1")
length(unique(predictGH_dgrp1$cow)) # 60 cows

predictGH_dgrp2 <- subset(predictGH, dgrp=="2")
length(unique(predictGH_dgrp2$cow)) #  56 cows

predictGH_dgrp3 <- subset(predictGH, dgrp=="3")
length(unique(predictGH_dgrp3$cow)) # 24 cows

predictGH_dgrp4 <- subset(predictGH, dgrp=="4")
length(unique(predictGH_dgrp4$cow)) #  6 cows

Accu_dgrp1 <-  round(cor(predictGH_dgrp1$milkZ,predictGH_dgrp1$mean),2)
Accu_dgrp1 # 0.36

Accu_dgrp2 <-  round(cor(predictGH_dgrp2$milkZ,predictGH_dgrp2$mean),2)
Accu_dgrp2 # 0.63

Accu_dgrp3 <-  round(cor(predictGH_dgrp3$milkZ,predictGH_dgrp3$mean),2)
Accu_dgrp3 # 0.66

Accu_dgrp4 <-  round(cor(predictGH_dgrp4$milkZ,predictGH_dgrp4$mean),2)
Accu_dgrp4 # 0.33

Accu_dgrp <- round((Accu_dgrp1 + Accu_dgrp2 + Accu_dgrp3 + Accu_dgrp4)/4,2)
Accu_dgrp # 0.50

#Fit GS dgrp

predictGS <- predictGS[, c("cow","region","dgrp", "milkZ", "mean")]
predictGS_dgrp1 <- subset(predictGS, dgrp=="1")
length(unique(predictGS_dgrp1$cow)) # 60 cows

predictGS_dgrp2 <- subset(predictGS, dgrp=="2")
length(unique(predictGS_dgrp2$cow)) #  56 cows

predictGS_dgrp3 <- subset(predictGS, dgrp=="3")
length(unique(predictGS_dgrp3$cow)) # 24 cows

predictGS_dgrp4 <- subset(predictGS, dgrp=="4")
length(unique(predictGS_dgrp4$cow)) #  6 cows

Accu_dgrp1 <-  round(cor(predictGS_dgrp1$milkZ,predictGS_dgrp1$mean),2)
Accu_dgrp1 # 0.37

Accu_dgrp2 <-  round(cor(predictGS_dgrp2$milkZ,predictGS_dgrp2$mean),2)
Accu_dgrp2 # 0.59
Accu_dgrp3 <-  round(cor(predictGS_dgrp3$milkZ,predictGS_dgrp3$mean),2)
Accu_dgrp3 # 0.73
Accu_dgrp4 <-  round(cor(predictGS_dgrp4$milkZ,predictGS_dgrp4$mean),2)
Accu_dgrp4 # 0.74
Accu_dgrp <- round((Accu_dgrp1 + Accu_dgrp2 + Accu_dgrp3 + Accu_dgrp4)/4,2)
Accu_dgrp # 0.61

#Fit GHS dgrp
predictGHS <- predictGHS[, c("cow","region","dgrp", "milkZ", "mean")]
predictGHS_dgrp1 <- subset(predictGHS, dgrp=="1")
length(unique(predictGHS_dgrp1$cow)) # 60 cows

predictGHS_dgrp2 <- subset(predictGHS, dgrp=="2")
length(unique(predictGHS_dgrp2$cow)) #  56 cows

predictGHS_dgrp3 <- subset(predictGHS, dgrp=="3")
length(unique(predictGHS_dgrp3$cow)) # 24 cows

predictGHS_dgrp4 <- subset(predictGHS, dgrp=="4")
length(unique(predictGHS_dgrp4$cow)) #  6 cows

Accu_dgrp1 <-  round(cor(predictGHS_dgrp1$milkZ,predictGHS_dgrp1$mean),2)
Accu_dgrp1 # 0.44

Accu_dgrp2 <-  round(cor(predictGHS_dgrp2$milkZ,predictGHS_dgrp2$mean),2)
Accu_dgrp2 # 0.67
Accu_dgrp3 <-  round(cor(predictGHS_dgrp3$milkZ,predictGHS_dgrp3$mean),2)
Accu_dgrp3 # 0.75

Accu_dgrp4 <-  round(cor(predictGHS_dgrp4$milkZ,predictGHS_dgrp4$mean),2)
Accu_dgrp4 # 0.66

Accu_dgrp <- round((Accu_dgrp1 + Accu_dgrp2 + Accu_dgrp3 + Accu_dgrp4)/4,2)
Accu_dgrp # 0.63
# region
#Fit G region
predictG <- predictG[, c("cow","region","dgrp", "milkZ", "mean")]
predictG_region1 <- subset(predictG, region=="1")
length(unique(predictG_region1$cow)) # 30 cows

predictG_region2 <- subset(predictG, region=="2")
length(unique(predictG_region2$cow)) #  60 cows

predictG_region3 <- subset(predictG, region=="3")
length(unique(predictG_region3$cow)) # 20 cows

predictG_region4 <- subset(predictG, region=="4")
length(unique(predictG_region4$cow)) #  36 cows

Accu_region1 <-  round(cor(predictG_region1$milkZ,predictG_region1$mean),2)
Accu_region1 # 0.38

Accu_region2 <-  round(cor(predictG_region2$milkZ,predictG_region2$mean),2)
Accu_region2 # 0.46

Accu_region3 <-  round(cor(predictG_region3$milkZ,predictG_region3$mean),2)
Accu_region3 # 0.11

Accu_region4 <-  round(cor(predictG_region4$milkZ,predictG_region4$mean),2)
Accu_region4 # 0.08

Accu_region <- round((Accu_region1 + Accu_region2 + Accu_region3 + Accu_region4)/4,2)
Accu_region # 0.26


#Fit GH region

predictGH <- predictGH[, c("cow","region","dgrp", "milkZ", "mean")]
predictGH_region1 <- subset(predictGH, region=="1")
length(unique(predictGH_region1$cow)) # 30 cows

predictGH_region2 <- subset(predictGH, region=="2")
length(unique(predictGH_region2$cow)) #  60 cows

predictGH_region3 <- subset(predictGH, region=="3")
length(unique(predictGH_region3$cow)) # 20 cows

predictGH_region4 <- subset(predictGH, region=="4")
length(unique(predictGH_region4$cow)) #  36 cows

Accu_region1 <-  round(cor(predictGH_region1$milkZ,predictGH_region1$mean),2)
Accu_region1 # 0.65

Accu_region2 <-  round(cor(predictGH_region2$milkZ,predictGH_region2$mean),2)
Accu_region2 # 0.57

Accu_region3 <-  round(cor(predictGH_region3$milkZ,predictGH_region3$mean),2)
Accu_region3 # 0.32

Accu_region4 <-  round(cor(predictGH_region4$milkZ,predictGH_region4$mean),2)
Accu_region4 # 0.35
Accu_region <- round((Accu_region1 + Accu_region2 + Accu_region3 + Accu_region4)/4,2)
Accu_region # 0.47
#Fit GS region
predictGS <- predictGS[, c("cow","region","region", "milkZ", "mean")]
predictGS_region1 <- subset(predictGS, region=="1")
length(unique(predictGS_region1$cow)) # 30 cows

predictGS_region2 <- subset(predictGS, region=="2")
length(unique(predictGS_region2$cow)) #  60 cows

predictGS_region3 <- subset(predictGS, region=="3")
length(unique(predictGS_region3$cow)) # 20 cows
predictGS_region4 <- subset(predictGS, region=="4")
length(unique(predictGS_region4$cow)) #  36 cows
Accu_region1 <-  round(cor(predictGS_region1$milkZ,predictGS_region1$mean),2)
Accu_region1 # 0.57
Accu_region2 <-  round(cor(predictGS_region2$milkZ,predictGS_region2$mean),2)
Accu_region2 # 0.53
Accu_region3 <-  round(cor(predictGS_region3$milkZ,predictGS_region3$mean),2)
Accu_region3 # 0.28
Accu_region4 <-  round(cor(predictGS_region4$milkZ,predictGS_region4$mean),2)
Accu_region4 # 0.37
Accu_region <- round((Accu_region1 + Accu_region2 + Accu_region3 + Accu_region4)/4,2)
Accu_region # 0.44

#Fit GHS region
predictGHS <- predictGHS[, c("cow","region","region", "milkZ", "mean")]
predictGHS_region1 <- subset(predictGHS, region=="1")
length(unique(predictGHS_region1$cow)) # 30 cows

predictGHS_region2 <- subset(predictGHS, region=="2")
length(unique(predictGHS_region2$cow)) #  60 cows

predictGHS_region3 <- subset(predictGHS, region=="3")
length(unique(predictGHS_region3$cow)) # 20 cows

predictGHS_region4 <- subset(predictGHS, region=="4")
length(unique(predictGHS_region4$cow)) #  36 cows

Accu_region1 <-  round(cor(predictGHS_region1$milkZ,predictGHS_region1$mean),2)
Accu_region1 #  0.65
Accu_region2 <-  round(cor(predictGHS_region2$milkZ,predictGHS_region2$mean),2)
Accu_region2 # 0.60
Accu_region3 <-  round(cor(predictGHS_region3$milkZ,predictGHS_region3$mean),2)
Accu_region3 # 0.37
Accu_region4 <-  round(cor(predictGHS_region4$milkZ,predictGHS_region4$mean),2)
Accu_region4 # 0.46
Accu_region <- round((Accu_region1 + Accu_region2 + Accu_region3 + Accu_region4)/4,2)
Accu_region #  0.52
