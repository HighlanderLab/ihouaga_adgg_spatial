# ---- Header ------------------------------------------------------------------
# Spatial modeling
# Trait: Test day milk yield
# Author: Isidore
# Version 1.0.0
# Date: 2023-07-20

# ---- Setup -------------------------------------------------------------------

# Working directory
# ... Isidore's window's laptop
getwd()
baseDir <- "C:/Users/Lenovo/OneDrive/Documents/ihouaga_adgg_spatial"

# ... Isidore's Mack laptop
baseDir <- "/Users/ihouaga2/ihouaga_adgg_spatial"

# ... Isidore's Eddie workspace
baseDir <- ""

# ... path on ILRI hpc
baseDir <- ""

# ... Gregor's office computer
baseDir <- "/Users/ggorjanc/Storages/GitBox/HighlanderLab/ihouaga_adgg_spatial_upstream/"

# Change working directory
setwd(dir = baseDir)
getwd()
dir()

# ---- Installing and load packages --------------------------------------------

if (FALSE) {
  requiredPackages <- c(
    "tidyverse", # for data manipulation
    "fmesher",
    "inlabru",
    "verification", # for Continuous Ranked Probability Score
    "irlba", # for fast PCA
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
library(fmesher)
library(gridExtra) # Visualize random field grid
library(verification) # Visualize random field grid
library(lattice)
library(raster) # GetData to plot country map
library(rgeos)
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

# Read in the phenotype data
data1 <- read.csv(file = "data/cleaned_data/milk_yield_pheno_cleaned.csv")
str(data1)

# Factors
data1 <- data1 %>%
  mutate_at(.vars = c("cow", "ward", "herd", "cyrsn", "tyrmn", "dgrp", "lac",
                      "lacgr", "ward_code", "region"),
            .funs = as.factor)
str(data1)
summary(data1$herd)
data1$herdI <- as.numeric(data1$herd) # these codes will now be 1:n
data1$ward_codeI <- as.numeric(data1$ward_code) # these codes will now be 1:n
data1$cowPe <- data1$cow # cow permanent environment
data1$cowPeI <- as.numeric(data1$cowPe) # these codes will now be 1:n
# data1$cowI <- as.numeric(data1$????) # see below!
data1$tyrmnI <- as.numeric(data1$tyrmn) # these codes will now be 1:n
data1$cyrsnI <- as.numeric(data1$cyrsn) # these codes will now be 1:n
data1$dgrpI <- as.numeric(data1$dgrp) # these codes will now be 1:n
# I made cyrn and tyrnm as numeric because for SPDE models I got error when they are factor (# I got: "Error in cyrsn + tyrmn : non-conformable arrays" and Error in cyrsn + tyrmn : non-numeric argument to binary operator"
#After making them as.numeric the spde model run well # Discuss this with @gg)
data1$regionI <- as.numeric(data1$region) #these codes will now be 1:n

# Read in the regional data for Besag model
nb.map <- inla.read.graph(filename = "data/cleaned_data/ward_neighbours.txt")
load(file = "data/cleaned_data/ward_neighbours_precision_matrix.RData")
# this will load nb.matrix and nb.matrixScaled

# Read in the genomic relationship matrix
load(file = "data/cleaned_data/GRMInv.RData")
str(GRMInv)
dim(GRMInv) # 1894 x 1894
class(GRMInv)
head(GRMInv)
#rm(nb.map,nb.matrix, nb.matrixScaled)
# Now we can code the cows in pheno data correctly
data1$cow <- factor(data1$cow, levels = 1:nrow(GRMInv))
summary(as.numeric(levels(data1$cow))) # we need 1:1911 numbers here!
data1$cowI <- as.numeric(as.character(data1$cow))
summary(data1$cowI) # we have 1:1894 ids here
head(data1)
tail(data1)
# TODO double check that the above cowI is correct (GG thinks it is)
# Standardise response variable and covariates
# ... INLA uses priors so best to be on the O(1) scale
data1$milkZ <- scale(data1$milk)
data1$ageZ <- scale(data1$age)

summary(data1$milk)
summary(data1$milkZ)

summary(data1$age)
summary(data1$ageZ)

summary(data1$leg0)
summary(data1$leg1)
summary(data1$leg2)

# -------Building mesh and prepare data for SPDE modelling---------------------------------------------------------
# Building the mesh
mapTZA <- getData('GADM', country = "TZA", level= 1) # Get map of Tanzania
# Extract the border from the map
TanzaniaBorder <- gUnaryUnion(mapTZA, rep(1, nrow(mapTZA)))
# Formatting boundary for r-INLA
TanzaniaBorder <- inla.sp2segment(TanzaniaBorder)
locations <- cbind(data1$long, data1$lat)
plot(locations)
bnd <- fm_nonconvex_hull_inla(locations, 0.5)
TanzaniaBorder <- fm_as_segm(TanzaniaBorder)

mesh <- fm_mesh_2d_inla(boundary = list(bnd, TanzaniaBorder),
                        max.edge = c(0.2, 1), cutoff = 0.08)

#mesh2 <- fm_mesh_2d_inla(boundary = list(bnd, TanzaniaBorder),
                        #max.edge = c(0.2, 1), cutoff = 0.3)

#mesh3 <- fm_mesh_2d_inla(boundary = list(bnd, TanzaniaBorder),
                        # max.edge = c(0.1, 1), cutoff = 0.3)

#mesh4 <- fm_mesh_2d_inla(boundary = list(bnd, TanzaniaBorder),
                         #max.edge = c(0.1, 0.5), cutoff = 0.3)


#mesh5 <- fm_mesh_2d_inla(boundary = list(bnd, TanzaniaBorder),
                         #max.edge = c(0.1, 0.4), cutoff = 0.3)


#mesh6 <- fm_mesh_2d_inla(boundary = list(bnd, TanzaniaBorder),
                        # max.edge = c(0.1, 0.3), cutoff = 0.3)

#mesh7 <- fm_mesh_2d_inla(boundary = list(bnd, TanzaniaBorder),
                        # max.edge = c(0.1, 0.2), cutoff = 0.3)



if (FALSE) {
  ggplot() +
    geom_fm(data = mesh) +
    geom_point(aes(locations[,1], locations[,2]))
}

#For INLABru define matern
matern <-
  inla.spde2.pcmatern(mesh,
                      prior.sigma = c(sqrt(0.25), 0.5),
                      prior.range = c(50, 0.8)
  )


data2 <- sf::st_as_sf(x = data1[, c("long", "lat")],
                      coords = c("long", "lat"))
data2$milkZ <- data1$milkZ
data2$milk <- data1$milk
data2$cyrsn <- data1$cyrsn
data2$tyrmn <- data1$tyrmn
data2$dgrp <- data1$dgrp
data2$lacgr <- data1$lacgr
data2$ageZ <- data1$ageZ
data2$age <- data1$age
data2$herd <- data1$herd
data2$ward_code <- data1$ward_code
data2$herdI <- data1$herdI
data2$cow <- data1$cow
data2$cowI <- data1$cowI
data2$cowPe <- data1$cowPe
# data2$cowPeI <- data1$cowPeI
data2$leg0 <- data1$leg0
data2$leg1 <- data1$leg1
data2$leg2 <- data1$leg2
data2$long <- data1$long
data2$lat <- data1$lat
data2$ward_codeI <- data1$ward_codeI
data2$cyrsnI <- data1$cyrsnI
data2$tyrmnI <- data1$tyrmnI
data2$dgrpI <- data1$dgrpI
length(unique(data2$geometry)) #1385 couples of GPS vs 1386 herds (2 herds at same location)

# Create a spatialpointdataframe needed for special effect prediction.
# Following instructions from  https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
# prepare coordinates, data, and proj4string
coords <- data1[ , c("long", "lat")]   # coordinates
data   <- data1          # data
crs    <- CRS("+init=epsg:28992") # proj4string of coords

# make the SpatialPointsDataFrame object
spdf <- SpatialPointsDataFrame(coords      = coords,
                               data        = data,
                               proj4string = crs)
class(spdf) # "SpatialPointsDataFrame"

# ---- Specify models R-INLAbru ----------------------------------------------------------
# Base model for milk - All effects without ward_code
modelBase <- ~ fixed_effects(main = cyrsn + tyrmn + dgrp +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + herd(main =herd, model = 'iid') +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv)

modelBaseFormula <- milkZ ~ .

# Adding ward_code as a fixed effect
modelWCF <- ~ fixed_effects(main = cyrsn + tyrmn + dgrp + ward_code +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + herd(main =herd, model = 'iid') +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv)

modelWCFFormula <- milkZ ~ .

# Adding ward_code as a Random independent
modelWCRI <- ~fixed_effects(main = cyrsn + tyrmn + dgrp +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + herd(main =herd, model = 'iid') +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) + ward(main=ward_code, model = 'iid')

modelWCRIFormula <- milkZ ~ .


# Adding ward_code as a Random Besag
modelWCRB <- ~ fixed_effects(main = cyrsn + tyrmn + dgrp +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + herd(main =herd, model = 'iid') +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) + ward(main=ward_code,  model = 'besag', graph = nb.map, scale.model = TRUE,
                                                                   mapper=bru_mapper_index(3644))

modelWCRBFormula <- milkZ ~ .

# Base model fitBase + Spatial effect (fitS)
modelfitS <- ~ fixed_effects(main = cyrsn + tyrmn + dgrp +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + herd(main =herd, model = 'iid') +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  field(geometry, model = matern)

modelfitSFormula <- milkZ ~ .

# Base model fitBase + Ward effect +   + Spatial effect  (fitWS)

modelWS <- ~
  fixed_effects(main = ~ 1 + cyrsn + tyrmn + dgrp +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + herd(main =herd, model = 'iid') +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  ward(main=ward_code, model = 'iid') +
  field(geometry, model = matern)

modelWSFormula <- milkZ ~ .


# ---- Specifying additional models_inlabru-------------------------------------
#-------------------------INLA_Bru-factors--------------------------------------
# Base model without herd effect  but with ward effect included (fitB)
# ModelBase without herd

modelB <- ~ fixed_effects(main = ~ cyrsn + tyrmn + dgrp +
                                ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + ward(main = ward_code, model = 'iid') +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv)

modelBFormula <- milkZ ~ .


# Base model fitB + herd  (fitBH)
modelBH <- ~ fixed_effects(main = ~ cyrsn + tyrmn + dgrp +
                                 ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + ward(main = ward_code, model = 'iid') +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  herd(main=herd, model = 'iid')

modelBHFormula <- milkZ ~ .

# Base model fitB + Spatial effect (fitBS)
modelBS <- ~ fixed_effects(main = ~ cyrsn + tyrmn + dgrp +
                                 ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + ward(main = ward_code, model = 'iid') +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  field(geometry, model = matern)

modelBSFormula <- milkZ ~ .

# Base model fitB + herd  + Spatial effect  (fitBHS)

modelBHS <- ~fixed_effects(main = ~ cyrsn + tyrmn + dgrp +
                                 ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + ward(main = ward_code, model = 'iid') +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  herd(main=herd, model = 'iid') +
  field(geometry, model = matern)

modelBHSFormula <- milkZ ~ .

#-------------------------INLA_Bru_numeric_factors---------------------------------------
# Base model without herd effect  but with ward effect included (fitB_bru)
# ModelBase without herd

modelB_bru <- ~ fixed_effects(main = ~ cyrsnI + tyrmnI + dgrpI +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + ward(main = ward_code, model = 'iid') +
  animal(main = cowI, model = 'generic', Cmatrix = GRMInv)

modelB_bruFormula <- milkZ ~ .


# Base model fitB + herd  (fitBH_bru)
modelBH_bru <- ~ fixed_effects(main = ~ cyrsnI + tyrmnI + dgrpI +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + ward(main = ward_code, model = 'iid') +
  animal(main = cowI, model = 'generic', Cmatrix = GRMInv) +
  herd(main=herd, model = 'iid')

modelBH_bruFormula <- milkZ ~ .

# Base model fitB + Spatial effect (fitBS_bru)
modelBS_bru <- ~ fixed_effects(main = ~ cyrsnI + tyrmnI + dgrpI +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + ward(main = ward_code, model = 'iid') +
  animal(main = cowI, model = 'generic', Cmatrix = GRMInv) +
  field(geometry, model = matern)

modelBS_bruFormula <- milkZ ~ .

# Base model fitB + herd  + Spatial effect  (fitBHS_bru)

modelBHS_bru <- ~fixed_effects(main = ~ cyrsnI + tyrmnI + dgrpI +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + ward(main = ward_code, model = 'iid') +
  animal(main = cowI, model = 'generic', Cmatrix = GRMInv) +
  herd(main=herd, model = 'iid') +
  field(geometry, model = matern)

modelBHS_bruFormula <- milkZ ~ .

# ---- Specifying additional models_inla----------------------------------------------------------
# Base model without herd effect  but with ward effect included (fitB_inla)
# ModelBase without herd

modelB_inla <- "milkZ ~ 1 + cyrsnI + tyrmnI + dgrpI +
                          (ageZ|lacgr) + (leg1|lacgr) + (leg2|lacgr) +
                          f(ward_code, model = 'iid') +
                          f(cowPe, model = 'iid') +
                         f(cowI, model = 'generic0', Cmatrix = GRMInv)"


# Base model fitB + herd  (fitBH_inla)
modelBH_inla <- as.formula(paste0(modelB_inla, " + f(herd, model = 'iid')"))


# Base model fitB + Spatial effect=fitBS_inla and model fitBH + Spatial effect=fitBHS_inla

# -------SPDE with R-INLA ------------------------------------------------------
# Make mesh and SPDE
# Priors
# SPDE
hyperRange  = c(50, 0.8)
hyperVarSpdeS = c(sqrt(0.25), 0.5)
hyperVarSpdeHS = c(sqrt(0.10), 0.5)
hyperResVarBHS = list(theta = list(prior="pc.prec", param=c(sqrt(0.15),0.5)))


A = inla.spde.make.A(mesh = mesh, loc = cbind(data2$long, data2$lat) )
spdeStatBS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeS)
meshIndexBS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatBS$n.spde)
spdeStatBHS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeHS)
meshIndexBHS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatBHS$n.spde)

# Make stack
# StackBS
stackBS = inla.stack(data = list(milkZ = data2$milkZ),
                    A = list(A,1),
                    effects = list(c(meshIndexBS, list(intercept = 1)),
                                   list(cowI = data2$cowI,ward_code = data2$ward_code, cowPe = data2$cowPe,
                                        cyrsnI=data2$cyrsnI, tyrmnI=data2$tyrmnI, dgrpI= data2$dgrpI, ageZ=data2$ageZ, lacgr=data2$lacgr, leg1=data2$leg1,leg2=data2$leg2)), tag = "data2.data")

# StackBHS (herd as random + Spatial effect)
stackBHS = inla.stack(data = list(milkZ = data2$milkZ),
                     A = list(A,1),
                     effects = list(c(meshIndexBHS, list(intercept = 1)),
                                    list(cowI = data2$cowI,herd = data2$herd, cowPe = data2$cowPe,
                                         cyrsnI=data2$cyrsnI, tyrmnI=data2$tyrmnI,ward_code= data2$ward_code, dgrpI= data2$dgrpI, ageZ=data2$ageZ, lacgr=data2$lacgr,leg1=data2$leg1,leg2=data2$leg2)), tag = "data2.data")


# ModelBS
formulaBS_inla <- as.formula(paste0(modelB_inla, " + f(fieldID, model = spdeStatBS)-1", collapse = " "))

# fitBHS_inla

formulaBHS_inla <- as.formula(paste0(modelB_inla, " + f(herd, model = 'iid') + f(fieldID, model = spdeStatBHS) -1", collapse = " "))



#----#Alternative Base model without herd effect  but with ward effect included (fitB)-------
# ModelBase without herd

modelB0 <- ~ Intercept(1) +
  fixed_effects(main = ~ cyrsn + tyrmn + dgrp +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv)

modelB0Formula <- milkZ ~ .


# Base model fitB0 + herd  (fitBH)
modelB0H <- ~ Intercept(1) +
  fixed_effects(main = ~ cyrsn + tyrmn + dgrp +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  herd(main=herd, model = 'iid')

modelB0HFormula <- milkZ ~ .

# Base model fitB0 + Spatial effect (fitBS)
modelB0S <- ~ Intercept(1) +
  fixed_effects(main = ~ cyrsn + tyrmn + dgrp +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  field(geometry, model = matern)

modelB0SFormula <- milkZ ~ .

# Base model fitB0 + herd  + Spatial effect  (fitBHS)

modelB0HS <- ~ Intercept(1) +
  fixed_effects(main = ~ cyrsn + tyrmn + dgrp +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  herd(main=herd, model = 'iid') +
  field(geometry, model = matern)

modelB0HSFormula <- milkZ ~ .



# ---- Run the models R-INLABru----------------------------------------------------------
#Fitbase
fitBase <- bru(modelBase,
           like(family = "Gaussian",
                modelBaseFormula,
                data = data2))
#save(fitBase,file = "data/cleaned_data/fitBase/fitBase.RData")

summary(fitBase) #DIC: 32632.72
fitBase$summary.random$fixed_effects
fitBase$summary.random$perm
fitBase$summary.random$animal

# Adding ward_code as a fixed effect (fitWCF)

fitWCF <- bru(modelWCF,
              like(family = "Gaussian",
                   modelWCFFormula,
                   data = data2))
#save(fitWCF,file = "data/cleaned_data/fitWCF/fitWCF.RData")
summary(fitWCF) #DIC: 32606.69
fitWCF$summary.random$fixed_effects
fitWCF$summary.random$perm
fitWCF$summary.random$animal


# Adding ward_code as random independent

fitWCRI <- bru(modelWCRI,
              like(family = "Gaussian",
                   modelWCRIFormula,
                   data = data2))
#save(fitWCRI,file = "data/cleaned_data/fitWCRI/fitWCRI.RData")
summary(fitWCRI) #DIC: 32606.18
fitWCRI$summary.random$fixed_effects
fitWCRI$summary.random$perm
fitWCRI$summary.random$animal
fitWCRI$summary.random$ward_code


# Adding ward_code as random Besag

fitWCRB <- bru(modelWCRB,
               like(family = "Gaussian",
                    modelWCRBFormula,
                    data = data2))

# ERROR *** 	Number of locations and N does not match: 72 != 3644
# I have added mapper=bru_mapper_index(3644) to the model components part
#to solve that.

save(fitWCRB,file = "data/cleaned_data/fitWCRB/fitWCRB.RData")
summary(fitWCRB) #DIC: 32633.57
fitWCRB$summary.random$fixed_effects
fitWCRB$summary.random$perm
fitWCRB$summary.random$animal
fitWCRB$summary.random$ward_code

# Base model fitBase + Spatail effect (fitS)

fitS <- bru(modelfitS,
             like(family = "Gaussian",
                  modelfitSFormula,
                  data = data2))
#save(fitS,file = "data/cleaned_data/fitS/fitS.RData")
summary(fitS) #DIC: 32563.33
fitS$summary.random$fixed_effects
fitS$summary.random$perm
fitS$summary.random$animal

#--------------------Predict spatial effect at my locations----------------------
# Spatial effect in model BS
fieldBS <- predict(fitBS,spdf, ~field)

fieldBS_df <- data.frame(fieldBS[,c(1 , 4,30, 31)])
fieldBS_df <- distinct(fieldBS_df)

# Spatial effect BS
Spatial_BS <- fieldBS_df[,c(1,3)]
names(Spatial_BS)[1]<- "ID"
EBV_BS<- fitBS$summary.random$animal[,c(1,2)]

Spatial_EBV_BS <- merge(EBV_BS,Spatial_BS, by="ID")
cor(Spatial_EBV_BS$mean.x,Spatial_EBV_BS$mean.y) # 0.05997991
plot(Spatial_EBV_BS$mean.x,Spatial_EBV_BS$mean.y)

EBV_B<- fitB$summary.random$animal[,c(1,2)]
Spatial_EBV_BS <- merge(EBV_B,Spatial_BS, by="ID")
cor(Spatial_EBV_BS$mean.x,Spatial_EBV_BOS$mean.y) # 0.5375262
#High correlations between EBV from B0 and B0S
summary(fitBH)

#ggplot() +
  #gg(mesh, color = fieldS$mean)

gps <- data1[,11:12]
vrt2 <- fm_vertices(mesh$g, format = "sf")
fieldS <- predict(fitS, vrt, ~field)

fieldS_df <- data.frame(fieldS[,1:2])
summary(mesh$ge)


# Base model fitBase + ward + Spatial effect (fitWS)

fitWS <- bru(modelWS,
              like(family = "Gaussian",
                   modelWSFormula,
                   data = data2))
save(fitWS,file = "data/cleaned_data/fitWS/fitWS.RData")
summary(fitWS) #DIC:32555.68
fitWS$summary.random$fixed_effects
fitWS$summary.random$perm
fitWS$summary.random$animal


#---------Run additional models-INLAbru--------------------------------------------------
#-------Run the additional models INLA_Bru-factors-----------------------------------
fitB <- bru(modelB,
                like(family = "Gaussian",
                     modelBFormula,
                     data = data2))
save(fitB,file = "data/cleaned_data/fitB/fitB.RData")
summary(fitB) #DIC: 32653.37 (fitB)

fitB$summary.random$fixed_effects
str(fitB$summary.random)
fitB$summary.random$perm
fitB$summary.random$animal

# Base model fitB + herd  (fitBH)
fitBH <- bru(modelBH,
                 like(family = "Gaussian",
                      modelBHFormula,
                      data = data2))
save(fitBH,file = "data/cleaned_data/fitBH/fitBH.RData")
summary(fitBH) #DIC:32606.31 (fitBH)
fitBH$summary.random$fixed_effects
fitBH$summary.random$perm
fitBH$summary.random$animal

# Base model fitB_bru + Spatial effect (fitBS)

fitBS <- bru(modelBS,
                 like(family = "Gaussian",
                      modelBSFormula,
                      data = data2))
save(fitBS,file = "data/cleaned_data/fitBS/fitBS.RData")
summary(fitBS) #DIC:32587.81 (fitBS)
fitBS$summary.random$fixed_effects
fitBS$summary.random$perm
fitBS$summary.random$animal

# Base model fitB + herd  + Spatial effect  (fitBHS)

fitBHS <- bru(modelBHS,
                  like(family = "Gaussian",
                       modelBHSFormula,
                       data = data2))
save(fitBHS,file = "data/cleaned_data/fitBHS/fitBHS.RData")
summary(fitBHS) #DIC: 32556.48 (fitBHS)
fitBHS$summary.random$fixed_effects
fitBHS$summary.random$perm
fitBHS$summary.random$animal



#-------Run the additional models INLA_Bru-numeric_factors------------------------------
# Base model fitB without herd but with ward included

fitB_bru <- bru(modelB_bru,
            like(family = "Gaussian",
                 modelB_bruFormula,
                 data = data2))
save(fitB_bru,file = "data/cleaned_data/fitB_bru/fitB_bru.RData")
summary(fitB_bru) #DIC: 32835.10(fitB_bru) vs 32653.37 (fitB) vs 35579.40 (fitB_inla)

fitB_bru$summary.random$fixed_effects
str(fitB_bru$summary.random)
fitB_bru$summary.random$perm
fitB_bru$summary.random$animal

# Base model fitB_bru + herd  (fitBH_bru)
fitBH_bru <- bru(modelBH_bru,
             like(family = "Gaussian",
                  modelBH_bruFormula,
                  data = data2))
save(fitBH_bru,file = "data/cleaned_data/fitBH_bru/fitBH_bru.RData")
summary(fitBH_bru) #DIC: 32783.17 (fitBH_bru) vs 32606.31 (fitBH) vs 35520.47 (fitBH_inla)
fitBH_bru$summary.random$fixed_effects
fitBH_bru$summary.random$perm
fitBH_bru$summary.random$animal

# Base model fitB_bru + Spatial effect (fitBS_bru)

fitBS_bru <- bru(modelBS_bru,
             like(family = "Gaussian",
                  modelBS_bruFormula,
                  data = data2))
save(fitBS_bru,file = "data/cleaned_data/fitBS_bru/fitBS_bru.RData")
summary(fitBS_bru) #DIC: 32751.03 (fitBS_bru) vs 32587.81 (fitB) vs 35489.22 (fitBS_inla)
fitBS_bru$summary.random$fixed_effects
fitBS_bru$summary.random$perm
fitBS_bru$summary.random$animal

# Base model fitB + herd  + Spatial effect  (fitBHS)

fitBHS_bru <- bru(modelBHS_bru,
              like(family = "Gaussian",
                   modelBHS_bruFormula,
                   data = data2))
save(fitBHS_bru,file = "data/cleaned_data/fitBHS_bru/fitBHS_bru.RData")
summary(fitBHS_bru) #DIC: 32721.15(fitBHS_bru) vs 32556.48 (fitB) vs 35455.30 (fitBHS_inla)
fitBHS_bru$summary.random$fixed_effects
fitBHS_bru$summary.random$perm
fitBHS_bru$summary.random$animal

#---------Run additional models-inla--------------------------------------------
# ----fitB_inla-----------------------------------------------------------------
modelB_inla <- as.formula(modelB_inla)
fitB_inla <- inla(formula = modelB_inla, data = data2,
                control.compute = list(dic = TRUE, config=TRUE))
summary(fitB_inla) # DIC= 35579.40
save(fitB_inla,file = "data/cleaned_data/fitB_inla/fitB_inla.RData")
# ----fitBH_inla----------------------------------------------------------------
fitBH_inla <- inla(formula = modelBH_inla, data = data2,
                  control.compute = list(dic = TRUE, config=TRUE))
summary(fitBH_inla) # DIC=35520.47
save(fitBH_inla,file = "data/cleaned_data/fitBH_inla/fitBH_inla.RData")
# ----fitBS_inla----------------------------------------------------------------
fitBS_inla= inla(formula = formulaBS_inla, data = inla.stack.data(stackBS),
                 family = "normal", control.predictor =list(A=inla.stack.A(stackBS),compute = T),
                 control.family=list(list(hyper=hyperResVarBHS)),
                 control.compute = list(dic=T,cpo=F, config=T),
                 control.fixed = list(expand.factor.strategy="inla"), verbose=T)
summary(fitBS_inla) # DIC= 35489.22
save(fitBS_inla,file = "data/cleaned_data/fitBS_inla/fitBS_inla.RData")

# ----fitBHS_inla---------------------------------------------
fitBHS_inla= inla(formula = formulaBHS_inla, data = inla.stack.data(stackBHS),
                  family = "normal", control.predictor =list(A=inla.stack.A(stackBHS),compute = T),
                  control.family=list(list(hyper=hyperResVarBHS)),
                  control.compute = list(dic=T,cpo=F, config=T),
                  control.fixed = list(expand.factor.strategy="inla"), verbose=T)
summary(fitBHS_inla) # DIC=35455.30
save(fitBHS_inla,file = "data/cleaned_data/fitBHS_inla/fitBHS_inla.RData")











#--------Run models Base without wards and herds--------------------------------

# Base model fitB0 without herd but with ward

fitB0 <- bru(modelB0,
            like(family = "Gaussian",
                 modelB0Formula,
                 data = data2))
save(fitB0,file = "data/cleaned_data/fitB/fitB0.RData")
summary(fitB0) #DIC: 32692.47
fitB0$summary.random$fixed_effects
fitB0$summary.random$perm
fitB0$summary.random$animal

# Base model fitB0 + herd  (fitB0H)
fitB0H <- bru(modelB0H,
             like(family = "Gaussian",
                  modelB0HFormula,
                  data = data2))
#save(fitB0H,file = "data/cleaned_data/fitBH/fitBH.RData")
summary(fitB0H) #DIC: 32632.69
fitB0H$summary.random$fixed_effects
fitB0H$summary.random$perm
fitB0H$summary.random$animal

# Base model fitB0 + Spatail effect (fitB0S)

fitB0S <- bru(modelB0S,
             like(family = "Gaussian",
                  modelB0SFormula,
                  data = data2))
#save(fitB0S,file = "data/cleaned_data/fitB0S/fitB0S.RData")
summary(fitB0S) #DIC:32598.42
fitB0S$summary.random$fixed_effects
fitB0S$summary.random$perm
fitB0S$summary.random$animal


# Base model fitB0 + herd  + Spatial effect  (fitB0HS)

fitB0HS <- bru(modelB0HS,
              like(family = "Gaussian",
                   modelB0HSFormula,
                   data = data2))
#save(fitB0HS,file = "data/cleaned_data/fitB0HS/fitB0HS.RData")
summary(fitB0HS) #DIC:32562.37
fitB0HS$summary.random$fixed_effects
fitB0HS$summary.random$perm
fitB0HS$summary.random$animal






# ---- Summarise functions----------------------------------------------------------
# Create a function to summarise precision to variance

SummarizeFun = function(x, Quantiles = c(0.025, 0.975)) {
  c(mean(x), sd(x), quantile(x, probs = Quantiles))
}

summarise_precision_to_variance = function(x, nSamples = 1000) {
  # Summarize INLA effects "precisions" in form of Standard deviations, Variances, and Proportions (var / sum(all vars))
  Terms = names(x$marginals.hyperpar)
  Terms = Terms[grepl(pattern = "Precision for ", x = Terms)]
  TermsShort = sub(pattern = "Precision for ", replacement = "", x = Terms)
  TermsShort = sub(pattern = "the ",           replacement = "", x = TermsShort)
  nTerms = length(Terms)
  Samples = matrix(data = numeric(), nrow = nSamples, ncol = nTerms + 1)
  Out = vector(mode = "list", length = 3)
  names(Out) = c("Sd", "Var", "Proportion")
  Out[[1]] = Out[[2]] = Out[[3]] = matrix(data = numeric(), nrow = nTerms + 1, ncol = 4)
  dimnames(Out[[1]]) = dimnames(Out[[2]]) = dimnames(Out[[3]]) = list(c(TermsShort, "Total"), c("Mean", "Sd", "Q0.025", "Q0.975"))
  for (Term in 1:nTerms) {
    # Term = 3
    Samples[, Term] = 1 / inla.rmarginal(n = nSamples, marginal = x$marginals.hyperpar[[Terms[Term]]])
  }
  Samples[, Term + 1] = rowSums(x = Samples[, 1:nTerms, drop = FALSE])
  Out$Var[]        = t(apply(X = Samples,                         MARGIN = 2, FUN = SummarizeFun))
  Out$Sd[]         = t(apply(X = sqrt(Samples),                   MARGIN = 2, FUN = SummarizeFun))
  Out$Proportion[] = t(apply(X = Samples / Samples[, nTerms + 1], MARGIN = 2, FUN = SummarizeFun))
  return(Out)
}

SummarizeInlaSpdeVars = function(x, nSamples = 1000, name = 'spatial', spde) {
  # Summarize INLA SPDE hyper-parameters
  SpdeParam = inla.spde2.result(inla = x, name = name, spde = spde)
  Samples = matrix(data = numeric(), nrow = nSamples, ncol = 2)
  Out = vector(mode = "list", length = 3)
  names(Out) = c("Sd", "Var", "Range")
  Out[[1]] = Out[[2]] = Out[[3]] = matrix(data = numeric(), nrow = 1, ncol = 4)
  colnames(Out[[1]]) = colnames(Out[[2]]) = colnames(Out[[3]]) = c("Mean", "Sd", "Q0.025", "Q0.975")
  Samples[, 1] = inla.rmarginal(n = nSamples, marginal = SpdeParam$marginals.variance.nominal[[1]])
  Samples[, 2] = inla.rmarginal(n = nSamples, marginal = SpdeParam$marginals.range.nominal[[1]])
  Out$Sd[]    = SummarizeFun(x = sqrt(Samples[, 1]))
  Out$Var[]   = SummarizeFun(x =      Samples[, 1])
  Out$Range[] = SummarizeFun(x =      Samples[, 2])
  return(Out)
}
SummarizeInlaSpdeVars(fitS) # Did not work. Error in SummarizeInlaSpdeVars(fitS) : object 'spdeStat' not found

# Summarize from inla.posterior.sample
SummarizeINLApostsample = function(x, nSamples = 1000) {
  # Summarize INLA effects "precisions" in form of Standard deviations, Variances, and Proportions (var / sum(all vars))
  Terms = names(x[[nSamples]]$hyperpar)
  Terms = Terms[grepl(pattern = "Precision for ", x = Terms)]
  TermsShort = sub(pattern = "Precision for ", replacement = "", x = Terms)
  TermsShort = sub(pattern = "the ",           replacement = "", x = TermsShort)
  nTerms = length(Terms)
  Samples = matrix(data = numeric(), nrow = nSamples, ncol = nTerms + 1)
  Out = vector(mode = "list", length = 3)
  names(Out) = c("Sd", "Var", "Proportion")
  Out[[1]] = Out[[2]] = Out[[3]] = matrix(data = numeric(), nrow = nTerms + 1, ncol = 4)
  dimnames(Out[[1]]) = dimnames(Out[[2]]) = dimnames(Out[[3]]) = list(c(TermsShort, "Total"), c("Mean", "Sd", "Q0.025", "Q0.975"))
  for (Term in 1:nTerms) {
    # Term = 3
    for (i in 1:nSamples) {
      Samples[i, Term] = 1 / x[[i]]$hyperpar[[Terms[Term]]]
    }
  }
  Samples[, Term + 1] = rowSums(x = Samples[, 1:nTerms, drop = FALSE])
  Out$Var[]        = t(apply(X = Samples,                         MARGIN = 2, FUN = SummarizeFun))
  Out$Sd[]         = t(apply(X = sqrt(Samples),                   MARGIN = 2, FUN = SummarizeFun))
  Out$Proportion[] = t(apply(X = Samples / Samples[, nTerms + 1], MARGIN = 2, FUN = SummarizeFun))
  return(Out)
}



summarise_precision_to_variance(fitWCRI)
#'summary fitWCF'
#summary(fitWCF)

#'Summarize Variances fitWCF'
#summarise_precision_to_variance(fitWCF)
#sink()

# fitWCRI

#sink('results/fitWCRI.txt') #Print outputs of fitWCRI
#'fitWCRI'
#fitWCRI <- inla(formula = modelWCRI, data = data1,
     #          control.compute = list(dic = TRUE, config=TRUE))
#save(fitWCRI,file = "data/cleaned_data/fitWCRI/fitWCRI.RData")
#load(file = "data/cleaned_data/fitWCRI/fitWCRI.RData") # to load the R object
#'Summary fitWCRI'
#summary(fitWCRI)
#'sumarize variance fitWCRI'
#summarise_precision_to_variance(fitWCRI)
#sink()
#summary(fitWCRI)

# fitWCRB

#sink('results/fitWCRB.txt') #Print outputs of fitWCRB
#fitWCRB <- inla(formula = modelWCRB, data = data1,
#               control.compute = list(dic = TRUE, config=TRUE))
#save(fitWCRB,file = "data/cleaned_data/fitWCRB/fitWCRB.RData")
#load(file = file = "data/cleaned_data/fitWCRB/fitWCRB.RData")

#'fitWCRB'
#'Summary fitWCRB'
#summary(fitWCRB)

#'Summarize variance fitWCRB'
#Summarise_precision_to_variance(fitWCRB)
#sink() # close sink fitWCRB





# ModelBase + herd


fit_herd <- bru(modelBaseherd,
           like(family = "Gaussian",
                modelBaseherdFormula,
                data = data2))

summary(fit_herd) #DIC: 32604.35
fit_herd$summary.random$fixed_effects
fit_herd$summary.random$perm
fit_herd$summary.random$animal
fit_herd$summary.random$herd


#Testing SPDE
#modelBaseS_test <- ~ Intercept(1) + fixed_effects(main = ~ 1 + cyrsn, model = "fixed") +
#  field(geometry, model = matern)

#modelBaseS_testFormula <- milkZ ~ .

#fitS_test <- bru(modelBaseS_test,
    #        like(family = "Gaussian",
  #               modelBaseS_testFormula,
  #               data = data2))
#summary(fitS_test)

# Real Model: fitS (Adding spatial effect to the base model without ward)
modelBaseS <- ~ Intercept(1) +
  fixed_effects(main = ~ cyrsn + tyrmn + dgrp + lacgr +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + herd(main=herd, model="iid")
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  field(geometry, model = matern)

modelBaseSFormula <- milkZ ~ .

fitS <- bru(modelBaseS,
             like(family = "Gaussian",
                  modelBaseSFormula,
                  data = data2))

summary(fitS) #DIC: 32597.74
fitS$summary.random$fixed_effects
fitS$summary.random$perm
fitS$summary.random$animal
fitS$summary.random$herd
fitS$summary.random$field


# Model fitWS (Adding ward and spatial effect to the base model with )
modelBaseWS <- ~ Intercept(1) +
  fixed_effects(main = ~ cyrsn + tyrmn + dgrp + lacgr +
                  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  perm(main = cowPe, model = 'iid') + ward(main = ward_code, model = 'iid') +
  herd(main=herd, model="iid") +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  field(geometry, model = matern)


modelBaseWSFormula <- milkZ ~ .

fitWS <- bru(modelBaseWS,
           like(family = "Gaussian",
                modelBaseWSFormula,
                data = data2))

summary(fitWS) #DIC: 32585.44
fitS$summary.random$fixed_effects
fitS$summary.random$perm
fitS$summary.random$animal
fitS$summary.random$herd
fitS$summary.random$field

#----------Plot posterior distributions of hyperparameters for all models--------

# Genetic variance
# * fitBase
ggplot(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                fitBase$marginals.hyperpar$`Precision for cowI`))), aes(x, y)) +
  geom_line() +
  theme_bw()
# * fitWCF
ggplot(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                fitWCF$marginals.hyperpar$`Precision for cowI`))), aes(x, y)) +
  geom_line() +
  theme_bw()
# * fitWCRI
ggplot(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                fitWCRI$marginals.hyperpar$`Precision for cowI`))), aes(x, y)) +
  geom_line() +
  theme_bw()
# * fitWCRB
ggplot(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                fitWCRB$marginals.hyperpar$`Precision for cowI`))), aes(x, y)) +
  geom_line() +
  theme_bw()
# * fitS
ggplot(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                fitS$marginals.hyperpar$`Precision for cowI`))), aes(x, y)) +
  geom_line() +
  theme_bw()
# * fitWS
ggplot(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                fitWS$marginals.hyperpar$`Precision for cowI`))), aes(x, y)) +
  geom_line() +
  theme_bw()
#Plot posterior distribution of all models in a single graph
# Plot genetic variance across models

Pg<-ggplot()  +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBase$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="B")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCF$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="BWF"),linetype="dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRI$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="BWRI")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRB$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="BWRB")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitS$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="BS"), linetype = "dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWS$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="BWRIS"), linetype = "dashed") +
  scale_color_manual(values = c("B" = "black",
                                "BWF" = "yellow",
                                "BWRI"="red",
                                "BWRB"="orange",
                                "BS"="blue",
                                "BWRIS"="green")) +
  labs(x="", y="", title ="Genetic variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model"))
# scale_y_continuous(breaks=NULL) to remove axis labels
Pg


# Residual variance

Pr <-ggplot()  +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBase$marginals.hyperpar$`Precision for the Gaussian observations`))), mapping = aes(x, y, colour="B")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCF$marginals.hyperpar$`Precision for the Gaussian observations`))), mapping = aes(x, y, colour="BWF"),linetype="dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRI$marginals.hyperpar$`Precision for the Gaussian observations`))), mapping = aes(x, y, colour="BWRI")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRB$marginals.hyperpar$`Precision for the Gaussian observations`))), mapping = aes(x, y, colour="BWRB")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitS$marginals.hyperpar$`Precision for the Gaussian observations`))), mapping = aes(x, y, colour="BS"), linetype = "dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWS$marginals.hyperpar$`Precision for the Gaussian observations`))), mapping = aes(x, y, colour="BWRIS"), linetype = "dashed") +
  scale_color_manual(values = c("B" = "black",
                                "BWF" = "yellow",
                                "BWRI"="red",
                                "BWRB"="orange",
                                "BS"="blue",
                                "BWRIS"="green")) +
  labs(x="", y="", title ="Residual variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model"))
# scale_y_continuous(breaks=NULL) to remove axis labels
Pr


# Herd_Variance
Ph<-ggplot()  +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBase$marginals.hyperpar$`Precision for herdI`))), mapping = aes(x, y, colour="B")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCF$marginals.hyperpar$`Precision for herdI`))), mapping = aes(x, y, colour="BWF"),linetype="dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRI$marginals.hyperpar$`Precision for herdI`))), mapping = aes(x, y, colour="BWRI")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRB$marginals.hyperpar$`Precision for herdI`))), mapping = aes(x, y, colour="BWRB")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitS$marginals.hyperpar$`Precision for herdI`))), mapping = aes(x, y, colour="BS"), linetype = "dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWS$marginals.hyperpar$`Precision for herdI`))), mapping = aes(x, y, colour="BWRIS"), linetype = "dashed") +
  scale_color_manual(values = c("B" = "black",
                                "BWF" = "yellow",
                                "BWRI"="red",
                                "BWRB"="orange",
                                "BS"="blue",
                                "BWRIS"="green")) +
  labs(x="", y="", title ="Herd variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model"))
# scale_y_continuous(breaks=NULL) to remove axis labels
Ph


# Ward_variance (fitWCRI, fitWCRB and fitWS)
Pw<-ggplot()  +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRI$marginals.hyperpar$`Precision for ward_codeI`))), mapping = aes(x, y, colour="BWRI")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRB$marginals.hyperpar$`Precision for ward_codeI`))), mapping = aes(x, y, colour="BWRB")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWS$marginals.hyperpar$`Precision for ward_codeI`))), mapping = aes(x, y, colour="BWRIS"), linetype = "dashed") +
  scale_color_manual(values = c("BWRI"="red",
                                "BWRB"="orange",
                                "BWRIS"="green")) +
  labs(x="", y="", title ="Ward variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model"))
# scale_y_continuous(breaks=NULL) to remove axis labels
Pw

# Spatial variance

#fitS
#Compute statistics in terms of range and variance
spde.est_fitS <- inla.spde2.result(inla = fitS, name = "fieldID",
                              spde = spdeStatS, do.transf = TRUE)

# Spatial Variance
#inla.zmarginal(spde.est_fitS$marginals.variance.nominal[[1]])

Sv<- ggplot()  +

  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) (exp(x))^2,
                                                     fitS$internal.marginals.hyperpar[[6]]))), mapping = aes(x, y, colour="BS"), linetype = "dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) (exp(x))^2,
                                                     fitWS$internal.marginals.hyperpar[[6]]))), mapping = aes(x, y, colour="BWRIS"), linetype = "dashed") +
  scale_color_manual(values = c("BS"="blue",
                                "BWRIS"="green")) +
  labs(x="", y="", title ="Spatial variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model"))
Sv




# Spatial range

Sr<- ggplot()  +

  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) (exp(x))^2,
                                                     fitS$internal.marginals.hyperpar[[5]]))), mapping = aes(x, y, colour="BS"), linetype = "dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) (exp(x))^2,
                                                     fitWS$internal.marginals.hyperpar[[5]]))), mapping = aes(x, y, colour="BWRIS"), linetype = "dashed") +
  scale_color_manual(values = c("BS"="blue",
                                "BWRIS"="green")) +
  labs(x="", y="", title ="Spatial range") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model"))
Sr

# Variance for spatial component:
teta1 = fitS$internal.marginals.hyperpar[[6]]
# Theta 1 = log(SD)
# To get VAR:
te1 = inla.tmarginal(function(x) (exp(x))^2, teta1)
Samples = matrix(data = numeric(), nrow = 1000, ncol = 1)
Samples[, 1] = inla.rmarginal(n = 1000, te1)
mean(Samples) #################### 17.36


# Range for spatial component:
teta2 = fitS$internal.marginals.hyperpar[[5]]
# Theta 2 = log(range)
# To get range:
te2 = inla.tmarginal(function(x) (exp(x)), teta2)

Samples = matrix(data = numeric(), nrow = 1000, ncol = 1)
Samples[, 1] = inla.rmarginal(n = 1000, te2)
mean(Samples) #################### 3.482276


#fitWS
#Compute statistics in terms of range and variance
spde.est_fitWS <- inla.spde2.result(inla = fitWS, name = "fieldID",
                                   spde = spdeStatWS, do.transf = TRUE)

# Variance for spatial component:
teta1 = fitWS$internal.marginals.hyperpar[[6]]
# Theta 1 = log(SD)
# To get VAR:
te1 = inla.tmarginal(function(x) (exp(x))^2, teta1)
Samples = matrix(data = numeric(), nrow = 1000, ncol = 1)
Samples[, 1] = inla.rmarginal(n = 1000, te1)
mean(Samples) #################### 19.6


# Range for spatial component:
teta2 = fitWS$internal.marginals.hyperpar[[5]]
# Theta 2 = log(range)
# To get range:
te2 = inla.tmarginal(function(x) (exp(x)), teta2)
Samples = matrix(data = numeric(), nrow = 1000, ncol = 1)
Samples[, 1] = inla.rmarginal(n = 1000, te2)
mean(Samples) # 38.9

#---------Combined Plot Variance components------------------------------------------------
(p <- ggarrange(Pg,Pr,Ph,Pw,Sv,Sr, ncol = 2, nrow = 3, common.legend = T, legend = "left", align = "h", widths = c(1,1,1)))

#-------------------Manuscript and Presentation_Plot----------------------------
# Plot genetic variance across models

PgE<-ggplot()  +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitB$marginals.hyperpar$`Precision for animal`))), mapping = aes(x, y, colour="B")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBH$marginals.hyperpar$`Precision for animal`))), mapping = aes(x, y, colour="BH")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBS$marginals.hyperpar$`Precision for animal`))), mapping = aes(x, y, colour="BS"), linetype = "dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBHS$marginals.hyperpar$`Precision for animal`))), mapping = aes(x, y, colour="BHS"), linetype = "dashed") +
  scale_color_manual(values = c("B" = "black",
                                "BH"="red",
                                "BS"="blue",
                                "BHS"="green")) +
  labs(x="", y="", title ="Genetic variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model"))
# scale_y_continuous(breaks=NULL) to remove axis labels
PgE<- PgE+ theme(plot.title = element_text(face = "bold"))
PgE


# Plot genetic variance
genvarB <- plot (fitB, "animal")
genvarB





# Plot spatial variance
varianceBS <- spde.posterior(fitBS, "field", "variance")
plot(varianceBS)
varianceBHS <- spde.posterior(fitBHS, "field", "variance")
plot(varianceBHS)

# Plot spatial range
rangeBS <- spde.posterior(fitBS, "field", "range")
plot(rangeBS)
rangeBHS <- spde.posterior(fitBHS, "field", "range")
plot(rangeBHS)


# Residual variance

PrE <-ggplot()  +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBase$marginals.hyperpar$`Precision for the Gaussian observations`))), mapping = aes(x, y, colour="B")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRI$marginals.hyperpar$`Precision for the Gaussian observations`))), mapping = aes(x, y, colour="BW")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitS$marginals.hyperpar$`Precision for the Gaussian observations`))), mapping = aes(x, y, colour="BS"), linetype = "dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWS$marginals.hyperpar$`Precision for the Gaussian observations`))), mapping = aes(x, y, colour="BWS"), linetype = "dashed") +
  scale_color_manual(values = c("B" = "black",
                                "BW"="red",
                                "BS"="blue",
                                "BWS"="green")) +
  labs(x="", y="", title ="Residual variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model"))
# scale_y_continuous(breaks=NULL) to remove axis labels
PrE<- PrE+ theme(plot.title = element_text(face = "bold"))


# Herd_Variance
PhE<-ggplot()  +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBase$marginals.hyperpar$`Precision for herdI`))), mapping = aes(x, y, colour="B")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRI$marginals.hyperpar$`Precision for herdI`))), mapping = aes(x, y, colour="BW")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitS$marginals.hyperpar$`Precision for herdI`))), mapping = aes(x, y, colour="BS"), linetype = "dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWS$marginals.hyperpar$`Precision for herdI`))), mapping = aes(x, y, colour="BWS"), linetype = "dashed") +
  scale_color_manual(values = c("B" = "black",
                                "BW"="red",
                                "BS"="blue",
                                "BWS"="green")) +
  labs(x="", y="", title ="Herd variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model"))
# scale_y_continuous(breaks=NULL) to remove axis labels
PhE<- PhE+ theme(plot.title = element_text(face = "bold"))


# Ward_variance (fitWCRI, fitWCRB and fitWS)
PwE<-ggplot()  +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRI$marginals.hyperpar$`Precision for ward_codeI`))), mapping = aes(x, y, colour="BW")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWS$marginals.hyperpar$`Precision for ward_codeI`))), mapping = aes(x, y, colour="BWS"), linetype = "dashed") +
  scale_color_manual(values = c("BW"="red",
                                "BWS"="green")) +
  labs(x="", y="", title ="Ward variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model"))
# scale_y_continuous(breaks=NULL) to remove axis labels
PwE <- PwE + theme(plot.title = element_text(face = "bold"))

# Permanent environmental variance

PeE<-ggplot()  +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBase$marginals.hyperpar$`Precision for cowPeI`))), mapping = aes(x, y, colour="B")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRI$marginals.hyperpar$`Precision for cowPeI`))), mapping = aes(x, y, colour="BW")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitS$marginals.hyperpar$`Precision for cowPeI`))), mapping = aes(x, y, colour="BS"), linetype = "dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWS$marginals.hyperpar$`Precision for cowPeI`))), mapping = aes(x, y, colour="BWS"), linetype = "dashed") +
  scale_color_manual(values = c("B" = "black",
                                "BW"="red",
                                "BS"="blue",
                                "BWS"="green")) +
  labs(x="", y="", title ="Permanent environment variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model"))
# scale_y_continuous(breaks=NULL) to remove axis labels
PeE<- PeE + theme(plot.title = element_text(face = "bold"))

PeE
  ggsave(plot = PeE + PreseTheme, filename = "Permanent_environmental_variance_presentation.png",
         height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation



# Spatial variance

#fitS


SvE<- ggplot()  +

  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) (exp(x))^2,
                                                     fitS$internal.marginals.hyperpar[[6]]))), mapping = aes(x, y, colour="BS"), linetype = "dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) (exp(x))^2,
                                                     fitWS$internal.marginals.hyperpar[[6]]))), mapping = aes(x, y, colour="BWS"), linetype = "dashed") +
  scale_color_manual(values = c("BS"="blue",
                                "BWS"="green")) +
  labs(x="", y="", title ="Spatial variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model"))
SvE <- SvE + theme(plot.title = element_text(face = "bold"))



# Spatial range

SrE<- ggplot()  +

  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) (exp(x))^2,
                                                     fitS$internal.marginals.hyperpar[[5]]))), mapping = aes(x, y, colour="BS"), linetype = "dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) (exp(x))^2,
                                                     fitWS$internal.marginals.hyperpar[[5]]))), mapping = aes(x, y, colour="BWS"), linetype = "dashed") +
  scale_color_manual(values = c("BS"="blue",
                                "BWS"="green")) +
  labs(x="", y="", title ="Spatial range") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Model")) +
  guides(colour=guide_legend(title="Model"))
SrE<- SrE + theme(plot.title = element_text(face = "bold")) # Making title Bold
SrE<- SrE + geom_line(size=3)
SrE
# Combined plot for EAAP
PaperTheme = theme_bw(base_size = 11, base_family = "serif") +
  theme(legend.position = "top")

PaperThemeLegendright = theme_bw(base_size = 11, base_family = "serif") +
  theme(legend.position = "right")

PaperThemeNoLegend = PaperTheme + theme(legend.position = "none")
PaperSize = 10
PreseTheme = theme_bw(base_size = 18, base_family = "sans") +
  theme(legend.position = "top")
PreseThemeLegendright = theme_bw(base_size = 18, base_family = "sans") +
  theme(legend.position = "right")
PreseThemeNoLegend = PreseTheme + theme(legend.position = "none")
PreseSize = 16

(pE <- ggarrange(PgE,PrE,PhE,PwE,SvE,SrE, ncol = 2, nrow = 3, common.legend = T, legend = "left", align = "h", widths = c(1,1,1))) + PreseTheme
pE

ggsave(plot = pE + PreseTheme, filename = "Posterior_distribution_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = pE + PaperTheme, filename = "Posterior_distribution_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper



# ---- Plot posterior mean (a) and standard deviation of the estimated spatial effect ----
#fitBHS

gproj <- inla.mesh.projector(mesh,  dims = c(300, 300))
g.mean <- inla.mesh.project(gproj, fitBHS$summary.random$field$mean)
g.sd <- inla.mesh.project(gproj, fitBHS$summary.random$field$sd)

plot_BHS<- grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='Easting', ylab='Northing', main='a',col.regions = viridis(16)),
             levelplot(g.sd, scal=list(draw=F), xla='Easting', yla='Northing', main='b' ,col.regions = viridis(16)), nrow=2) + PaperTheme







# a: posterior mean
#b: uncertainty (posterior standard deviation)

#FitBS
gproj <- inla.mesh.projector(mesh,  dims = c(300, 300))
g.mean <- inla.mesh.project(gproj, fitBS$summary.random$field$mean)
g.sd <- inla.mesh.project(gproj, fitBS$summary.random$field$sd)

plot_BS<- grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='Easting', ylab='Northing', main='a',col.regions = viridis(16)),
             levelplot(g.sd, scal=list(draw=F), xla='Easting', yla='Northing', main='b' ,col.regions = viridis(16)), nrow=2) + PaperTheme

# a: posterior mean
#b: uncertainty (posterior standard deviation)
#Title: "Posterior mean (a) and standard deviation (b) of the estimated spatial effect (in units of posterior spatial standard deviation) from model BHS .png"

# Saving plot for paper





#------------Boxplot-----------------------------------------------------------

# Correlation between EBVs in models with spatial effects and in models without spatial effect


EBV_Spde <- cbind(EBV_B,EBV_S, EBV_BW, EBV_BWS, SPDE_cow)

EBV_Spde<- subset(EBV_Spde, select = c(ID, EBV_B,EBV_S, EBV_BW, EBV_BWS, spdeS_mean, spdeWS_mean) )
EBV_Spde$dBS <-EBV_Spde$EBV_B-EBV_Spde$EBV_S
round(summary(EBV_Spde$spdeS_mean),2)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-1.59   -0.63   -0.16   -0.04    0.48    2.70
EBV_Spde <- EBV_Spde %>%
  mutate(group=
           case_when((EBV_Spde$spdeS_mean > -1.59 & EBV_Spde$spdeS_mean <= -0.63) ~"(-1.59,-0.63]",
                     (EBV_Spde$spdeS_mean > -0.63 & EBV_Spde$spdeS_mean<= -0.16) ~  "(-0.63,-0.16]",
                     (EBV_Spde$spdeS_mean > -0.16 & EBV_Spde$spdeS_mean<= -0.04) ~ "(-0.16, -0.04]",
                     (EBV_Spde$spdeS_mean> -0.04  & EBV_Spde$spdeS_mean<= 0.48) ~ "(-0.04,0.48]",
                     (EBV_Spde$spdeS_mean> 0.48  & EBV_Spde$spdeS_mean<= 2.70) ~  "(0.48,2.70]"))


# Set the default theme to theme_classic() with the legend at the right of the plot:
theme_set(
  theme_classic() +
    theme(legend.position = "right")
)
# Change the default order of items
EBV_Spde$group <- as.factor(EBV_Spde$group)
e <- ggplot(EBV_Spde, aes(x =group , y = dBS)) +
  geom_boxplot()
e
e<- e + geom_boxplot() +
  scale_x_discrete(limits=c("(-1.59,-0.63]","(-0.63,-0.16]","(-0.16, -0.04]","(-0.04,0.48]","(0.48,2.70]")) # Reorder x axis items.
e <- e + geom_boxplot(lwd=0.2, fatten=5) # make line width thicker

# Now let's make the median line less thicker
e
# Changing axis names
e<- e + labs(x = "Spatial effect", y = "Difference in breeding value (BH-BHS)")
e<- e+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=20,face="bold"))



ggsave(plot = e + PreseTheme, filename = "Boxplot_difference_EBVs_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm")


ggsave(plot = e + PaperTheme, filename = "Boxplot_difference_EBVs_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm")


summary(sd_EBV_Spde$dBS)
sd_EBV_Spde <- sd_EBV_Spde %>%
  mutate(group=
           case_when((sd_EBV_Spde$sd_spdeS_mean > 0.09 & sd_EBV_Spde$sd_spdeS_mean <= 0.15) ~"(0.09,0.15]",
                     (sd_EBV_Spde$sd_spdeS_mean > 0.15 & sd_EBV_Spde$sd_spdeS_mean<= 0.18) ~  "(0.15,0.18]",
                     (sd_EBV_Spde$sd_spdeS_mean > 0.18 & sd_EBV_Spde$sd_spdeS_mean<= 0.19) ~ "(0.18, 0.19]",
                     (sd_EBV_Spde$sd_spdeS_mean> 0.19  & sd_EBV_Spde$sd_spdeS_mean<= 0.22) ~ "(0.19,0.22]",
                     (sd_EBV_Spde$sd_spdeS_mean> 0.22  & sd_EBV_Spde$sd_spdeS_mean<= 0.29) ~  "(0.22,0.29]"))


# Set the default theme to theme_classic() with the legend at the right of the plot:
theme_set(
  theme_classic() +
    theme(legend.position = "right")
)
# Change the default order of items
sd_EBV_Spde$group <- as.factor(sd_EBV_Spde$group)
e <- ggplot(sd_EBV_Spde, aes(x =group , y = dBS)) +
  geom_boxplot()
e
e<- e + geom_boxplot() +
  scale_x_discrete(limits=c("(0.09,0.15]", "(0.15,0.18]", "(0.18, 0.19]", "(0.19,0.22]","(0.22,0.29]" )) # Reorder x axis items.
e

# Changing axis names
e<- e + labs(x = "Spatial effect", y = "Difference in breeding value")
e<- e+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=20,face="bold"))

#-------------------Accuracy of Prediction: Cross validation--------------------

#----------- Accuracy model fitBase----------------------------------------------

# Masking phenotypes of cows in dgrp 1 (Breed proportion class 1)
#dgrp1 masked
data1_1_NA<- data1 %>% mutate(milkZ = ifelse(dgrp == "1", NA, milkZ))
sum(is.na(data1$milkZ)) # 0
sum(is.na(data1_1_NA$milkZ)) #7466
length(unique(data1_1_NA$cowI)) # 1894
length(unique(data1$cowI)) # 1894
fitBase_NA1 <- inla(formula = modelfitBase, data = data1_1_NA,
               control.compute = list(dic = TRUE,config=TRUE))

# save fitBase_NA1 as R object
save(fitBase_NA1,file = "data/cleaned_data/fitBase_pred/fitBase_NA1.RData")
 #load(file = "data/cleaned_data/FitBase_pred/FitBase_NA1.RData") # to load the R object
pheno_pred1 <- fitBase_NA1$summary.linear.predictor
colnames(pheno_pred1)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred1$mean
data1$milkZ_pred_sd <- pheno_pred1$sd
pheno1 <-   subset(data1, dgrp==1)
pheno1<- subset(pheno1,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno1$cowI))# 731 Cows in breed composition class 1

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotyped
accuracy_fitBase_1 <- round(cor(pheno1$milkZ,pheno1$milkZ_pred),2)
Coef1<- lm (pheno1$milkZ~pheno1$milkZ_pred)
accuracy_fitBase_1 #  0.24
R2_fitBase_1<- summary(Coef1)
R2_fitBase_1 #  0.05932
#CRPS
obs <- pheno1$milkZ
pred<- subset(pheno1,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBase_1 <- crps(obs,pred)
round((crps_fitBase_1$CRPS),2) # 0.6

# Masking phenotypes of cows in dgrp 2 (Breed proportion class 2)
#dgrp2 masked

data1_2_NA<- data1 %>% mutate(milkZ = ifelse(dgrp == "2", NA, milkZ))
sum(is.na(data1_2_NA$milkZ)) # 8149
length(unique(data1_2_NA$cowI)) # 1894

fitBase_NA2 <- inla(formula = modelfitBase, data = data1_2_NA,
                   control.compute = list(dic = TRUE,config=TRUE))
save(fitBase_NA2,file = "data/cleaned_data/fitBase_pred/fitBase_NA2.RData")
#load(file = "data/cleaned_data/fitBase_pred/fitBase_NA2.RData") # to load the R object

pheno_pred2 <- fitBase_NA2$summary.linear.predictor
colnames(pheno_pred2)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred2$mean
data1$milkZ_pred_sd <- pheno_pred2$sd

pheno2 <-   subset(data1, dgrp==2)
pheno2<- subset(pheno2,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno2$cowI))#  770 Cows in breed composition class 2

# Correlating observed phenotype of cows in dgrp 2 with predicted phenotypes
accuracy_fitBase_2 <- round(cor(pheno2$milkZ,pheno2$milkZ_pred),2)
Coef2<- lm (pheno2$milkZ~pheno2$milkZ_pred)
accuracy_fitBase_2# 0.34
R2_fitBase_2<- summary(Coef2)
R2_fitBase_2 # 0.1123
#CRPS
obs <- pheno2$milkZ
pred<- subset(pheno2,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBase_2 <- crps(obs,pred)
round((crps_fitBase_2$CRPS),2) # 0.53

# Masking phenotypes of cows in dgrp 3 (Breed proportion class 3)
#dgrp3 masked
data1_3_NA<- data1 %>% mutate(milkZ = ifelse(dgrp == "3", NA, milkZ))
sum(is.na(data1_3_NA$milkZ)) # 3058 records
length(unique(data1_3_NA$cowI)) # 1894 cows in data1

fitBase_NA3 <- inla(formula = modelfitBase, data = data1_3_NA,
                   control.compute = list(dic = TRUE,config=TRUE))

save(fitBase_NA3,file = "data/cleaned_data/fitBase_pred/fitBase_NA3.RData")
 #load(file = "data/cleaned_data/fitBase_pred/fitBase_NA3.RData") # to load the R object

pheno_pred3 <- fitBase_NA3$summary.linear.predictor
sum(is.na(pheno_pred3$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred3$mean
data1$milkZ_pred_sd <- pheno_pred3$sd

pheno3 <-   subset(data1, dgrp==3)
pheno3<- subset(pheno3,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno3$cowI))#  309 Cows in breed composition class 3

# Correlating observed phenotype of cows in dgrp 3 with predicted phenotypes
accuracy_fitBase_3 <- round(cor(pheno3$milkZ,pheno3$milkZ_pred),2)
Coef3<- lm (pheno3$milkZ~pheno3$milkZ_pred)
accuracy_fitBase_3# 0.28
R2_fitBase_3<- summary(Coef3)
R2_fitBase_3 #0.08083
#CRPS
obs <- pheno3$milkZ
pred<- subset(pheno3,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBase_3 <- crps(obs,pred)
round((crps_fitBase_3$CRPS),2) # 0.52


# Masking phenotypes of cows in dgrp 4 (Breed proportion class 4)
#dgrp4 masked
data1_4_NA<- data1 %>% mutate(milkZ = ifelse(dgrp == "4", NA, milkZ))
sum(is.na(data1_4_NA$milkZ)) # 702 records
length(unique(data1_4_NA$cowI)) # 1894

fitBase_NA4 <- inla(formula = modelfitBase, data = data1_4_NA,
                   control.compute = list(dic = TRUE,config=TRUE))
 save(fitBase_NA4,file = "data/cleaned_data/fitBase_pred/fitBase_NA4.RData")
 #load(file = "data/cleaned_data/fitBase_pred/fitBase_NA4.RData") # to load the R object

pheno_pred4 <- fitBase_NA4$summary.linear.predictor
sum(is.na(pheno_pred4$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred4$mean
data1$milkZ_pred_sd <- pheno_pred4$sd

pheno4 <-   subset(data1, dgrp==4)
pheno4<- subset(pheno4,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno4$cowI))#  84 Cows in breed composition class 4

# Correlating observed phenotype of cows in dgrp 4 with predicted phenotypes
accuracy_fitBase_4 <- round(cor(pheno4$milkZ,pheno4$milkZ_pred),2)
Coef4<- lm (pheno4$milkZ~pheno4$milkZ_pred)
accuracy_fitBase_4 # -0.05
R2_fitBase_4<- summary(Coef4)
R2_fitBase_4 # 0.0008917
#CRPS
obs <- pheno4$milkZ
pred<- subset(pheno4,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBase_4 <- crps(obs,pred)
round((crps_fitBase_4$CRPS),2) # 0.58


# Accuracy of prediction and degree under/overprediction of model FitBase

accuracy_fitBase = (accuracy_fitBase_1 + accuracy_fitBase_2 + accuracy_fitBase_3 + accuracy_fitBase_4)/4
round((accuracy_fitBase),2) # 0.2

R2_fitBase = (0.05932+0.1123+0.08083+0.0008917)/4
round((R2_fitBase),2) # 0.06

# crps_fitBase= (crps_fitBase_1+crps_fitBase_2+crps_fitBase_3+crps_fitBase_4)/4
crps_fitBase = (0.6+ 0.53 + 0.52 + 0.58)/4
crps_fitBase=round((crps_fitBase),2)
crps_fitBase #0.56

# Another method to calculate Accuracy, CRPS and R2 (degree of over/underprediction)
# Let's combine the predicted and observed phenotypes
pheno_fitBase <- rbind(pheno1,pheno2,pheno3,pheno4)
accuracy_fitBase<- round(cor(pheno_fitBase$milkZ,pheno_fitBase$milkZ_pred),2)
Coef_fitBase <- lm (pheno_fitBase$milkZ~pheno_fitBase$milkZ_pred)
accuracy_fitBase # 0.24
R2_fitBase<- summary(Coef_fitBase)
R2_fitBase# 0.05613
#CRPS_fitBase
obs <- pheno_fitBase$milkZ
pred<- subset(pheno_fitBase,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBase <- crps(obs,pred)
round((crps_fitBase$CRPS),2) # 0.55

#-------------------------------------------------------------------------------
#----------- Accuracy model FitWCF----------------------------------------------

# Masking phenotypes of cows in dgrp 1 (Breed proportion class 1)
#dgrp1 masked

fitWCF_NA1 <- inla(formula = modeltWCF, data = data1_1_NA,
                    control.compute = list(dic = TRUE,config=TRUE))
# save fitWCF_NA1 as R object
save(fitWCF_NA1,file = "data/cleaned_data/fitWCF_pred/fitWCF_NA1.RData")
#load(file = "data/cleaned_data/fitWCF_pred/fitWCF_NA1.RData") # to load the R object

pheno_pred1 <- fitWCF_NA1$summary.linear.predictor
colnames(pheno_pred1)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred1$mean
data1$milkZ_pred_sd <- pheno_pred1$sd
pheno1 <-   subset(data1, dgrp==1)
pheno1<- subset(pheno1,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno1$cowI))# 731 Cows in breed composition class 1

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotyped
accuracy_fitWCF_1 <- round(cor(pheno1$milkZ,pheno1$milkZ_pred),2)
Coef1<- lm (pheno1$milkZ~pheno1$milkZ_pred)
accuracy_fitWCF_1 #  0.24
R2_fitWCF_1<- summary(Coef1)
R2_fitWCF_1 #  0.0593
#CRPS
obs <- pheno1$milkZ
pred<- subset(pheno1,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCF_1 <- crps(obs,pred)
round((crps_fitWCF_1$CRPS),2) # 0.6

# Masking phenotypes of cows in dgrp 2 (Breed proportion class 2)
#dgrp2 masked

fitWCF_NA2 <- inla(formula = modelWCF, data = data1_2_NA,
                    control.compute = list(dic = TRUE,config=TRUE))
save(fitWCF_NA2,file = "data/cleaned_data/fitWCF_pred/fitWCF_NA2.RData")
#load(file = "data/cleaned_data/fitWCF_pred/fitWCF_NA2.RData") # to load the R object

pheno_pred2 <- fitWCF_NA2$summary.linear.predictor
colnames(pheno_pred2)
sum(is.na(pheno_pred2$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred2$mean
data1$milkZ_pred_sd <- pheno_pred2$sd

pheno2 <-   subset(data1, dgrp==2)
pheno2<- subset(pheno2,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno2$cowI))#  770 Cows in breed composition class 2

# Correlating observed phenotype of cows in dgrp 2 with predicted phenotypes
accuracy_fitWCF_2 <- round(cor(pheno2$milkZ,pheno2$milkZ_pred),2)
Coef2<- lm (pheno2$milkZ~pheno2$milkZ_pred)
accuracy_fitWCF_2# 0.34
R2_fitWCF_2<- summary(Coef2)
R2_fitWCF_2 # 0.1124
#CRPS
obs <- pheno2$milkZ
pred<- subset(pheno2,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCF_2 <- crps(obs,pred)
round((crps_fitWCF_2$CRPS),2) # 0.53

# Masking phenotypes of cows in dgrp 3 (Breed proportion class 3)
#dgrp3 masked

fitWCF_NA3 <- inla(formula = modelWCF, data = data1_3_NA,
                    control.compute = list(dic = TRUE,config=TRUE))

save(fitWCF_NA3,file = "data/cleaned_data/fitWCF_pred/fitWCF_NA3.RData")
#load(file = "data/cleaned_data/fitWCF_pred/fitWCF_NA3.RData") # to load the R object

pheno_pred3 <- fitWCF_NA3$summary.linear.predictor
sum(is.na(pheno_pred3$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred3$mean
data1$milkZ_pred_sd <- pheno_pred3$sd

pheno3 <-   subset(data1, dgrp==3)
pheno3<- subset(pheno3,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno3$cowI))#  309 Cows in breed composition class 3

# Correlating observed phenotype of cows in dgrp 3 with predicted phenotypes
accuracy_fitWCF_3 <- round(cor(pheno3$milkZ,pheno3$milkZ_pred),2)
Coef3<- lm (pheno3$milkZ~pheno3$milkZ_pred)
accuracy_fitWCF_3# 0.28
R2_fitWCF_3<- summary(Coef3)
R2_fitWCF_3 #0.08083
#CRPS
obs <- pheno3$milkZ
pred<- subset(pheno3,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCF_3 <- crps(obs,pred)
round((crps_fitWCF_3$CRPS),2) # 0.52


# Masking phenotypes of cows in dgrp 4 (Breed proportion class 4)
#dgrp4 masked

fitWCF_NA4 <- inla(formula = modelWCF, data = data1_4_NA,
                    control.compute = list(dic = TRUE,config=TRUE))
save(fitWCF_NA4,file = "data/cleaned_data/fitWCF_pred/fitWCF_NA4.RData")
#load(file = "data/cleaned_data/fitWCF_pred/fitWCF_NA4.RData") # to load the R object

pheno_pred4 <- fitWCF_NA4$summary.linear.predictor
sum(is.na(pheno_pred4$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred4$mean
data1$milkZ_pred_sd <- pheno_pred4$sd

pheno4 <-   subset(data1, dgrp==4)
pheno4<- subset(pheno4,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno4$cowI))#  84 Cows in breed composition class 4

# Correlating observed phenotype of cows in dgrp 4 with predicted phenotypes
accuracy_fitWCF_4 <- round(cor(pheno4$milkZ,pheno4$milkZ_pred),2)
Coef4<- lm (pheno4$milkZ~pheno4$milkZ_pred)
accuracy_fitWCF_4 # -0.05
R2_fitWCF_4<- summary(Coef4)
R2_fitWCF_4 # 0.0006388
#CRPS
obs <- pheno4$milkZ
pred<- subset(pheno4,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCF_4 <- crps(obs,pred)
round((crps_fitWCF_4$CRPS),2) # 0.58


# Accuracy of prediction and degree under/overprediction of model fitWCF

accuracy_fitWCF = (accuracy_fitWCF_1 + accuracy_fitWCF_2 + accuracy_fitWCF_3 + accuracy_fitWCF_4)/4
round((accuracy_fitWCF),2) # 0.2

R2_fitWCF = (0.05932+0.1123+0.08083+0.0008917)/4
round((R2_fitWCF),2) # 0.06

# crps_fitWCF= (crps_fitWCF_1+crps_fitWCF_2+crps_fitWCF_3+crps_fitWCF_4)/4
crps_fitWCF = (0.6+ 0.53 + 0.52 + 0.58)/4
crps_fitWCF=round((crps_fitWCF),2)
crps_fitWCF #0.56

# Another method to calculate Accuracy, CRPS and R2 (degree of over/underprediction)
# Let's combine the predicted and observed phenotypes
pheno_fitWCF <- rbind(pheno1,pheno2,pheno3,pheno4)
accuracy_fitWCF<- round(cor(pheno_fitWCF$milkZ,pheno_fitWCF$milkZ_pred),2)
Coef_fitWCF <- lm (pheno_fitWCF$milkZ~pheno_fitWCF$milkZ_pred)
accuracy_fitWCF # 0.24
R2_fitWCF<- summary(Coef_fitWCF)
R2_fitWCF# 0.05623
#CRPS_fitWCF
obs <- pheno_fitWCF$milkZ
pred<- subset(pheno_fitWCF,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCF <- crps(obs,pred)
round((crps_fitWCF$CRPS),2) # 0.55


#Saving prediction models on External Drive
#save(fitWCF_NA1,file = "D:/Results_ADGG_Spatial/fitWCF_pred/fitWCF_NA1.RData")
#save(fitWCF_NA2,file = "D:/Results_ADGG_Spatial/fitWCF_pred/fitWCF_NA2.RData")
#save(fitWCF_NA3,file = "D:/Results_ADGG_Spatial/fitWCF_pred/fitWCF_NA3.RData")
#save(fitWCF_NA4,file = "D:/Results_ADGG_Spatial/fitWCF_pred/fitWCF_NA4.RData")

#-------------------------------------------------------------------------------

#----------- Accuracy model fitWCRI---------------------------------------------
# Masking phenotypes of cows in dgrp 1 (Breed proportion class 1)
#dgrp1 masked


fitWCRI_NA1 <- inla(formula = modelWCRI, data = data1_1_NA,
                  control.compute = list(dic = TRUE,config=TRUE))
# save fitWCRI_NA1 as R object
save(fitWCRI_NA1,file = "data/cleaned_data/fitWCRI_pred/fitWCRI_NA1.RData")
#load(file = "data/cleaned_data/fitWCRI_pred/fitWCRI_NA1.RData") # to load the R object
pheno_pred1 <- fitWCRI_NA1$summary.linear.predictor
colnames(pheno_pred1)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred1$mean
data1$milkZ_pred_sd <- pheno_pred1$sd
pheno1 <-   subset(data1, dgrp==1)
pheno1<- subset(pheno1,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno1$cowI))# 731 Cows in breed composition class 1

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotyped
accuracy_fitWCRI_1 <- round(cor(pheno1$milkZ,pheno1$milkZ_pred),2)
Coef1<- lm (pheno1$milkZ~pheno1$milkZ_pred)
accuracy_fitWCRI_1 #  0.42
R2_fitWCRI_1<- summary(Coef1)
R2_fitWCRI_1 =   0.1728
#CRPS
obs <- pheno1$milkZ
pred<- subset(pheno1,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCRI_1 <- crps(obs,pred)
round((crps_fitWCRI_1$CRPS),2) # 0.52

# Masking phenotypes of cows in dgrp 2 (Breed proportion class 2)
#dgrp2 masked

fitWCRI_NA2 <- inla(formula = modelWCRI, data = data1_2_NA,
                  control.compute = list(dic = TRUE,config=TRUE))
save(fitWCRI_NA2,file = "data/cleaned_data/fitWCRI_pred/fitWCRI_NA2.RData")
#load(file = "data/cleaned_data/fitWCRI_pred/fitWCRI_NA2.RData") # to load the R object

pheno_pred2 <- fitWCRI_NA2$summary.linear.predictor
colnames(pheno_pred2)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred2$mean
data1$milkZ_pred_sd <- pheno_pred2$sd

pheno2 <-   subset(data1, dgrp==2)
pheno2<- subset(pheno2,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno2$cowI))#  770 Cows in breed composition class 2

# Correlating observed phenotype of cows in dgrp 2 with predicted phenotypes
accuracy_fitWCRI_2 <- round(cor(pheno2$milkZ,pheno2$milkZ_pred),2)
Coef2<- lm (pheno2$milkZ~pheno2$milkZ_pred)
accuracy_fitWCRI_2# 0.51
R2_fitWCRI_2<- summary(Coef2)
R2_fitWCRI_2 = 0.2555
#CRPS
obs <- pheno2$milkZ
pred<- subset(pheno2,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCRI_2 <- crps(obs,pred)
round((crps_fitWCRI_2$CRPS),2) # 0.48

# Masking phenotypes of cows in dgrp 3 (Breed proportion class 3)
#dgrp3 masked

fitWCRI_NA3 <- inla(formula = modelWCRI, data = data1_3_NA,
                  control.compute = list(dic = TRUE,config=TRUE))

save(fitWCRI_NA3,file = "data/cleaned_data/fitWCRI_pred/fitWCRI_NA3.RData")
#load(file = "data/cleaned_data/fitWCRI_pred/fitWCRI_NA3.RData") # to load the R object

pheno_pred3 <- fitWCRI_NA3$summary.linear.predictor
sum(is.na(pheno_pred3$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred3$mean
data1$milkZ_pred_sd <- pheno_pred3$sd

pheno3 <-   subset(data1, dgrp==3)
pheno3<- subset(pheno3,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno3$cowI))#  309 Cows in breed composition class 3

# Correlating observed phenotype of cows in dgrp 3 with predicted phenotypes
accuracy_fitWCRI_3 <- round(cor(pheno3$milkZ,pheno3$milkZ_pred),2)
Coef3<- lm (pheno3$milkZ~pheno3$milkZ_pred)
accuracy_fitWCRI_3# 0.39
R2_fitWCRI_3<- summary(Coef3)
R2_fitWCRI_3 =0.1557
#CRPS
obs <- pheno3$milkZ
pred<- subset(pheno3,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCRI_3 <- crps(obs,pred)
round((crps_fitWCRI_3$CRPS),2) # 0.45


# Masking phenotypes of cows in dgrp 4 (Breed proportion class 4)
#dgrp4 masked

fitWCRI_NA4 <- inla(formula = modelWCRI, data = data1_4_NA,
                  control.compute = list(dic = TRUE,config=TRUE))
save(fitWCRI_NA4,file = "data/cleaned_data/fitWCRI_pred/fitWCRI_NA4.RData")
#load(file = "data/cleaned_data/fitWCRI_pred/fitWCRI_NA4.RData") # to load the R object

pheno_pred4 <- fitWCRI_NA4$summary.linear.predictor
sum(is.na(pheno_pred4$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred4$mean
data1$milkZ_pred_sd <- pheno_pred4$sd

pheno4 <-   subset(data1, dgrp==4)
pheno4<- subset(pheno4,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno4$cowI))#  84 Cows in breed composition class 4

# Correlating observed phenotype of cows in dgrp 4 with predicted phenotypes
accuracy_fitWCRI_4 <- round(cor(pheno4$milkZ,pheno4$milkZ_pred),2)
Coef4<- lm (pheno4$milkZ~pheno4$milkZ_pred)
accuracy_fitWCRI_4 # -0.05
R2_fitWCRI_4<- summary(Coef4)
R2_fitWCRI_4 = 0.04624
#CRPS
obs <- pheno4$milkZ
pred<- subset(pheno4,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCRI_4 <- crps(obs,pred)
round((crps_fitWCRI_4$CRPS),2) # 0.44


# Accuracy of prediction and degree under/overprediction of model fitWCRI

accuracy_fitWCRI = (accuracy_fitWCRI_1 + accuracy_fitWCRI_2 + accuracy_fitWCRI_3 + accuracy_fitWCRI_4)/4
round((accuracy_fitWCRI),2) # 0.38

R2_fitWCRI = (R2_fitWCRI_1 +R2_fitWCRI_2+R2_fitWCRI_3+R2_fitWCRI_4)/4
round((R2_fitWCRI),2) #0.16

#crps_fitWCRI= (crps_fitWCRI_1+crps_fitWCRI_2+crps_fitWCRI_3+crps_fitWCRI_4)/4
crps_fitWCRI =
crps_fitWCRI=round((crps_fitWCRI),2)
crps_fitWCRI #0.56

# Another method to calculate Accuracy, CRPS and R2 (degree of over/underprediction)
# Let's combine the predicted and observed phenotypes
pheno_fitWCRI <- rbind(pheno1,pheno2,pheno3,pheno4)
accuracy_fitWCRI<- round(cor(pheno_fitWCRI$milkZ,pheno_fitWCRI$milkZ_pred),2)
Coef_fitWCRI <- lm (pheno_fitWCRI$milkZ~pheno_fitWCRI$milkZ_pred)
accuracy_fitWCRI # 0.48
R2_fitWCRI<- summary(Coef_fitWCRI)
R2_fitWCRI# 0.2298
#CRPS_fitWCRI
obs <- pheno_fitWCRI$milkZ
pred<- subset(pheno_fitWCRI,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCRI <- crps(obs,pred)
round((crps_fitWCRI$CRPS),2) # 0.49

# Saving the prediction models on external Drive
save(fitWCRI_NA1,file = "D:/Results_ADGG_Spatial/fitWCRI_pred/fitWCRI_NA1.RData")
save(fitWCRI_NA2,file = "D:/Results_ADGG_Spatial/fitWCRI_pred/fitWCRI_NA2.RData")
save(fitWCRI_NA3,file = "D:/Results_ADGG_Spatial/fitWCRI_pred/fitWCRI_NA3.RData")
save(fitWCRI_NA4,file = "D:/Results_ADGG_Spatial/fitWCRI_pred/fitWCRI_NA4.RData")

#-------------------------------------------------------------------------------
#----------- Accuracy model fitWCRB---------------------------------------------

# Masking phenotypes of cows in dgrp 1 (Breed proportion class 1)
#dgrp1 masked

fitWCRB_NA1 <- inla(formula = modelWCRB, data = data1_1_NA,
                    control.compute = list(dic = TRUE,config=TRUE))
# save fitWCRB_NA1 as R object
save(fitWCRB_NA1,file = "data/cleaned_data/fitWCRB_pred/fitWCRB_NA1.RData")
#load(file = "data/cleaned_data/fitWCRB_pred/fitWCRB_NA1.RData") # to load the R object
pheno_pred1 <- fitWCRB_NA1$summary.linear.predictor
colnames(pheno_pred1)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred1$mean
data1$milkZ_pred_sd <- pheno_pred1$sd
pheno1 <-   subset(data1, dgrp==1)
pheno1<- subset(pheno1,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno1$cowI))# 731 Cows in breed composition class 1

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotyped
accuracy_fitWCRB_1 <- round(cor(pheno1$milkZ,pheno1$milkZ_pred),2)
Coef1<- lm (pheno1$milkZ~pheno1$milkZ_pred)
accuracy_fitWCRB_1 #  0.42
R2_fitWCRB_1<- summary(Coef1)
R2_fitWCRB_1 #   0.1739
#CRPS
obs <- pheno1$milkZ
pred<- subset(pheno1,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCRB_1 <- crps(obs,pred)
round((crps_fitWCRB_1$CRPS),2) # 0.52

# Masking phenotypes of cows in dgrp 2 (Breed proportion class 2)
#dgrp2 masked

fitWCRB_NA2 <- inla(formula = modelWCRB, data = data1_2_NA,
                    control.compute = list(dic = TRUE,config=TRUE))
save(fitWCRB_NA2,file = "data/cleaned_data/fitWCRB_pred/fitWCRB_NA2.RData")
#load(file = "data/cleaned_data/fitWCRB_pred/fitWCRB_NA2.RData") # to load the R object

pheno_pred2 <- fitWCRB_NA2$summary.linear.predictor
colnames(pheno_pred2)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred2$mean
data1$milkZ_pred_sd <- pheno_pred2$sd

pheno2 <-   subset(data1, dgrp==2)
pheno2<- subset(pheno2,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno2$cowI))#  770 Cows in breed composition class 2

# Correlating observed phenotype of cows in dgrp 2 with predicted phenotypes
accuracy_fitWCRB_2 <- round(cor(pheno2$milkZ,pheno2$milkZ_pred),2)
Coef2<- lm (pheno2$milkZ~pheno2$milkZ_pred)
accuracy_fitWCRB_2# 0.51
R2_fitWCRB_2<- summary(Coef2)
R2_fitWCRB_2 #  0.2564
#CRPS
obs <- pheno2$milkZ
pred<- subset(pheno2,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCRB_2 <- crps(obs,pred)
round((crps_fitWCRB_2$CRPS),2) # 0.48

# Masking phenotypes of cows in dgrp 3 (Breed proportion class 3)
#dgrp3 masked

fitWCRB_NA3 <- inla(formula = modelWCRB, data = data1_3_NA,
                    control.compute = list(dic = TRUE,config=TRUE))

save(fitWCRB_NA3,file = "data/cleaned_data/fitWCRB_pred/fitWCRB_NA3.RData")
#load(file = "data/cleaned_data/fitWCRB_pred/fitWCRB_NA3.RData") # to load the R object

pheno_pred3 <- fitWCRB_NA3$summary.linear.predictor
sum(is.na(pheno_pred3$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred3$mean
data1$milkZ_pred_sd <- pheno_pred3$sd

pheno3 <-   subset(data1, dgrp==3)
pheno3<- subset(pheno3,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno3$cowI))#  309 Cows in breed composition class 3

# Correlating observed phenotype of cows in dgrp 3 with predicted phenotypes
accuracy_fitWCRB_3 <- round(cor(pheno3$milkZ,pheno3$milkZ_pred),2)
Coef3<- lm (pheno3$milkZ~pheno3$milkZ_pred)
accuracy_fitWCRB_3# 0.39
R2_fitWCRB_3<- summary(Coef3)
R2_fitWCRB_3 # 0.1508
#CRPS
obs <- pheno3$milkZ
pred<- subset(pheno3,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCRB_3 <- crps(obs,pred)
round((crps_fitWCRB_3$CRPS),2) # 0.46


# Masking phenotypes of cows in dgrp 4 (Breed proportion class 4)
#dgrp4 masked

fitWCRB_NA4 <- inla(formula = modelWCRB, data = data1_4_NA,
                    control.compute = list(dic = TRUE,config=TRUE))
save(fitWCRB_NA4,file = "data/cleaned_data/fitWCRB_pred/fitWCRB_NA4.RData")
#load(file = "data/cleaned_data/fitWCRB_pred/fitWCRB_NA4.RData") # to load the R object

pheno_pred4 <- fitWCRB_NA4$summary.linear.predictor
sum(is.na(pheno_pred4$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred4$mean
data1$milkZ_pred_sd <- pheno_pred4$sd


pheno4 <-   subset(data1, dgrp==4)
pheno4<- subset(pheno4,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno4$cowI))#  84 Cows in breed composition class 4

# Correlating observed phenotype of cows in dgrp 4 with predicted phenotypes
accuracy_fitWCRB_4 <- round(cor(pheno4$milkZ,pheno4$milkZ_pred),2)
Coef4<- lm (pheno4$milkZ~pheno4$milkZ_pred)
accuracy_fitWCRB_4 # 0.23
R2_fitWCRB_4<- summary(Coef4)
R2_fitWCRB_4 # 0.05038
#CRPS
obs <- pheno4$milkZ
pred<- subset(pheno4,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCRB_4 <- crps(obs,pred)
round((crps_fitWCRB_4$CRPS),2) # 0.44


# Accuracy of prediction and degree under/overprediction of model fitWCRB

accuracy_fitWCRB = (accuracy_fitWCRB_1 + accuracy_fitWCRB_2 + accuracy_fitWCRB_3 + accuracy_fitWCRB_4)/4
round((accuracy_fitWCRB),2) # 0.39

R2_fitWCRB = ( 0.1739 +0.2564+0.1508+0.05038 )/4
round((R2_fitWCRB),2) # 0.16

# crps_fitWCRB= (crps_fitWCRB_1+crps_fitWCRB_2+crps_fitWCRB_3+crps_fitWCRB_4)/4
crps_fitWCRB = ( 0.52+ 0.48 +  0.46+ 0.44)/4
crps_fitWCRB=round((crps_fitWCRB),2)
crps_fitWCRB #0.48

# Another method to calculate Accuracy, CRPS and R2 (degree of over/underprediction)
# Let's combine the predicted and observed phenotypes
pheno_fitWCRB <- rbind(pheno1,pheno2,pheno3,pheno4)
accuracy_fitWCRB<- round(cor(pheno_fitWCRB$milkZ,pheno_fitWCRB$milkZ_pred),2)
Coef_fitWCRB <- lm (pheno_fitWCRB$milkZ~pheno_fitWCRB$milkZ_pred)
accuracy_fitWCRB # 0.48
R2_fitWCRB<- summary(Coef_fitWCRB)
R2_fitWCRB#  0.232
#CRPS_fitWCRB
obs <- pheno_fitWCRB$milkZ
pred<- subset(pheno_fitWCRB,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCRB <- crps(obs,pred)
round((crps_fitWCRB$CRPS),2) # 0.49
#-------------------------------------------------------------------------------

#----------- Accuracy model fitS------------------------------------------------
# Masking phenotypes of cows in dgrp 1 (Breed proportion class 1)
#dgrp1 masked

mesh = inla.mesh.2d(cbind(data1$long, data1$lat), max.edge=c(10, 20), cutoff = 2.5, offset = 30)
A.pred = inla.spde.make.A(mesh = mesh, loc = cbind(data1$long, data1$lat) )
spdeStatS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeS)
meshIndexS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatS$n.spde)

# Make stack_pred
# StackS_NA1
stackS_pred_NA1 = inla.stack(data = list(milkZ = NA),
                              A = list(A.pred,1),
                              effects = list(c(meshIndexS, list(intercept = 1)),
                                             list(cowI = data1_1_NA$cowI,herdI = data1_1_NA$herdI, cowPeI = data1_1_NA$cowPeI,
                                                  cyrsnI=data1_1_NA$cyrsnI, tyrmnI=data1_1_NA$tyrmnI,dgrpI= data1_1_NA$dgrpI, ageZ=data1_1_NA$ageZ, lacgr=data1_1_NA$lacgr, leg0=data1_1_NA$leg0, leg1=data1_1_NA$leg1,leg2=data1_1_NA$leg2)), tag = "data1_1_NA.data")

# Create joint stack
join.stack_NA1 <- inla.stack(stackS, stackS_pred_NA1)

# ModelS
formulaS <- as.formula(paste0(modelBase, " + f(fieldID, model = spdeStatS)-1"))


fitS_NA1= inla(formula = formulaS, data = inla.stack.data(join.stack_NA1),
           family = "normal", control.predictor =list(A=inla.stack.A(join.stack_NA1),compute = T),
           control.family=list(list(hyper=hyperResVarGWS)),
           control.compute = list(dic=T,cpo=F, config=T),
           control.fixed = list(expand.factor.strategy="inla"), verbose=T)

# save fitS_NA1 as R object
save(fitS_NA1,file = "data/cleaned_data/fitS_pred/fitS_NA1.RData")
#load(file = "data/cleaned_data/fitS_pred/fitS_NA1.RData") # to load the R object

#Get the prediction index

pred.ind_NA1 <- inla.stack.index(join.stack_NA1, tag='data1_1_NA.data')$data
pheno_pred1 <- fitS_NA1$summary.linear.predictor [pred.ind_NA1,]

sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred1$mean
data1$milkZ_pred_sd <- pheno_pred1$sd
pheno1 <-   subset(data1, dgrp==1)
pheno1<- subset(pheno1,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno1$cowI))# 731 Cows in breed composition class 1

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotyped
accuracy_fitS_1 <- round(cor(pheno1$milkZ,pheno1$milkZ_pred),2)
Coef1<- lm (pheno1$milkZ~pheno1$milkZ_pred)
accuracy_fitS_1 #  0.81
R2_fitS_1<- summary(Coef1)
R2_fitS_1 #  0.6536
#CRPS
obs <- pheno1$milkZ
pred<- subset(pheno1,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitS_1 <- crps(obs,pred)
round((crps_fitS_1$CRPS),2) # 0.35

# Masking phenotypes of cows in dgrp 2 (Breed proportion class 2)
#dgrp2 masked

# Make stack_pred
# StackS_NA2
stackS_pred_NA2 = inla.stack(data = list(milkZ = NA),
                             A = list(A.pred,1),
                             effects = list(c(meshIndexS, list(intercept = 1)),
                                            list(cowI = data1_2_NA$cowI,herdI = data1_2_NA$herdI, cowPeI = data1_2_NA$cowPeI,
                                                 cyrsnI=data1_2_NA$cyrsnI, tyrmnI=data1_2_NA$tyrmnI,dgrpI= data1_2_NA$dgrpI, ageZ=data1_2_NA$ageZ, lacgr=data1_2_NA$lacgr, leg0=data1_2_NA$leg0, leg1=data1_2_NA$leg1,leg2=data1_2_NA$leg2)), tag = "data1_2_NA.data")

# Create joint stack
join.stack_NA2 <- inla.stack(stackS, stackS_pred_NA2)

# ModelS
formulaS <- as.formula(paste0(modelBase, " + f(fieldID, model = spdeStatS)-1"))


fitS_NA2= inla(formula = formulaS, data = inla.stack.data(join.stack_NA2),
               family = "normal", control.predictor =list(A=inla.stack.A(join.stack_NA2),compute = T),
               control.family=list(list(hyper=hyperResVarGWS)),
               control.compute = list(dic=T,cpo=F, config=T),
               control.fixed = list(expand.factor.strategy="inla"), verbose=T)




save(fitS_NA2,file = "data/cleaned_data/fitS_pred/fitS_NA2.RData")
#load(file = "data/cleaned_data/fitS_pred/fitS_NA2.RData") # to load the R object

#Get the prediction index

pred.ind_NA2 <- inla.stack.index(join.stack_NA2, tag='data1_2_NA.data')$data
pheno_pred2 <- fitS_NA2$summary.linear.predictor [pred.ind_NA2,]


sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred2$mean
data1$milkZ_pred_sd <- pheno_pred2$sd

pheno2 <-   subset(data1, dgrp==2)
pheno2<- subset(pheno2,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno2$cowI))#  770 Cows in breed composition class 2

# Correlating observed phenotype of cows in dgrp 2 with predicted phenotypes
accuracy_fitS_2 <- round(cor(pheno2$milkZ,pheno2$milkZ_pred),2)
Coef2<- lm (pheno2$milkZ~pheno2$milkZ_pred)
accuracy_fitS_2# 0.82
R2_fitS_2<- summary(Coef2)
R2_fitS_2 #  0.6749
#CRPS
obs <- pheno2$milkZ
pred<- subset(pheno2,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitS_2 <- crps(obs,pred)
round((crps_fitS_2$CRPS),2) # 0.34

# Masking phenotypes of cows in dgrp 3 (Breed proportion class 3)
#dgrp3 masked
# Make stack_pred
# StackS_NA3
stackS_pred_NA3 = inla.stack(data = list(milkZ = NA),
                             A = list(A.pred,1),
                             effects = list(c(meshIndexS, list(intercept = 1)),
                                            list(cowI = data1_3_NA$cowI,herdI = data1_3_NA$herdI, cowPeI = data1_3_NA$cowPeI,
                                                 cyrsnI=data1_3_NA$cyrsnI, tyrmnI=data1_3_NA$tyrmnI,dgrpI= data1_3_NA$dgrpI, ageZ=data1_3_NA$ageZ, lacgr=data1_3_NA$lacgr, leg0=data1_3_NA$leg0, leg1=data1_3_NA$leg1,leg2=data1_3_NA$leg2)), tag = "data1_3_NA.data")

# Create joint stack
join.stack_NA3 <- inla.stack(stackS, stackS_pred_NA3)

# ModelS
formulaS <- as.formula(paste0(modelBase, " + f(fieldID, model = spdeStatS)-1"))


fitS_NA3= inla(formula = formulaS, data = inla.stack.data(join.stack_NA3),
               family = "normal", control.predictor =list(A=inla.stack.A(join.stack_NA3),compute = T),
               control.family=list(list(hyper=hyperResVarGWS)),
               control.compute = list(dic=T,cpo=F, config=T),
               control.fixed = list(expand.factor.strategy="inla"), verbose=T)

save(fitS_NA3,file = "data/cleaned_data/fitS_pred/fitS_NA3.RData")
#load(file = "data/cleaned_data/fitS_pred/fitS_NA3.RData") # to load the R object

#Get the prediction index

pred.ind_NA3 <- inla.stack.index(join.stack_NA3, tag='data1_3_NA.data')$data
pheno_pred3 <- fitS_NA3$summary.linear.predictor [pred.ind_NA3,]


sum(is.na(pheno_pred3$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred3$mean
data1$milkZ_pred_sd <- pheno_pred3$sd

pheno3 <-   subset(data1, dgrp==3)
pheno3<- subset(pheno3,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno3$cowI))#  309 Cows in breed composition class 3

# Correlating observed phenotype of cows in dgrp 3 with predicted phenotypes
accuracy_fitS_3 <- round(cor(pheno3$milkZ,pheno3$milkZ_pred),2)
Coef3<- lm (pheno3$milkZ~pheno3$milkZ_pred)
accuracy_fitS_3# 0.78
R2_fitS_3<- summary(Coef3)
R2_fitS_3 #0.6118
#CRPS
obs <- pheno3$milkZ
pred<- subset(pheno3,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitS_3 <- crps(obs,pred)
round((crps_fitS_3$CRPS),2) # 0.31


# Masking phenotypes of cows in dgrp 4 (Breed proportion class 4)
#dgrp4 masked
A = inla.spde.make.A(mesh = mesh, loc = cbind(data1_4_NA$long, data1_4_NA$lat) )
spdeStatS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeS)
meshIndexS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatS$n.spde)

# Make stack
# StackS_NA4
# Make stack_pred
stackS_pred_NA4 = inla.stack(data = list(milkZ = NA),
                             A = list(A.pred,1),
                             effects = list(c(meshIndexS, list(intercept = 1)),
                                            list(cowI = data1_4_NA$cowI,herdI = data1_4_NA$herdI, cowPeI = data1_4_NA$cowPeI,
                                                 cyrsnI=data1_4_NA$cyrsnI, tyrmnI=data1_4_NA$tyrmnI,dgrpI= data1_4_NA$dgrpI, ageZ=data1_4_NA$ageZ, lacgr=data1_4_NA$lacgr, leg0=data1_4_NA$leg0, leg1=data1_4_NA$leg1,leg2=data1_4_NA$leg2)), tag = "data1_4_NA.data")

# Create joint stack
join.stack_NA4 <- inla.stack(stackS, stackS_pred_NA4)

# ModelS
formulaS <- as.formula(paste0(modelBase, " + f(fieldID, model = spdeStatS)-1"))


fitS_NA4= inla(formula = formulaS, data = inla.stack.data(join.stack_NA4),
               family = "normal", control.predictor =list(A=inla.stack.A(join.stack_NA4),compute = T),
               control.family=list(list(hyper=hyperResVarGWS)),
               control.compute = list(dic=T,cpo=F, config=T),
               control.fixed = list(expand.factor.strategy="inla"), verbose=T)

save(fitS_NA4,file = "data/cleaned_data/fitS_pred/fitS_NA4.RData")
#load(file = "data/cleaned_data/fitS_pred/fitS_NA4.RData") # to load the R object

#Get the prediction index

pred.ind_NA4 <- inla.stack.index(join.stack_NA4, tag='data1_4_NA.data')$data
pheno_pred4 <- fitS_NA4$summary.linear.predictor [pred.ind_NA4,]

sum(is.na(pheno_pred4$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred4$mean
data1$milkZ_pred_sd <- pheno_pred4$sd

pheno4 <-   subset(data1, dgrp==4)
pheno4<- subset(pheno4,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno4$cowI))#  84 Cows in breed composition class 4

# Correlating observed phenotype of cows in dgrp 4 with predicted phenotypes
accuracy_fitS_4 <- round(cor(pheno4$milkZ,pheno4$milkZ_pred),2)
Coef4<- lm (pheno4$milkZ~pheno4$milkZ_pred)
accuracy_fitS_4 # 0.75
R2_fitS_4<- summary(Coef4)
R2_fitS_4 # # 0.5626
#CRPS
obs <- pheno4$milkZ
pred<- subset(pheno4,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitS_4 <- crps(obs,pred)
round((crps_fitS_4$CRPS),2) # 0.25


# Accuracy of prediction and degree under/overprediction of model fitS

accuracy_fitS = (accuracy_fitS_1 + accuracy_fitS_2 + accuracy_fitS_3 + accuracy_fitS_4)/4
round((accuracy_fitS),2) # 0.79

R2_fitS = ( 0.6536+ 0.6749+ 0.6118 +0.5626 )/4
round((R2_fitS),2) # 0.63

# crps_fitS= (crps_fitS_1+crps_fitS_2+crps_fitS_3+crps_fitS_4)/4
crps_fitS = (0.35+ 0.34 + 0.31 + 0.25)/4
crps_fitS=round((crps_fitS),2)
crps_fitS #0.31

# Another method to calculate Accuracy, CRPS and R2 (degree of over/underprediction)
# Let's combine the predicted and observed phenotypes
pheno_fitS <- rbind(pheno1,pheno2,pheno3,pheno4)
accuracy_fitS<- round(cor(pheno_fitS$milkZ,pheno_fitS$milkZ_pred),2)
Coef_fitS <- lm (pheno_fitS$milkZ~pheno_fitS$milkZ_pred)
accuracy_fitS # 0.83
R2_fitS<- summary(Coef_fitS)
R2_fitS# 0.6907
#CRPS_fitS
obs <- pheno_fitS$milkZ
pred<- subset(pheno_fitS,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitS <- crps(obs,pred)
round((crps_fitS$CRPS),2) # 0.34

#Saving prediction models on External Drive
save(fitS_NA1,file = "D:/Results_ADGG_Spatial/fitS_pred/fitS_NA1.RData")
save(fitS_NA2,file = "D:/Results_ADGG_Spatial/fitS_pred/fitS_NA2.RData")
save(fitS_NA3,file = "D:/Results_ADGG_Spatial/fitS_pred/fitS_NA3.RData")
save(fitS_NA4,file = "D:/Results_ADGG_Spatial/fitS_pred/fitS_NA4.RData")


#-------------------------------------------------------------------------------
#----------- Accuracy model fitWS-----------------------------------------------
# Masking phenotypes of cows in dgrp 1 (Breed proportion class 1)
#dgrp1 masked

#Create data structure for prediction
mesh = inla.mesh.2d(cbind(data1$long, data1$lat), max.edge=c(10, 20), cutoff = 2.5, offset = 30)
A.pred = inla.spde.make.A(mesh = mesh, loc = cbind(data1$long, data1$lat) )
spdeStatWS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeWS)
meshIndexWS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatWS$n.spde)

# Make stack_pred
# StackWS_NA1
stackWS_pred_NA1 = inla.stack(data = list(milkZ = NA),
                        A = list(A.pred,1),
                        effects = list(c(meshIndexWS, list(intercept = 1)),
                                       list(cowI = data1_1_NA$cowI,herdI = data1_1_NA$herdI, cowPeI = data1_1_NA$cowPeI,
                                            cyrsnI=data1_1_NA$cyrsnI, tyrmnI=data1_1_NA$tyrmnI,dgrpI= data1_1_NA$dgrpI, ageZ=data1_1_NA$ageZ, lacgr=data1_1_NA$lacgr, leg0=data1_1_NA$leg0, leg1=data1_1_NA$leg1,leg2=data1_1_NA$leg2)), tag = "data1_1_NA.data")

# Create joint stack
join.stack_NA1 <- inla.stack(stack, stackWS_pred_NA1)

# ModelWS
formulaWS <- as.formula(paste0(modelBase, " + f(fieldID, model = spdeStatWS)-1"))


fitWS_NA1= inla(formula = formulaWS, data = inla.stack.data(join.stack_NA1),
               family = "normal", control.predictor =list(A=inla.stack.A(join.stack_NA1),compute = T),
               control.family=list(list(hyper=hyperResVarGWS)),
               control.compute = list(dic=T,cpo=F, config=T),
               control.fixed = list(expand.factor.strategy="inla"), verbose=T)

# save fitWS_NA1 as R object
save(fitWS_NA1,file = "data/cleaned_data/fitWS_pred/fitWS_NA1.RData")
#load(file = "data/cleaned_data/fitWS_pred/fitWS_NA1.RData") # to load the R object

#Get the prediction index

pred.ind_NA1 <- inla.stack.index(join.stack_NA1, tag='data1_1_NA.data')$data
pheno_pred1 <- fitWS_NA1$summary.linear.predictor [pred.ind_NA1,]
colnames(pheno_pred1)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred1$mean
data1$milkZ_pred_sd <- pheno_pred1$sd
pheno1 <-   subset(data1, dgrp==1)
pheno1<- subset(pheno1,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno1$cowI))# 731 Cows in breed composition class 1


# Correlating observed phenotype of cows in dgrp 1 with predicted phenotyped
accuracy_fitWS_1 <- round(cor(pheno1$milkZ,pheno1$milkZ_pred),2)
Coef1<- lm (pheno1$milkZ~pheno1$milkZ_pred)
accuracy_fitWS_1 #  0.81
R2_fitWS_1<- summary(Coef1)
R2_fitWS_1 #  0.6536
#CRPS
obs <- pheno1$milkZ
pred<- subset(pheno1,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWS_1 <- crps(obs,pred)
round((crps_fitWS_1$CRPS),2) # 0.35

# Masking phenotypes of cows in dgrp 2 (Breed proportion class 2)
#dgrp2 masked
#Create data structure for prediction
#idem as in dgrp1 masked

# Make stack_pred
# StackS_NA2
stackWS_pred_NA2 = inla.stack(data = list(milkZ = NA),
                              A = list(A.pred,1),
                              effects = list(c(meshIndexWS, list(intercept = 1)),
                                             list(cowI = data1_2_NA$cowI,herdI = data1_2_NA$herdI, cowPeI = data1_2_NA$cowPeI,
                                                  cyrsnI=data1_2_NA$cyrsnI, tyrmnI=data1_2_NA$tyrmnI,dgrpI= data1_2_NA$dgrpI, ageZ=data1_2_NA$ageZ, lacgr=data1_2_NA$lacgr, leg0=data1_2_NA$leg0, leg1=data1_2_NA$leg1,leg2=data1_2_NA$leg2)), tag = "data1_2_NA.data")

# Create joint stack
join.stack_NA2 <- inla.stack(stackWS, stackWS_pred_NA2)


# ModelWS
formulaWS <- as.formula(paste0(modelBase, " + f(fieldID, model = spdeStatWS)-1"))


fitWS_NA2= inla(formula = formulaWS, data = inla.stack.data(join.stack_NA2),
               family = "normal", control.predictor =list(A=inla.stack.A(join.stack_NA2),compute = T),
               control.family=list(list(hyper=hyperResVarGWS)),
               control.compute = list(dic=T,cpo=F, config=T),
               control.fixed = list(expand.factor.strategy="inla"), verbose=T)


save(fitWS_NA2,file = "data/cleaned_data/fitWS_pred/fitWS_NA2.RData")
#load(file = "data/cleaned_data/fitWS_pred/fitWS_NA2.RData") # to load the R object

#Get the prediction index

pred.ind_NA2 <- inla.stack.index(join.stack_NA2, tag='data1_2_NA.data')$data
pheno_pred2 <- fitWS_NA2$summary.linear.predictor [pred.ind_NA2,]
colnames(pheno_pred2)
sum(is.na(pheno_pred2$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred2$mean
data1$milkZ_pred_sd <- pheno_pred2$sd

pheno2 <-   subset(data1, dgrp==2)
pheno2<- subset(pheno2,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno2$cowI))#  770 Cows in breed composition class 2

# Correlating observed phenotype of cows in dgrp 2 with predicted phenotypes
accuracy_fitWS_2 <- round(cor(pheno2$milkZ,pheno2$milkZ_pred),2)
Coef2<- lm (pheno2$milkZ~pheno2$milkZ_pred)
accuracy_fitWS_2# 0.82
R2_fitWS_2<- summary(Coef2)
R2_fitWS_2 # 0.675
#CRPS
obs <- pheno2$milkZ
pred<- subset(pheno2,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWS_2 <- crps(obs,pred)
round((crps_fitWS_2$CRPS),2) # 0.34

# Masking phenotypes of cows in dgrp 3 (Breed proportion class 3)
#dgrp3 masked
# Make stack_pred
# StackS_NA3
stackWS_pred_NA3 = inla.stack(data = list(milkZ = NA),
                              A = list(A.pred,1),
                              effects = list(c(meshIndexWS, list(intercept = 1)),
                                             list(cowI = data1_3_NA$cowI,herdI = data1_3_NA$herdI, cowPeI = data1_3_NA$cowPeI,
                                                  cyrsnI=data1_3_NA$cyrsnI, tyrmnI=data1_3_NA$tyrmnI,dgrpI= data1_3_NA$dgrpI, ageZ=data1_3_NA$ageZ, lacgr=data1_3_NA$lacgr, leg0=data1_3_NA$leg0, leg1=data1_3_NA$leg1,leg2=data1_3_NA$leg2)), tag = "data1_3_NA.data")

# Create joint stack
join.stack_NA3 <- inla.stack(stackWS, stackWS_pred_NA3)

# ModelWS
formulaWS <- as.formula(paste0(modelBase, " + f(fieldID, model = spdeStatWS)-1"))


fitWS_NA3= inla(formula = formulaWS, data = inla.stack.data(join.stack_NA3),
               family = "normal", control.predictor =list(A=inla.stack.A(join.stack_NA3),compute = T),
               control.family=list(list(hyper=hyperResVarGWS)),
               control.compute = list(dic=T,cpo=F, config=T),
               control.fixed = list(expand.factor.strategy="inla"), verbose=T)

save(fitWS_NA3,file = "data/cleaned_data/fitWS_pred/fitWS_NA3.RData")
#load(file = "data/cleaned_data/fitWS_pred/fitWS_NA3.RData") # to load the R object

#Get the prediction index
pred.ind_NA3 <- inla.stack.index(join.stack_NA3, tag='data1_3_NA.data')$data
pheno_pred3 <- fitWS_NA3$summary.linear.predictor [pred.ind_NA3,]
colnames(pheno_pred3)
sum(is.na(pheno_pred3$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred3$mean
data1$milkZ_pred_sd <- pheno_pred3$sd

pheno3 <-   subset(data1, dgrp==3)
pheno3<- subset(pheno3,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno3$cowI))#  309 Cows in breed composition class 3

# Correlating observed phenotype of cows in dgrp 3 with predicted phenotypes
accuracy_fitWS_3 <- round(cor(pheno3$milkZ,pheno3$milkZ_pred),2)
Coef3<- lm (pheno3$milkZ~pheno3$milkZ_pred)
accuracy_fitWS_3 # 0.78
R2_fitWS_3<- summary(Coef3)
R2_fitWS_3 #  0.6118
#CRPS
obs <- pheno3$milkZ
pred<- subset(pheno3,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWS_3 <- crps(obs,pred)
round((crps_fitWS_3$CRPS),2) #  0.31


# Masking phenotypes of cows in dgrp 4 (Breed proportion class 4)
#dgrp4 masked
# Make stack_pred
# StackS_NA3
stackWS_pred_NA4 = inla.stack(data = list(milkZ = NA),
                              A = list(A.pred,1),
                              effects = list(c(meshIndexWS, list(intercept = 1)),
                                             list(cowI = data1_4_NA$cowI,herdI = data1_4_NA$herdI, cowPeI = data1_4_NA$cowPeI,
                                                  cyrsnI=data1_4_NA$cyrsnI, tyrmnI=data1_4_NA$tyrmnI,dgrpI= data1_4_NA$dgrpI, ageZ=data1_4_NA$ageZ, lacgr=data1_4_NA$lacgr, leg0=data1_4_NA$leg0, leg1=data1_4_NA$leg1,leg2=data1_4_NA$leg2)), tag = "data1_4_NA.data")

# Create joint stack
join.stack_NA4 <- inla.stack(stackWS, stackWS_pred_NA4)



# ModelWS
formulaWS <- as.formula(paste0(modelBase, " + f(fieldID, model = spdeStatWS)-1"))


fitWS_NA4= inla(formula = formulaWS, data = inla.stack.data(join.stack_NA4),
               family = "normal", control.predictor =list(A=inla.stack.A(join.stack_NA4),compute = T),
               control.family=list(list(hyper=hyperResVarGWS)),
               control.compute = list(dic=T,cpo=F, config=T),
               control.fixed = list(expand.factor.strategy="inla"), verbose=T)

save(fitWS_NA4,file = "data/cleaned_data/fitWS_pred/fitWS_NA4.RData")
#load(file = "data/cleaned_data/fitWS_pred/fitWS_NA4.RData") # to load the R object

#Get the prediction index
pred.ind_NA4 <- inla.stack.index(join.stack_NA4, tag='data1_4_NA.data')$data
pheno_pred4 <- fitWS_NA4$summary.linear.predictor [pred.ind_NA4,]
colnames(pheno_pred4)
sum(is.na(pheno_pred4$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred4$mean
data1$milkZ_pred_sd <- pheno_pred4$sd

pheno4 <-   subset(data1, dgrp==4)
pheno4<- subset(pheno4,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno4$cowI))#  84 Cows in breed composition class 4

# Correlating observed phenotype of cows in dgrp 4 with predicted phenotypes
accuracy_fitWS_4 <- round(cor(pheno4$milkZ,pheno4$milkZ_pred),2)
Coef4<- lm (pheno4$milkZ~pheno4$milkZ_pred)
accuracy_fitWS_4 # 0.75
R2_fitWS_4<- summary(Coef4)
R2_fitWS_4 # 0.5626
#CRPS
obs <- pheno4$milkZ
pred<- subset(pheno4,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWS_4 <- crps(obs,pred)
round((crps_fitWS_4$CRPS),2) # 0.25


# Accuracy of prediction and degree under/overprediction of model fitWS

accuracy_fitWS = (accuracy_fitWS_1 + accuracy_fitWS_2 + accuracy_fitWS_3 + accuracy_fitWS_4)/4
accuracy_fitWS ## 0.79

R2_fitWS =  (0.6536 + 0.675+ 0.6118  +0.5626)/4
R2_fitWS
round((R2_fitWS),2) #0.63

crps_fitWS= (crps_fitWS_1$CRPS+crps_fitWS_2$CRPS+crps_fitWS_3$CRPS+crps_fitWS_4$CRPS)/4
crps_fitWS
crps_fitWS=round((crps_fitWS),2)
crps_fitWS #0.31

# Another method to calculate Accuracy, CRPS and R2 (degree of over/underprediction)
# Let's combine the predicted and observed phenotypes
pheno_fitWS <- rbind(pheno1,pheno2,pheno3,pheno4)
accuracy_fitWS<- round(cor(pheno_fitWS$milkZ,pheno_fitWS$milkZ_pred),2)
Coef_fitWS <- lm (pheno_fitWS$milkZ~pheno_fitWS$milkZ_pred)
accuracy_fitWS # 0.83
R2_fitWS<- summary(Coef_fitWS)
R2_fitWS# 0.6907
#CRPS_fitWS
obs <- pheno_fitWS$milkZ
pred<- subset(pheno_fitWS,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWS <- crps(obs,pred)
round((crps_fitWS$CRPS),2) # 0.34

#Saving prediction models on External Drive (TO DO)
save(fitWS_NA1,file = "D:/Results_ADGG_Spatial/fitWS_pred/fitWS_NA1.RData")
save(fitWS_NA2,file = "D:/Results_ADGG_Spatial/fitWS_pred/fitWS_NA2.RData")
save(fitWS_NA3,file = "D:/Results_ADGG_Spatial/fitWS_pred/fitWS_NA3.RData")
save(fitWS_NA4,file = "D:/Results_ADGG_Spatial/fitWS_pred/fitWS_NA4.RData")

#------Accuracy of prediction: Foward validation--------------------------------

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


#-----------------------Foward validation FitB-------------------------------
#Let's make NA milkz of cows born in 2016 and 2017
data2_NA<- data2%>% mutate(milkZ = ifelse(birthyear=="2016" | birthyear=="2017", NA, milkZ))
sum(is.na(data2_NA$milkZ))   #1424 Expected
length(unique(data2_NA$cow)) #1894 Expected
# Yield deviation of all cows

predictB <- predict(fitB, newdata=data2,
                   formula= ~ fixed_effects + perm)
predictB2 <- predictB
colnames(predictB)
summary(data2$milkZ)
summary(predictB$milkZ)
summary(predictB$mean)
round(cor(predictB$milkZ,predictB$mean),2) # 0.77
predictB <- predictB[,c("cow","milkZ","mean")]
predictB <- predictB %>% mutate(YDB=milkZ-mean) %>%
 group_by(cow) %>%
  summarise(YDB_cow=mean(YDB))

# Breeding values and Yield deviation model fitB

# Predict Breeding value of young cows
fitB_NA <- bru(modelB,
                       like(family = "Gaussian",
                            modelBFormula,
                            data = data2_NA))

EBV_B_NA <- fitB_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_B_NA)[1] <- "cow"
colnames(EBV_B_NA)[2] <- "EBV_B"

# data2 for 2016 and 2017
data2_2016_2017 <- subset(data2, birthyear==2016 | birthyear==2017 )
length(unique(data2_2016_2017$cow)) # 146 cows
# Young cows IDS  (born in 2016 and 2017)
young_cow_ID <- subset(data2_2016_2017, select=c(cow))
young_cow_ID <- unique(young_cow_ID)
summary(young_cow_ID)
rownames(young_cow_ID) <- NULL

sel <- predictB$cow %in% young_cow_ID$cow
predictB <- predictB[sel,]
predictB<- subset(predictB,select= c("cow","YDB_cow"))
predictB <- data.frame(predictB)
predictB <- predictB[,c("cow","YDB_cow")]

sel <- EBV_B_NA$cow %in% young_cow_ID$cow
EBV_B_NA_young<- EBV_B_NA[sel,]

EBV_YDB <- merge(predictB,EBV_B_NA_young, by="cow", all.x=TRUE)

# Forward Validation for accuracy for model fitB
round(cor(EBV_YDB$YDB_cow,EBV_YDB$EBV_B),2) # 0.13

#------ Forward Validation model fitBH------------------------------------------
predictBH <- predict(fitBH, newdata=data2,
                    formula= ~fixed_effects + perm)
predictBH2 <- predictBH # Saving a copy of predictBH

summary(fitBH)
colnames(predictBH)
summary(data2$milkZ)
summary(predictBH$milkZ)
summary(predictBH$mean)
plot(predictBH$mean,predictBH$milkZ)
round(cor(predictBH$mean,predictBH$milkZ),2) # 0.64
predictBH <- predictBH[,c("cow","milkZ","mean")]
predictBH <- predictBH %>% mutate(YDBH=milkZ-mean) %>%
  group_by(cow) %>%
  summarise(YDBH_cow=mean(YDBH))

# Predict Breeding value of young cows
fitBH_NA <- bru(modelBH,
               like(family = "Gaussian",
                    modelBHFormula,
                    data = data2_NA))

EBV_BH_NA <- fitBH_NA$summary.random$animal[,c("ID", "mean")]

colnames(EBV_BH_NA)[1] <- "cow"
colnames(EBV_BH_NA)[2] <- "EBV_BH"

sel <- predictBH$cow %in% young_cow_ID$cow
predictBH <- predictBH[sel,]
predictBH<- subset(predictBH,select= c("cow","YDBH_cow"))
predictBH <- data.frame(predictBH)
predictBH <- predictBH[,c("cow", "YDBH_cow")]

sel <- EBV_BH_NA$cow %in% young_cow_ID$cow
EBV_BH_NA_young<- EBV_BH_NA[sel,]

EBV_YDBH <- merge(predictBH,EBV_BH_NA_young, by="cow", all.x=TRUE)

# Forward Validation for accuracy for model fitBH
round(cor(EBV_YDBH$YDBH_cow,EBV_YDBH$EBV_BH),2) # -0.03

#------ Forward Validation model fitBS------------------------------------------
predictBS <- predict(fitBS, newdata=data2,
                     formula= ~fixed_effects + perm)
summary(fitBS) # 89222.12 (non scaled)

predictBS2 <- predictBS # Saving a copy of predictBS

colnames(predictBS)
summary(data2$milkZ)
summary(predictBS$milkZ)
summary(predictBS$mean)
plot(predictBS$mean,predictBS$milk)
round(cor(predictBS$mean,predictBS$milk),2) #0.67
predictBS <- predictBS[,c("cow","milkZ","mean")]
predictBS <- predictBS %>% mutate(YDBS=milkZ-mean) %>%
  group_by(cow) %>%
  summarise(YDBS_cow=mean(YDBS))


# Predict Breeding value of young cows
fitBS_NA <- bru(modelBS,
                like(family = "Gaussian",
                     modelBSFormula,
                     data = data2_NA))

EBV_BS_NA <- fitBS_NA$summary.random$animal[,c("ID", "mean")]

colnames(EBV_BS_NA)[1] <- "cow"
colnames(EBV_BS_NA)[2] <- "EBV_BS"

sel <- predictBS$cow %in% young_cow_ID$cow
predictBS <- predictBS[sel,]
predictBS<- subset(predictBS,select= c("cow","YDBS_cow"))
predictBS <- data.frame(predictBS)
predictBS <- predictBS[,c("cow", "YDBS_cow")]

sel <- EBV_BS_NA$cow %in% young_cow_ID$cow
EBV_BS_NA_young<- EBV_BS_NA[sel,]

EBV_YDBS <- merge(predictBS,EBV_BS_NA_young, by="cow", all.x=TRUE)

# Forward Validation for accuracy for model fitBS
round(cor(EBV_YDBS$YDBS_cow,EBV_YDBS$EBV_BS),2) # 0.09

#------ Forward Validation model fitBHS------------------------------------------
predictBHS <- predict(fitBHS, newdata=data2,
                     formula= ~fixed_effects + perm)

predictBHS2 <- predictBHS # Saving a copy of predictBS

colnames(predictBHS)
summary(data2$milkZ)
summary(predictBHS$milkZ)
summary(predictBHS$mean)
plot(predictBHS$mean,predictBHS$milkZ)
round(cor(predictBHS$mean,predictBHS$milk),2) # 0.63
predictBHS <- predictBHS[,c("cow","milkZ","mean")]
predictBHS <- predictBHS %>% mutate(YDBHS=milkZ-mean) %>%
  group_by(cow) %>%
  summarise(YDBHS_cow=mean(YDBHS))

# Predict Breeding value of young cows
fitBHS_NA <- bru(modelBHS,
                like(family = "Gaussian",
                     modelBHSFormula,
                     data = data2_NA))

EBV_BHS_NA <- fitBHS_NA$summary.random$animal[,c("ID", "mean")]

colnames(EBV_BHS_NA)[1] <- "cow"
colnames(EBV_BHS_NA)[2] <- "EBV_BHS"

sel <- predictBHS$cow %in% young_cow_ID$cow
predictBHS <- predictBHS[sel,]
predictBHS<- subset(predictBHS,select= c("cow","YDBHS_cow"))
predictBHS <- data.frame(predictBHS)
predictBHS <- predictBHS[,c("cow", "YDBHS_cow")]

sel <- EBV_BHS_NA$cow %in% young_cow_ID$cow
EBV_BHS_NA_young<- EBV_BHS_NA[sel,]

EBV_YDBHS <- merge(predictBHS,EBV_BHS_NA_young, by="cow", all.x=TRUE)

# Forward Validation for accuracy for model fitBS
round(cor(EBV_YDBHS$YDBHS_cow,EBV_YDBHS$EBV_BHS),2) # -0.01

#------------Foward validation:Correcting for all non genetic factors-----------
#FitB
predictB2 <- predict(fitB, newdata=data2,
                    formula= ~ fixed_effects + perm + ward)
colnames(predictB2)
summary(data2$milkZ)
summary(predictB2$milkZ)
summary(predictB2$mean)
round(cor(predictB2$milkZ,predictB2$mean),2) # 0.85
predictB2 <- predictB2[,c("cow","milkZ","mean")]
predictB2 <- predictB2 %>% mutate(YDB=milkZ-mean) %>%
  group_by(cow) %>%
  summarise(YDB_cow=mean(YDB))

sel <- predictB2$cow %in% young_cow_ID$cow
predictB2 <- predictB2[sel,]
predictB2<- subset(predictB2,select= c("cow","YDB_cow"))
predictB2 <- data.frame(predictB2)
predictB2 <- predictB2[,c("cow","YDB_cow")]

sel <- EBV_B_NA$cow %in% young_cow_ID$cow
EBV_B_NA_young<- EBV_B_NA[sel,]

EBV_YDB <- merge(predictB2,EBV_B_NA_young, by="cow", all.x=TRUE)

# Forward Validation for accuracy for model fitB
round(cor(EBV_YDB$YDB_cow,EBV_YDB$EBV_B),2) # 0.26


#------ Forward Validation model fitBH------------------------------------------
predictBH2 <- predict(fitBH, newdata=data2,
                     formula= ~fixed_effects + perm + ward + herd)
colnames(predictBH2)
summary(data2$milkZ)
summary(predictBH2$milkZ)
summary(predictBH2$mean)
plot(predictBH2$mean,predictBH2$milkZ)
round(cor(predictBH2$mean,predictBH2$milkZ),2) # 0.85
predictBH2 <- predictBH2[,c("cow","milkZ","mean")]
predictBH2 <- predictBH2 %>% mutate(YDBH=milkZ-mean) %>%
  group_by(cow) %>%
  summarise(YDBH_cow=mean(YDBH))

sel <- predictBH2$cow %in% young_cow_ID$cow
predictBH2 <- predictBH2[sel,]
predictBH2<- subset(predictBH2,select= c("cow","YDBH_cow"))
predictBH2 <- data.frame(predictBH2)
predictBH2 <- predictBH2[,c("cow", "YDBH_cow")]

sel <- EBV_BH_NA$cow %in% young_cow_ID$cow
EBV_BH_NA_young<- EBV_BH_NA[sel,]

EBV_YDBH <- merge(predictBH2,EBV_BH_NA_young, by="cow", all.x=TRUE)

# Forward Validation for accuracy for model fitBH
round(cor(EBV_YDBH$YDBH_cow,EBV_YDBH$EBV_BH),2) # 0.20

#------ Forward Validation model fitBS------------------------------------------
predictBS2 <- predict(fitBS, newdata=data2,
                     formula= ~fixed_effects + perm + ward + field)

colnames(predictBS2)
summary(data2$milkZ)
summary(predictBS2$milkZ)
summary(predictBS2$mean)
plot(predictBS2$mean,predictBS2$milk)
round(cor(predictBS2$mean,predictBS2$milk),2) #0.86
predictBS2 <- predictBS2[,c("cow","milkZ","mean")]
predictBS2 <- predictBS2 %>% mutate(YDBS=milkZ-mean) %>%
  group_by(cow) %>%
  summarise(YDBS_cow=mean(YDBS))

sel <- predictBS2$cow %in% young_cow_ID$cow
predictBS2 <- predictBS2[sel,]
predictBS2<- subset(predictBS2,select= c("cow","YDBS_cow"))
predictBS2 <- data.frame(predictBS2)
predictBS2 <- predictBS2[,c("cow", "YDBS_cow")]

sel <- EBV_BS_NA$cow %in% young_cow_ID$cow
EBV_BS_NA_young<- EBV_BS_NA[sel,]

EBV_YDBS <- merge(predictBS2,EBV_BS_NA_young, by="cow", all.x=TRUE)

# Forward Validation for accuracy for model fitBS
round(cor(EBV_YDBS$YDBS_cow,EBV_YDBS$EBV_BS),2) # 0.06


#------ Forward Validation model fitBHS------------------------------------------
predictBHS2 <- predict(fitBHS, newdata=data2,
                      formula= ~fixed_effects + perm + ward + herd + field)

colnames(predictBHS2)
summary(data2$milkZ)
summary(predictBHS2$milkZ)
summary(predictBHS2$mean)
plot(predictBHS2$mean,predictBHS2$milk)
round(cor(predictBHS2$mean,predictBHS2$milk),2) #0.86
predictBHS2 <- predictBHS2[,c("cow","milkZ","mean")]
predictBHS2 <- predictBHS2 %>% mutate(YDBHS=milkZ-mean) %>%
  group_by(cow) %>%
  summarise(YDBHS_cow=mean(YDBHS))

sel <- predictBHS2$cow %in% young_cow_ID$cow
predictBHS2 <- predictBHS2[sel,]
predictBHS2<- subset(predictBHS2,select= c("cow","YDBHS_cow"))
predictBHS2 <- data.frame(predictBHS2)
predictBHS2 <- predictBHS2[,c("cow", "YDBHS_cow")]

sel <- EBV_BHS_NA$cow %in% young_cow_ID$cow
EBV_BHS_NA_young<- EBV_BHS_NA[sel,]

EBV_YDBHS <- merge(predictBHS2,EBV_BHS_NA_young, by="cow", all.x=TRUE)

# Forward Validation for accuracy for model fitBHS
round(cor(EBV_YDBHS$YDBHS_cow,EBV_YDBHS$EBV_BHS),2) # 0.01



#------------Foward validation:Predictive ability (EAAP)-----------------------
#FitB
predictB3 <- predict(fitB, newdata=data2_NA,
                     formula= ~fixed_effects + perm + ward + animal)
colnames(predictB3)
summary(data2_NA$milkZ)
summary(predictB3$milkZ)
summary(predictB3$mean)

sum(is.na(predictB3$mean)) # 0 expected
data2$milkZ_pred <- predictB3$mean
data2$milkZ_pred_sd <- predictB3$sd
predictB3 <-   subset(data2, birthyear=="2016" | birthyear=="2017")
predictB3<- subset(predictB3,select= c(cow,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(predictB3$cow))# 146 Cows expected to be born between 2016 (141 cows) and 2017 (5cows)
summary(predictB3)
summary(data2)

# Correlating observed phenotype of young cows  with predicted phenotype
accuracy_fitB <- round(cor(predictB3$milkZ,predictB3$milkZ_pred),2)
Coef<- lm (predictB3$milkZ~predictB3$milkZ_pred)
accuracy_fitB # 0.84
R2_fitBase<- summary(Coef)
R2_fitBase #
#CRPS
predictB3<- data.frame(predictB3)
obs <- predictB3$milkZ
pred<- subset(predictB3,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitB <- crps(obs,pred)
round((crps_fitB$CRPS),2) # 5.84

#FitBH
predictBH3 <- predict(fitBH, newdata=data2_NA,
                     formula= ~fixed_effects + perm + ward + herd + animal)
colnames(predictBH3)
summary(data2_NA$milkZ)
summary(predictBH3$milkZ)
summary(predictBH3$mean)

sum(is.na(predictBH3$mean)) # 0 expected
data2$milkZ_pred <- predictBH3$mean
data2$milkZ_pred_sd <- predictBH3$sd
predictBH3 <-   subset(data2, birthyear=="2016" | birthyear=="2017")
predictBH3<- subset(predictBH3,select= c(cow,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(predictBH3$cow))# 146 Cows expected to be born between 2016 (141 cows) and 2017 (5cows)
summary(predictBH3)
summary(data2)

# Correlating observed phenotype of young cows  with predicted phenotype
accuracy_fitBH <- round(cor(predictBH3$milkZ,predictBH3$milkZ_pred),2)
Coef<- lm (predictBH3$milkZ~predictBH3$milkZ_pred)
accuracy_fitBH # 0.84
R2_fitBH<- summary(Coef)
R2_fitBH #  0.70
#CRPS
predictBH3<- data.frame(predictBH3)
obs <- predictBH3$milkZ
pred<- subset(predictBH3,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBH <- crps(obs,pred)
round((crps_fitBH$CRPS),2) # 5.38


#FitBS
predictBS3 <- predict(fitBS, newdata=data2_NA,
                      formula= ~fixed_effects + perm + ward + field + animal)
colnames(predictBS3)
summary(data2_NA$milkZ)
summary(predictBS3$milkZ)
summary(predictBS3$mean)

sum(is.na(predictBS3$mean)) # 0 expected
data2$milkZ_pred <- predictBS3$mean
data2$milkZ_pred_sd <- predictBS3$sd
predictBS3 <-   subset(data2, birthyear=="2016" | birthyear=="2017")
predictBS3<- subset(predictBS3,select= c(cow,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(predictBS3$cow))# 146 Cows expected to be born between 2016 (141 cows) and 2017 (5cows)
summary(predictBS3)
summary(data2)

# Correlating observed phenotype of young cows  with predicted phenotype
accuracy_fitBS <- round(cor(predictBS3$milkZ,predictBS3$milkZ_pred),2)
Coef<- lm (predictBS3$milkZ~predictBS3$milkZ_pred)
accuracy_fitBS # 0.84
R2_fitBS<- summary(Coef)
R2_fitBS # 0.70
#CRPS
predictBS3<- data.frame(predictBS3)
obs <- predictBS3$milkZ
pred<- subset(predictBS3,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBS <- crps(obs,pred)
round((crps_fitBS$CRPS),2) # 5.77

#FitBHS
predictBHS3 <- predict(fitBHS, newdata=data2_NA,
                      formula= ~fixed_effects + perm + ward + herd + field + animal)
colnames(predictBHS3)
summary(data2_NA$milkZ)
summary(predictBHS3$milkZ)
summary(predictBHS3$mean)

sum(is.na(predictBHS3$mean)) # 0 expected
data2$milkZ_pred <- predictBHS3$mean
data2$milkZ_pred_sd <- predictBHS3$sd
predictBHS3 <-   subset(data2, birthyear=="2016" | birthyear=="2017")
predictBHS3<- subset(predictBHS3,select= c(cow,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(predictBHS3$cow))# 146 Cows expected to be born between 2016 (141 cows) and 2017 (5cows)
summary(predictBHS3)
summary(data2)

# Correlating observed phenotype of young cows  with predicted phenotype
accuracy_fitBHS <- round(cor(predictBHS3$milkZ,predictBHS3$milkZ_pred),2)
Coef<- lm (predictBHS3$milkZ~predictBHS3$milkZ_pred)
accuracy_fitBHS # 0.84
R2_fitBHS<- summary(Coef)
R2_fitBHS #0.7
#CRPS
predictBHS3<- data.frame(predictBHS3)
obs <- predictBHS3$milkZ
pred<- subset(predictBS3,select= c(milkZ_pred,milkZ_pred_sd)) # 0.84
crps_fitBHS <- crps(obs,pred)
round((crps_fitBHS$CRPS),2) # 5.77

#------------------------------Cross-validation_Corrected Phenotype---------------------------------
#-----------------Cross-validation fitB-----------------------------------------
# Masking phenotypes of cows in dgrp 1 (Breed proportion class 1)
#dgrp1 masked
data2_1_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "1", NA, milkZ))
sum(is.na(data2$milkZ)) # 0
sum(is.na(data2_1_NA$milkZ)) #7466
length(unique(data2_1_NA$cow)) # 1894
length(unique(data2$cow)) # 1894

predictB <- predict(fitB, newdata=data2,
                    formula= ~ fixed_effects + perm)
colnames(predictB)
summary(data2$milkZ)
summary(predictB$milkZ)
summary(predictB$mean)
round(cor(predictB$milkZ,predictB$mean),2) # 0.77
predictB <- predictB[,c("cow","milkZ","mean")]
predictB <- predictB %>% mutate(YDB=milkZ-mean) %>%
  group_by(cow) %>%
  summarise(YDB_cow=mean(YDB))

predictB2 <- predictB # Saving predictB

# Breeding values and Yield deviation model fitB

# Predict Breeding value of cows in dgrp1
fitB_1_NA <- bru(modelB,
               like(family = "Gaussian",
                    modelBFormula,
                    data = data2_1_NA))

EBV_B_1_NA <- fitB_1_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_B_1_NA)[1] <- "cow"
colnames(EBV_B_1_NA)[2] <- "EBV_B"

# data2 for dgrp1
data2_1 <- subset(data2, dgrp==1)
length(unique(data2_1$cow)) # 731 cows
# cows with dgrp1 IDS  (dgrp=1)
cow_ID_1 <- subset(data2_1, select=c(cow))
cow_ID_1 <- unique(cow_ID_1)
summary(cow_ID_1)
rownames(cow_ID_1) <- NULL
length(unique(cow_ID_1$cow)) # 731 cows
sel <- predictB$cow %in% cow_ID_1$cow
predictB <- predictB[sel,]
predictB<- subset(predictB,select= c("cow","YDB_cow"))
predictB <- data.frame(predictB)
predictB <- predictB[,c("cow","YDB_cow")]

sel <- EBV_B_1_NA$cow %in% cow_ID_1$cow
EBV_B_1_NA <- EBV_B_1_NA[sel,]

EBV_YDB <- merge(predictB,EBV_B_1_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitB for cows in dgrp1
round(cor(EBV_YDB$YDB_cow,EBV_YDB$EBV_B),2) # 0.05

#dgrp2 masked
data2_2_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "2", NA, milkZ))
sum(is.na(data2$milkZ)) # 0
sum(is.na(data2_2_NA$milkZ)) #8149
length(unique(data2_2_NA$cow)) # 1894 cows
length(unique(data2$cow)) # 1894

# Predict Breeding value of cows in dgrp2
fitB_2_NA <- bru(modelB,
               like(family = "Gaussian",
                    modelBFormula,
                    data = data2_2_NA))

EBV_B_2_NA <- fitB_2_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_B_2_NA)[1] <- "cow"
colnames(EBV_B_2_NA)[2] <- "EBV_B"

# data2 for dgrp2
data2_2 <- subset(data2, dgrp==2)
length(unique(data2_2$cow)) #  cows
# cows with dgrp2 IDS  (dgrp=2)
cow_ID_2 <- subset(data2_2, select=c(cow))
cow_ID_2 <- unique(cow_ID_2)
summary(cow_ID_2)
rownames(cow_ID_2) <- NULL
length(unique(cow_ID_2$cow)) # 770 cows
predictB<- predictB2 # Recalling predictB
sel <- predictB$cow %in% cow_ID_2$cow
predictB <- predictB[sel,]
predictB<- subset(predictB,select= c("cow","YDB_cow"))
predictB <- data.frame(predictB)
predictB <- predictB[,c("cow","YDB_cow")]

sel <- EBV_B_2_NA$cow %in% cow_ID_2$cow
EBV_B_2_NA <- EBV_B_2_NA[sel,]

EBV_YDB <- merge(predictB,EBV_B_2_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitB for cows in dgrp2
round(cor(EBV_YDB$YDB_cow,EBV_YDB$EBV_B),2) # 0. 14

#dgrp3 masked
data2_3_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "3", NA, milkZ))
sum(is.na(data2$milkZ)) # 0
sum(is.na(data2_3_NA$milkZ)) #3058
length(unique(data2_3_NA$cow)) # 1894cows
length(unique(data2$cow)) # 1894

# Predict Breeding value of cows in dgrp2
fitB_3_NA <- bru(modelB,
                 like(family = "Gaussian",
                      modelBFormula,
                      data = data2_3_NA))

EBV_B_3_NA <- fitB_3_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_B_3_NA)[1] <- "cow"
colnames(EBV_B_3_NA)[2] <- "EBV_B"

# data2 for dgrp3
data2_3 <- subset(data2, dgrp==3)
length(unique(data2_3$cow)) # 309 cows
# cows with dgrp2 IDS  (dgrp=3)
cow_ID_3 <- subset(data2_3, select=c(cow))
cow_ID_3 <- unique(cow_ID_3)
summary(cow_ID_3)
rownames(cow_ID_3) <- NULL
length(unique(cow_ID_3$cow)) # 309 cows
predictB<- predictB2 # Recalling predictB
sel <- predictB$cow %in% cow_ID_3$cow
predictB <- predictB[sel,]
predictB<- subset(predictB,select= c("cow","YDB_cow"))
predictB <- data.frame(predictB)
predictB <- predictB[,c("cow","YDB_cow")]

sel <- EBV_B_3_NA$cow %in% cow_ID_3$cow
EBV_B_3_NA <- EBV_B_3_NA[sel,]

EBV_YDB <- merge(predictB,EBV_B_3_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitB for cows in dgrp3
round(cor(EBV_YDB$YDB_cow,EBV_YDB$EBV_B),2) # 0.14


#dgrp4 masked
data2_4_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "4", NA, milkZ))
sum(is.na(data2$milkZ)) # 0
sum(is.na(data2_4_NA$milkZ)) #702
length(unique(data2_4_NA$cow)) # 1894cows
length(unique(data2$cow)) # 1894

# Predict Breeding value of cows in dgrp2
fitB_4_NA <- bru(modelB,
                 like(family = "Gaussian",
                      modelBFormula,
                      data = data2_4_NA))

EBV_B_4_NA <- fitB_4_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_B_4_NA)[1] <- "cow"
colnames(EBV_B_4_NA)[2] <- "EBV_B"

# data2 for dgrp4
data2_4 <- subset(data2, dgrp==4)
length(unique(data2_4$cow)) #  84 cows
# cows with dgrp4 IDS  (dgrp=4)
cow_ID_4 <- subset(data2_4, select=c(cow))
cow_ID_4 <- unique(cow_ID_4)
summary(cow_ID_4)
rownames(cow_ID_4) <- NULL
length(unique(cow_ID_4$cow)) #  84cows
predictB<- predictB2 # Recalling predictB
sel <- predictB$cow %in% cow_ID_4$cow
predictB <- predictB[sel,]
predictB<- subset(predictB,select= c("cow","YDB_cow"))
predictB <- data.frame(predictB)
predictB <- predictB[,c("cow","YDB_cow")]

sel <- EBV_B_4_NA$cow %in% cow_ID_4$cow
EBV_B_4_NA <- EBV_B_4_NA[sel,]

EBV_YDB <- merge(predictB,EBV_B_4_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitB for cows in dgrp4
round(cor(EBV_YDB$YDB_cow,EBV_YDB$EBV_B),2) # -0.04

# The Average accuracy for model fitB
Acc_B <- round((0.05+0.14 + 0.14-0.04)/4,2)
Acc_B # 0.07

#------------Cross-validation fitBH---------------------------------------------
# Masking phenotypes of cows in dgrp 1 (Breed proportion class 1)
#dgrp1 masked

predictBH <- predict(fitBH, newdata=data2,
                    formula= ~ fixed_effects + perm)
colnames(predictBH)
summary(data2$milkZ)
summary(predictBH$milkZ)
summary(predictBH$mean)
round(cor(predictBH$milkZ,predictBH$mean),2) # 0.65
predictBH <- predictBH[,c("cow","milkZ","mean")]
predictBH <- predictBH %>% mutate(YDBH=milkZ-mean) %>%
  group_by(cow) %>%
  summarise(YDBH_cow=mean(YDBH))

predictBH2 <- predictBH # Saving predictBH

# Breeding values and Yield deviation model fitBH

# Predict Breeding value of cows in dgrp1
fitBH_1_NA <- bru(modelBH,
                 like(family = "Gaussian",
                      modelBHFormula,
                      data = data2_1_NA))

EBV_BH_1_NA <- fitBH_1_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_BH_1_NA)[1] <- "cow"
colnames(EBV_BH_1_NA)[2] <- "EBV_BH"
length(unique(cow_ID_1$cow)) # 731 cows
sel <- predictBH$cow %in% cow_ID_1$cow
predictBH <- predictBH[sel,]
predictBH<- subset(predictBH,select= c("cow","YDBH_cow"))
predictBH <- data.frame(predictBH)
predictBH <- predictBH[,c("cow","YDBH_cow")]

sel <- EBV_BH_1_NA$cow %in% cow_ID_1$cow
EBV_BH_1_NA <- EBV_BH_1_NA[sel,]

EBV_YDBH <- merge(predictBH,EBV_BH_1_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitBH for cows in dgrp1
round(cor(EBV_YDBH$YDBH_cow,EBV_YDBH$EBV_BH),2) # 0.03

#dgrp2 masked

# Predict Breeding value of cows in dgrp2
fitBH_2_NA <- bru(modelBH,
                 like(family = "Gaussian",
                      modelBHFormula,
                      data = data2_2_NA))

EBV_BH_2_NA <- fitBH_2_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_BH_2_NA)[1] <- "cow"
colnames(EBV_BH_2_NA)[2] <- "EBV_BH"

# data2 for dgrp2
data2_2 <- subset(data2, dgrp==2)
length(unique(data2_2$cow)) #  cows
# cows with dgrp2 IDS  (dgrp=2)
cow_ID_2 <- subset(data2_2, select=c(cow))
cow_ID_2 <- unique(cow_ID_2)
summary(cow_ID_2)
rownames(cow_ID_2) <- NULL
length(unique(cow_ID_2$cow)) # 770 cows
predictBH<- predictBH2 # Recalling predictBH
sel <- predictBH$cow %in% cow_ID_2$cow
predictBH <- predictBH[sel,]
predictBH<- subset(predictBH,select= c("cow","YDBH_cow"))
predictBH <- data.frame(predictBH)
predictBH <- predictBH[,c("cow","YDBH_cow")]

sel <- EBV_BH_2_NA$cow %in% cow_ID_2$cow
EBV_BH_2_NA <- EBV_BH_2_NA[sel,]

EBV_YDBH <- merge(predictBH,EBV_BH_2_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitB for cows in dgrp2
round(cor(EBV_YDBH$YDBH_cow,EBV_YDBH$EBV_BH),2) # 0.07

#dgrp3 masked

# Predict Breeding value of cows in dgrp2
fitBH_3_NA <- bru(modelBH,
                 like(family = "Gaussian",
                      modelBHFormula,
                      data = data2_3_NA))

EBV_BH_3_NA <- fitBH_3_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_BH_3_NA)[1] <- "cow"
colnames(EBV_BH_3_NA)[2] <- "EBV_BH"

# data2 for dgrp3
data2_3 <- subset(data2, dgrp==3)
length(unique(data2_3$cow)) #  cows
# cows with dgrp3 IDS  (dgrp=3)
cow_ID_3 <- subset(data2_3, select=c(cow))
cow_ID_3 <- unique(cow_ID_3)
summary(cow_ID_3)
rownames(cow_ID_3) <- NULL
length(unique(cow_ID_3$cow)) # 770 cows
predictBH<- predictBH2 # Recalling predictBH
sel <- predictBH$cow %in% cow_ID_3$cow
predictBH <- predictBH[sel,]
predictBH<- subset(predictBH,select= c("cow","YDBH_cow"))
predictBH <- data.frame(predictBH)
predictBH <- predictBH[,c("cow","YDBH_cow")]

sel <- EBV_BH_3_NA$cow %in% cow_ID_3$cow
EBV_BH_3_NA <- EBV_BH_3_NA[sel,]

EBV_YDBH <- merge(predictBH,EBV_BH_3_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitB for cows in dgrp3
round(cor(EBV_YDBH$YDBH_cow,EBV_YDBH$EBV_BH),2) # 0.08


#dgrp4 masked

# Predict Breeding value of cows in dgrp4
fitBH_4_NA <- bru(modelBH,
                 like(family = "Gaussian",
                      modelBHFormula,
                      data = data2_4_NA))

EBV_BH_4_NA <- fitBH_4_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_BH_4_NA)[1] <- "cow"
colnames(EBV_BH_4_NA)[2] <- "EBV_BH"

# data2 for dgrp4
data2_4 <- subset(data2, dgrp==4)
length(unique(data2_4$cow)) #  84 cows
# cows with dgrp4 IDS  (dgrp=4)
cow_ID_4 <- subset(data2_4, select=c(cow))
cow_ID_4 <- unique(cow_ID_4)
summary(cow_ID_4)
rownames(cow_ID_4) <- NULL
length(unique(cow_ID_4$cow)) #  cows
predictBH<- predictBH2 # Recalling predictB
sel <- predictBH$cow %in% cow_ID_4$cow
predictBH <- predictBH[sel,]
predictBH<- subset(predictBH,select= c("cow","YDBH_cow"))
predictBH <- data.frame(predictBH)
predictBH <- predictBH[,c("cow","YDBH_cow")]

sel <- EBV_BH_4_NA$cow %in% cow_ID_4$cow
EBV_BH_4_NA <- EBV_BH_4_NA[sel,]

EBV_YDBH <- merge(predictBH,EBV_BH_4_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitB for cows in dgrp4
round(cor(EBV_YDBH$YDBH_cow,EBV_YDBH$EBV_BH),2) # 0.01

# The Average accuracy for model fitBH
Acc_BH <- round((0.03 + 0.07+ 0.08+0.01)/4,2)
Acc_BH # 0.05


#------------Cross-validation fitBS---------------------------------------------
# Masking phenotypes of cows in dgrp 1 (Breed proportion class 1)
#dgrp1 masked

predictBS <- predict(fitBS, newdata=data2,
                     formula= ~ fixed_effects + perm)
colnames(predictBS)
summary(data2$milkZ)
summary(predictBS$milkZ)
summary(predictBS$mean)
round(cor(predictBS$milkZ,predictBS$mean),2) # 0.65
predictBS <- predictBS[,c("cow","milkZ","mean")]
predictBS <- predictBS %>% mutate(YDBS=milkZ-mean) %>%
  group_by(cow) %>%
  summarise(YDBS_cow=mean(YDBS))

predictBS2 <- predictBS # Saving predictBS

# Breeding values and Yield deviation model fitBS

# Predict Breeding value of cows in dgrp1
fitBS_1_NA <- bru(modelBS,
                  like(family = "Gaussian",
                       modelBSFormula,
                       data = data2_1_NA))

EBV_BS_1_NA <- fitBS_1_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_BS_1_NA)[1] <- "cow"
colnames(EBV_BS_1_NA)[2] <- "EBV_BS"
length(unique(cow_ID_1$cow)) # 731 cows
sel <- predictBS$cow %in% cow_ID_1$cow
predictBS <- predictBS[sel,]
predictBS<- subset(predictBS,select= c("cow","YDBS_cow"))
predictBS <- data.frame(predictBS)
predictBS <- predictBS[,c("cow","YDBS_cow")]

sel <- EBV_BS_1_NA$cow %in% cow_ID_1$cow
EBV_BS_1_NA <- EBV_BS_1_NA[sel,]

EBV_YDBS <- merge(predictBS,EBV_BS_1_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitBS for cows in dgrp1
round(cor(EBV_YDBS$YDBS_cow,EBV_YDBS$EBV_BS),2) # 0.07

#dgrp2 masked

# Predict Breeding value of cows in dgrp2
fitBS_2_NA <- bru(modelBS,
                  like(family = "Gaussian",
                       modelBSFormula,
                       data = data2_2_NA))

EBV_BS_2_NA <- fitBS_2_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_BS_2_NA)[1] <- "cow"
colnames(EBV_BS_2_NA)[2] <- "EBV_BS"

# data2 for dgrp2
data2_2 <- subset(data2, dgrp==2)
length(unique(data2_2$cow)) #  cows
# cows with dgrp2 IDS  (dgrp=2)
cow_ID_2 <- subset(data2_2, select=c(cow))
cow_ID_2 <- unique(cow_ID_2)
summary(cow_ID_2)
rownames(cow_ID_2) <- NULL
length(unique(cow_ID_2$cow)) # 770 cows
predictBS<- predictBS2 # Recalling predictBS
sel <- predictBS$cow %in% cow_ID_2$cow
predictBS <- predictBS[sel,]
predictBS<- subset(predictBS,select= c("cow","YDBS_cow"))
predictBS <- data.frame(predictBS)
predictBS <- predictBS[,c("cow","YDBS_cow")]

sel <- EBV_BS_2_NA$cow %in% cow_ID_2$cow
EBV_BS_2_NA <- EBV_BS_2_NA[sel,]

EBV_YDBS <- merge(predictBS,EBV_BS_2_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitB for cows in dgrp2
round(cor(EBV_YDBS$YDBS_cow,EBV_YDBS$EBV_BS),2) # 0.07

#dgrp3 masked

# Predict Breeding value of cows in dgrp2
fitBS_3_NA <- bru(modelBS,
                  like(family = "Gaussian",
                       modelBSFormula,
                       data = data2_3_NA))

EBV_BS_3_NA <- fitBS_3_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_BS_3_NA)[1] <- "cow"
colnames(EBV_BS_3_NA)[2] <- "EBV_BS"

# data2 for dgrp3
data2_3 <- subset(data2, dgrp==3)
length(unique(data2_3$cow)) #  cows
# cows with dgrp3 IDS  (dgrp=3)
cow_ID_3 <- subset(data2_3, select=c(cow))
cow_ID_3 <- unique(cow_ID_3)
summary(cow_ID_3)
rownames(cow_ID_3) <- NULL
length(unique(cow_ID_3$cow)) # 770 cows
predictBS<- predictBS2 # Recalling predictBS
sel <- predictBS$cow %in% cow_ID_3$cow
predictBS <- predictBS[sel,]
predictBS<- subset(predictBS,select= c("cow","YDBS_cow"))
predictBS <- data.frame(predictBS)
predictBS <- predictBS[,c("cow","YDBS_cow")]

sel <- EBV_BS_3_NA$cow %in% cow_ID_3$cow
EBV_BS_3_NA <- EBV_BS_3_NA[sel,]

EBV_YDBS <- merge(predictBS,EBV_BS_3_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitB for cows in dgrp3
round(cor(EBV_YDBS$YDBS_cow,EBV_YDBS$EBV_BS),2) # 0.08


#dgrp4 masked

# Predict Breeding value of cows in dgrp4
fitBS_4_NA <- bru(modelBS,
                  like(family = "Gaussian",
                       modelBSFormula,
                       data = data2_4_NA))

EBV_BS_4_NA <- fitBS_4_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_BS_4_NA)[1] <- "cow"
colnames(EBV_BS_4_NA)[2] <- "EBV_BS"

# data2 for dgrp4
data2_4 <- subset(data2, dgrp==4)
length(unique(data2_4$cow)) #  84 cows
# cows with dgrp4 IDS  (dgrp=4)
cow_ID_4 <- subset(data2_4, select=c(cow))
cow_ID_4 <- unique(cow_ID_4)
summary(cow_ID_4)
rownames(cow_ID_4) <- NULL
length(unique(cow_ID_4$cow)) #  cows
predictBS<- predictBS2 # Recalling predictB
sel <- predictBS$cow %in% cow_ID_4$cow
predictBS <- predictBS[sel,]
predictBS<- subset(predictBS,select= c("cow","YDBS_cow"))
predictBS <- data.frame(predictBS)
predictBS <- predictBS[,c("cow","YDBS_cow")]

sel <- EBV_BS_4_NA$cow %in% cow_ID_4$cow
EBV_BS_4_NA <- EBV_BS_4_NA[sel,]

EBV_YDBS <- merge(predictBS,EBV_BS_4_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitB for cows in dgrp4
round(cor(EBV_YDBS$YDBS_cow,EBV_YDBS$EBV_BS),2) # 0.01

# The Average accuracy for model fitBS
Acc_BS <- round(( + + +)/4,2)
Acc_BS #

#------------------------------Cross-validation_Predicted Pheno Vs Obs_Pheno----
#-----------------Cross-validation fitB-----------------------------------------
# Masking phenotypes of cows in dgrp (Breed proportion class)
#dgrp1 masked
#Let's create unique ID (merge_ID) columns for data merging
data2_merge_ID <- data2
data2_merge_ID$merge_ID <- 1:nrow(data2_merge_ID)
data2_1_NA<- data2_merge_ID %>% mutate(milkZ = ifelse(dgrp == "1", NA, milkZ))
sum(is.na(data2$milkZ)) # 0
sum(is.na(data2_1_NA$milkZ)) #7466
length(unique(data2_1_NA$cow)) # 1894
length(unique(data2$cow)) # 1894

predictB_2_1_NA <- predict(fitB, newdata=data2_1_NA,
                    formula= ~ fixed_effects + perm + ward + animal)
colnames(predictB_2_1_NA)
pheno1 <- data.frame(predictB_2_1_NA)
pheno1 <- pheno1[,c("merge_ID","mean")]
dim(pheno1)
class(pheno1)
data2_pheno1 <- data.frame(data2_merge_ID)
data2_pheno1 <- data2_pheno1[,c("merge_ID","milkZ","dgrp")]
summary(data2_pheno1)
colnames(data2_pheno1)
summary(data2_merge_ID$milkZ)
class(data2_pheno1)
class(pheno1)
dim(data2_pheno1)
dim(pheno1)
length(unique(data2_pheno1$merge_ID))
length(unique(pheno1$merge_ID))

pheno1_final <- merge(data2_pheno1,pheno1, by="merge_ID", all.x=TRUE)

dim(pheno1_final)
Pheno_1_NA <- subset(pheno1_final, dgrp=="1")
summary(Pheno_1_NA$dgrp)
Accu1 <- round(cor(Pheno_1_NA$milkZ,Pheno_1_NA$mean),2) # 0.83
Accu1 # 0.83

#dgrp2 masked
data2_2_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "2", NA, milkZ))
sum(is.na(data2$milkZ)) # 0
sum(is.na(data2_2_NA$milkZ)) #8149
length(unique(data2_2_NA$cow)) # 1894
length(unique(data2$cow)) # 1894

predictB_2_2_NA <- predict(fitB, newdata=data2_2_NA,
                           formula= ~ fixed_effects + perm + ward + animal)
colnames(predictB_2_2_NA)

summary(data2$milkZ)
summary(predictB_2_2_NA$mean)
predictB_2_2_NA$milkZ <- data2$milkZ

Pheno_2_NA <- subset(predictB_2_2_NA, dgrp=="2")
summary(Pheno_2_NA$dgrp)
Accu2 <- round(cor(Pheno_2_NA$milkZ,Pheno_2_NA$mean),2)
Accu2 # 0.85

#dgrp3 masked
data2_3_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "3", NA, milkZ))
sum(is.na(data2$milkZ)) # 0
sum(is.na(data2_3_NA$milkZ)) # 3058
length(unique(data2_3_NA$cow)) # 1894
length(unique(data2$cow)) # 1894

predictB_2_3_NA <- predict(fitB, newdata=data2_3_NA,
                           formula= ~ fixed_effects + perm + ward + animal)
colnames(predictB_2_3_NA)

summary(data2$milkZ)
summary(predictB_2_3_NA$mean)
predictB_2_3_NA$milkZ <- data2$milkZ

Pheno_3_NA <- subset(predictB_2_3_NA, dgrp=="3")
summary(Pheno_3_NA$dgrp)
Accu3 <- round(cor(Pheno_3_NA$milkZ,Pheno_3_NA$mean),2) # 0.
Accu3 # 0.83


#dgrp4 masked
data2_4_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "4", NA, milkZ))
sum(is.na(data2$milkZ)) # 0
sum(is.na(data2_4_NA$milkZ)) # 702
length(unique(data2_4_NA$cow)) # 1894
length(unique(data2$cow)) # 1894

predictB_2_4_NA <- predict(fitB, newdata=data2_4_NA,
                           formula= ~ fixed_effects + perm + ward + animal)
colnames(predictB_2_4_NA)

summary(data2$milkZ)
summary(predictB_2_4_NA$mean)
predictB_2_4_NA$milkZ <- data2$milkZ

Pheno_4_NA <- subset(predictB_2_4_NA, dgrp=="4")
summary(Pheno_4_NA$dgrp)
Accu4 <- round(cor(Pheno_4_NA$milkZ,Pheno_4_NA$mean),2) # 0.
Accu4 # 0.82

Accu_B <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_B # 0.83


#-----------------Cross-validation fitBH-----------------------------------------
# Masking phenotypes of cows in dgrp (Breed proportion class)
#dgrp1 masked
predictBH_2_1_NA <- predict(fitBH, newdata=data2_1_NA,
                           formula= ~ fixed_effects + perm + ward + animal + herd)
colnames(predictBH_2_1_NA)

summary(data2$milkZ)
summary(predictBH_2_1_NA$mean)
predictBH_2_1_NA$milkZ <- data2$milkZ

Pheno_1_NA <- subset(predictBH_2_1_NA, dgrp=="1")
summary(Pheno_1_NA$dgrp)
Accu1 <- round(cor(Pheno_1_NA$milkZ,Pheno_1_NA$mean),2) #
Accu1 # 0.83

#dgrp2 masked
data2_2_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "2", NA, milkZ))
sum(is.na(data2$milkZ)) # 0
sum(is.na(data2_2_NA$milkZ)) #8149
length(unique(data2_2_NA$cow)) # 1894
length(unique(data2$cow)) # 1894

predictBH_2_2_NA <- predict(fitBH, newdata=data2_2_NA,
                           formula= ~ fixed_effects + perm + ward + animal + herd)
colnames(predictBH_2_2_NA)

summary(data2$milkZ)
summary(predictBH_2_2_NA$mean)
predictBH_2_2_NA$milkZ <- data2$milkZ

Pheno_2_NA <- subset(predictBH_2_2_NA, dgrp=="2")
summary(Pheno_2_NA$dgrp)
Accu2 <- round(cor(Pheno_2_NA$milkZ,Pheno_2_NA$mean),2)
Accu2 # 0.85

#dgrp3 masked

predictBH_2_3_NA <- predict(fitBH, newdata=data2_3_NA,
                           formula= ~ fixed_effects + perm + ward + animal+herd)
colnames(predictBH_2_3_NA)

summary(data2$milkZ)
summary(predictBH_2_3_NA$mean)
predictBH_2_3_NA$milkZ <- data2$milkZ

Pheno_3_NA <- subset(predictBH_2_3_NA, dgrp=="3")
summary(Pheno_3_NA$dgrp)
Accu3 <- round(cor(Pheno_3_NA$milkZ,Pheno_3_NA$mean),2)
Accu3 # 0.83


#dgrp4 masked
predictBH_2_4_NA <- predict(fitBH, newdata=data2_4_NA,
                           formula= ~ fixed_effects + perm + ward + animal +herd)
colnames(predictBH_2_4_NA)

summary(data2$milkZ)
summary(predictBH_2_4_NA$mean)
predictBH_2_4_NA$milkZ <- data2$milkZ

Pheno_4_NA <- subset(predictBH_2_4_NA, dgrp=="4")
summary(Pheno_4_NA$dgrp)
Accu4 <- round(cor(Pheno_4_NA$milkZ,Pheno_4_NA$mean),2) # 0.
Accu4 # 0.82

Accu_BH <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_BH # 0.83


#-----------------Cross-validation fitBS-----------------------------------------
# Masking phenotypes of cows in dgrp (Breed proportion class)
#dgrp1 masked
predictBS_2_1_NA <- predict(fitBS, newdata=data2_1_NA,
                            formula= ~ fixed_effects + perm + ward + animal + field)
colnames(predictBS_2_1_NA)

summary(data2$milkZ)
summary(predictBS_2_1_NA$mean)
predictBS_2_1_NA$milkZ <- data2$milkZ

Pheno_1_NA <- subset(predictBS_2_1_NA, dgrp=="1")
summary(Pheno_1_NA$dgrp)
Accu1 <- round(cor(Pheno_1_NA$milkZ,Pheno_1_NA$mean),2) #
Accu1 #

#dgrp2 masked
data2_2_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "2", NA, milkZ))
sum(is.na(data2$milkZ)) # 0
sum(is.na(data2_2_NA$milkZ)) #8149
length(unique(data2_2_NA$cow)) # 1894
length(unique(data2$cow)) # 1894

predictBH_2_2_NA <- predict(fitBH, newdata=data2_2_NA,
                            formula= ~ fixed_effects + perm + ward + animal + herd)
colnames(predictBH_2_2_NA)

summary(data2$milkZ)
summary(predictBH_2_2_NA$mean)
predictBH_2_2_NA$milkZ <- data2$milkZ

Pheno_2_NA <- subset(predictBH_2_2_NA, dgrp=="2")
summary(Pheno_2_NA$dgrp)
Accu2 <- round(cor(Pheno_2_NA$milkZ,Pheno_2_NA$mean),2)
Accu2 # 0.85

#dgrp3 masked

predictBH_2_3_NA <- predict(fitBH, newdata=data2_3_NA,
                            formula= ~ fixed_effects + perm + ward + animal+herd)
colnames(predictBH_2_3_NA)

summary(data2$milkZ)
summary(predictBH_2_3_NA$mean)
predictBH_2_3_NA$milkZ <- data2$milkZ

Pheno_3_NA <- subset(predictBH_2_3_NA, dgrp=="3")
summary(Pheno_3_NA$dgrp)
Accu3 <- round(cor(Pheno_3_NA$milkZ,Pheno_3_NA$mean),2)
Accu3 # 0.83


#dgrp4 masked
predictBH_2_4_NA <- predict(fitBH, newdata=data2_4_NA,
                            formula= ~ fixed_effects + perm + ward + animal +herd)
colnames(predictBH_2_4_NA)

summary(data2$milkZ)
summary(predictBH_2_4_NA$mean)
predictBH_2_4_NA$milkZ <- data2$milkZ

Pheno_4_NA <- subset(predictBH_2_4_NA, dgrp=="4")
summary(Pheno_4_NA$dgrp)
Accu4 <- round(cor(Pheno_4_NA$milkZ,Pheno_4_NA$mean),2) # 0.
Accu4 # 0.82

Accu_BH <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_BH # 0.83
























predictB2 <- predictB # Saving predictB

# Breeding values and Yield deviation model fitB

# Predict Breeding value of cows in dgrp1
fitB_1_NA <- bru(modelB,
                 like(family = "Gaussian",
                      modelBFormula,
                      data = data2_1_NA))

EBV_B_1_NA <- fitB_1_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_B_1_NA)[1] <- "cow"
colnames(EBV_B_1_NA)[2] <- "EBV_B"

# data2 for dgrp1
data2_1 <- subset(data2, dgrp==1)
length(unique(data2_1$cow)) # 731 cows
# cows with dgrp1 IDS  (dgrp=1)
cow_ID_1 <- subset(data2_1, select=c(cow))
cow_ID_1 <- unique(cow_ID_1)
summary(cow_ID_1)
rownames(cow_ID_1) <- NULL
length(unique(cow_ID_1$cow)) # 731 cows
sel <- predictB$cow %in% cow_ID_1$cow
predictB <- predictB[sel,]
predictB<- subset(predictB,select= c("cow","YDB_cow"))
predictB <- data.frame(predictB)
predictB <- predictB[,c("cow","YDB_cow")]

sel <- EBV_B_1_NA$cow %in% cow_ID_1$cow
EBV_B_1_NA <- EBV_B_1_NA[sel,]

EBV_YDB <- merge(predictB,EBV_B_1_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitB for cows in dgrp1
round(cor(EBV_YDB$YDB_cow,EBV_YDB$EBV_B),2) # 0.05

#dgrp2 masked
data2_2_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "2", NA, milkZ))
sum(is.na(data2$milkZ)) # 0
sum(is.na(data2_2_NA$milkZ)) #8149
length(unique(data2_2_NA$cow)) # 1894 cows
length(unique(data2$cow)) # 1894

# Predict Breeding value of cows in dgrp2
fitB_2_NA <- bru(modelB,
                 like(family = "Gaussian",
                      modelBFormula,
                      data = data2_2_NA))

EBV_B_2_NA <- fitB_2_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_B_2_NA)[1] <- "cow"
colnames(EBV_B_2_NA)[2] <- "EBV_B"

# data2 for dgrp2
data2_2 <- subset(data2, dgrp==2)
length(unique(data2_2$cow)) #  cows
# cows with dgrp2 IDS  (dgrp=2)
cow_ID_2 <- subset(data2_2, select=c(cow))
cow_ID_2 <- unique(cow_ID_2)
summary(cow_ID_2)
rownames(cow_ID_2) <- NULL
length(unique(cow_ID_2$cow)) # 770 cows
predictB<- predictB2 # Recalling predictB
sel <- predictB$cow %in% cow_ID_2$cow
predictB <- predictB[sel,]
predictB<- subset(predictB,select= c("cow","YDB_cow"))
predictB <- data.frame(predictB)
predictB <- predictB[,c("cow","YDB_cow")]

sel <- EBV_B_2_NA$cow %in% cow_ID_2$cow
EBV_B_2_NA <- EBV_B_2_NA[sel,]

EBV_YDB <- merge(predictB,EBV_B_2_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitB for cows in dgrp2
round(cor(EBV_YDB$YDB_cow,EBV_YDB$EBV_B),2) # 0. 14

#dgrp3 masked
data2_3_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "3", NA, milkZ))
sum(is.na(data2$milkZ)) # 0
sum(is.na(data2_3_NA$milkZ)) #3058
length(unique(data2_3_NA$cow)) # 1894cows
length(unique(data2$cow)) # 1894

# Predict Breeding value of cows in dgrp2
fitB_3_NA <- bru(modelB,
                 like(family = "Gaussian",
                      modelBFormula,
                      data = data2_3_NA))

EBV_B_3_NA <- fitB_3_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_B_3_NA)[1] <- "cow"
colnames(EBV_B_3_NA)[2] <- "EBV_B"

# data2 for dgrp3
data2_3 <- subset(data2, dgrp==3)
length(unique(data2_3$cow)) # 309 cows
# cows with dgrp2 IDS  (dgrp=3)
cow_ID_3 <- subset(data2_3, select=c(cow))
cow_ID_3 <- unique(cow_ID_3)
summary(cow_ID_3)
rownames(cow_ID_3) <- NULL
length(unique(cow_ID_3$cow)) # 309 cows
predictB<- predictB2 # Recalling predictB
sel <- predictB$cow %in% cow_ID_3$cow
predictB <- predictB[sel,]
predictB<- subset(predictB,select= c("cow","YDB_cow"))
predictB <- data.frame(predictB)
predictB <- predictB[,c("cow","YDB_cow")]

sel <- EBV_B_3_NA$cow %in% cow_ID_3$cow
EBV_B_3_NA <- EBV_B_3_NA[sel,]

EBV_YDB <- merge(predictB,EBV_B_3_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitB for cows in dgrp3
round(cor(EBV_YDB$YDB_cow,EBV_YDB$EBV_B),2) # 0.14


#dgrp4 masked
data2_4_NA<- data2 %>% mutate(milkZ = ifelse(dgrp == "4", NA, milkZ))
sum(is.na(data2$milkZ)) # 0
sum(is.na(data2_4_NA$milkZ)) #702
length(unique(data2_4_NA$cow)) # 1894cows
length(unique(data2$cow)) # 1894

# Predict Breeding value of cows in dgrp2
fitB_4_NA <- bru(modelB,
                 like(family = "Gaussian",
                      modelBFormula,
                      data = data2_4_NA))

EBV_B_4_NA <- fitB_4_NA$summary.random$animal[,c("ID","mean")]

colnames(EBV_B_4_NA)[1] <- "cow"
colnames(EBV_B_4_NA)[2] <- "EBV_B"

# data2 for dgrp4
data2_4 <- subset(data2, dgrp==4)
length(unique(data2_4$cow)) #  84 cows
# cows with dgrp4 IDS  (dgrp=4)
cow_ID_4 <- subset(data2_4, select=c(cow))
cow_ID_4 <- unique(cow_ID_4)
summary(cow_ID_4)
rownames(cow_ID_4) <- NULL
length(unique(cow_ID_4$cow)) #  84cows
predictB<- predictB2 # Recalling predictB
sel <- predictB$cow %in% cow_ID_4$cow
predictB <- predictB[sel,]
predictB<- subset(predictB,select= c("cow","YDB_cow"))
predictB <- data.frame(predictB)
predictB <- predictB[,c("cow","YDB_cow")]

sel <- EBV_B_4_NA$cow %in% cow_ID_4$cow
EBV_B_4_NA <- EBV_B_4_NA[sel,]

EBV_YDB <- merge(predictB,EBV_B_4_NA, by="cow", all.x=TRUE)

# Cross Validation for accuracy for model fitB for cows in dgrp4
round(cor(EBV_YDB$YDB_cow,EBV_YDB$EBV_B),2) # -0.04

# The Average accuracy for model fitB
Acc_B <- round((0.05+0.14 + 0.14-0.04)/4,2)
Acc_B # 0.07





























#----Old (EAAP) cross_validation------------------------------------------------
# INLA Software
fitBase_NA1 <- inla(formula = modelfitBase, data = data1_1_NA,
                    control.compute = list(dic = TRUE,config=TRUE))

# save fitBase_NA1 as R object
save(fitBase_NA1,file = "data/cleaned_data/fitBase_pred/fitBase_NA1.RData")
#load(file = "data/cleaned_data/FitBase_pred/FitBase_NA1.RData") # to load the R object
pheno_pred1 <- fitBase_NA1$summary.linear.predictor
colnames(pheno_pred1)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred1$mean
data1$milkZ_pred_sd <- pheno_pred1$sd
pheno1 <-   subset(data1, dgrp==1)
pheno1<- subset(pheno1,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno1$cowI))# 731 Cows in breed composition class 1

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotyped
accuracy_fitBase_1 <- round(cor(pheno1$milkZ,pheno1$milkZ_pred),2)
Coef1<- lm (pheno1$milkZ~pheno1$milkZ_pred)
accuracy_fitBase_1 #  0.24
R2_fitBase_1<- summary(Coef1)
R2_fitBase_1 #  0.05932
#CRPS
obs <- pheno1$milkZ
pred<- subset(pheno1,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBase_1 <- crps(obs,pred)
round((crps_fitBase_1$CRPS),2) # 0.6

# Masking phenotypes of cows in dgrp 2 (Breed proportion class 2)
#dgrp2 masked

data1_2_NA<- data1 %>% mutate(milkZ = ifelse(dgrp == "2", NA, milkZ))
sum(is.na(data1_2_NA$milkZ)) # 8149
length(unique(data1_2_NA$cowI)) # 1894

fitBase_NA2 <- inla(formula = modelfitBase, data = data1_2_NA,
                    control.compute = list(dic = TRUE,config=TRUE))
save(fitBase_NA2,file = "data/cleaned_data/fitBase_pred/fitBase_NA2.RData")
#load(file = "data/cleaned_data/fitBase_pred/fitBase_NA2.RData") # to load the R object

pheno_pred2 <- fitBase_NA2$summary.linear.predictor
colnames(pheno_pred2)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred2$mean
data1$milkZ_pred_sd <- pheno_pred2$sd

pheno2 <-   subset(data1, dgrp==2)
pheno2<- subset(pheno2,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno2$cowI))#  770 Cows in breed composition class 2

# Correlating observed phenotype of cows in dgrp 2 with predicted phenotypes
accuracy_fitBase_2 <- round(cor(pheno2$milkZ,pheno2$milkZ_pred),2)
Coef2<- lm (pheno2$milkZ~pheno2$milkZ_pred)
accuracy_fitBase_2# 0.34
R2_fitBase_2<- summary(Coef2)
R2_fitBase_2 # 0.1123
#CRPS
obs <- pheno2$milkZ
pred<- subset(pheno2,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBase_2 <- crps(obs,pred)
round((crps_fitBase_2$CRPS),2) # 0.53

# Masking phenotypes of cows in dgrp 3 (Breed proportion class 3)
#dgrp3 masked
data1_3_NA<- data1 %>% mutate(milkZ = ifelse(dgrp == "3", NA, milkZ))
sum(is.na(data1_3_NA$milkZ)) # 3058 records
length(unique(data1_3_NA$cowI)) # 1894 cows in data1

fitBase_NA3 <- inla(formula = modelfitBase, data = data1_3_NA,
                    control.compute = list(dic = TRUE,config=TRUE))

save(fitBase_NA3,file = "data/cleaned_data/fitBase_pred/fitBase_NA3.RData")
#load(file = "data/cleaned_data/fitBase_pred/fitBase_NA3.RData") # to load the R object

pheno_pred3 <- fitBase_NA3$summary.linear.predictor
sum(is.na(pheno_pred3$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred3$mean
data1$milkZ_pred_sd <- pheno_pred3$sd

pheno3 <-   subset(data1, dgrp==3)
pheno3<- subset(pheno3,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno3$cowI))#  309 Cows in breed composition class 3

# Correlating observed phenotype of cows in dgrp 3 with predicted phenotypes
accuracy_fitBase_3 <- round(cor(pheno3$milkZ,pheno3$milkZ_pred),2)
Coef3<- lm (pheno3$milkZ~pheno3$milkZ_pred)
accuracy_fitBase_3# 0.28
R2_fitBase_3<- summary(Coef3)
R2_fitBase_3 #0.08083
#CRPS
obs <- pheno3$milkZ
pred<- subset(pheno3,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBase_3 <- crps(obs,pred)
round((crps_fitBase_3$CRPS),2) # 0.52


# Masking phenotypes of cows in dgrp 4 (Breed proportion class 4)
#dgrp4 masked
data1_4_NA<- data1 %>% mutate(milkZ = ifelse(dgrp == "4", NA, milkZ))
sum(is.na(data1_4_NA$milkZ)) # 702 records
length(unique(data1_4_NA$cowI)) # 1894

fitBase_NA4 <- inla(formula = modelfitBase, data = data1_4_NA,
                    control.compute = list(dic = TRUE,config=TRUE))
save(fitBase_NA4,file = "data/cleaned_data/fitBase_pred/fitBase_NA4.RData")
#load(file = "data/cleaned_data/fitBase_pred/fitBase_NA4.RData") # to load the R object

pheno_pred4 <- fitBase_NA4$summary.linear.predictor
sum(is.na(pheno_pred4$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred4$mean
data1$milkZ_pred_sd <- pheno_pred4$sd

pheno4 <-   subset(data1, dgrp==4)
pheno4<- subset(pheno4,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno4$cowI))#  84 Cows in breed composition class 4

# Correlating observed phenotype of cows in dgrp 4 with predicted phenotypes
accuracy_fitBase_4 <- round(cor(pheno4$milkZ,pheno4$milkZ_pred),2)
Coef4<- lm (pheno4$milkZ~pheno4$milkZ_pred)
accuracy_fitBase_4 # -0.05
R2_fitBase_4<- summary(Coef4)
R2_fitBase_4 # 0.0008917
#CRPS
obs <- pheno4$milkZ
pred<- subset(pheno4,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBase_4 <- crps(obs,pred)
round((crps_fitBase_4$CRPS),2) # 0.58


# Accuracy of prediction and degree under/overprediction of model FitBase

accuracy_fitBase = (accuracy_fitBase_1 + accuracy_fitBase_2 + accuracy_fitBase_3 + accuracy_fitBase_4)/4
round((accuracy_fitBase),2) # 0.2

R2_fitBase = (0.05932+0.1123+0.08083+0.0008917)/4
round((R2_fitBase),2) # 0.06

# crps_fitBase= (crps_fitBase_1+crps_fitBase_2+crps_fitBase_3+crps_fitBase_4)/4
crps_fitBase = (0.6+ 0.53 + 0.52 + 0.58)/4
crps_fitBase=round((crps_fitBase),2)
crps_fitBase #0.56

# Another method to calculate Accuracy, CRPS and R2 (degree of over/underprediction)
# Let's combine the predicted and observed phenotypes
pheno_fitBase <- rbind(pheno1,pheno2,pheno3,pheno4)
accuracy_fitBase<- round(cor(pheno_fitBase$milkZ,pheno_fitBase$milkZ_pred),2)
Coef_fitBase <- lm (pheno_fitBase$milkZ~pheno_fitBase$milkZ_pred)
accuracy_fitBase # 0.24
R2_fitBase<- summary(Coef_fitBase)
R2_fitBase# 0.05613
#CRPS_fitBase
obs <- pheno_fitBase$milkZ
pred<- subset(pheno_fitBase,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBase <- crps(obs,pred)
round((crps_fitBase$CRPS),2) # 0.55



































#Saving forward validation prediction fitBase External Drive
save(fitB_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitBase_NA.RData")
pheno_pred <- fitB_NA[,21:22]
summary(pheno_pred)

colnames(pheno_pred)
sum(is.na(pheno_pred$mean)) # 0 expected
data2_New$milkZ_pred <- pheno_pred$mean
data2_New$milkZ_pred_sd <- pheno_pred$sd
pheno <-   subset(data2_New, birthyear=="2016" | birthyear=="2017")
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno$cow))# 146 Cows expected to be born between 2016 (141 cows) and 2017 (5cows)
summary(pheno)
summary(data1)
# Correlating observed phenotype of young cows  with predicted phenotype
accuracy_fitB <- round(cor(pheno$milkZ,pheno$milkZ_pred),2)
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitB # 0.43
R2_fitBase<- summary(Coef)
R2_fitBase #   0.05162
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitB <- crps(obs,pred)
round((crps_fitBase$CRPS),2) # 0.56

#-----------------------Foward validation FitBH-------------------------------
modelBH_INLA <-   "milkZ ~ 1 + cyrsnI + tyrmnI + dgrpI + (ageZ|lacgr) +
  (leg1|lacgr) + (leg2|lacgr) +
  f(herd, model = 'iid') +
  f(cowPe, model = 'iid') +
  f(ward_code, model = 'iid') +
  f(cowI, model = 'generic0', Cmatrix = GRMInv)"

fitBH_NA <- inla(formula = modelBH_INLA, data = data2_New_NA,
                  control.compute = list(dic = TRUE,config=TRUE))
#Saving forward validation prediction fitWCF External Drive
save(fitBH_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitWCF_NA.RData")
pheno_pred <- fitBH_NA$summary.linear.predictor
colnames(pheno_pred)
sum(is.na(pheno_pred$mean)) # 0 expected
data2_New$milkZ_pred <- pheno_pred$mean
data2_New$milkZ_pred_sd <- pheno_pred$sd
pheno <-   subset(data2_New, birthyear=="2016" | birthyear=="2017")
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno$cowI))# 146 Cows expected to be born between 2016 (141 cows) and 2017 (5cows)

# Correlating observed phenotype of young cows  with predicted phenotype
accuracy_fitBH <- round(cor(pheno$milkZ,pheno$milkZ_pred),2)
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitBH #  0.54
R2_fitBH<- summary(Coef)
R2_fitBH #  0.2878
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBH <- crps(obs,pred)
round((crps_fitWCF$CRPS),2) #

#-----------------------Foward validation fitBS-------------------------------
data2_New$cyrsnI <- data1$cyrsnI
data2_New_NA$cyrsnI <- data1$cyrsnI
data2_New$dgrpI <- data1$dgrpI
data2_New_NA$dgrpI <- data1$dgrpI
data2_New$tyrmnI <- data1$tyrmnI
data2_New_NA$tyrmnI <- data1$tyrmnI

#Priors
hyperRange <- c(50, 0.8)
hyperVarSpdeS <- c(sqrt(0.25), 0.5)
hyperVarSpdeWS <- c(sqrt(0.10), 0.5)
hyperResVarGWS <- list(theta = list(prior = "pc.prec", param = c(sqrt(0.15),0.5)))

spdeStatS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeS)
spdeStatWS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeWS)
meshIndexS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatS$n.spde)
meshIndexWS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatWS$n.spde)
A = inla.spde.make.A(mesh = mesh, loc = cbind(data2_New$long, data2_New$lat) )
A.pred = inla.spde.make.A(mesh = mesh, loc = cbind(data2_New$long, data2_New$lat))

# Make stack
# StackS
stackS = inla.stack(data = list(milkZ = data2_New$milkZ),
                    A = list(A,1),
                    effects = list(c(meshIndexS, list(intercept = 1)),
                                   list(cowI = data2_New$cowI,cowPe = data2_New$cowPe,
                                        cyrsnI=data2_New$cyrsnI, tyrmnI=data2_New$tyrmnI,dgrpI= data2_New$dgrpI, ward_code=data2_New$ward_code, ageZ=data2_New$ageZ, lacgr=data2_New$lacgr,leg1=data2_New$leg1,leg2=data2_New$leg2)), tag = "data2_New.data")


# StackWS (herd as random + Spatial effect)
stackWS = inla.stack(data = list(milkZ = data2_New$milkZ),
                     A = list(A,1),
                     effects = list(c(meshIndexWS, list(intercept = 1)),
                                    list(cowI = data2_New$cowI, cowPe = data2_New$cowPe,
                                         cyrsnI=data2_New$cyrsnI,herd=data2_New$herd, tyrmnI=data2_New$tyrmnI,dgrpI= data2_New$dgrpI,ward_code=data2_New$ward_code, ageZ=data2_New$ageZ, lacgr=data2_New$lacgr,leg1=data2_New$leg1,leg2=data2_New$leg2)), tag = "data2_New.data")

# Make stack_pred
# StackS_NA
stackS_pred_NA = inla.stack(data = list(milkZ = NA),
                             A = list(A.pred,1),
                             effects = list(c(meshIndexS, list(intercept = 1)),
                                            list(cowI = data2_New_NA$cowI, cowPe = data2_New_NA$cowPe,
                                                 cyrsnI=data2_New_NA$cyrsnI,tyrmnI=data2_New_NA$tyrmnI,dgrpI= data2_New_NA$dgrpI, ward_code=data2_New_NA$ward_code,ageZ=data2_New_NA$ageZ, lacgr=data2_New_NA$lacgr,leg1=data2_New_NA$leg1,leg2=data2_New_NA$leg2)), tag = "data2_New_NA.data")


# Create joint stack
join.stack_NA <- inla.stack(stackS, stackS_pred_NA)

# ModelBS
formulaBS <- as.formula(paste0(modelB_INLA, " + f(fieldID, model = spdeStatS)-1"))


fitBS_NA= inla(formula = formulaBS, data = inla.stack.data(join.stack_NA),
               family = "normal", control.predictor =list(A=inla.stack.A(join.stack_NA),compute = T),
               control.family=list(list(hyper=hyperResVarGWS)),
               control.compute = list(dic=T,cpo=F, config=T),
               control.fixed = list(expand.factor.strategy="inla"), verbose=T)

# save fitS_NA as R object
#save(fitBS_NA,file ="data/cleaned_data/foward_valid/fitBS_NA.RData")

#save(fitBS_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitS_NA.RData")
#Get the prediction index

pred.ind_NA <- inla.stack.index(join.stack_NA, tag='data2_New_NA.data')$data
pheno_pred <- fitBS_NA$summary.linear.predictor [pred.ind_NA,]

sum(is.na(pheno_pred$mean)) # 0 expected
data2_New$milkZ_pred <- pheno_pred$mean
data2_New$milkZ_pred_sd <- pheno_pred$sd
pheno <-   subset(data2_New, birthyear=="2016" | birthyear=="2017")
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno$cowI))# 146 Cows expected

# Correlating observed phenotype of young cows  with predicted phenotyped
accuracy_fitBS <- round(cor(pheno$milkZ,pheno$milkZ_pred),2)
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitBS # 0.81
R2_fitS<- summary(Coef)
R2_fitS #
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitS <- crps(obs,pred)
round((crps_fitS$CRPS),2) # 0.31

#-----------------------Foward validation model fitBHS-------------------------------

# Make stack_pred
# StackWS_NA
stackWS_pred_NA = inla.stack(data = list(milkZ = NA),
                             A = list(A.pred,1),
                             effects = list(c(meshIndexWS, list(intercept = 1)),
                                            list(cowI = data2_New_NA$cowI, cowPe = data2_New_NA$cowPe,
                                                 cyrsnI=data2_New_NA$cyrsnI,tyrmnI=data2_New_NA$tyrmnI,dgrpI= data2_New_NA$dgrpI, ward_code=data2_New_NA$ward_code,ageZ=data2_New_NA$ageZ, lacgr=data2_New_NA$lacgr,leg1=data2_New_NA$leg1,leg2=data2_New_NA$leg2)), tag = "data2_New_NA.data")


# Create joint stack
join.stack_NA <- inla.stack(stackWS, stackWS_pred_NA)

# ModelS
formulaWS <- as.formula(paste0(modelBH_INLA, " + f(fieldID, model = spdeStatWS)-1"))


fitBHS_NA= inla(formula = formulaWS, data = inla.stack.data(join.stack_NA),
              family = "normal", control.predictor =list(A=inla.stack.A(join.stack_NA),compute = T),
              control.family=list(list(hyper=hyperResVarGWS)),
              control.compute = list(dic=T,cpo=F, config=T),
              control.fixed = list(expand.factor.strategy="inla"), verbose=T)

# save fitBHS_NA as R object
save(fitBHS_NA,file ="data/cleaned_data/foward_valid/fitBHS_NA.RData")
#save(fitWS_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitWS_NA.RData")

#Get the prediction index

pred.ind_NA <- inla.stack.index(join.stack_NA, tag='data2_New_NA.data')$data
pheno_pred <- fitBHS_NA$summary.linear.predictor [pred.ind_NA,]

sum(is.na(pheno_pred$mean)) # 0 expected
data2_New$milkZ_pred <- pheno_pred$mean
data2_New$milkZ_pred_sd <- pheno_pred$sd
pheno <-   subset(data2_New, birthyear=="2016" | birthyear=="2017")
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno$cowI))# 146 Cows expected

# Correlating observed phenotype of young cows  with predicted phenotyped
accuracy_fitBHS <- round(cor(pheno$milkZ,pheno$milkZ_pred),2)
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitBHS # 0.74
R2_fitBHS<- summary(Coef)
R2_fitBHS #
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWS <- crps(obs,pred)
round((crps_fitWS$CRPS),2) # 0.31
#-------------------------------------------------------------------------------




#------Accuracy of prediction: No_herd Foward validation--------------------------------

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
length(unique(cowbdate$cow))

# Create birthyear column in data1 by merging with cowbdate

data1 <-  merge(data1, cowbdate, by= "cow", all.x = TRUE)
data1$birthyear <- as.factor(data1$birthyear)
summary(data1$birthyear)
#2004  2005  2006  2007  2008  2009  2010  2011  2012  2013  2014  2015  2016  2017
#17    12     7    38    41     7   122   232  3044  1066 12017  1298  1444    30

#data1_2016_2017 <- subset(data1, birthyear==2016 | birthyear==2017 )

#-----------------------Foward validation FitbaseNoherd-------------------------------
#Let's make NA milkz of cows born in 2016 and 2017
data1_NA<- data1 %>% mutate(milkZ = ifelse(birthyear=="2016" | birthyear=="2017", NA, milkZ))
sum(is.na(data1_NA$milkZ)) # 1474 Expected
length(unique(data1_NA$cowI)) # 1894 expected
fitBaseNoherd_NA <- inla(formula = modelBaseNoherd, data = data1_NA,
                   control.compute = list(dic = TRUE,config=TRUE))

#Saving forward validation prediction fitBase External Drive
save(fitBaseNoherd_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitBaseNoherd_NA.RData")
pheno_pred <- fitBaseNoherd_NA$summary.linear.predictor
colnames(pheno_pred)
sum(is.na(pheno_pred$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred$mean
data1$milkZ_pred_sd <- pheno_pred$sd
pheno <-   subset(data1, birthyear=="2016" | birthyear=="2017")
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno$cowI))# 146 Cows expected to be born between 2016 (141 cows) and 2017 (5cows)

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotype
accuracy_fitBaseNoherd <- round(cor(pheno$milkZ,pheno$milkZ_pred),2)
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitBaseNoherd # 0.19
R2_fitBaseNoherd<- summary(Coef)
R2_fitBaseNoherd # 0.03427
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBaseNoherd <- crps(obs,pred)
round((crps_fitBaseNoherd$CRPS),2) #  0.6

#-----------------------Foward validation FitWCF-------------------------------
fitWCFNoherd_NA <- inla(formula = modelWCFNoherd, data = data1_NA,
                  control.compute = list(dic = TRUE,config=TRUE))


#Saving forward validation prediction fitWCF External Drive
save(fitWCFNoherd_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitWCFNoherd_NA.RData")
pheno_pred <- fitWCFNoherd_NA$summary.linear.predictor
colnames(pheno_pred)
sum(is.na(pheno_pred$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred$mean
data1$milkZ_pred_sd <- pheno_pred$sd
pheno <-   subset(data1, birthyear=="2016" | birthyear=="2017")
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno$cowI))# 146 Cows expected to be born between 2016 (141 cows) and 2017 (5cows)

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotype
accuracy_fitWCFNoherd <- round(cor(pheno$milkZ,pheno$milkZ_pred),2)
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitWCFNoherd #  0.19
R2_fitWCFNoherd<- summary(Coef)
R2_fitWCFNoherd # 0.0342
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCFNoherd <- crps(obs,pred)
round((crps_fitWCFNoherd$CRPS),2) #0.6

#-----------------------Foward validation fitWCRI-------------------------------
fitWCRINoherd_NA <- inla(formula = modelWCRINoherd, data = data1_NA,
                   control.compute = list(dic = TRUE,config=TRUE))


#Saving forward validation prediction fitWCRI External Drive
#save(fitWCRI_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitWCRI_NA.RData")
pheno_pred <- fitWCRINoherd_NA$summary.linear.predictor
colnames(pheno_pred)
sum(is.na(pheno_pred$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred$mean
data1$milkZ_pred_sd <- pheno_pred$sd
pheno <-   subset(data1, birthyear=="2016" | birthyear=="2017")
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno$cowI))# 146 Cows expected to be born between 2016 (141 cows) and 2017 (5cows)

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotype
accuracy_fitWCRINoherd <- round(cor(pheno$milkZ,pheno$milkZ_pred),2)
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitWCRINoherd #  0.52
R2_fitWCRINoherd<- summary(Coef)
R2_fitWCRINoherd # 0.2688
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCRINoherd <- crps(obs,pred)
round((crps_fitWCRINoherd$CRPS),2) # 0.5

#-----------------------Foward validation fitWCRB-------------------------------
fitWCRBNoherd_NA <- inla(formula = modelWCRBNoherd, data = data1_NA,
                   control.compute = list(dic = TRUE,config=TRUE))


#Saving forward validation prediction fitWCRB External Drive
save(fitWCRB_NANoherd,file = "D:/Results_ADGG_Spatial/forward_valid/fitWCRB_NA.RData")
pheno_pred <- fitWCRBNoherd_NA$summary.linear.predictor
colnames(pheno_pred)
sum(is.na(pheno_pred$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred$mean
data1$milkZ_pred_sd <- pheno_pred$sd
pheno <-   subset(data1, birthyear=="2016" | birthyear=="2017")
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno$cowI))# 146 Cows expected to be born between 2016 (141 cows) and 2017 (5cows)

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotype
accuracy_fitWCRBNoherd <- round(cor(pheno$milkZ,pheno$milkZ_pred),2)
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitWCRBNoherd # 0.52
R2_fitWCRBNoherd<- summary(Coef)
R2_fitWCRBNoherd #0.2677
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCRBNoherd <- crps(obs,pred)
round((crps_fitWCRBNoherd$CRPS),2) # 0.5

#-----------------------Foward validation model fitSNoherd-------------------------------

A.pred = inla.spde.make.A(mesh = mesh, loc = cbind(data1$long, data1$lat) )

# Make stack_pred
# StackSNoherd_NA
stackSNoherd_pred_NA = inla.stack(data = list(milkZ = NA),
                            A = list(A.pred,1),
                            effects = list(c(meshIndexS, list(intercept = 1)),
                                           list(cowI = data1_NA$cowI, cowPeI = data1_NA$cowPeI,
                                                cyrsnI=data1_NA$cyrsnI, tyrmnI=data1_NA$tyrmnI,dgrpI= data1_NA$dgrpI, ageZ=data1_NA$ageZ, lacgr=data1_NA$lacgr, leg0=data1_NA$leg0, leg1=data1_NA$leg1,leg2=data1_NA$leg2)), tag = "data1Noherd_NA.data")

# Create joint stack
join.stackSNoherd_NA <- inla.stack(stackSNoherd, stackSNoherd_pred_NA)

# ModelS


fitSNoherd_NA= inla(formula = formulaSNoherd, data = inla.stack.data(join.stackNoherd_NA),
              family = "normal", control.predictor =list(A=inla.stack.A(join.stackNoherd_NA),compute = T),
              control.family=list(list(hyper=hyperResVarGWS)),
              control.compute = list(dic=T,cpo=F, config=T),
              control.fixed = list(expand.factor.strategy="inla"), verbose=T)

# save fitS_NA as R object
save(fitSNoherd_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitSNoherd_NA.RData")

#Get the prediction index

pred.ind_NA <- inla.stack.index(join.stackNoherd_NA, tag='data1Noherd_NA.data')$data
pheno_pred <- fitSNoherd_NA$summary.linear.predictor [pred.ind_NA,]

sum(is.na(pheno_pred$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred$mean
data1$milkZ_pred_sd <- pheno_pred$sd
pheno <-   subset(data1, birthyear=="2016" | birthyear=="2017")
pheno<- subset(pheno,ect= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno$cowI))# 146 Cows expected

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotyped
accuracy_fitSNoherd <- round(cor(pheno$milkZ,pheno$milkZ_pred),2)
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitSNoherd # 0.86
R2_fitSNoherd<- summary(Coef)
R2_fitSNoherd # 0.7445
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitSNoherd <- crps(obs,pred)
round((crps_fitSNoherd$CRPS),2) # 0.31

#-----------------------Foward validation model fitWS-------------------------------

# Make stack_pred
# StackWS_NA
stackWSNoherd_pred_NA = inla.stack(data = list(milkZ = NA),
                             A = list(A.pred,1),
                             effects = list(c(meshIndexWS, list(intercept = 1)),
                                            list(cowI = data1_NA$cowI, cowPeI = data1_NA$cowPeI,
                                                 cyrsnI=data1_NA$cyrsnI, tyrmnI=data1_NA$tyrmnI,dgrpI= data1_NA$dgrpI, ageZ=data1_NA$ageZ, lacgr=data1_NA$lacgr, leg0=data1_NA$leg0, leg1=data1_NA$leg1,leg2=data1_NA$leg2)), tag = "data1Noherd_NA.data")

# Create joint stack
join.stackWSNoherd_NA <- inla.stack(stackWSNoherd, stackWSNoherd_pred_NA)

# ModelS

fitWSNoherd_NA= inla(formula = formulaWSNoherd, data = inla.stack.data(join.stackWSNoherd_NA),
               family = "normal", control.predictor =list(A=inla.stack.A(join.stackWSNoherd_NA),compute = T),
               control.family=list(list(hyper=hyperResVarGWS)),
               control.compute = list(dic=T,cpo=F, config=T),
               control.fixed = list(expand.factor.strategy="inla"), verbose=T)

# save fitWS_NA as R object
save(fitWS_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitWS_NA.RData")

#Get the prediction index

pred.ind_NA <- inla.stack.index(join.stackNoherd_NA, tag='data1Noherd_NA.data')$data
pheno_pred <- fitWSNoherd_NA$summary.linear.predictor [pred.ind_NA,]

sum(is.na(pheno_pred$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred$mean
data1$milkZ_pred_sd <- pheno_pred$sd
pheno <-   subset(data1, birthyear=="2016" | birthyear=="2017")
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))
length(unique(pheno$cowI))# 146 Cows expected

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotyped
accuracy_fitWSNoherd <- round(cor(pheno$milkZ,pheno$milkZ_pred),2)
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitWSNoherd # 0.84
R2_fitWSNoherd<- summary(Coef)
R2_fitWSNoherd #  0.7112
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWSNoherd <- crps(obs,pred)
round((crps_fitWSNoherd$CRPS),2) # 0.36



hist(data1$milk)
summary(data1$milk)
mean(data1$milk)
sd(data1$milk)
8.3-c(1,2,3)*4.3
8.3+c(1,2,3)*4.3

sel <- data1$milk < 20
hist(data1$milk[sel])
summary(data1$milk[sel])
mean(data1$milk[sel])
sd(data1$milk[sel])

#-----------------------Principal Components analysis---------------------------
# Import SNP data (already QCed by Raphael so we will just take it as it is!)
geno <- read.table(file = "data/original_data/snpref-isi.txt", header = FALSE)

geno[1:10, 1:10]
dim(geno) # 1911 664823
colnames(geno)
summary(geno$V1) # IDs from 1 to 1916
colnames(geno)[1] <- "cow"
#Save geno as an R object
#save(geno,file = "data/cleaned_data/geno.RData")
load(file = "data/cleaned_data/geno.RData")
#install.packages(pkg = "irlba")
library(package = "irlba")

# stats::prcomp(x = geno[1:10, c(2:20)])
pcasAll <- prcomp_irlba(x = geno[, -1], n = 3)
summary(pcasAll)
summary(pcas)
par(mfrow = c(2, 2))
plot(pcasAll$x[, 1], pcasAll$x[, 2], pch = 19, cex = 0.1)
plot(pcasAll$x[, 1], pcasAll$x[, 3], pch = 19, cex = 0.1)
plot(pcasAll$x[, 2], pcasAll$x[, 3], pch = 19, cex = 0.1)
par(mfrow = c(1, 1))
plot(pcasAll$x[, 1], pcasAll$x[, 2], pch = 19, cex = 0.5)

pcas <- as.data.frame(cbind(geno[, 1], pcasAll$x))
colnames(pcas)[1] <- "cow"
head(pcas)
sel <- pcas$cow %in% data1$cow
pcas <- pcas[sel, ]
dim(pcas)
tmp <- data1[, c("cow", "dgrp", "region")]
tmp <- tmp[!duplicated(tmp), ]
dim(tmp)
head(tmp)
head(pcas)
pcas <- merge(x = pcas, y = tmp, by = "cow", all.x = TRUE)
dim(pcas)
str(pcas)
plot()

# Visualize PCA by breed Proportion dgrp
# PC1 vs PC2

dgrp_pca1 <- ggplot(data = pcas) +
  geom_point(mapping = aes(x = `PC1`, y = `PC2`,color = dgrp), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "PCA of SNP genotypes by breed proportion",
       x = paste0("PC1(3.1%)"),
       y = paste0("PC2(0.8%)")) +
  theme_minimal()


# PC1 vs PC3
dgrp_pca2<- ggplot(data = pcas) +
  geom_point(mapping = aes(x = `PC1`, y = `PC3`,color = dgrp), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x = paste0("PC1(3.1%)"),
       y = paste0("PC3(0.5%)")) +
  theme_minimal()

# PC2 vs PC3
dgrp_pca3<- ggplot(data = pcas) +
  geom_point(mapping = aes(x = `PC2`, y = `PC3`,color = dgrp), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x = paste0("PC2(0.8%)"),
       y = paste0("PC3(0.5%)")) +
  theme_minimal()

dgrp_pca <- ggarrange(dgrp_pca1,dgrp_pca2,dgrp_pca3, ncol = 3, nrow = 1, common.legend = T, legend = "right", align = "h", widths = c(1,1,1))
dgrp_pca

ggsave(plot = dgrp_pca + PreseTheme, filename = "PCA_breed_proportion_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = dgrp_pca + PaperTheme, filename = "PCA_breed_proportion_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper


# Visualize PCA by region
# PC1 vs PC2
region_pca1<- ggplot(data = pcas) +
  geom_point(mapping = aes(x = `PC1`, y = `PC2`,color = region), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "PCA of SNP genotypes across regions",
       x = paste0("PC1(3.1%)"),
       y = paste0("PC2(0.8%)")) +
  theme_minimal()


# PC1 vs PC3
region_pca2<- ggplot(data = pcas) +
  geom_point(mapping = aes(x = `PC1`, y = `PC3`,color = region), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x = paste0("PC1(3.1%)"),
       y = paste0("PC3(0.5%)")) +
  theme_minimal()

# PC2 vs PC3
region_pca3<- ggplot(data = pcas) +
  geom_point(mapping = aes(x = `PC2`, y = `PC3`,color = dgrp), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x = paste0("PC2(0.8%)"),
       y = paste0("PC3(0.5%)")) +
  theme_minimal()


region_pca <- ggarrange(region_pca1,region_pca2,region_pca3, ncol = 3, nrow = 1, common.legend = T, legend = "right", align = "h", widths = c(1,1,1))
region_pca

ggsave(plot = region_pca + PreseTheme, filename = "PCA_region_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = region_pca + PaperTheme, filename = "PCA_region_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper



table(data1$dgrp,data1$region)
#Breed proportion in First column
#    1    2    3    4
#1  497 4663 1076 1230
#2 2248 3012 1704 1185
#3 1495  704  622  237
#4  307  103  272   20


#-----------------------GRM vs Spatial Relationship Matrix---------------------------
# LOAD GRM
load(file = "data/cleaned_data/GRM.RData")
dim(GRM)
colnames(GRM) <- 1:ncol(GRM)
rownames(GRM) <- 1:nrow(GRM)
class(GRM) # "matrix" "array"
GRM_df <- data.frame(cow_i=rep(row.names(GRM),ncol(GRM)),
                cow_j=rep(colnames(GRM),each=nrow(GRM)),
                animal_cov=as.vector(GRM))

cow_herd <- data1[,c(1,4)]
length(unique(cow_herd$cow)) #1894
length(unique(cow_herd$herd)) #1386


# Add herd id for cow i
cow_herd_i <- cow_herd
names(cow_herd_i)[1] <- "cow_i"
names(cow_herd_i)[2] <- "herd_i"

cow_herd_j <- cow_herd
names(cow_herd_j)[1] <- "cow_j"
names(cow_herd_j)[2] <- "herd_j"

dim(GRM_df)
#Adding gps coordinates to herds
herd_GPS <- data1[,c(4,11:12)]
herd_GPS <- distinct(herd_GPS)
dim(herd_GPS)


# Add GPS for herd i
herd_GPS_i <- herd_GPS
names(herd_GPS_i)[1] <- "herd_i"
names(herd_GPS_i)[2] <- "long_herd_i"
names(herd_GPS_i)[3] <- "lat_herd_i"
# Add GPS for herd j
herd_GPS_j <- herd_GPS
names(herd_GPS_j)[1] <- "herd_j"
names(herd_GPS_j)[2] <- "long_herd_j"
names(herd_GPS_j)[3] <- "lat_herd_j"

#Merge Cow-herdi with herd_i GPS
head(cow_herd_i)
head(herd_GPS_i)
Merge1 <- merge(herd_GPS_i, cow_herd_i, by="herd_i", all.x = TRUE)
Merge1<- distinct(Merge1)

#Merge Cow-herdj with herd_j GPS
head(cow_herd_j)
head(herd_GPS_j)
Merge2 <- merge(herd_GPS_j, cow_herd_j, by="herd_j", all.x = TRUE)
Merge2<- distinct(Merge2)
dim(distinct(GRM_df))

Merge3<- merge(GRM_df, Merge1, by="cow_i", all.x = TRUE)
dim(distinct(Merge3))
head(Merge3)
head(Merge2)

Merge4 <- merge(Merge2, Merge3, by="cow_j", all.x = TRUE)

GRM_distance_herd<- Merge4 %>% mutate(distance_herd_i_herd_j=sqrt((long_herd_i-long_herd_j)^2 + (lat_herd_i-lat_herd_j)^2))

summary(GRM_df_final$dij)
# Let's group relationship (score) by herd_i and herd_j
GRM_df_final$herd_ij <- with(GRM_df_final, paste0(herd_i, "-", herd_j))

GRM_SRM <- GRM_df_final %>%  group_by(herd_ij) %>%
  e(mean_score = mean(score), mean_dij =mean(dij), n = n())

length(unique(GRM_SRM$herd_ij))
save(GRM_SRM,file = "data/cleaned_data/GRM_SRM.RData")

summary(GRM_SRM)
head(GRM_SRM)
#plot(GRM_df_final$score~ GRM_df_final$dij)

GRM_Plot<- ggplot(data = GRM_SRM) +
  geom_point(mapping = aes(x = `mean_dij`, y = `mean_score`)) +
  geom_hline(yintercept = 0, linetype="dotted", linewidth=0.001) +
  geom_vline(xintercept = 0, linetype="dotted", linewidth=0.001) +
  labs(title = "Genomic kinship versus distance between herds",
       x = paste0("Euclidian distance between herds"),
       y = paste0("Genomic relationship coefficient")) +
  theme_minimal()

GRM_Plot

ggsave(plot = GRM_Plot + PreseTheme, filename = "GRM_SRM_Presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = GRM_Plot + PaperTheme, filename = "GRM_SRM__paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper


Same_herd_to_herd_distance <- subset(GRM_SRM, mean_dij=="0")
head(GRM_SRM)




GRM_Plot2 <- smoothScatter(GRM_SRM$mean_dij,GRM_SRM$mean_dij,
                     ## pch=NA: do not draw them
                     nrpoints = 2000, ret.selection=TRUE)


d=densCols(GRM_SRM$mean_dij, GRM_SRM$mean_score)
GRM_Plot2  <- ggplot(GRM_SRM) +
  geom_hexbin(aes(mean_dij, mean_score, col = d), size = 0.1) +
  scale_color_identity() +
  theme_bw()
print(GRM_Plot2)

x= GRM_SRM$mean_dij
y=GRM_SRM$mean_score
df <- data.frame(x = x, y = y,
                 d = densCols(x, y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
p <- ggplot(df) +
  geom_point(aes(x, y, col = d), size = 0.1) +
  scale_color_identity() +
  labs(title = "Genomic kinship versus distance between herds",
       x = paste0("Euclidian distance between herds"),
       y = paste0("Genomic relationship coefficient")) +
  theme_bw()
print(p)

cor.test(GRM_SRM$mean_dij,GRM_SRM$mean_score)



GRM_df_final_Reduced <- GRM_df_final[,c(5,6,10)]

dim(GRM_df_final_Reduced)

x= GRM_df_final$dij
y=GRM_df_final$score

df <- data.frame(x = x, y = y,
                 d = densCols(x, y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
p <- ggplot(df) +
  geom_point(aes(x, y, col = d), size = 0.1) +
  scale_color_identity() +
  labs(title = "Individual Genomic kinship versus distance between herds",
       x = paste0("Euclidian distance between herds"),
       y = paste0("Genomic relationship coefficient")) +
  theme_bw()
print(p)




x= GRM_df_final$dij
y=GRM_df_final$score

df <- data.frame(x = x, y = y,
                 d = densCols(x, y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
p <- ggplot(df) +
  geom_point(aes(x, y, col = d), size = 0.1) +
  scale_color_identity() +
  labs(title = "Individual Genomic kinship versus distance between herds",
       x = paste0("Euclidian distance between herds"),
       y = paste0("Genomic relationship coefficient")) +
  theme_bw()
print(p)



GRM_df_final_i <- GRM_df_final %>% group_by(i) %>% summarise(score_i=mean(score), mean_dij=mean(dij))


x= GRM_df_final_i$mean_dij
y=GRM_df_final_i$score_i


df <- data.frame(x = x, y = y,
                 d = densCols(x, y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
p <- ggplot(df) +
  geom_point(aes(x, y, col = d), size = 1) +
  scale_color_identity() +
  labs(title = "Individual Genomic kinship versus distance between herds",
       x = paste0("Euclidian distance between herds"),
       y = paste0("Genomic relationship coefficient")) +
  theme_bw()
print(p)

cor.test(GRM_df_final_i$mean_dij,GRM_df_final_i$score_i)

str(GRM_df_final_i)

summary(GRM_df_final_i$score_i)

#SmoothScatter

smooth_plot<- smoothScatter(GRM_distance_herd$animal_cov ~ GRM_distance_herd$distance_herd_i_herd_j,
              nrpoints = 1000, pch = 1,bandwidth = 0.05, col = "black")


smoothScatter(GRM_distance_herd$animal_cov ~ GRM_distance_herd$distance_herd_i_herd_j,
             pch = 1,bandwidth = 0.05, col = "black")




smooth_plot + labs(y="y", x="x")





























### Extract breed names
#fam <- data.frame(famids=read.table("dataForPCA.mdist.id")[,1])
### Extract individual names
#famInd <- data.frame(IID=read.table("dataForPCA.mdist.id")[,2])

## Perform PCA using the cmdscale function
# Time intensive step - takes a few minutes with the 4.5K animals
mds_populations <- cmdscale(dataForPCA,eig=T,5)

## Extract the eigen vectors
#eigenvec_populations <- cbind(fam,famInd,mds_populations$points)

eigenvec_populations <- data.frame(mds_populations$points)

## Proportion of variation captured by each eigen vector
eigen_percent <- round(((mds_populations$eig)/sum(mds_populations$eig))*100,2)

eigen_percent <- round(((mds_populations$eig)/sum(mds_populations$eig))*100,2)
eigen_percent <- data.frame(eigen_percent)
# Visualize PCA
ggplot(data = eigenvec_populations) +
  geom_point(mapping = aes(x = `X1`, y = `X2`), show.legend = FALSE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "PCA of wordwide goat populations",
       x = paste0("Principal component 1 (",eigen_percent[1,]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,]," %)")) +
  theme_minimal()



#-----------------------Manuscript- Correlations and Animal Ranking------------------------------------------
#fitB
EBV_B <- fitB$summary.random$animal[,1:3]
names(EBV_B)[2] <- "EBV_B"
names(EBV_B)[3] <- "sd_B"


# fitBH
EBV_BH <- fitBH$summary.random$animal[,1:3]
names(EBV_BH)[2] <- "EBV_BH"
names(EBV_BH)[3] <- "sd_BH"

# fitBS

EBV_BS <- fitBS$summary.random$animal[,1:3]
names(EBV_BS)[2] <- "EBV_BS"
names(EBV_BS)[3] <- "sd_BS"


# fitBHS

EBV_BHS <- fitBHS$summary.random$animal[,1:3]
names(EBV_BHS)[2] <- "EBV_BHS"
names(EBV_BHS)[3] <- "sd_BHS"

#----------Combined EBVs manuscripts models-------------------------------------
EBV_sd<- merge(EBV_B,EBV_BH, by="ID")
EBV_sd<- merge(EBV_sd, EBV_BS,by="ID")
EBV_sd<- merge(EBV_sd, EBV_BHS, by="ID")

EBV <- EBV_sd[,c("ID","EBV_B","EBV_BH","EBV_BS","EBV_BHS")]

#--------------------Predict spatial effect at my locations (Models BS and BHS)-
# Spatial effect in model BS
fieldBS <- predict(fitBS,spdf, ~field)

fieldBS_df <- data.frame(fieldBS[,c("cow" , "herd","mean", "sd")])
fieldBS_df <- distinct(fieldBS_df)
length(levels(fieldBS_df$herd)) # 1386 herds
# Spatial effect BS
Spatial_BS <- fieldBS_df[,c("cow","herd", "mean")]
names(Spatial_BS)[1]<- "ID"
names(Spatial_BS)[3]<- "Spatial_effect_BS"

# Spatial effect in model BHS
fieldBHS <- predict(fitBHS,spdf, ~field)

fieldBHS_df <- data.frame(fieldBHS[,c("cow" , "herd","mean", "sd")])
fieldBHS_df <- distinct(fieldBHS_df)

# Spatial effect BHS
Spatial_BHS <- fieldBHS_df[,c("cow","herd", "mean")]
names(Spatial_BHS)[1]<- "ID"
names(Spatial_BHS)[3]<- "Spatial_effect_BHS"

# ------------EBVs and Spatial effect-------------------------------------------

#EBV_Spatial<- data.frame(EBV_B$ID,EBV_B$EBV_B,EBV_WCF$EBV_WCF, EBV_WCRI$EBV_WCRI,EBV_WCRB$EBV_WCRB,EBV_S$EBV_S, EBV_WS$EBV_WS,SPDE_cow_S$spdeS_mean)
#names(EBV_Spatial)[1:8] <- c("ID","EBV_B","EBV_WCF","EBV_WCRI","EBV_WCRB","EBV_S","EBV_WS","Spatial_effect_S")
# Remove herd from Spatial_BS
Spatial_BS <- Spatial_BS[,c("ID", "Spatial_effect_BS")]
EBVs_Spatial_effect <- merge(EBV,Spatial_BS, by="ID", all.x=TRUE)
EBVs_Spatial_effect <- merge(EBVs_Spatial_effect, Spatial_BHS, by="ID", all.x=TRUE )
plot(EBVs_Spatial_effect$Spatial_effect_BS,EBVs_Spatial_effect$EBV_B)
EBVs_Spatial_Matrix <- EBVs_Spatial_effect[, -c(1)]

cor(EBVs_Spatial_effect$EBV_BH,EBVs_Spatial_effect$EBV_BHS)

#create matrix of correlation coefficients and p-values for EBVs between models
# How to Create a Correlation Matrix in R (4 Examples) - Statology
# https://www.statology.org/correlation-matrix-in-r/

  # Pearson Correlation
round(cor(EBVs_Spatial_Matrix, method = "pearson"),2)

# Spearman's Rank Correlations
round(cor(EBVs_Spatial_Matrix, method = "spearman"),2)



#Adding difference EBV_BH - EBV_BHS

EBVs_Spatial_effect <- EBVs_Spatial_effect %>% mutate(dEBV= EBV_BH - EBV_BHS)
EBVs_Spatial_effect
round(cor(EBVs_Spatial_effect$dEBV,EBVs_Spatial_effect$Spatial_effect_BHS),2)
#0.29

#--------Adding average milk yield and herd effect
#Milk yield
colnames(EBVs_Spatial_effect)
data2_MY <- data.frame(data2[, c("cow","milkZ")])
data2_MY <- data2_MY[, c("cow","milkZ")]
data2_MY <- data2_MY %>% group_by(cow)  %>%
  summarise(mean_milkZ= mean(milkZ))

EBVs_Spatial_effect <- merge(EBVs_Spatial_effect, data2_MY, by= "cow", all.x = TRUE)

#Adding herd effects from model BHS
herd_effect_BHS <- data.frame(fitBHS$summary.random$herd[,c("ID","mean")])
names(herd_effect_BHS)[1]<- "herd"
names(herd_effect_BHS)[2]<- "herd_effect"

#Merging herd_effect_BHS with EBVs_Spatial_effect
EBVs_Spatial_effect <- merge(EBVs_Spatial_effect, herd_effect_BHS, by= "herd", all.x = TRUE)


#Adding Permanent environmental effect from model BHS
perm_effect_BHS <- data.frame(fitBHS$summary.random$perm[,c("ID","mean")])
names(perm_effect_BHS)[1]<- "cow"
names(perm_effect_BHS)[2]<- "perm_effect"

#Merging herd_effect_BHS with EBVs_Spatial_effect
EBVs_Spatial_effect <- merge(EBVs_Spatial_effect, perm_effect_BHS, by= "cow", all.x = TRUE)

# Save EBvs and spatial effects
colnames(EBVs_Spatial_effect)
#"cow"                "herd"               "EBV_B"
#"EBV_BH"             "EBV_BS"             "EBV_BHS"
#"Spatial_effect_BS"  "Spatial_effect_BHS" "dEBV"
#"mean_milkZ"         "herd_effect"        "perm_effect"

names(EBVs_Spatial_effect)[1]<- "cow"

save(EBVs_Spatial_effect,file = "data/cleaned_data/EBVs_Spatial_effect.RData")

#-----PCA by spatial effect-----------------------------------
pcas_2 <- merge(x = pcas, y = EBVs_Spatial_effect, by = "cow", all.x = TRUE)
dim(pcas_2)
str(pcas_2)


# Visualize PCA by spatial effects_BHS
# PC1 vs PC2

spatial_pca1 <- ggplot(data = pcas_2) +
  geom_point(mapping = aes(x = `PC1`, y = `PC2`,color = Spatial_effect_BHS), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "PCA of SNP genotypes by spatial effect BHS",
       x = paste0("PC1(3.1%)"),
       y = paste0("PC2(0.8%)")) +
  theme_minimal()
spatial_pca1

# Visualize PCA by spatial EBV_BH
# PC1 vs PC2

EBV_BH_pca1 <- ggplot(data = pcas_2) +
  geom_point(mapping = aes(x = `PC1`, y = `PC2`,color = EBV_BH), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "PCA of SNP genotypes by EBV_BH",
       x = paste0("PC1(3.1%)"),
       y = paste0("PC2(0.8%)")) +
  theme_minimal()
EBV_BH_pca1


# Visualize PCA by spatial EBV_BHS
# PC1 vs PC2

EBV_BHS_pca1 <- ggplot(data = pcas_2) +
  geom_point(mapping = aes(x = `PC1`, y = `PC2`,color = EBV_BHS), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "PCA of SNP genotypes by EBV_BHS",
       x = paste0("PC1(3.1%)"),
       y = paste0("PC2(0.8%)")) +
  theme_minimal()
EBV_BHS_pca1


ggsave(plot = spatial_pca1 + PreseTheme, filename = "spatial_pca1_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = EBV_BHS_pca1 + PreseTheme, filename = "EBV_BHS_pca1_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = EBV_BH_pca1 + PreseTheme, filename = "EBV_BH_pca1_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation







# PC1 vs PC3
dgrp_pca2<- ggplot(data = pcas) +
  geom_point(mapping = aes(x = `PC1`, y = `PC3`,color = dgrp), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x = paste0("PC1(3.1%)"),
       y = paste0("PC3(0.5%)")) +
  theme_minimal()

# PC2 vs PC3
dgrp_pca3<- ggplot(data = pcas) +
  geom_point(mapping = aes(x = `PC2`, y = `PC3`,color = dgrp), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x = paste0("PC2(0.8%)"),
       y = paste0("PC3(0.5%)")) +
  theme_minimal()

dgrp_pca <- ggarrange(dgrp_pca1,dgrp_pca2,dgrp_pca3, ncol = 3, nrow = 1, common.legend = T, legend = "right", align = "h", widths = c(1,1,1))
dgrp_pca

ggsave(plot = dgrp_pca + PreseTheme, filename = "PCA_breed_proportion_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = dgrp_pca + PaperTheme, filename = "PCA_breed_proportion_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper

#--------MANOVA----------------------------------

with(pcas_3, plot(PC1 ~ PC2, col = c("blue", "red", "black", "yellow")[dgrp]))
fit <- manova(cbind(PC1, PC2, PC3) ~ dgrp, data = pcas_3)
summary(fit) # p=1.911e-10 --> highly significant!
summary.aov(fit)
# PC1: p=9.645e-11
# PC2: p=0.04549
# PC3: p=0.07128
fit <- lm(PC1 ~ dgrp, data = pcas_3)
summary(fit)
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  -15.194      3.193  -4.758 2.10e-06 ***
# dgrp2         18.850      4.459   4.228 2.47e-05 ***
# dgrp3         32.516      5.858   5.550 3.25e-08 ***
# dgrp4         50.191      9.947   5.046 4.95e-07 ***
fit <- lm(PC2 ~ dgrp, data = pcas_3)
summary(fit)
# --> nothing is really significant

with(pcas_3, plot(PC1 ~ PC2, col = c("blue", "red", "black", "yellow")[region]))
fit <- manova(cbind(PC1, PC2, PC3) ~ region, data = pcas_3)
summary(fit) # p=2.2e-16 --> highly significant
summary.aov(fit)
# PC1: p=2.2e-16
# PC2: p=2.2e-16
# PC3: p=6.568e-10
fit <- lm(PC1 ~ region, data = pcas_3)
summary(fit)
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)   32.539      4.090   7.956 3.03e-15 ***
#   region2    -46.455      5.115  -9.083  < 2e-16 ***
#   region3    -33.487      5.864  -5.711 1.30e-08 ***
#   region4    -43.890      6.633  -6.617 4.77e-11 ***
plot(y = data1$lat, x = data1$long, col = c("blue", "red", "black", "yellow")[data1$region])
plot(y = data2$lat, x = data2$long)
nrow(pcas_3)

fit <- lm(PC2 ~ region, data = pcas_3)
summary(fit)

fit <- lm(PC3 ~ region, data = pcas_3)
summary(fit)

tmp <- data1[!duplicated(data1[, c("cow")]), ]
with(tmp, table(dgrp, region))
summary(lm(as.numeric(dgrp) ~ region, data = tmp))

par(mfrow = c(2, 2))
lims <- range(pcas_3$PC1)
tmp <- pcas_3[pcas_3$dgrp == "1", ]
hist(tmp$PC1, xlim = lims, main = "1")
tmp <- pcas_3[pcas_3$dgrp == "2", ]
hist(tmp$PC1, xlim = lims, main = "2")
tmp <- pcas_3[pcas_3$dgrp == "3", ]
hist(tmp$PC1, xlim = lims, main = "3")
tmp <- pcas_3[pcas_3$dgrp == "4", ]
hist(tmp$PC1, xlim = lims, main = "4")

plot(y = pcas_3$EBV_BHS, x = pcas_3$Spatial_effect_BHS,
     col = c("blue", "red", "black", "yellow")[pcas_3$region])

lims <- range(c(pcas_3$EBV_B, pcas_3$EBV_BHS))
range(pcas_3$EBV_B)
range(pcas_3$EBV_BH)
range(pcas_3$EBV_BS)
range(pcas_3$EBV_BHS)
plot(y = pcas_3$EBV_B, x = pcas_3$EBV_BH,
     col = c("blue", "red", "black", "yellow")[pcas_3$region],
     pch = 19, cex = 0.1)
plot(y = pcas_3$EBV_B, x = pcas_3$EBV_BS,
     col = c("blue", "red", "black", "yellow")[pcas_3$region],
     pch = 19, cex = 0.1)
plot(y = pcas_3$EBV_BH, x = pcas_3$EBV_BS,
     col = c("blue", "red", "black", "yellow")[pcas_3$region],
     pch = 19, cex = 0.1)

summary(lm(EBV_B ~ region, data = pcas_3))
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.001682   0.003877  -0.434   0.6644
# region2     -0.001412   0.004849  -0.291   0.7709
# region3     -0.009324   0.005559  -1.677   0.0937 .
# region4      0.030381   0.006289   4.831 1.47e-06 ***

summary(lm(EBV_BH ~ region, data = pcas_3))
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.0013604  0.0018458  -0.737 0.461217
# region2      0.0006607  0.0023083   0.286 0.774718
# region3     -0.0015272  0.0026464  -0.577 0.563937
# region4      0.0100666  0.0029937   3.363 0.000788 ***

summary(lm(EBV_BS ~ region, data = pcas_3))
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.0001187  0.0002865  -0.414    0.679
# region2     -0.0001831  0.0003583  -0.511    0.609
# region3      0.0005166  0.0004108   1.258    0.209
# region4      0.0005951  0.0004647   1.280    0.201
# --> not significant anymore, good or bad? need validation

summary(lm(EBV_BHS ~ region, data = pcas_3))
# --> not significant anymore, good or bad? need validation

#-----------------Box plot dEBV vs Spatial effect BHS------------------

round(summary(EBVs_Spatial_Matrix$Spatial_effect_BHS),2)
str(EBVs_Spatial_Matrix$Spatial_effect_BHS)
# Duplicate Spatial_effect_BHS rename spatial_Box
EBVs_Spatial_Matrix$spatial_Box<- EBVs_Spatial_Matrix$Spatial_effect_BHS
summary(EBVs_Spatial_Matrix$spatial_Box)
#Min. 1st Qu.  Median   Mean 3rd Qu. Max.
# -1.15 -0.6 -0.23 -0.22  0.11        0.8
# I will use -1.15 -0.6 -0.23 0.11 0.8 for grouping (mean removed)

EBVs_Spatial_Matrix <- EBVs_Spatial_Matrix %>%
  mutate(group=
           case_when((EBVs_Spatial_Matrix$spatial_Box > -1.15 & EBVs_Spatial_Matrix$spatial_Box <= -0.6) ~"(-1.15,-0.6]",
                     (EBVs_Spatial_Matrix$spatial_Box > -0.6 & EBVs_Spatial_Matrix$spatial_Box<= -0.23) ~"(-0.6,-0.23]",
                     (EBVs_Spatial_Matrix$spatial_Box> -0.23  & EBVs_Spatial_Matrix$spatial_Box<= 0.11) ~"(-0.23,0.11]",
                     (EBVs_Spatial_Matrix$spatial_Box> 0.11  & EBVs_Spatial_Matrix$spatial_Box<= 0.8) ~"(0.11,0.8]"  ))

table(EBVs_Spatial_Matrix$group)
# Set the default theme to theme_classic() with the legend at the right of the plot:
theme_set(
  theme_classic() +
    theme(legend.position = "right")
)
# Change the default order of items
EBVs_Spatial_Matrix$group <- as.factor(EBVs_Spatial_Matrix$group)
#sd_EBV_Spde$group <- as.factor(sd_EBV_Spde$group)

Box <- ggplot(data=subset(EBVs_Spatial_Matrix, !is.na(group)), aes(x=reorder(group,dEBV), y=dEBV, fill=group)) +
  geom_boxplot() +
  theme(legend.position="none")
Box

# Changing axis names
Box<- Box + labs(x = "Spatial effect", y = "Difference in estimated breeding values (BH-BHS)")

ggsave(plot = Box + PaperThemeNoLegend, filename = "Boxplot_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper
table(EBVs_Spatial_Matrix$group)
#(-0.23,0.11] (-0.6,-0.23] (-1.15,-0.6]   (0.11,0.8]
#477          475          471          470
#----------------Scatter plot dEBV vs Spatial effect BHS------------------
library(devtools)
devtools::install_github("kassambara/easyGgplot2")
library(easyGgplot2)

scatter_plot<- ggplot(EBVs_Spatial_effect, aes(x=Spatial_effect_BHS, y=dEBV, colour=Spatial_effect_BHS)) +
  geom_point(size=0.5) +
  geom_smooth(method=lm) + labs(x="Spatial effects (BHS)", y="Differences between breeding values(BH-BHS)") + labs(color = "Spatial effects") +
  theme_bw()

scatter_plot

ggsave(plot = scatter_plot + PreseThemeLegendright, filename = "spatial_dEBV_scatter_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = scatter_plot + PaperThemeLegendright, filename = "spatial_dEBV_scatter_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper


scatter_plot_Paper<- ggplot(EBVs_Spatial_effect, aes(x=Spatial_effect_BHS, y=dEBV)) +
  geom_point(size=0.5) +
  geom_smooth(method=lm) + labs(x="Spatial effects (BHS)", y="Differences between breeding values(BH-BHS)") +
  theme_bw()

scatter_plot_Paper

ggsave(plot = scatter_plot_Paper + PreseThemeNoLegend, filename = "spatial_dEBV_scatter_Paper_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = scatter_plot_Paper + PaperThemeNoLegend, filename = "spatial_dEBV_scatter_Paper_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper

#------------Adding average cow milk yield percow----------------------------
scatter_plot_MY<- ggplot(EBVs_Spatial_effect, aes(x=Spatial_effect_BHS, y=mean_milkZ)) +
  geom_point(size=0.5) +
  geom_smooth(method=lm) + labs(x="Spatial effects (BHS)", y="Cow's average milk yield") +
  theme_bw()

scatter_plot_MY

ggsave(plot = scatter_plot_MY + PreseThemeLegendright, filename = "spatial_MY_scatter_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation



#------------Animal Ranking----------------------------------------------------

# Ranking EBV_B
EBV_B_Rank <- EBV_B[order(EBV_B$EBV_B,decreasing = TRUE),]

# Ranking EBV_fitBH
EBV_BH_Rank <- EBV_BH[order(EBV_BH$EBV_BH,decreasing = TRUE),]

# Ranking EBV_fitBS
EBV_BS_Rank <- EBV_BS[order(EBV_BS$EBV_BS,decreasing = TRUE),]

# Ranking EBV_fitBHS
EBV_BHS_Rank <- EBV_BHS[order(EBV_BHS$EBV_BHS,decreasing = TRUE),]

# Top 10 EBV_BH
EBV_BH_top10 <- EBV_BH_Rank[1:10,c("ID","EBV_BH")]
head(EBV_BH_top10, n=10)
ID    EBV_BH
#536 0.1429519 OK
#272 0.1350679 OK
#764 0.1328426 OK
#241 0.1304377 OK
#734 0.1244745 NO
#1818 0.1228799NO
#621 0.1205020  OK
#1670 0.1146873 NO
#750  1673 0.1120357 NO
#965  1867 0.1114881 NO
#NO BH:734,1818,1670,750,965
# Top 10 EBV_BHS
EBV_BHS_top10 <- EBV_BHS_Rank[1:10,c("ID","EBV_BHS")]
head(EBV_BHS_top10, n=10)
ID     EBV_BHS
#536 0.002140894 OK
#764 0.001890106 OK
#272 0.001639369 OK
#1867 0.001622924 NO
#241 0.001622191 OK
#621 0.001597591 OK
#607 0.001587894 NO
#561 0.001582671 NO
#1600 0.001535894 NO
#51 0.001530524 NO

#Comparing Top 10 cows (BH vs BHS) #5/10 overlaps
#NO BH:734,1818,1670,750,965
#NO BHS: 1867,607,561,1600,51
#Overlaps:536,764,272,241,621
ID_top_10_BH_BHS <- data.frame(c(734,1818,1670,750,965,1867,607,561,1600,51,536,764,272,241,621))
names(ID_top_10_BH_BHS)[1]<- "cow"

# Top 20 EBV_BH
EBV_BH_top20 <- EBV_BH_Rank[1:20,c("ID","EBV_BH")]
head(EBV_BH_top20, n=20)

ID     EBV_BH
#536 0.14295193  OK
#272 0.13506792  OK
#764 0.13284255  OK
#241 0.13043770  OK
#734 0.12447452  OK
#1818 0.12287991 NO
#621 0.12050204  OK
#1670 0.11468730 OK
#1673 0.11203567 NO
#1867 0.11148815 OK
#1600 0.11003804 OK
#955 0.10859908  OK
#738 0.10821744  NO
#52 0.10283803   OK
#828 0.10263070 NO
#394 0.10203878 NO
#607 0.10104092 OK
#561 0.09966824 OK
#1766 0.09770793 NO
#1765 0.09762628 OK

# Top 20 EBV_BHS
EBV_BHS_top20 <- EBV_BHS_Rank[1:20,c("ID","EBV_BHS")]
head(EBV_BHS_top20, n=20)
ID     EBV_BHS
#536 0.002140894 OK
#764 0.001890106 OK
#272 0.001639369 OK
#1867 0.001622924 OK
#241 0.001622191 OK
#621 0.001597591 OK
#607 0.001587894 OK
#561 0.001582671 OK
#1600 0.001535894 OK
#51 0.001530524 NO
#734 0.001500917 OK
#955 0.001494856 OK
#994 0.001432891  NO
#1148 0.001417028 NO
#560 0.001408543 NO
#1670 0.001371440 OK
#1419 0.001367967 NO
#889 0.001366775 NO
#1765 0.001359953 OK
#52 0.001312502 OK

14 overlaps/20

#---------------Schematic illustration of spatial effect------------------------
# PCA by phenotype, EBV_BHS, herd_effect, spatial effect and Permanent effect_top10-
pcas_3 <- merge(x = pcas, y = EBVs_Spatial_effect, by = "cow", all.x = TRUE)
dim(pcas_3)
str(pcas_3)
# PCA of top 10 ranked cows from BH and BHS (5 overlaps and 5 mismatches resulting in 15 animals in total)
sel <-pcas_3$cow %in% ID_top_10_BH_BHS$cow
pcas_top_10 <- pcas_3[sel,]

# PCA by phenotype

PCA_phenotype <-ggplot(data = pcas_top_10 ) +
  geom_point(mapping = aes(x = `PC1`, y = `PC2`,color = mean_milkZ), show.legend = TRUE ) +
  labs(title = "Phenotype",
       x = paste0(""),
       y = paste0("")) + labs(color="") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

PCA_phenotype

# PCA by breeding values
PCA_EBV_BHS <- ggplot(data=pcas_top_10, aes(PC1,PC2,color = EBV_BHS)) +
  geom_point() + labs(title = "Genetics",
       x = paste0(""),
       y = paste0("")) + labs(color="") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
PCA_EBV_BHS
# PCA by herd efect
PCA_herd_effect <- ggplot(data=pcas_top_10, aes(PC1,PC2,color = herd_effect)) +
  geom_point() + labs(title = "Herd",
                      x = paste0(""),
                      y = paste0("")) + labs(color="") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
PCA_herd_effect

# PCA by spatial effect
PCA_spatial_effect <- ggplot(data=pcas_top_10, aes(PC1,PC2,color = Spatial_effect_BHS)) +
  geom_point() + labs(title = "Spatial",
                      x = paste0(""),
                      y = paste0("")) + labs(color="") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
PCA_spatial_effect


# PCA by permanent environmental effect
PCA_perm_effect <- ggplot(data=pcas_top_10, aes(PC1,PC2,color = perm_effect)) +
  geom_point() + labs(title = "Permanent",
                      x = paste0(""),
                      y = paste0("")) + labs(color="") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
PCA_perm_effect


schematic_pca <- ggarrange(PCA_phenotype,PCA_EBV_BHS,PCA_herd_effect,PCA_spatial_effect,PCA_perm_effect, ncol = 3, nrow = 2, common.legend = F, legend = "right", align = "h", widths = c(1,1,1))
schematic_pca

ggsave(plot = schematic_pca + PreseThemeLegendright, filename = "schematic_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation

ggsave(plot = schematic_pca + PaperThemeLegendright, filename = "schematic_paper.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation

#--PCA by phenotype, EBV_BHS, herd_effect, spatial effect and Permanent effect_top10-missmatches-

# PCA of top 10 ranked_mismatches (10) cows from BH and BHS (5 overlaps and 5 mismatches resulting in 15 animals in total)
sel <-pcas_3$cow %in% ID_top_10_BH_BHS[1:10,]
pcas_top_10_miss <- pcas_3[sel,]

# PCA by phenotype

PCA_phenotype <-ggplot(data = pcas_top_10_miss ) +
  geom_point(mapping = aes(x = `PC1`, y = `PC2`,color = mean_milkZ), show.legend = TRUE ) +
  labs(title = "Phenotype",
       x = paste0(""),
       y = paste0("")) + labs(color="") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

PCA_phenotype

# PCA by breeding values
PCA_EBV_BHS <- ggplot(data=pcas_top_10_miss, aes(PC1,PC2,color = EBV_BHS)) +
  geom_point() + labs(title = "Genetics",
                      x = paste0(""),
                      y = paste0("")) + labs(color="") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
PCA_EBV_BHS
# PCA by herd efect
PCA_herd_effect <- ggplot(data=pcas_top_10_miss, aes(PC1,PC2,color = herd_effect)) +
  geom_point() + labs(title = "Herd",
                      x = paste0(""),
                      y = paste0("")) + labs(color="") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
PCA_herd_effect

# PCA by spatial effect
PCA_spatial_effect <- ggplot(data=pcas_top_10_miss, aes(PC1,PC2,color = Spatial_effect_BHS)) +
  geom_point() + labs(title = "Spatial",
                      x = paste0(""),
                      y = paste0("")) + labs(color="") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
PCA_spatial_effect


# PCA by permanent environmental effect
PCA_perm_effect <- ggplot(data=pcas_top_10_miss, aes(PC1,PC2,color = perm_effect)) +
  geom_point() + labs(title = "Permanent",
                      x = paste0(""),
                      y = paste0("")) + labs(color="") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
PCA_perm_effect


schematic_pca_miss <- ggarrange(PCA_phenotype,PCA_EBV_BHS,PCA_herd_effect,PCA_spatial_effect,PCA_perm_effect, ncol = 3, nrow = 2, common.legend = F, legend = "right", align = "h", widths = c(1,1,1))
schematic_pca_miss

ggsave(plot = schematic_pca_miss + PreseThemeLegendright, filename = "schematic_miss_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation

ggsave(plot = schematic_pca_miss + PaperThemeLegendright, filename = "schematic_miss_Paper.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation

#-----Scatter_plot_Matrix-------------------------------------------------
colnames(pcas_top_10)
Scatter_Matrix <- pairs.panels(EBVs_Spatial_effect[,c("mean_milkZ","EBV_BHS","herd_effect","Spatial_effect_BHS", "perm_effect")], main = "Scatter Plot Matrix")

ggsave(plot = Scatter_Matrix + PreseThemeLegendright, filename = "Scatter_Matrix_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation
