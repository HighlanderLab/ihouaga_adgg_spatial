# ---- Header ------------------------------------------------------------------
# Spatial modelling
# Trait: Test day milk yield
# Author: Isidore Houaga, Ivan Pocrnic and Gregor Gorjanc.
# Version 1.0.0
# Date: 2023-07-20
# ---- Setup--------------------------------------------------------------------

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
# Subset of data1
#data1_1_2 <- subset(data1, lac=="1" | lac=="2")

#rm(nb.map,nb.matrix, nb.matrixScaled)
# Now we can code the cows in pheno data correctly
data1$cow <- factor(data1$cow, levels = 1:nrow(GRMInv))
summary(as.numeric(levels(data1$cow))) # we need 1:1894 numbers here!
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

#--------Building mesh and prepare data for SPDE modelling----------------------
# Create dataframe of coordinates
Mydata_gps <- data1[, c("long", "lat")]
#save(Mydata_gps, file = "Mydata_gps.RData")

dat <- sf::st_as_sf(Mydata_gps, coords = c("long", "lat"), crs = 4326)
crs_tz <- fm_crs(21035)
fm_length_unit(crs_tz) <- "km"
dat <- fm_transform(dat, crs_tz)

ggplot() + geom_sf(data = dat)
ggplot() + geom_sf(data = dat) + coord_sf(datum = crs_tz)
ggplot() + geom_sf(data = dat) + coord_sf(datum = crs_tz) +
  geom_hline(yintercept = 9500) +
  geom_vline(xintercept = 1500) +
  geom_vline(xintercept = 1300) +
  geom_hline(yintercept = 9200)

bnd_inner <- fm_nonconvex_hull_inla(dat, 50)

# Need to extend the mesh a bit beyond the artificial country borders to avoid
# edge and "island" effects:
# Load Tanzania map (level 1: regional boundaries)
#mapTZA <- geodata::gadm(country = "TZA", level = 1, path = tempdir())

#Alternative to load map data

# Retrieve all countries as an sf object
countries <- ne_countries(returnclass = "sf")
str(countries)
# Filter for Tanzania
mapTZA <- countries[countries$admin == "United Republic of Tanzania", ]

# Plot Tanzania boundary
plot(st_geometry(mapTZA), main = "Tanzania Boundary")

# Inspect mapTZA
class(mapTZA) # "SpatVector"

# Convert mapTZA from "SpatVector" to sf format
#mapTZA<- st_as_sf(mapTZA)
#class(mapTZA) # "sf"         "data.frame"

#Inspect the CRS of the mapTZA
st_crs(mapTZA) #  the CRS is (longlat: EPSG:4326) ie WGS84

data_border <-  st_union(mapTZA) # whatever code needed to load the border data and sets its original crs
plot(data_border)
data_border <- fm_transform(data_border, crs_tz)

dat_border_points <- fm_as_segm(data_border)$loc

bnd_outer <- fm_nonconvex_hull_inla(dat_border_points, 100, convex = 325)
plot(data_border)
plot(bnd_inner)
plot(bnd_outer)
mesh <- fm_mesh_2d_inla(boundary = list(bnd_inner, bnd_outer),
                        max.edge = c(10, 40),
                        cutoff = 10,
                        crs = crs_tz)

mesh2 <- fm_mesh_2d_inla(boundary = list(bnd_inner, data_border),
                        max.edge = c(10, 40),
                        cutoff = 10,
                        crs = crs_tz)
#save(mesh, file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/mesh.RData")
#save(mapTZA, file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/mapTZA.RData")
#save(data_border, file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/data_border.RData")
#save(bnd_inner, file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/bnd_inner.RData")
#save(bnd_outer, file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/bnd_outer.RData")

meshplot<- ggplot() +
  geom_fm(data = mesh) +
  geom_sf(data = dat, size = 0.1) + geom_sf(data=data_border,alpha=0.2) +
  coord_sf(datum = crs_tz) +
  labs(x = "Longitude (km)", y = "Latitude (km)")
meshplot

#meshplot2<- ggplot() +
  geom_fm(data = mesh2) +
  geom_sf(data = dat, size = 0.1) + geom_sf(data=data_border,alpha=0.2) +
  coord_sf(datum = crs_tz) +
  labs(x = "Longitude (km)", y = "Latitude (km)")
#meshplot2
#--------Formating of plot------------------------------------------------------
PaperTheme = theme_bw(base_size = 12, base_family = "serif") +
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
#-------------------------------------------------------------------------------

# Export the mesh in Paper_plot format.
ggsave(plot = meshplot + PaperTheme, filename = "Figure1New.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper

#For INLABru define matern
matern <-
  inla.spde2.pcmatern(mesh,
                      prior.sigma = c(sqrt(0.25), 0.5),
                      prior.range = c(50, 0.8)
  )

#save(matern, file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/matern.RData")
matern2 <-
  inla.spde2.pcmatern(mesh2,
                      prior.sigma = c(sqrt(0.25), 0.5),
                      prior.range = c(50, 0.8)
  )

save(matern2, file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/matern2.RData")


#data2 <- sf::st_as_sf(x = locations[, c("long", "lat")],
                      #coords = c("long", "lat")) # X and Y?

data2 <- dat[,"geometry"]

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
data2$ward_raph <- data1$ward
data2$lac <- data1$lac
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
data2$region <- data1$region
# Create age groups
data2$age_group <- cut(data2$age, breaks = c(0, 36, 60, Inf), labels = c("18-36 months", "36-60 months", "> 60 months"))
colnames(data2)
length(unique(data2$geometry)) #1385 couples of GPS vs 1386 herds (2 herds at same location)
# Create HYS (Herd year season)
data2$hys <- with(data2, paste0(herd, "-", cyrsn))
table(table(data2$hys))
# Create a spatialpointdataframe needed for special effect prediction.
# Following instructions from  https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
# prepare coordinates, data, and proj4string
#coords <- data1[ , c("long", "lat")]   # coordinates
#data   <- data1          # data
#crs    <- crs_tz # proj4string of coords. This was a wrong one.
#28992 is for Netherlands.
#crs   <- ??????  Use local system
#coords <- locations[ , c("X", "Y")]   # coordinates
#data   <- data1          # data
#crs    <- crs_tz # proj4string of coords. This was a wrong one.

class(dat)

# make the SpatialPointsDataFrame object
#spdf <- SpatialPointsDataFrame(coords      = coords,
                              # data        = data,
                               #proj4string = crs)

spdf <- as(dat, "Spatial")
class(spdf) # "SpatialPointsDataFrame"

#------Specify models R-INLAbru-------------------------------------------------
#-------Define G models G (without permanent environmental effect)--------------
# ModelG
modelG <- ~ fixed_effects(main = ~ tyrmn + cyrsn + dgrp + ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv)
modelGFormula <- milkZ ~ .

#ModelGH
modelGH <- ~ fixed_effects(main = ~ tyrmn + cyrsn + dgrp + ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  herd(main=herd, model = 'iid')
modelGHFormula <- milkZ ~ .

# ModelGS
modelGS <- ~ fixed_effects(main = ~ tyrmn + cyrsn + dgrp + ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  field(geometry, model = matern)

modelGSFormula <- milkZ ~ .

# ModelGS2 # Spatial plot in paper-Supplementary
modelGS2 <- ~ fixed_effects(main = ~ tyrmn + cyrsn + dgrp + ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  field(geometry, model = matern2)

modelGSFormula2 <- milkZ ~ .

# ModelGHS
modelGHS <- ~ fixed_effects(main = ~ tyrmn + cyrsn + dgrp + ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  herd(main=herd, model = 'iid') +
  field(geometry, model = matern)
modelGHSFormula <- milkZ ~ .


# ModelGHS2
modelGHS2 <- ~ fixed_effects(main = ~ tyrmn + cyrsn + dgrp + ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  herd(main=herd, model = 'iid') +
  field(geometry, model = matern2)
modelGHSFormula2 <- milkZ ~ .

#-------Define GP models with permanent environmental effect--------------------
#---------------------Define GP models------------------------------------------
# ModelGP
modelGP <- ~ fixed_effects(main = ~ tyrmn + cyrsn + dgrp + ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  perm(main = cowPe, model = 'iid')
modelGPFormula <- milkZ ~ .

# ModelGPH
modelGPH <- ~ fixed_effects(main = ~ tyrmn + cyrsn + dgrp +  ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  perm(main = cowPe, model = 'iid') +
  herd(main=herd, model = 'iid')

modelGPHFormula <- milkZ ~ .


# ModelGPS
modelGPS <- ~ fixed_effects(main = ~ tyrmn + cyrsn + dgrp + ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  perm(main = cowPe, model = 'iid') +
  field(geometry, model = matern)
modelGPSFormula <- milkZ ~ .

# ModelGPHS
modelGPHS <- ~ fixed_effects(main = ~ tyrmn + cyrsn + dgrp + ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  perm(main = cowPe, model = 'iid') +
  herd(main=herd, model = 'iid') +
  field(geometry, model = matern)

modelGPHSFormula <- milkZ ~ .

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


#---Run the models--------------------------------------------------------------
#-----------Run GP models-----------------------------------------------------
# Run fitGP
fitGP <- bru(modelGP,
             like(family = "Gaussian",
                  modelGPFormula,
                  data = data2))

summary(fitGP) #32692.42
summarise_precision_to_variance(fitGP)

# Run fitGPH
fitGPH <- bru(modelGPH,
             like(family = "Gaussian",
                  modelGPHFormula,
                  data = data2))

summary(fitGPH) #DIC: 32633.02
summarise_precision_to_variance(fitGPH)


# Run fitGPS
fitGPS <- bru(modelGPS,
             like(family = "Gaussian",
                  modelGPSFormula,
                  data = data2))

summary(fitGPS) #DIC:32620.02
summarise_precision_to_variance(fitGPS)

# Predict spatial effects at our locations and calculate Spatial variance for GPS
#fitGPS_sample <- inla.posterior.sample(1000, fitGPS) # Sampling posterior means

fieldGPS <- predict(fitGPS,spdf, ~field)
fieldGPS_df <- data.frame(fieldGPS[,c("cow" , "mean")])
fieldGPS_df <- distinct(data.frame(fieldGPS_df))

# Run fitGPHS
fitGPHS <- bru(modelGPHS,
               like(family = "Gaussian",
                    modelGPHSFormula,
                    data = data2))

summary(fitGPHS) #DIC:  32557.85 vs old 32554.88
summarise_precision_to_variance(fitGPHS)

#--------Run G models-----------------------------------------------------------
# Run fitG
fitG<- bru(modelG,
           like(family = "Gaussian",
                modelGFormula,
                data = data2))

summary(fitG) # 32715.99
summarise_precision_to_variance(fitG)


# Run fitGH
fitGH<- bru(modelGH,
            like(family = "Gaussian",
                 modelGHFormula,
                 data = data2))

summary(fitGH) #
summarise_precision_to_variance(fitGH)


# Run fitGS
fitGS<- bru(modelGS,
            like(family = "Gaussian",
                 modelGSFormula,
                 data = data2))

summary(fitGS) # 32735.66
summarise_precision_to_variance(fitGS)

# Run fitGS2 # for the plot in supplementary
fitGS2<- bru(modelGS2,
            like(family = "Gaussian",
                 modelGSFormula2,
                 data = data2))

summary(fitGS2) # 32736.75
summarise_precision_to_variance(fitGS2)



# Run fitGHS
fitGHS<- bru(modelGHS,
             like(family = "Gaussian",
                  modelGHSFormula,
                  data = data2))

summary(fitGHS) # 32596.25
summarise_precision_to_variance(fitGHS)

# Run fitGHS2 with mesh2
fitGHS2<- bru(modelGHS2,
             like(family = "Gaussian",
                  modelGHSFormula2,
                  data = data2))

summary(fitGHS2) # 32597.73
summarise_precision_to_variance(fitGHS2)

#--------Save the models on External drive-------------------------------------
save(fitG,file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/Gmodels/fitG.RData")
save(fitGH,file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/Gmodels/fitGH.RData")
save(fitGS,file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/Gmodels/fitGS.RData")
save(fitGHS,file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/Gmodels/fitGHS.RData")
save(fitGS2,file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/Gmodels/fitGS2.RData")
save(fitGHS2,file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/Gmodels/fitGHS2.RData")
save(fitGP,file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/GPmodels/fitGP.RData")
save(fitGPH,file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/GPmodels/fitGPH.RData")
save(fitGPS,file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/GPmodels/fitGPS.RData")
save(fitGPHS,file = "/Users/ihouaga2/OneDrive - University of Edinburgh/Spatial_results/GPmodels/fitGPHS.RData")

#------------------------ Variance components Paper-----------------------------
#GP
summary(fitGP)
summarise_precision_to_variance(fitGP)
#GPH
summary(fitGPH)
summarise_precision_to_variance(fitGPH)
#GPS
summary(fitGPS)
summarise_precision_to_variance(fitGPS)

#--------------------Spatial variance-------------------------------------------
#td_dev1 = 7.24e-01
#std_dev_sd1 = 8.00e-02
#variance1 = std_dev1 ** 2
#variance_sd1 = 2 * std_dev1 * std_dev_sd1
#variance1
#variance_sd1

#GP models
summary(fitGPS)
summarise_precision_to_variance(fitGPS)
#Spatial variance
# Set the seed for reproducibility
set.seed(123)
fieldGPS <- generate(fitGPS,spdf, ~field, 1000)
str(fieldGPS)
# Reporting Spatial variance for our herd (locations)
fieldGPS <- as.data.frame(fieldGPS)
fieldGPS$herd <- data1$herd
fieldGPS_herd <- fieldGPS[, c(ncol(fieldGPS), 1:(ncol(fieldGPS) - 1))]
fieldGPS_herd <- distinct(fieldGPS_herd)
fieldGPS_herd <- fieldGPS_herd[, -1]
column_variances_GPS <- data.frame(apply(fieldGPS_herd, 2, var)) #1 means row; 2 means cols
colnames(column_variances_GPS)[1] <- "spatial_var"
round(mean(column_variances_GPS$spatial_var),3) # 0.34
round(sd(column_variances_GPS$spatial_var),3) # 0.07

# To find out which elements are directly available to use in the prediction/formula expression in generate()/predict(), call
#generate(the_fit, n.samples = 1)
#and look at the names in the output. This should include all of the covariance parameters.
#fieldGS1 <- generate(fitGS, n.samples = 1)

Range_fieldGPS <- generate(fitGPS,spdf, ~Range_for_field, 1000)
str(Range_fieldGPS)
Range_fieldGPS_df <- data.frame(Range_fieldGPS)
round(rowMeans(Range_fieldGPS_df),3) #33.5
# Define a function to calculate row-wise standard deviation
row_sd <- function(data) {
  apply(data, 1, sd, na.rm = TRUE)
}
sd_range_GPS <- round(row_sd(Range_fieldGPS_df),3)
sd_range_GPS # 6.6

#GPHS
summary(fitGPHS)
summarise_precision_to_variance(fitGPHS)
# Spatial variance

# Reporting Spatial variance for fitGPHS at our herd (locations)
fieldGPHS <- generate(fitGPHS,spdf, ~field, 1000)
fieldGPHS <- as.data.frame(fieldGPHS)
fieldGPHS$herd <- data1$herd
fieldGPHS_herd <- fieldGPHS[, c(ncol(fieldGPHS), 1:(ncol(fieldGPHS) - 1))]
fieldGPHS_herd <- distinct(fieldGPHS_herd)
fieldGPHS_herd <- fieldGPHS_herd[, -1]
column_variances_GPHS <- data.frame(apply(fieldGPHS_herd, 2, var)) #1 means row; 2 means cols
colnames(column_variances_GPHS)[1] <- "spatial_var"
round(mean(column_variances_GPHS$spatial_var),3) #0.33
round(sd(column_variances_GPHS$spatial_var),3) # 0.06


#Range of spatial effect fitGPHS
Range_fieldGPHS <- generate(fitGPHS,spdf, ~Range_for_field, 1000)

Range_fieldGPHS_df <- data.frame(Range_fieldGPHS)
round(rowMeans(Range_fieldGPHS_df),3) #32.3
sd_range_GPHS <- round(row_sd(Range_fieldGPHS_df),3)
sd_range_GPHS # 7.8


#GS
summary(fitGS)
summarise_precision_to_variance(fitGS)
#Spatial variance
# Set the seed for reproducibility
set.seed(123)
fieldGS <- generate(fitGS,spdf, ~field, 1000)
str(fieldGS)
# Reporting Spatial variance for our herd (locations)
fieldGS <- as.data.frame(fieldGS)
fieldGS$herd <- data1$herd
fieldGS_herd <- fieldGS[, c(ncol(fieldGS), 1:(ncol(fieldGS) - 1))]
fieldGS_herd <- distinct(fieldGS_herd)
fieldGS_herd <- fieldGS_herd[, -1]
column_variances_GS <- data.frame(apply(fieldGS_herd, 2, var)) #1 means row; 2 means cols
colnames(column_variances_GS)[1] <- "spatial_var"
round(mean(column_variances_GS$spatial_var),3) # 0.32
round(sd(column_variances_GS$spatial_var),3) # 0.06

# To find out which elements are directly available to use in the prediction/formula expression in generate()/predict(), call
#generate(the_fit, n.samples = 1)
#and look at the names in the output. This should include all of the covariance parameters.
#fieldGS1 <- generate(fitGS, n.samples = 1)

Range_fieldGS <- generate(fitGS,spdf, ~Range_for_field, 1000)
str(Range_fieldGS)
Range_fieldGS_df <- data.frame(Range_fieldGS)
round(rowMeans(Range_fieldGS_df),3) #32.4
# Define a function to calculate row-wise standard deviation
row_sd <- function(data) {
  apply(data, 1, sd, na.rm = TRUE)
}
sd_range_GS <- round(row_sd(Range_fieldGS_df),3)
sd_range_GS #7.5


#GHS
summary(fitGHS)
summarise_precision_to_variance(fitGHS)
# Spatial variance

# Reporting Spatial variance for fitGHS at our herd (locations)
fieldGHS <- generate(fitGHS,spdf, ~field, 1000)
fieldGHS <- as.data.frame(fieldGHS)
fieldGHS$herd <- data1$herd
fieldGHS_herd <- fieldGHS[, c(ncol(fieldGHS), 1:(ncol(fieldGHS) - 1))]
fieldGHS_herd <- distinct(fieldGHS_herd)
fieldGHS_herd <- fieldGHS_herd[, -1]
column_variances_GHS <- data.frame(apply(fieldGHS_herd, 2, var)) #1 means row; 2 means cols
colnames(column_variances_GHS)[1] <- "spatial_var"
round(mean(column_variances_GHS$spatial_var),3) # 0.321
round(sd(column_variances_GHS$spatial_var),3) #  0.054


#Range of spatial effect fitGHS
Range_fieldGHS <- generate(fitGHS,spdf, ~Range_for_field, 1000)
Range_fieldGHS_df <- data.frame(Range_fieldGHS)
round(rowMeans(Range_fieldGHS_df),3) #33.4
sd_range_GHS <- round(row_sd(Range_fieldGHS_df),3)
sd_range_GHS #8.4


#----------------Check distribution of phenotypes-------------------------------

#Distribution of ages by lactation group - for example look at
#histogram of ages by lactation group and
#mean of yield by age and lactation group


# Histogram of age by lactation group

ggplot(data2, aes(x = age)) +
  geom_histogram(binwidth = 1, color = "black", fill = "skyblue", alpha = 0.7) +
  facet_wrap(~ lacgr, scales = "free") +
  labs(title = "Histogram of Age by Lactation Group",
       x = "Age",
       y = "Frequency")


# Histogram of milk yield by lactation group
ggplot(data2, aes(x = milk)) +
  geom_histogram(binwidth = 1, color = "black", fill = "skyblue", alpha = 0.7) +
  facet_wrap(~ lacgr, scales = "free") +
  labs(title = "Histogram of milk yield by Lactation Group",
       x = "Milk yield",
       y = "Frequency")


# Histogram of milk yield by age group

ggplot(data2, aes(x = milk)) +
  geom_histogram(binwidth = 1, color = "black", fill = "skyblue", alpha = 0.7) +
  facet_wrap(~ age_group, scales = "free") +
  labs(title = "Histogram of milk yield by age Group",
       x = "Milk yield",
       y = "Frequency")

# Mean of milk yield by lactation group
plot(data2$milk~data2$lacgr,
     main = "  Mean of milk yield by lactation group",
     xlab = "Lactation group",
     ylab = "Milk yield")
tmp <- tapply(X=data2$milk, INDEX = data2$lacgr, FUN=mean)
points(tmp~ names(tmp), pch=19, col="red", cex=2)

tapply(X=data2$milk, INDEX = data2$lacgr, FUN=mean)

# 1                2             3
# 7.383087    8.267883         9.233086

# Scatter plot of milk yield by lactation group
ggplot(data = data2, aes(x = lacgr, y = milk)) +
  geom_point() +  # Scatter plot
  labs(title = "Scatter Plot of Milk Yield by lactation group", x = "Lactation Group", y = "Milk Yield")



# Mean of milk yield by age_group

plot(data2$milk~data2$age_group,
     main = "  Mean of milk yield by age group",
     xlab = "Age group",
     ylab = "Milk yield")
tmp <- tapply(X=data2$milk, INDEX = data2$age_group, FUN=mean)
points(tmp~ names(tmp), pch=19, col="red", cex=2)

tapply(X=data2$milk, INDEX = data2$age_group, FUN=mean)

# 18-36 months  36-60 months    > 60 months
# 7.979894      8.080243           9.204795


# Scatter plot of milk yield by age
P<- ggplot(data = data2, aes(x = age, y = milk)) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) + # Scatter plot + smooth line
  labs(title = "Scatter Plot of Milk Yield by Age across parities", x = "Age", y = "Milk Yield")
P

# Scatter plot of milk yield by age in each parity
table(data1$lac)
data1_1 <- subset(data1, lac=="1" )

# Scatter plot of milk yield by age
P1<- ggplot(data = data1_1, aes(x = age, y = milk)) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) + # Scatter plot + smooth line
  labs(title = "Scatter Plot of Milk Yield by Age in Parity 1", x = "Age", y = "Milk Yield")
P1
data1_2 <- subset(data1, lac=="2" )
P2<- ggplot(data = data1_2, aes(x = age, y = milk)) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) + # Scatter plot + smooth line
  labs(title = "Scatter Plot of Milk Yield by Age in Parity 2", x = "Age", y = "Milk Yield")
P2
data1_3 <- subset(data1, lac=="3" )
P3<- ggplot(data = data1_3, aes(x = age, y = milk)) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) + # Scatter plot + smooth line
  labs(title = "Scatter Plot of Milk Yield by Age in Parity 3", x = "Age", y = "Milk Yield")
P3

data1_4 <- subset(data1, lac=="4" )
P4<- ggplot(data = data1_4, aes(x = age, y = milk)) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) + # Scatter plot + smooth line
  labs(title = "Scatter Plot of Milk Yield by Age in Parity 4", x = "Age", y = "Milk Yield")
P4

data1_5 <- subset(data1, lac=="5" )
P5<- ggplot(data = data1_5, aes(x = age, y = milk)) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) + # Scatter plot + smooth line
  labs(title = "Scatter Plot of Milk Yield by Age in Parity 5", x = "Age", y = "Milk Yield")
P5


data1_6 <- subset(data1, lac=="6" )
P6<- ggplot(data = data1_6, aes(x = age, y = milk)) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) + # Scatter plot + smooth line
  labs(title = "Scatter Plot of Milk Yield by Age in Parity 6", x = "Age", y = "Milk Yield")
P6

data1_7 <- subset(data1, lac=="7" )
P7<- ggplot(data = data1_7, aes(x = age, y = milk)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", se = FALSE) + # Scatter plot + smooth line
  labs(title = "Scatter Plot of Milk Yield by Age in Parity 7", x = "Age", y = "Milk Yield")
P7
data1_8 <- subset(data1, lac=="8" )

P8<- ggplot(data = data1_8, aes(x = age, y = milk)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", se = FALSE) + # Scatter plot + smooth line
  labs(title = "Scatter Plot of Milk Yield by Age in Parity 8", x = "Age", y = "Milk Yield")
P8

(P1_8 <- ggarrange(P1,P2,P3,P4,P5,P6,P7,P8, ncol = 2, nrow = 4, common.legend = T, legend = "left", align = "h", widths = c(1,1,1)))
# Scatter plot of milk yield by age group
ggplot(data = data2, aes(x = age_group, y = milk)) +
  geom_point() +  # Scatter plot
  labs(title = "Scatter Plot of Milk Yield by Age group", x = "Age Group", y = "Milk Yield")

# Look at how repeated is our data
#library(dplyr)
data2_by_cow <- data2 %>%
 group_by(data2$cow) %>%
 summarise(Total_Count = n())

summary(data2_by_cow$Total_Count)

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

#Summary statistics for the paper--------------------------------------------
# Average Milk yield
round(mean(data1$milk),2)
#8.32
round(sd(data1$milk),2)
# 4.31

#Milk yield by breed proportion

round(tapply(data1$milk, data1$dgrp,mean),2)
#       1        2        3        4
# 9.61 8.24 6.26 4.51

round(tapply(data1$milk, data1$dgrp,sd),2)
#1    2    3    4
#4.24 4.24 3.59 2.64

#Milk yield by region
round(tapply(data1$milk, data1$region,mean),2)
# 1     2     3     4
# 6.70  8.04  8.00 12.41

round(tapply(data1$milk, data1$region,sd),2)
#1    2    3    4
# 3.93 3.95 3.63 4.40


#Milk yield by parity group
round(tapply(data1$milk, data1$lacgr,mean),2)
#1    2    3
#7.38 8.27 9.23

round(tapply(data1$milk, data1$lacgr,sd),2)

#1    2    3
#4.12 4.19 4.42

# Number of records
dim(data1)
#19375    29
table(data1$dgrp)
#1    2    3    4
#7466 8149 3058  702

table(data1$lacgr)
#1    2    3
#5883 7018 6474

table(data1$region)
#1    2    3    4
#4547 8482 3674 2672

# Number of cows
#Creating data for each parity group
data2_P1 <- subset(data2, lacgr=="1")# cows
data2_P2 <- subset(data2, lacgr=="2")# cows
data2_P3<- subset(data2, lacgr=="3")# cows

length(unique(data2_P1$cow)) # 1024
length(unique(data2_P2$cow)) #1206
length(unique(data2_P3$cow)) # 777

#Create box plots of Milk Yield variation by breed proportion across regions
# Define new labels for regions and dgrp
region_levels <- c(2, 1, 4, 3)  # Custom order: North-Central first
region_labels <- c("North-Central", "North-East", "South-West", "South-Central")
dgrp_labels <- c("[100, 87.5]", "(87.5, 60]", "(60, 36]", "(36, 0]")

# Create the box plots

p <- ggplot(data1, aes(x = factor(dgrp), y = milk, fill = factor(dgrp))) +
  geom_boxplot(outlier.size = 0.3, size=0.3) +
  facet_wrap(~ factor(region, levels = region_levels, labels = region_labels)) +
  scale_x_discrete(labels = dgrp_labels) +
  labs(title = "", x = "Percentage of exotic genes (%)", y = "Milk Yield (L)") +
  theme_minimal()

# Title: Boxplot of milk yield  by percentage of exotic genes for each region

# Customize the plot
p <- p +
  theme(
    legend.position = "none",  # Remove legend
    strip.text = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 16),  # Set axis title font size
    axis.text = element_text(face = "bold", size = 16) # Make region labels bold
  )

p
ggsave(plot = p + PaperThemeNoLegend, filename = "Figure2_MilkBoxplot.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper


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
# Plot genomic variance across models

Pg<-ggplot()  +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitB_inla$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="B")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBH_inla$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="BH")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBS_inla$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="BS"), linetype = "dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBHS_inla$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="BHS"), linetype = "dashed") +
  scale_color_manual(values = c("B" = "black",
                                "BH"="red",
                                "BS"="blue",
                                "BHS"="green")) +
  labs(x="", y="", title ="Genomic variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model"))

# scale_y_continuous(breaks=NULL) to remove axis labels
Pg<- Pg+ theme(plot.title = element_text(face = "bold"))
Pg

# Individual plot

#PgEB<- ggplot()  +
# geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                   #fitB_inla$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="B")) +
  # scale_color_manual(values = c("B" = "black",
  #                             "BH"="red",
  #                              "BS"="blue",
  #                             "BHS"="green")) + labs(x="", y="", title ="Genomic variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  #theme(legend.position = "left",
  #   legend.box.background = element_blank(),
  #     legend.title = element_text("Models")) +
  #guides(colour=guide_legend(title="Model")) + theme_bw()

  #PgEB

#PgEBH<- ggplot()  +
#geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBH_inla$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="BH")) +
  #scale_color_manual(values = c("B" = "black",
                                "BH"="red",
                                "BS"="blue",
                                "BHS"="green")) + labs(x="", y="", title ="Genomic variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  #theme(legend.position = "left",
  #    legend.box.background = element_blank(),
  #    legend.title = element_text("Models")) +
  #guides(colour=guide_legend(title="Model")) + theme_bw()

  #PgEBH



  #PgEBS<- ggplot()  +
  # geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBS_inla$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="BS")) +
  #scale_color_manual(values = c("B" = "black",
  #                             "BH"="red",
  #                              "BS"="blue",
  #                              "BHS"="green")) + labs(x="", y="", title ="Genomic variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  #theme(legend.position = "left",
  #     legend.box.background = element_blank(),
  #      legend.title = element_text("Models")) +
  #guides(colour=guide_legend(title="Model")) + theme_bw()

  #PgEBS



#PgEBHS<- ggplot()  +
# geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
        #                                             fitBHS_inla$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="BHS")) +
  # scale_color_manual(values = c("B" = "black",
  #                              "BH"="red",
  #                              "BS"="blue",
  #                              "BHS"="green")) + labs(x="", y="", title ="Genomic variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  #theme(legend.position = "left",
  #     legend.box.background = element_blank(),
  #      legend.title = element_text("Models")) +
  # guides(colour=guide_legend(title="Model")) + theme_bw()

  #PgEBHS

# Combine the posterior distribution of genomic variance

##(Pg <- ggarrange(PgEB,PgEBH,PgEBS,PgEBHS,PgE, ncol = 3, nrow = 2, common.legend = F, legend = "left", align = "h", widths = c(1,1,1))) + PreseTheme



# Spatial range
Sr<- ggplot()  +

  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) (exp(x))^2,
                                                     fitGPS$internal.marginals.hyperpar[[4]]))), mapping = aes(x, y, colour="GPS")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) (exp(x))^2,
                                                     fitGPHS$internal.marginals.hyperpar[[5]]))), mapping = aes(x, y, colour="GPHS")) +
  scale_color_manual(values = c("GPS"="blue",
                                "GPHS"="green")) +
  labs(x="", y="", title ="Spatial range") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(),
        legend.title = element_text("Model")) +
  guides(colour=guide_legend(title="Model"))
Sr<- Sr + theme(plot.title = element_text(face = "bold")) # Making title Bold
Sr

#-------------------- GP Models_Final Variance components Paper and presentation-----------
(pE <- ggarrange(Pg,Pr,Ph,Pe,Sv,Sr, ncol = 2, nrow = 4, common.legend = T, legend = "left", align = "h", widths = c(1,1,1)))
pE

ggsave(plot = pE + PreseTheme, filename = "FigGPmodels_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = pE + PaperTheme, filename = "FigGPmodels_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper


#-------------------------Plot genomic variance for each model separately------

# Load required libraries
library(ggplot2)
library(INLA)

# Ensure that your fitted models (fitGP, fitGPH, fitGPS, fitGPHS) are loaded into the environment

# Step 1: Transform data
data_GP <- data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x, fitGP$marginals.hyperpar$`Precision for animal`)))
data_GP$Source <- "GP"

data_GPH <- data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x, fitGPH$marginals.hyperpar$`Precision for animal`)))
data_GPH$Source <- "GPH"

data_GPS <- data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x, fitGPS$marginals.hyperpar$`Precision for animal`)))
data_GPS$Source <- "GPS"

data_GPHS <- data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x, fitGPHS$marginals.hyperpar$`Precision for animal`)))
data_GPHS$Source <- "GPHS"

# Step 2: Plot each dataset individually

# Plot for data_GP
plot_GP <- ggplot(data_GP, aes(x = x, y = y)) +
  geom_line(colour = "black") +
  labs(x = "", y = "", title = "Genomic variance for GP") +
  theme_minimal() +
  theme(legend.position = "none")

# Plot for data_GPH
plot_GPH <- ggplot(data_GPH, aes(x = x, y = y)) +
  geom_line(colour = "red") +
  labs(x = "", y = "", title = "Genomic variance for GPH") +
  theme_minimal() +
  theme(legend.position = "none")

# Plot for data_GPS
plot_GPS <- ggplot(data_GPS, aes(x = x, y = y)) +
  geom_line(colour = "blue") +
  labs(x = "", y = "", title = "Genomic variance for GPS") +
  theme_minimal() +
  theme(legend.position = "none")

# Plot for data_GPHS
plot_GPHS <- ggplot(data_GPHS, aes(x = x, y = y)) +
  geom_line(colour = "green") +
  labs(x = "", y = "", title = "Genomic variance for GPHS") +
  theme_minimal() +
  theme(legend.position = "none")

# Display the plots
print(plot_GP)
print(plot_GPH)
print(plot_GPS)
print(plot_GPHS)

#-------------------------------------------------------------------------------







#-------------------------------------------------------------------------------
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
#genvarB <- plot (fitB, "animal")
#genvarB



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

(pE <- ggarrange(PgE,PrE,PhE,PwE,SvE,SrE, ncol = 2, nrow = 3, common.legend = T, legend = "left", align = "h", widths = c(1,1,1))) + PreseTheme
pE

ggsave(plot = pE + PreseTheme, filename = "Posterior_distribution_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = pE + PaperTheme, filename = "Posterior_distribution_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper

# ---- Plot posterior mean and Posterior standard deviation of the estimated spatial effects----

#fitGHS
#Adapted code from INLABru for Paper
# Generate prediction mesh (including uncertainty)

# Create a finer grid over the extent of data_border.
pred_GHS <- predict(
  fitGHS,
  fm_pixels(mesh, mask = data_border),
  ~ data.frame(
    mean_field = field,
    log_field = log(field)
  )
)

# Convert predictions into a usable data frame
plot_data <- pred_GHS$mean_field  # Extract the mean field

# Plot for Posterior Mean (Response Scale)
pl_mean <- ggplot(plot_data) +
  geom_tile(aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = mean)) +
  geom_sf(data = data_border, alpha = 0.1) +
  scale_fill_viridis_c(option = "D", direction = 1,name = NULL) +  # Higher values = lighter, lower = darker
  ggtitle("Posterior mean ") +
  xlab("Longitude (km)") + ylab("Latitude (km)") +
  theme_minimal() + theme(
    plot.title = element_text(face = "plain", size = 15,hjust = 0.5),   # Title
    axis.title = element_text(face = "plain", size = 15),   # Axis labels
    legend.text = element_text(face = "plain", size = 15),  # Legend text
    legend.title = element_text(face = "plain", size = 15),
    panel.grid = element_blank(),
    axis.text = element_blank()
  )
pl_mean
# Plot for Posterior Uncertainty (Standard Deviation)
pl_sd <- ggplot(plot_data) +
  geom_tile(aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = sd)) +
  geom_sf(data = data_border, alpha = 0.1) +
  scale_fill_viridis_c(option = "D", direction = 1,name = NULL) +  # Different color scale for variation
  ggtitle("Posterior standard deviation") +
  xlab("Longitude (km)") + ylab("Latitude (km)") +
  theme_minimal() + theme(
    plot.title = element_text(face = "plain", size = 15, hjust = 0.5),   # Title
    axis.title = element_text(face = "plain", size = 15),   # Axis labels
    legend.text = element_text(face = "plain", size = 15),
    panel.grid = element_blank(),# Legend text
    legend.title = element_text(face = "plain", size = 15),
      # Remove axis lines inside the plot area
    axis.text = element_blank()# Legend title
  )
pl_sd
# Arrange the plots in a single figure
plot_GHS<- grid.arrange(pl_mean, pl_sd, ncol = 2,
             top = NULL, # Optionally set top title
             padding = unit(c(0, 0, 0, 0), "cm"))


# Save the combined plot
PaperSize=10
ggsave("Figure4_GHS.png", plot = plot_GHS, height = PaperSize, width = PaperSize * 1.5, unit = "cm")

ggsave(plot = plot_GHS, filename = "Figure4_GHS.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper

#---------------------------------trying---------------------------------------
#fitGHS
# Load required packages
library(ggplot2)
library(sf)
library(gridExtra)
library(viridis)

# Convert predictions into a usable data frame
plot_data <- pred_GHS$mean_field  # Extract the mean field

# Plot for Posterior Mean (Response Scale)
pl_mean <- ggplot(plot_data) +
  geom_tile(aes(x = st_coordinates(geometry)[,1],
                y = st_coordinates(geometry)[,2],
                fill = mean)) +
  geom_sf(data = data_border, alpha = 0.1, size = 0.3) +  # Thinner border
  scale_fill_viridis_c(option = "D", direction = 1, name = NULL) +
  ggtitle("Posterior Mean") +
  xlab("Longitude (km)") + ylab("Latitude (km)") +
  theme_minimal(base_size = 18) +  # Increase base font size
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(face = "plain", size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank()  # Remove background grid
  )

# Plot for Posterior Uncertainty (Standard Deviation)
pl_sd <- ggplot(plot_data) +
  geom_tile(aes(x = st_coordinates(geometry)[,1],
                y = st_coordinates(geometry)[,2],
                fill = sd)) +
  geom_sf(data = data_border, alpha = 0.1, size = 0.3) +  # Thinner border
  scale_fill_viridis_c(option = "D", direction = 1, name = NULL) +
  ggtitle("Posterior Standard Deviation") +
  xlab("Longitude (km)") + ylab("Latitude (km)") +
  theme_minimal(base_size = 18) +  # Increase base font size
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(face = "plain", size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank()  # Remove background grid
  )

# Arrange the plots side by side
plot_GHS <- grid.arrange(pl_mean, pl_sd, ncol = 2)

# Save the figure with high resolution
ggsave("high_res_plot.png", plot_GHS, width = PaperSize * 1.5, height = PaperSize, dpi = 600)

#fitGHS------> Paper

gproj <- inla.mesh.projector(mesh2,  dims = c(300,300))

g.mean <- inla.mesh.project(gproj, fitGHS2$summary.random$field$mean)
g.sd <- inla.mesh.project(gproj, fitGHS2$summary.random$field$sd)
summary(fitGHS2)
str(g.sd)
hist(c(g.sd))
g.mean[g.sd > 5] = NA
g.sd[g.sd > 5] = NA



mean_plot <- levelplot(g.mean, scales=list(draw=F),
                       xlab='Longitude (km)', ylab='Latitude (km)',
                       main=list(label='Posterior mean', fontface = "plain"),
                       col.regions = viridis(16))

sd_plot <- levelplot(g.sd, scales=list(draw=F),
                     xlab='Longitude (km)', ylab='Latitude (km)',
                     main=list(label='Posterior standard deviation', fontface = "plain"),
                     col.regions = viridis(16))

# Combine the plots with grid.arrange and reduce margins
plot_GHS <- grid.arrange(mean_plot, sd_plot, ncol=2,
                         top = NULL, # Optionally set top title
                         padding = unit(c(0, 0, 0, 0), "cm")) # Set padding to zero


# Save the combined plot
ggsave("Figure4_GHS.png", plot = plot_GHS, height = PaperSize, width = PaperSize * 1.5, unit = "cm")


#fitGS2------> Paper

gproj <- inla.mesh.projector(mesh2,  dims = c(300,300))

g.mean <- inla.mesh.project(gproj, fitGS2$summary.random$field$mean)
g.sd <- inla.mesh.project(gproj, fitGS2$summary.random$field$sd)
summary(fitGS2)
str(g.sd)
hist(c(g.sd))
g.mean[g.sd > 5] = NA
g.sd[g.sd > 5] = NA



mean_plot <- levelplot(g.mean, scales=list(draw=F),
                       xlab='Longitude (km)', ylab='Latitude (km)',
                       main=list(label='Posterior mean', fontface = "plain"),
                       col.regions = viridis(16))

sd_plot <- levelplot(g.sd, scales=list(draw=F),
                     xlab='Longitude (km)', ylab='Latitude (km)',
                     main=list(label='Posterior standard deviation', fontface = "plain"),
                     col.regions = viridis(16))

# Combine the plots with grid.arrange and reduce margins
plot_GS <- grid.arrange(mean_plot, sd_plot, ncol=2,
                         top = NULL, # Optionally set top title
                         padding = unit(c(0, 0, 0, 0), "cm")) # Set padding to zero


# Save the combined plot
ggsave("Supplementary_Figure_GS.png", plot = plot_GS, height = PaperSize, width = PaperSize * 1.5, unit = "cm")



#----------------------------------#-----------------------#--------------------
#GS

pred_GS <- predict(
  fitGS,
  fm_pixels(mesh, mask = data_border),
  ~ data.frame(
    mean_field = field,
    log_field = log(field)
  )
)

# Convert predictions into a usable data frame
plot_data <- pred_GS$mean_field  # Extract the mean field

# Plot for Posterior Mean (Response Scale)
pl_mean <- ggplot(plot_data) +
  geom_tile(aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = mean)) +
  geom_sf(data = data_border, alpha = 0.1) +
  scale_fill_viridis_c(option = "D", direction = 1,name = NULL) +  # Higher values = lighter, lower = darker
  ggtitle("Posterior mean ") +
  xlab("Longitude (km)") + ylab("Latitude (km)") +
  theme_minimal() + theme(
    plot.title = element_text(face = "plain", size = 15,hjust = 0.5),   # Title
    axis.title = element_text(face = "plain", size = 15),   # Axis labels
    legend.text = element_text(face = "plain", size = 15),  # Legend text
    legend.title = element_text(face = "plain", size = 15),
    panel.grid = element_blank(),
    axis.text = element_blank()
  )
pl_mean
# Plot for Posterior Uncertainty (Standard Deviation)
pl_sd <- ggplot(plot_data) +
  geom_tile(aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = sd)) +
  geom_sf(data = data_border, alpha = 0.1) +
  scale_fill_viridis_c(option = "D", direction = 1,name = NULL) +  # Different color scale for variation
  ggtitle("Posterior standard deviation") +
  xlab("Longitude (km)") + ylab("Latitude (km)") +
  theme_minimal() + theme(
    plot.title = element_text(face = "plain", size = 15, hjust = 0.5),   # Title
    axis.title = element_text(face = "plain", size = 15),   # Axis labels
    legend.text = element_text(face = "plain", size = 15),
    panel.grid = element_blank(),# Legend text
    legend.title = element_text(face = "plain", size = 15),
    # Remove axis lines inside the plot area
    axis.text = element_blank()# Legend title
  )

# Arrange the plots in a single figure
plot_GS<- grid.arrange(pl_mean, pl_sd, ncol = 2,
                        top = NULL, # Optionally set top title
                        padding = unit(c(0, 0, 0, 0), "cm"))

# Save the combined plot
PaperSize=10
ggsave("Supplementary_Figure_GS.png", plot = plot_GS, height = PaperSize, width = PaperSize * 1.5, unit = "cm")




ggsave("FigureGHS_papertheme.png", plot = plot_GHS, height = PaperSize, width = PaperSize * 1.5, unit = "cm")
#--------------------------Number with herds with single cows (997/1386=72%)----
df <- data2[,c("herd","cow")]
df <- distinct(df)
single_cow_herds <- df %>%
  group_by(herd) %>%
  summarise(cow_count = n())

summary(single_cow_herds$cow_count)

table(single_cow_herds$cow_count)

prop.table(table(single_cow_herds$cow_count))

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
save(fitGP_4_NA,file = "data/cleaned_data/cross_validation_breed_proportion/fitGP/fitGP_4_NA.RData")


# Predicting with-held phenotypes
predictGP_1 <- predict(fitGP_1_NA, newdata=data2_1,
                      formula= ~ fixed_effects + perm  + animal)

predictGP_2 <- predict(fitGP_2_NA, newdata=data2_2,
                      formula= ~ fixed_effects + perm  + animal)

predictGP_3 <- predict(fitGP_3_NA, newdata=data2_3,
                      formula= ~ fixed_effects + perm  + animal)

predictGP_4 <- predict(fitGP_4_NA, newdata=data2_4,
                      formula= ~ fixed_effects + perm  + animal)

# Accuracy of fitGP

Accu1 <- round(cor(predictGP_1$milkZ,predictGP_1$mean),2)
Accu1 #0.29

Accu2 <- round(cor(predictGP_2$milkZ,predictGP_2$mean),2)
Accu2 #0.35

Accu3 <- round(cor(predictGP_3$milkZ,predictGP_3$mean),2)
Accu3 # 0.39

Accu4 <- round(cor(predictGP_4$milkZ,predictGP_4$mean),2)
Accu4 # 0.4
# Cross-validation accuracy of fitGP
Accu_GP <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GP # 0.37

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
Accu1 # 0.56 (old 0.55)

Accu2 <- round(cor(predictGPS_2$milkZ,predictGPS_2$mean),2)
Accu2 # 0.65 (old0.63)

Accu3 <- round(cor(predictGPS_3$milkZ,predictGPS_3$mean),2)
Accu3 # 0.64 (old 0.63)

Accu4 <- round(cor(predictGPS_4$milkZ,predictGPS_4$mean),2)
Accu4 # 0.57 (Old 0.56)

# Cross-validation accuracy of fitGPS
Accu_GPS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GPS #  0.60 (old 0.59)


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
Accu1 #0.58 (0.57)

Accu2 <- round(cor(predictGPHS_2$milkZ,predictGPHS_2$mean),2)
Accu2 # 0.67 (0.66)

Accu3 <- round(cor(predictGPHS_3$milkZ,predictGPHS_3$mean),2)
Accu3 # 0.66 (old 0.65)

Accu4 <- round(cor(predictGPHS_4$milkZ,predictGPHS_4$mean),2)
Accu4 # 0.57(0.54)
# Cross-validation accuracy of fitGPH
Accu_GPHS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GPHS # 0.62  (0.60)

#-----------------Cross_validation_Region_Paper-----------------------------------------
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
Accu1 # 0.10 (0.02)

Accu2 <- round(cor(predictGPS_2$milkZ,predictGPS_2$mean),2)
Accu2 # -0.08 (0.04)

Accu3 <- round(cor(predictGPS_3$milkZ,predictGPS_3$mean),2)
Accu3 # 0.27 (0.11)

Accu4 <- round(cor(predictGPS_4$milkZ,predictGPS_4$mean),2)
Accu4 # 0.17 (0.19)
# Cross-validation accuracy of fitGPS
Accu_GPS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GPS # 0.12 (0.09)

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
Accu1 # 0.04 (0.05)

Accu2 <- round(cor(predictGPHS_2$milkZ,predictGPHS_2$mean),2)
Accu2 #-0.03 (0.03)

Accu3 <- round(cor(predictGPHS_3$milkZ,predictGPHS_3$mean),2)
Accu3 # 0.16 (0.07)

Accu4 <- round(cor(predictGPHS_4$milkZ,predictGPHS_4$mean),2)
Accu4 # 0.13 (0.54)
# Cross-validation accuracy of fitGPHS
Accu_GPHS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GPHS # 0.08 (0.17)


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
AccuGP #

# Accuracy forward validation GPH
AccuGPH <- round(cor(predictGPH$milkZ,predictGPH$mean),2)
AccuGPH #

# Accuracy forward validation GPS
AccuGPS <- round(cor(predictGPS$milkZ,predictGPS$mean),2)
AccuGPS# 0.62

# Accuracy forward validation GPHS
AccuGPHS <- round(cor(predictGPHS$milkZ,predictGPHS$mean),2)
AccuGPHS #0.67

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
Accu_dgrp1 #0.4

Accu_dgrp2 <-  round(cor(predictGPS_dgrp2$milkZ,predictGPS_dgrp2$mean),2)
Accu_dgrp2 #0.64

Accu_dgrp3 <-  round(cor(predictGPS_dgrp3$milkZ,predictGPS_dgrp3$mean),2)
Accu_dgrp3 #0.74

Accu_dgrp4 <-  round(cor(predictGPS_dgrp4$milkZ,predictGPS_dgrp4$mean),2)
Accu_dgrp4 #0.72

Accu_dgrp <- round((Accu_dgrp1 + Accu_dgrp2 + Accu_dgrp3 + Accu_dgrp4)/4,2)
Accu_dgrp #0.62

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
Accu4 # 0.4
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
Accu_GH #0.4


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
Accu1 #  0.53  (old 0.54)

Accu2 <- round(cor(predictGS_2$milkZ,predictGS_2$mean),2)
Accu2 # 0.62 (old 0.60)

Accu3 <- round(cor(predictGS_3$milkZ,predictGS_3$mean),2)
Accu3 # 0.60 (Old 0.58)

Accu4 <- round(cor(predictGS_4$milkZ,predictGS_4$mean),2)
Accu4 # 0.48 (Old0.50)

# Cross-validation accuracy of fitGS
Accu_GS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GS # 0.56 (old 0.57)
#rm(fitPS_4_NA,fitGPS_3_NA,fitGPS_2_NA,fitGPS_1_NA)


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
Accu1 # 0.57 (old0.58)

Accu2 <- round(cor(predictGHS_2$milkZ,predictGHS_2$mean),2)
Accu2 # 0.65 (old0.66)

Accu3 <- round(cor(predictGHS_3$milkZ,predictGHS_3$mean),2)
Accu3 # 0.64 (old0.65)

Accu4 <- round(cor(predictGHS_4$milkZ,predictGHS_4$mean),2)
Accu4 # 0.53 (old0.53)
# Cross-validation accuracy of fitGH
Accu_GHS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GHS #  0.60 (old0.60)

#-----------------Cross_validation_Region_Paper-----------------------------------------
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
Accu1 #-0.06

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
Accu1 # 0.1

Accu2 <- round(cor(predictGH_2$milkZ,predictGH_2$mean),2)
Accu2 #0.1

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
Accu1 # 0.02 (old 0.05)

Accu2 <- round(cor(predictGS_2$milkZ,predictGS_2$mean),2)
Accu2 # 0.011 ( old 0.0)

Accu3 <- round(cor(predictGS_3$milkZ,predictGS_3$mean),2)
Accu3 # 0.24 (old 0.24)

Accu4 <- round(cor(predictGS_4$milkZ,predictGS_4$mean),2)
Accu4 # 0.15(old 0.10)
# Cross-validation region accuracy of fitGS
Accu_GS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GS # 0.13 (old0.2)

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
Accu1 # 0.11 (old 0.06)

Accu2 <- round(cor(predictGHS_2$milkZ,predictGHS_2$mean),2)
Accu2 # -0.02 (Old -0.01)

Accu3 <- round(cor(predictGHS_3$milkZ,predictGHS_3$mean),2)
Accu3 # 0.18 (old 0.19)

Accu4 <- round(cor(predictGHS_4$milkZ,predictGHS_4$mean),2)
Accu4 # 0.14 (Old 0.13)
# Cross-validation accuracy of fitGHS
Accu_GHS <- round((Accu1 + Accu2 + Accu3 + Accu4)/4,2)
Accu_GHS # 0.10 (Old 0.09)


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
AccuG #0.44

# Accuracy forward validation GH
AccuGH <- round(cor(predictGH$milkZ,predictGH$mean),2)
AccuGH #0.59

# Accuracy forward validation GS
AccuGS <- round(cor(predictGS$milkZ,predictGS$mean),2)
AccuGS# 0.59 (old0.57)

# Accuracy forward validation GHS
AccuGHS <- round(cor(predictGHS$milkZ,predictGHS$mean),2)
AccuGHS #0.65 (old 0.63)


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
Accu_dgrp # 0.48


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
Accu_dgrp1 # 0.37 (0.34)

Accu_dgrp2 <-  round(cor(predictGS_dgrp2$milkZ,predictGS_dgrp2$mean),2)
Accu_dgrp2 # 0.59 (old 0.57)

Accu_dgrp3 <-  round(cor(predictGS_dgrp3$milkZ,predictGS_dgrp3$mean),2)
Accu_dgrp3 # 0.73 (old 0.71)

Accu_dgrp4 <-  round(cor(predictGS_dgrp4$milkZ,predictGS_dgrp4$mean),2)
Accu_dgrp4 # 0.74 (0.77)

Accu_dgrp <- round((Accu_dgrp1 + Accu_dgrp2 + Accu_dgrp3 + Accu_dgrp4)/4,2)
Accu_dgrp # 0.61 (old 0.6).

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
Accu_dgrp1 # 0.44 (old 0.43)

Accu_dgrp2 <-  round(cor(predictGHS_dgrp2$milkZ,predictGHS_dgrp2$mean),2)
Accu_dgrp2 # 0.67 (old 0.65)

Accu_dgrp3 <-  round(cor(predictGHS_dgrp3$milkZ,predictGHS_dgrp3$mean),2)
Accu_dgrp3 # 0.75 (old 0.73)

Accu_dgrp4 <-  round(cor(predictGHS_dgrp4$milkZ,predictGHS_dgrp4$mean),2)
Accu_dgrp4 # 0.66 (old 0.60)

Accu_dgrp <- round((Accu_dgrp1 + Accu_dgrp2 + Accu_dgrp3 + Accu_dgrp4)/4,2)
Accu_dgrp # 0.63 (old 0.60)


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
Accu_region1 # 0.57 (0.56)

Accu_region2 <-  round(cor(predictGS_region2$milkZ,predictGS_region2$mean),2)
Accu_region2 # 0.53 (same)

Accu_region3 <-  round(cor(predictGS_region3$milkZ,predictGS_region3$mean),2)
Accu_region3 # 0.28 (old 0.21)

Accu_region4 <-  round(cor(predictGS_region4$milkZ,predictGS_region4$mean),2)
Accu_region4 # 0.37 (old 0.34)

Accu_region <- round((Accu_region1 + Accu_region2 + Accu_region3 + Accu_region4)/4,2)
Accu_region # 0.44 (old 0.41)



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
Accu_region1 #  0.65 (0.64)

Accu_region2 <-  round(cor(predictGHS_region2$milkZ,predictGHS_region2$mean),2)
Accu_region2 # 0.60 (old 0.61)

Accu_region3 <-  round(cor(predictGHS_region3$milkZ,predictGHS_region3$mean),2)
Accu_region3 # 0.37 (old 0.33)

Accu_region4 <-  round(cor(predictGHS_region4$milkZ,predictGHS_region4$mean),2)
Accu_region4 # 0.46(old 0.42)

Accu_region <- round((Accu_region1 + Accu_region2 + Accu_region3 + Accu_region4)/4,2)
Accu_region #  0.52 (old 0.50)

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
#load(file = "data/cleaned_data/geno.RData")
#install.packages(pkg = "irlba")
library(package = "irlba")
# stats::prcomp(x = geno[1:10, c(2:20)])
pcasAll <- prcomp_irlba(x = geno[, -1], n = 3)
summary(pcasAll)
#summary(pcas)
par(mfrow = c(2, 2))
plot(pcasAll$x[, 1], pcasAll$x[, 2], pch = 19, cex = 0.1)
plot(pcasAll$x[, 1], pcasAll$x[, 3], pch = 19, cex = 0.1)
plot(pcasAll$x[, 2], pcasAll$x[, 3], pch = 19, cex = 0.1)
par(mfrow = c(1, 1))
plot(pcasAll$x[, 1], pcasAll$x[, 2], pch = 19, cex = 0.5)

pcas <- as.data.frame(cbind(geno[, 1], pcasAll$x))
colnames(pcas)[1] <- "cow"
head(pcas)
tmp <- data1[, c("cow", "dgrp", "region")]
tmp <- tmp[!duplicated(tmp), ]
dim(tmp)
head(tmp)
head(pcas)
pcas <- merge(x =pcas , y =tmp , by = "cow", all.x = TRUE)
dim(pcas)
str(pcas)


# Visualize PCA by breed Proportion dgrp

#Proportion of exotic genes, as follows:
# 4 =[0, 36)\%,
# 3=[36, 60)\%,
# 2=[60, 87.5)\%, and
#1=[87.5, 100)\%.

# PC1 vs PC2
dgrp_pca1<- ggplot(data = pcas) +
  geom_point(mapping = aes(x = `PC1`, y = `PC2`,color = dgrp, size=I(0.2), alpha=I(0.5)), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x = paste0("PC1(3.1%)"),
       y = paste0("PC2(0.8%)"), color="Proportion of exotic genes(%)") +
  scale_color_manual(values=c("1"="green", "2"="red", "3"="orange","4"="blue"),
                     labels=c("1"="[100, 87.5]", "2"="(87.5, 60]", "3"="(60, 36]","4"="(36, 0]"),
                     name= "Proportion of exotic genes(%)") +
  theme_minimal()

dgrp_pca1

# PC1 vs PC3

dgrp_pca2<- ggplot(data = pcas) +
  geom_point(mapping = aes(x = `PC1`, y = `PC3`,color = dgrp, size=I(0.2), alpha=I(0.5)), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x = paste0("PC1(3.1%)"),
       y = paste0("PC3(0.5%)"), color="Proportion of exotic genes(%)") +
  scale_color_manual(values=c("1"="green", "2"="red", "3"="orange","4"="blue"),
                     labels=c("1"="[87.5, 100)", "2"="[60, 87.5)", "3"="[36, 60)","4"="[0, 36)"),
                     name= "Proportion of exotic genes(%)") +
  theme_minimal()

dgrp_pca2


# PC2 vs PC3
dgrp_pca3<- ggplot(data = pcas) +
  geom_point(mapping = aes(x = `PC2`, y = `PC3`,color = dgrp, size=I(0.2), alpha=I(0.5)), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x = paste0("PC2(0.8%)"),
       y = paste0("PC3(0.5%)"), color="Proportion of exotic genes(%)") +
  scale_color_manual(values=c("1"="green", "2"="red", "3"="orange","4"="blue"),
                     labels=c("1"="[87.5, 100)", "2"="[60, 87.5)", "3"="[36, 60)","4"="[0, 36)"),
                     name= "Proportion of exotic genes(%)") +
  theme_minimal()

dgrp_pca3


dgrp_pca <- ggarrange(dgrp_pca1,dgrp_pca2,dgrp_pca3, ncol = 3, nrow = 1, common.legend = T, legend = "right", align = "h", widths = c(1,1,1))
dgrp_pca

ggsave(plot = dgrp_pca + PreseTheme, filename = "aloha0.1PCA_breed_proportion_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = dgrp_pca + PaperTheme, filename = "PCA_breed_proportion_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper


#Visualize PCA by region
#North-Central (latitude $>$ -4),
#North-East (latitude $<$ -4 \& longitude  $>$ 36),
#South-Central (latitude $<$ -7 \& longitude $>$ 34), and
#South-West (latitude $<$ -7 \& longitude $<$ 34).

#North-Central (latitude $>$ -4): ----> Region 2 (light)

#North-East (latitude $<$ -4 & longitude > 36)--->  Region 1 (dark)

#South-Central (latitude < -7 & longitude > 34)-Region 3 (dark)

#South-West (latitude < -7 & longitude < 34)----> Region 4 light

# PC1 vs PC2
region_pca1<- ggplot(data = pcas) +
  geom_point(mapping = aes(x = `PC1`, y = `PC2`,color = region, size=I(0.2), alpha=I(0.5)), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x = paste0("PC1(3.1%)"),
       y = paste0("PC2(0.8%)"), color="Regions") +
  scale_color_manual(values=c("1"="blue", "2"="green", "3"="black","4"="orange"),
                     labels=c("1"="North-East", "2"="North-Central", "3"="South-Central","4"="South-West"),
                     name= "Regions") +
  theme_minimal()

region_pca1



# region_pca1 (Polygone with transparency)
#ggplot(pcas, aes(x = PC1, y = PC2, color = region)) +
 # geom_point(size=0.2) +
 # stat_ellipse(geom = "polygon",
   #            aes(fill = region),
    #           alpha = 0.1)

# PC1 vs PC3
region_pca2<- ggplot(data = pcas) +
  geom_point(mapping = aes(x = `PC1`, y = `PC3`,color = region, size=I(0.2), alpha=I(0.5)), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x = paste0("PC1(3.1%)"),
       y = paste0("PC3(0.5%)"), color="Regions") +
  scale_color_manual(values=c("1"="blue", "2"="green", "3"="black","4"="orange"),
                     labels=c("1"="North-East", "2"="North-Central", "3"="South-Central","4"="South-West"),
                     name= "Regions") +
  theme_minimal()

region_pca2


# PC2 vs PC3

region_pca3<- ggplot(data = pcas) +
  geom_point(mapping = aes(x = `PC2`, y = `PC3`,color = region, size=I(0.2), alpha=I(0.5)), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x = paste0("PC2(0.8%)"),
       y = paste0("PC3(0.5%)"), color="Regions") +
  scale_color_manual(values=c("1"="blue", "2"="green", "3"="black","4"="orange"),
                     labels=c("1"="North-East", "2"="North-Central", "3"="South-Central","4"="South-West"),
                     name= "Regions") +
  theme_minimal()

region_pca3

# region_pca2 (Polygone with transparency)
#ggplot(pcas, aes(x = PC1, y = PC3, color = region)) +
#  geom_point(size=0.2) +
#  stat_ellipse(geom = "polygon",
         #      aes(fill = region),
             #  alpha = 0.1)


# region_pca3 (Polygone with transparency)
#ggplot(pcas, aes(x = PC2, y = PC3, color = region)) +
#  geom_point(size=0.2) +
 # stat_ellipse(geom = "polygon",
 #              aes(fill = region),
  #             alpha = 0.1)

region_pca <- ggarrange(region_pca1,region_pca2,region_pca3, ncol = 3, nrow = 1, common.legend = T, legend = "right", align = "h", widths = c(1,1,1))
region_pca

ggsave(plot = region_pca + PreseTheme, filename = "alpha0.1PCA_region_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = region_pca + PaperTheme, filename = "PCA_region_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper

#-----Final Plot paper----------------------------------------------------------

pca <-  ggarrange(dgrp_pca, region_pca, ncol = 1, nrow = 2, common.legend = T, legend = "right", align = "h", widths = c(1,1,1))
pca
ggsave(plot = pca + PreseTheme, filename = "pca_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot =  pca + PaperTheme, filename = "pca_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper



table(data1$dgrp,data1$region)
#Breed proportion in First column
#    1    2    3    4
#1  497 4663 1076 1230
#2 2248 3012 1704 1185
#3 1495  704  622  237
#4  307  103  272   20

data3 <- data1[, c("cow", "dgrp", "region")]
data3<- distinct(data3)

table(data3$dgrp,data3$region)

#   NE  NC  SC  SW
#1  49 417 134 131
#2 198 275 183 114
#3 154  71  64  20
#4  36  12  33   3

#-----------------------Manuscript- Correlations and Animal Ranking------------------------------------------
#fitG
EBV_G <- fitG$summary.random$animal[,1:3]
names(EBV_G)[2] <- "EBV_G"
names(EBV_G)[3] <- "sd_G"


# fitGH
EBV_GH <- fitGH$summary.random$animal[,1:3]
names(EBV_GH)[2] <- "EBV_GH"
names(EBV_GH)[3] <- "sd_GH"

# fitGS

EBV_GS <- fitGS$summary.random$animal[,1:3]
names(EBV_GS)[2] <- "EBV_GS"
names(EBV_GS)[3] <- "sd_GS"


# fitGHS

EBV_GHS <- fitGHS$summary.random$animal[,1:3]
names(EBV_GHS)[2] <- "EBV_GHS"
names(EBV_GHS)[3] <- "sd_GHS"

#----------Combined EBVs manuscripts models-------------------------------------
EBV_sd<- merge(EBV_G,EBV_GH, by="ID")
EBV_sd<- merge(EBV_sd, EBV_GS,by="ID")
EBV_sd<- merge(EBV_sd, EBV_GHS, by="ID")

EBV <- EBV_sd[,c("ID","EBV_G","EBV_GH","EBV_GS","EBV_GHS")]

#--------------------Predict spatial effect at my locations (Models GS and GHS)-
# Spatial effect in model GS
fieldGS <- predict(fitGS,spdf, ~field)
#fieldGS_df <- data.frame(fieldGS$[,c("cow" , "herd","mean", "sd")])
fieldGS_df <- data.frame(fieldGS[,c("mean","sd")])
fieldGS_df$cow <- data2$cow
fieldGS_df$herd <- data2$herd
fieldGS_df <- distinct(fieldGS_df)
length(levels(fieldGS_df$herd)) # 1386 herds
# Spatial effect GS
Spatial_GS <- fieldGS_df[,c("cow","herd", "mean")]
names(Spatial_GS)[1]<- "ID"
names(Spatial_GS)[3]<- "Spatial_effect_GS"

# Spatial effect in model GHS
fieldGHS <- predict(fitGHS,spdf, ~field)
fieldGHS_df <- data.frame(fieldGHS[,c("mean", "sd")])
fieldGHS_df$cow <- data2$cow
fieldGHS_df$herd <- data2$herd
fieldGHS_df <- distinct(fieldGHS_df)

# Spatial effect GHS
Spatial_GHS <- fieldGHS_df[,c("cow", "mean")]
names(Spatial_GHS)[1]<- "ID"
names(Spatial_GHS)[2]<- "Spatial_effect_GHS"

# Herd effect in model GHS
spdf_herd <- spdf
spdf_herd$herd <- data2$herd
herdGHS <- predict(fitGHS,spdf_herd, ~herd)

herdGHS_df <- data.frame(herdGHS[,c("mean", "sd")])
herdGHS_df$cow <- data2$cow
herdGHS_df <- distinct(herdGHS_df)

# Herd effect GHS
Herd_GHS <- herdGHS_df[,c("cow", "mean")]
names(Herd_GHS)[1]<- "ID"
names(Herd_GHS)[2]<- "Herd_effect_GHS"

# Herd effect in model GH
herdGH <- predict(fitGH,spdf_herd, ~herd)

herdGH_df <- data.frame(herdGH[,c("mean", "sd")])
herdGH_df$cow <- data2$cow
herdGH_df <- distinct(herdGH_df)

# Herd effect GH
Herd_GH <- herdGH_df[,c("cow", "mean")]
names(Herd_GH)[1]<- "ID"
names(Herd_GH)[2]<- "Herd_effect_GH"


Herd_effect <- merge(Herd_GH, Herd_GHS, by="ID", all.x=TRUE )

# ------------EBVs, Spatial effect and herd effects ----------------------------
EBV_herd_effect <- merge(EBV, Herd_effect, by="ID", all.x=TRUE )
EBVs_Spatial_effect <- merge(EBV_herd_effect, Spatial_GHS, by="ID", all.x=TRUE )
#plot(EBVs_Spatial_effect$Spatial_effect_GHS,EBVs_Spatial_effect$EBV_G)
EBVs_Spatial_Matrix <- EBVs_Spatial_effect[, -c(1)]

round(cor(EBVs_Spatial_effect$EBV_GH,EBVs_Spatial_effect$EBV_GHS),2) #0.87 vs (0.89 old)

#create matrix of correlation coefficients and p-values for EBVs between models
# How to Create a Correlation Matrix in R (4 Examples)-Statology
# https://www.statology.org/correlation-matrix-in-r/

str(EBVs_Spatial_Matrix)
  # Pearson Correlation
round(cor(EBVs_Spatial_Matrix, method = "pearson"),2)
#                 EBV_G EBV_GH EBV_GS EBV_GHS Herd_effect_GH Herd_effect_GHS Spatial_effect_GHS
#EBV_G               1.00   0.78   0.76    0.66           0.92            0.66               0.64
#EBV_GH              0.78   1.00   0.68    0.87           0.52            0.33               0.43
#EBV_GS              0.76   0.68   1.00    0.87           0.65            0.87               0.06
#EBV_GHS             0.66   0.87   0.87    1.00           0.42            0.54               0.06
#Herd_effect_GH      0.92   0.52   0.65    0.42           1.00            0.72               0.71
#Herd_effect_GHS     0.66   0.33   0.87    0.54           0.72            1.00               0.06
#Spatial_effect_GHS  0.64   0.43   0.06    0.06           0.71            0.06               1.00

# Spearman's Rank Correlations
round(cor(EBVs_Spatial_Matrix, method = "spearman"),2)
#                   EBV_G EBV_GH EBV_GS EBV_GHS Herd_effect_GH Herd_effect_GHS Spatial_effect_GHS
#EBV_G               1.00   0.77   0.74    0.63           0.92            0.64               0.65
#EBV_GH              0.77   1.00   0.65    0.84           0.53            0.32               0.44
#EBV_GS              0.74   0.65   1.00    0.86           0.62            0.85               0.06
#EBV_GHS             0.63   0.84   0.86    1.00           0.40            0.53               0.06
#Herd_effect_GH      0.92   0.53   0.62    0.40           1.00            0.69               0.71
#Herd_effect_GHS     0.64   0.32   0.85    0.53           0.69            1.00               0.06
#Spatial_effect_GHS  0.65   0.44   0.06    0.06           0.71            0.06               1.00

#Adding difference EBV_GH - EBV_GHS

EBVs_Spatial_effect <- EBVs_Spatial_effect %>% mutate(dEBV= EBV_GH - EBV_GHS)
EBVs_Spatial_effect
round(cor(EBVs_Spatial_effect$dEBV,EBVs_Spatial_effect$Spatial_effect_GHS),2)
#0.75

#--------MANOVA----------------------------------------------------------------
# Multivariate analysis of variances: Genotype Principal components (PC) by breed proportion (dgrp, with 4 levels) and by region (with 4 levels)
#Pcas is a dataframe with Principal components and factors (dgrp and region)
with(pcas, plot(PC1 ~ PC2, col = c("blue", "red", "black", "yellow")[dgrp]))
fit <- manova(cbind(PC1, PC2, PC3) ~ dgrp, data = pcas)
summary(fit) # p= 0.0069 ** --> 0.0069 **!
summary.aov(fit)
# PC1: p=0.0001885 ***
# PC2: p=0.6658
# PC3: p=0.7028

#PC1 ~ dgrp
fit <- lm(PC1 ~ dgrp, data = pcas)
summary(fit)
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)   10.961      3.225   3.399 0.000691 ***
#  dgrp2        -15.191      4.495  -3.379 0.000741 ***
#  dgrp3        -19.947      5.902  -3.379 0.000741 ***
#  dgrp4        -27.966     10.015  -2.792 0.005283 **

## PC2 ~ dgrp
fit <- lm(PC2 ~ dgrp, data = pcas)
summary(fit)
# --> nothing is really significant
#Coefficients:
 # Estimate Std. Error t value Pr(>|t|)
#(Intercept)  -0.7973     1.7024  -0.468    0.640
#dgrp2         1.9755     2.3729   0.832    0.405
#dgrp3         0.7710     3.1157   0.247    0.805
#dgrp4        -3.8048     5.2865  -0.720    0.472

## PC3 ~ dgrp
fit <- lm(PC3 ~ dgrp, data = pcas)
summary(fit)
# --> nothing is really significant
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.8370     1.2728   0.658    0.511
#dgrp2        -1.9389     1.7741  -1.093    0.275
#dgrp3        -0.0158     2.3294  -0.007    0.995
#dgrp4        -1.1569     3.9524  -0.293    0.770

#MANOVA Region
with(pcas, plot(PC1 ~ PC2, col = c("blue", "red", "black", "yellow")[region]))
fit <- manova(cbind(PC1, PC2, PC3) ~ region, data = pcas)
summary(fit) # p=5.991e-08 *** --> highly significant
summary.aov(fit)
# PC1: p=0.1287
# PC2: p=4.448e-07 ***
# PC3: p=0.003763 **

#PC1 ~ region
fit <- lm(PC1 ~ region, data = pcas)
summary(fit)
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)   -7.820      4.172  -1.874   0.0611 .
#region2       12.402      5.224   2.374   0.0177 *
#region3        7.920      5.982   1.324   0.1857
#region4        9.104      6.767   1.345   0.1787

plot(y = data1$lat, x = data1$long, col = c("blue", "red", "black", "yellow")[data1$region])
plot(y = data2$lat, x = data2$long)
nrow(pcas)
#PC2 ~ region
fit <- lm(PC2 ~ region, data = pcas)
summary(fit)
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)    1.336      2.176   0.614 0.539365
#region2       -6.299      2.725  -2.312 0.020899 *
# region3       -2.283      3.120  -0.732 0.464556
#region4       11.960      3.530   3.388 0.000718 ***

#PC3 ~ region
fit <- lm(PC3 ~ region, data = pcas)
summary(fit)
#Coefficients:
 # Estimate Std. Error t value Pr(>|t|)
#(Intercept)  -0.8764     1.6353  -0.536   0.5921
#region2       2.6398     2.0474   1.289   0.1974
#region3       2.7453     2.3446   1.171   0.2418
#region4      -5.6999     2.6523  -2.149   0.0318 *

  tmp <- data1[!duplicated(data1[, c("cow")]), ]
with(tmp, table(dgrp, region))
summary(lm(as.numeric(dgrp) ~ region, data = tmp))

par(mfrow = c(2, 2))
lims <- range(pcas$PC1)
tmp <- pcas[pcas$dgrp == "1", ]
hist(tmp$PC1, xlim = lims, main = "1")
tmp <- pcas[pcas$dgrp == "2", ]
hist(tmp$PC1, xlim = lims, main = "2")
tmp <- pcas[pcas$dgrp == "3", ]
hist(tmp$PC1, xlim = lims, main = "3")
tmp <- pcas[pcas$dgrp == "4", ]
hist(tmp$PC1, xlim = lims, main = "4")

# Interaction dgrp x region

fit <- manova(cbind(PC1, PC2, PC3) ~ dgrp * region, data = pcas)
summary(fit) # p=5.991e-08 *** --> highly significant
#Df   Pillai approx F num Df den Df    Pr(>F)
#dgrp           3 0.012180   2.5451      9   5619 0.0064856 **
#region         3 0.024588   5.1593      9   5619 5.361e-07 ***
#dgrp:region    9 0.032912   2.3084     27   5619 0.0001381 ***

summary.aov(fit)
# PC1: dgrp:region p=0.0489639 *
# PC2: dgrp:region p=0.01582 *
# PC3: p= 0.003458 **

#PC1 ~ dgrp*region
fit <- lm(PC1 ~ dgrp * region, data = pcas)
summary(fit)
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)     17.426     12.395   1.406   0.1599
#dgrp2          -28.037     13.844  -2.025   0.0430 *
#  dgrp3          -34.314     14.231  -2.411   0.0160 *
#  dgrp4           -5.463     19.046  -0.287   0.7743
#region2         -6.328     13.111  -0.483   0.6294
#region3         -2.556     14.485  -0.176   0.8599
#region4        -13.312     14.529  -0.916   0.3597
#dgrp2:region2   17.051     15.404   1.107   0.2685
#dgrp3:region2   16.321     18.078   0.903   0.3667
#dgrp4:region2  -54.432     31.754  -1.714   0.0867 .
#dgrp2:region3   11.249     16.999   0.662   0.5082
#dgrp3:region3   16.037     19.399   0.827   0.4085
#dgrp4:region3  -51.288     25.437  -2.016   0.0439 *
  #dgrp2:region4   16.595     17.752   0.935   0.3500
#dgrp3:region4   56.797     25.226   2.252   0.0245 *
 # dgrp4:region4   37.537     54.125   0.694   0.4881

#PC2 ~ dgrp*region
fit <- lm(PC2 ~ dgrp * region, data = pcas)
summary(fit)
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)    -12.913      6.483  -1.992 0.046535 *
#  dgrp2           17.581      7.241   2.428 0.015278 *
#  dgrp3           15.581      7.443   2.093 0.036465 *
 # dgrp4            9.627      9.962   0.966 0.333987
#region2          7.052      6.858   1.028 0.303903
#region3         17.107      7.576   2.258 0.024064 *
#  region4         27.468      7.600   3.614 0.000309 ***
#  dgrp2:region2  -17.466      8.057  -2.168 0.030303 *
 # dgrp3:region2   -7.724      9.456  -0.817 0.414121
#dgrp4:region2   -1.140     16.609  -0.069 0.945307
#dgrp2:region3  -19.760      8.891  -2.222 0.026376 *
  #dgrp3:region3  -33.265     10.147  -3.278 0.001063 **
  #dgrp4:region3  -27.730     13.305  -2.084 0.037279 *
  #dgrp2:region4  -21.656      9.286  -2.332 0.019794 *
 # dgrp3:region4  -14.996     13.195  -1.137 0.255890
#dgrp4:region4   28.894     28.310   1.021 0.307560

#PC3 ~ dgrp*region
fit <- lm(PC3 ~ dgrp * region, data = pcas)
summary(fit)
#Coefficients:
 # Estimate Std. Error t value Pr(>|t|)
#(Intercept)     12.993      4.866   2.670 0.007644 **
#  dgrp2          -12.385      5.434  -2.279 0.022784 *
 # dgrp3          -19.908      5.586  -3.564 0.000375 ***
 # dgrp4          -15.075      7.476  -2.016 0.043911 *
 # region2        -10.925      5.147  -2.123 0.033913 *
 # region3        -11.056      5.686  -1.944 0.051997 .
#region4        -21.697      5.703  -3.804 0.000147 ***
#  dgrp2:region2    9.945      6.047   1.645 0.100228
#dgrp3:region2   27.405      7.097   3.862 0.000116 ***
#  dgrp4:region2    7.120     12.465   0.571 0.567938
#dgrp2:region3    9.535      6.673   1.429 0.153188
#dgrp3:region3   27.367      7.615   3.594 0.000334 ***
#  dgrp4:region3   15.561      9.985   1.558 0.119314
#dgrp2:region4   14.954      6.969   2.146 0.032018 *
 # dgrp3:region4   30.535      9.903   3.084 0.002075 **
 # dgrp4:region4   36.712     21.247   1.728 0.084178 .


tmp <- data1[!duplicated(data1[, c("cow")]), ]
with(tmp, table(dgrp, region))
  #        region
#dgrp   1   2   3   4
#1  49 417 134 131
#2 198 275 183 114
#3 154  71  64  20
#4  36  12  33   3
tmp2 <- prop.table(with(tmp, table(dgrp, region)))*100
tmp2

#----------------Scatter plot dEBV (EBV_GH-EBV_GHS) vs Spatial effect GHS------------------
#library(ggpmisc)
# Example data (assuming EBVs_Spatial_effect is your dataframe)
# Fit the linear model
model <- lm(dEBV ~ Spatial_effect_GHS, data = EBVs_Spatial_effect)

# Extract R-squared value
r_squared <- summary(model)$r.squared

# Create the scatter plot with regression line and R-squared annotation
scatter_plot <- ggplot(EBVs_Spatial_effect, aes(x = Spatial_effect_GHS, y = dEBV)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line without confidence interval
  labs(x = "Spatial effects (GHS)",
       y = "Difference between EBVs (GH-GHS)",
       color = "Spatial effects") +
  theme_bw() +
  annotate("text", x = Inf, y = Inf,
           label = paste("R² =", round(r_squared, 2)),
           hjust = 1.1, vjust = 2, size = 5, color = "black") +
  theme(
    axis.text.x = element_text(size = 15.3),  # Increase X-axis text size
    axis.text.y = element_text(size = 15.3),
    axis.title.x = element_text(size = 15.3),  # Increase X-axis label size
    axis.title.y = element_text(size = 15.3)   # Increase Y-axis label size
  ) +
  ylim(-0.5, 0.5)

print(scatter_plot)

ggsave(plot = scatter_plot, filename = "2_spatial_dEBV_scatter_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper

ggsave(plot = scatter_plot, filename = "2_spatial_dEBV_scatter_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


#------------Animal Ranking----------------------------------------------------

# Ranking EBV_G
EBV_G_Rank <- EBV_G[order(EBV_G$EBV_G,decreasing = TRUE),]

# Ranking EBV_fitGH
EBV_GH_Rank <- EBV_GH[order(EBV_GH$EBV_GH,decreasing = TRUE),]

# Ranking EBV_fitGS
EBV_GS_Rank <- EBV_GS[order(EBV_GS$EBV_GS,decreasing = TRUE),]

# Ranking EBV_fitGHS
EBV_GHS_Rank <- EBV_GHS[order(EBV_GHS$EBV_GHS,decreasing = TRUE),]

# Top 10 EBV_GH
EBV_GH_top10 <- EBV_GH_Rank[1:10,c("ID","EBV_GH")]
head(EBV_GH_top10, n=10)

# Top 10 EBV_GHS
EBV_GHS_top10 <- EBV_GHS_Rank[1:10,c("ID","EBV_GHS")]
head(EBV_GHS_top10, n=10)

sel_top10 <- EBV_GH_top10$ID %in% EBV_GHS_top10$ID
table(sel_top10)

#FALSE  TRUE
#4     6

# Top 20 EBV_GH
EBV_GH_top20 <- EBV_GH_Rank[1:20,c("ID","EBV_GH")]
head(EBV_GH_top20, n=20)

# Top 20 EBV_GHS
EBV_GHS_top20 <- EBV_GHS_Rank[1:20,c("ID","EBV_GHS")]
head(EBV_GHS_top20, n=20)

sel_top20 <- EBV_GH_top20$ID %in% EBV_GHS_top20$ID
table(sel_top20)
#FALSE  TRUE
#5    15


# Top 50 EBV_GH
EBV_GH_top50 <- EBV_GH_Rank[1:50,c("ID","EBV_GH")]

# Top 50 EBV_GHS
EBV_GHS_top50 <- EBV_GHS_Rank[1:50,c("ID","EBV_GHS")]


sel_top50 <- EBV_GH_top50$ID %in% EBV_GHS_top50$ID
table(sel_top50)
#FALSE  TRUE
#15      35


# Top 100 EBV_GH
EBV_GH_top100 <- EBV_GH_Rank[1:100,c("ID","EBV_GH")]

# Top 100 EBV_GHS
EBV_GHS_top100 <- EBV_GHS_Rank[1:100,c("ID","EBV_GHS")]


sel_top100 <- EBV_GH_top100$ID %in% EBV_GHS_top100$ID
table(sel_top100)
#FALSE  TRUE
30    70


# Top 1000 EBV_GH
EBV_GH_top1000 <- EBV_GH_Rank[1:1000,c("ID","EBV_GH")]

# Top 1000 EBV_GHS
EBV_GHS_top1000 <- EBV_GHS_Rank[1:1000,c("ID","EBV_GHS")]


sel_top1000 <- EBV_GH_top1000$ID %in% EBV_GHS_top1000$ID
table(sel_top1000)
#FALSE  TRUE
164   836

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



summary(fitG) 32716
summary(fitGH) # 32650
summary(fitGS) #32729
summary(fitGHS)# 32604
#------------GS_BOAT-------------
# Phenotype
pheno_boat <- read.csv(file= "data/original_data/pheno_boat.csv")
colnames(pheno_boat)
#[1] "region"               "district"             "ward"                 "village"              "Farm_id"
#[6] "farmergender"         "cattletotalowned"     "tag_id"               "animalid"             "closest_calvdate"
#[11] "milkdate"             "MilkAM"               "MilkMidDay"           "MilkPM"               "TotalMilk"
#[16] "Days.In.Milk"         "MilkFat"              "MilkProt"             "SCC"                  "Heartgirth"
#[21] "BodyLength"           "Weight"               "EstimatedWt"          "Bodyscore"            "parity"
#[26] "testdaynumber"        "latitude"             "longitude"            "original_tag_id"      "event_id"
#[31] "farmer_name"          "farm_id"              "project"              "birthdate"            "farmtype"
#[36] "Extract_Datetime_UTC" "Version"

length(pheno_boat$animalid) #319525 records
length(unique(pheno_boat$animalid)) #19,991 cows
length(unique(pheno_boat$farm_id))  #10,736 herds
table(pheno_boat$region)
#Arusha Dar es Salaam        Iringa        Kagera   Kilimanjaro          Mara         Mbeya      Morogoro
#61932            24         18308             8         81036             5         55790             3
#Njombe         Pwani        Simiyu         Tanga
#31262             3             9         71145
hist(pheno_boat$TotalMilk)

#Pedigree
pedigree_boat <- read.csv(file= "data/original_data/pedigree_boat.csv")
colnames(pedigree_boat)

length(pedigree_boat$animal_id) # 69901 in pedigree records
length(unique(pedigree_boat$animal_id)) # 69,901 animals

length(pedigree_boat$sire_id) # 69901 in pedigree records
length(unique(pedigree_boat$sire_id)) # 652 unique sires ID
length(unique(pedigree_boat$dam_id)) # 11162 dams ID

#--------------------Priors-for-Penalized-complexity----------------------------

#Help me understand by sharing the following:
# fitGS model
# Extract the posterior marginals for hyperparameters
hyperpar <- INLA::inla.hyperpar(fitGS)

# The list of hyperparameter names and marginals
hyperpar_names <- names(hyperpar$marginals.hyperpar)

# View the available hyperparameters
print(hyperpar_names)

#[1] "Precision for the Gaussian observations"
#[2] "Precision for animal"
#[3] "Range for field"
#[4] "Stdev for field"
#------plot of posterior distributions of hyperparameters from fitGS------------
# Genetic variance for G models
#Pg
# Step 1: Transform and combine data

data_GS <- data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x, fitGS$marginals.hyperpar$`Precision for animal`)))
data_GS$Source <- "GS"



# Step 2: Plot the combined data
Pg<- ggplot(data_GS, aes(x = x, y = y, colour = Source)) +
  geom_line() +
  labs(colour = "Model") +  scale_color_manual(values = c("GS"="black")) +

  theme_minimal() + labs(x="", y="", title ="Genomic variance (GS)") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "",
        legend.box.background = element_blank(),
        legend.title = element_text("")) +
  guides(colour=guide_legend(title="Model"))

# scale_y_continuous(breaks=NULL) to remove axis labels
Pg<- Pg+ theme(plot.title = element_text(face = "bold"))
Pg

# Herd variance
# Ph
# Step 1: Transform and combine data

data_GH <- data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x, fitGH$marginals.hyperpar$`Precision for herd`)))
data_GH$Source <- "GH"

data_GHS <- data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x, fitGHS$marginals.hyperpar$`Precision for herd`)))
data_GHS$Source <- "GHS"

combined_data <- rbind( data_GH, data_GHS)

# Residual variance
#Pr

data_GS <- data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x, fitGS$marginals.hyperpar$`Precision for the Gaussian observations`)))
data_GS$Source <- "GS"

# Step 2: Plot the combined data
Pr<- ggplot(data_GS, aes(x = x, y = y, colour = Source)) +
  geom_line() +
  theme_minimal() + labs(x="", y="", title ="Residual variance (GS)") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "",
        legend.box.background = element_blank(),
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model")) + scale_color_manual(values = c("GS"="black"))

# scale_y_continuous(breaks=NULL) to remove axis labels
Pr<- Pr+ theme(plot.title = element_text(face = "bold"))
Pr


# Plot spatial variance
varianceGS <- spde.posterior(fitGS, "field", "variance")

P_Svar<-  plot(varianceGS)

P_Svar <-P_Svar + ggtitle("Spatial Variance (GS)") +
  labs(x = "", y = "") + theme(plot.title = element_text(face = "bold")) + scale_color_manual(values = c("GS"="black"))
P_Svar

# Spatial range
# Plot spatial range
rangeGS <- spde.posterior(fitGS, "field", "range")
P_range <- plot(rangeGS)
P_range <- P_range + ggtitle("Spatial range (GS)") +
  labs(x = "", y = "") + theme(plot.title = element_text(face = "bold")) + geom_line(color="black")

P_range



(p <- ggarrange(Pg,Pr, P_Svar, P_range, ncol = 2, nrow = 3, common.legend = T, legend = "", align = "h", widths = c(1,1,1))) + PreseThemeNoLegend
p

#structure of the object fieldGS which I got from fieldGPS <- generate(fitGPS, spdf, ~field, 1000) where spdf is our data of class SpatialPointsDataFrame

str(fieldGS)












