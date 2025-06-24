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

# Read in the genomic relationship matrix
load(file = "data/cleaned_data/GRMInv.RData")
str(GRMInv)
dim(GRMInv) # 1894 x 1894
class(GRMInv)
head(GRMInv)
# Now we can code the cows in pheno data correctly
data1$cow <- factor(data1$cow, levels = 1:nrow(GRMInv))
summary(as.numeric(levels(data1$cow))) # we need 1:1894 numbers here!
data1$cowI <- as.numeric(as.character(data1$cow))
summary(data1$cowI) # we have 1:1894 ids here
head(data1)
tail(data1)
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

#Mesh plot for manuscript
meshplot<- ggplot() +
  geom_fm(data = mesh) +
  geom_sf(data = dat, size = 0.1) + geom_sf(data=data_border,alpha=0.2) +
  coord_sf(datum = crs_tz) +
  labs(x = "Longitude (km)", y = "Latitude (km)")
meshplot

#Formating of plot
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
# Export the mesh in Paper_plot format.
ggsave(plot = meshplot + PaperTheme, filename = "Figure1New.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper

#For INLABru define matern
matern <-
  inla.spde2.pcmatern(mesh,
                      prior.sigma = c(sqrt(0.25), 0.5),
                      prior.range = c(50, 0.8)
  )

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
#-------Define G models (without permanent environmental effect)--------------
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

# ModelGHS
modelGHS <- ~ fixed_effects(main = ~ tyrmn + cyrsn + dgrp + ageZ:lacgr + leg1:lacgr + leg2:lacgr, model = "fixed") +
  animal(main = cow, model = 'generic', Cmatrix = GRMInv) +
  herd(main=herd, model = 'iid') +
  field(geometry, model = matern)
modelGHSFormula <- milkZ ~ .

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
# Run fitGPH
fitGPH <- bru(modelGPH,
             like(family = "Gaussian",
                  modelGPHFormula,
                  data = data2))

# Run fitGPS
fitGPS <- bru(modelGPS,
             like(family = "Gaussian",
                  modelGPSFormula,
                  data = data2))

#--------Run G models-----------------------------------------------------------
# Run fitG
fitG<- bru(modelG,
           like(family = "Gaussian",
                modelGFormula,
                data = data2))
summary(fitG)
summarise_precision_to_variance(fitG)

# Run fitGH
fitGH<- bru(modelGH,
            like(family = "Gaussian",
                 modelGHFormula,
                 data = data2))
# Run fitGS
fitGS<- bru(modelGS,
            like(family = "Gaussian",
                 modelGSFormula,
                 data = data2))
# Run fitGHS
fitGHS<- bru(modelGHS,
             like(family = "Gaussian",
                  modelGHSFormula,
                  data = data2))

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
round(mean(column_variances_GPS$spatial_var),3)
round(sd(column_variances_GPS$spatial_var),3) 

Range_fieldGPS <- generate(fitGPS,spdf, ~Range_for_field, 1000)
str(Range_fieldGPS)
Range_fieldGPS_df <- data.frame(Range_fieldGPS)
round(rowMeans(Range_fieldGPS_df),3) 
# Define a function to calculate row-wise standard deviation
row_sd <- function(data) {
  apply(data, 1, sd, na.rm = TRUE)
}
sd_range_GPS <- round(row_sd(Range_fieldGPS_df),3)
sd_range_GPS #

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
round(mean(column_variances_GPHS$spatial_var),3) 
round(sd(column_variances_GPHS$spatial_var),3)


#Range of spatial effect fitGPHS
Range_fieldGPHS <- generate(fitGPHS,spdf, ~Range_for_field, 1000)

Range_fieldGPHS_df <- data.frame(Range_fieldGPHS)
round(rowMeans(Range_fieldGPHS_df),3)
sd_range_GPHS <- round(row_sd(Range_fieldGPHS_df),3)
sd_range_GPHS

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
round(mean(column_variances_GS$spatial_var),3) 
round(sd(column_variances_GS$spatial_var),3) 

Range_fieldGS <- generate(fitGS,spdf, ~Range_for_field, 1000)
str(Range_fieldGS)
Range_fieldGS_df <- data.frame(Range_fieldGS)
round(rowMeans(Range_fieldGS_df),3) 
# Define a function to calculate row-wise standard deviation
row_sd <- function(data) {
  apply(data, 1, sd, na.rm = TRUE)
}
sd_range_GS <- round(row_sd(Range_fieldGS_df),3)
sd_range_GS 

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
round(mean(column_variances_GHS$spatial_var),3) 
round(sd(column_variances_GHS$spatial_var),3)  

#Range of spatial effect fitGHS
Range_fieldGHS <- generate(fitGHS,spdf, ~Range_for_field, 1000)
Range_fieldGHS_df <- data.frame(Range_fieldGHS)
round(rowMeans(Range_fieldGHS_df),3)
sd_range_GHS <- round(row_sd(Range_fieldGHS_df),3)
sd_range_GHS

# ----------------Plotting Spatial effects from GHS and GS models------
#fitGHS
# Load required packages
library(ggplot2)
library(sf)
library(gridExtra)
library(viridis)
#fitGHS (Manuscript Figure 4)

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

#fitGS (Manuscript Additional file 2)

gproj <- inla.mesh.projector(mesh,  dims = c(300,300))

g.mean <- inla.mesh.project(gproj, fitGS$summary.random$field$mean)
g.sd <- inla.mesh.project(gproj, fitGS$summary.random$field$sd)
summary(fitGS)
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
