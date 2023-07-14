#=======================================================================
# ADGG data analysis
# Trait: Test day milk yield
# Author: Isidore                                   #
# Version 1.0.0                                                        #
# Date: 26/10/2022                                                    #
# Descriptive analysis and data manipulation                           #                                                             #
#=======================================================================
#-----------------------------------------------------------------------
# Installing and Uploading packages   
#-----------------------------------------------------------------------
library(dplyr)
library(stringr)
library(ggplot2)     
library(ggExtra)
library(EnvStats)
library(Matrix)
library(lme4)     #for mixed model
library(afex)     #for p-values
library(readr)
library(GGally)
library(sf) # import shapefile
library(sp) # poly2nb for graph
library(foreach)
library(parallel)
library(INLA) # inla.read.graph
library(spData)
library(spdep)
library(tidyr)
library(rgdal)
library(terra)
library(devtools)
library(cowplot)
install.packages('spDataLarge', repos='https://nowosad.github.io/drat/',
                 type='source')
(.packages()) # Checking loaded packages

rm(list = ls())
# Setting working directory
setwd("C:/Users/Lenovo/OneDrive/Documents/adgg")
# Importing ADGG data
data<- read.table(file="data2.dat", header = TRUE) 
colnames(data2) <- c("cow","milk","ward","herd","cyrsn","tyrmn","dgrp","lac","lacgr","age","long","lat", "intercept", "leg1","leg2", "mean")
head(data2)

