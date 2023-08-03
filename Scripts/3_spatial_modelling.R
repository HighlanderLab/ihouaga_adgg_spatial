# ---- Header ------------------------------------------------------------------

# Spatial modeling
# Trait: Test day milk yield
# Author: Isidore
# Version 1.0.0
# Date: 2023-07-20

# ---- Setup -------------------------------------------------------------------

# Working directory
# ... Isidore's laptop
baseDir <- "C:/Users/Lenovo/OneDrive/Documents/ihouaga_adgg_spatial"
getwd()
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
  )
  install.packages(pkgs = requiredPackages)
  install.packages(pkgs = "INLA",
                   repos=c(getOption("repos"),
                           INLA = "https://inla.r-inla-download.org/R/stable"),
                   dep = TRUE)
  inla.upgrade(testing = TRUE)
}
library(tidyverse)
library(INLA) # TODO change to inlabru and add inlabrue estimation above

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
summary(data1$herdI)
data1$herdI <- as.numeric(data1$herd) # these codes will now be 1:n
data1$ward_codeI <- as.numeric(data1$ward_code) # these codes will now be 1:n
data1$cowPe <- data1$cow # cow permanent environment
data1$cowPeI <- as.numeric(data1$cowPe) # these codes will now be 1:n
# data1$cowI <- as.numeric(data1$????) # see below!

# Read in the regional data for Besag model
nb.map <- inla.read.graph(filename = "data/cleaned_data/ward_neighbours.txt")
load(file = "data/cleaned_data/ward_neighbours_precision_matrix.RData")
# this will load nb.matrix and nb.matrixScaled

# Read in the genomic relationship matrix
load(file = "data/cleaned_data/GRMInv.RData")
str(GRMInv)
dim(GRMInv) # 1911 x 1911
class(GRMInv)
# Now we can code the cows in pheno data correctly
data1$cow <- factor(data1$cow, levels = 1:nrow(GRMInv))
summary(as.numeric(levels(data1$cow))) # we need 1:1911 numbers here!
data1$cowI <- as.numeric(as.character(data1$cow))
summary(data1$cowI) # we have 1:1906 ids here
head(data1)
tail(data1)
# TODO double check that the above cowI is correct (GG thinks it is )

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

# ---- Specify models ----------------------------------------------------------

# Base model for milk - All effects without ward_codeI 
modelBase <- "milkZ ~ 1 + cyrsn + tyrmn + dgrp + (ageZ|lacgr) + (leg0|lacgr) + (leg1|lacgr) + (leg2|lacgr) + f(herdI, model = 'iid') + f(cowPeI, model = 'iid') + f(cowI, model = 'generic0', Cmatrix = GRMInv)"

# Adding ward_codeI as a fixed effect
modelWCF <- as.formula(paste0(modelBase, " + ward_codeI"))

# Adding ward_codeI as a random IID effect
modelWCRI <- as.formula(paste0(modelBase, " + f(ward_codeI, model = 'iid')"))

# Adding ward_codeI as a random Besag effect
modelWCRB <- as.formula(paste0(modelBase, " + f(ward_codeI, model = 'besag', graph = nb.map, scale.model = TRUE)"))

# * base model for milk_corrected
modelBase_milk_corrected <- "milkZ ~ 1 + cyrsn + tyrmn + dgrp + (ageZ|lacgr) + (leg0|lacgr) + (leg1|lacgr) + (leg2|lacgr) + f(cowPeI, model = 'iid')" # Milk corrected for all the fixed effects + cowPeI as in Mrode et al.2021

# * model for milk_corrected: Adding ward_code as a fixed effect
modelWCF_milk_corrected <- as.formula(paste0(modelBase_milk_corrected, " + ward_codeI"))

# * model for milk_corrected: Adding ward_code as a random IID effect
modelWCRI_milk_corrected <- as.formula(paste0(modelBase_milk_corrected, " + f(ward_codeI, model = 'iid')"))

# * model for milk_corrected: Adding ward_code as a random Besag effect
modelWCRB_milk_corrected <- as.formula(paste0(modelBase_milk_corrected, " + f(ward_codeI, model = 'besag', graph = nb.map, scale.model = TRUE)"))

# ---- Run the models ----------------------------------------------------------
# *fitWCF
fitWCF <- inla(formula = modelWCF, data = data1,
               control.compute = list(dic = TRUE))
(tmp <- summary(fitWCF))
tmp2 <- 1 / tmp$hyperpar$mean
names(tmp2) <- sub(x = rownames(tmp$hyperpar),
                   pattern = "Precision",
                   replacement = "Variance")
tmp2

# fitWCRI
fitWCRI <- inla(formula = modelWCRI, data = data1,
               control.compute = list(dic = TRUE))
(tmp <- summary(fitWCRI))
tmp2 <- 1 / tmp$hyperpar$mean
names(tmp2) <- sub(x = rownames(tmp$hyperpar),
                   pattern = "Precision",
                   replacement = "Variance")
tmp2

# fitWCRB
fitWCRB <- inla(formula = modelWCRB, data = data1,
               control.compute = list(dic = TRUE))
(tmp <- summary(fitWCRB))
tmp2 <- 1 / tmp$hyperpar$mean
names(tmp2) <- sub(x = rownames(tmp$hyperpar),
                   pattern = "Precision",
                   replacement = "Variance")
tmp2

# *fitWCF_milk_corrected
fitWCF_milk_corrected<- inla(formula = modelWCF_milk_corrected, data = data1,
               control.compute = list(dic = TRUE))
(tmp <- summary(fitWCF_milk_corrected))
tmp2 <- 1 / tmp$hyperpar$mean
names(tmp2) <- sub(x = rownames(tmp$hyperpar),
                   pattern = "Precision",
                   replacement = "Variance")
tmp2

# fitWCRI_milk_corrected
fitWCRI_milk_corrected <- inla(formula = modelWCRI_milk_corrected, data = data1,
                control.compute = list(dic = TRUE))
(tmp <- summary(fitWCRI_milk_corrected))
tmp2 <- 1 / tmp$hyperpar$mean
names(tmp2) <- sub(x = rownames(tmp$hyperpar),
                   pattern = "Precision",
                   replacement = "Variance")
tmp2

# fitWCRB_milk_corrected
fitWCRB_milk_corrected <- inla(formula = modelWCRB_milk_corrected, data = data1,
                control.compute = list(dic = TRUE))
(tmp <- summary(fitWCRB_milk_corrected))
tmp2 <- 1 / tmp$hyperpar$mean
names(tmp2) <- sub(x = rownames(tmp$hyperpar),
                   pattern = "Precision",
                   replacement = "Variance")
tmp2

# -------SPDE-------------------------------------------------------------------
























#===========================================================================================================================================================
# ADGG-Modelling
# Trait: Test day milk yield
# Author: Isidore Houaga                                  #
# Version 1.0.0                                                        #
# Date: 12/07/2023                                                    #
#===========================================================================================================================================================################################################################################
#Loading Libraries

library(INLA)
library(sf) #spatial
library(sp)
library(spData)
library(spdep) #spatial
library(rgdal) #spatial
library(tidyverse)
library(Matrix)
library(foreach)
library(parallel)
library(kableExtra)
library(gdata) # drop unused factor levels
library(patchwork) # combine plots
library(naniar) # missing values
library(ggpubr)
library(ggplot2)
library(usethis)
library(devtools)
library(RcppArmadillo)
library(raster)
library(dplyr)
find_rtools() # TRUE
has_devel()

# sample dataset
sample <- data5[sample(nrow(data5), 3000), ]
dim(sample)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Base model G
G <- milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (leg0|lacgr) + (leg1|lacgr) + (leg2|lacgr) + f(herd, model = 'iid') + f(cowPe, model = 'iid') + f(IId, model = "generic0", Cmatrix = Ginv2)

YDG <- milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (leg0|lacgr) + (leg1|lacgr) +(leg2|lacgr) + f(cowPe, model = 'iid')
YDG# Corrected phenotype (Yield deviation) for model G
# Alternative models

# G + ward_code as fixed
G1 <- milk ~ 1 + ward_code + cyrsn + tyrmn + dgrp + (age|lacgr) + (leg0|lacgr) + (leg1|lacgr) +(leg2|lacgr)+ f(herd, model = 'iid') + f(cowPe, model = 'iid') + f(IId, model = "generic0", Cmatrix = Ginv2)

YDG1 <- milk ~ 1 + ward_code +cyrsn + tyrmn + dgrp + (age|lacgr) + (leg0|lacgr) + (leg1|lacgr) +(leg2|lacgr) + f(cowPe, model = 'iid')
 YDG1 # Corrected phenotype (Yield deviation) for model G1

# G + ward_code as random iid
G2 <- milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (leg0|lacgr) + (leg1|lacgr) +(leg2|lacgr)+ f(herd, model = 'iid') + f(ward_code, model = 'iid') + f(cowPe, model = 'iid') + f(IId, model = "generic0", Cmatrix = Ginv2)

YDG2 <-milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (leg0|lacgr) + (leg1|lacgr) +(leg2|lacgr) + f(ward_code, model = 'iid') + f(cowPe, model = 'iid')
YDG2 # Corrected phenotype (Yield deviation) for model G2

# G + ward_code as random besag
G3<- milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (leg0|lacgr) + (leg1|lacgr) +(leg2|lacgr)+ f(herd, model = 'iid') + f(cowPe, model = 'iid') + f(IId, model = "generic0", Cmatrix = Ginv2) +  f(ward_code, model = "besag", graph = gwa, scale.model = TRUE)
# G3 is the Besag model
YDG3<- milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (leg0|lacgr) + (leg1|lacgr) +(leg2|lacgr) + f(cowPe, model = 'iid') + f(ward_code, model = "besag", graph = gwa, scale.model = TRUE)
YDG3 # Corrected phenotype (Yield deviation) for model G3

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Model G
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variance components
modelG<- inla(formula = G,  control.compute = list(dic = TRUE), data=pheno)# Model G
summary(modelG)
modelG$summary.random
modelG$summary.hyperpar
modelG$summary.random$herd
#herd model G
marg.varianceG_herd <- inla.tmarginal(function(x) 1/x,
                                       modelG$marginals.hyperpar$`Precision for herd`)
inla.zmarginal(marg.varianceG_herd)
#Mean            9.54961
#Stdev           0.207374

#additive (cow) ie IId in model G
marg.varianceG_IId <- inla.tmarginal(function(x) 1/x,
                                      modelG$marginals.hyperpar$`Precision for IId`)
inla.zmarginal(marg.varianceG_IId)
#Mean            1.09332
#Stdev           0.190655

#Permanent env. (cow)
marg.varianceG_cowPe <- inla.tmarginal(function(x) 1/x,
                                        modelG$marginals.hyperpar$`Precision for cowPe`)
inla.zmarginal(marg.varianceG_cowPe)

#Mean            1.20414
#Stdev           0.15874

#Residual variance

marg.varianceG_residual <- inla.tmarginal(function(x) 1/x,
                                           modelG$marginals.hyperpar$`Precision for the Gaussian observations`)
inla.zmarginal(marg.varianceG_residual)

#Mean            6.32345
#Stdev           0.0257729

# Phenotype variance

PhenovarG <- marg.varianceG_herd + marg.varianceG_IId + marg.varianceG_cowPe + marg.varianceG_residual
inla.zmarginal(PhenovarG)

#Mean            18.2094
#Stdev           0.587476

#HeritabilityG

h2G<-  1.09332/18.2387
round(h2G,2) #   0.06
sdh2G<- 0.190655/ 0.587476
round(sdh2G,2) # 0.32

# Ploting  posterior of the variance of the random effect herd in model G

ggplot(data.frame(inla.smarginal(marg.varianceG_herd)), aes(x, y)) +
  geom_line() +
  theme_bw()

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Model G1
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
modelG1<- inla(formula = G1,  control.compute = list(dic = TRUE), data=pheno)#Model G1
summary(modelG1)
modelG1$summary.random
modelG1$summary.hyperpar

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Model G2
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
modelG2<- inla(formula = G2, control.compute = list(dic = TRUE), data=pheno) #model G2
summary(modelG2)
modelG2$summary.random
modelG2$summary.hyperpar

#herd model G2
marg.varianceG2_herd <- inla.tmarginal(function(x) 1/x,
                                       modelG2$marginals.hyperpar$`Precision for herd`)
inla.zmarginal(marg.varianceG2_herd)
#Mean           5.83766
#Stdev          0.358565

#additive (cow) model G2
marg.varianceG2_IId <- inla.tmarginal(function(x) 1/x,
                                      modelG2$marginals.hyperpar$`Precision for IId`)
inla.zmarginal(marg.varianceG2_IId)
#Mean           0.795031
#Stdev           0.209466

#Permanent env. (cow)
marg.varianceG2_cowPe <- inla.tmarginal(function(x) 1/x,
                                        modelG2$marginals.hyperpar$`Precision for cowPe`)
inla.zmarginal(marg.varianceG2_cowPe)

#Mean            1.52253
#Stdev          0.283257

#Residual variance

marg.varianceG2_residual <- inla.tmarginal(function(x) 1/x,
                                           modelG2$marginals.hyperpar$`Precision for the Gaussian observations`)
inla.zmarginal(marg.varianceG2_residual)
#Mean            6.30954
#Stdev           0.0674283

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Model G3
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
modelG3<- inla(formula = G3, control.compute = list(dic = TRUE), data=pheno) #Besag
summary(modelG3)
modelG3$summary.random
modelG3$summary.hyperpar

# Model G3 Mean posterior of variance and sd of herd


modelG3$summary.hyperpar

#herd model G3
marg.varianceG3_herd <- inla.tmarginal(function(x) 1/x,
                                       modelG3$marginals.hyperpar$`Precision for herd`)
inla.zmarginal(marg.varianceG3_herd)
# Mean            5.91992
# Stdev           0.331839


#additive (cow) model G3
marg.varianceG3_IId <- inla.tmarginal(function(x) 1/x,
                                      modelG3$marginals.hyperpar$`Precision for IId`)
inla.zmarginal(marg.varianceG3_IId)

#Mean            0.878264
#Stdev           0.25309

#ward_code
marg.varianceG3_ward_code <- inla.tmarginal(function(x) 1/x,
                                       modelG3$marginals.hyperpar$`Precision for ward_code`)
inla.zmarginal(marg.varianceG3_ward_code)
#Mean            10.2872
#Stdev           1.65742

#Residual variance model G3

marg.varianceG3_residual <- inla.tmarginal(function(x) 1/x,
                                           modelG3$marginals.hyperpar$`Precision for the Gaussian observations`)
inla.zmarginal(marg.varianceG3_residual)

#Mean            6.32046
#Stdev           0.0604557
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#MODEL ACCURACY (ALL ANIMALS)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#ACCURACY MODEL G (ALL ANIMALS)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
YDG <- milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (leg0|lacgr) + (leg1|lacgr) +(leg2|lacgr) + f(cowPe, model = 'iid')
modelYDG<- inla(formula = YDG,  control.compute = list(dic = TRUE), data=pheno)
YDG <- data.frame(modelYDG$summary.random$cowPe)
YDG <- subset(YDG, select= c(ID,mean))
names(YDG)[1]<- "INLA_cowID"
names(YDG)[2]<- "PhenoCor"

EBVG <- modelG$summary.random$cow

# Remove non phenotyped animals from EBV
head(EBVG)
head(YDG)

sel2 <- EBVG$ID %in% YDG$INLA_cowID
EBVG <- EBVG[sel2, ]
dim(EBVG)
dim(YDG)
round(cor(EBVG$mean,YDG$PhenoCor),2) # 0.76
Coef<- lm (YDG$PhenoCor~EBVG$mean)
summary(Coef) # R2=0.572 (overstimation)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#ACCURACY MODEL G1 (ALL ANIMALS)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
YDG1<- milk ~ 1 + ward_code +cyrsn + tyrmn + dgrp + (age|lacgr) + (leg0|lacgr) + (leg1|lacgr) +(leg2|lacgr) + f(cowPe, model = 'iid')
modelYDG1<- inla(formula = YDG1,  control.compute = list(dic = TRUE), data=pheno)
YDG1 <- data.frame(modelYDG1$summary.random$cowPe)
YDG1 <- subset(YDG1, select= c(ID,mean))
names(YDG1)[1]<- "INLA_cowID"
names(YDG1)[2]<- "PhenoCor"


EBVG1 <- modelG1$summary.random$cow


# Remove non phenotyped animals from EBV

head(YDG1)
sel3 <- EBVG1$ID %in% YDG1$INLA_cowID
EBVG1 <- EBVG1[sel3, ]
dim(EBVG1)
dim(YDG1)
round(cor(EBVG1$mean, YDG1$PhenoCor),2)# 0.76
Coef1<- lm (YDG1$PhenoCor~EBVG1$mean)
summary(Coef1) # R2= 0.572 (Overstimation)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#ACCURACY MODEL G2 (ALL ANIMALS)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
YDG2 <-milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (leg0|lacgr) + (leg1|lacgr) +(leg2|lacgr) + f(ward_code, model = 'iid') + f(cowPe, model = 'iid')
modelYDG2<- inla(formula = YDG2,  control.compute = list(dic = TRUE), data=pheno)
YDG2 <- data.frame(modelYDG2$summary.random$cowPe)
YDG2 <- subset(YDG2, select= c(ID,mean))
names(YDG2)[1]<- "INLA_cowID"
names(YDG2)[2]<- "PhenoCor"


EBVG2 <- modelG2$summary.random$cow


# Remove non phenotyped animals from EBV

head(YDG2)
sel4 <- EBVG2$ID %in% YDG2$INLA_cowID
EBVG2 <- EBVG2[sel4, ]
dim(EBVG2)
dim(YDG2)

round(cor(EBVG2$mean, YDG2$PhenoCor),2)#  0.84

Coef2<- lm (YDG2$PhenoCor~EBVG2$mean)
summary(Coef2) # R2=0.6978 (Overstimation)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#ACCURACY MODEL G3 (ALL ANIMALS)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
YDG3<- milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (leg0|lacgr) + (leg1|lacgr) + (leg2|lacgr) + f(cowPe, model = 'iid') + f(ward_code, model = "besag", graph = gwa, scale.model = TRUE)
modelYDG3<- inla(formula = YDG3,  control.compute = list(dic = TRUE), data=pheno)
YDG3 <- data.frame(modelYDG3$summary.random$cowPe)
YDG3 <- subset(YDG3, select= c(ID,mean))
names(YDG3)[1]<- "INLA_cowID"
names(YDG3)[2]<- "PhenoCor"

EBVG3 <- modelG3$summary.random$cow


# Remove non phenotyped animals from EBV

head(YDG3)
sel5 <- EBVG3$ID %in% YDG3$INLA_cowID
EBVG3 <- EBVG3[sel5, ]
dim(EBVG3)
dim(YDG3)
round(cor(EBVG3$mean, YDG3$PhenoCor),2)# 0.84
Coef3<- lm (YDG3$PhenoCor~EBVG3$mean)
summary(Coef3)# R2model3= 0.6995 Overestimation


 #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 #Prepare data for Raphael
 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

 pheno_isi <- pheno %>% relocate("IId", .before=cow)
 dim(pheno_isi)
 pheno_isi<- subset(pheno_isi,select= -c(cow))
 names(pheno_isi)[1]<- "cow"
 dim(pheno_isi)

write.table(pheno_isi, "pheno_isi.txt",
             quote = FALSE, row.names = FALSE, col.names = TRUE, sep = " ", na = "0")

# pheno_isi and GLinv.txt (G inverse (1911 x 1911)) were shared with Raphael for ASReml test.
