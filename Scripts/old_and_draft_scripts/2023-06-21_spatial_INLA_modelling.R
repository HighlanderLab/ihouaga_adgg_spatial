################################################################################
# Geostatistical Models ADGG
################################################################################

# RUN MODELS (exp)

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
find_rtools() # TRUE
has_devel()
# Import data
rm(list = ls())
setwd("C:/Users/Lenovo/OneDrive/Documents/adgg")
getwd()
dir()
# Import data
data5<- read.table(file = "data5.dat", header = FALSE)
colnames(data5) <- c("cow","milk","ward","herd","cyrsn","tyrmn","dgrp","lac","lacgr","age","long","lat", "intercept", "leg1","leg2", "mean", "ward_code", "Region_Cod")
head(data5)
dim(data5)
# sample dataset
sample <- data5[sample(nrow(data5), 3000), ]
dim(sample)
# Import SNP data
SNP<- read.table(file = "snpref-isi.txt", header = FALSE)
head(SNP)

# write.csv(SNP_ID,"/Users/Lenovo/OneDrive/Documents/adgg/SNP.csv")
#########################################################################################################
# Import G matrix and make it a full matrix (Ginv2)
Ginv = data.frame(read.table("GLinv.txt", header = FALSE))
names(Ginv) = list("var1", "var2", "value")
dimension = max(Ginv[,1])   
Ginv2 = matrix(0, nrow = dimension, ncol = dimension)   

Ginv2[as.matrix(Ginv[c("var1", "var2")])] = Ginv[["value"]]  
Ginv2[as.matrix(Ginv[c("var2", "var1")])] = Ginv[["value"]] 
dim(Ginv2)

# Remove rows and columns with zeros in GinV2 (These are rows/columns numbers 18, 151,294,849 and 1428 )
# removing the 18th row and 18th colummn of the matrix
Ginv2_18<- Ginv2[-c(18), -c(18)]
Ginv2_18_151<- Ginv2_18[-c(150), -c(150)]
Ginv2_18_151_294<- Ginv2_18_151[-c(292), -c(292)]
Ginv2_18_151_294_849<-Ginv2_18_151_294[-c(846), -c(846)]
Ginv2_18_151_294_849_1428<- Ginv2_18_151_294_849[-c(1424), -c(1424)]
Ginv3 <-Ginv2_18_151_294_849_1428 # Matrix without missed animals
dim(Ginv3)
#--------------------------------------------------------------------------------
#Reorganising Cow'IDs'
#-------------------------------------------------------------------------------
# Cow_ID in dataset
cow_ID <- unique(data5$cow)
length(cow_ID) # 
cow_ID<- as.data.frame(cow_ID)
#SNP_ID in genotype data
SNP<- read.table(file = "snpref-isi.txt", header = FALSE)
head(SNP)
SNP_ID<- data.frame(unique(SNP[,1]))
length(SNP_ID) #
#IDs in G inverse (Ginv3)
ID_Ginv<-unique(Ginv$var1)
ID_Ginv<- as.data.frame(ID_Ginv)
length(ID_Ginv)
#  Let's add INLA_ID numbered from 1 to 1911 (number of genotyped cows) to Gin_IV_ID
INLA_Ginv_ID<-ID_Ginv 
INLA_Ginv_ID$INLA_ID <-row_number(INLA_Ginv_ID)
# Merging data5 with INLA_Ginv_ID
#rename IDG_inv to cow
names(INLA_Ginv_ID)[1] <- "cow"
# Merging IDs from Ginv with cow in my data (maintaining all my cows)
data6<- merge(x=data5, y=INLA_Ginv_ID, by="cow", all.x=TRUE)
length(unique(data6$cow)) # 1988
max(data6$cow)# 1911
max(data6$INLA_ID) # 1906

#-------------------------------------------------------------------------------

# Base model G


G1<- milk ~ 1 + ward +cyrsn + tyrmn + dgrp + (age|lacgr) + (intercept|lacgr) + (leg1|lacgr) +(leg2|lacgr)+ f(herd, model = 'iid') + f(cow, model = "generic0", Cmatrix = Ginv3) 

data5$cowPe <- data5$cow
G2<- milk ~ 1 + ward +cyrsn + tyrmn + dgrp + (age|lacgr) + (intercept|lacgr) + (leg1|lacgr) +(leg2|lacgr)+ f(herd, model = 'iid') + f(cowPe, model = 'iid') + f(cow, model = "generic0", Cmatrix = Ginv3) 

G3<- milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (intercept|lacgr) + (leg1|lacgr) +(leg2|lacgr)+ f(herd, model = 'iid') + f(ward, model = 'iid')+ f(cow, model = "generic0", Cmatrix = Ginv3) 

G4 <- milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (intercept|lacgr) + (leg1|lacgr) +(leg2|lacgr)+ f(herd, model = 'iid') + f(ward, model = 'iid')+ + f(cowPe, model = 'iid') + f(cow, model = "generic0", Cmatrix = Ginv3) 


modelG1<- inla(formula = G1,  control.compute = list(dic = TRUE), data=data5)
summary(modelG1)
modelG1$summary.random
modelG1$summary.hyperpar
  
modelG2<- inla(formula = G2, control.compute = list(dic = TRUE), data=data5)
summary(modelG2)
modelG2$summary.random
modelG2$summary.hyperpar


modelG3<- inla(formula = G3, control.compute = list(dic = TRUE), data=data5)
summary(modelG3)
modelG3$summary.random
modelG3$summary.hyperpar


modelG4<- inla(formula = G4, control.compute = list(dic = TRUE), data=data5)
summary(modelG4)
modelG4$summary.random
 modelG4$summary.hyperpar
 
################################################################################
 # Marginal variances This is described here (after figure 4.3) 
 #https://www.paulamoraga.com/book-geospatial/sec-inla.html
 
 # 1. Model G1 Mean posterior of variance and sd of herd
 modelG1$summary.hyperpar  
 
 # Mean posterior of Variance and sd of herd 
 marg.varianceG1_herd <- inla.tmarginal(function(x) 1/x,
                                   modelG1$marginals.hyperpar$`Precision for herd`)
 mG1_herd <- inla.emarginal(function(x) x, marg.varianceG1_herd)
 mG1_herd # 4.928808
 
 #sd_herd
 mG1_herdmG1_herd <- inla.emarginal(function(x) x^2, marg.varianceG1_herd)
 sd_herd=sqrt( mG1_herdmG1_herd - mG1_herd^2)
 sd_herd # 0.3987141
 
 ## Mean posterior of Variance and sd of cow (additive genetic variance and sd)
 
 marg.varianceG1_cow <- inla.tmarginal(function(x) 1/x,
                                        modelG1$marginals.hyperpar$`Precision for cow`)
 mG1_cow <- inla.emarginal(function(x) x, marg.varianceG1_cow)
 mG1_cow # 4.14846
 
 #sd_cow
 mG1_cowmG1_cow <- inla.emarginal(function(x) x^2, marg.varianceG1_cow)
 sd_cow=sqrt( mG1_cowmG1_cow - mG1_cow^2)
 sd_cow # 0.2457276 
 
 ## Mean posterior of Variance and sd of residual 
 
 marg.varianceG1_residual <- inla.tmarginal(function(x) 1/x,
                                       modelG1$marginals.hyperpar$`Precision for the Gaussian observations`)
 mG1_residual <- inla.emarginal(function(x) x, marg.varianceG1_residual)
 mG1_residual # 6.303613
 
 #sd_residual
 mG1_residualmG1_residual <- inla.emarginal(function(x) x^2, marg.varianceG1_residual)
 sd_residual=sqrt(mG1_residualmG1_residual - mG1_residual^2) 
 sd_residual #  
 
# We can also compute mean and sd of mean posterior of variance using  inla.zmarginal
# eg: residual
 inla.zmarginal(marg.varianceG1_residual) 
 #Mean            6.30361 
 #Stdev           0.0679027 
 
 
# A posterior of the variance of the random effect herd in G1

 ggplot(data.frame(inla.smarginal(marg.varianceG1_herd)), aes(x, y)) +
   geom_line() +
   theme_bw()

 # A posterior of the variance of the random effect herd in G4

 
 
 # Combine plot herd G1 G4 side by side

 
 # 2.Model G2 Mean posterior of variance and sd of herd
 
 modelG2$summary.hyperpar

 #herd model G2
marg.varianceG2_herd <- inla.tmarginal(function(x) 1/x,
                                            modelG2$marginals.hyperpar$`Precision for herd`)
 inla.zmarginal(marg.varianceG2_herd)
 #Mean            5.10843 
 #Stdev           0.48924 
 
 #additive (cow) model G2
 marg.varianceG2_cow <- inla.tmarginal(function(x) 1/x,
                                        modelG2$marginals.hyperpar$`Precision for cow`)
 inla.zmarginal(marg.varianceG2_cow)
 #Mean            3.19519 
 #Stdev           0.319685  
 
 #Permanent env. (cow)
 marg.varianceG2_cowPe <- inla.tmarginal(function(x) 1/x,
                                       modelG2$marginals.hyperpar$`Precision for cowPe`)
 inla.zmarginal(marg.varianceG2_cowPe)
 
#Mean            0.827498 
#Stdev           0.0843645 
 
#Residual variance
 
marg.varianceG2_residual <- inla.tmarginal(function(x) 1/x,
                                         modelG2$marginals.hyperpar$`Precision for the Gaussian observations`)
 inla.zmarginal(marg.varianceG2_residual)
 
 
# 3.Model G3 Mean posterior of variance and sd of herd
 
 
 modelG3$summary.hyperpar
 
 #herd model G3
 marg.varianceG3_herd <- inla.tmarginal(function(x) 1/x,
                                        modelG3$marginals.hyperpar$`Precision for herd`)
 inla.zmarginal(marg.varianceG3_herd)
 # Mean            2.18306 
 # Stdev           0.198157 
 
 
 #additive (cow) model G3
 marg.varianceG3_cow <- inla.tmarginal(function(x) 1/x,
                                       modelG3$marginals.hyperpar$`Precision for cow`)
 inla.zmarginal(marg.varianceG3_cow)
 
 
 #ward 
 marg.varianceG3_ward <- inla.tmarginal(function(x) 1/x,
                                         modelG3$marginals.hyperpar$`Precision for ward`)
 inla.zmarginal(marg.varianceG3_ward)
 
 
  #Residual variance model G3
 
 marg.varianceG3_residual <- inla.tmarginal(function(x) 1/x,
                                            modelG3$marginals.hyperpar$`Precision for the Gaussian observations`)
 inla.zmarginal(marg.varianceG3_residual)
 
 #Mean            6.3187 
 #Stdev           0.0673344
 
 
 # 4.Model G3 Mean posterior of variance and sd of herd
 
 modelG4$summary.hyperpar
 
 #herd model G4
 marg.varianceG4_herd <- inla.tmarginal(function(x) 1/x,
                                        modelG4$marginals.hyperpar$`Precision for herd`)
 inla.zmarginal(marg.varianceG4_herd)
 #Mean            2.15752 
 #Stdev           0.191723 
 
 
 #additive (cow) model G4
 marg.varianceG4_cow <- inla.tmarginal(function(x) 1/x,
                                       modelG4$marginals.hyperpar$`Precision for cow`)
 inla.zmarginal(marg.varianceG4_cow)
 
 #Mean            2.39907 
 #Stdev           0.200297 
 
 #Permanent env. (cow) model G4
 marg.varianceG4_cowPe <- inla.tmarginal(function(x) 1/x,
                                       modelG4$marginals.hyperpar$`Precision for cowPe`)
 inla.zmarginal(marg.varianceG4_cowPe)
 #Mean            0.0103343 
# Stdev           0.00314107 
 
 #ward G4
 marg.varianceG4_ward <- inla.tmarginal(function(x) 1/x,
                                        modelG4$marginals.hyperpar$`Precision for ward`)
 inla.zmarginal(marg.varianceG4_ward)
# Mean            6.387 
# Stdev           0.948287 
 
#Residual variance model G4
 
marg.varianceG4_residual <- inla.tmarginal(function(x) 1/x,
                                            modelG4$marginals.hyperpar$`Precision for the Gaussian observations`)
 inla.zmarginal(marg.varianceG4_residual)
 #Mean            6.31695 
 #Stdev           0.0683571 
 #
 
 ###############################################################################
 #Fitting G1, G2, G3, G4 and Besag with ward_code
 ###############################################################################
 G1.1<- milk ~ 1 + ward_code +cyrsn + tyrmn + dgrp + (age|lacgr) + (intercept|lacgr) + (leg1|lacgr) +(leg2|lacgr)+ f(herd, model = 'iid') + f(cow, model = "generic0", Cmatrix = Ginv3) 
 
 G2.1<- milk ~ 1 + ward_code +cyrsn + tyrmn + dgrp + (age|lacgr) + (intercept|lacgr) + (leg1|lacgr) +(leg2|lacgr)+ f(herd, model = 'iid') + f(cowPe, model = 'iid') + f(cow, model = "generic0", Cmatrix = Ginv3) 
 
 G3.1<- milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (intercept|lacgr) + (leg1|lacgr) +(leg2|lacgr)+ f(herd, model = 'iid') + f(ward_code, model = 'iid')+ f(cow, model = "generic0", Cmatrix = Ginv3) 
 
 G4.1 <- milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (intercept|lacgr) + (leg1|lacgr) +(leg2|lacgr)+ f(herd, model = 'iid') + f(ward_code, model = 'iid')+ + f(cowPe, model = 'iid') + f(cow, model = "generic0", Cmatrix = Ginv3) 

 ###
 modelG1.1<- inla(formula = G1.1,  control.compute = list(dic = TRUE), data=data5)
 summary(modelG1.1)
 modelG1.1$summary.random
 modelG1.1$summary.hyperpar
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 modelG2.1<- inla(formula = G2.1, control.compute = list(dic = TRUE), data=data5)
 summary(modelG2.1)
 modelG2.1$summary.random
 modelG2.1$summary.hyperpar
 
 
 modelG3.1<- inla(formula = G3.1, control.compute = list(dic = TRUE), data=data5)
 summary(modelG3.1)
 modelG3.1$summary.random
 modelG3.1$summary.hyperpar
 
 
 modelG4.1<- inla(formula = G4.1, control.compute = list(dic = TRUE), data=data5)
 summary(modelG4.1)
 modelG4.1$summary.random
 modelG4.1$summary.hyperpar
 
 length(unique(data5$cow))
 
################################################################################
#Besag model
###############################################################################
 # Importing map ward Tanzania
 mapwa2 <- st_read("/Users/Lenovo/OneDrive/Documents/adgg/TZwards.shp", quiet = TRUE)
 class(mapwa2) # class sf
 head(mapwa2)
 plot(mapwa)
 ggplot(mapwa2) + geom_sf()+ geom_sf_text(data = mapwa2, aes(label = Ward_Code), size = 3, colour = "black")
 
 # construct graph of neighbours (Library INLA)
 library(spdep) # To identify if two regions are neighbours
 #remotes::install_github("r-spatial/s2")
 library(s2)
 sf_use_s2(FALSE)
 
 nb.mapwa2 <- poly2nb(mapwa2)# Look at map and tell neighbours
 nb2INLA("map.graphwa2",nb.mapwa2)
 gwa <- inla.read.graph(filename = "map.graphwa2") # Lines 322 and 323 run together
 
#Formula Ward Besag
#GBesag 1
 
G4_Besag <- milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (intercept|lacgr) + (leg1|lacgr) +(leg2|lacgr)+ f(herd, model = 'iid') + f(cowPe, model = 'iid') + f(cow, model = "generic0", Cmatrix = Ginv3) +  f(ward_code, model = "besag", graph = gwa, scale.model = TRUE)
 
 
modelG4_Besag <- inla(formula = G4_Besag, control.compute = list(dic = TRUE), data=data5)
summary(modelG4_Besag)
modelG4_Besag$summary.random
modelG4_Besag$summary.hyperpar
 
 
############################################
# calculate distances between locations
############################################
#Rcpp::sourceCpp("data/other/euclideandistance.cpp")
# DISTANCES==
rhoHS <- 27.5
rhoS <- 10
# Exponential distances
# small data
fun_distances <- function(data5, rho, model = "S"){
  data5 <- data5[order(data5$herd),] # 
  coordinate <- unique(paste(data5$long, data$lat))
  coordinate <- matrix(as.numeric(unlist(strsplit(coordinate, split=" "))), ncol=2, byrow = T)
  paired.distances <- fastPairDist(coordinate, coordinate)
  paired.distances.exp <- exp(-paired.distances/rho)
  Exp.inv <- solve(paired.distances.exp) # inverz matrike
  if (model == "S"){
    data5$idExpS <- as.numeric(as.factor(data5$herd))
  }
  if (model == "HS"){
    data5$idExpHS <- as.numeric(as.factor(data5$herd))
  }
  
  list(data5 = data5, dist = paired.distances.exp, dist.inv = Exp.inv)
}

dist.small <- fun_distances(phenoData.small, rhoS, "S")
phenoData.small <- dist.small[[1]]
Exp.inv.small.S <- dist.small[[3]]

dist.small <- fun_distances(phenoData.small, rhoS, "HS")
phenoData.small <- dist.small[[1]]
Exp.inv.small.HS <- dist.small[[3]]


phenoData.small <- phenoData.small[order(phenoData.small$index), ]


############################################
# define priors
############################################
hyperVarPed = list(theta = list(prior="pc.prec", param=c(sqrt(0.1),0.5)))
hyperVarIdlok.GH = list(theta = list(prior="pc.prec", param=c(sqrt(0.25),0.5)))
hyperVarIdlok.GHS = list(theta = list(prior="pc.prec", param=c(sqrt(0.15),0.5)))
hyperResVarG = list(theta = list(prior="pc.prec", param=c(sqrt(0.3),0.5)))
hyperResVarGHS = list(theta = list(prior="pc.prec", param=c(sqrt(0.15),0.5)))

# Exp
hyperVarExp.S = list(theta = list(prior="pc.prec", param=c(sqrt(0.25),0.5)))
hyperVarExp.HS = list(theta = list(prior="pc.prec", param=c(sqrt(0.10),0.5)))



formulaB.S.Exp.prior.small = last004scaledNA ~ 1 + vpliv1 + vpliv2 + vpliv3  +  
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed) + 
  f(idExpS, model = "generic0", Cmatrix = Exp.inv.small.S, hyper = hyperVarExp.S)
formulaB.HS.Exp.prior.small <- last004scaledNA ~ 1 + vpliv1 + vpliv2 + vpliv3 + 
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed) + 
  f(idlok, model = "iid", hyper = hyperVarIdlok.GHS) + 
  f(idExpHS, model = "generic0", Cmatrix = Exp.inv.small.HS, hyper = hyperVarExp.HS)




model.Exp.S.prior.small.NA <- inla(formulaB.S.Exp.prior.small,  family = "gaussian", 
                                   data = phenoData.small, control.compute = list(dic=T,cpo=F,config=T), verbose=T)
model.Exp.HS.prior.small.NA <- inla(formulaB.HS.Exp.prior.small,  family = "gaussian", 
                                    data = phenoData.small, control.compute = list(dic=T,cpo=F,config=T), verbose=T)








########################################################################################
#G0 Base model
#####################################################################################
G0 <- milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (intercept|lacgr) + (leg1|lacgr) +(leg2|lacgr) + f(ward, model = 'iid')+ + f(cowPe, model = 'iid') + f(cow, model = "generic0", Cmatrix = Ginv3)

modelG0<- inla(formula = G0, control.compute = list(dic = TRUE), data=data5)
summary(modelG0)
modelG0$summary.random
modelG0$summary.hyperpar



