# ---- Header ------------------------------------------------------------------
# Spatial modeling
# Trait: Test day milk yield
# Author: Isidore
# Version 1.0.0
# Date: 2023-07-20

# ---- Setup -------------------------------------------------------------------

# Working directory
# ... Isidore's window's laptop
baseDir <- "C:/Users/Lenovo/OneDrive/Documents/ihouaga_adgg_spatial"
getwd()

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
dir()

# ---- Installing and load packages --------------------------------------------

if (FALSE) {
  requiredPackages <- c(
    "tidyverse",  # tidyverse for data manipulation and verification for Continuous Ranked Probability Score 
  )
  install.packages(pkgs = requiredPackages)
  install.packages(pkgs = "INLA",
                   repos=c(getOption("repos"),
                           INLA = "https://inla.r-inla-download.org/R/stable"),
                   dep = TRUE)
  inla.upgrade(testing = TRUE)
}
install.packages('verification')

library(tidyverse)
library(INLA)# TODO change to inlabru and add inlabrue estimation above
library(gridExtra)
library(verification)
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

# Adding ward_code as a fixed effect
modelWCF <- as.formula(paste0(modelBase, " + ward_code"))

# Adding ward_codeI as a random IID effect
modelWCRI <- as.formula(paste0(modelBase, " + f(ward_codeI, model = 'iid')"))

# Adding ward_codeI as a random Besag effect
modelWCRB <- as.formula(paste0(modelBase, " + f(ward_codeI, model = 'besag', graph = nb.map, scale.model = TRUE)"))

# ---- Run the models ----------------------------------------------------------

# *fitWCF
sink('results/fitWCF.txt') # Print outputs of fitWCF
'fitWCF'
fitWCF <- inla(formula = modelWCF, data = data1,
               control.compute = list(dic = TRUE,config=TRUE))
# save fitWCF as R object
#save(fitWCF,file = "data/cleaned_data/fitWCF.RData") 
load(file = "data/cleaned_data/fitWCF.RData") # to load the R object 

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
'summary fitWCF'
summary(fitWCF)

'Summarize Variances fitWCF'
summarise_precision_to_variance(fitWCF)
sink()

# fitWCRI

sink('results/fitWCRI.txt') #Print outputs of fitWCRI
'fitWCRI'
fitWCRI <- inla(formula = modelWCRI, data = data1,
               control.compute = list(dic = TRUE, config=TRUE))
#save(fitWCRI,file = "data/cleaned_data/fitWCRI.RData") 
load(file = "data/cleaned_data/fitWCRI.RData") # to load the R object 
'Summary fitWCRI'
summary(fitWCRI)
'sumarize variance fitWCRI'
summarise_precision_to_variance(fitWCRI)
sink() 

# fitWCRB

sink('results/fitWCRB.txt') #Print outputs of fitWCRB
fitWCRB <- inla(formula = modelWCRB, data = data1,
               control.compute = list(dic = TRUE, config=TRUE))
#save(fitWCRB,file = "data/cleaned_data/fitWCRB.RData") 
load(file = "data/cleaned_data/fitWCRB.RData")

'fitWCRB'
'Summary fitWCRB'
summary(fitWCRB)

'Summarize variance fitWCRB'
summarise_precision_to_variance(fitWCRB)
sink() # close sink fitWCRB

#----------Plot posterior distributions of hyperparameters for all models--------

# Genetic variance
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

#TO DO: Plot posterior distribution of all models in a single graph for genetic effect

# Plot genetic variance across models
Pg <- ggplot()  +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCF$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y), color="black") + 
  theme_bw() +
  
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRI$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y), color="blue", linetype = "dashed") +
  theme_bw() +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
      
                                              fitWCRB$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y), color="red", linetype = "dotted") +
  theme_bw() +
  labs(x="", y="", title ="Genetic variance")

Pg <- Pg + scale_y_continuous(breaks=NULL)
  
Pg 
# Residual variance 
# * fitWCF
ggplot(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                fitWCF$marginals.hyperpar$`Precision for the Gaussian observations`))), aes(x, y)) +
  geom_line() +
  theme_bw()

# * fitWCRI
ggplot(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                fitWCRI$marginals.hyperpar$`Precision for the Gaussian observations`))), aes(x, y)) +
  geom_line() +
  theme_bw()

# * fitWCRB
ggplot(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                fitWCRB$marginals.hyperpar$`Precision for the Gaussian observations`))), aes(x, y)) +
  geom_line() +
  theme_bw()
#TO DO: Plot posterior distributions of all models in a single graph for residual variance

# Plotting Residual variance posterior distribution across models

Pr <- ggplot()  +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCF$marginals.hyperpar$`Precision for the Gaussian observations`))), mapping = aes(x, y), color="black") + 
  theme_bw() +
  
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRI$marginals.hyperpar$`Precision for the Gaussian observations`))), mapping = aes(x, y), color="blue", linetype = "dashed") +
  theme_bw() +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     
                                                     fitWCRB$marginals.hyperpar$`Precision for the Gaussian observations`))), mapping = aes(x, y), color="red", linetype = "dotted") +
  theme_bw() +
  labs(x="", y="", title ="Residual variance")
Pr<- Pr + scale_y_continuous(breaks=NULL)
Pr


# Herd_Variance 

# * fitWCF
ggplot(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                fitWCF$marginals.hyperpar$`Precision for herdI`))), aes(x, y)) +
  geom_line() +
  theme_bw()

# * fitWCRI
ggplot(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                fitWCRI$marginals.hyperpar$`Precision for herdI`))), aes(x, y)) +
  geom_line() +
  theme_bw()

# * fitWCRB
ggplot(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                fitWCRB$marginals.hyperpar$`Precision for herdI`))), aes(x, y)) +
  geom_line() +
  theme_bw()

#TO DO: Plot posterior distributions of all models in a single graph for herd effect

Ph <- ggplot()  +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCF$marginals.hyperpar$`Precision for herdI`))), mapping = aes(x, y), color="black") + 
  theme_bw() +
  
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRI$marginals.hyperpar$`Precision for herdI`))), mapping = aes(x, y), color="blue", linetype = "dashed") +
  theme_bw() +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     
                                                     fitWCRB$marginals.hyperpar$`Precision for herdI`))), mapping = aes(x, y), color="red", linetype = "dotted") +
  theme_bw() +
  labs(x="", y="", title ="Herd effect variance")

Ph<- Ph + scale_y_continuous(breaks=NULL)
Ph

# Ward_variance (fitWCRI and fitWCRB)

# * fitWCRI
ggplot(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                fitWCRI$marginals.hyperpar$`Precision for ward_codeI`))), aes(x, y)) +
  geom_line() +
  theme_bw()

# * fitWCRB
ggplot(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                fitWCRB$marginals.hyperpar$`Precision for ward_codeI`))), aes(x, y)) +
  geom_line() +
  theme_bw()
#TO DO: Plot posterior distributions of all models in a single graph for ward effect

Pw <- ggplot()  +
  
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRI$marginals.hyperpar$`Precision for ward_codeI`))), mapping = aes(x, y), color="blue", linetype = "dashed") +
  theme_bw() +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     
                                                     fitWCRB$marginals.hyperpar$`Precision for ward_codeI`))), mapping = aes(x, y), color="red", linetype = "dotted") +
  theme_bw() +
  labs(x="", y="", title ="Ward effect variance")


Pw<- Pw + scale_y_continuous(breaks=NULL)
Pw

#---------Combined Plot Variance components------------------------------------------------
grid.arrange(Pg, Pr, Ph, Pw) # Black= fitWCF, Blue=fitWCRI and Red=fitWCRB

# TO DO: Add proper legend.
 

# -------SPDE-------------------------------------------------------------------
# Make mesh and SPDE 
# Priors
# SPDE
hyperRange  = c(50, 0.8)
hyperVarSpdeS = c(sqrt(0.25), 0.5)
hyperVarSpdeHS = c(sqrt(0.10), 0.5)
hyperResVarGHS = list(theta = list(prior="pc.prec", param=c(sqrt(0.15),0.5)))

mesh = inla.mesh.2d(cbind(data1$long, data1$lat), max.edge=c(10, 20), cutoff = 2.5, offset = 30)
A = inla.spde.make.A(mesh = mesh, loc = cbind(data1$long, data1$lat) )
spdeStatS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeS)
meshIndexS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatS$n.spde) 
spdeStatHS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeHS)
meshIndexHS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatHS$n.spde) 

# Make stack

stack = inla.stack(data = list(milkZ = data1$milkZ),
                   A = list(A,1),
                   effects = list(c(meshIndexS, list(intercept = 1)),
                                  list(cowI = data1$cowI,herdI = data1$herdI, cowPeI = data1$cowPeI,
                                       cyrsn=data1$cyrsn, tyrmn=data1$tyrmn,dgrp= data1$dgrp, ageZ=data1$ageZ, lacgr=data1$lacgr, leg0=data1$leg0, leg1=data1$leg1,leg2=data1$leg2)), tag = "data1.data") 


formulaspde <- as.formula(paste0(modelBase, " + f(fieldID, model = spdeStatS)"))


modelSPDE= inla(formula = formulaspde, data = inla.stack.data(stack),
                family = "normal", control.predictor =list(A=inla.stack.A(stack),compute = T),
                control.family=list(list(hyper=hyperResVarGHS)),
                control.compute = list(dic=T,cpo=F, config=T), 
                control.fixed = list(expand.factor.strategy="inla"), verbose=T)



# I got: "Error in cyrsn + tyrmn : non-conformable arrays" and Error in cyrsn + tyrmn : non-numeric argument to binary operator"

# After making them as.numeric the spde model run # Discuss this with @gg

data1$tyrmn <- as.numeric(data1$tyrmn) 
data1$cyrsn <- as.numeric(data1$cyrsn)
data1$dgrp <- as.numeric(data1$dgrp) 


#-------------------Accuracy of Prediction--------------------------------------
#----------- Accuracy model fitWCF----------------------------------------------

# Masking phenotypes of cows in dgrp 1 (Breed proportion class 1)
#dgrp1 masked
data1_1_NA<- data1 %>% mutate(milkZ = ifelse(dgrp == "1", NA, milkZ))
sum(is.na(data1$milkZ)) # 0 
sum(is.na(data1_1_NA$milkZ)) #7466 
length(unique(data1_1_NA$cowI)) # 1894
length(unique(data1$cowI)) # 1894
fitWCF_NA1 <- inla(formula = modelWCF, data = data1_1_NA,
               control.compute = list(dic = TRUE,config=TRUE))

# save fitWCF_NA1 as R object
save(fitWCF_NA1,file = "data/cleaned_data/fitWCF_NA1.RData") 
 #load(file = "data/cleaned_data/fitWCF_NA1.RData") # to load the R object 

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
accuracy_fitWCF_1 = 0.24
R2_fitWCF_1<- summary(Coef1)
R2_fitWCF_1 = 0.05932
#CRPS
obs <- pheno1$milkZ
pred<- subset(pheno1,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCF_1 <- crps(obs,pred)
round((crps_fitWCF_1$CRPS),2) # 0.6
crps_fitWCF_1= 0.6
# Masking phenotypes of cows in dgrp 2 (Breed proportion class 2)
#dgrp2 masked
data1_2_NA<- data1 %>% mutate(milkZ = ifelse(dgrp == "2", NA, milkZ))
 sum(is.na(data1_2_NA$milkZ)) # 8149
length(unique(data1_2_NA$cowI)) # 1894

fitWCF_NA2 <- inla(formula = modelWCF, data = data1_2_NA,
                   control.compute = list(dic = TRUE,config=TRUE))
save(fitWCF_NA2,file = "data/cleaned_data/fitWCF_NA2.RData") 
#load(file = "data/cleaned_data/fitWCF_NA2.RData") # to load the R object 

pheno_pred2 <- fitWCF_NA2$summary.linear.predictor 
colnames(pheno_pred2)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred2$mean
data1$milkZ_pred_sd <- pheno_pred2$sd

pheno2 <-   subset(data1, dgrp==2) 
pheno2<- subset(pheno2,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))  
length(unique(pheno2$cowI))#  Cows in breed composition class 2

# Correlating observed phenotype of cows in dgrp 2 with predicted phenotypes 
accuracy_fitWCF_2 <- round(cor(pheno2$milkZ,pheno2$milkZ_pred),2) 
Coef2<- lm (pheno2$milkZ~pheno2$milkZ_pred)
accuracy_fitWCF_2 = 0.34
R2_fitWCF_2<- summary(Coef2)
R2_fitWCF_2 = 0.1123
#CRPS
obs <- pheno2$milkZ
pred<- subset(pheno2,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCF_2 <- crps(obs,pred)
round((crps_fitWCF_2$CRPS),2)
crps_fitWCF_2= 0.53

# Masking phenotypes of cows in dgrp 3 (Breed proportion class 3)
#dgrp3 masked
data1_3_NA<- data1 %>% mutate(milkZ = ifelse(dgrp == "3", NA, milkZ))
sum(is.na(data1_3_NA$milkZ)) # 3058 
length(unique(data1_3_NA$cowI)) # 1894

fitWCF_NA3 <- inla(formula = modelWCF, data = data1_3_NA,
                   control.compute = list(dic = TRUE,config=TRUE))

save(fitWCF_NA3,file = "data/cleaned_data/fitWCF_NA3.RData") 
 #load(file = "data/cleaned_data/fitWCF_NA3.RData") # to load the R object 

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
accuracy_fitWCF_3 = 0.28
R2_fitWCF_3<- summary(Coef3)
R2_fitWCF_3= 0.08083
#CRPS
obs <- pheno3$milkZ
pred<- subset(pheno3,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCF_3 <- crps(obs,pred)
round((crps_fitWCF_3$CRPS),2) # 
crps_fitWCF_3= 0.52

# Masking phenotypes of cows in dgrp 4 (Breed proportion class 4)
#dgrp4 masked
data1_4_NA<- data1 %>% mutate(milkZ = ifelse(dgrp == "4", NA, milkZ))
sum(is.na(data1_4_NA$milkZ)) # 702
length(unique(data1_4_NA$cowI)) # 1894

fitWCF_NA4 <- inla(formula = modelWCF, data = data1_4_NA,
                   control.compute = list(dic = TRUE,config=TRUE))
save(fitWCF_NA4,file = "data/cleaned_data/fitWCF_NA4.RData") 
 #load(file = "data/cleaned_data/fitWCF_NA4.RData") # to load the R object 

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
R2_fitWCF_4=0.0008917 
#CRPS
obs <- pheno4$milkZ
pred<- subset(pheno4,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCF_4 <- crps(obs,pred)
round((crps_fitWCF_4$CRPS),2) # 0.58
crps_fitWCF_4= 0.58

# Accuracy of prediction and degree under/overprediction of model fitWCF
  
accuracy_fitWCF = (accuracy_fitWCF_1 + accuracy_fitWCF_2 + accuracy_fitWCF_3 + accuracy_fitWCF_4)/4
round((accuracy_fitWCF),2) # 0.2 

R2_fitWCF = R2_fitWCF_1 + R2_fitWCF_2 + R2_fitWCF_3 + R2_fitWCF_4
round((R2_fitWCF),2) # 

R2_fitWCF=0.25 # TO DO: will edit code to avoid assigning object to numeric value as value may change

crps_fitWCF =(crps_fitWCF_1+crps_fitWCF_2+crps_fitWCF_3+crps_fitWCF_4)/4
round((crps_fitWCF),2) 

crps_fitWCF # 0.56

#----------- Accuracy model fitWCRI (TO Edit as in fitWCF)----------------------------------------------

fitWCRI_NA1 <- inla(formula = modelWCRI, data = data1_1_NA,
                   control.compute = list(dic = TRUE,config=TRUE))

# save fitWCF_NA1 as R object
save(fitWCRI_NA1,file = "data/cleaned_data/fitWCRI_NA1.RData") 
# load(file = "data/cleaned_data/fitWCRI_NA1.RData") # to load the R object 

pheno_pred1 <- fitWCRI_NA1$summary.linear.predictor 
colnames(pheno_pred1)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred1$mean

pheno1 <-   subset(data1, dgrp==1) 
pheno1<- subset(pheno1,select= c(cowI,milkZ,milkZ_pred))  
length(unique(pheno1$cowI))# 731 Cows in breed composition class 1

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotyped 
accuracy_fitWCRI_1 <- round(cor(pheno1$milkZ,pheno1$milkZ_pred),2) 
Coef1<- lm (pheno1$milkZ~pheno1$milkZ_pred)
accuracy_fitWCRI_1 #
R2_fitWCRI_1<- summary(Coef1)

fitWCRI_NA2 <- inla(formula = modelWCRI, data = data1_2_NA,
                   control.compute = list(dic = TRUE,config=TRUE))
save(fitWCRI_NA2,file = "data/cleaned_data/fitWCRI_NA2.RData") 
# load(file = "data/cleaned_data/fitWCRI_NA2.RData") # to load the R object 

pheno_pred2 <- fitWCRI_NA2$summary.linear.predictor 
colnames(pheno_pred2)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred2$mean

pheno2 <-   subset(data1, dgrp==2) 
pheno2<- subset(pheno2,select= c(cowI,milkZ,milkZ_pred))  
length(unique(pheno2$cowI))#  Cows in breed composition class 2

# Correlating observed phenotype of cows in dgrp 2 with predicted phenotypes 
accuracy_fitWCRI_2 <- round(cor(pheno2$milkZ,pheno2$milkZ_pred),2) 
Coef2<- lm (pheno2$milkZ~pheno2$milkZ_pred)
accuracy_fitWCRI_2 #
R2_fitWCRI_2<- summary(Coef2)


#dgrp3 masked

fitWCRI_NA3 <- inla(formula = modelWCRI, data = data1_3_NA,
                   control.compute = list(dic = TRUE,config=TRUE))

save(fitWCRI_NA3,file = "data/cleaned_data/fitWCRI_NA3.RData") 
# load(file = "data/cleaned_data/fitWCRI_NA3.RData") # to load the R object 

pheno_pred3 <- fitWCRI_NA3$summary.linear.predictor 
colnames(pheno_pred3)
sum(is.na(pheno_pred3$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred3$mean

pheno3 <-   subset(data1, dgrp==3) 
pheno3<- subset(pheno3,select= c(cowI,milkZ,milkZ_pred))  
length(unique(pheno3$cowI))#  Cows in breed composition class 3

# Correlating observed phenotype of cows in dgrp 3 with predicted phenotypes 
accuracy_fitWCRI_3 <- round(cor(pheno3$milkZ,pheno3$milkZ_pred),2) 
Coef3<- lm (pheno3$milkZ~pheno3$milkZ_pred)
accuracy_fitWCRI_3 #
R2_fitWCRI_3 <-summary(Coef3)

fitWCRI_NA4 <- inla(formula = modelWCRI, data = data1_4_NA,
                   control.compute = list(dic = TRUE,config=TRUE))
save(fitWCRI_NA4,file = "data/cleaned_data/fitWCRI_NA4.RData") 
# load(file = "data/cleaned_data/fitWCRI_NA4.RData") # to load the R object

pheno_pred4 <- fitWCRI_NA4$summary.linear.predictor 
colnames(pheno_pred4)
sum(is.na(pheno_pred4$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred4$mean

pheno4 <-   subset(data1, dgrp==4) 
pheno4<- subset(pheno4,select= c(cowI,milkZ,milkZ_pred))  
length(unique(pheno4$cowI))#  Cows in breed composition class 4

# Correlating observed phenotype of cows in dgrp 4 with predicted phenotypes 
accuracy_fitWCRI_4 <- round(cor(pheno4$milkZ,pheno4$milkZ_pred),2) 
Coef4<- lm (pheno4$milkZ~pheno4$milkZ_pred)
accuracy_fitWCRI_4 #
R2_fitWCRI_4 <- summary(Coef4)

# Accuracy of prediction and degree under/overprediction of model fitWCF

accuracy_fitWCRI = (accuracy_fitWCRI_1 + accuracy_fitWCRI_2 + accuracy_fitWCRI_3 + accuracy_fitWCRI_4)/4
accuracy_fitWCRI # 

R2_fitWCRI = R2_fitWCRI_1 + R2_fitWCRI_2 + R2_fitWCRI_3 + R2_fitWCRI_4
R2_fitWCRI # 


#----------- Accuracy model fitWCRB---------------------------------------------

fitWCRB_NA1 <- inla(formula = modelWCRB, data = data1_1_NA,
                    control.compute = list(dic = TRUE,config=TRUE))

# save fitWCRB_NA1 as R object
save(fitWCRB_NA1,file = "data/cleaned_data/fitWCRB_NA1.RData") 
# load(file = "data/cleaned_data/fitWCRB_NA1.RData") # to load the R object 

pheno_pred1 <- fitWCRB_NA1$summary.linear.predictor 
colnames(pheno_pred1)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred1$mean

pheno1 <-   subset(data1, dgrp==1) 
pheno1<- subset(pheno1,select= c(cowI,milkZ,milkZ_pred))  
length(unique(pheno1$cowI))# 731 Cows in breed composition class 1

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotyped 
accuracy_fitWCRB_1 <- round(cor(pheno1$milkZ,pheno1$milkZ_pred),2) 
Coef1<- lm (pheno1$milkZ~pheno1$milkZ_pred)
accuracy_fitWCRB_1 #
R2_fitWCRB_1<- summary(Coef1)

fitWCRB_NA2 <- inla(formula = modelWCRB, data = data1_2_NA,
                    control.compute = list(dic = TRUE,config=TRUE))
save(fitWCRB_NA2,file = "data/cleaned_data/fitWCRB_NA2.RData") 
# load(file = "data/cleaned_data/fitWCRB_NA2.RData") # to load the R object 

pheno_pred2 <- fitWCRB_NA2$summary.linear.predictor 
colnames(pheno_pred2)
sum(is.na(pheno_pred1$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred2$mean

pheno2 <-   subset(data1, dgrp==2) 
pheno2<- subset(pheno2,select= c(cowI,milkZ,milkZ_pred))  
length(unique(pheno2$cowI))#  Cows in breed composition class 2

# Correlating observed phenotype of cows in dgrp 2 with predicted phenotypes 
accuracy_fitWCRB_2 <- round(cor(pheno2$milkZ,pheno2$milkZ_pred),2) 
Coef2<- lm (pheno2$milkZ~pheno2$milkZ_pred)
accuracy_fitWCRB_2 #
R2_fitWCRB_2<- summary(Coef2)


#dgrp3 masked

fitWCRB_NA3 <- inla(formula = modelWCRB, data = data1_3_NA,
                    control.compute = list(dic = TRUE,config=TRUE))

save(fitWCRB_NA3,file = "data/cleaned_data/fitWCRB_NA3.RData") 
# load(file = "data/cleaned_data/fitWCRB_NA3.RData") # to load the R object 

pheno_pred3 <- fitWCRB_NA3$summary.linear.predictor 
colnames(pheno_pred3)
sum(is.na(pheno_pred3$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred3$mean

pheno3 <-   subset(data1, dgrp==3) 
pheno3<- subset(pheno3,select= c(cowI,milkZ,milkZ_pred))  
length(unique(pheno3$cowI))#  Cows in breed composition class 3

# Correlating observed phenotype of cows in dgrp 3 with predicted phenotypes 
accuracy_fitWCRB_3 <- round(cor(pheno3$milkZ,pheno3$milkZ_pred),2) 
Coef3<- lm (pheno3$milkZ~pheno3$milkZ_pred)
accuracy_fitWCRB_3 #
R2_fitWCRB_3 <-summary(Coef3)

fitWCRB_NA4 <- inla(formula = modelWCRB, data = data1_4_NA,
                    control.compute = list(dic = TRUE,config=TRUE))
save(fitWCRB_NA4,file = "data/cleaned_data/fitWCRB_NA4.RData") 
# load(file = "data/cleaned_data/fitWCRB_NA4.RData") # to load the R object

pheno_pred4 <- fitWCRB_NA4$summary.linear.predictor 
colnames(pheno_pred4)
sum(is.na(pheno_pred4$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred4$mean

pheno4 <-   subset(data1, dgrp==4) 
pheno4<- subset(pheno4,select= c(cowI,milkZ,milkZ_pred))  
length(unique(pheno4$cowI))#  Cows in breed composition class 4

# Correlating observed phenotype of cows in dgrp 4 with predicted phenotypes 
accuracy_fitWCRB_4 <- round(cor(pheno4$milkZ,pheno4$milkZ_pred),2) 
Coef4<- lm (pheno4$milkZ~pheno4$milkZ_pred)
accuracy_fitWCRB_4 #
R2_fitWCRB_4 <- summary(Coef4)

# Accuracy of prediction and degree of under/overprediction of model fitWCRB

accuracy_fitWCRB = (accuracy_fitWCRB_1 + accuracy_fitWCRB_2 + accuracy_fitWCRB_3 + accuracy_fitWCRB_4)/4
accuracy_fitWCRB # 

R2_fitWCRB = (R2_fitWCRB_1 + R2_fitWCRB_2 + R2_fitWCRB_3 + R2_fitWCRB_4)/4
R2_fitWCRB # 











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

# Explore the the model outputs







