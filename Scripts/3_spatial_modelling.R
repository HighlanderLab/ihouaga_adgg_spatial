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
    "tidyverse", "verification", "MASS", # tidyverse for data manipulation and verification for Continuous Ranked Probability Score 
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
library(gridExtra) # Visualize random field grid
library(verification) # Visualize random field grid
library(lattice)
library(raster) # GetData to plot country map
library(rgeos)
library(viridis) # Plot mean and sd spatial effect
library(fields)# Plot mean and sd spatial effect
library(ggpubr) # Plot mean and sd spatial effect
library("MASS") # Pca
library("factoextra") # Pca
library(Hmisc) # Correlation matrix with P-values

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

# Base model for milk - All effects without ward_code/ward_codeI 
modelBase <- "milkZ ~ 1 + cyrsnI + tyrmnI + dgrpI + (ageZ|lacgr) +
                          (leg0|lacgr) + (leg1|lacgr) + (leg2|lacgr) +
                          f(herdI, model = 'iid') +
                          f(cowPeI, model = 'iid') +
                          f(cowI, model = 'generic0', Cmatrix = GRMInv)"

# Adding ward_code as a fixed effect
modelWCF <- as.formula(paste0(modelBase, " + ward_code"))

# Adding ward_codeI as a random IID effect
modelWCRI <- as.formula(paste0(modelBase, " + f(ward_codeI, model = 'iid')"))

# Adding ward_codeI as a random Besag effect
modelWCRB <- as.formula(paste0(modelBase, " + f(ward_codeI, model = 'besag', graph = nb.map, scale.model = TRUE)"))

modelBase <- as.formula(modelBase)

# ---- Specifying additional models----------------------------------------------------------
#Base model without herd effect and without word effect

# Base model for milk - All effects without ward_code/ward_codeI 
modelBaseNoherd <- " milkZ ~ 1 + cyrsnI + tyrmnI + dgrpI + (ageZ|lacgr) +
                             (leg0|lacgr) + (leg1|lacgr) + (leg2|lacgr) +
                             f(cowPeI, model = 'iid') +
                             f(cowI, model = 'generic0', Cmatrix = GRMInv)"

# Adding ward_code as a fixed effect
modelWCFNoherd <- as.formula(paste0(modelBaseNoherd, " + ward_code"))

# Adding ward_codeI as a random IID effect
modelWCRINoherd <- as.formula(paste0(modelBaseNoherd, " + f(ward_codeI, model = 'iid')"))

# Adding ward_codeI as a random Besag effect
modelWCRBNoherd <- as.formula(paste0(modelBaseNoherd, " + f(ward_codeI, model = 'besag', graph = nb.map, scale.model = TRUE)"))

modelBase <- as.formula(modelBaseNoherd)

# Base model + Region (fixed)
# modelBaseRegion <- as.formula(paste0(modelBase, " + regionI"))
# Adding region fixed to WCRI
# modelWCRIRegion <- as.formula(paste0(modelWCRI, " + regionI"))

# ---- Run the models ----------------------------------------------------------
# *fitBase

sink('results/fitBase.txt') # Print outputs of fitBase
cat('fitBase\n')
fitBase <- inla(formula = modelBase, data = data1,
               control.compute = list(dic = TRUE, config=TRUE))

save(fitBase,file = "data/cleaned_data/fitBase/fitBase.RData") 
load(file = "data/cleaned_data/fitBase/fitBase.RData") # to load the R object 

# *fitWCF
sink('results/fitWCF.txt') # Print outputs of fitWCF
cat('fitWCF\n')
fitWCF <- inla(formula = modelWCF, data = data1,
               control.compute = list(dic = TRUE, config=TRUE))
# save fitWCF as R object
save(fitWCF,file = "data/cleaned_data/fitWCF/fitWCF.RData") 
load(file = "data/cleaned_data/fitWCF/fitWCF.RData") # to load the R object 

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
save(fitWCRI,file = "data/cleaned_data/fitWCRI/fitWCRI.RData") 
load(file = "data/cleaned_data/fitWCRI/fitWCRI.RData") # to load the R object 
'Summary fitWCRI'
summary(fitWCRI)
'sumarize variance fitWCRI'
summarise_precision_to_variance(fitWCRI)
sink() 
summary(fitWCRI)

# fitWCRB

sink('results/fitWCRB.txt') #Print outputs of fitWCRB
fitWCRB <- inla(formula = modelWCRB, data = data1,
               control.compute = list(dic = TRUE, config=TRUE))
save(fitWCRB,file = "data/cleaned_data/fitWCRB/fitWCRB.RData") 
load(file = file = "data/cleaned_data/fitWCRB/fitWCRB.RData")

'fitWCRB'
'Summary fitWCRB'
summary(fitWCRB)

'Summarize variance fitWCRB'
summarise_precision_to_variance(fitWCRB)
sink() # close sink fitWCRB

# -------SPDE-------------------------------------------------------------------
# Make mesh and SPDE 
# Priors
# SPDE
hyperRange  = c(50, 0.8)
hyperVarSpdeS = c(sqrt(0.25), 0.5)
hyperVarSpdeWS = c(sqrt(0.10), 0.5)
hyperResVarGWS = list(theta = list(prior="pc.prec", param=c(sqrt(0.15),0.5)))

mapTZA<-getData('GADM', country="TZA", level= 1) #Get map of Tanzania
# Extract the border from the map
hr.border <- gUnaryUnion(mapTZA, rep(1, nrow(mapTZA)))
# Formatting boundary for r-INLA
hr.bdry <- inla.sp2segment(hr.border)
locations = cbind(data1$long, data1$lat)
mesh = inla.mesh.2d(loc=locations, boundary  = hr.bdry, max.edge = c(0.4,0.8), cutoff = 0.08, offset = c(0.7, 0.7))
plot(mesh, main="3rd mesh") ; points(data1$long, data1$lat, col="red") # Plotting mesh

A = inla.spde.make.A(mesh = mesh, loc = cbind(data1$long, data1$lat) )
spdeStatS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeS)
meshIndexS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatS$n.spde) 
spdeStatWS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeWS)
meshIndexWS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatWS$n.spde) 

# Make stack
# StackS
stackS = inla.stack(data = list(milkZ = data1$milkZ),
                    A = list(A,1),
                    effects = list(c(meshIndexS, list(intercept = 1)),
                                   list(cowI = data1$cowI,herdI = data1$herdI, cowPeI = data1$cowPeI,
                                        cyrsnI=data1$cyrsnI, tyrmnI=data1$tyrmnI,dgrpI= data1$dgrpI, ageZ=data1$ageZ, lacgr=data1$lacgr, leg0=data1$leg0, leg1=data1$leg1,leg2=data1$leg2)), tag = "data1.data") 



# StackWS (Ward as random + Spatial effect)
stackWS = inla.stack(data = list(milkZ = data1$milkZ),
                     A = list(A,1),
                     effects = list(c(meshIndexWS, list(intercept = 1)),
                                    list(cowI = data1$cowI,herdI = data1$herdI, cowPeI = data1$cowPeI,
                                         cyrsnI=data1$cyrsnI, tyrmnI=data1$tyrmnI,ward_codeI= data1$ward_codeI, dgrpI= data1$dgrpI, ageZ=data1$ageZ, lacgr=data1$lacgr, leg0=data1$leg0, leg1=data1$leg1,leg2=data1$leg2)), tag = "data1.data") 


# ModelS
formulaS <- as.formula(paste0(modelBase, " + f(fieldID, model = spdeStatS)-1"))


fitS= inla(formula = formulaS, data = inla.stack.data(stackS),
           family = "normal", control.predictor =list(A=inla.stack.A(stackS),compute = T),
           control.family=list(list(hyper=hyperResVarGWS)),
           control.compute = list(dic=T,cpo=F, config=T), 
           control.fixed = list(expand.factor.strategy="inla"), verbose=T)
save(fitS,file = "D:/Results_ADGG_Spatial/fitS/fitS.RData")


# fitWS

modelWCRI <- "milkZ ~ 1 + cyrsnI + tyrmnI + dgrpI + (ageZ|lacgr) + (leg0|lacgr) + (leg1|lacgr) + (leg2|lacgr) + f(herdI, model = 'iid') + f(cowPeI, model = 'iid') + f(ward_codeI, model = 'iid') +  f(cowI, model = 'generic0', Cmatrix = GRMInv)"

formulaWS <- as.formula(paste0(modelWCRI, " + f(fieldID, model = spdeStatWS)-1"))

fitWS= inla(formula = formulaWS, data = inla.stack.data(stackWS),
            family = "normal", control.predictor =list(A=inla.stack.A(stackWS),compute = T),
            control.family=list(list(hyper=hyperResVarGWS)),
            control.compute = list(dic=T,cpo=F, config=T), 
            control.fixed = list(expand.factor.strategy="inla"), verbose=T)

save(fitWS,file = "D:/Results_ADGG_Spatial/fitWS/fitWS.RData") 
# loading from different directory
getwd()
summary(fitBase)
summary(fitWCF)
summary(fitWCRI)
summary(fitWCRB)
summary(fitS)
summary(fitWS)
SummarizeInlaSpdeVars(fitS) #Error "object 'SpdeStat' not found"
SummarizeInlaSpdeVars(fitWS) #Error "object 'SpdeStat' not found"

#---------Run additional models--------------------------------------------------
#FitBaseNoherd (B no herd)
fitBaseNoherd <- inla(formula = modelBaseNoherd, data = data1,
                control.compute = list(dic = TRUE,config=TRUE))
summary(fitBaseNoherd) # DIC=35546.55

#fitWCFNoherd 
# *fitWCFNoherd 

fitWCFNoherd  <- inla(formula = modelWCFNoherd, data = data1,
               control.compute = list(dic = TRUE,config=TRUE))


#FitBaseNoherdWCRI (WCRI no herd)

fitWCRINoherd  <- inla(formula = modelWCRINoherd, data = data1,
                           control.compute = list(dic = TRUE,config=TRUE)) 
summary(fitWCRINoherd)




fitWCRBNoherd <- inla(formula = modelWCRBNoherd, data = data1,
                control.compute = list(dic = TRUE, config=TRUE))
summary(fitWCRBNoherd) 

# Base+Region (fixed) B+Region (BR)
#fitBaseregion <- inla(formula = modelBaseRegion, data = data1,
                      #control.compute = list(dic = TRUE,config=TRUE)) 

#summary(fitBaseregion)

# WCRI +Region (fixed) WCRI + Region (Fixed)
#fitWCRIregion <- inla(formula = modelWCRIRegion, data = data1,
                   #control.compute = list(dic = TRUE,config=TRUE)) 
#table(data1$regionI)

#summary(fitWCRIregion)

# -------SPDE_No Herd-------------------------------------------------------------------
# Make mesh and SPDE 
# Priors
# SPDE
hyperRange  = c(50, 0.8)
hyperVarSpdeS = c(sqrt(0.25), 0.5)
hyperVarSpdeWS = c(sqrt(0.10), 0.5)
hyperResVarGWS = list(theta = list(prior="pc.prec", param=c(sqrt(0.15),0.5)))

mapTZA<-getData('GADM', country="TZA", level= 1) #Get map of Tanzania
# Extract the border from the map
hr.border <- gUnaryUnion(mapTZA, rep(1, nrow(mapTZA)))
# Formatting boundary for r-INLA
hr.bdry <- inla.sp2segment(hr.border)
locations = cbind(data1$long, data1$lat)
mesh = inla.mesh.2d(loc=locations, boundary  = hr.bdry, max.edge = c(0.4,0.8), cutoff = 0.08, offset = c(0.7, 0.7))
plot(mesh, main="3rd mesh") ; points(data1$long, data1$lat, col="red") # Plotting mesh

A = inla.spde.make.A(mesh = mesh, loc = cbind(data1$long, data1$lat) )
spdeStatS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeS)
meshIndexS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatS$n.spde) 
spdeStatWS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeWS)
meshIndexWS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatWS$n.spde) 

# Make stack
# StackSNoherd
stackSNoherd = inla.stack(data = list(milkZ = data1$milkZ),
                    A = list(A,1),
                    effects = list(c(meshIndexS, list(intercept = 1)),
                                   list(cowI = data1$cowI, cowPeI = data1$cowPeI,
                                        cyrsnI=data1$cyrsnI, tyrmnI=data1$tyrmnI,dgrpI= data1$dgrpI, ageZ=data1$ageZ, lacgr=data1$lacgr, leg0=data1$leg0, leg1=data1$leg1,leg2=data1$leg2)), tag = "data1Noherd.data") 



# StackWS (Ward as random + Spatial effect)
stackWSNoherd = inla.stack(data = list(milkZ = data1$milkZ),
                     A = list(A,1),
                     effects = list(c(meshIndexWS, list(intercept = 1)),
                                    list(cowI = data1$cowI, cowPeI = data1$cowPeI,
                                         cyrsnI=data1$cyrsnI, tyrmnI=data1$tyrmnI,ward_codeI= data1$ward_codeI, dgrpI= data1$dgrpI, ageZ=data1$ageZ, lacgr=data1$lacgr, leg0=data1$leg0, leg1=data1$leg1,leg2=data1$leg2)), tag = "data1Noherd.data") 


# ModelS
formulaSNoherd <- as.formula(paste0(modelBaseNoherd, " + f(fieldID, model = spdeStatS)-1"))


fitSNoherd= inla(formula = formulaSNoherd, data = inla.stack.data(stackSNoherd),
           family = "normal", control.predictor =list(A=inla.stack.A(stackSNoherd),compute = T),
           control.family=list(list(hyper=hyperResVarGWS)),
           control.compute = list(dic=T,cpo=F, config=T), 
           control.fixed = list(expand.factor.strategy="inla"), verbose=T)
save(fitS,file = "D:/Results_ADGG_Spatial/fitS/fitS.RData")


# fitWS

modelWCRINoherd <- "milkZ ~ 1 + cyrsnI + tyrmnI + dgrpI + (ageZ|lacgr) + (leg0|lacgr) + (leg1|lacgr) + (leg2|lacgr) + f(cowPeI, model = 'iid') + f(ward_codeI, model = 'iid') +  f(cowI, model = 'generic0', Cmatrix = GRMInv)"

formulaWSNoherd <- as.formula(paste0(modelWCRINoherd, " + f(fieldID, model = spdeStatWS)-1"))

fitWSNoherd= inla(formula = formulaWSNoherd, data = inla.stack.data(stackWSNoherd),
            family = "normal", control.predictor =list(A=inla.stack.A(stackWSNoherd),compute = T),
            control.family=list(list(hyper=hyperResVarGWS)),
            control.compute = list(dic=T,cpo=F, config=T), 
            control.fixed = list(expand.factor.strategy="inla"), verbose=T)

save(fitWS,file = "D:/Results_ADGG_Spatial/fitWS/fitWS.RData") 
# loading from different directory
getwd()
summary(fitBaseNoherd)
summary(fitWCFNoherd)
summary(fitWCRINoherd)
summary(fitWCRBNoherd)
summary(fitSNoherd)
summary(fitWSNoherd)

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

#-------------------EAAP_Plot---------------------------------------------------
# Plot genetic variance across models

PgE<-ggplot()  +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitBase$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="B")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWCRI$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="BW")) +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitS$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="BS"), linetype = "dashed") +
  geom_line(data.frame(inla.smarginal(inla.tmarginal(function(x) 1/x,
                                                     fitWS$marginals.hyperpar$`Precision for cowI`))), mapping = aes(x, y, colour="BWS"), linetype = "dashed") +
  scale_color_manual(values = c("B" = "black", 
                                "BW"="red",
                                "BS"="blue",
                                "BWS"="green")) +
  labs(x="", y="", title ="Genetic variance") + scale_y_continuous(breaks=NULL) + theme_bw()  +
  theme(legend.position = "left",
        legend.box.background = element_blank(), 
        legend.title = element_text("Models")) +
  guides(colour=guide_legend(title="Model"))
# scale_y_continuous(breaks=NULL) to remove axis labels
PgE<- PgE+ theme(plot.title = element_text(face = "bold"))


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

# Permanent environmental effect

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
PaperThemeNoLegend = PaperTheme + theme(legend.position = "none")
PaperSize = 10
PreseTheme = theme_bw(base_size = 18, base_family = "sans")  
  theme(legend.position = "top")
PreseThemeNoLegend = PreseTheme + theme(legend.position = "none")
PreseSize = 16

(pE <- ggarrange(PgE,PrE,PhE,PwE,SvE,SrE, ncol = 2, nrow = 3, common.legend = T, legend = "left", align = "h", widths = c(1,1,1))) + PreseTheme
pE

ggsave(plot = pE + PreseTheme, filename = "Posterior_distribution_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = pE + PaperTheme, filename = "Posterior_distribution_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper



# ---- Plot posterior mean (a) and standard deviation of the estimated spatial effect ----
#fitWS

# spatial_est=fitWS$summary.random$fieldID[,1:3]
#fitWS
gproj <- inla.mesh.projector(mesh,  dims = c(300, 300))
g.mean <- inla.mesh.project(gproj, fitWS$summary.random$fieldID$mean)
g.sd <- inla.mesh.project(gproj, fitWS$summary.random$fieldID$sd)

grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='Easting', ylab='Northing', main='mean',col.regions = viridis(16)),
             levelplot(g.sd, scal=list(draw=F), xla='Easting', yla='Northing', main='sd' ,col.regions = viridis(16)), nrow=2)


#FitS
gproj <- inla.mesh.projector(mesh,  dims = c(300, 300))
g.mean <- inla.mesh.project(gproj, fitS$summary.random$fieldID$mean)
g.sd <- inla.mesh.project(gproj, fitS$summary.random$fieldID$sd)

grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='Easting', ylab='Northing', main='mean',col.regions = viridis(16)),
             levelplot(g.sd, scal=list(draw=F), xla='Easting', yla='Northing', main='sd' ,col.regions = viridis(16)), nrow=2)


#------------Boxplot------------------------------------------------------------

# Correlation between EBVs in models with spatial effects and in models without spatial effect

EBV_BW <- data.frame(fitWCRI$summary.random$cowI)
EBV_BWS <- data.frame(fitWS$summary.random$cowI)
EBV_B <- data.frame(fitBase$summary.random$cowI)
EBV_S <- data.frame(fitS$summary.random$cowI)
# Remove EBVs of non-phenotyped animals
sel <- EBV_BW$ID %in% data1$cowI
EBV_BW <- EBV_BW[sel, ]
dim(EBV_BW)

sel <- EBV_BWS$ID %in% data1$cowI
EBV_BWS <- EBV_BWS[sel, ]
dim(EBV_BWS)

sel <- EBV_B$ID %in% data1$cowI
EBV_B <- EBV_B[sel, ]
dim(EBV_B)

sel <- EBV_S$ID %in% data1$cowI
EBV_S <- EBV_S[sel, ]
dim(EBV_S)

# Correlation between EBV_BW and EBV_BWS
round(cor(EBV_BW$mean,EBV_BWS$mean), 2) # 0.79

# Correlation between EBV_B and EBV_BS
round(cor(EBV_B$mean,EBV_S$mean),2) # 0.67

spde.result = inla.spde.result(fitS, "fieldID", spdeStatS, do.transform = TRUE)

fitWS$marginals.linear.predictor$fieldID

spatial<- data.frame(fitWS$marginals.random)

# Get index for random spatial field (fitS and fitWS)

#fitS
index_fitS <- inla.stack.index(stackS,'data1.data')$data

spatial_effect_S <- data.frame(fitS$summary.linear.pred[index_fitS,1])
colnames(spatial_effect_S)[1] <- "spdeS"
data1$spdeS<- spatial_effect_S$spdeS
#fitWS
index_fitWS <- inla.stack.index(stackWS,'data1.data')$data

spatial_effect_WS <- data.frame(fitWS$summary.linear.pred[index_fitWS,1])
colnames(spatial_effect_WS)[1] <- "spdeWS"

data1$spdeWS<- spatial_effect_WS$spdeWS

spatial_effect<- subset(data1, select = c(cowI, spdeS, spdeWS))

SPDE_cow <- spatial_effect %>% group_by(cowI) %>% 
                            summarise(spdeS_mean= mean(spdeS), spdeWS_mean= mean(spdeWS)) 

# Correlating spatial effect with EBV_BW

round(cor(EBV_BW$mean,SPDE_cow$spdeWS_mean),2) # 0.63

round(cor(EBV_BWS$mean,SPDE_cow$spdeWS_mean),2) # 0.51

round(cor(EBV_B$mean,SPDE_cow$spdeS_mean),2) # 0.66

colnames(EBV_B)[2]<- "EBV_B"
colnames(EBV_S)[2]<- "EBV_S"

colnames(EBV_BW)[2]<- "EBV_BW"
colnames(EBV_BWS)[2]<- "EBV_BWS"
colnames(SPDE_cow)[1] <- "ID"

EBV_B<- EBV_B[,1:2]
EBV_S<- EBV_S[,1:2]
EBV_BW<- EBV_BW[,1:2]
EBV_BWS<- EBV_BWS[,1:2]

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
e<- e + labs(x = "Spatial effect", y = "Difference in breeding value (B-BS)")
e<- e+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=20,face="bold")) 


  
ggsave(plot = e + PreseTheme, filename = "Boxplot_difference_EBVs_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm")


ggsave(plot = e + PaperTheme, filename = "Boxplot_difference_EBVs_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm")


# Let's Express difference in EBv as sd
#fitS
index_fitS <- inla.stack.index(stackS,'data1.data')$data

spatial_effect_sd_S <- data.frame(fitS$summary.linear.pred[index_fitS,2])
colnames(spatial_effect_sd_S)[1] <- "sd_spdeS"
data1$sd_spdeS<- spatial_effect_sd_S$sd_spdeS
#fitWS
index_fitWS <- inla.stack.index(stackWS,'data1.data')$data

spatial_effect_sd_WS <- data.frame(fitWS$summary.linear.pred[index_fitWS,2])
colnames(spatial_effect_sd_WS)[1] <- "sd_spdeWS"

data1$sd_spdeWS<- spatial_effect_sd_WS$sd_spdeWS

spatial_effect<- subset(data1, select = c(cowI, sd_spdeS, sd_spdeWS))

sd_SPDE_cow <- spatial_effect %>% group_by(cowI) %>% 
  summarise(sd_spdeS_mean= mean(sd_spdeS), sd_spdeWS_mean= mean(sd_spdeWS)) 


colnames(EBV_B)[3]<- "sd_EBV_B"
colnames(EBV_S)[3]<- "sd_EBV_S"

colnames(EBV_BW)[3]<- "sd_EBV_BW"
colnames(EBV_BWS)[3]<- "sd_EBV_BWS"
colnames(SPDE_cow)[1] <- "ID"

sd_EBV_B<- EBV_B[,1:3]
sd_EBV_S<- EBV_S[,1:3]
sd_EBV_BW<- EBV_BW[,1:3]
sd_EBV_BWS<- EBV_BWS[,1:3]

sd_EBV_Spde <- cbind(sd_EBV_B,sd_EBV_S, sd_EBV_BW, sd_EBV_BWS, sd_SPDE_cow)

sd_EBV_Spde<- subset(sd_EBV_Spde, select = c(ID, sd_EBV_B,sd_EBV_S, sd_EBV_BW, sd_EBV_BWS, sd_spdeS_mean, sd_spdeWS_mean) )
sd_EBV_Spde$dBS <- sd_EBV_Spde$sd_EBV_B-sd_EBV_Spde$sd_EBV_S
round(summary(sd_EBV_Spde$sd_spdeS_mean),2)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.09    0.15    0.18    0.19    0.22    0.29 
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
length(unique(cowbdate$cow))

# Create birthyear column in data1 by merging with cowbdate

data1 <-  merge(data1, cowbdate, by= "cow", all.x = TRUE)
summary(data1$birthyear)
#2004  2005  2006  2007  2008  2009  2010  2011  2012  2013  2014  2015  2016  2017 
#17    12     7    38    41     7   122   232  3044  1066 12017  1298  1444    30 
data1$birthyear <- as.factor(data1$birthyear)
summary(data1$birthyear)
#data1_2016_2017 <- subset(data1_birthdate, birthyear==2016 | birthyear==2017 ) 

#-----------------------Foward validation Fitbase-------------------------------
#Let's make NA milkz of cows born in 2016 and 2017
data1_NA<- data1 %>% mutate(milkZ = ifelse(birthyear=="2016" | birthyear=="2017", NA, milkZ))
sum(is.na(data1_NA$milkZ)) # 1474 Expected
length(unique(data1_NA$cowI)) # 1894 expected
fitBase_NA <- inla(formula = modelfitBase, data = data1_NA,
                    control.compute = list(dic = TRUE,config=TRUE))

#Saving forward validation prediction fitBase External Drive
save(fitBase_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitBase_NA.RData") 
pheno_pred <- fitBase_NA$summary.linear.predictor 
colnames(pheno_pred)
sum(is.na(pheno_pred$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred$mean
data1$milkZ_pred_sd <- pheno_pred$sd 
pheno <-   subset(data1, birthyear=="2016" | birthyear=="2017") 
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))  
length(unique(pheno$cowI))# 146 Cows expected to be born between 2016 (141 cows) and 2017 (5cows)

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotype 
accuracy_fitBase <- round(cor(pheno$milkZ,pheno$milkZ_pred),2) 
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitBase #  0.4
R2_fitBase<- summary(Coef)
R2_fitBase #  0.1582 
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitBase <- crps(obs,pred)
round((crps_fitBase$CRPS),2) # 0.56

#-----------------------Foward validation FitWCF-------------------------------
fitWCF_NA <- inla(formula = modelWCF, data = data1_NA,
                   control.compute = list(dic = TRUE,config=TRUE))


#Saving forward validation prediction fitWCF External Drive
save(fitWCF_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitWCF_NA.RData") 
pheno_pred <- fitWCF_NA$summary.linear.predictor 
colnames(pheno_pred)
sum(is.na(pheno_pred$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred$mean
data1$milkZ_pred_sd <- pheno_pred$sd 
pheno <-   subset(data1, birthyear=="2016" | birthyear=="2017") 
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))  
length(unique(pheno$cowI))# 146 Cows expected to be born between 2016 (141 cows) and 2017 (5cows)

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotype 
accuracy_fitWCF <- round(cor(pheno$milkZ,pheno$milkZ_pred),2) 
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitWCF #  0.4
R2_fitWCF<- summary(Coef)
R2_fitWCF #  0.1582 
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCF <- crps(obs,pred)
round((crps_fitWCF$CRPS),2) #0.56

#-----------------------Foward validation fitWCRI-------------------------------
fitWCRI_NA <- inla(formula = modelWCRI, data = data1_NA,
                   control.compute = list(dic = TRUE,config=TRUE))


#Saving forward validation prediction fitWCRI External Drive
#save(fitWCRI_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitWCRI_NA.RData") 
pheno_pred <- fitWCRI_NA$summary.linear.predictor 
colnames(pheno_pred)
sum(is.na(pheno_pred$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred$mean
data1$milkZ_pred_sd <- pheno_pred$sd 
pheno <-   subset(data1, birthyear=="2016" | birthyear=="2017") 
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))  
length(unique(pheno$cowI))# 146 Cows expected to be born between 2016 (141 cows) and 2017 (5cows)

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotype 
accuracy_fitWCRI <- round(cor(pheno$milkZ,pheno$milkZ_pred),2) 
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitWCRI #  0.59
R2_fitWCRI<- summary(Coef)
R2_fitWCRI #  0.3529 
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCRI <- crps(obs,pred)
round((crps_fitWCRI$CRPS),2) #0.47

#-----------------------Foward validation fitWCRB-------------------------------
fitWCRB_NA <- inla(formula = modelWCRB, data = data1_NA,
                  control.compute = list(dic = TRUE,config=TRUE))


#Saving forward validation prediction fitWCRB External Drive
save(fitWCRB_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitWCRB_NA.RData") 
pheno_pred <- fitWCRB_NA$summary.linear.predictor 
colnames(pheno_pred)
sum(is.na(pheno_pred$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred$mean
data1$milkZ_pred_sd <- pheno_pred$sd 
pheno <-   subset(data1, birthyear=="2016" | birthyear=="2017") 
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))  
length(unique(pheno$cowI))# 146 Cows expected to be born between 2016 (141 cows) and 2017 (5cows)

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotype 
accuracy_fitWCRB <- round(cor(pheno$milkZ,pheno$milkZ_pred),2) 
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitWCRB #  0.59
R2_fitWCRB<- summary(Coef)
R2_fitWCRB #  0.3504
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitWCRB <- crps(obs,pred)
round((crps_fitWCRB$CRPS),2) # 0.47 

#-----------------------Foward validation model fitS-------------------------------

A.pred = inla.spde.make.A(mesh = mesh, loc = cbind(data1$long, data1$lat) )

# Make stack_pred
# StackS_NA
stackS_pred_NA = inla.stack(data = list(milkZ = NA),
                             A = list(A.pred,1),
                             effects = list(c(meshIndexS, list(intercept = 1)),
                                            list(cowI = data1_NA$cowI,herdI = data1_NA$herdI, cowPeI = data1_NA$cowPeI,
                                                 cyrsnI=data1_NA$cyrsnI, tyrmnI=data1_NA$tyrmnI,dgrpI= data1_NA$dgrpI, ageZ=data1_NA$ageZ, lacgr=data1_NA$lacgr, leg0=data1_NA$leg0, leg1=data1_NA$leg1,leg2=data1_NA$leg2)), tag = "data1_NA.data") 

# Create joint stack
join.stack_NA <- inla.stack(stackS, stackS_pred_NA)

# ModelS
formulaS <- as.formula(paste0(modelBase, " + f(fieldID, model = spdeStatS)-1"))


fitS_NA= inla(formula = formulaS, data = inla.stack.data(join.stack_NA),
               family = "normal", control.predictor =list(A=inla.stack.A(join.stack_NA),compute = T),
               control.family=list(list(hyper=hyperResVarGWS)),
               control.compute = list(dic=T,cpo=F, config=T), 
               control.fixed = list(expand.factor.strategy="inla"), verbose=T)

# save fitS_NA as R object
save(fitS_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitS_NA.RData")

#Get the prediction index

pred.ind_NA <- inla.stack.index(join.stack_NA, tag='data1_NA.data')$data
pheno_pred <- fitS_NA$summary.linear.predictor [pred.ind_NA,]

sum(is.na(pheno_pred$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred$mean
data1$milkZ_pred_sd <- pheno_pred$sd 
pheno <-   subset(data1, birthyear=="2016" | birthyear=="2017") 
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))  
length(unique(pheno$cowI))# 146 Cows expected

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotyped 
accuracy_fitS <- round(cor(pheno$milkZ,pheno$milkZ_pred),2) 
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitS # 0.86
R2_fitS<- summary(Coef)
R2_fitS #  0.744 
#CRPS
obs <- pheno$milkZ
pred<- subset(pheno,select= c(milkZ_pred,milkZ_pred_sd))
crps_fitS <- crps(obs,pred)
round((crps_fitS$CRPS),2) # 0.31

#-----------------------Foward validation model fitWS-------------------------------

# Make stack_pred
# StackWS_NA
stackWS_pred_NA = inla.stack(data = list(milkZ = NA),
                            A = list(A.pred,1),
                            effects = list(c(meshIndexWS, list(intercept = 1)),
                                           list(cowI = data1_NA$cowI,herdI = data1_NA$herdI, cowPeI = data1_NA$cowPeI,
                                                cyrsnI=data1_NA$cyrsnI, tyrmnI=data1_NA$tyrmnI,dgrpI= data1_NA$dgrpI, ageZ=data1_NA$ageZ, lacgr=data1_NA$lacgr, leg0=data1_NA$leg0, leg1=data1_NA$leg1,leg2=data1_NA$leg2)), tag = "data1_NA.data") 

# Create joint stack
join.stack_NA <- inla.stack(stackWS, stackWS_pred_NA)

# ModelS
formulaWS <- as.formula(paste0(modelBase, " + f(fieldID, model = spdeStatWS)-1"))


fitWS_NA= inla(formula = formulaWS, data = inla.stack.data(join.stack_NA),
              family = "normal", control.predictor =list(A=inla.stack.A(join.stack_NA),compute = T),
              control.family=list(list(hyper=hyperResVarGWS)),
              control.compute = list(dic=T,cpo=F, config=T), 
              control.fixed = list(expand.factor.strategy="inla"), verbose=T)

# save fitWS_NA as R object
save(fitWS_NA,file = "D:/Results_ADGG_Spatial/forward_valid/fitWS_NA.RData")

#Get the prediction index

pred.ind_NA <- inla.stack.index(join.stack_NA, tag='data1_NA.data')$data
pheno_pred <- fitWS_NA$summary.linear.predictor [pred.ind_NA,]

sum(is.na(pheno_pred$mean)) # 0 expected
data1$milkZ_pred <- pheno_pred$mean
data1$milkZ_pred_sd <- pheno_pred$sd 
pheno <-   subset(data1, birthyear=="2016" | birthyear=="2017") 
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))  
length(unique(pheno$cowI))# 146 Cows expected

# Correlating observed phenotype of cows in dgrp 1 with predicted phenotyped 
accuracy_fitWS <- round(cor(pheno$milkZ,pheno$milkZ_pred),2) 
Coef<- lm (pheno$milkZ~pheno$milkZ_pred)
accuracy_fitWS # 0.86
R2_fitWS<- summary(Coef)
R2_fitWS #  0.7451 
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

#data1_2016_2017 <- subset(data1_birthdate, birthyear==2016 | birthyear==2017 ) 

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
pheno<- subset(pheno,select= c(cowI,milkZ,milkZ_pred,milkZ_pred_sd))  
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

data1_variable_pca <- data1[,c(1,3,7,18)]
data1_variable_pca


dim(geno)
# 1911 664823
str(geno)
summary(geno)
# geno$ID as character
geno$cow <- as.character(geno$cow) 
str(geno)

# Delete Cases with Missingness  
geno_nomiss <-  na.omit(geno)
dim(geno_nomiss)
# 1911 664823
#Exclude Categorical Data
geno_sample <- geno_nomiss[,-c(1)]

#Run PCA
geno_pca <- prcomp(geno_sample, 
                     scale = TRUE)

#Summary of Analysis 
summary(geno_pca)

#Elements of PCA object 
names(geno_pca)

#Std Dev of Components 
geno_pca$sdev

#Eigenvectors 
geno_pca$rotation

#Std Dev and Mean of Variables 
geno_pca$center
geno_pca$scale

#Principal Component Scores
geno_pca$x


#Scree Plot of Variance 
fviz_eig(geno_pca, 
         addlabels = TRUE,
         ylim = c(0, 70))

#Biplot with Default Settings
fviz_pca_biplot(geno_pca)

#Biplot without Labeled Variables
#fviz_pca_biplot(geno_pca,
                label="var")

#Biplot with Colored Groups
#fviz_pca_biplot(geno_pca,
               # label="var",
               # habillage = geno_variable$dgrp)

# Biplot with Customized Colored Groups and Variables
#fviz_pca_biplot(biopsy_pca,
               # label="var",
               # habillage = biopsy_nomiss$class, 
                col.var = "black") +
 # scale_color_manual(values=c("orange", "purple"))



# Biplot with Customized Colored Groups and Variables
#fviz_pca_biplot(biopsy_pca,
               # label="var",
                #habillage = biopsy_nomiss$class, 
               # col.var = "black") +
  #scale_color_manual(values=c("orange", "purple"))

#-----------------------Animal Ranking------------------------------------------
# Baseline Model (B)
EBV_B <- data.frame(fitBase$summary.random$cowI)
EBV_B <- EBV_B[,1:2]
names(EBV_B)[2] <- "EBV_B"
# Remove EBVS of non-phenotyped animals
sel <- EBV_B$ID %in% data1$cowI
EBV_B <- EBV_B[sel, ]
dim(EBV_B)

# WCF
EBV_WCF <- data.frame(fitWCF$summary.random$cowI)
EBV_WCF <- EBV_WCF[,1:2]
names(EBV_WCF)[2] <- "EBV_WCF"
sel <- EBV_WCF$ID %in% data1$cowI
EBV_WCF <- EBV_WCF[sel, ]
dim(EBV_WCF)

#WCRI
EBV_WCRI <- data.frame(fitWCRI$summary.random$cowI)
EBV_WCRI <- EBV_WCRI[,1:2]
names(EBV_WCRI)[2] <- "EBV_WCRI"
sel <- EBV_WCRI$ID %in% data1$cowI
EBV_WCRI <- EBV_WCRI[sel, ]
dim(EBV_WCRI)

#WCRB
EBV_WCRB <- data.frame(fitWCRB$summary.random$cowI)
EBV_WCRB <- EBV_WCRB[,1:2]
names(EBV_WCRB)[2] <- "EBV_WCRB"
sel <- EBV_WCRB$ID %in% data1$cowI
EBV_WCRB <- EBV_WCRB[sel, ]
dim(EBV_WCRB)


# S

EBV_S <- data.frame(fitS$summary.random$cowI)
EBV_S <- EBV_S[,1:2]
names(EBV_S)[2] <- "EBV_S"
sel <- EBV_S$ID %in% data1$cowI
EBV_S <- EBV_S[sel, ]
dim(EBV_S)

# WS

EBV_WS <- data.frame(fitWS$summary.random$cowI)
EBV_WS <- EBV_WS[,1:2]
names(EBV_WS)[2] <- "EBV_WS" 
sel <- EBV_WS$ID %in% data1$cowI
EBV_WS <- EBV_WS[sel, ]
dim(EBV_WS)

#Spatial effect by Cow for model S

SPDE_cow_S <- SPDE_cow[,1:2]

# EBVs and Spatial effect

EBV_Spatial<- data.frame(EBV_B$ID,EBV_B$EBV_B,EBV_WCF$EBV_WCF, EBV_WCRI$EBV_WCRI,EBV_WCRB$EBV_WCRB,EBV_S$EBV_S, EBV_WS$EBV_WS,SPDE_cow_S$spdeS_mean)
names(EBV_Spatial)[1:8] <- c("ID","EBV_B","EBV_WCF","EBV_WCRI","EBV_WCRB","EBV_S","EBV_WS","Spatial_effect_S")

#Spatial effect by  herd for model S



#create matrix of correlation coefficients and p-values for EBVs between models
# How to Create a Correlation Matrix in R (4 Examples) - Statology
# https://www.statology.org/correlation-matrix-in-r/
  
  # Pearson Correlation
rcorr(as.matrix(EBV_Spatial))

# Spearman's Rank  Correlation

cor(EBV_Spatial, method = "spearman")



#Animal Ranking
  
# Ranking EBV_B
  EBV_B_Rank <- EBV_B[order(EBV_B$EBV_B,decreasing = TRUE),] 

# Ranking EBV_WCF
EBV_WCF_Rank <- EBV_WCF[order(EBV_WCF$EBV_WCF,decreasing = TRUE),] 

# Ranking EBV_WCRI
EBV_WCRI_Rank <- EBV_WCRI[order(EBV_WCRI$EBV_WCRI,decreasing = TRUE),] 

# Ranking EBV_WCRB
EBV_WCRB_Rank <- EBV_WCRB[order(EBV_WCRB$EBV_WCRB,decreasing = TRUE),]

# Ranking EBV_S
EBV_S_Rank <- EBV_S[order(EBV_S$EBV_S,decreasing = TRUE),] 

# Ranking EBV_WS
EBV_WS_Rank <- EBV_WS[order(EBV_WS$EBV_WS,decreasing = TRUE),] 








corr1 <- cor.test(x=EBV_B$mean, y=EBV_WCF$mean, method = 'spearman')
corr1
corr2 <- cor.test(x=EBV_B$mean, y=EBV_WCRI$mean, method = 'spearman')
corr2

corr1 <- cor.test(x=EBV_B$mean, y=EBV_S$mean, method = 'spearman')
corr1
corr2 <- cor.test(x=EBV_B$mean, y=EBV_S$mean, method = 'spearman')
corr2
corr3 <- cor.test(x=EBV_S$mean, y=EBV_WS$mean, method = 'spearman')
corr3
