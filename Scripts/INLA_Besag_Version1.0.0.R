## Isidore,  Eva, Thiago,  Ivan and Gregor                                   #
# Version 1.0.0                                                        #
# Date: 09/02/2022 
#================================================================================================
# Running the models using INLA
#================================================================================================
##Things to ask:
## how to compute Ainv with  1.000.000 animals (Ivan)
## how to fit AgeNestedParity (Age as covariate nested within parity)
## no prior values included (discuss this with Thiago and Ivan)
## check with Gregor for the INLA code if things are ok (I am not sure everything is correct)
## Don't run INLA on full dataset when checking how things work in INLA may crash and takes long


# BASE MODEL # Keep HERD_ID or HYS in the model and not both
#Import the Holstein pedigree data of South Africa
HST_PED_IH <- read.csv(file = "HST_PED.csv", stringsAsFactors = FALSE)
#Extract ANI_ID, SIRE_ID, DAM_ID 
HST_PED_HST_SA_IH<- HST_PED_IH %>%
  select(ANIM_ID,SIRE_ID, DAM_ID)
HST_PED_HST_SA_IH
HST_PED_HST_SA_IH<- HST_PED_HST_SA_IH %>% distinct()
colnames(HST_PED_HST_SA_IH)
# Rename  HST_PED_HST_SA$ANIM_ID to HST_PED_HST_SA$ANIMAL_NUMBER
names(HST_PED_HST_SA_IH)[names(HST_PED_HST_SA_IH) == 'ANIM_ID'] <- 'COW_ID'
HST_PED_HST_SA_IH 

#Importing data_clean_final
data_clean_final <- read.csv(file = "data_clean_final.csv", stringsAsFactors = FALSE)
colnames(data_clean_final)
dataset1<- subset(data_clean_final, PARITY==1 |PARITY==2)

# Build relationship matrix from pedigree
#ped = read.csv("sorod.csv") #### Line 17 t0 45 to get Ainv
#head(ped)


# Relationship matrix from pedigree should be built from idziva (COW_ID), idocea (Dam_ID) and idmata (Sire_ID)

pedHST_SA = data.frame(label =HST_PED_HST_SA_IH$COW_ID, dam = HST_PED_HST_SA_IH$DAM_ID, sire = HST_PED_HST_SA_IH$SIRE_ID)
library("pedigree")
library(package="pedigreemm")
pedHST_SA[pedHST_SA == 0] <- NA

pedHST_SA3 = editPed(sire = pedHST_SA$sire, dam = pedHST_SA$dam, label = pedHST_SA$label) ## Stopped here

# Pedigree object for pedigreemm functions
pedHSTMM = pedigree(label = pedHST_SA3$label,dam=pedHST_SA3$dam,pedHST_SA3$sire)


# Precision matrix (A inverse)
Tinv    = as(pedHSTMM, "sparseMatrix") ## T^{-1} in A^{-1} = (T^{-1})' D^{-1} T^{-1}
DmatTemp = pedigreemm::Dmat(pedHSTMM)
D       = Diagonal(x=DmatTemp)   ## D      in A = TDT'
Dinv    = solve(D)                   ## ...
Ainv = t(Tinv) %*% Dinv %*% Tinv  ## ... Ainv to be save 

# Labelling of the A matrix is according to pedPMM
labelAinv = as.numeric(pedHSTMM@label)
mappingAinv = data.frame(COW_ID = labelAinv, rowNumberAinv = 1:length(labelAinv))

data_clean_final = merge(data_clean_final, mappingAinv)


# BASE MODEL # Keep HERD_ID or HYS in the model and not both
nb.map <- poly2nb(map)
nb2INLA("map.graph",nb.map)
gREGION_ID <- inla.read.graph(filename = "map.graph")

# Herd as random (G1 models)
formulabase1 = MILK_YLD_305D ~ 1 + AgeNestedParity  + 
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv) + 
  f(code_HERD_ID, model = "iid") +
  f(COW_ID, model = "iid")

formulabase1Alt1 = MILK_YLD_305D ~ 1 + AgeNestedParity  + REGION_ID +
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv) + 
  f(code_HERD_ID, model = "iid") +
  f(COW_ID, model = "iid")

formulabase1Alt2 = MILK_YLD_305D ~ 1 + AgeNestedParity  + 
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv) + 
  f(code_HERD_ID, model = "iid") +
  f(COW_ID, model = "iid")+
  f(region, model = "iid")

formulabase1Alt3 = MILK_YLD_305D ~ 1 + AgeNestedParity  + 
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv) + 
  f(code_HERD_ID, model = "iid") +
  f(COW_ID, model = "iid") +
  f(region, model = "besag", graph = gREGION_ID, scale.model = FALSE)


### HYS as random (G2 models)

formulabase2 = MILK_YLD_305D ~ 1 + AgeNestedParity  + 
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv) + 
  f(HYS, model = "iid") +
  f(COW_ID, model = "iid")

formulabase2Alt1 = MILK_YLD_305D ~ 1 + AgeNestedParity  + REGION_ID +
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv) + 
  f(HYS, model = "iid") +
  f(COW_ID, model = "iid")

formulabase2Alt2 = MILK_YLD_305D ~ 1 + AgeNestedParity  + 
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv) + 
  f(HYS, model = "iid") +
  f(COW_ID, model = "iid")+
  f(region, model = "iid")

formulabase2Alt3 = MILK_YLD_305D ~ 1 + AgeNestedParity  + 
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv) + 
  f(HYS, model = "iid") +
  f(COW_ID, model = "iid") +
  f(region, model = "besag", graph = gREGION_ID, scale.model = FALSE)


### Code to run each model
model_G1 <- inla(formulabase1,  family = "gaussian", data = data_clean_final, control.compute = list(dic=T,cpo=F), verbose=T)

# dic = T: computes DIC value, which we can use for comparing the fit of the model
# verbose = T: gives you real time log on what is happening, you may not need it when running on a server
