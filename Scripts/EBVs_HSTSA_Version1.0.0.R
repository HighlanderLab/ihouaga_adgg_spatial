#=======================================================================
# EBVS_Holstein (South Africa) data analysis                                #
# Isidore,  Thiago,  Ivan and Gregor                                   #
# Version 1.0.0                                                        #
# Date: 09/02/2022                                                     #
# EBs manipulation                                                         #
#=======================================================================
#-----------------------------------------------------------------------
# Uploading packages   
#-----------------------------------------------------------------------
library(dplyr)
library(stringr)
library(ggplot2)     
library(lme4)     #for mixed model
library(spData)
library(tidyr)
library(devtools)
library(cowplot)
library(psych)
(.packages()) # Checking loaded packages
#=======================================================================
# Uploading datasets
#=======================================================================
# Importing EBVs for each model

EBV_G1 <- read.csv(file = "estimates_G1.csv", stringsAsFactors = FALSE)
EBV_G2 <- read.csv(file = "estimates_G2.csv", stringsAsFactors = FALSE)
EBV_GR1 <- read.csv(file = "estimates_GR1.csv", stringsAsFactors = FALSE)
EBV_GR2 <- read.csv(file = "estimates_GR2.csv", stringsAsFactors = FALSE)
EBV_GRU1 <- read.csv(file = "estimates_GRU1.csv", stringsAsFactors = FALSE)
EBV_GRU2 <- read.csv(file = "estimates_GRU2.csv", stringsAsFactors = FALSE)
EBV_GRC1 <- read.csv(file = "estimates_GRC1.csv", stringsAsFactors = FALSE)
EBV_GRC2 <- read.csv(file = "estimates_GRC2.csv", stringsAsFactors = FALSE)
EBV_GRC1_NS <- read.csv(file = "estimates_GRC1_NS.csv", stringsAsFactors = FALSE)# file to be added
EBV_GRC2_NS <- read.csv(file = "estimates_GRC2_NS.csv", stringsAsFactors = FALSE)# file to be added
#=======================================================================
# Data organization
#=======================================================================
## EBV_G1_anim
EBV_G1$factor_G1<- as.factor(EBV_G1$factor_G1)
EBV_G1$level_G1<- as.factor(EBV_G1$level_G1)
str(EBV_G1)
EBV_G1_anim <- EBV_G1 %>% 
  filter(factor_G1 == 4)
summary(EBV_G1_anim)
var(EBV_G1_anim$EBV_G1)
## EBV_G2_anim
EBV_G2$factor_G2<- as.factor(EBV_G2$factor_G2)
EBV_G2$level_G2<- as.factor(EBV_G2$level_G2)
str(EBV_G2)
EBV_G2_anim <- EBV_G2 %>% 
  filter(factor_G2 == 4)
summary(EBV_G2_anim)

## EBV_GR1_anim
EBV_GR1$factor_GR1<- as.factor(EBV_GR1$factor_GR1)
EBV_GR1$level_GR1<- as.factor(EBV_GR1$level_GR1)
str(EBV_GR1)
summary(EBV_GR1$factor_GR1)
EBV_GR1_anim <- EBV_GR1 %>% 
  filter(factor_GR1 == 5)
summary(EBV_GR1_anim)  

## EBV_GR2_anim
EBV_GR2$factor_GR2<- as.factor(EBV_GR2$factor_GR2)
EBV_GR2$level_GR2<- as.factor(EBV_GR2$level_GR2)
str(EBV_GR2)
summary(EBV_GR2$factor_GR2)
EBV_GR2_anim <- EBV_GR2 %>% 
  filter(factor_GR2 == 5)
summary(EBV_GR2_anim)  

## EBV_GRU1_anim
EBV_GRU1$factor_GRU1<- as.factor(EBV_GRU1$factor_GRU1)
EBV_GRU1$level_GRU1<- as.factor(EBV_GRU1$level_GRU1)
str(EBV_GRU1)
summary(EBV_GRU1$factor_GRU1)
EBV_GRU1_anim <- EBV_GRU1 %>% 
  filter(factor_GRU1 == 5)
summary(EBV_GRU1_anim) 

## EBV_GRU2_anim
EBV_GRU2$factor_GRU2<- as.factor(EBV_GRU2$factor_GRU2)
EBV_GRU2$level_GRU2<- as.factor(EBV_GRU2$level_GRU2)
str(EBV_GRU2)
summary(EBV_GRU2$factor_GRU2)
EBV_GRU2_anim <- EBV_GRU2 %>% 
  filter(factor_GRU2 == 5)
summary(EBV_GRU2_anim)

## EBV_GRC1_anim
EBV_GRC1$factor_GRC1<- as.factor(EBV_GRC1$factor_GRC1)
EBV_GRC1$level_GRC1<- as.factor(EBV_GRC1$level_GRC1)
str(EBV_GRC1)
summary(EBV_GRC1$factor_GRC1)
EBV_GRC1_anim <- EBV_GRC1 %>% 
  filter(factor_GRC1 == 4)
summary(EBV_GRC1_anim) 

## EBV_GRC2_anim
EBV_GRC2$factor_GRC2<- as.factor(EBV_GRC2$factor_GRC2)
EBV_GRC2$level_GRC2<- as.factor(EBV_GRC2$level_GRC2)
str(EBV_GRC2)
summary(EBV_GRC2$factor_GRC2)
EBV_GRC2_anim <- EBV_GRC2 %>% 
  filter(factor_GRC2 == 4)
summary(EBV_GRC2_anim) 

EBV_anim<- cbind(EBV_G1_anim,EBV_G2_anim,EBV_GR1_anim, EBV_GR2_anim,EBV_GRU1_anim, EBV_GRU2_anim, EBV_GRC1_anim, EBV_GRC2_anim) 

EBVs<- EBV_anim %>% select(level_G1, EBV_G1, EBV_GR1,EBV_GRU1,EBV_GRC1,EBV_G2,EBV_GR2,EBV_GRU2,EBV_GRC2)
colnames(EBVs)[1] <- "renum_ID"  
write.csv(EBVs,"/Users/ihouaga/Documents/Project1/EBVs.csv") #Save EBvs

##Estimate Region GRC1 and GRC2

region_estimates_GRC1 <- EBV_GRC1 %>% 
  filter(factor_GRC1 == 6) %>%
  select(level_GRC1,EBV_GRC1) 

region_estimates_GRC2 <- EBV_GRC2 %>% 
  filter(factor_GRC2 == 6) %>%
  select(level_GRC2,EBV_GRC2) 
region_estimates_GRC2
#Rename level_GRC2 to level_GRC1 before merge
colnames(region_estimates_GRC2)[1] <- "level_GRC1"
RC_E <- merge(region_estimates_GRC1,region_estimates_GRC2, by= "level_GRC1", all.x = TRUE)
colnames(RC_E)[2] <- "E_GRC1"
colnames(RC_E)[3] <- "E_GRC2"
write.csv(RC_E,"/Users/ihouaga/Documents/Project1/region_estimates.csv" )

#=======================================================================
# Exporting Pedigre_renum IDs and merging with EBVs
#=======================================================================
renum_ori_ID<- read.csv(file = "renum_ori_ID.csv", stringsAsFactors = FALSE)
### EBVs with animal ID
EBVs_ID<- merge(EBVs,renum_ori_ID, by="renum_ID", all.x = TRUE)
##Region estimates
region_estimates<- merge(region_estimates_GRC1,region_estimates_GRC2, by="level_GRC1", all.x=TRUE) 
colnames(region_estimates)[1]<- "region"
### correlations (EBVs)
cor.test(EBVs_ID$EBV_G1,EBVs$EBV_GR1)
cor.test(EBVs_ID$EBV_G1,EBVs$EBV_GRU1)
cor.test(EBVs_ID$EBV_G1,EBVs$EBV_GRC1)
cor.test(EBVs_ID$EBV_GR1,EBVs$EBV_GRU1)
cor.test(EBVs_ID$EBV_GR1,EBVs$EBV_GRC1)
cor.test(EBVs_ID$EBV_GRU1,EBVs$EBV_GRC1)
EBVs_cor <- EBVs_ID %>%  select(-c(renum_ID))
corPlot(EBVs_cor,cex = 1.2) ## Correlation between EBVs
### Density plot of animals grouped by Region
data_region<- data_clean_test_3 %>% select(COW_ID,region)
colnames(data_region)[1] <- "ori_anim_ID"

length(unique(data_clean_test_3$COW_ID))

EBVs_ID <-  EBVs_ID %>% relocate(ori_anim_ID)
EBVs_region_ID<- merge(data_region,EBVs_ID, by="ori_anim_ID", all.x= TRUE) ## cows with orignal anima ID and region ID

# Ploting EBvs cows phenotyped by region, GRC1 model
EBVs_region_ID %>% ggplot(aes(x = EBV_GRC1, colour=region)) +
  geom_density()

# Ploting EBvs cows phenotyped by region, GRC2 model
EBVs_region_ID %>% ggplot(aes(x = EBV_GRC2, colour=region)) +
  geom_density()

# Ploting EBvs cows phenotyped by region, G1 model
EBVs_region_ID %>% ggplot(aes(x = EBV_G1, colour=region)) +
  geom_density()
# Ploting EBvs cows phenotyped by region, G1 model
EBVs_region_ID %>% ggplot(aes(x = EBV_G2, colour=region)) +
  geom_density()

# Ploting EBvs cows phenotyped by region, GR1 model
EBVs_region_ID %>% ggplot(aes(x = EBV_GR1, colour=region)) +
  geom_density()
# Ploting EBvs cows phenotyped by region, GR2 model
EBVs_region_ID %>% ggplot(aes(x = EBV_GR2, colour=region)) +
  geom_density()

# Ploting EBvs cows phenotyped by region, GRU1 model
EBVs_region_ID %>% ggplot(aes(x = EBV_GRU1, colour=region)) +
  geom_density()
# Ploting EBvs cows phenotyped by region, GRU2 model
EBVs_region_ID %>% ggplot(aes(x = EBV_GRU2, colour=region)) +
  geom_density()

#Correlation between regions_estimates
sol_region_fixed_random<- read.csv(file = "sol_region_fixed_random.csv", stringsAsFactors = FALSE)
sol_region_GRC<- read.csv(file = "sol_region_GRC.csv", stringsAsFactors = FALSE)

sol_regions<- merge(sol_region_GRC,sol_region_fixed_random, by="region_ID", all.x=TRUE)
sol_regions<- sol_regions %>% select(region_ID, sol_region_GR1,sol_region_GRU1,sol_region_GRC1,sol_region_GRC1_NS,sol_region_GR2,sol_region_GRU2,sol_region_GRC2,sol_region_GRC2_NS)
sol_regions_cor <- sol_regions %>%  select(-c(region_ID))
corPlot(sol_regions_cor,cex = 0.5) ## Correlation between regions estimates 
  cor.test(sol_regions_cor$sol_region_GR1,sol_regions_cor$sol_region_GRU1)
  cor.test(sol_regions_cor$sol_region_GR1,sol_regions_cor$sol_region_GRC1)
  cor.test(sol_regions_cor$sol_region_GRU1,sol_regions_cor$sol_region_GRC1)
  
  
  #=======================================================================
  # Prediction Error Variance (PEV) of each model 
  #=======================================================================
  # Importing SE of EBVs for each model
  
  SE_G1 <- read.csv(file = "SE_G1.csv", stringsAsFactors = FALSE)
  SE_G2 <- read.csv(file = "SE_G2.csv", stringsAsFactors = FALSE)
  SE_GR1 <- read.csv(file = "SE_GR1.csv", stringsAsFactors = FALSE)
  SE_GR2 <- read.csv(file = "SE_GR2.csv", stringsAsFactors = FALSE)
  SE_GRU1 <- read.csv(file = "SE_GRU1.csv", stringsAsFactors = FALSE)
  SE_GRU2 <- read.csv(file = "SE_GRU2.csv", stringsAsFactors = FALSE)
  SE_GRC1 <- read.csv(file = "SE_GRC1.csv", stringsAsFactors = FALSE)
  SE_GRC2 <- read.csv(file = "SE_GRC2.csv", stringsAsFactors = FALSE)
  SE_GRC1_NS <- read.csv(file = "SE_GRC1_NS.csv", stringsAsFactors = FALSE)
  SE_GRC2_NS <- read.csv(file = "SE_GRC2_NS.csv", stringsAsFactors = FALSE)
  #=======================================================================
  # PEV Data organization
  #=======================================================================
  ##SE_G1_anim
  SE_G1$factor_G1<- as.factor(SE_G1$factor_G1)
  SE_G1$level_G1<- as.factor(SE_G1$level_G1)
  str(SE_G1)
 SE_G1_anim <- SE_G1 %>% 
    filter(factor_G1 == 4)
  summary(SE_G1_anim)
  
  ##SE_G2_anim
SE_G2$factor_G2<- as.factor(SE_G2$factor_G2)
SE_G2$level_G2<- as.factor(SE_G2$level_G2)
  str(SE_G2)
  SE_G2_anim <- SE_G2 %>% 
    filter(factor_G2 == 4)
  summary(SE_G2_anim)
  
  ## SE_GR1_anim
  SE_GR1$factor_GR1<- as.factor(SE_GR1$factor_GR1)
  SE_GR1$level_GR1<- as.factor(SE_GR1$level_GR1)
  str(SE_GR1)
  summary(SE_GR1$factor_GR1)
  SE_GR1_anim <- SE_GR1 %>% 
    filter(factor_GR1 == 5)
  summary(SE_GR1_anim)  
  
  ## SE_GR2_anim
  SE_GR2$factor_GR2<- as.factor(SE_GR2$factor_GR2)
  SE_GR2$level_GR2<- as.factor(SE_GR2$level_GR2)
  str(SE_GR2)
  summary(SE_GR2$factor_GR2)
  SE_GR2_anim <- SE_GR2 %>% 
    filter(factor_GR2 == 5)
  summary(SE_GR2_anim)  
  
  ##SE_GRU1_anim
  SE_GRU1$factor_GRU1<- as.factor(SE_GRU1$factor_GRU1)
  SE_GRU1$level_GRU1<- as.factor(SE_GRU1$level_GRU1)
  str(SE_GRU1)
  summary(SE_GRU1$factor_GRU1)
  SE_GRU1_anim <- SE_GRU1 %>% 
    filter(factor_GRU1 == 5)
  summary(SE_GRU1_anim) 
  
  ## SE_GRU2_anim
  SE_GRU2$factor_GRU2<- as.factor(SE_GRU2$factor_GRU2)
  SE_GRU2$level_GRU2<- as.factor(SE_GRU2$level_GRU2)
  str(SE_GRU2)
  summary(SE_GRU2$factor_GRU2)
  SE_GRU2_anim <-SE_GRU2 %>% 
    filter(factor_GRU2 == 5)
  summary(SE_GRU2_anim)
  
  ## SE_GRC1_anim
  SE_GRC1$factor_GRC1<- as.factor(SE_GRC1$factor_GRC1)
  SE_GRC1$level_GRC1<- as.factor(SE_GRC1$level_GRC1)
  str(SE_GRC1)
  summary(SE_GRC1$factor_GRC1)
  SE_GRC1_anim <- SE_GRC1 %>% 
    filter(factor_GRC1 == 4)
  summary(SE_GRC1_anim) 
  
  ## SE_GRC1_NS_anim
  SE_GRC1_NS$factor_GRC1_NS<- as.factor(SE_GRC1_NS$factor_GRC1_NS)
  SE_GRC1_NS$level_GRC1_NS<- as.factor(SE_GRC1_NS$level_GRC1_NS)
  str(SE_GRC1_NS)
  summary(SE_GRC1$factor_GRC1)
  SE_GRC1_NS_anim <- SE_GRC1_NS %>% 
    filter(factor_GRC1_NS == 4)
  summary(SE_GRC1_NS_anim) 
  
  
  ## SE_GRC2_anim
  SE_GRC2$factor_GRC2<- as.factor(SE_GRC2$factor_GRC2)
  SE_GRC2$level_GRC2<- as.factor(SE_GRC2$level_GRC2)
  str(SE_GRC2)
  summary(SE_GRC2$factor_GRC2)
  SE_GRC2_anim <- SE_GRC2 %>% 
    filter(factor_GRC2 == 4)
  summary(SE_GRC2_anim) 
  
  ## SE_GRC2_NS_anim
  SE_GRC2_NS$factor_GRC2_NS<- as.factor(SE_GRC2_NS$factor_GRC2_NS)
  SE_GRC2_NS$level_GRC2_NS<- as.factor(SE_GRC2_NS$level_GRC2_NS)
  str(SE_GRC2_NS)
  summary(SE_GRC2_NS$factor_GRC2_NS)
  SE_GRC2_NS_anim <- SE_GRC2_NS %>% 
    filter(factor_GRC2_NS == 4)
  summary(SE_GRC2_NS_anim) 
  
  SE_anim<- cbind(SE_G1_anim,SE_G2_anim,SE_GR1_anim, SE_GR2_anim,SE_GRU1_anim, SE_GRU2_anim, SE_GRC1_anim, SE_GRC2_anim,SE_GRC1_NS_anim, SE_GRC2_NS_anim) 
  
  SE<- SE_anim %>% select(level_G1, SE_G1, SE_GR1,SE_GRU1,SE_GRC1,SE_GRC1_NS, SE_G2,SE_GR2,SE_GRU2,SE_GRC2,SE_GRC2_NS)
  colnames(SE)[1] <- "renum_ID"  
  
  ## Calculate PEV
  PEV<- SE  %>%
    mutate(PEV_G1=SE_G1^2, PEV_GR1=SE_GR1^2,PEV_GRU1=SE_GRU1^2,PEV_GRC1=SE_GRC1^2,PEV_GRC1_NS=SE_GRC1_NS^2,PEV_G2=SE_G2^2,PEV_GR2=SE_GR2^2,PEV_GRU2=SE_GRU2^2,PEV_GRC2=SE_GRC2^2,PEV_GRC2_NS=SE_GRC2_NS^2)
  
  
  write.csv(PEV,"/Users/ihouaga/Documents/Project1/PEV.csv") #Save EBvs
  
  summary(PEV)
  
#================================================================================================
# Grouping EBV by region
#================================================================================================ 
  EBVs_region<- EBVs_region_ID %>% 
    group_by(region) %>%  
  summarise(mean_EBV_GRC1=mean(EBV_GRC1),sd_EBV_GRC1=sd(EBV_GRC1), mean_EBV_GRC2=mean(EBV_GRC2),sd_EBV_GRC2=sd(EBV_GRC2))
  write.csv(EBVs_region,"/Users/ihouaga/Documents/Project1/EBVs_region.csv") #Save EBVs_region
 
  #================================================================================================
  # Grouping PEV by region
  #================================================================================================  
  EBVs_region_ID_2<- EBVs_region_ID %>% relocate(renum_ID)
EBVs_PEV_region <- merge(EBVs_region_ID_2,PEV, by="renum_ID", all.x=TRUE)  

PEV_region<- EBVs_PEV_region %>% 
  group_by(region)%>%
  summarise(mean_PEV_GRC1=mean(PEV_GRC1),mean_PEV_GRC2=mean(PEV_GRC2))
write.csv(PEV_region,"/Users/ihouaga/Documents/Project1/PEV_region.csv")

# Scatterplot (x-y) between sd of EBV and average PEV per region 
EBVs_PEV_region_sd_mean <- merge(EBVs_region, PEV_region, by="region", all.x=TRUE)

ggplot(EBVs_PEV_region_sd_mean, aes(x = sd_EBV_GRC1, y = mean_PEV_GRC1, colour=region)) +
  geom_point()
#================================================================================================
# Grouping PEV by SIRE, COW_ID and non phenotype ID
#================================================================================================  
## Import PEV and data_clean_final
PEV<- read.csv(file = "PEV.csv", stringsAsFactors = FALSE)
data_clean_final<- read.csv(file = "data_clean_final.csv", stringsAsFactors = FALSE)
## Import animal IDs (renum+Original)
ani_ID<- read.csv(file = "renum_ori_ID.csv", stringsAsFactors = FALSE)
PEV_ani_ID<- merge(PEV,ani_ID, by="renum_ID", all.x=TRUE)

PEV_ani_ID<- PEV_ani_ID %>% select(ori_anim_ID, everything())
## Sire ID
sire_ID <- data_clean_final %>%
  select(SIRE_ID)
sire_ID_uniq<-  sire_ID %>% distinct()  # 5978 Sire

sire_ID_uniq<- sire_ID_uniq %>%
filter(SIRE_ID >0) ## Remove missing sire

 COW_ID<- data_clean_final %>%
  select(COW_ID)
COW_ID_uniq<-  COW_ID %>% distinct()

ani_ID_ori <- ani_ID %>%
  select(ori_anim_ID) 
# Save COW_ID_uniqu and sire_ID_uniq

write.csv(COW_ID_uniq,"/Users/ihouaga/Documents/Project1/COW_ID_uniq.csv")
write.csv(sire_ID_uniq,"/Users/ihouaga/Documents/Project1/sire_ID_uniq.csv")

#rename Colname in COW_ID_unique 

colnames(COW_ID_uniq)[1] <- "ori_anim_ID"
COW_ID_uniq
COW_ID_uniq$mean <- rep(1, nrow(COW_ID_uniq))
Non_pheno_anim_ID<- merge (ani_ID_ori,COW_ID_uniq, by="ori_anim_ID", all.x=TRUE) ## NA= non pheno animal
Non_pheno_anim_ID<- Non_pheno_anim_ID %>% 
  filter(is.na(mean)) %>% 
  select(-c(mean))
Non_pheno_anim_ID## 98146
colnames(sire_ID_uniq)[1] <- "ori_anim_ID" 
sire_ID_uniq$mean <- rep(1, nrow(sire_ID_uniq))
Non_pheno_COW_ID<- merge (Non_pheno_anim_ID,sire_ID_uniq, by="ori_anim_ID", all.x=TRUE) ## NA= non pheno animal

Non_pheno_COW_ID<- Non_pheno_COW_ID  %>% 
  filter(is.na(mean)) %>% 
  select(-c(mean))
Non_pheno_COW_ID# 92169 non phenotyped cows

#Phenotyped cows=98632
 # Sire: 5977
# 5977+98632+92169= 196778
# Save non phenotyped animal IDs  
write.csv(Non_pheno_anim_ID,"/Users/ihouaga/Documents/Project1/Non_pheno_anim_ID.csv")

## PEV sire
PEV_sire<- merge(sire_ID_uniq, PEV_ani_ID, by="ori_anim_ID", all.x=TRUE)
#PEV_pheno_cow
PEV_pheno_cow<- merge(COW_ID_uniq, PEV_ani_ID, by="ori_anim_ID", all.x=TRUE)

#PEV_non_pheno_anim 
PEV_non_pheno_anim<- merge(Non_pheno_anim_ID, PEV_ani_ID, by="ori_anim_ID", all.x=TRUE)

#PEV_non_pheno_COW 
PEV_non_pheno_cow<- merge(Non_pheno_COW_ID, PEV_ani_ID, by="ori_anim_ID", all.x=TRUE)

# Correlation Sire_PEV
cor.test(PEV_sire$PEV_G1,PEV_sire$PEV_GR1)
cor.test(PEV_sire$PEV_G1,PEV_sire$PEV_GRU1)
cor.test(PEV_sire$PEV_G1,PEV_sire$PEV_GRC1)
cor.test(PEV_sire$PEV_GR1,PEV_sire$PEV_GRU1)
cor.test(PEV_sire$PEV_GR1,PEV_sire$PEV_GRC1)
cor.test(PEV_sire$PEV_GRU1,PEV_sire$PEV_GRC1)

#Correlation PEV_pheno_cow
cor.test(PEV_pheno_cow$PEV_G1,PEV_pheno_cow$PEV_GR1)
cor.test(PEV_pheno_cow$PEV_G1,PEV_pheno_cow$PEV_GRU1)
cor.test(PEV_pheno_cow$PEV_G1,PEV_pheno_cow$PEV_GRC1)
cor.test(PEV_pheno_cow$PEV_GR1,PEV_pheno_cow$PEV_GRU1)
cor.test(PEV_pheno_cow$PEV_GR1,PEV_pheno_cow$PEV_GRC1)
cor.test(PEV_pheno_cow$PEV_GRU1,PEV_pheno_cow$PEV_GRC1)

#Correlation PEV_non_pheno_anim
cor.test(PEV_non_pheno_anim$PEV_G1,PEV_non_pheno_anim$PEV_GR1)
cor.test(PEV_non_pheno_anim$PEV_G1,PEV_non_pheno_anim$PEV_GRU1)
cor.test(PEV_non_pheno_anim$PEV_G1,PEV_non_pheno_anim$PEV_GRC1)
cor.test(PEV_non_pheno_anim$PEV_GR1,PEV_non_pheno_anim$PEV_GRU1)
cor.test(PEV_non_pheno_anim$PEV_GR1,PEV_non_pheno_anim$PEV_GRC1)
cor.test(PEV_non_pheno_anim$PEV_GRU1,PEV_non_pheno_anim$PEV_GRC1)

#Correlation PEV_non_pheno_anim
cor.test(PEV_non_pheno_cow$PEV_G1,PEV_non_pheno_cow$PEV_GR1)
cor.test(PEV_non_pheno_cow$PEV_G1,PEV_non_pheno_cow$PEV_GRU1)
cor.test(PEV_non_pheno_cow$PEV_G1,PEV_non_pheno_cow$PEV_GRC1)
cor.test(PEV_non_pheno_cow$PEV_GR1,PEV_non_pheno_cow$PEV_GRU1)
cor.test(PEV_non_pheno_cow$PEV_GR1,PEV_non_pheno_cow$PEV_GRC1)
cor.test(PEV_non_pheno_cow$PEV_GRU1,PEV_non_pheno_cow$PEV_GRC1)

# HYS models Sire_PEV
cor.test(PEV_sire$PEV_G2,PEV_sire$PEV_GR2)
cor.test(PEV_sire$PEV_G2,PEV_sire$PEV_GRU2)
cor.test(PEV_sire$PEV_G2,PEV_sire$PEV_GRC2)
cor.test(PEV_sire$PEV_GR2,PEV_sire$PEV_GRU2)
cor.test(PEV_sire$PEV_GR2,PEV_sire$PEV_GRC2)
cor.test(PEV_sire$PEV_GRU2,PEV_sire$PEV_GRC2)

#================================================================================================
# Reliability and accuracy  SIRE_EBV
#================================================================================================  
summary(PEV_sire$PEV_G1)
summary(PEV_sire$PEV_GR1)
summary(PEV_sire$PEV_GRU1)
summary(PEV_sire$PEV_GRC1)
summary(PEV_sire$PEV_GRC1_NS)


summary(PEV_sire$PEV_G2)
summary(PEV_sire$PEV_GR2)
summary(PEV_sire$PEV_GRU2)
summary(PEV_sire$PEV_GRC2)
summary(PEV_sire$PEV_GRC2_NS)
#================================================================================================
# Reliability and accuracy Phenoyped Cow
#================================================================================================  
summary(PEV_pheno_cow$PEV_G1)
summary(PEV_pheno_cow$PEV_GR1)
summary(PEV_pheno_cow$PEV_GRU1)
summary(PEV_pheno_cow$PEV_GRC1)
summary(PEV_pheno_cow$PEV_GRC1_NS)


summary(PEV_pheno_cow$PEV_G2)
summary(PEV_pheno_cow$PEV_GR2)
summary(PEV_pheno_cow$PEV_GRU2)
summary(PEV_pheno_cow$PEV_GRC2)
summary(PEV_pheno_cow$PEV_GRC2_NS)


#================================================================================================
# Reliability and accuracy non pheno-anim_ID
#================================================================================================  
summary(PEV_non_pheno_anim$PEV_G1)
summary(PEV_non_pheno_anim$PEV_GR1)
summary(PEV_non_pheno_anim$PEV_GRU1)
summary(PEV_non_pheno_anim$PEV_GRC1)
summary(PEV_non_pheno_anim$PEV_GRC1_NS)

summary(PEV_non_pheno_anim$PEV_G2)
summary(PEV_non_pheno_anim$PEV_GR2)
summary(PEV_non_pheno_anim$PEV_GRU2)
summary(PEV_non_pheno_anim$PEV_GRC2)
summary(PEV_non_pheno_anim$PEV_GRC2_NS)

#================================================================================================
# Reliability and accuracy non pheno-cow
#================================================================================================  
summary(PEV_non_pheno_cow$PEV_G1)
summary(PEV_non_pheno_cow$PEV_GR1)
summary(PEV_non_pheno_cow$PEV_GRU1)
summary(PEV_non_pheno_cow$PEV_GRC1)
summary(PEV_non_pheno_cow$PEV_GRC1_NS)

summary(PEV_non_pheno_cow$PEV_G2)
summary(PEV_non_pheno_cow$PEV_GR2)
summary(PEV_non_pheno_cow$PEV_GRU2)
summary(PEV_non_pheno_cow$PEV_GRC2)
summary(PEV_non_pheno_cow$PEV_GRC2_NS)
#=======================================================================
# EBVs and prediction Error Variance (PEV) of each model (G1Mu,GR1Mu,GRU1Mu and GRC1Mu)
#=======================================================================
# Importing SE of EBVs for each model

SE_G1Mu <- read.csv(file = "SE_G1Mu.csv", stringsAsFactors = FALSE)
#SE_G2Mu <- read.csv(file = "SE_G2Mu.csv", stringsAsFactors = FALSE)
SE_GR1Mu <- read.csv(file = "SE_GR1Mu.csv", stringsAsFactors = FALSE)
#SE_GR2Mu <- read.csv(file = "SE_GR2Mu.csv", stringsAsFactors = FALSE)
SE_GRU1Mu <- read.csv(file = "SE_GRU1Mu.csv", stringsAsFactors = FALSE)
#SE_GRU2Mu <- read.csv(file = "SE_GRU2Mu.csv", stringsAsFactors = FALSE)
SE_GRC1Mu <- read.csv(file = "SE_GRC1Mu.csv", stringsAsFactors = FALSE)
#SE_GRC2Mu <- read.csv(file = "SE_GRC2Mu.csv", stringsAsFactors = FALSE)
#SE_GRC1_NS <- read.csv(file = "SE_GRC1_NSMu.csv", stringsAsFactors = FALSE)
#SE_GRC2_NS <- read.csv(file = "SE_GRC2_NSMu.csv", stringsAsFactors = FALSE)
#=======================================================================
# PEVMu Data organization (Municipality)
#=======================================================================
##SE_G1Mu_anim
SE_G1Mu$factor_G1<- as.factor(SE_G1Mu$factor_G1)
SE_G1Mu$level_G1<- as.factor(SE_G1Mu$level_G1)
str(SE_G1Mu)
SE_G1Mu_anim <- SE_G1Mu %>% 
  filter(factor_G1 == 4)
summary(SE_G1Mu_anim)

## SE_GR1Mu_anim
SE_GR1Mu$factor_GR1<- as.factor(SE_GR1Mu$factor_GR1)
SE_GR1Mu$level_GR1<- as.factor(SE_GR1Mu$level_GR1)
str(SE_GR1Mu)
summary(SE_GR1Mu$factor_GR1)
SE_GR1Mu_anim <- SE_GR1Mu %>% 
  filter(factor_GR1 == 5)
summary(SE_GR1Mu_anim)  

##SE_GRU1Mu_anim
SE_GRU1Mu$factor_GRU1<- as.factor(SE_GRU1Mu$factor_GRU1)
SE_GRU1$level_GRU1<- as.factor(SE_GRU1Mu$level_GRU1)
str(SE_GRU1Mu)
summary(SE_GRU1Mu$factor_GRU1)
SE_GRU1Mu_anim <- SE_GRU1Mu %>% 
  filter(factor_GRU1 == 5)
summary(SE_GRU1Mu_anim) 

## SE_GRC1_anim
SE_GRC1Mu$factor_GRC1<- as.factor(SE_GRC1Mu$factor_GRC1)
SE_GRC1Mu$level_GRC1<- as.factor(SE_GRC1Mu$level_GRC1)
str(SE_GRC1Mu)
summary(SE_GRC1Mu$factor_GRC1)
SE_GRC1Mu_anim <- SE_GRC1Mu %>% 
  filter(factor_GRC1 == 4)
summary(SE_GRC1Mu_anim) 

SE_EBVMu_anim<- cbind(SE_G1Mu_anim,SE_GR1Mu_anim, SE_GRU1Mu_anim, SE_GRC1Mu_anim) # EBVs and SE for models with municipality

write.csv(SE_EBVMu_anim,"/Users/ihouaga/Documents/Project1/SE_EBVMu_anim.csv") ## Save EBVs and SE for models with municipality 

SEMu<- SE_EBVMu_anim %>% select(level_G1, SE_G1Mu, SE_GR1Mu,SE_GRU1Mu,SE_GRC1Mu)
colnames(SEMu)[1] <- "renum_ID"  

## Calculate PEV
PEVMu<- SEMu  %>%
  mutate(PEV_G1Mu=SE_G1Mu^2, PEV_GR1Mu=SE_GR1Mu^2,PEV_GRU1Mu=SE_GRU1Mu^2,PEV_GRC1Mu=SE_GRC1Mu^2)

write.csv(PEVMu,"/Users/ihouaga/Documents/Project1/PEV.csv") #Save EBvs

summary(PEVMu)
### Municipality Estimates
E_Mu<- SE_GRC1Mu %>% 
filter(factor_GRC1 == 6)
colnames(E_Mu)[3] <- "E_Mu"
## Region solutions in GRC1Mu follow same as region_ID in neighbours_scaledMu.txt
Mu_Scaled_ID<- read.csv(file = "Mu_Scaled_ID.csv", stringsAsFactors = FALSE)

Mu_Scaled_ID<- Mu_Scaled_ID %>% distinct()

Mu_Scaled_ID$level_GRC1= 1:213 
Mu_Scaled_ID <- Mu_Scaled_ID %>% select(level_GRC1,Mu_ID)
E_Mu<- E_Mu %>% select(level_GRC1,E_Mu)
E_Mu_final<- merge(Mu_Scaled_ID,E_Mu, by="level_GRC1", all.x= TRUE)

E_Mu_final<- E_Mu_final %>% select(Mu_ID,E_Mu)


## Colour the map with region estimates from besag1

colnames(E_Mu_final)[1] <- "OBJECTID"
(.packages())
## Reading mapMu2 South Africa for ploting as sf . (municipalities)
mapMu2 <- st_read("/Users/ihouaga/Documents/Project1/1cecc5ab-304d-494a-a29d-1b33611768702020329-1-ai5f0p.uax88.shp", quiet = TRUE)          
class(mapMu2)
mapMu2_estimates<- merge(mapMu2,E_Mu_final,duplicateGeoms = T)
head(mapMu2_estimates)
E_Mu_final$E_Mu<- as.numeric(E_Mu_final$E_Mu)
summary(E_Mu_final)

## Plot one with scale_fill_gradients
ggplot(data = mapMu2_estimates, aes(fill = E_Mu)) + geom_sf() +
  scale_fill_gradient(low = "blue", high = "red") + 
  labs(fill='Municipality estimates') 

## Plot two without scale_fill_gradients
ggplot(data = mapMu2_estimates, aes(fill = E_Mu)) + geom_sf() +
  labs(fill='Municipality estimates') 

#### Scripts below need to be updated
library(devtools)
install_github("trimped")

(.packages())
#================================================================================================
# Grouping EBV by region
#================================================================================================ 
EBVs_region<- EBVs_region_ID %>% 
  group_by(region) %>%  
  summarise(mean_EBV_GRC1=mean(EBV_GRC1),sd_EBV_GRC1=sd(EBV_GRC1), mean_EBV_GRC2=mean(EBV_GRC2),sd_EBV_GRC2=sd(EBV_GRC2))
write.csv(EBVs_region,"/Users/ihouaga/Documents/Project1/EBVs_region.csv") #Save EBVs_region

#================================================================================================
# Grouping PEV by region
#================================================================================================  
EBVs_region_ID_2<- EBVs_region_ID %>% relocate(renum_ID)
EBVs_PEV_region <- merge(EBVs_region_ID_2,PEV, by="renum_ID", all.x=TRUE)  

PEV_region<- EBVs_PEV_region %>% 
  group_by(region)%>%
  summarise(mean_PEV_GRC1=mean(PEV_GRC1),mean_PEV_GRC2=mean(PEV_GRC2))
write.csv(PEV_region,"/Users/ihouaga/Documents/Project1/PEV_region.csv")

# Scatterplot (x-y) between sd of EBV and average PEV per region 
EBVs_PEV_region_sd_mean <- merge(EBVs_region, PEV_region, by="region", all.x=TRUE)

ggplot(EBVs_PEV_region_sd_mean, aes(x = sd_EBV_GRC1, y = mean_PEV_GRC1, colour=region)) +
  geom_point()
#================================================================================================
# Grouping PEV by SIRE, COW_ID and non phenotype ID
#================================================================================================  
## Import PEV and data_clean_final
PEV<- read.csv(file = "PEV.csv", stringsAsFactors = FALSE)
data_clean_final<- read.csv(file = "data_clean_final.csv", stringsAsFactors = FALSE)
## Import animal IDs (renum+Original)
ani_ID<- read.csv(file = "renum_ori_ID.csv", stringsAsFactors = FALSE)
PEV_ani_ID<- merge(PEV,ani_ID, by="renum_ID", all.x=TRUE)

PEV_ani_ID<- PEV_ani_ID %>% select(ori_anim_ID, everything())
## Sire ID
sire_ID <- data_clean_final %>%
  select(SIRE_ID)
sire_ID_uniq<-  sire_ID %>% distinct()  # 5978 Sire

sire_ID_uniq<- sire_ID_uniq %>%
  filter(SIRE_ID >0) ## Remove missing sire

COW_ID<- data_clean_final %>%
  select(COW_ID)
COW_ID_uniq<-  COW_ID %>% distinct()

ani_ID_ori <- ani_ID %>%
  select(ori_anim_ID) 
# Save COW_ID_uniqu and sire_ID_uniq

write.csv(COW_ID_uniq,"/Users/ihouaga/Documents/Project1/COW_ID_uniq.csv")
write.csv(sire_ID_uniq,"/Users/ihouaga/Documents/Project1/sire_ID_uniq.csv")

#rename Colname in COW_ID_unique 

colnames(COW_ID_uniq)[1] <- "ori_anim_ID"
COW_ID_uniq
COW_ID_uniq$mean <- rep(1, nrow(COW_ID_uniq))
Non_pheno_anim_ID<- merge (ani_ID_ori,COW_ID_uniq, by="ori_anim_ID", all.x=TRUE) ## NA= non pheno animal
Non_pheno_anim_ID<- Non_pheno_anim_ID %>% 
  filter(is.na(mean)) %>% 
  select(-c(mean))
Non_pheno_anim_ID## 98146
colnames(sire_ID_uniq)[1] <- "ori_anim_ID" 
sire_ID_uniq$mean <- rep(1, nrow(sire_ID_uniq))
Non_pheno_COW_ID<- merge (Non_pheno_anim_ID,sire_ID_uniq, by="ori_anim_ID", all.x=TRUE) ## NA= non pheno animal

Non_pheno_COW_ID<- Non_pheno_COW_ID  %>% 
  filter(is.na(mean)) %>% 
  select(-c(mean))
Non_pheno_COW_ID# 92169 non phenotyped cows

#Phenotyped cows=98632
# Sire: 5977
# 5977+98632+92169= 196778
# Save non phenotyped animal IDs  
write.csv(Non_pheno_anim_ID,"/Users/ihouaga/Documents/Project1/Non_pheno_anim_ID.csv")

## PEV sire
PEV_sire<- merge(sire_ID_uniq, PEV_ani_ID, by="ori_anim_ID", all.x=TRUE)
#PEV_pheno_cow
PEV_pheno_cow<- merge(COW_ID_uniq, PEV_ani_ID, by="ori_anim_ID", all.x=TRUE)

#PEV_non_pheno_anim 
PEV_non_pheno_anim<- merge(Non_pheno_anim_ID, PEV_ani_ID, by="ori_anim_ID", all.x=TRUE)

#PEV_non_pheno_COW 
PEV_non_pheno_cow<- merge(Non_pheno_COW_ID, PEV_ani_ID, by="ori_anim_ID", all.x=TRUE)

# Correlation Sire_PEV
cor.test(PEV_sire$PEV_G1,PEV_sire$PEV_GR1)
cor.test(PEV_sire$PEV_G1,PEV_sire$PEV_GRU1)
cor.test(PEV_sire$PEV_G1,PEV_sire$PEV_GRC1)
cor.test(PEV_sire$PEV_GR1,PEV_sire$PEV_GRU1)
cor.test(PEV_sire$PEV_GR1,PEV_sire$PEV_GRC1)
cor.test(PEV_sire$PEV_GRU1,PEV_sire$PEV_GRC1)

#Correlation PEV_pheno_cow
cor.test(PEV_pheno_cow$PEV_G1,PEV_pheno_cow$PEV_GR1)
cor.test(PEV_pheno_cow$PEV_G1,PEV_pheno_cow$PEV_GRU1)
cor.test(PEV_pheno_cow$PEV_G1,PEV_pheno_cow$PEV_GRC1)
cor.test(PEV_pheno_cow$PEV_GR1,PEV_pheno_cow$PEV_GRU1)
cor.test(PEV_pheno_cow$PEV_GR1,PEV_pheno_cow$PEV_GRC1)
cor.test(PEV_pheno_cow$PEV_GRU1,PEV_pheno_cow$PEV_GRC1)

#Correlation PEV_non_pheno_anim
cor.test(PEV_non_pheno_anim$PEV_G1,PEV_non_pheno_anim$PEV_GR1)
cor.test(PEV_non_pheno_anim$PEV_G1,PEV_non_pheno_anim$PEV_GRU1)
cor.test(PEV_non_pheno_anim$PEV_G1,PEV_non_pheno_anim$PEV_GRC1)
cor.test(PEV_non_pheno_anim$PEV_GR1,PEV_non_pheno_anim$PEV_GRU1)
cor.test(PEV_non_pheno_anim$PEV_GR1,PEV_non_pheno_anim$PEV_GRC1)
cor.test(PEV_non_pheno_anim$PEV_GRU1,PEV_non_pheno_anim$PEV_GRC1)

#Correlation PEV_non_pheno_anim
cor.test(PEV_non_pheno_cow$PEV_G1,PEV_non_pheno_cow$PEV_GR1)
cor.test(PEV_non_pheno_cow$PEV_G1,PEV_non_pheno_cow$PEV_GRU1)
cor.test(PEV_non_pheno_cow$PEV_G1,PEV_non_pheno_cow$PEV_GRC1)
cor.test(PEV_non_pheno_cow$PEV_GR1,PEV_non_pheno_cow$PEV_GRU1)
cor.test(PEV_non_pheno_cow$PEV_GR1,PEV_non_pheno_cow$PEV_GRC1)
cor.test(PEV_non_pheno_cow$PEV_GRU1,PEV_non_pheno_cow$PEV_GRC1)

# HYS models Sire_PEV
cor.test(PEV_sire$PEV_G2,PEV_sire$PEV_GR2)
cor.test(PEV_sire$PEV_G2,PEV_sire$PEV_GRU2)
cor.test(PEV_sire$PEV_G2,PEV_sire$PEV_GRC2)
cor.test(PEV_sire$PEV_GR2,PEV_sire$PEV_GRU2)
cor.test(PEV_sire$PEV_GR2,PEV_sire$PEV_GRC2)
cor.test(PEV_sire$PEV_GRU2,PEV_sire$PEV_GRC2)

#================================================================================================
# Reliability and accuracy  SIRE_EBV
#================================================================================================  
summary(PEV_sire$PEV_G1)
summary(PEV_sire$PEV_GR1)
summary(PEV_sire$PEV_GRU1)
summary(PEV_sire$PEV_GRC1)
summary(PEV_sire$PEV_GRC1_NS)


summary(PEV_sire$PEV_G2)
summary(PEV_sire$PEV_GR2)
summary(PEV_sire$PEV_GRU2)
summary(PEV_sire$PEV_GRC2)
summary(PEV_sire$PEV_GRC2_NS)





