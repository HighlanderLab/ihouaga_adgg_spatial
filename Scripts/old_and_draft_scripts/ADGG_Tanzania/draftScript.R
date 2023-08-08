#Drafted Script

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





#Yield deviation model G ie YDG
YDG <- milk ~ 1 + cyrsn + tyrmn + dgrp + (age|lacgr) + (intercept|lacgr) + (leg1|lacgr) +(leg2|lacgr) + f(herd, model = 'iid') + f(cowPe, model = 'iid') 

modelYDG<- inla(formula = YDG,  control.compute = list(dic = TRUE), data=pheno)
YDG <- data.frame(modelYDG$summary.random$cowPe)

# Predict EBVs for cow with breed composition dgrp==1 (Lets set milk=NA for dgrp==1)

table(data6$dgrp) # dgrp1 (7436 records),  dgrp2 (8114 records), dgrp3 (3058records) and dgrp4(702 records)

# Mask phenotype for dgrp1

data_dgrp2_3_4 <- data6 %>% mutate(milk = ifelse(dgrp == 1, NA, milk))
sum(is.na(data_dgrp2_3_4))# 7436 records
modelG_dgrp_1 <- inla(formula = G,  control.compute = list(dic = TRUE), data=data_dgrp2_3_4) 
PGBV2_3_4 <-data.frame(modelG_dgrp_1$summary.random$INLA_ID)
PGBV2_3_4 <- subset(PGBV2_3_4, select= c(ID,mean))
head(PGBV2_3_4)
names(PGBV2_3_4)[1]<- "INLA_cowID"
names(PGBV2_3_4)[2]<- "PGBV"
# Cow_ID_dgrp1
datadgrp1 <- data6 %>% subset(dgrp==1)
Cow_ID_dgrp1 <- unique(datadgrp1$INLA_ID)
length(Cow_ID_dgrp1)# 728 cows in dgrp1
Cow_ID_dgrp1<- as.data.frame(Cow_ID_dgrp1)
names(Cow_ID_dgrp1)[1]<- "INLA_cowID"

#YDG corrected Phenotype 

YDG <- subset(YDG, select= c(ID,mean)) 
names(YDG)[1]<- "INLA_cowID"
names(YDG)[2]<- "PhenoCor"

## Merging PGBV2_3_4, YDG and Cow_ID_dgrp1 by "INLA_cowID
# Merging Cow_ID_dgrp1, YDG
PhenoCor1G <- merge(x=Cow_ID_dgrp1,y= YDG, by="INLA_cowID",all.x=TRUE)
length(PhenoCor1G$INLA_cowID)
EBV1G <- merge(x=PhenoCor1G,y= PGBV2_3_4, by="INLA_cowID",all.x=TRUE)
length(EBV1G$INLA_cowID)
#Option 1: 
Accuracy1.1G = cor(EBV1G$PhenoCor,EBV1G$PGBV)
Accuracy1.1G # -0.01185069

#Option 2
Accuracy1.2G = cor(EBV1G$PhenoCor,EBV1G$PGBV)/sqrt(0.03743606)
Accuracy1.2G #


##### Mask phenotype for dgrp2

data_dgrp1_3_4 <- data6 %>% mutate(milk = ifelse(dgrp == 2, NA, milk))
sum(is.na(data_dgrp1_3_4))# 8114 records
modelG_dgrp_2 <- inla(formula = G,  control.compute = list(dic = TRUE), data=data_dgrp1_3_4) 
PGBV1_3_4 <-data.frame(modelG_dgrp_2$summary.random$INLA_ID)
PGBV1_3_4 <- subset(PGBV1_3_4, select= c(ID,mean))
head(PGBV1_3_4)
names(PGBV1_3_4)[1]<- "INLA_cowID"
names(PGBV1_3_4)[2]<- "PGBV"
# Cow_ID_dgrp2
datadgrp2 <- data6 %>% subset(dgrp==2)
Cow_ID_dgrp2 <- unique(datadgrp2$INLA_ID)
length(Cow_ID_dgrp2)#  cows in dgrp2
Cow_ID_dgrp2<- as.data.frame(Cow_ID_dgrp2)
names(Cow_ID_dgrp2)[1]<- "INLA_cowID"
