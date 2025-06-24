# ---- Header ------------------------------------------------------------------
# Spatial modelling
# Trait: Test day milk yield
# Author: Isidore Houaga, Ivan Pocrnic and Gregor Gorjanc.
# Version 1.0.0
# Date: 2023-07-20
# ---- Setup--------------------------------------------------------------------

# Working directory
baseDir <- "/Users/ihouaga2/ihouaga_adgg_spatial"
# Change working directory
setwd(dir = baseDir)
getwd()
# ---- Installing and loading packages--------------------------------------------
#Load required libraries
library(tidyverse)
library(INLA)
library(inlabru)
library(irlba)
library(Hmisc) # Correlation matrix with P-values
library(easyGgplot2)
library(psych) # scatter-plot matrix
(.packages()) # Check loaded packages
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
# Figure 3 Manuscript
pca <-  ggarrange(dgrp_pca, region_pca, ncol = 1, nrow = 2, common.legend = T, legend = "right", align = "h", widths = c(1,1,1))
pca
ggsave(plot = pca + PreseTheme, filename = "pca_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot =  pca + PaperTheme, filename = "pca_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper


#-----------------------Manuscript-EBVs Correlations with model components and Animal Ranking------------------------------------------
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

# Spearman's Rank Correlations
round(cor(EBVs_Spatial_Matrix, method = "spearman"),2)
#Adding difference EBV_GH - EBV_GHS
EBVs_Spatial_effect <- EBVs_Spatial_effect %>% mutate(dEBV= EBV_GH - EBV_GHS)

round(cor(EBVs_Spatial_effect$dEBV,EBVs_Spatial_effect$Spatial_effect_GHS),2)

#--------Multivariate analysis of variances (MANOVA)----------------------------------------------------------------
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
#Figure 5 manuscript
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
