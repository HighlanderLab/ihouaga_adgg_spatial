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
library(Hmisc) 
library(easyGgplot2)
library(psych) 
(.packages()) # Check loaded packages
#-----------------------Principal Components analysis---------------------------
# Import SNP data (QC done)
geno <- read.table(file = "data/original_data/snpref-isi.txt", header = FALSE)

geno[1:10, 1:10]
dim(geno) # 1911 664823
colnames(geno)
summary(geno$V1) 
colnames(geno)[1] <- "cow"
library(package = "irlba")

pcasAll <- prcomp_irlba(x = geno[, -1], n = 3)
summary(pcasAll)
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

# Visualize PCA by the proportion of exotic genes

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

# Visualize PCA by region

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

region_pca <- ggarrange(region_pca1,region_pca2,region_pca3, ncol = 3, nrow = 1, common.legend = T, legend = "right", align = "h", widths = c(1,1,1))
region_pca

ggsave(plot = region_pca + PreseTheme, filename = "alpha0.1PCA_region_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") #Presentation


ggsave(plot = region_pca + PaperTheme, filename = "PCA_region_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") # Paper

#-----Final Plot paper----------------------------------------------------------
# Figure 3 Manuscript
pca <-  ggarrange(dgrp_pca, region_pca, ncol = 1, nrow = 2, common.legend = T, legend = "right", align = "h", widths = c(1,1,1))

ggsave(plot = pca + PreseTheme, filename = "pca_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") 

ggsave(plot =  pca + PaperTheme, filename = "pca_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") 

#-----------------------Manuscript-EBVs Correlations with model components and Animal Ranking------------------------------------------
# fitG
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
EBVs_Spatial_Matrix <- EBVs_Spatial_effect[, -c(1)]

round(cor(EBVs_Spatial_effect$EBV_GH,EBVs_Spatial_effect$EBV_GHS),2) 

str(EBVs_Spatial_Matrix)
# Pearson's Correlation
round(cor(EBVs_Spatial_Matrix, method = "pearson"),2)

# Spearman's Rank Correlations
round(cor(EBVs_Spatial_Matrix, method = "spearman"),2)

# Difference EBV_GH - EBV_GHS
EBVs_Spatial_effect <- EBVs_Spatial_effect %>% mutate(dEBV= EBV_GH - EBV_GHS)

round(cor(EBVs_Spatial_effect$dEBV,EBVs_Spatial_effect$Spatial_effect_GHS),2)

#--------Multivariate analysis of variances (MANOVA)----------------------------------------------------------------
with(pcas, plot(PC1 ~ PC2, col = c("blue", "red", "black", "yellow")[dgrp]))
fit <- manova(cbind(PC1, PC2, PC3) ~ dgrp, data = pcas)
summary(fit) 
summary.aov(fit)

fit <- lm(PC1 ~ dgrp, data = pcas)
summary(fit)

fit <- lm(PC2 ~ dgrp, data = pcas)
summary(fit)

fit <- lm(PC3 ~ dgrp, data = pcas)
summary(fit)

# MANOVA Region
with(pcas, plot(PC1 ~ PC2, col = c("blue", "red", "black", "yellow")[region]))
fit <- manova(cbind(PC1, PC2, PC3) ~ region, data = pcas)
summary(fit) 
summary.aov(fit)

fit <- lm(PC1 ~ region, data = pcas)
summary(fit)

plot(y = data1$lat, x = data1$long, col = c("blue", "red", "black", "yellow")[data1$region])
plot(y = data2$lat, x = data2$long)
nrow(pcas)

fit <- lm(PC2 ~ region, data = pcas)
summary(fit)

fit <- lm(PC3 ~ region, data = pcas)
summary(fit)

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
summary(fit) 

summary.aov(fit)

fit <- lm(PC1 ~ dgrp * region, data = pcas)
summary(fit)

fit <- lm(PC2 ~ dgrp * region, data = pcas)
summary(fit)

fit <- lm(PC3 ~ dgrp * region, data = pcas)
summary(fit)

with(tmp, table(dgrp, region))

tmp2 <- prop.table(with(tmp, table(dgrp, region)))*100
tmp2

#----------------Scatter plot dEBV (EBV_GH-EBV_GHS) vs Spatial effect GHS------------------
# Figure 5 manuscript

# Fit the linear model
model <- lm(dEBV ~ Spatial_effect_GHS, data = EBVs_Spatial_effect)

# Extract R-squared value
r_squared <- summary(model)$r.squared

# Create the scatter plot with regression line and R-squared annotation

scatter_plot <- ggplot(EBVs_Spatial_effect, aes(x = Spatial_effect_GHS, y = dEBV)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +  
  labs(x = "Spatial effects (GHS)",
       y = "Difference between EBVs (GH-GHS)",
       color = "Spatial effects") +
  theme_bw() +
  annotate("text", x = Inf, y = Inf,
           label = paste("R² =", round(r_squared, 2)),
           hjust = 1.1, vjust = 2, size = 5, color = "black") +
  theme(
    axis.text.x = element_text(size = 15.3),  
    axis.text.y = element_text(size = 15.3),
    axis.title.x = element_text(size = 15.3),  
    axis.title.y = element_text(size = 15.3)  
  ) +
  ylim(-0.5, 0.5)

print(scatter_plot)

ggsave(plot = scatter_plot, filename = "2_spatial_dEBV_scatter_paper.png",
       height = PaperSize, width = PaperSize * 1.5, unit = "cm") 

ggsave(plot = scatter_plot, filename = "2_spatial_dEBV_scatter_presentation.png",
       height = PreseSize, width = PreseSize * 1.5, unit = "cm") 

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

# Top 20 EBV_GH
EBV_GH_top20 <- EBV_GH_Rank[1:20,c("ID","EBV_GH")]
head(EBV_GH_top20, n=20)

# Top 20 EBV_GHS
EBV_GHS_top20 <- EBV_GHS_Rank[1:20,c("ID","EBV_GHS")]
head(EBV_GHS_top20, n=20)

sel_top20 <- EBV_GH_top20$ID %in% EBV_GHS_top20$ID
table(sel_top20)

# Top 50 EBV_GH
EBV_GH_top50 <- EBV_GH_Rank[1:50,c("ID","EBV_GH")]

# Top 50 EBV_GHS
EBV_GHS_top50 <- EBV_GHS_Rank[1:50,c("ID","EBV_GHS")]

sel_top50 <- EBV_GH_top50$ID %in% EBV_GHS_top50$ID
table(sel_top50)

# Top 100 EBV_GH
EBV_GH_top100 <- EBV_GH_Rank[1:100,c("ID","EBV_GH")]

# Top 100 EBV_GHS
EBV_GHS_top100 <- EBV_GHS_Rank[1:100,c("ID","EBV_GHS")]

sel_top100 <- EBV_GH_top100$ID %in% EBV_GHS_top100$ID
table(sel_top100)

# Top 1000 EBV_GH
EBV_GH_top1000 <- EBV_GH_Rank[1:1000,c("ID","EBV_GH")]

# Top 1000 EBV_GHS
EBV_GHS_top1000 <- EBV_GHS_Rank[1:1000,c("ID","EBV_GHS")]

sel_top1000 <- EBV_GH_top1000$ID %in% EBV_GHS_top1000$ID
table(sel_top1000)
