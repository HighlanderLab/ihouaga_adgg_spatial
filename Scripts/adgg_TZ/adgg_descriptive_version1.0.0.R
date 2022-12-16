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
#install.packages(c("dplyr","stringr","ggplot2","ggExtra","EnvStats","lme4","afex","readr","GGally","sf","sp","INLA","spData","spdep","tidyr","rgdal","terra","devtools","cowplot"))
library(dplyr)
#library(stringr)
library(ggplot2) 
library("orthopolynom", lib.loc = "/home/ihouaga/Rpackages")# Legendre polynomial
#library(ggExtra)
#library(EnvStats)
#library(lme4)     #for mixed model
#library(afex)     #for p-values
#library(readr)
#library(GGally)
#library(sf) # import shapefile
#library(sp) # poly2nb for graph
#library(INLA) # inla.read.graph
#library(spData)
#library(spdep)
#library(tidyr)
#library(rgdal)
#library(terra)
#library(devtools)
#library(cowplot)
(.packages()) # Checking loaded packages
#=======================================================================
# Uploading datasets
#=======================================================================
# Importing ADGG data
rm(list = ls() ) # Clear the environment
data1 <- read.table(file = "/home/rmrode/data/milk4d-gps-pre.txt", header = FALSE)# path on ILRI hpc
colnames(data1) <- c("cow","ward","herd","cyrsn","tyrmn","dgrp","lacest","lac","dim","htd","hcyr","htyr","pym","ksea","bprop","age","milk","hgirth","bodywt","bcs","long","lat")
head(data1)
#=======================================================================
# Variables definition
#=======================================================================
#cow
#ward
#herd
#cyrsn: calving-year-season
#tyrmn-test-year-month
#dgrp-breed class (4)
#lacest--lactation 1  and 2+later
#lac:lactations (uncombined)
#dim: days in milk
#htd: herd-test-day
#hcyr:herd calving-year
#htyr: herd-test-year
#pym: milk-year-month
#ksea: calving season
#bprop: breed proportion (exotic)
#age (month)
#milk
#hgirth: hearth-girth
#bodywt: body weight
#bcs:Body condition Score
#long: Longitude
#lat: Latitude
#=======================================================================
# Data organization
#=======================================================================
data2 <- data1 %>%
  mutate_at(c("cow", "ward","herd", "cyrsn","tyrmn","dgrp","lacest","lac","htd","hcyr","htyr","pym","ksea","bprop"),factor)
data2$dim <- as.numeric(data2$dim)
data2$age <- as.numeric(data2$age)
data2$milk <- as.numeric(data2$milk)
data2$hgirth <- as.numeric(data2$hgirth)
data2$bodywt <- as.numeric(data2$bodywt)
##=======================================================================
# Check for duplicates and NA
#=======================================================================
#Remove potential duplicated rows
data <- dplyr::distinct(data2) ###Efficient duplicate removal
# Check for missing data (NA)
sum(is.na(data))
#=======================================================================
# Summary 
#=======================================================================
#Summary data
summary(data)
# number of records 
nrecords<- length(data$cow)
nrecords
#number of cows
ncows<- length(unique(data$cow))
ncows
#number of herds
nherds<- length(unique(data$herd))
nherds
# ward
nwards<- length(unique(data$ward))
nwards

#number of herds per ward
nherdward <- data %>%
  group_by(ward, herd)%>%
  summarise(nherd_byward= n())
nherdward_2<- nherdward %>%
  select(ward) %>%
  summarise(nherdinward=n())
summary(nherdward_2$nherdinward)

# Number of records per cow per herd
ncowherd<- data %>%
  group_by(herd,cow)%>%
  summarise(ncow_byherd= n())
summary(ncowherd$ncow_byherd)
# number of cows per herd
ncowherd_2<- ncowherd %>%
  select(herd) %>%
  summarise(ncowherd=n())
summary(ncowherd_2$ncowherd)

# number of records per herd
nrecords_byherd<- data %>%
  group_by(herd) %>%
  summarise(recordherd=n()) 
summary(nrecords_byherd$recordherd)

#number of records per ward
nrecords_byward <- data %>%
  group_by(ward) %>%
  summarise(nrecord=n())
summary(nrecords_byward$nrecord)

# Number of records per cow per ward
ncowward<- data %>%
  group_by(ward,cow)%>%
  summarise(ncow_byward= n())
summary(ncowward$ncow_byward)
# number of cows per ward
ncowward_2<- ncowward %>%
  select(ward) %>%
  summarise(ncowward=n())
summary(ncowward_2$ncowward)
#-----------------------------------------------------------------------
# Number of cows in Parity/Summary age   
#-----------------------------------------------------------------------
#Number of parity
summary(data$lac)

#Number of cows in Parity 1
p1<- subset(data, data$lac=="1")
ncows1<- length(unique(p1$cow))
ncows1

#Number of cows in Parity 2
p2<- subset(data, data$lac=="2")
ncows2<- length(unique(p2$cow))
ncows2

#Number of cows in Parity 3
p3<- subset(data, data$lac=="3")
ncows3<- length(unique(p3$cow))
ncows3

#Number of cows in Parity 4
p4<- subset(data, data$lac=="4")
ncows4<- length(unique(p4$cow))
ncows4

#Number of cows in Parity 5
p5<- subset(data, data$lac=="5")
ncows5<- length(unique(p5$cow))
ncows5

#Number of cows in Parity 6
p6<- subset(data, data$lac=="6")
ncows6<- length(unique(p6$cow))
ncows6

#Number of cows in Parity 7
p7<- subset(data, data$lac=="7")
ncows7<- length(unique(p7$cow))
ncows7

#Number of cows in Parity 8
p8<- subset(data, data$lac=="8")
ncows8<- length(unique(p8$cow))
ncows8

#Number of cows in Parity 9
p9<- subset(data, data$lac=="9")
ncows9<- length(unique(p9$cow))
ncows9

# Number of cows with records in grouped parities
summary(data$lacest)
# Number of cows with records in Parity 1
p1_1<- subset(data, data$lacest=="1")
ncows1_1<- length(unique(p1_1$cow))
ncows1_1
ncows1==ncows1_1 # TRUE expected
# Number of Cows with records in Parity 2+
p2_2<- subset(data, data$lacest=="2")
ncows2_2<- length(unique(p2_2$cow))
ncows2_2
# Summary age
summary(data$age)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#18.05   38.79   47.83   50.28   58.88  162.10 
# Summary age at first calving (AFC)
summary(p1_1$age)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#18.05   34.19   39.91   42.09   47.60  117.19
#=======================================================================
# Descriptive statistics of test-day milk yield across 9 lactations
#=======================================================================
#write.csv(dataHST, "/Users/ihouaga/Documents/Project1/dataHSTRaw.csv") ## save .csv copy of raw data
# My Toy data to test scripts
#cow<- c(1,1,1,1,1,2,2,2,3,3,4,4,5,5,5)
#herd<- c(1,1,1,1,1,3,3,3,2,2,3,3, 2,2,2)
#ward <- c(1,1,1,1,1,2,2,2,2,2,2,2,2,2,2)
#lac <- c(1,2,3,4,5,1,2,3,1,2,1,2,1,2,3)
#milk<- c(5,6,7,8,9,4,5,6,3,4,5,6,5,6,7)
#toy<- data.frame(cow,herd,ward,lac,milk)
#toy
#toy$milk <- as.numeric(toy$milk)
#toy$cow <- as.factor(toy$cow)
#toy$herd <- as.factor(toy$herd)
#toy$ward <- as.factor(toy$ward)
#toy$lac <- as.factor(toy$lac)
#-----------------------------------------------------------------------
#Number of levels of calving season,calving year season, test-year-month 
#and breed class (breed proportion of exotic)  
#------------------------------------------------------------------------
#number of levels of calving seasons 
nlevels(data$ksea)
#number of levels of calving year season
nlevels(data$cyrsn)
#levels of test-year-month
nlevels(data$tyrmn)
# number of levels of breed proportion (exotic)
nlevels(data$dgrp)
#=======================================================================
# Descriptive Statistics
#=======================================================================
summary(data1)
#=======================================================================
# Means, sd, histogram
#=======================================================================
#Milk test-day milk yield across 9 lactations 
min(data1$milk)
max(data1$milk)
mean(data1$milk,na.rm = TRUE)
sd(data1$milk,na.rm = TRUE)
hist(na.omit(data1$milk))
#=======================================================================
# Mean and variance
#=======================================================================
# Mean and Variance - grouped by parity (lactation)
SMilk <- data1%>%
  group_by(lac) %>%
  summarise(
    MeanMilk = mean(milk, na.rm = TRUE),
    MedianMilk = median(milk, na.rm = TRUE),
    VarMilk = var(milk, na.rm = TRUE),
    SdMilk =  sqrt(VarMilk),
    nMilk = n()
  )
SMilk
#-----------------------------------------------------------------------
# Mean vs. Variance
SMilk %>%
  ggplot(aes(y = MeanMilk, x = SdMilk,  colour = lac)) +
  geom_text(aes(label = lac)) +
  xlab("SD(Milk)") +
  ylab("E(Milk)") +
  theme_bw()
# Linear relationship between mean and variance????????
# Mean vs. Variance grouped by herd and parity
SMilk %>%
  ggplot(aes(y = MeanMilk, x = SdMilk,  label = lac)) +
  geom_label(aes(fill = factor(lac)), colour = "white",
             fontface = "bold") +
  xlab("SD(Milk)") +
  ylab("E(Milk)") +
  theme_bw()
ggsave("parity.png")
#-----------------------------------------------------------------------
# Mean vs. number of observations
SMilk %>%
  ggplot(aes(y = MeanMilk, x = nMilk,  label = lac)) +
  geom_label(aes(fill = factor(lac)), colour = "white",
             fontface = "bold") +
  xlab("Number of Observations") +
  ylab("E(Milk)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# We can see the milk production increases up to parity xxxx 3-4 and then
# starts to decrease. This is related to the combination of effect of
# parity and number of observations????

# SD vs. number of observations
SMilk %>%
  ggplot(aes(y = SdMilk, x = nMilk,  label = lac)) +
  geom_label(aes(fill = factor(lac)), colour = "white",
             fontface = "bold") +
  xlab("Number of Observations") +
  ylab("SD(Milk)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Same for Standard deviation???????
#-----------------------------------------------------------------------
# Mean and Median relationship
SMilk %>%
  ggplot(aes(y = MeanMilk,  x = MedianMilk)) +
  geom_abline(intercept =0, slope=1) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("Median(Milk)") +
  ylab("E(Milk)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Ok -> linear relationship?????
#=======================================================================
# Centering the variable Milk
#=======================================================================
data1 %>%
  mutate_at("lac", factor) %>%
  mutate(MYc = milk-mean(milk)) %>%
  ggplot(aes(y = MYc, x = lac)) +
  geom_boxplot() +
  xlab("Levels of Parity") +
  ylab("Test-day Milk Yield") +
  theme_bw()
# Distribution is assymetric for test day milk yield??????
# data<- read.csv(file = "testdata.csv", stringsAsFactors = FALSE)

#=======================================================================
# #Data preparation for Mrode et al. 2021 GBLUP in blupf90
#=======================================================================
#Create column for grouped parity (1,2 and 3+)
data$lac <- as.numeric(data$lac)
data3 <- data %>% mutate(lacgr=lac) %>% 
    mutate(lacgr = ifelse(lacgr >= "3", "3", lacgr))
head(data3)
data3$lacgr <- as.factor(data3$lacgr)
summary(data3$lacgr)
#Create columns for intercept and legendre polynomial

m1 = legendre.polynomials(n=2, normalized=TRUE)
# u + ((v - u) * ((x - min(x))/(max(x) - min(x)))
# if u=-1 and v=1 than -1 + 2*((x - min(x))/(max(x) - min(x)))
dataLegendre <- as.matrix(as.data.frame(polynomial.values(polynomials=m1,
                                            x=scaleX(data3$dim, u=-1, v=1))))
colnames(dataLegendre) <- c("intercept", "leg1", "leg2")
data4 <- cbind(data3,dataLegendre)
head(data4)
# Exporting data for blupf90
#Selecting variable of interest
data5<- data4 %>% select(cow,milk,ward,herd,cyrsn,tyrmn,dgrp,lac,lacgr,age,long,lat,intercept,leg1,leg2)
#1. cow
#2. milk
#3. ward
#4. herd
#5. cyrsn
#6. tyrmn
#7. dgrp
#8.lac
#9.lacgr
#10.age
#11.long
#12.lat
#13.intercept
#14.leg1
#15.leg2
#Exporting the data for blupf90
write.table(data5, "data.dat",
           quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")


