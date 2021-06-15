#Isidore Houaga, 25/05/2021
# MACE Project
#South Africa Data
#Holstein data
# Scripts Preliminary data cleaning and descriptive statistics
#Importing South Africa (SA) holstein data (two datasets HST_01 and HST_02)

# Importing HSTSA_01 dataset
HSTSA1 <- read.csv(file = "HSTSA_01.csv", stringsAsFactors = FALSE)
summary(HSTSA1)
# Importing HST_02 dataset
HSTSA2 <- read.csv(file = "HSTSA_02.csv", stringsAsFactors = FALSE)
summary(HSTSA2)
library(dplyr)
glimpse(HSTSA1)
glimpse(HSTSA2)
# check Columns names of the two datasets before merging
colnames(HSTSA1)==colnames(HSTSA2)

# Merging the two datasets by multiple columns into HSTSA
HSTSA<- merge(HSTSA1, HSTSA2, by = c('ANIMAL_NUMBER','CALVING_DTM', 'LACTATION_DTM','PARTICIPANT', 'PARITY', 'END_LACTATION_DTM', 'MILK_YLD','FAT_YLD','PROTEIN_YLD','MILK_YLD_305D','FAT_YLD_305D','PROTEIN_YLD_305D','LACT_LNGTH','LACT_INDEX', 'FARMERS_NO'))
summary(HSTSA)

HSTSA$PARTICIPANT<- as.factor(HSTSA$PARTICIPANT)
HSTSA$PARITY<- as.factor(HSTSA$PARITY)

#as.Date('1/15/2001',format='%m/%d/%Y')
HSTSA$MILK_YLD_305D<- as.numeric(HSTSA$MILK_YLD_305D)
summary(HSTSA)
# Minimum Milk yield=0 found. These are treated as missing data as discussed with Yinka
#Maximum of 23300 litres seems high, minimum and maximum sets need to be discussed
## Checking missing data
is.na(HSTSA)

#Use the "apply" function in combination with the "gsub" function to remove all whitespace from each data frame column.
HSTSA_new1 <- as.data.frame(apply(HSTSA,              # Remove blanks
                    2,
           function(x) gsub("\\s+", "", x)))
summary(HSTSA_new1)

#Set Blank cells & cells with Space & "Unknown" to NA in Data (HSTSA)
HSTSA_new2 <- HSTSA_new1                                    
HSTSA_new2[HSTSA_new2 == ""|HSTSA_new2 == " "|HSTSA_new2 == "0"|HSTSA_new2 == "0.00"|HSTSA_new2== "Unknown"] <- NA 
summary(HSTSA_new2)        
#Summary HSTSA_new2 (Raw data Holstein South Africa)
library(dplyr)
library(ggplot2)
HSTSA_new2$PARTICIPANT<- as.factor(HSTSA_new2$PARTICIPANT)
HSTSA_new2$PARITY<- as.factor(HSTSA_new2$PARITY)
HSTSA_new2$MILK_YLD <- as.numeric(HSTSA_new2$MILK_YLD)
HSTSA_new2$FAT_YLD  <- as.numeric(HSTSA_new2$FAT_YLD) 
HSTSA_new2$PROTEIN_YLD <- as.numeric(HSTSA_new2$PROTEIN_YLD)
HSTSA_new2$MILK_YLD_305D<- as.numeric(HSTSA_new2$MILK_YLD_305D)
HSTSA_new2$FAT_YLD_305D <- as.numeric(HSTSA_new2$FAT_YLD_305D)
HSTSA_new2$PROTEIN_YLD_305D <- as.numeric(HSTSA_new2$PROTEIN_YLD_305D)
HSTSA_new2$LACT_LNGTH <- as.numeric(HSTSA_new2$LACT_LNGTH)
#as.Date('1/15/2001',format='%m/%d/%Y')
#Histogram of 305 Milk yield
hist(na.omit(HSTSA_new2$MILK_YLD_305D))

################################################################################################################################################################################
###Counting obrservation in R
#######sum(with(data,gender == "M" & stream == "Commerce"))
#Alternatively
#library(dplyr)
##nrow(filter(data,gender == "M" & stream == "Commerce"))

###############################
#Number of lactation milk yield records
HSTSA_new2$MILK_YLD
length(HSTSA_new2$MILK_YLD)#282440
#Number of 305days milk yield (MILK_YLD_305D) records
HSTSA_new2$MILK_YLD_305D
length(HSTSA_new2$MILK_YLD_305D)##282440

## Number of missing lactation milk yield and 305days milk yield
sum(is.na(HSTSA_new2$MILK_YLD))###40317
sum(is.na(HSTSA_new2$MILK_YLD_305D)) ##128925
glimpse(HSTSA_new2)

#Number of cows
ncows<- length(unique(HSTSA_new2$ANIMAL_NUMBER))
ncows ##96307

# Number of herds Holstein South Africa
nherds<- length(unique(HSTSA_new2$PARTICIPANT))
nherds# 5855
### Clear the environment
### rm(list=ls())

#Total number of missing values
sum(is.na(HSTSA_new2))### NA=1100110
#Percentage of missing data
colMeans(is.na(HSTSA_new2))*100
#Total number of missing values of MILK_YLD_305D 
sum(is.na(HSTSA_new2$MILK_YLD_305D))###128925
# Number of observed parities
table(HSTSA_new2$PARITY)

#MILK_YLD_305D, FAT_YLD_305D and PROTEIN_YLD_305D across all lactation 

mean(HSTSA_new2$MILK_YLD_305D,na.rm = TRUE)### 7227.577+2241.825
sd(HSTSA_new2$MILK_YLD_305D,na.rm = TRUE)

mean(HSTSA_new2$FAT_YLD_305D,na.rm = TRUE)### 3266.78+3658.95
sd(HSTSA_new2$FAT_YLD_305D,na.rm = TRUE)

mean(HSTSA_new2$PROTEIN_YLD_305D,na.rm = TRUE)### 3974.056+3551.132
sd(HSTSA_new2$PROTEIN_YLD_305D,na.rm = TRUE)
## Lactation length
mean(HSTSA_new2$LACT_LNGTH,na.rm = TRUE) ###229.2397+100.8147
sd(HSTSA_new2$LACT_LNGTH,na.rm = TRUE)
 
 # 305D Milk yield by Parity ##### @Thiago, please check the following two codes for plots
 library(tidyr)
 newdata <- HSTSA_new2 %>%
 ggplot(aes(y = MILK_YLD_305D_1, x = PARITY)) +
   geom_point() +
   ylab("305D Milk Yield") +
     xlab("PARITY") +
      theme_bw()
 table(HSTSA_new2$PARITY)
 
 # Number of observation by Parity
 library(tidyr)
 newdata%>%
   ggplot(aes(y = MILK_YLD_305D_3, x = PARITY)) +
   geom_point() +
     ylab("Number of observation") +
   xlab("PARITY") +
   theme_bw()
 
 # Milk yield by Parity, Boxplot to show the distribution
 library(tidyr)
 HSTSA_new2 %>%
   drop_na() %>%
   ggplot(aes(y = MILK_YLD_305D, x = PARITY)) +
   geom_boxplot() +
   ylab("305D Milk Yield") +
   xlab("PARITY") +
   theme_bw()
 table(HSTSA_new2$PARITY)
 
##########################################################################################################
##Data editing ##
###########################################################################################################
# Remove rows with parity more or equal to 10
 HSTSA_new2$PARITY<- as.numeric(HSTSA_new2$PARITY)
 HSTSA_new3<- subset(HSTSA_new2, HSTSA_new2$PARITY < 9 | HSTSA_new2$PARITY == "9")
 table(HSTSA_new3$PARITY)
 # Duplicate participant column and rename it "TYPE"
 HSTSA_new3$TYPE<- HSTSA_new3$PARTICIPANT
 HSTSA_new3
 ## Relocate TYPE near participant (herd)
 HSTSA_new4<- subset(HSTSA_new3, select = c(1,2,3,4,16,5,6,7,8,9,10,11,12,13,14,15))
 head(HSTSA_new4)
 
 ##### Extract the strings from numeric in the HSTSA_new4$TYPE 
 
 