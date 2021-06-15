#Isidore Houaga, 25/05/2021
# MACE Project
#South Africa Data
#Holstein data
# Scripts Preliminary data cleaning and descriptive statistics
#Importing South Africa (SA) holstein data (two datasets HST_01 and HST_02)
#Loading Libraries

library(data.table)
library(bit64)
# Importing HST_01 dataset
HSTSA1 <- fread("C:\\Users\\Lenovo\\OneDrive\\Documents\\MACE\\MACE\\HSTSA_01.csv")
# Importing HST_02 dataset
HSTSA2 <- fread("C:\\Users\\Lenovo\\OneDrive\\Documents\\MACE\\MACE\\HSTSA_02.csv")
library(dplyr)
glimpse(HSTSA1)
glimpse(HSTSA2)
# Columns names
colnames(HSTSA1)==colnames(HSTSA2)

# Merging the two datasets by multiple columns into HSTSA
HSTSA<- merge(HSTSA1, HSTSA2, by = c('ANIMAL_NUMBER','CALVING_DTM', 'LACTATION_DTM','PARTICIPANT', 'PARITY', 'END_LACTATION_DTM', 'MILK_YLD','FAT_YLD','PROTEIN_YLD','MILK_YLD_305D','FAT_YLD_305D','PROTEIN_YLD_305D','LACT_LNGTH','LACT_INDEX', 'FARMERS_NO'))
summary(HSTSA)
# Set Blank & Space & "Unknown" to NA in Data (HSTSA)

#Use the "apply" function in combination with the "gsub" function to remove all whitespace from each data frame column.
HSTSA_new1 <- as.data.frame(apply(HSTSA,              # Remove blanks
                    2,
           function(x) gsub("\\s+", "", x)))
summary(HSTSA_new1)

#Set Blank cells in dataset to missing data "NA"
HSTSA_new2 <- HSTSA_new1                                    
HSTSA_new2[HSTSA_new2 == ""|HSTSA_new2 == " "|HSTSA_new2== "Unknown"] <- NA 
summary(HSTSA_new2)        

#Summary HSTSA_new2
library(dplyr)
library(ggplot2)
HSTSA_new2$PARTICIPANT<- as.factor(HSTSA_new2$PARTICIPANT)
HSTSA_new2$PARITY<- as.factor(HSTSA_new2$PARITY)
HSTSA_new2$FARMERS_NO<- as.factor(HSTSA_new2$FARMERS_NO)
################################################################################################################################################################################
###Counting obrservation in R
#######sum(with(data,gender == "M" & stream == "Commerce"))
#Alternatively
#library(dplyr)
##nrow(filter(data,gender == "M" & stream == "Commerce"))

###############################
#Number of milk yield records
HSTSA_new2$MILK_YLD
length(HSTSA_new2$MILK_YLD)#282440

#Number of cows ### 53919
HSTSA_new2$ANIMAL_NUMBER<- as.factor(HSTSA_new2$ANIMAL_NUMBER)
ncowsHSTSA<- levels(HSTSA_new2$ANIMAL_NUMBER)
summary(ncowsHSTSA) ### Number of cows 53919
length(ncowsHSTSA)  ### Number of cows 53919

#Number of herds Holstein South Africa

nherdsHSTSA<- levels(HSTSA_new2$PARTICIPANT)
length(nherdsHSTSA)  ### Number of herds=5855

# Number of farmers ###36005
nfarmersHSTSA<- levels(HSTSA_new2$FARMERS_NO)
length(nfarmersHSTSA)

#Total number of missing values
sum(is.na(HSTSA_new2))### NA=994399
#otal number of missing values of ILK_YLD_305D ##
sum(is.na(HSTSA_new1$MILK_YLD_305D))

# Checking total missing values in a dataset
#sum(is.na(df)) #  For entire dataset
#for a particular column in a dataset

#sum(is.na(df$col1)) 
#Or to check for all the columns as mentioned by @nicola

#colSums(is.na(df))

#Check percentages and counts of missing values in columns:
   
   #colMeans(is.na(data))*100
#colSums(is.na(data))

head(HSTSA_new1)
#MILK_YLD_305D FAT_YLD_305D PROTEIN_YLD_305D
mean(HSTSA_new1$MILK_YLD_305D,na.rm = TRUE)### 6093.334+2181.314
sd(HSTSA_new1$MILK_YLD_305D,na.rm = TRUE)

mean(HSTSA_new1$FAT_YLD_305D,na.rm = TRUE)### 217.5417+79.74909
sd(HSTSA_new1$FAT_YLD_305D,na.rm = TRUE)

mean(HSTSA_new1$PROTEIN_YLD_305D,na.rm = TRUE)### 178.0691+91.35244
sd(HSTSA_new1$PROTEIN_YLD_305D,na.rm = TRUE)

#Compute the average 305D Milk yield by PARTICIPANT
MY305DPA<- HSTSA_new1 %>%
  group_by(PARTICIPANT) %>%
  summarise(mean305DMY_Participant = mean(MILK_YLD_305D))
head(MY305DPA)

#Compute the average 305D Milk yield by PARITY
HSTSA_new1 %>%
  group_by(PARITY) %>%
  summarise(mean305DMY_Parity = mean(MILK_YLD_305D,na.rm=TRUE)) %>%
ggplot(aes(x = PARITY, y = mean305DMY_Parity, fill = PARITY)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(
    x = "Parity",
    y = "Average 305D Milk yield (Kg)",
    title = paste(
      "Average 305 days milk yield by parity"
    )
  )

# LinePlot#### No line on output
HSTSA_new1 %>%
  group_by(PARITY) %>%
  summarise(mean305DMY_Parity = mean(MILK_YLD_305D,na.rm=TRUE))
 ggplot(aes(x = PARITY, y = mean305DMY_Parity)) +
  geom_line() +
  theme_classic() +
  labs(
    x = "Parity",
    y = "Average 305D Milk yield (Kg)",
    title = paste(
      "Average 305 days milk yield by parity"
    )
  )
 #Compute the average 305D fat yield by PARITY
 HSTSA_new1 %>%
   group_by(PARITY) %>%
   summarise(mean305DFY_Parity = mean(FAT_YLD_305D,na.rm=TRUE)) %>%
   ggplot(aes(x = PARITY, y = mean305DFY_Parity, fill = PARITY)) +
   geom_bar(stat = "identity") +
   theme_classic() +
   labs(
     x = "Parity",
     y = "Average 305D Fat yield (Kg)",
     title = paste(
       "Average 305 days fat yield by parity"
     )
   )
 #Compute the average 305D Milk yield by Participant
 HSTSA_new1 %>%
   group_by(PARTICIPANT) %>%
   summarise(mean305DMY_Participant = mean(MILK_YLD_305D,na.rm=TRUE)) %>%
   ggplot(aes(x = PARTICIPANT, y = mean305DMY_Participant, fill = PARTICIPANT)) +
   geom_bar(stat = "identity") +
   theme_classic() +
   labs(
     x = "Participant",
     y = "Average 305D Milk yield (Kg)",
     title = paste(
       "Average 305 days milk yield by participant"
     )
   )
 str(HSTSA_new1$PARTICIPANT)
 summary(HSTSA_new1$FARMERS_NO) 
 
 #Scripts from Thiago
 HSTSA_new1 %>%
  group_by(PARTICIPANT, PARITY, FARMERS_NO) %>%
 summarise(across(everything(), list(mean, sd))) %>%
 mutate_at(vars(PARTICIPANT, PARITY, FARMERS_NO), factor) %>%
   ggplot(aes(y = MILK_YLD_305D, x = PARITY,  linetype = Type)) +
    geom_line() +
  ylab("305D Mean Yield") +
  xlab("PARITY") +
  labs(fill = "Mean + 1.96 SD") +
   geom_vline(xintercept = 20, alpha = 0.3,  linetype = 2) +
      theme_bw()
   ##########
 newdata <- HSTSA_new1 %>%
      group_by(PARITY) %>%
   drop_na()%>%
   summarise(across(MILK_YLD_305D:PROTEIN_YLD_305D, list(mean, sd, length)))
 summary(newdata)
 #head(newdata)
 # Milk yield by Parity
 library(tidyr)
 newdata%>%
 ggplot(aes(y = MILK_YLD_305D_1, x = PARITY)) +
   geom_point() +
   ylab("305D Milk Yield") +
     xlab("PARITY") +
      theme_bw()
 table(HSTSA_new1$PARITY)
 
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
 HSTSA_new1 %>%
   drop_na() %>%
   ggplot(aes(y = MILK_YLD_305D, x = PARITY)) +
   geom_boxplot() +
   ylab("305D Milk Yield") +
   xlab("PARITY") +
   theme_bw()
 table(HSTSA_new1$PARITY)
