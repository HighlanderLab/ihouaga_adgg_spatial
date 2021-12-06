#=======================================================================
# Holstein (South Africa) data analysis                                #
# Isidore,  Thiago,  Ivan and Gregor                                   #
# Version 1.0.1                                                        #
# Date: 12/12/2021                                                     #
# Descriptive analysis and data manipulation                           #
# Importing South Africa (SA) holstein data (two datasets HST_01 and   #
# HST_02)                                                              #
#=======================================================================
#-----------------------------------------------------------------------
# Uploading packages   
#-----------------------------------------------------------------------
library(dplyr)
library(stringr)
library(ggplot2)     
library(ggExtra)
library(EnvStats)
library(lme4)     #for mixed model
library(afex)     #for p-values
library(readr)
library(GGally)
library(sf) # import shapefile
library(sp) # poly2nb for graph
library(INLA) # inla.read.graph
library(spData)
library(spdep)
(.packages()) # Checking loaded packages
#=======================================================================
# Uploading datasets
#=======================================================================
# Importing HST_01 dataset

HSTSA1 <- read.csv(file = "HSTSA_01.csv", stringsAsFactors = FALSE)

# Importing HST_02 dataset
HSTSA2 <- read.csv(file = "HSTSA_02.csv", stringsAsFactors = FALSE)
#=======================================================================
# Data organization
#=======================================================================
dataHST1 <- HSTSA1 %>%
  mutate_at(c("ANIMAL_NUMBER", "PARTICIPANT","PARITY", "FARMERS_NO"),  factor)

  dataHST2 <- HSTSA2 %>%
  mutate_at(c("ANIMAL_NUMBER", "PARTICIPANT","PARITY", "FARMERS_NO"),  factor)
#=======================================================================
# Stacking the data
#=======================================================================
  glimpse(dataHST1)
  glimpse(dataHST2)
  dataHST <- rbind(dataHST1, dataHST2)
dim(dataHST)[1] == dim(dataHST1)[1] + dim(dataHST2)[1] # Expect TRUE
#=======================================================================
# Organizing date variables
#=======================================================================
dataHST$CALVING_DTM <- as.Date(dataHST$CALVING_DTM,
                               format = "%d/%m/%Y")
dataHST$LACTATION_DTM <- as.Date(dataHST$LACTATION_DTM,
                                 format = "%d/%m/%Y")
dataHST$END_LACTATION_DTM <- as.Date(dataHST$END_LACTATION_DTM,
                                     format = "%d/%m/%Y")
summary(dataHST)
### 144,774 records from Unknown herd (PARTICIPANT)

#=======================================================================
# Summary of Raw data (number of records, cows and herds)
#=======================================================================
## Periods of calving: between 1949-05-18 and 2014-09-18, lactation date between 1949-05-18 and 2014-09-26
## Number of records Raw data
length(dataHST$ANIMAL_NUMBER) ## 1331015
## Number of cows raw data
length(unique(dataHST$ANIMAL_NUMBER))## 377921
## Number of herds Raw data
length(unique(dataHST$PARTICIPANT)) ## 7121
#### Summarry of Unknown herds
DataUnknownherd<- dataHST %>%
  filter(dataHST$PARTICIPANT=="Unknown")
## Number of records from Unknown herd
length(DataUnknownherd$ANIMAL_NUMBER) ## 144,774 records.
## Number of cows from Unknown herd
length(unique(DataUnknownherd$ANIMAL_NUMBER)) ## 64465
#=======================================================================
# Factor PARTICIPANT -> Unknown level to NA
#=======================================================================
dataHST$PARTICIPANT[dataHST$PARTICIPANT == "Unknown"] <- NA
#=======================================================================
# Two columns to identify commercial and Stud animals
#=======================================================================
dataHST$Type <- str_extract_all(as.character(dataHST$PARTICIPANT),
                                "[:upper:]", simplify = TRUE)
dataHST$Class <- paste0(dataHST$Type[, 1], dataHST$Type[, 2],
                        dataHST$Type[, 3])
summary(dataHST$Class)
summary(dataHST$Type)
dataHST$Class <- as.factor(dataHST$Class)
dataHST <- subset(dataHST,  select = -c(Type))
#-----------------------------------------------------------------------
dataHST$herd <- readr::parse_number(as.character(dataHST$PARTICIPANT))
summary(dataHST)
summary(dataHST$Class)
#=======================================================================
# One column to identify breed (SAHOL: South African Holstein)
#=======================================================================
dataHST<- data.frame(append(dataHST, c(BREED='SAHOL'), after=1))
glimpse(dataHST)
summary(dataHST$Class)
write.csv(dataHST, "/Users/ihouaga/Documents/Project1/dataHSTRaw.csv") ## save .csv copy of raw data
#=======================================================================
# Check for duplicates
#=======================================================================
#Remove potential duplicated rows
dataHST2 <- dplyr::distinct(dataHST) ###Efficient duplicate removal
length(dataHST2$ANIMAL_NUMBER) # New number of records 1048575
length(unique(dataHST2$ANIMAL_NUMBER))### New number of cows is 377921
length(unique(dataHST2$PARTICIPANT)) ### 7121; no herd removed. 
## Number of removed duplicated records (rows) = 1331015-1048575=282440
###There was 282440 duplicated rows

#=======================================================================
# Organizing the whole data set - removing NA and Zeros
dataHST3<- dataHST2 %>%
  filter(Class %in% c("HST", "MCM", "NA", "UNK")) %>% #Remove animals with class JSE
  droplevels()
  length(dataHST3$ANIMAL_NUMBER) #1048568
  ####Number of records removed ie records of JSE class =1048568-1048575 =7
  ### Remove NAs in 305D Milk yield
  dataHST4<- dataHST3 %>%
  filter(is.na(MILK_YLD_305D) == FALSE)
  droplevels() ###Error message "no applicable method for 'droplevels' applied to an object of class "NULL"
  #Number of records after removing MILK_YLD_305D missing data 
  length(dataHST4$ANIMAL_NUMBER) # 606954
  # 1048568-606954= 441614 records removed.
  #Number of cows dataHST4
  length(unique(dataHST4$ANIMAL_NUMBER)) #254314
  ##Number of cows removed = 377921-254314=123607
  #123607 cows removed
  #Number of herds dataHST4
  length(unique(dataHST4$PARTICIPANT)) # 5364 herds
  # 7121-5364=1757 herds removed
  ##  Remove rows with MILK_YLD_305D < 305 kg (Makaglela et al., 2007; Banga, 2009)
  dataHST5<- dataHST4 %>%
  filter(MILK_YLD > 305) %>% # Remove cows with milk yield < 305 litres (Makaglela et al., 2007; Banga, 2009)
    droplevels()
     #Number of records after filtering MILK_YLD > 305
  length(dataHST5$ANIMAL_NUMBER) ##606382
  ### Number of records removed 606954-606382=572 records
 #### Number of removed cows
  length(unique(dataHST4$ANIMAL_NUMBER))-length(unique(dataHST5$ANIMAL_NUMBER)) ### 214 cows removed
  length(unique(dataHST4$ANIMAL_NUMBER)) ##254314
  length(unique(dataHST5$ANIMAL_NUMBER))## 254100
  levels(dataHST5$Class)
  dataHSTUNK<- dataHST5 %>%
    filter(Class=="UNK")
  droplevels() # 59 records of UNK. Error message droplevels' not applied to an object of class "NULL
  summary(dataHST5$Class)
levels(dataHST5$Class)[3] <- "Unclassified"
summary(dataHST5$Class)
levels(dataHST5$Class)
levels(dataHST5$Class)[4] <- "Unclassified" ##UNK herds treated as unclassifided
summary(dataHST5$Class)
summary(dataHST5)
#=======================================================================
# Let's create Class2 with levels Classified (HST and HCM) and Unclassified
#=======================================================================
dataHST5$Class2<- dataHST5$Class
dataHST5
levels(dataHST5$Class)
summary(dataHST5$Class2)
levels(dataHST5$Class2)[1:2]<- "Classified"
summary(dataHST5$Class2)
summary(dataHST5$Class)
#=======================================================================
# Let's create new dataset by removing Unclassified animal
#=======================================================================
dataHSTnew<- dataHST5 %>%
  filter(Class %in% c("HST", "MCM")) %>%
  droplevels()
summary(dataHSTnew)
summary(dataHSTnew$Class)
#Number of records removed
length(dataHSTnew$ANIMAL_NUMBER) ### 563765
### 606382-563765=42617 records removed
#Number of cows removed
length(unique(dataHSTnew$ANIMAL_NUMBER)) #232516
##254100-232516=21584 cows removed
## Number of herds removed
length(unique(dataHST5$PARTICIPANT))-length(unique(dataHSTnew$PARTICIPANT))
### 37 herds removed
### Filtering PARITY
dataHSTnew$PARITY<- as.numeric(dataHSTnew$PARITY)
summary(dataHSTnew$PARITY)
dataHSTnew1<- dataHSTnew %>%
  filter(PARITY<10) %>%  
  droplevels()
summary(dataHSTnew1$PARITY)
dataHSTnew1$PARITY <- as.factor(dataHSTnew1$PARITY)
### Number of records, animals and herds removed 
#Records removed
length(dataHSTnew$ANIMAL_NUMBER)-length(dataHSTnew1$ANIMAL_NUMBER)# 1616
# Cows removed
length(unique(dataHSTnew$ANIMAL_NUMBER))-length(unique(dataHSTnew1$ANIMAL_NUMBER))##104
# Herds removed 
length(unique(dataHSTnew$PARTICIPANT))-length(unique(dataHSTnew1$PARTICIPANT)) # 1 herd removed
#=======================================================================
# Summary
#=======================================================================
summary(dataHSTnew1)
#=======================================================================
# Descriptive Statistics
#=======================================================================
#Number of lactation milk yield records
length(dataHSTnew1$ANIMAL_NUMBER) #562149
#Number of cows
length(unique(dataHSTnew1$ANIMAL_NUMBER))
ncows<- length(unique(dataHSTnew1$ANIMAL_NUMBER))
ncows ##232412
#Number of cows in Parity 1
PARITY1<- subset(dataHSTnew1, dataHST$PARITY=="1")
ncows1<- length(unique(PARITY1$ANIMAL_NUMBER))
ncows1 ##113007

#Number of cows in Parity 2
PARITY2<- subset(dataHSTnew1, dataHST$PARITY=="2")
ncows2<- length(unique(PARITY2$ANIMAL_NUMBER))
ncows2 ##113174

#Number of cows in Parity 3
PARITY3<- subset(dataHSTnew1, dataHST$PARITY=="3")
ncows3<- length(unique(PARITY3$ANIMAL_NUMBER))
ncows3 ## 97875

#Number of cows in Parity 4
PARITY4<- subset(dataHSTnew1, dataHST$PARITY=="4")
ncows4<- length(unique(PARITY4$ANIMAL_NUMBER))
ncows4 ##65844

#Number of cows in Parity 5
PARITY5<- subset(dataHSTnew1, dataHST$PARITY=="5")
ncows5<- length(unique(PARITY5$ANIMAL_NUMBER))
ncows5 ##43851

#Number of cows in Parity 6
PARITY6<- subset(dataHSTnew1, dataHST$PARITY=="6")
ncows6<- length(unique(PARITY6$ANIMAL_NUMBER))
ncows6 ##26855

#Number of cows in Parity 7
PARITY7<- subset(dataHSTnew1, dataHST$PARITY=="7")
ncows7<- length(unique(PARITY7$ANIMAL_NUMBER))
ncows7 ##14976

#Number of cows in Parity 8
PARITY8<- subset(dataHSTnew1, dataHST$PARITY=="8")
ncows8<- length(unique(PARITY8$ANIMAL_NUMBER))
ncows8 ## 7421

#Number of cows in Parity 9
PARITY9<- subset(dataHSTnew1, dataHST$PARITY=="9")
ncows9<- length(unique(PARITY9$ANIMAL_NUMBER))
ncows9 ## 3309

# Number of herds Holstein South Africa dataHSTnew1
nherds<- length(unique(dataHSTnew1$PARTICIPANT))
nherds# 5323
#=======================================================================
# Clear the environment
# rm(list=ls())
#=======================================================================
#MILK_YLD_305D across 9 lactations 

mean(dataHSTnew1$MILK_YLD_305D,na.rm = TRUE)### 6658.969±2385.502
sd(dataHSTnew1$MILK_YLD_305D,na.rm = TRUE)
hist(na.omit(dataHSTnew1$MILK_YLD_305D)) ###Error that figure margins too large

mean(dataHSTnew1$FAT_YLD_305D,na.rm = TRUE)### 240.5781±89.77458
sd(dataHSTnew1$FAT_YLD_305D,na.rm = TRUE)

mean(dataHSTnew1$PROTEIN_YLD_305D,na.rm = TRUE)### 203.7155±88.65431
sd(dataHSTnew1$PROTEIN_YLD_305D,na.rm = TRUE)

#=======================================================================
# MILK 305
#=======================================================================
# Mean and Variance - grouped by PARITY
dataHSTnew1 %>%
  group_by(PARITY) %>%
  summarise(
    MeanMilk = mean(MILK_YLD_305D, na.rm = TRUE),
    MedianMilk = median(MILK_YLD_305D, na.rm = TRUE),
    VarMilk = var(MILK_YLD_305D, na.rm = TRUE),
    SdMilk =  sqrt(VarMilk),
    nMilk = n()
  )
#-----------------------------------------------------------------------
# Mean and Variance - grouped by Class
dataHSTnew1$herd<- as.factor(dataHSTnew1$herd)
dataHSTnew1 %>%
  group_by(Class) %>%
  summarise(
    MeanMilk = mean(MILK_YLD_305D, na.rm = TRUE),
    MedianMilk = median(MILK_YLD_305D, na.rm = TRUE),
    VarMilk = var(MILK_YLD_305D, na.rm = TRUE),
    SdMilk =  sqrt(VarMilk),
    nMilk = n()
  )
summary(dataHSTnew1$Class)
dataHSTnew1 %>%
  filter(Class != "NA") %>%
  group_by(herd, Class) %>%
  summarise(
    MeanMilk = mean(MILK_YLD_305D, na.rm = TRUE),
    MedianMilk = median(MILK_YLD_305D, na.rm = TRUE),
    VarMilk = var(MILK_YLD_305D, na.rm = TRUE),
    SdMilk =  sqrt(VarMilk),
    nMilk = n()
  ) %>%
  ggplot(aes(y = MeanMilk, x = SdMilk)) +
  geom_point(shape = 1) +
  facet_grid(~Class) +
  xlab("SD(Milk)") +
  ylab("Milk Average") +
  theme_bw()
ggsave("herd_mean_sd.png")

# We have to investigate why there are herds with zero SD
# You can talk on this in the lab meeting -> Same as in MILK_YIELD
#-----------------------------------------------------------------------
# Filtering herds with variance zero for 305D milk yield
#-----------------------------------------------------------------------
sd305 <- dataHSTnew1 %>%
  filter(Class != "NA") %>%
  group_by(herd, Class) %>%
  summarise(
    Mean305DMilk = mean(MILK_YLD_305D, na.rm = TRUE),
    Median305DMilk = median(MILK_YLD_305D, na.rm = TRUE),
    Var305DMilk = var(MILK_YLD_305D, na.rm = TRUE),
    Sd305DMilk =  sqrt(Var305DMilk),
    n305DMilk = n()
  ) %>%
  filter(Var305DMilk <= 50) %>%
  droplevels()

sd305 %>%
  ggplot(aes(y = Mean305DMilk,  x = n305DMilk,  group = Var305DMilk)) +
  geom_point(aes(colour = Var305DMilk)) +
  xlab("Number of Observations") +
  ylab("305DMilk Average") +
  theme_bw()
# All herds with variance zero has 2 identical observations
# We have to discover the reason why they have two identical
# observations!!
# Lets check for duplicate
duplicated(dataHSTnew1) 
# Lets check for duplicated 305D MILK_YLD
duplicated(dataHSTnew1$MILK_YLD_305D)#### There was duplicated 305d milk yield
# Lets extract duplicated MILK_YLD
dataHSTnew1$MILK_YLD_305D[duplicated(dataHST$MILK_YLD_305D)]
dataHSTnew1
### Check for what is happening with duplicated 305 milk yield 
# eg: MILK_YLD== 11621 Litre or 
MY9152<- subset(dataHSTnew1, dataHSTnew1$MILK_YLD_305D == "9152")
MY9152
MY9152_distinct<- MY9152 %>% distinct()
MY9152_unique<- unique(MY9152)
MY9152==MY9152_distinct### True. There is no duplicated rows. The duplicated MILK_YLD_305D are likely to be 
# different animal from same herd having same performance or herd with single cow
#-----------------------------------------------------------------------
# Mean and Variance - groupd by Class and PARITY
SMilk <- dataHSTnew1 %>%
  group_by(PARITY, Class) %>%
  summarise(
    MeanMilk = mean(MILK_YLD_305D, na.rm = TRUE),
    MedianMilk = median(MILK_YLD_305D, na.rm = TRUE),
    VarMilk = var(MILK_YLD_305D, na.rm = TRUE),
    SdMilk =  sqrt(VarMilk),
    nMilk = n()
  )

# Mean vs. Variance
SMilk %>%
  ggplot(aes(y = MeanMilk, x = SdMilk,  colour = PARITY)) +
  geom_text(aes(label = PARITY)) +
  xlab("SD(Milk)") +
  ylab("E(Milk)") +
  theme_bw()
# Linear relationship between mean and variance
# Mean vs. Variance grouped by Class
SMilk %>%
  ggplot(aes(y = MeanMilk, x = SdMilk,  label = PARITY)) +
  geom_label(aes(fill = factor(PARITY)), colour = "white",
             fontface = "bold") +
  facet_wrap(~Class) +
  xlab("SD(Milk 305)") +
  ylab("E(Milk 305)") +
  theme_bw()
ggsave("parity.png")
# # Here we can see:
## HST: High mean and high variance compared to MCM and NA HST: Mean and
## SD increases from Parity 1 to parity 3 and starts to decrease from
## parity 4

## MCM: Mean lower than HST but higher than NA
## Distance between Parity 1 to 2 is higher than observed to HSTSA
## Mean decrease is slower in MCM compared to HST -> you can see that
## numbers 9 and 8 are close on MCM compared to HSTSA

# NA: The lowest mean and SD. Why this is happeneing here? Which kind
# of herd/animals we have as NA?? Seems a different group. If NA were
# composed by animals from HST and MCM I expected a Mean-Variance
# relationship between HST and MCM -> something to think!
#-----------------------------------------------------------------------
# Mean vs. number of observations
SMilk %>%
  ggplot(aes(y = MeanMilk, x = nMilk,  label = PARITY)) +
  geom_label(aes(fill = factor(PARITY)), colour = "white",
             fontface = "bold") +
  facet_wrap(~Class) +
  xlab("Number of Observations") +
  ylab("E(Milk 305)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# We can see the Milk 305 production increases up to PARITY 3-4 and then
# starts to decrease. This is related to the combination of effect of
# PARITY and number of observations

# SD vs. number of observations
SMilk %>%
  ggplot(aes(y = SdMilk, x = nMilk,  label = PARITY)) +
  geom_label(aes(fill = factor(PARITY)), colour = "white",
             fontface = "bold") +
  facet_wrap(~Class) +
  xlab("Number of Observations") +
  ylab("SD(Milk 305)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Same for Standard deviation
#-----------------------------------------------------------------------
# Mean and Median relationship
SMilk %>%
  ggplot(aes(y = MeanMilk,  x = MedianMilk)) +
  geom_abline(intercept =0, slope=1) + 
  facet_wrap(~Class) +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("Median(Milk 305)") +
  ylab("E(Milk 305)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Ok -> linear relationship
#=======================================================================
# Centering the variable Milk 305
#=======================================================================
dataHSTnew1 %>%
  mutate_at("PARITY", factor) %>%
  mutate(MYc = MILK_YLD_305D-mean(MILK_YLD_305D)) %>%
  ggplot(aes(y = MYc, x = PARITY)) +
  geom_boxplot() +
  facet_grid(~ Class) +
  xlab("Levels of Parity") +
  ylab("Centred Milk 305 Yield") +
  theme_bw()

# Distribution is assymetric for both HST and MCM class

#=======================================================================
# Pairs
#=======================================================================
ggpairs(dataHSTnew1[, c(5, 7:12)], aes(colour = factor(PARITY),
                                   alpha = 0.4))
ggsave("pairsPlot.png",  width = 14, height = 14)
ggsave("pairsPlot.png",  width = 16, height = 16) ## @Thiago:Check output please

## @Thiago, please see error message below
"Column 'PARTICIPANT' has more levels (5323) than the threshold (15) allowed
Please remove the column or increase the 'cardinality_threshold' parameter 
Increasing the cardinality_threshold may produce long processing times"
#=======================================================================
# Extract birth date,sire and dam's IDs from Holstein pedigree data 
#=======================================================================
#Import the Holstein pedigree data of South Africa
HST_PED <- read.csv(file = "HST_PED.csv", stringsAsFactors = FALSE)
#Extract ANI_ID, SIRE_ID, DAM_ID and BIRTH_DATE
HST_PED_HST_SA<- HST_PED %>%
  select(ANIM_ID,SIRE_ID, DAM_ID,BIRTH_DATE)
HST_PED_HST_SA
# Rename  HST_PED_HST_SA$ANIM_ID to HST_PED_HST_SA$ANIMAL_NUMBER
names(HST_PED_HST_SA)[names(HST_PED_HST_SA) == 'ANIM_ID'] <- 'ANIMAL_NUMBER'
HST_PED_HST_SA 
# Bring the Cow birth date, sire and dam-IDs to the peformance data
dataHSTnew2 <- merge(dataHSTnew1, HST_PED_HST_SA, by= "ANIMAL_NUMBER", all.x = TRUE)
dataHSTnew2
head(dataHSTnew2)
### Sire ID and Dam ID NAs to 0
dataHSTnew2 <- mutate_at(dataHSTnew2, c("SIRE_ID", "DAM_ID"), ~replace(., is.na(.), 0))
dataHSTnew2
# Create Age columns
dataHSTnew2$BIRTH_DATE <- as.Date(dataHSTnew2$BIRTH_DATE,
                                  format = "%d/%m/%Y")
class(dataHSTnew2$BIRTH_DATE)

dataHSTnew2$Age <- as.numeric(difftime(dataHSTnew2$CALVING_DTM, dataHSTnew2$BIRTH_DATE,
                                       units ="days")) 
class(dataHSTnew2$Age)
# Rename  dataHSTNew2$ANIM_NUMBER to dataHSTNew2$COW_ID and dataHSTNew2$PARTICIPANT to dataHSTNew2$HERD_ID
names(dataHSTnew2)[names(dataHSTnew2) == 'ANIMAL_NUMBER'] <- 'COW_ID'
names(dataHSTnew2)[names(dataHSTnew2) == 'PARTICIPANT'] <- 'HERD_ID'
names(dataHSTnew2)[names(dataHSTnew2) == 'Age'] <- 'AGE'
glimpse(dataHSTnew2)
#=======================================================================
# Create Columns for additional traits 
#=======================================================================
# Remove cows without BIRTH_DATE (dataHSTnew4) or cows without CALVING_DTM 
# Or cows without LACTATION_DTM (test-day date) 
is.na(dataHSTnew2$BIRTH_DATE)
sum(is.na(dataHSTnew2$BIRTH_DATE)) ####### 132813 records from cows without birth date
summary(dataHSTnew3$BIRTH_DATE)
is.na(dataHSTnew2$CALVING_DTM) # Seems no missing calving date
sum(is.na(dataHSTnew2$CALVING_DTM))## 0 NA
is.na(dataHSTnew2$LACTATION_DTM) #Seems no missing lactation date
sum(is.na(dataHSTnew2$LACTATION_DTM)) ## 0NA
dataHSTnew3 <- dataHSTnew2[complete.cases(dataHSTnew2$BIRTH_DATE),]
sum(is.na(dataHSTnew3$BIRTH_DATE))
sum(is.na(dataHSTnew3$CALVING_DTM))
sum(is.na(dataHSTnew3$LACTATION_DTM))
###Comments @Thiago: "dataHSTnew3 <- dataHSTnew2[complete.cases(dataHSTnew2$BIRTH_DATE),]" different from
###                  "dataHSTnew3 <- na.omit(dataHSTnew2, cols=c("BIRTH_DATE"))"
##### NUmber of records from Cows without birthdate removed
length(dataHSTnew2$COW_ID)-length(dataHSTnew3$COW_ID) ##### 132813
## Number of cows removed
length(unique(dataHSTnew2$COW_ID))-length(unique(dataHSTnew3$COW_ID))
### 60038 Cows without Birthdate removed.
## Number of herds removed due to  cows without Birthdate 
length(unique(dataHSTnew2$HERD_ID))-length(unique(dataHSTnew3$HERD_ID))
## 723 herds removed due to Cows without Birthdae
#-----------------------------------------------------------------------
#Create  year and season columns based on CALVING_DTM
#-----------------------------------------------------------------------
#YEAR
dataHSTnew3$CALVING_DTM <- as.POSIXct(dataHSTnew3$CALVING_DTM, format = "%m/%d/%Y %H:%M:%S") 
dataHSTnew3$YEAR<-   format(dataHSTnew3$CALVING_DTM, format="%Y") # Extract YEAR of calving
dataHSTnew3$YEAR<- as.factor(dataHSTnew3$YEAR)
# Seasons in South Africa: April to September= Dry Season (D)
dataHSTnew3$SEASON<-   format(dataHSTnew3$CALVING_DTM, format="%m") # Extract Month of calving
colnames(dataHSTnew3)
head(dataHSTnew3)
dataHSTnew3$SEASON<- as.factor(dataHSTnew3$SEASON) # SEASON variable being created  as factor
levels(dataHSTnew3$SEASON)
levels(dataHSTnew3$SEASON)[4:9] <- "Dry" # April (4) to September (9)
levels(dataHSTnew3$SEASON)
levels(dataHSTnew3$SEASON)[1:3] <- "Wet"
levels(dataHSTnew3$SEASON)
levels(dataHSTnew3$SEASON)[3:5] <- "Wet"
levels(dataHSTnew3$SEASON)
dataHSTnew3$SEASON
#Create  HYS (Herd-Year-Season interaction) column based on HERD_ID, YEAR (Calving year) and SEASON (calving season)
dataHSTnew3$HYS <- paste(dataHSTnew3$HERD_ID, dataHSTnew3$YEAR,dataHSTnew3$SEASON, sep= "-")
dataHSTnew3$HYS
dataHSTnew3$HYS<- as.factor(dataHSTnew3$HYS)
levels(dataHSTnew3$HYS)
head(dataHSTnew3)
#### Filtering HYS: Remove levels HYS< 3 ????????????????
dataHSTnew3<- transform(dataHSTnew3, count = table(HYS)[HYS])
summary(dataHSTnew3$count.Freq)
## How many records HYS
K<- table(dataHSTnew3$count.HYS)
table(K)
### @ Ivan @GG @Thiago. We have 15708 records with HYS frequency of 1 and 9658 records with HYS frequency of 2
#Should we remove them and keep HYS group of more than 2 as in most publications???
#-----------------------------------------------------------------------
#Create column AFC (Age at first calving) = Age of each cow when PARITY=1
#-----------------------------------------------------------------------
dataHSTnew3P1<- subset(dataHSTnew3, dataHSTnew3$PARITY=="1")
dataHSTnew3P1
dataHSTnew3P1AFC<- dataHSTnew3P1 %>%
  rename(AFC=AGE) %>%
  select(COW_ID, AFC)

dataHSTnew4<- merge(dataHSTnew3, dataHSTnew3P1AFC, by= "COW_ID", all.x = TRUE)
dataHSTnew4
colnames(dataHSTnew4)
#-----------------------------------------------------------------------
#CREATE COLUMNS FOR CALVING INTERVAL
#Create columns for calving dates for each parity
#-----------------------------------------------------------------------
# Subsetting the dataset to get calving date for each parity
dataHSTnew4_P1<- subset(dataHSTnew4, dataHSTnew4$PARITY=="1")
dataHSTnew4_P2<- subset(dataHSTnew4, dataHSTnew4$PARITY=="2")
dataHSTnew4_P3<- subset(dataHSTnew4, dataHSTnew4$PARITY=="3")
dataHSTnew4_P4<- subset(dataHSTnew4, dataHSTnew4$PARITY=="4")
dataHSTnew4_P5<- subset(dataHSTnew4, dataHSTnew4$PARITY=="5")
dataHSTnew4_P6<- subset(dataHSTnew4, dataHSTnew4$PARITY=="6")
dataHSTnew4_P7<- subset(dataHSTnew4, dataHSTnew4$PARITY=="7")
dataHSTnew4_P8<- subset(dataHSTnew4, dataHSTnew4$PARITY=="8")
dataHSTnew4_P9<- subset(dataHSTnew4, dataHSTnew4$PARITY=="9")

# Calving date Parity one
CD1<- dataHSTnew4_P1 %>%
  rename(CALVING_DTM_1=CALVING_DTM) %>%
  select(COW_ID, CALVING_DTM_1)
# Calving date Parity 2
CD2<- dataHSTnew4_P2 %>%
  rename(CALVING_DTM_2=CALVING_DTM) %>%
  select(COW_ID, CALVING_DTM_2)
# Calving date Parity 3
CD3<- dataHSTnew4_P3 %>%
  rename(CALVING_DTM_3=CALVING_DTM) %>%
  select(COW_ID, CALVING_DTM_3)
# Calving date Parity 4
CD4<- dataHSTnew4_P4 %>%
  rename(CALVING_DTM_4=CALVING_DTM) %>%
  select(COW_ID, CALVING_DTM_4)
# Calving date Parity 5
CD5<- dataHSTnew4_P5 %>%
  rename(CALVING_DTM_5=CALVING_DTM) %>%
  select(COW_ID, CALVING_DTM_5)
# Calving date Parity 6
CD6<- dataHSTnew4_P6 %>%
  rename(CALVING_DTM_6=CALVING_DTM) %>%
  select(COW_ID, CALVING_DTM_6)
# Calving date Parity 7
CD7<- dataHSTnew4_P7 %>%
  rename(CALVING_DTM_7=CALVING_DTM) %>%
  select(COW_ID, CALVING_DTM_7)
# Calving date Parity 8
CD8<- dataHSTnew4_P8 %>%
  rename(CALVING_DTM_8=CALVING_DTM) %>%
  select(COW_ID, CALVING_DTM_8)
# Calving date Parity 9
CD9<- dataHSTnew4_P9 %>%
  rename(CALVING_DTM_9=CALVING_DTM) %>%
  select(COW_ID, CALVING_DTM_9)
## Create new dataset dataHSTnew5 by adding calving date columns to dataHSTnew4
dataHSTnew5_1<- merge(dataHSTnew4, CD1, by= "COW_ID", all.x = TRUE)
dataHSTnew5_2<- merge(dataHSTnew5_1,CD2, by= "COW_ID", all.x = TRUE)
dataHSTnew5_3<-merge(dataHSTnew5_2, CD3, by= "COW_ID", all.x = TRUE)
dataHSTnew5_4<- merge(dataHSTnew5_3, CD4, by= "COW_ID", all.x = TRUE)
dataHSTnew5_5<-merge(dataHSTnew5_4, CD5, by= "COW_ID", all.x = TRUE)
dataHSTnew5_6<-merge(dataHSTnew5_5, CD6, by= "COW_ID", all.x = TRUE)   
dataHSTnew5_7<-merge(dataHSTnew5_6, CD7, by= "COW_ID", all.x = TRUE)
dataHSTnew5_8<-merge(dataHSTnew5_7, CD8, by= "COW_ID", all.x = TRUE) 
dataHSTnew5_9<-merge(dataHSTnew5_8, CD9, by= "COW_ID", all.x = TRUE)
dataHSTnew5<- dataHSTnew5_9 # New dataset with calving dates at each parity
colnames(dataHSTnew5)
# Create columns for calving intervals, eg: CI1 (Interval between first and Second calving)= Claving date for  cow_i in Pariry 2 minus (-) Calving date for cow_i in Parity 1
dataHSTnew6<- dataHSTnew5 %>%
  mutate(CI1=CALVING_DTM_2-CALVING_DTM_1, CI2=CALVING_DTM_3-CALVING_DTM_2,CI3=CALVING_DTM_4-CALVING_DTM_3, CI4=CALVING_DTM_5-CALVING_DTM_4, CI5=CALVING_DTM_6-CALVING_DTM_5,CI6=CALVING_DTM_7-CALVING_DTM_6, CI7=CALVING_DTM_8-CALVING_DTM_7, CI8=CALVING_DTM_9-CALVING_DTM_8)
## Create Column for average calving interval 

  ## Average calving interval, let's add Average calving interval to the data 
  dataHSTnew6$CI1<-as.numeric(dataHSTnew6$CI1)
dataHSTnew6$CI2<-as.numeric(dataHSTnew6$CI2)
dataHSTnew6$CI3<-as.numeric(dataHSTnew6$CI3)
dataHSTnew6$CI4<-as.numeric(dataHSTnew6$CI4)
dataHSTnew6$CI5<-as.numeric(dataHSTnew6$CI5)
dataHSTnew6$CI6<-as.numeric(dataHSTnew6$CI6)
dataHSTnew6$CI7<-as.numeric(dataHSTnew6$CI7)
dataHSTnew6$CI8<-as.numeric(dataHSTnew6$CI8)
colnames(dataHSTnew6)
AveragegCI<- rowMeans(dataHSTnew6[38:45], na.rm =TRUE)
AveragegCI
dataHSTnew6$AvgCI <- AveragegCI
dataHSTnew6 #  dataset with Average calving interval column
colnames(dataHSTnew6)
#=======================================================================
# Adding Postcode, Region and GPS Coordinates
#=======================================================================
## Adding postcodes
# Importing herd location data HST and JSE South Africa dataset
HERD_SA1 <- read.csv(file = "HST_JSE_DR_BANGA.csv", stringsAsFactors = FALSE)
glimpse(HERD_SA1)
HERD_SA1<- HERD_SA1 %>%
  rename(HERD_ID=PARTICIPANT_CODE)
HERD_SA1<- HERD_SA1[,c(1,3,2)] # Reorder
class(HERD_SA1$POSTCODE)
sum(is.na(HERD_SA1$POSTCODE))## 0 missing postcodes from raw 
HERD_SA1$POSTCODE<- as.numeric(HERD_SA1$POSTCODE)
HERD_SA1$POSTCODE[HERD_SA1$POSTCODE == ""] <- NA
HERD_SA1$POSTCODE[HERD_SA1$POSTCODE == "-"]<- NA

##Creating a new dataset dataHSTnew7 by adding the postcodes and districts to dataHSTnew6
dataHSTnew6 <- dataHSTnew6 %>%   select(HERD_ID, COW_ID,SIRE_ID, DAM_ID, BIRTH_DATE, everything())
dataHSTnew7 <- merge(dataHSTnew6, HERD_SA1, by= "HERD_ID", all.x = TRUE)
#Before merging make sure merging column are at the beginining of datasets?????
summary(dataHSTnew7$POSTCODE) ### Postcode in South Africa is between 1 and 9999. 
#  Create new dataset by removing records with postcode zero (Zimbabwe,Swaziland)

### @Thiago: I tried dataHSTnew7[!(dataHSTnew7$POSTCODE==0),]### This removed zero but 
#The change in Herd or Cow size is huge. We need to figure out. or do we set Postcode 0 as missing data??? 
summary(dataHSTnew7$POSTCODE)
#### dataHSTnew8 <- dataHSTnew7[!(dataHSTnew7$POSTCODE=="0"), ] ### I used this code to remove
# Postcode =0 but this gave me COW_ID with NAs
# and Filtering POSTCODE> 0 removes NAs. CAUTIONS!!!!

# Create a new dataset dataHSTnew8 by adding Region column based on Postcode
dataHSTnew8 <- dataHSTnew7 %>% 
  mutate(REGION=
           case_when((POSTCODE<=2899) ~"NorthenRegion",
                     (POSTCODE>=2900 & POSTCODE<4731) ~  "EasternRegion",
                     (POSTCODE>=4731 & POSTCODE<6500) ~ "SouthernRegion",
                     (POSTCODE>=6500 & POSTCODE<8300) ~ "WesternRegion",
                     (POSTCODE>=8300 & POSTCODE<10000) ~  "CentralRegion"))
summary(dataHSTnew8$MILK_YLD_305D)
length(unique(dataHSTnew8$POSTCODE)) ###590 different Postcodes
length(unique(dataHSTnew8$COW_ID))## 172374 Cows
length(unique(dataHSTnew8$HERD_ID))### 4600 herds
sum(is.na(dataHSTnew8$COW_ID))
sum(is.na(dataHSTnew8$HERD_ID))
length(dataHSTnew8$COW_ID)
#=======================================================================
# Filtering MILK_YLD_305D 
#=======================================================================
dataHSTnew8 <- dataHSTnew8 %>%   select(COW_ID,HERD_ID,SIRE_ID, DAM_ID, BIRTH_DATE, everything())
# Yinka's Method of filtering: MILK_YLD_305D betwen Mean±3SD
dataHSTnew9<- dataHSTnew8 %>%
  filter(MILK_YLD_305D >= (mean(MILK_YLD_305D)-3*sd(MILK_YLD_305D)), MILK_YLD_305D <= (mean(MILK_YLD_305D)+3*sd(MILK_YLD_305D)))
##Number of lost records
length(dataHSTnew8$COW_ID)-length(dataHSTnew9$COW_ID) ### 3541 records
##Number of lost cows 
length(unique(dataHSTnew8$COW_ID))-length(unique(dataHSTnew9$COW_ID)) ### 117cows
##Number of lost herds
length(unique(dataHSTnew8$HERD_ID))-length(unique(dataHSTnew9$HERD_ID)) ### 1 herd

# Summarry of MILK_YLD_305D
summary(dataHSTnew8$MILK_YLD_305D) 
hist(na.omit(dataHSTnew8$MILK_YLD_305D)) ##Error message: Figure margins too large
### I believed no need to apply mean±3SD. Let's proceed with dataHSTnew8.

#=======================================================================
# New Descriptive Statistics: dataHSTnew8 
#=======================================================================
#-----------------------------------------------------------------------
#Summary dataHSTnew8
#-----------------------------------------------------------------------
mean(dataHSTnew8$MILK_YLD_305D)
sd(dataHSTnew8$MILK_YLD_305D)
hist(na.omit(dataHSTnew8$MILK_YLD_305D))
## Number of records
length(dataHSTnew8$COW_ID) ##429336
## Number of cows
length(unique(dataHSTnew8$COW_ID)) #172374
## Number of herds
length(unique(dataHSTnew8$HERD_ID))# 4600
## Number of records per herd 
nrecords_byherd8<- dataHSTnew8 %>%
  group_by(HERD_ID) %>%
  summarise(nCOW_ID=n()) %>%
  mutate(cumCOW= cumsum(nCOW_ID))
summary(nrecords_byherd8$nCOW_ID)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    4.00   21.00   93.33   88.00  4921.00 
# Number of records per cow per herd
ncowherd8<- dataHSTnew8 %>%
  group_by(HERD_ID, COW_ID)%>%
  summarise(ncow_byherd= n())
summary(ncowherd8$ncow_byherd)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   2.000   2.345   3.000   9.000 

# Number of  cow per herd
ncowherd8_2<- ncowherd8 %>%
  select(HERD_ID) %>%
  summarise(ncowherd=n())
summary(ncowherd8_2$ncowherd)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.0     3.0    12.0    39.8    40.0  2233.0
### Number of postcodes
summary(dataHSTnew8$POSTCODE)
N8Postcode<- length(unique(dataHSTnew8$POSTCODE))
N8Postcode ###589, NA removed as level 
# Number of records per cow per Postcode # dataHSTnew8
ncowpostcode8<- dataHSTnew8 %>%
  group_by(POSTCODE, COW_ID)%>%
  summarise(ncow_bypostcode8= n())
summary(ncowpostcode8$ncow_bypostcode8)
# Number of  cow per postcode # dataHSTnew8
ncowpostcode8_2<- ncowpostcode8 %>%
  select(POSTCODE) %>%
  summarise(ncowbypostcode8=n()) %>%
filter(ncowbypostcode8 < 65642)  # don't count NA level
  summary(ncowpostcode8_2$ncowbypostcode8) 
  #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 #1.0     9.0    52.0   192.2   204.0  4291.0  

  summary(ncowpostcode8_2$ncowbypostcode8) 

## 65642 cows without Postcode

# Number of  herds per postcode # dataHSTnew8
nherdpostcode8<- (dataHSTnew8) %>%
  group_by(POSTCODE, HERD_ID)%>%
  summarise(nherd_bypostcode8= n())

nherdpostcode8_2<- nherdpostcode8 %>%
  select(POSTCODE) %>%
  summarise(nherdinpostcode8=n())%>%
filter(nherdinpostcode8< 2065)
summary(nherdpostcode8_2$nherdinpostcode8)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   2.000   4.304   5.000  42.000 
#-----------------------------------------------------------------------
#305D traits 
#-----------------------------------------------------------------------
#MILK_YLD_305D
# Mean, SD, median, min, max
dataHSTnew8 %>%
  summarise(Milk305D=mean(MILK_YLD_305D, na.rm = TRUE), SDMilk305D=sd(MILK_YLD_305D, na.rm = TRUE), MedianMilk305D=median(MILK_YLD_305D,na.rm = TRUE), MinMilk305D=min(MILK_YLD_305D,na.rm = TRUE), MaxMilk305D=max(MILK_YLD_305D,na.rm = TRUE)) 
# Histogram of MILK_YLD_305D
hist(na.omit(dataHSTnew8$MILK_YLD_305D))

## FAT_YLD_305D
# Mean, SD, median, min, max
dataHSTnew8 %>%
  summarise(Fat305D=mean(FAT_YLD_305D, na.rm = TRUE), SDfat305D=sd(FAT_YLD_305D, na.rm = TRUE), Medianfat305D=median(FAT_YLD_305D,na.rm = TRUE), Minfat305D=min(FAT_YLD_305D,na.rm = TRUE), Maxfat305D=max(FAT_YLD_305D,na.rm = TRUE)) 
# Histogram of FAT_YLD_305D
hist(na.omit(dataHSTnew8$FAT_YLD_305D))

## PROTEIN_YLD_305D
# Mean, SD, median, min, max
dataHSTnew8 %>%
  summarise(Prot305D=mean(PROTEIN_YLD_305D, na.rm = TRUE), SDProt305D=sd(PROTEIN_YLD_305D, na.rm = TRUE), MedianProt305D=median(PROTEIN_YLD_305D,na.rm = TRUE), MinProt305D=min(PROTEIN_YLD_305D,na.rm = TRUE), MaxProt305D=max(PROTEIN_YLD_305D,na.rm = TRUE)) 
# Histogram of PROTEIN_YLD_305D
hist(na.omit(dataHSTnew8$PROTEIN_YLD_305D))
#-----------------------------------------------------------------------  
# Reproductive traits
#-----------------------------------------------------------------------  
### Age at First Calving
# Mean, SD, median, min, max
dataHSTnew8 %>%
  summarise(MeanAFC=mean(AFC, na.rm = TRUE), SD_AFC=sd(AFC, na.rm = TRUE), MedianAFC=median(AFC,na.rm = TRUE), MinAFC=min(AFC,na.rm = TRUE), MaxAFC=max(AFC,na.rm = TRUE)) 
# Histogram of AFC
hist(na.omit(dataHSTnew8$AFC)) 

### ### Calving Interval
# Mean, SD, median, min, max,  calving interval 1 
dataHSTnew8 %>%
  summarise(MeanCI1=mean(CI1, na.rm = TRUE), SD_CI1=sd(CI1, na.rm = TRUE), MedianCI1=median(CI1,na.rm = TRUE), MinCI1=min(CI1,na.rm = TRUE), MaxCI1=max(CI1,na.rm = TRUE)) 
# Histogram of CI1
hist(na.omit(dataHSTnew8$CI1)) 

#### CI1 and AFC and AvgCI needed to be cleaned.

#=======================================================================
# Filtering Age  
# 1day= 0,0328767 month
#=======================================================================
### The minimum age/age at first calving 18 months (547days) and maximum age at Parity 9 is Mean+3SD months (3859.815days)
dataHSTnew10<-dataHSTnew8 %>% 
filter(AGE>547, AGE<= (mean(dataHSTnew8$AGE, na.rm=TRUE)+3*sd(dataHSTnew8$AGE,na.rm=TRUE))) 

##Number of lost records
length(dataHSTnew8$COW_ID)-length(dataHSTnew10$COW_ID) ###  3998 records
##Number of lost cows 
length(unique(dataHSTnew8$COW_ID))-length(unique(dataHSTnew10$COW_ID))  ### 224 cows
## Number of lost herds
length(unique(dataHSTnew8$HERD_ID))-length(unique(dataHSTnew10$HERD_ID)) ###  2 herds

#=======================================================================
# Prepare data for Blupf90
#=======================================================================
colnames(dataHSTnew10)
dataHSTnew10$REGION<- as.factor(dataHSTnew10$REGION)
dataHSTnew10$SIRE_ID<- as.factor(dataHSTnew10$SIRE_ID)
dataHSTnew10$DAM_ID<- as.factor(dataHSTnew10$DAM_ID)
class(dataHSTnew10$HERD_ID)
class(dataHSTnew10$PARITY)
class(dataHSTnew10$AGE)

datahtsablup2<- dataHSTnew10 %>%
  select(COW_ID,SIRE_ID,DAM_ID,HERD_ID,HYS,AGE,PARITY,REGION,MILK_YLD_305D)
datahtsablup2

#### datahtsablup2$code_HERD_ID = datahtsablup2$HERD_ID %>% factor() %>% as.numeric()  
colnames(datahtsablup2)
summary(datahtsablup2$MILK_YLD_305D)
summary(datahtsablup2$AGE)
datahtsablup2<- datahtsablup2 %>%
  select(COW_ID,SIRE_ID,DAM_ID,HERD_ID,HYS, AGE,PARITY,REGION, MILK_YLD_305D)
datahtsablup2
dataHSTnew10

write.table(datahtsablup2, "rawdatahstsanew.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
length(unique(dataHSTnew10$SIRE_ID))### 6816 sires

#=======================================================================
# Postcode statistics
#=======================================================================
HerdPostcode<- dataHSTnew10 %>%
  select(HERD_ID,POSTCODE)

summary(HerdPostcode)
Distinctherdpostcode<- HerdPostcode %>% distinct()
summary(Distinctherdpostcode$POSTCODE)

###Removed herd Postcode 0
HerdPostcode2<- dataHSTnew10 %>%
  select(COW_ID, HERD_ID,POSTCODE,PARITY) %>%
  filter(POSTCODE>0)
summary(HerdPostcode2)
Distinctherdpostcode2<- HerdPostcode2 %>% distinct()
summary(Distinctherdpostcode2$POSTCODE)
length(unique(HerdPostcode2$POSTCODE))
### New Number of records
length(HerdPostcode2$COW_ID)
### New Number cows
length(unique(HerdPostcode2$COW_ID))
### New Number herds
length(unique(HerdPostcode2$HERD_ID))
summary(HerdPostcode2$PARITY)
### Postcode statistics

# Number of records per cow per Postcode # HerdPostcode2
ncowpostcode10<- HerdPostcode2 %>%
  group_by(POSTCODE, COW_ID)%>%
  summarise(ncow_bypostcode10= n())
summary(ncowpostcode10$ncow_bypostcode10)

  
# Number of  cow per postcode # dataHSTnew10
ncowpostcode10_2<- ncowpostcode10 %>%
  select(POSTCODE) %>%
  summarise(ncowbypostcode10=n()) 
  summary(ncowpostcode10_2$ncowbypostcode10) 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.0     9.0    52.0   192.2   204.0  4291.0  

# Number of  herds per postcode # dataHSTnew10
nherdpostcode10<- HerdPostcode2 %>%
  group_by(POSTCODE, HERD_ID)%>%
  summarise(nherd_bypostcode10= n())

nherdpostcode10_2<- nherdpostcode10 %>%
  select(POSTCODE) %>%
  summarise(nherdinpostcode10=n())
summary(nherdpostcode10_2$nherdinpostcode10)

# Number of records per cow per herd
ncowherd11<- HerdPostcode2 %>%
  group_by(HERD_ID, COW_ID)%>%
  summarise(ncow_byherd= n())
summary(ncowherd11$ncow_byherd)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   2.000   2.345   3.000   9.000 

# Number of  cow per herd
ncowherd11_2<- ncowherd11 %>%
  select(HERD_ID) %>%
  summarise(ncowherd=n())
summary(ncowherd11_2$ncowherd)
colnames(dataHSTnew10)
summary(dataHSTnew10$DISTRICT)
sum(is.na(dataHSTnew10$DISTRICT))
dataHSTnew10$DISTRICT<-as.factor(dataHSTnew10$DISTRICT)
summary(dataHSTnew10$YEAR)

#=======================================================================
# Postcodedata. Only herds with postcode
#=======================================================================

datahtsapostcode<- dataHSTnew10 %>%
  select(COW_ID,SIRE_ID,DAM_ID,HERD_ID,HYS,AGE,PARITY,POSTCODE,MILK_YLD_305D) %>%
  filter(POSTCODE>0)
  colnames(datahtsapostcode)
  datahtsapostcode$POSTCODE<- as.factor(datahtsapostcode$POSTCODE)
  length(datahtsapostcode$COW_ID)
  length(unique(datahtsapostcode$COW_ID))
  length(unique(datahtsapostcode$HERD_ID))
  ###Import data in blupf9o format
  
  #datahtsapostcode$code_HERD_ID = datahtsapostcode$HERD_ID %>% factor() %>% as.numeric()  
  colnames(datahtsapostcode)
  summary(datahtsapostcode$MILK_YLD_305D)
  summary(datahtsapostcode$AGE)
  datahtsapostcode2<- datahtsapostcode %>%
    select(COW_ID,SIRE_ID,DAM_ID,code_HERD_ID,HYS, AGE,PARITY,POSTCODE, MILK_YLD_305D)
  datahtsapostcode2
  write.table(datahtsapostcode2, "data_base.dat",
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
  
  ## Save data datahtsapostcode and datahtsapostcode2 as .csv file
  write.csv(datahtsapostcode, "/Users/ihouaga/Documents/Project1/datahtsapostcode.csv")
  write.csv(datahtsapostcode2, "/Users/ihouaga/Documents/Project1/datahtsapostcode2.csv")
  
  # Extract the Postcodes
  
 postcodesa<- datahtsapostcode2 %>%
   select(POSTCODE)%>%
   distinct()
   postcodesa$country<- "South Africa"
  
   ### Let's split the postcodesa to size 
 
   postcodesa_1<- postcodesa[1:100,]
   postcodesa_2<- postcodesa[101:201,]
   postcodesa_3<- postcodesa[202:302,]
   postcodesa_4<- postcodesa[303:404,]
   postcodesa_5<- postcodesa[405:505,]
   postcodesa_6<- postcodesa[506:588,]
   
   write.csv(postcodesa,"/Users/ihouaga/Documents/Project1/postcodesa.csv")
   write.csv(postcodesa_1,"/Users/ihouaga/Documents/Project1/postcodesa_1.csv")
   write.csv(postcodesa_2,"/Users/ihouaga/Documents/Project1/postcodesa_2.csv")
   write.csv(postcodesa_3,"/Users/ihouaga/Documents/Project1/postcodesa_3.csv")
   write.csv(postcodesa_4,"/Users/ihouaga/Documents/Project1/postcodesa_4.csv")
   write.csv(postcodesa_5,"/Users/ihouaga/Documents/Project1/postcodesa_5.csv")
   write.csv(postcodesa_6,"/Users/ihouaga/Documents/Project1/postcodesa_6.csv")
   
   # Raw postcode with GPS
   Rawgps <- read.csv(file = "postcodesa_gps.csv", stringsAsFactors = FALSE)
   Rawgps1<- Rawgps %>%
     select(POSTCODE,Country, Latitude,Longitude)%>%
     filter(Country=="ZAF") ### 221 POSTCODE with GPS
  
    za_postal_codes<- read.csv(file = "za_postal_codes.csv", stringsAsFactors = FALSE)
   colnames(za_postal_codes)
   
   postecogps <- merge(postcodesa,  za_postal_codes, by= "POSTCODE", all.x = TRUE)
   sum(is.na(postecogps$Latitude))
   sum(is.na(postecogps$Longitude))
   
   poscodemissing<- postecogps %>%
     filter(is.na(Latitude))
  
poscoderecover<- merge(poscodemissing,  Rawgps1, by= "POSTCODE", all.x = TRUE)
poscoderecover2<- rename(poscoderecover, Latitude= Latitude.y, Longitude=Longitude.y) %>%  
  select(POSTCODE, Latitude,Longitude) 
   sum(is.na(poscoderecover2$Latitude))
   poscoderecover3<- na.omit(poscoderecover2)
   postecogps<- postecogps %>% select(POSTCODE,Latitude,Longitude) 
   postecogps2<- rbind(postecogps,poscoderecover3)
  dim(postecogps2)[1] == dim(postecogps)[1] + dim(poscoderecover3)[1]
  postecogps3 <- na.omit(postecogps2) ## 483
  allpostcode<- postecogps %>% select(POSTCODE)
  postecogps4<- merge(allpostcode,  postecogps3, by= "POSTCODE", all.x = TRUE)# Postcode dataset with missing 105 missing Latitude and Longitude
  sum(is.na(postecogps4$Latitude))
  finalmissingpostcode <-postecogps4 %>% 
  filter(is.na(Latitude)) %>% 
  select(POSTCODE)
  write.csv(postecogps3, "/Users/ihouaga/Documents/Project1/postecogps3.csv")
  write.csv(finalmissingpostcode, "/Users/ihouaga/Documents/Project1/finalmissingpostcode.csv")
  ### Import POSTCODE_ZAF_IH.csv
  oldpostcode <- read.csv(file = "POSTCODE_ZAF_IH.csv", stringsAsFactors = FALSE)
  oldpostcode<- oldpostcode %>% select(POSTCODE, X,Y)
  newrecoverpostcode<- merge(finalmissingpostcode, oldpostcode, by= "POSTCODE", all.x = TRUE) %>%
   na.omit(newrecoverpostcode) # 17 postcode recover.
  newrecoverpostcode<- newrecoverpostcode %>%
    rename(Latitude=X, Longitude=Y)

  postcode500<- rbind(postecogps3,newrecoverpostcode) 
dim(postcode500)[1]==dim(postecogps3)[1]+dim(newrecoverpostcode)[1]
write.csv(postcode500, "/Users/ihouaga/Documents/Project1/postcode500.csv")  
 ### 88 missing postcodes
postcodemissing88<- merge(allpostcode,postcode500, by="POSTCODE", all.x= TRUE)
sum(is.na(postcodemissing88))
postcodemissing88_2<- postcodemissing88 %>%
  filter(is.na(Latitude)) %>%
  select(POSTCODE)
write.csv(postcodemissing88_2, "/Users/ihouaga/Documents/Project1/postcodemissing88_2.csv") 

#=======================================================================
# datahtsapostcodegps herds with POSTCODE and GPS (503 Postcodes with GPS)
#=======================================================================

datahtsapostcodegps <- merge(postcode500,datahtsapostcode2, by="POSTCODE", all.x=TRUE) 

length(datahtsapostcodegps$COW_ID) #256208 records

length(unique(datahtsapostcodegps$COW_ID))##102676 Cows

length(unique(datahtsapostcodegps$code_HERD_ID))## 2286 herds

datahtsapostcodegps2<- datahtsapostcodegps %>% select(COW_ID,SIRE_ID,DAM_ID,code_HERD_ID,HYS, AGE,PARITY,POSTCODE, Latitude, Longitude, MILK_YLD_305D)
datahtsapostcodegps2
summary(datahtsapostcodegps2$PARITY)
### Export data for blupf90 (preclean data)

write.table(datahtsapostcodegps2, "preclean_data_base.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
####Save datahtsapostcodegps2 as .csv

write.csv(datahtsapostcodegps2,"/Users/ihouaga/Documents/Project1/datahtsapostcodegps2.csv")

###Finding the unique combinations of COW_ID and HERD_ID
  cowherd<- unique(datahtsapostcode[c("COW_ID","HERD_ID")])
  length(cowherd$COW_ID)###116337 unique cow herd. Cow that did not change herd
  length(unique(cowherd$COW_ID))
  
  ### ## Connecting region IDs to my data by using the Postcode
  
  library(sp)
  # Make POSTCODE first column
  datahtsapostcodegps3<- datahtsapostcodegps2 %>% select(POSTCODE,everything())
  df_coordinates <- datahtsapostcodegps3[,c("POSTCODE", "Latitude", "Longitude")]
  df_coordinates<- df_coordinates[!duplicated(df_coordinates$POSTCODE),]
  df_coordinates<- df_coordinates[order(df_coordinates$POSTCODE),]
  df_coordinates$post_id_nr= 1:500
  df_coordinates$Latitude <- as.numeric(df_coordinates$Latitude)
  df_coordinates$Longitude <- as.numeric(df_coordinates$Longitude)
  write.csv(df_coordinates,"/Users/ihouaga/Documents/Project1/df_coordinates.csv")
 
   df <- read.csv("df_coordinates.csv")
  df1 <- read.csv("df_coordinates.csv")
  library(sf)
  library(rgdal)
  map <- readOGR("/Users/ihouaga/Documents/Project1/MDB_District_Municipal_Boundary_2011.shp")
  coordinates(df) = ~Longitude+Latitude
  # over
  proj4string(df) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  return_dfcoordinates<- over(map, geometry(df), returnList = TRUE)
  library(reshape2)
  region_codes <- melt(return_dfcoordinates)
  colnames(region_codes) = c("post_id_nr", "region")
  df1 <- df1[,c("Latitude","Longitude", "POSTCODE", "post_id_nr")]
  df_regions <- merge(df1, region_codes, by = "post_id_nr", all.x = TRUE)
  df_regions$region <- as.numeric(df_regions$region)
  df_regions_merge <- df_regions[,c("POSTCODE", "region")]
  summary(df_regions_merge)
  # merge with data by POSTCODE
  # find missing regions (8)
  df_na <- df_regions[is.na(df_regions$region),]
  coordinates(df_na) = ~Longitude+Latitude
  proj4string(df_na) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  plot(map)
  points(df_na, col = "red")
# checked coordinates for those regions!! they are on the coast!
# They are 8 postcode with missing regions all at the Coast. Five of them (5205, 5206, 5208,5241,5259) have same coordinates.
# Adding missing region
  df_regions_merge<- df_regions_merge  %>% mutate(region = ifelse(POSTCODE == "5205" |POSTCODE == "5206"|POSTCODE == "5208"|POSTCODE == "5241"|POSTCODE == "5259", "1", region)) %>%
   mutate(region = ifelse(POSTCODE == "8130", "3", region)) %>%
   mutate(region = ifelse(POSTCODE == "7441", "2", region)) %>%
   mutate(region = ifelse(POSTCODE == "6191", "4", region))
   
   # Merge df_regions_merge with my data by postcode 
  df_regions_merge  
data_clean<- merge(datahtsapostcodegps3,df_regions_merge, by="POSTCODE", all.x=TRUE)
 sum(is.na(data_clean$region))  # 0 missing region
 sum(is.na(data_clean)) ### Zero missing data
 #=======================================================================
 # Quick summarry of data_clean
 #=======================================================================
 #Create column HY and HS in data_clean
 library(tidyr)
 data_clean$HYS2<- data_clean$HYS
 data_clean<- data_clean %>% 
   separate(HYS2, c("H", "Y", "S"))
 data_clean$HY <- paste(data_clean$H, data_clean$Y, sep= "-")
 data_clean<- data_clean %>%
   select(-c(H,Y,S))
 summary(data_clean)
 data_clean$HY <- as.factor( data_clean$HY)
 data_clean$region <- as.factor(data_clean$region)
 data_clean$code_HERD_ID<- as.factor(data_clean$code_HERD_ID)
#### Number  herd 
 sum(is.na(data_clean$region))
 n_clean_herd<- length(unique(data_clean$code_HERD_ID))
 n_clean_herd# 2286
 #### Number of records per herd 
 nherdrecords_clean<- data_clean %>%
   group_by(code_HERD_ID)%>%
   summarise(nrecordsherd=n()) 
 summary(nherdrecords_clean$nrecordsherd)
 region_records<-data_clean %>%
   select(COW_ID, region)%>%
   group_by(region)%>%
 summarise(nrecords_region=n()) 
 summary(region_records$nrecords_region)
 ### Save data_clean 
 write.csv(data_clean, "/Users/ihouaga/Documents/Project1/data_clean.csv")
data_clean_records<- merge(data_clean, nherdrecords_clean, by="code_HERD_ID", all.x=TRUE) 
data_clean_records_2<- merge(data_clean_records,region_records, by="region", all.x=TRUE)
 ##Ploting number of recordsby region as x and by herd as Y.
ggplot(aes(nrecords_region,nrecordsherd), data = data_clean_records_2) +
  geom_line()

### Filtering herd with at least Median records(24) to test the model
  data_clean_test_2 <- subset(data_clean_records_2, nrecordsherd >=24)
  summary(data_clean_test_2)
  length(unique(data_clean_test_2$code_HERD_ID))
  256208-248597 #7611 records lost after this filtering
 # Export data_clean_test_2 for blupf90
  data_clean_test_2 <- data_clean_test_2 %>% select(COW_ID,SIRE_ID,DAM_ID,code_HERD_ID,HYS, AGE,PARITY,POSTCODE, Latitude, Longitude, region, MILK_YLD_305D,HY) 
  write.table(data_clean_test_2, "data_clean_test_2.dat",
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
  summary(data_clean_test_2$MILK_YLD_305D)
 ## Subsetting data (1/3) as test sample for Thaigo
 data_test_3<- data_clean_test_2[sample(nrow(data_clean_test_2), 82865), ]
 
 save(data_test_3, file = "table.RData")

 #### Export dataset to run herd as random effect plus mean (1) as fixed effect as a test
 data_clean_test_3<- data_clean_test_2
 ## Save data clean_test 3 as .csv
 write.csv(data_clean_test_3, "/Users/ihouaga/Documents/Project1/data_clean_test_3.csv")
 
 summary(is.na(data_clean_test_2))
 # Export data clean_test 3 to bupf90
 sum(is.na(data_clean_test_3$region))
 write.table(data_clean_test_3, "data_clean_test_3.dat",
             quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
 
 
# Export final data_clean for blupf90
 summary(data_clean$region)
#data_clean_final<- data_clean%>% select(COW_ID,SIRE_ID,DAM_ID,code_HERD_ID,HYS, AGE,PARITY,POSTCODE, Latitude, Longitude, region, MILK_YLD_305D,HY, mean) 
#write.table(data_clean, "data_clean.dat",
            #quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
  
 # construct graph of neighbours
 nb.map <- poly2nb(map)
  nb2INLA("map.graph",nb.map)
  gREGION <- inla.read.graph(filename = "map.graph")
  
  # then: construct (inverse) of neighbour matrix
  R = -inla.graph2matrix(gREGION) # 52 * 52
  Rinv<- solve(R)
  # class(R) #dgTMatrix
  diag(R) = gREGION$nnbs # diagonal = nr of regions
  
  # triangular form
  Qlow <- tril(R)
  trio <- summary(Qlow) # export that - neighbours.txt file -> blupf90
  
  write.table(trio, "neighbours.txt",
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
  
  ### Data check
  
  C = inla.graph2matrix(gREGION) 
  Cinv = solve(C)  
  identical(Cinv,R) # False
  
  ### Check summary of region
  summary (data_clean_test_3$region)
  length(unique(data_clean_test_3$region)) ###46 different regions out of 52 reprensented in our data.
  sum(is.na(data_clean_test_3$region))
  ### Order of region in data
  crossref_region_id <- data_clean_test_3 %>%
    select(region)
   crossref_region_id$nr_order= 1:248597
   write.csv(crossref_region, "/Users/ihouaga/Documents/Project1/crossref_region_id.csv") 
  #=======================================================================
  # Running the models using INLA
  #=======================================================================
  
   # BASE MODEL # Keep HERD_ID or HYS in the model and not both
  
  formulabase = MILK_YLD_305D ~ 1 + HYS + AgeNestedParity  + 
    f(rowNumberAinv, model = "generic0", Cmatrix = Ainv) + 
    f(HERD_ID, model = "iid") +
    f(COW_ID, model = "iid")
  
  formulaAlt1 = MILK_YLD_305D ~ 1 + HYS + AgeNestedParity  + REGION_ID +
    f(rowNumberAinv, model = "generic0", Cmatrix = Ainv) + 
    f(HERD_ID, model = "iid") +
    f(COW_ID, model = "iid")
  
    formulaAlt2 = MILK_YLD_305D ~ 1 + HYS + AgeNestedParity  + 
    f(rowNumberAinv, model = "generic0", Cmatrix = Ainv) + 
    f(HERD_ID, model = "iid") +
    f(COW_ID, model = "iid")+
    f(REGION_ID, model = "iid")
  
    formulaAlt3 = MILK_YLD_305D ~ 1 + HYS + AgeNestedParity  + 
      f(rowNumberAinv, model = "generic0", Cmatrix = Ainv) + 
      f(HERD_ID, model = "iid") +
      f(COW_ID, model = "iid") +
      f(REGION_ID, model = "besag", graph = gREGION_ID, scale.model = FALSE)
    install.packages("spdep")
    install.packages("terra")
    library(spdep)
    library(terra)
     nb.map <- poly2nb(map)
    nb2INLA("map.graph",nb.map)
    gREGION_ID <- inla.read.graph(filename = "map.graph")
    
    load("/Users/ihouaga/Documents/Project1/region_graph.RData")
head(map)    

## Extracting region IDs
#library(raster)
#pols <- shapefile(system.file("external/lux.shp", package="raster")) 
#set.seed(20180121) 
#pts <- data.frame(coordinates(spsample(pols, 5, 'random')), name=letters[1:5])
#plot(pols); points(pts)





  
  