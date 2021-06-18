#=======================================================================
# Holstein data analysis                                               #
# Isidore,  Thiago,  Ivan and Gregor                                   #
# Version 1.0.0                                                        #
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
library(lme4)
#=======================================================================
# Uploading datasets
#=======================================================================
HSTSA1 <- read.csv(
  file = "./MACE-Data/Holstein_South Africa/HSTSA_01.csv",
  stringsAsFactors = FALSE)

# Importing HST_02 dataset
HSTSA2 <- read.csv(
  file = "./MACE-Data/Holstein_South Africa/HSTSA_02.csv",
  stringsAsFactors = FALSE)

#=======================================================================
# Data organization
#=======================================================================
dataHST1 <- HSTSA1 %>%
  mutate_at(c("ANIMAL_NUMBER", "PARTICIPANT", "FARMERS_NO"),  factor) %>%
  filter(PARITY < 10) %>%
  droplevels()

dataHST2 <- HSTSA2 %>%
  mutate_at(c("ANIMAL_NUMBER", "PARTICIPANT", "FARMERS_NO"),  factor) %>%
  filter(PARITY < 10) %>%
  droplevels()

glimpse(dataHST1)
glimpse(dataHST2)
#=======================================================================
# Stacking the data
#=======================================================================
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
#=======================================================================
# Factor PARTICIPANT -> Unknown level to NA
#=======================================================================
dataHST$PARTICIPANT[dataHST$PARTICIPANT == "Unknown"] <- NA
#=======================================================================
# Two columns to identify commercial and under study animals
#=======================================================================
dataHST$Type <- str_extract_all(as.character(dataHST$PARTICIPANT),
                                 "[:upper:]", simplify = TRUE)
dataHST$Class <- paste0(dataHST$Type[, 1], dataHST$Type[, 2],
                         dataHST$Type[, 3])
dataHST$Class <- as.factor(dataHST$Class)
dataHST <- subset(dataHST,  select = -c(Type))
#=======================================================================
# Summary
#=======================================================================
summary(dataHST)
#=======================================================================
# Descriptive Statistics
#=======================================================================
# Mean and Variance - groupd by PARITY
dataHST %>%
  group_by(PARITY) %>%
  summarise(
    MeanMilk = mean(MILK_YLD, na.rm = TRUE),
    MedianMilk = median(MILK_YLD, na.rm = TRUE),
    VarMilk = var(MILK_YLD, na.rm = TRUE),
    SdMilk =  sqrt(VarMilk),
    nMilk = n()
    )
#-----------------------------------------------------------------------
# Mean and Variance - groupd by Class
dataHST %>%
  group_by(Class) %>%
  summarise(
    MeanMilk = mean(MILK_YLD, na.rm = TRUE),
    MedianMilk = median(MILK_YLD, na.rm = TRUE),
    VarMilk = var(MILK_YLD, na.rm = TRUE),
    SdMilk =  sqrt(VarMilk),
    nMilk = n()
  )
#-----------------------------------------------------------------------
# Mean and Variance - groupd by Class and PARITY
SMilk <- dataHST %>%
  group_by(PARITY, Class) %>%
  summarise(
    MeanMilk = mean(MILK_YLD, na.rm = TRUE),
    MedianMilk = median(MILK_YLD, na.rm = TRUE),
    VarMilk = var(MILK_YLD, na.rm = TRUE),
    SdMilk =  sqrt(VarMilk),
    nMilk = n()
  )

# Mean vs. Variance
SMilk %>%
  ggplot(aes(y = MeanMilk, x = SdMilk,  colour = PARITY,
             size = PARITY)) +
  geom_point(shape = 1) +
  xlab("SD(Milk)") +
  ylab("E(Milk)") +
  theme_bw()
# Linear relationship between mean and variance

Mean vs. Variance grouped by Class
SMilk %>%
  ggplot(aes(y = MeanMilk, x = SdMilk,  colour = PARITY,
             size = PARITY)) +
  geom_point(shape = 1) +
  facet_wrap(~Class) +
  xlab("SD(Milk)") +
  ylab("E(Milk)") +
  theme_bw()
# Linear relationship between mean and variance -> HST and MCM
# NAs group has variance more concentrate around 1900 and mean 3800
# Few observations for JSE and UNK - remove from the dataset?
#-----------------------------------------------------------------------
# Mean vs. number of observations
SMilk %>%
  ggplot(aes(y = MeanMilk, x = nMilk,  colour = PARITY,
             size = PARITY)) +
  geom_point(shape = 1) +
  facet_wrap(~Class) +
  xlab("Number of Observations") +
  ylab("E(Milk)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# We can see the Milk production increases up to PARITY 3-4 and then
# starts to decrese. This is related to the combination of effect of
# PARITY and number of observations

# SD vs. number of observations
SMilk %>%
  ggplot(aes(y = SdMilk, x = nMilk,  colour = PARITY,
             size = PARITY)) +
  geom_point(shape = 1) +
  facet_wrap(~Class) +
  xlab("Number of Observations") +
  ylab("SD(Milk)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Same for Standard deviation
#-----------------------------------------------------------------------
# Mean and Median relationship
SMilk %>%
  ggplot(aes(y = MeanMilk,  x = MedianMilk)) +
  geom_point(shape = 1) +
  facet_wrap(~Class) +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("Median(Milk)") +
  ylab("E(Milk)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Ok -> linear relationship
#=======================================================================
# BoxPlots
#=======================================================================
dataHST %>%
  mutate_at("PARITY", factor) %>%
  ggplot(aes(y = MILK_YLD, x = PARITY)) +
  geom_boxplot(na.rm = TRUE) +
  facet_grid(~ Class) +
  theme_bw()

# Distribution is very assymetric for both HST and MCM class
# We probably have to handle with estimation problem (high leverage)

#NA has a mixture of animals from HST, MCM, UNK and JSE -> also
# assymetric
#=======================================================================
# Testing basic initial models

# Linear regression
m1 <- lmer(MILK_YLD ~ PARITY + Class + (1 | PARTICIPANT) +
             (1 | PARTICIPANT:ANIMAL_NUMBER), data = dataHST)

# Gamma distribution
m1 <- glmer(MILK_YLD ~ PARITY + Class + (1 | PARTICIPANT) +
              (1 | PARTICIPANT:ANIMAL_NUMBER),  family = "Gamma",
            data = dataHST)
#=======================================================================
