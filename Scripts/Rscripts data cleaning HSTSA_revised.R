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
library(readr)
library(GGally)
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
# Two columns to identify commercial and under stud animals
#=======================================================================
dataHST$Type <- str_extract_all(as.character(dataHST$PARTICIPANT),
                                 "[:upper:]", simplify = TRUE)
dataHST$Class <- paste0(dataHST$Type[, 1], dataHST$Type[, 2],
                         dataHST$Type[, 3])
dataHST$Class <- as.factor(dataHST$Class)
dataHST <- subset(dataHST,  select = -c(Type))
#-----------------------------------------------------------------------
dataHST$Herd <- readr::parse_number(as.character(dataHST$PARTICIPANT))
#=======================================================================
# Organizing the whole data set - removing NA and Zeros
dataHST  <- dataHST %>%
  filter(Class %in% c("HST", "MCM", "NA")) %>%
  filter(is.na(MILK_YLD) == FALSE) %>%
  filter(MILK_YLD > 0) %>%
  droplevels()
levels(dataHST$Class)[3] <- "Unclassified"
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

dataHST %>%
  filter(Class != "NA") %>%
  group_by(Herd, Class) %>%
  summarise(
    MeanMilk = mean(MILK_YLD, na.rm = TRUE),
    MedianMilk = median(MILK_YLD, na.rm = TRUE),
    VarMilk = var(MILK_YLD, na.rm = TRUE),
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
# You can talk on this in the lab meeting

# There three observation in the HST to be inbestigated
# observations above a Mean of 15000
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
  ggplot(aes(y = MeanMilk, x = SdMilk,  colour = PARITY)) +
  geom_text(aes(label = PARITY)) +
  xlab("SD(Milk)") +
  ylab("E(Milk)") +
  theme_bw()
# Linear relationship between mean and variance

Mean vs. Variance grouped by Class
SMilk %>%
  ggplot(aes(y = MeanMilk, x = SdMilk,  label = PARITY)) +
  geom_label(aes(fill = factor(PARITY)), colour = "white",
             fontface = "bold") +
  facet_wrap(~Class) +
  xlab("SD(Milk)") +
  ylab("E(Milk)") +
  theme_bw()
ggsave("parity.png")
# Here we can see:
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
# Centering the variable Milk
#=======================================================================
dataHST %>%
  mutate_at("PARITY", factor) %>%
  mutate(MYc = MILK_YLD-mean(MILK_YLD)) %>%
  ggplot(aes(y = MYc, x = PARITY)) +
  geom_boxplot() +
  facet_grid(~ Class) +
  xlab("Levels of Parity") +
  ylab("Centred Milk Yield") +
  theme_bw()

# Distribution is very assymetric for both HST and MCM class
# We probably have to handle with estimation problem (high leverage)
#=======================================================================
#=======================================================================
# MILK 305
#=======================================================================
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
# Two columns to identify commercial and under stud animals
#=======================================================================
dataHST$Type <- str_extract_all(as.character(dataHST$PARTICIPANT),
                                 "[:upper:]", simplify = TRUE)
dataHST$Class <- paste0(dataHST$Type[, 1], dataHST$Type[, 2],
                         dataHST$Type[, 3])
dataHST$Class <- as.factor(dataHST$Class)
dataHST <- subset(dataHST,  select = -c(Type))
#-----------------------------------------------------------------------
dataHST$Herd <- readr::parse_number(as.character(dataHST$PARTICIPANT))
#=======================================================================
# Organizing the whole data set - removing NA and Zeros
dataHST  <- dataHST %>%
  filter(Class %in% c("HST", "MCM", "NA")) %>%
  filter(is.na(MILK_YLD_305D) == FALSE) %>%
  filter(MILK_YLD_305D > 0) %>%
  droplevels()
levels(dataHST$Class)[3] <- "Unclassified"
#=======================================================================
# Descriptive Statistics
#=======================================================================
# Mean and Variance - groupd by PARITY
dataHST %>%
  group_by(PARITY) %>%
  summarise(
    MeanMilk = mean(MILK_YLD_305D, na.rm = TRUE),
    MedianMilk = median(MILK_YLD_305D, na.rm = TRUE),
    VarMilk = var(MILK_YLD_305D, na.rm = TRUE),
    SdMilk =  sqrt(VarMilk),
    nMilk = n()
    )
#-----------------------------------------------------------------------
# Mean and Variance - groupd by Class
dataHST %>%
  group_by(Class) %>%
  summarise(
    MeanMilk = mean(MILK_YLD_305D, na.rm = TRUE),
    MedianMilk = median(MILK_YLD_305D, na.rm = TRUE),
    VarMilk = var(MILK_YLD_305D, na.rm = TRUE),
    SdMilk =  sqrt(VarMilk),
    nMilk = n()
  )

dataHST %>%
  filter(Class != "NA") %>%
  group_by(Herd, Class) %>%
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
# Mean and Variance - groupd by Class and PARITY
SMilk <- dataHST %>%
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
# Almost the same as discussed to Milk_YIELD
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
# starts to decrese. This is related to the combination of effect of
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
  geom_point(shape = 1) +
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
dataHST %>%
  mutate_at("PARITY", factor) %>%
  mutate(MYc = MILK_YLD_305D-mean(MILK_YLD_305D)) %>%
  ggplot(aes(y = MYc, x = PARITY)) +
  geom_boxplot() +
  facet_grid(~ Class) +
  xlab("Levels of Parity") +
  ylab("Centred Milk 305 Yield") +
  theme_bw()

# Distribution is very assymetric for both HST and MCM class

#=======================================================================
# Pairs
#=======================================================================
ggpairs(dataHST[, c(5, 7:12)], aes(colour = factor(PARITY),
                                   alpha = 0.4))
ggsave("pairsPlot.png",  width = 14, height = 14)
ggsave("pairsPlot.png",  width = 16, height = 16)
#=======================================================================
