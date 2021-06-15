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
# Organizing date variables
#=======================================================================
dataHST1$CALVING_DTM <- as.Date(dataHST1$CALVING_DTM,
                                format = "%d/%m/%Y")
dataHST1$LACTATION_DTM <- as.Date(dataHST1$LACTATION_DTM,
                                  format = "%d/%m/%Y")
dataHST1$END_LACTATION_DTM <- as.Date(dataHST1$END_LACTATION_DTM,
                                      format = "%d/%m/%Y")
#=======================================================================
# Factor PARTICIPANT -> Unknown level to NA
#=======================================================================
dataHST1$PARTICIPANT[dataHST1$PARTICIPANT == "Unknown"] <- NA

#=======================================================================
# Two columns to identify commercial and under study animals
#=======================================================================
dataHST1$Type <- str_extract_all(as.character(dataHST1$PARTICIPANT),
                                 "[:upper:]", simplify = TRUE)
dataHST1$Class <- paste0(dataHST1$Type[, 1], dataHST1$Type[, 2],
                         dataHST1$Type[, 3])
dataHST1$Class <- as.factor(dataHST1$Class)
dataHST1 <- subset(dataHST1,  select = -c(Type))
#=======================================================================
