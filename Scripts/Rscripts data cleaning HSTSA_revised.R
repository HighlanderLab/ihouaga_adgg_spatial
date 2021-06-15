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
