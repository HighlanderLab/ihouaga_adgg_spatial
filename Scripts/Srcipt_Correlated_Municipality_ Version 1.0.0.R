#=======================================================================
# Holstein (South Africa) data analysis                                #
# Isidore,  Thiago,  Ivan and Gregor                                   #
# Version 1.0.0                                                        #
# Date: 09/02/2022 
#  Municipality Spatial modelling                                                       #
#=======================================================================

#### Connect my data to municipality

### ## Connecting municipality ID to my data by using the Postcode

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

dfMu <- read.csv("df_coordinates.csv")
df1Mu <- read.csv("df_coordinates.csv")
library(sf)
library(rgdal)
mapMu <- readOGR("/Users/ihouaga/Documents/Project1/1cecc5ab-304d-494a-a29d-1b33611768702020329-1-ai5f0p.uax88.shp")
coordinates(dfMu) = ~Longitude+Latitude
# over
proj4string(dfMu) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
return_dfcoordinatesMu<- over(mapMu, geometry(dfMu), returnList = TRUE)
library(reshape2)
region_codesMu <- melt(return_dfcoordinatesMu)
colnames(region_codesMu) = c("post_id_nr", "Mu")
df1Mu <- df1Mu[,c("Latitude","Longitude", "POSTCODE", "post_id_nr")]
df_regionsMu <- merge(df1Mu, region_codesMu, by = "post_id_nr", all.x = TRUE)
df_regionsMu$Mu <- as.numeric(df_regionsMu$Mu)
df_regions_mergeMu <- df_regionsMu[,c("POSTCODE", "Mu")]
summary(df_regions_mergeMu$Mu)
# merge with data by POSTCODE
# find missing Mu (8)
df_naMu <- df_regionsMu[is.na(df_regionsMu$Mu),]
coordinates(df_naMu) = ~Longitude+Latitude
proj4string(df_naMu) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
plot(mapMu)
points(df_naMu, col = "red")


mapMu2 <- st_read("/Users/ihouaga/Documents/Project1/1cecc5ab-304d-494a-a29d-1b33611768702020329-1-ai5f0p.uax88.shp", quiet = TRUE)
class(mapMu2) # class sf
ggplot(mapMu2) + geom_sf()+ geom_sf_text(data = mapMu2, aes(label = OBJECTID), colour = "black")
ggplot(mapMu2) + geom_sf()+ geom_sf_text(data = mapMu2, aes(label = OBJECTID), size = 3, colour = "black")
ggplot(mapMu2) + geom_sf()+ geom_sf_text(data = mapMu2, aes(label = OBJECTID), size = 1, colour = "black")
ggplot(mapMu2) + geom_sf()+ geom_sf_text(data = mapMu2, aes(label = OBJECTID), size = 0.5, colour = "black")

mapMu2 <- st_read("/Users/ihouaga/Documents/Project1/1cecc5ab-304d-494a-a29d-1b33611768702020329-1-ai5f0p.uax88.shp", quiet = TRUE)

# checked coordinates for those Municipalities!! they are on the coast!
# They are 8 postcode with missing municipality all at the Coast. Five of them (5205, 5206, 5208,5241,5259) have same coordinates.
# Adding missing Municipality
df_regions_mergeMu<- df_regions_mergeMu  %>% mutate(Mu = ifelse(POSTCODE == "5205" |POSTCODE == "5206"|POSTCODE == "5208"|POSTCODE == "5241"|POSTCODE == "5259", "180", Mu)) %>%
  mutate(Mu = ifelse(POSTCODE == "8130", "3", Mu)) %>%
  mutate(Mu = ifelse(POSTCODE == "7441", "198", Mu)) %>%
  mutate(Mu = ifelse(POSTCODE == "6191", "132", Mu))

# Merge df_regions_mergeMu with my data by postcode 
df_regions_mergeMu  
summary (df_regions_mergeMu$Mu)
df_regions_mergeMu$Mu <- as.factor(df_regions_mergeMu$Mu)

data_cleanMu<- merge(datahtsapostcodegps3,df_regions_mergeMu, by="POSTCODE", all.x=TRUE)
sum(is.na(data_cleanMu$Mu))  # 0 missing region
sum(is.na(data_cleanMu)) ### Zero missing data
summary(data_cleanMu$Mu)
length(unique(data_cleanMu$Mu)) ##146 Local Municipalities
length(unique(mapMu2$OBJECTID)) ## 213 Local Municipalities in South Africa
#=======================================================================
# Quick summarry of data_cleanMu
#=======================================================================
#Create column HY and HS in data_clean
data_cleanMu$HYS2<- data_clean$HYS
data_cleanMu<- data_cleanMu %>% 
  separate(HYS2, c("H", "Y", "S"))
data_cleanMu$HY <- paste(data_cleanMu$H, data_cleanMu$Y, sep= "-")
data_cleanMu<- data_cleanMu %>%
  select(-c(H,Y,S))
summary(data_cleanMu)
data_cleanMu$HY <- as.factor( data_cleanMu$HY)
data_cleanMu$Mu <- as.factor(data_cleanMu$Mu)
data_cleanMu$code_HERD_ID<- as.factor(data_cleanMu$code_HERD_ID)

#### Number  herd 
sum(is.na(data_cleanMu$Mu))
n_clean_herdMu<- length(unique(data_cleanMu$code_HERD_ID))
n_clean_herdMu# 2286
#### Number of records per herd 
nherdrecords_cleanMu<- data_cleanMu %>%
  group_by(code_HERD_ID)%>%
  summarise(nrecordsherdMu=n()) 
summary(nherdrecords_cleanMu$nrecordsherdMu)# Records by herd
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.0     4.0    24.0   112.1    96.0  4914.0 
Mu_records<-data_cleanMu %>%
  select(COW_ID, Mu)%>%
  group_by(Mu)%>%
  summarise(nrecords_Mu=n()) 
summary(Mu_records$nrecords_Mu) #Records by municipality
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.0   108.8   827.0  1754.8  1731.5 20222.0
### Save data_cleanMu 
write.csv(data_cleanMu, "/Users/ihouaga/Documents/Project1/data_cleanMu.csv")

##Make code_herd_ID first column in data_cleanMu 
data_cleanMu<- data_cleanMu  %>% select(code_HERD_ID,everything())
data_clean_recordsMu<- merge(data_cleanMu, nherdrecords_cleanMu, by="code_HERD_ID", all.x=TRUE) 
data_clean_records_2Mu<- merge(data_clean_recordsMu,Mu_records, by="Mu", all.x=TRUE)

##Ploting number of recordsby region as x and by herd as Y.
ggplot(aes(nrecords_Mu,nrecordsherdMu), data = data_clean_records_2Mu) +
  geom_line()

### Filtering herd with at least Median records(24) to test the model
data_clean_test_2Mu <- subset(data_clean_recordsMu, nrecordsherdMu >=24)
summary(data_clean_test_2Mu)
length(unique(data_clean_test_2Mu$code_HERD_ID)) #1145
length(data_clean_test_2Mu$COW_ID)#248597
256208-248597 #7611 records lost after this filtering
# Export data_clean_test_2Mu for blupf90

data_clean_test_2Mu <- data_clean_test_2Mu %>% select(COW_ID,SIRE_ID,DAM_ID,code_HERD_ID,HYS, AGE,PARITY,POSTCODE, Latitude, Longitude, Mu, MILK_YLD_305D,HY) 
write.table(data_clean_test_2Mu, "data_clean_test_2Mu.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
summary(data_clean_test_2Mu$MILK_YLD_305D)
summary(data_clean_test_2Mu$Mu)
#### Export dataset to run herd as random effect plus mean (1) as fixed effect as a test
#It was advised to add a mean (constant to the model to capture intercept)
data_clean_finalMu <- data_clean_test_2Mu 
data_clean_finalMu$mean <- rep(1, nrow(data_clean_test_2Mu))
summary(data_clean_finalMu)

## Save data clean_test_3Mu and data_clean_finalMu as .csv
write.csv(data_clean_finalMu, "/Users/ihouaga/Documents/Project1/data_clean_finalMu.csv")

sum(is.na(data_clean_finalMu))
# Export data clean_test3 to bupf90
sum(is.na(data_clean_test_3$Mu))
write.table(data_clean_test_3Mu, "data_clean_test_3Mu.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")


# Export final data_clean_finalMu for blupf90
summary(data_clean_finalMu$Mu)
sum(is.na(data_clean_finalMu$Mu))# 0 missing data for municipality
sum(is.na(data_clean_finalMu)) #0 missing data in the dataset
data_clean_finalMu$Mu<- as.factor(data_clean_finalMu$Mu)

#Export_data_cleanMu used for analysis 
write.table(data_clean_finalMu, "data_clean_finalMu.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")

# construct graph of neighbours
nb.mapMu2 <- poly2nb(mapMu2)
nb2INLA("map.graphMu2",nb.mapMu2)
gMu <- inla.read.graph(filename = "map.graphMu2")

# then: construct (inverse) of neighbour matrix
RMu = -inla.graph2matrix(gMu) # 213 * 213
RMuinv<- solve(RMu)
RMuinv
RMu
diag(RMu) = gMu$nnbs # diagonal = nr of municipalities
diag(RMu)
# reinstall and upgrade inla
devtools::install_github(repo = "https://github.com/hrue/r-inla", ref = "stable", subdir = "rinla", build = FALSE)
(.packages())

#detach(package:INLA)# to unload a package
#Scaling RMu matrix
save(RMu, file = "RMu.RData")
scaledRMu = inla.scale.model(RMu, constr = list(A = matrix(1, 1, 213), e=0))
save(scaledRMu, file = "scaledRMu.RData")

# triangular form
QlowMu <- tril(scaledRMu)
trioMu <- summary(QlowMu) # export that - neighbours_scaledMu.txt file -> blupf90
summary(trioMu$i)
write.table(trioMu, "neighbours_scaledMu.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")


### Check summary of Municipality
summary (data_clean_test_3Mu$Mu)
length(unique(data_clean_test_3$region)) ###126 different municipality out of 213 reprensented in our data. Use the crossref of region in the renf90.table in blupf90 to match the each region with their renumbered id.
sum(is.na(data_clean_test_3$region))


#=======================================================================
# data_clean_fullMu Without filtering
#=======================================================================

data_clean_fullMu <- data_clean_recordsMu %>% select(COW_ID,SIRE_ID,DAM_ID,code_HERD_ID,HYS, AGE,PARITY,POSTCODE, Latitude, Longitude, Mu, MILK_YLD_305D,HY) 
### Add the mean=1 (intercept)
data_clean_fullMu$mean <- rep(1, nrow(data_clean_fullMu))
summary(data_clean_fullMu)
write.table(data_clean_fullMu, "data_clean_fullMu.dat",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")

