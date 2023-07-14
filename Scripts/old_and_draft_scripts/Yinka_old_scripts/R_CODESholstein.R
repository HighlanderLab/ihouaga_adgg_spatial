total4 <- read.csv("D:\\Users\\OAbejide\\305D\\total4.csv")
total4[1:2,]
is.factor(total4$HerdID)
total4$HerdID <- factor(total4$HerdID)
total4 <- read.csv("D:\\Users\\OAbejide\\305D\\total4.csv")
is.factor(total4$Herd_ID)
total4$Herd_ID[1:5]
total4$Herd_ID <- factor(as.character(total4$Herd_ID))
summary(total4$Herd_ID, maxsum=100000)
summary(total4$Herd_ID, maxsum=100000)
total4[1,
]
ncows <- tapply(total4$CowID, total4$Herd_ID, function(x){length(unique(x))})
ncows
ncows <- tapply(total4$CowID, total4$Herd_ID, function(x){length(unique(x))})
names(ncows)[ncows > 50]
herdkeep <- names(ncows)[ncows > 50]
length(herdkeep)
herdkeep <- names(ncows)[ncows > -1]
length(herdkeep)
herdkeep <- names(ncows)[ncows > 50]
herdkeep
match(herdkeep, total4$Herd_ID)
match(total4$Herd_ID, herdkeep)
recordkeep <- ! is.na(match(total4$Herd_ID, herdkeep))
recordkeep
sum(recordkeep)
sum(! recordkeep)
total4 <- total4[recordkeep,]
dim(total4)
total4[1:5,]
wetmonths <- c(4,5,6,11)
total$Month <- total$LSDseason
total4$Month <- total4$LSDseason
total4$WetSeason <- ! is.na(match(total4$Month, wetmonths))
total4[1:5,]
total4[100:105,]
total4$WetSeason <- 2 - (! is.na(match(total4$Month, wetmonths)))
total4[100:105,]
total4$Season <- 2 - (! is.na(match(total4$Month, wetmonths)))
total4$WetSeason <- NULL
total4$Season[1:5,]
total4[1:5,]
total4[100:105,]
write.csv(total4, file="D:\\Users\\OAbejide\\305D\\total4-modified.csv")
total4[1:2,]
total4$Combi <- factor(paste(total4$Herd_ID, total4$LSDyear, total4$Season))
ncow.combi <- tapply(total4$CowID, total4$Combi, function(x){length(unique(x))})
ncow.combi
ncow.combi <- tapply(total4$CowID, total4$Combi, function(x){length(unique(x))})
ncow.combi[1:5,]
ncow.dat <- data.frame(Combi = names(ncow.combi), Ncow = as.numeric(ncow.combi))
write.csv(ncow.dat, file="D:\\Users\\OAbejide\\305D\\ncow.combi.csv")
combikeep <- names(ncow.combi)[ncow.combi > 3]
sum(combikeep)
recordkeep.combi <- ! is.na(match(total4$Combi, combikeep))
sum(recordkeep.combi)
sum(! recordkeep.combi)
total4 <- total4[recordkeep.combi,]
write.csv(total4, file="D:\\Users\\OAbejide\\305D\\total4-modified2.csv")
head(total4)
newest <- data.frame(CowID = tapply(total4$CowID, total4$CowID, function(x){x[1]}))
newest$Herd_ID <- tapply(total4$Season, total4$Herd_ID, function(x){x[1]})
newest$Herd_ID <- tapply(total4$Herd_ID, total4$Cow_ID, function(x){x[1]})
newest <- data.frame(CowID = tapply(total4$CowID, total4$CowID, function(x){x[1]}))
newest$Herd_ID <- tapply(total4$Herd_ID, total4$Cow_ID, function(x){x[1]})
newest <- data.frame(CowID = tapply(total4$CowID, total4$CowID, function(x){x[1]}))
newest$Herd_ID <- tapply(total4$Herd_ID, total4$Cow_ID, function(x){x[1]})
newest <- data.frame(CowID = tapply(total4$CowID, total4$CowID, function(x){x[1]}))
newest$Herd_ID <- tapply(total4$Herd_ID, total4$CowID, function(x){x[1]})
newest <- data.frame(CowID = tapply(total4$CowID, total4$CowID, function(x){x[1]}))
newest$Herd_ID <- tapply(total4$Herd_ID, total4$CowID, function(x){x[1]})
newest$Season <- tapply(total4$Season, total4$CowID, function(x){x[1]})
newest$LSDyear <- tapply(total4$LSDyear, total4$CowID, function(x){x[1]})
newest[1:5,]
dim(newest)
newest$Combi <- tapply(total4$Combi, total4$CowID, function(x){x[1]})
newest[1:2,]
newest <- newest[order(newest$Combi),]
newest[1:2,]
newest <- newest[,c(2,3,4,1)]
newest[1:2,]
newest[1:10,]
row.names(newest)
row.names(newest) <- NULL
newest[1:5,]
newest <- newest[,c(4,1:3)]
newest[1:5,]
write.csv(newest, file="D:\\Users\\OAbejide\\305D\\cow-id-combination.csv")
clean <- subset(total4-modified2, !duplicated(total4-modified2$CowID))
clean <- subset(total4-modified2, !duplicated(total4-modified$CowID))
clean <- subset(total4-modified, !duplicated(total4-modified$CowID))
clean <- subset(total4-modified, !duplicated(total4))
clean <- subset(total4, !duplicated(total4$CowID))
head(clean)
write.csv(clean,"clean.csv")
write.csv(clean,file="D:\\Users\\OAbejide\\305D\\clean.csv")
length(unique(clean$CowID))
length(unique(clean$Combi))
HYS <- subset(clean, K305D > 0) # showed all +ve values
HY <- subset(clean, K305D < 0) # showed all +ve values
head(HY)
write.csv(total4, file="D:\\Users\\OAbejide\\305D\\total4-modified2.csv")
write.csv(HYS, file="D:\\Users\\OAbejide\\305D\\HYS.csv")
length(unique(HYS$CowID))
length(unique(HY$CowID))
write.csv(HY, file="D:\\Users\\OAbejide\\305D\\HY.csv")
length(unique(HYS$combi))
new<- HYS(combi,sep='')
new<- paste(HYS, combi,sep='')
new<- paste(HYS, HYS$combi,sep='')
head(new)
new<- paste(HYS, HYS$combi,sep='')
head(new)
head(HYS)
concatenated<-read.csv("Concatenated-HYS.csv")
concatenated<-read.csv("Concatenated-HYS.csv")
head(concatenated)
summary(concatenated)
length(unique(Contenated$CowID))
length(unique(Concatenated$CowID))
length(unique(concatenated$CowID))
length(unique(concatenated$Herd_ID))
str(concatenated)
summary(concatenated$HYS)
range(concatenated$HYS)
head(d2)
head(LactSeason)
findunique <- function(x){length(unique(x))}
y <- tapply(concatenated$CowID, concatenated$HYS, findunique)
summary(factor(y))
y <- tapply(concatenated$HYS, concatenated$CowID, findunique)
summary(factor(y))
NGTH(Y)LE
length(y)
y <- tapply(concatenated$CowID, concatenated$HYS, findunique)
length(y)
y <- tapply(concatenated$CowID, concatenated$Herd_ID, findunique)
summary(factor(y))
y <- tapply(concatenated$Herd_ID, concatenated$CowID, findunique)
summary(factor(y))
n<-concatenated[,CowID, HYS]
n<-concatenated[CowID, HYS]
new <- concatenated[! is.na(match(concatenated$CowID, names(q[q == 117]))),] 
new <- concatenated[! is.na(match(concatenated$CowID, names(q[q == 1]))),] 
concatenated$Herd_ID <- factor(as.character(concatenated$Herd_ID))
concatenated$HYS <- factor(as.character(concatenated$HYS))
findunique <- function(x){length(unique(x))}
y <- tapply(concatenated$CowID, concatenated$HYS, findunique)
summary(factor(y))
data<-read.txt("HST.SA.txt",sep='\t',header=T)
data<-read.table("HST.SA.txt",sep='\t',header=T)
data<-read.table("HST.SA.txt",sep='\t',header=T)
data<-read.table("HST_SA.txt",sep='\t',header=T)
head(data)
data<-read.csv("HST_SA.csv")
head(data)
data<-read.csv("HST_SA.csv")
head(data)
is.na(data)
holstein<-na.omit(data)
write.csv("holstein.csv")
write.csv(holstein,"holstein.csv")
length(unique(holstein$FARMER_NO)
length(unique(holstein$FARMER_NO))
head(holstein)
length(unique(holstein$FARMERS_NO))
length(unique(holstein$ANIMAL_NUMBER))
range(holstein$PARITY)
dd<-subset(holstein, holstein$PARITY=="30")
head(dd)
hist(holstein$MILK_YLD_305)
summary(holstein$MILK_YLD_305)
summary(holstein$MILK_YLD)
hist(holstein$MILK_YLD)
hist(holstein$PARITY)
summary(holstein$PARITY)
length(unique(holstein$PARTICIPANT))
sd(holstein$MILK_YLD)
sd(holstein$MILK_YLD_305)
sd(holstein$PARITY)
summary(holstein$MILK_YLD)
summary(holstein$PARITY)
sd(holstein$PARTICIPANT)
summary(holstein$PARTICIPANT)
summary(holstein$FARMERS_NO)
sd(holstein$MILK_YLD_305)
hist(holstein$MILK_YLD_305)
PARITY1<-subset(holstein,holstein$PARITY=="1")
write.csv(PARITY1,"PARITY1.csv")
PARITY2<-subset(holstein,holstein$PARITY=="2")
write.csv(PARITY2,"PARITY2.csv")
PARITY3<-subset(holstein,holstein$PARITY=="3")
write.csv(PARITY3,"PARITY3.csv")
PARITY4<-subset(holstein,holstein$PARITY=="4")
write.csv(PARITY4,"PARITY4.csv")
PARITY5<-subset(holstein,holstein$PARITY=="5")
write.csv(PARITY5,"PARITY5.csv")
length(unique(PARITY1$ANIMAL_NUMBER))
hist(PARITY1$MILK_YLD_305)
length(unique(PARITY2$ANIMAL_NUMBER))
hist(PARITY2$MILK_YLD_305)
length(unique(PARITY3$ANIMAL_NUMBER))
hist(PARITY3$MILK_YLD_305)
length(unique(PARITY4$ANIMAL_NUMBER))
length(unique(PARITY5$ANIMAL_NUMBER))
hist(PARITY5$MILK_YLD_305)
hist(PARITY4$MILK_YLD_305)
hist(PARITY3$MILK_YLD_305)
Holstein<-complete.cases(holstein)
write.csv(Holstein,"Holstein.csv")
PARITY1(Holstein,"Holstein$Parity=="1")
PARITY1(Holstein,Holstein$Parity=="1")
PARITY1<-(Holstein,Holstein$Parity=="1")
PARITY1<-subset(Holstein,Holstein$Parity=="1")
PARITY1<-subset(Holstein,Parity=="1")
PARITY1<-subset(Holstein,Holstein$Parity=="1")
PARITY1<-subset(Holstein,Holstein$PARITY=="1")
head(Holstein)
g<-holstein[!complete.cases(holstein),]
head(g)
is,na(holstein)
is.na(holstein)
g<-holstein[!complete.cases(holstein)]
head(g)
g<-holstein[!complete.cases(holstein),]
head(g)
g<-complete.cases(holstein) 
head(g)
stopifnot(complete.cases(y) != is.na(y))
g<-stopifnot(complete.cases(holstein) != is.na(holstein))
head(g)
c<-holstein[complete.cases(holstein), ]
head(c)
Holstein<-holstein[complete.cases(holstein), ]
write.csv(Holstein,"Holstein.csv")
PARITY1<-subset(PARITY, Holstein$PARITY==”1”)
PARITY1<-subset(Holstein, Holstein$PARITY==”1”)
PARITY1<-subset(Holstein, Holstein$PARITY=="1")
write.csv(PARITY1,"PARITY1.csv")
PARITY2<-subset(Holstein, Holstein$PARITY=="2")
write.csv(PARITY2,"PARITY2.csv")
PARITY3<-subset(Holstein, Holstein$PARITY=="3")
write.csv(PARITY3,"PARITY3.csv")
PARITY4<-subset(Holstein, Holstein$PARITY=="4")
write.csv(PARITY4,"PARITY4.csv")
PARITY5<-subset(Holstein, Holstein$PARITY=="5")
write.csv(PARITY5,"PARITY5.csv")
write.csv(Holstein,"Holstein.csv")
write.csv(Holstein,"Holstein.csv")
write.csv(Holstein,"Y://Country Files//SOUTH AFRICA//HOLSTEIN//Holstein.csv")
write.csv(Holstein,"Y:\\Country Files\\SOUTH AFRICA\\HOLSTEIN\\Holstein.csv")
length(Holstein$CowID)
length(Holstein$ANIMAL_NUMBER)

####################################################
SORTING PEDIGREE( original script from Gorjanc)
The sequence of commands should be:
pedOrig <- read.table("blabla.txt")
ord <- order(orderPed(pedOrig))
ped <- pedOrig[ord, ]
################################
ord <- order(orderPed(hst2))
ped <- hst2[ord, ]
hist(ped$gen, main="South African Holstein-Friesian generations",col="orange",ylab="Holstein-Friesians",xlab="Generations (24)")
###### find all males (complete cases) with international ids
IDS<-hstmales[complete.cases(hstmales$INTNL_ID),]
#########################
ord1 <- order(orderPed(jer))
ped1 <- jer[ord1, ]
hist(ped1$gen, main="South African Holstein-Friesian generations",col="orange",ylab="Holstein-Friesians",xlab="Generations (28)")
################################
KEEP HERDS WITH 50 OR MORE ANIMALS
dat$Herd_ID <- factor(as.character(total4$Herd_ID))
ncows <- tapply(dat$ANIMAL_NUMBER, dat$PARTICIPANT, function(x){length(unique(x))})
herdkeep <- names(ncows)[ncows >= 50]
***   new<-dat[which(dat$PARTICIPANT %in% herdkeep),] *****
new<-dat[which(dat$PARTICIPANT %in% herdkeep),] ######## TO CHECK WHICH MATCH AND TO PULL OUT RECORDS FROM ORIGINAL FILE FOR THE HERDS LISTED BY NAMES
##################################
new$CALVING_DTM <- as.Date(as.character(new$CALVING_DTM), "%d/%m/%Y")
 new$LACTATION_DTM <- as.Date(as.character(Data$LACTATION_DTM), "%d/%m/%Y")

 HOLSTEIN_50<- new_dates[order(new_dates$ANIMAL_NUMBER, new_dates$BIRTH_DATE),]
 ########################################### BASED ON 27TH JULY, 2016 (Mostert et al, 2006, Kgole et al, 2012)
 tmp$C_SEASON <- factor(as.character(tmp$CALVING_MTH), levels=1:12) season of birth for using calving date
levels(tmp$C_SEASON) <- c(2,2,2,1,1,1,1,1,1,2,2,2)
season 1 (winter): Apr, May, Jun, Jul, Aug, Sep
season 2 (dry): Oct, Nov, Dec, Jan, Feb, March 
###############################
EBVs and reliability for sires in South Africa
SIRE821_MYEBV$Rel_MY <- 1 - (SIRE821_MYEBV$V4^2 / (186531)) #305 day MY
SIRE821_MYEBV$Rel_L1 <- 1 - (SIRE821_MYEBV$V4^2 / (178970)) #1st Lactation MY 178970
sire$Rel_AFC <- 1 - (sire$AFC_se^2 / (1154.91 )) #AFC 1154.91 
sire$Rel_CI1 <- 1 - (sire$V4^2 / (151.350)) #CI1 151.350

sire$Rel_CI1[sire$Rel_CI1 <= 0] <- 0.01 #as min was -0.08760
sire$Rel_CI1 <- 0.01 



sln<-read.table("D:\\Users\\OAbejide\\HOLSTEIN\HF_myL1st.sln", fill=TRUE, head=FALSE) #  exporting large file in .sln to R



################################# 10/more HYS combinations
tmp$C_HYS<-as.factor(tmp$C_HYS)                                    #Convert as factor
niv<-levels(tmp$C_HYS)                                           #List of levels of the factor
extracted<-as.data.frame(which(table(tmp$C_HYS)>=10))[,1]         #List of levels that appear three or more times
rep_HYS<-vector(length=length(extracted))                      #Create vector to store HYS levels that appear three or more times
for (i in 1:length(extracted)){                                 #Loop to go through the entire list of repeated levels                       
  rep_HYS[i]<-niv[extracted[i]]                                #The levels we want is the one which is: niv[listed level in extracted]
}

#rep_HYS                                                        #List of herds that appear three or more times

HFSA_HYS10<-merge(tmp,rep_HYS,by="C_HYS",by.y=TRUE)	
#data4 <-match(khf,cow by.x=c("LactNo"),by.y=c("LactNo"))  #match function
 ########################################################################################################################
P1$CALVING_AGE[P1$CALVING_AGE > 1188 & ! is.na(P1$CALVING_AGE)] <- NA
P1$CALVING_AGE[P1$CALVING_AGE < 466 & ! is.na(P1$CALVING_AGE)] <- NA
summary(P1$CALVING_AGE)
p1<-na.omit(P1)

P2$CALVING_AGE[P2$CALVING_AGE > 1760 & ! is.na(P2$CALVING_AGE)] <- NA
P2$CALVING_AGE[P2$CALVING_AGE < 762 & ! is.na(P2$CALVING_AGE)] <- NA
summary(P2$CALVING_AGE)

P3$CALVING_AGE[P3$CALVING_AGE > 2278 & ! is.na(P3$CALVING_AGE)] <- NA
P3$CALVING_AGE[P3$CALVING_AGE < 1096 & ! is.na(P3$CALVING_AGE)] <- NA
summary(P3$CALVING_AGE)

P4$CALVING_AGE[P4$CALVING_AGE > 2771 & ! is.na(P4$CALVING_AGE)] <- NA
P4$CALVING_AGE[P4$CALVING_AGE < 1437 & ! is.na(P4$CALVING_AGE)] <- NA
summary(P4$CALVING_AGE)

P5$CALVING_AGE[P5$CALVING_AGE > 3248 & ! is.na(P5$CALVING_AGE)] <- NA
P5$CALVING_AGE[P5$CALVING_AGE < 1772 & ! is.na(P5$CALVING_AGE)] <- NA
summary(P5$CALVING_AGE)


combi.per.herd <- tapply(HF4$ANIMAL_NUMBER, HF4$HYS, function(x){length(unique(x))})
 
herd.names <- levels(HF4$HYS)

for(i in 1:length(herd.names)){
 
 if(combi.per.herd[i] < 3){
 
    which <- (herd.names[i] == HF4$HYS)

    HF4 <- HF4[! which,]
 }
}

animals.per.herd <- tapply(dat$ANIMAL_NUMBER, dat$PARTICIPANT, function(x){length(unique(x))})
 
herd.names <- levels(dat$PARTICIPANT)

for(i in 1:length(herd.names)){
 
 if(animals.per.herd[i] <= 3){
 
    which <- (herd.names[i] == dat$PARTICIPANT)

    dat <- dat[! which,]
 }
}


hst$CALVING_DTM <-as.Date(hst$CALVING_DTM,"%d/%m/%Y")
hst$BIRTH_DATE<-as.Date(hst$BIRTH_DATE,"%d/%m/%Y")

ordering animals by parity to get differences between calving dates ####
hos<-hst[order(hst$ANIMAL_NUMBER,hst$PARITY),]
head(hos)
###################################
**************animals across parity 1 and 5.
this was done to determine which animals appeared  from 1  to 5 so as to check intervals between calving.
pa<-tapply(p4$PARITY, p4$ANIMAL_NUMBER, function(x) length(table(x))) ##hos is HF_SA HYSlactationfile
table(pa)# list of ani_id with 1, 2, 3, or more lact per animal
table(pa==5)
pa5<-names(pa[pa==5]) ##to pull out names of  cows from pa to pa5; NAMES of all ani_id in lact1-5
p5<-hos[hos$ANIMAL_NUMBER%in%pa5, ]
p5<-hos[hos$ANIMAL_NUMBER%in%as.numeric(pa5), ]; ##because the ani_ids are numeric; maybe i should change to factor or do as.numeric. 
					p5 has the RECORD LIST of all ani_id in all parities
length(pa5)*5  # multiplied by 5 since there are 5 lactations per animal

CALVING INTERVAL FOR 2,521 HOLSTEIN-FRIESIAN COWS IN SOUTH AFRICA
k<-1
for(i in ANI){
  extract<-subset(ci,ci$ANI_ID==i)
  date1<-extract$CALVING_DATE[2]-extract$CALVING_DATE[1]
  date2<-extract$CALVING_DATE[3]-extract$CALVING_DATE[2]
  date3<-extract$CALVING_DATE[4]-extract$CALVING_DATE[3]
  date4<-extract$CALVING_DATE[5]-extract$CALVING_DATE[4]
  write(paste(ANI[k]," ",as.numeric(date1),as.numeric(date2),as.numeric(date3),as.numeric(date4)),file="ci1",append=T)
  k<-k+1
 ###############
 #####################NUMBER OF DAUGHTERS PER SIRES
HST <- HST[order(HST$ANIMAL_NUMBER,HST$PARITY),]
HST$SIRE_ID<-as.factor(HST$SIRE_ID)
listlevels<-levels(HST$SIRE_ID)
listlevels[1:5]
#HST_DAUGHTERS<-matrix(,nrow=0, ncol=length(HST[1,]),colnames(HST))                             #Convert SIRES as a factor
singlecounts<-subset(HST, !duplicated(ANIMAL_NUMBER))          #Removes duplicates ANIMAL_IDs (creates a file containing one row per ANIMAL_ID)
count<-as.data.frame(table(singlecounts$SIRE_ID))                      #Counts in how many rows every SIRE appears (thus, equal to the number of daugthers each sire has)
HST_extract5more<-subset(count,count$Freq>=5)                                  #Extracts a dataset with only those sires with 5 or more daughters
HST_DAUGHTERS<-matrix(, nrow = 0, ncol = length(HST[1,]))             #Create void array to add results
colnames(HST_DAUGHTERS)=colnames(HST)                                             #Give the void array the same column names as your original array
for (i in HST_extract5more$Var1){                                                                       #Loop for all sires with 5 or more daughters
  extractedinformation<-subset(HST,HST$SIRE_ID==i)  #extract subset of daugther information corresponding to the sire i
  HST_DAUGHTERS<-rbind(HST_DAUGHTERS,extractedinformation)                            #Append the extracted subset to the array containing all the information
}

#####################NUMBER OF DAUGHTERS PER SIRES
siredota <- siredota[order(siredota$ANI_ID,siredota$PARITY),]
siredota$SIRE_ID<-as.factor(siredota$SIRE_ID)
listlevels<-levels(siredota$SIRE_ID)
listlevels[1:5]
#finaldata<-matrix(,nrow=0, ncol=length(siredota[1,]),colnames(siredota))                             #Convert SIRES as a factor
singlecounts<-subset(siredota, !duplicated(ANI_ID))          #Removes duplicates ANIMAL_IDs (creates a file containing one row per ANIMAL_ID)
count<-as.data.frame(table(singlecounts$SIRE_ID))                      #Counts in how many rows every SIRE appears (thus, equal to the number of daugthers each sire has)
extract5more<-subset(count,count$Freq>=5)                                  #Extracts a dataset with only those sires with 5 or more daughters
finaldata<-matrix(, nrow = 0, ncol = length(siredota[1,]))             #Create void array to add results
colnames(finaldata)=colnames(siredota)                                             #Give the void array the same column names as your original array
for (i in extract5more$Var1){                                                                       #Loop for all sires with 5 or more daughters
  extractedinformation<-subset(siredota,siredota$SIRE_ID==i)  #extract subset of daugther information corresponding to the sire i
  finaldata<-rbind(finaldata,extractedinformation)                            #Append the extracted subset to the array containing all the information
}

##############################################################
k<-1
HST_DAUGHTERS$ANIMAL_NUMBER<-as.factor(HST_DAUGHTERS$ANIMAL_NUMBER)
HST_DAUGHTERS$CI <- NA
ANI<-levels(HST_DAUGHTERS$ANIMAL_NUMBER)
for(i in ANI){
  #print.noquote(i)
  ikeep <- (HST_DAUGHTERS$ANIMAL_NUMBER==i)
  ikeep1 <- ikeep & (HST_DAUGHTERS$PARITY == 1)
  ikeep2 <- ikeep & (HST_DAUGHTERS$PARITY == 2)
  ikeep3 <- ikeep & (HST_DAUGHTERS$PARITY == 3)
  ikeep4 <- ikeep & (HST_DAUGHTERS$PARITY == 4)
  #ikeep5 <- ikeep & (HST_DAUGHTERS$PARITY == 5)
  if(sum(ikeep1) > 0){
    if(sum(ikeep2) > 0){
      HST_DAUGHTERS$CI[ikeep1] <- as.numeric(HST_DAUGHTERS$CALVING_DTM[ikeep2] - HST_DAUGHTERS$CALVING_DTM[ikeep1])
    }
  }
  if(sum(ikeep2) > 0){
    if(sum(ikeep3) > 0){
      HST_DAUGHTERS$CI[ikeep2] <- as.numeric(HST_DAUGHTERS$CALVING_DTM[ikeep3] - HST_DAUGHTERS$CALVING_DTM[ikeep2])
    }
  }
  if(sum(ikeep3) > 0){
    if(sum(ikeep4) > 0){
      HST_DAUGHTERS$CI[ikeep3] <- as.numeric(HST_DAUGHTERS$CALVING_DTM[ikeep4] - HST_DAUGHTERS$CALVING_DTM[ikeep3])
    }
  }
  if(sum(ikeep4) > 0){
    if(sum(ikeep5) > 0){
      HST_DAUGHTERS$CI[ikeep4] <- as.numeric(HST_DAUGHTERS$CALVING_DTM[ikeep5] - HST_DAUGHTERS$CALVING_DTM[ikeep4])
    }
  }
}
##########################################
dat3$CowID<-as.factor(dat3$CowID)

dat3$CI <- NA
dat3$CI1 <- NA
dat3$CI2 <- NA
dat3$CI3 <- NA
dat3$CI4 <- NA
ANI<-levels(JJ$CowID)

for(i in ANI){

  #print.noquote(i)

  ikeep <- (JJ$CowID==i)

  ikeep1 <- ikeep & (JJ$LactNo == 1)
  ikeep2 <- ikeep & (JJ$LactNo == 2)
  ikeep3 <- ikeep & (JJ$LactNo == 3)
  ikeep4 <- ikeep & (JJ$LactNo == 4)
  ikeep5 <- ikeep & (JJ$LactNo == 5)

  if(sum(ikeep1) > 0){

    if(sum(ikeep2) > 0){

      JJ$CI[ikeep1] <- as.numeric(JJ$LactStartDate[ikeep2] - JJ$LactStartDate[ikeep1])
    }
  }

  if(sum(ikeep2) > 0){

    if(sum(ikeep3) > 0){

      JJ$CI[ikeep2] <- as.numeric(JJ$LactStartDate[ikeep3] - JJ$LactStartDate[ikeep2])
    }
  }

  if(sum(ikeep3) > 0){

    if(sum(ikeep4) > 0){

     JJ$CI[ikeep3] <- as.numeric(JJ$LactStartDate[ikeep4] - JJ$LactStartDate[ikeep3])
    }
  }

  if(sum(ikeep4) > 0){

    if(sum(ikeep5) > 0){

      JJ$CI[ikeep4] <- as.numeric(JJ$LactStartDate[ikeep5] - JJ$LactStartDate[ikeep4])
    }
  }
}
