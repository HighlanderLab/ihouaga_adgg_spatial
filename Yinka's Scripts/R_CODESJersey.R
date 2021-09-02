MILK_305D<-jer1$MILK_YLD_305D 
PROTEIN_305D<-jer1$PROTEIN_YLD_305D
FAT_305D<-jer1$FAT_YLD_305D                             #Save a vector with ages
IDs<-jer1$ANIMAL_NUMBER 
Parity<-jer1$PARITY                   
MILK_305D<-vector(length=length(IDs))		#Create  a vector for all animals  that will store the mean of MILK_YLD_305D
PROTEIN_305D<-vector(length=length(IDs))	#Create  a vector for all animals  that will store the mean of PROTEIN_YLD_305D
FAT_305D<-vector(length=length(IDs))	#Create  a vector for all animals  that will store the mean of FAT_YLD_305D
SD_M<-vector(length=length(IDs))
SD_P<-vector(length=length(IDs))
SD_F<-vector(length=length(IDs))
Parity<-vector(length=length(IDs))
for(i in 1:length(IDs)){                            #Loop for all animals in ID vector
     
MILK_305D[i]<-mean(jer1$MILK_YLD_305D[which(jer1$ANIMAL_NUMBER==IDs[i])]) 
PROTEIN_305D[i]<-mean(jer1$PROTEIN_YLD_305D[which(jer1$ANIMAL_NUMBER==IDs[i])])
FAT_305D[i]<-mean(jer1$FAT_YLD_305D[which(jer1$ANIMAL_NUMBER==IDs[i])])       #Compute total as the sum of NEWYIELD and removing NAs
SD_M[i]<-sd(jer1$MILK_YLD_305D[which(jer1$ANIMAL_NUMBER==IDs[i])])           #Compute total as the sum of NEWYIELD and removing NAs
SD_P[i]<-sd(jer1$PROTEIN_YLD_305D[which(jer1$ANIMAL_NUMBER==IDs[i])])
SD_F[i]<-sd(jer1$FAT_YLD_305D[which(jer1$ANIMAL_NUMBER==IDs[i])])
DIM[i]<-mean(jer1$LACTATION_LENGTH[which(jer1$ANIMAL_NUMBER==IDs[i])])
}

Jersey_PARITY5_305<-array(0,dim=c(length(IDs),7))
Jersey_PARITY5_305[,1]<-IDs
Jersey_PARITY5_305[,2]<-MILK_305D
Jersey_PARITY5_305[,3]<-SD_M
Jersey_PARITY5_305[,4]<-PROTEIN_305D
Jersey_PARITY5_305[,5]<-SD_P
Jersey_PARITY5_305[,6]<-FAT_305D
Jersey_PARITY5_305[,7]<-SD_F
head(Jersey_PARITY5_305)

That creates the new column named "MY_NEW_COLUMN" filled with "NA"
JERSEY_LACT1["DIM"] <- NA 
JERSEY_LACT1$DIM <- LACTATION_DTM - CALVING_DTM 

Jersey_PARITY1$LACTATION_DTM<-as.Date(as.character(Jersey_PARITY1$LACTATION_DTM),"%d/%m/%Y")
data$CALVING_DATE<-as.Date(as.character(data$CALVING_DTM),"%d/%m/%Y")

##################################################
meansDIM<-aggregate(MILK_YLD_305~LACT_LENGTH, data=JERSEY_LACT1, FUN=mean)
plot(meansDIM$LACT_LENGTH, meansDIM$MILK_YLD_305, type="l")

##################################################
CONVERT MONTH IN DATES FROM LETTERS  TO NUMBER FORMAT
p4["LACT_DATE"]<-NA
p4$LACT_DATE<-as.Date(p4$LACTATION_DTM,"%d-%b-%y")
p4$LACT_DATE<-format(p4$LACT_DATE,"%d/%m/%Y")
p4["CALVING_DATE"]<-NA
p4$CALVING_DATE<-as.Date(p4$CALVING_DTM,"%d-%b-%y")
p4$LCALVING_DATE<-format(p4$CALVING_DATE,"%d/%m/%Y")

CONVERTING NUMERIC DATES TO DATES ("%d/%m/%Y")
P4$test<-format(as.Date(as.character(P4$BIRTH_DATE), format="%Y%m%d"),"%d/%m/%Y") ####use as.character 1st, then as.Date 2nd; do not join together
jersey1$DOB<-as.Date(jersey1$DOB,"%d/%m/%Y")
p4$DOB<-as.Date(p4$BIRTH_DATE, format= "%Y%m%d")
p4$DOB<-as.numeric(p4$BIRTH_DATE, format= "%Y%m%d")
prodped$SEASON <- factor(as.character(prodped$LACT_MTH), levels=1:12)  #### SET LEVELS FOR COMPUTING SEASONS IN SOUTH AFRICA
jersey1$CALVING_DATE<-as.Date(jersey1$CALVING_DATE,"%d/%m/%Y") # treat calving dates as.Date in each parity
jersey2<-jersey2[order(as.Date(jersey2$LACT_DATE, format="%d/%m/%Y")),]
jersey3<-jersey3[order(as.Date(jersey3$LACT_DATE, format="%d/%m/%Y")),]
jersey4<-jersey4[order(as.Date(jersey4$LACT_DATE, format="%d/%m/%Y")),]
jersey5<-jersey5[order(as.Date(jersey5$LACT_DATE, format="%d/%m/%Y")),]
jerseyy$HYS<-as.factor(jerseyy$HYS)                                    #Convert as factor
niv<-levels(jerseyy$HYS)                                           #List of levels of the factor
extracted<-as.data.frame(which(table(jerseyy$HYS)>=3))[,1]         #List of levels that appear three or more times
rep_HYS<-vector(length=length(extracted))                      #Create vector to store HYS levels that appear three or more times
for (i in 1:length(extracted)){                                 #Loop to go through the entire list of repeated levels                       
  rep_HYS[i]<-niv[extracted[i]]                                #The levels we want is the one which is: niv[listed level in extracted]
}
head(rep_HYS)
prodped$test<-format(as.character(prodped$BIRTH_DATE, format="%Y%m%d"),"%d/%m/%Y")
SETTING CALVING_AGE AT MEAN3SD
p1$CALVING_AGE[p1$CALVING_AGE > 1225 & ! is.na(p1$CALVING_AGE)] <- NA
p1$CALVING_AGE[p1$CALVING_AGE < 365 & ! is.na(p1$CALVING_AGE)] <- NA
summary(p1$CALVING_AGE)

p2$CALVING_AGE[p2$CALVING_AGE > 1775 & ! is.na(p2$CALVING_AGE)] <- NA
p2$CALVING_AGE[p2$CALVING_AGE < 640 & ! is.na(p2$CALVING_AGE)] <- NA
summary(P2$CALVING_AGE)

p3$CALVING_AGE[p3$CALVING_AGE > 2251 & ! is.na(p3$CALVING_AGE)] <- NA
p3$CALVING_AGE[p3$CALVING_AGE < 967 & ! is.na(p3$CALVING_AGE)] <- NA
summary(p3$CALVING_AGE)

p4$CALVING_AGE[p4$CALVING_AGE > 2738 & ! is.na(p4$CALVING_AGE)] <- NA
p4$CALVING_AGE[p4$CALVING_AGE < 1282 & ! is.na(p4$CALVING_AGE)] <- NA
summary(P4$CALVING_AGE)

p5$CALVING_AGE[p5$CALVING_AGE > 3191 & ! is.na(p5$CALVING_AGE)] <- NA
p5$CALVING_AGE[p5$CALVING_AGE < 1619 & ! is.na(p5$CALVING_AGE)] <- NA
summary(p5$CALVING_AGE)

NUMBER OF ANIMALS ACROSS PARITY animals across parity 1 and 5.
this was done to determine which animals appeared  from 1  to 5 so as to check intervals between calving.
pa<-tapply(p4$PARITY, p4$ANI_ID, function(x) length(table(x))) ##p4 is mainunediitedJersey_lactationfile
table(pa)# list of ani_id with 1, 2, 3, or more lact per animal
table(pa==5)
pa5<-names(pa[pa==5]) ##to pull out names of  cows from pa to pa5; NAMES of all ani_id in lact1-5
p5<-p4[p4$ANI_ID%in%pa5, ]
p5<-p4[p4$ANI_ID%in%as.numeric(pa5), ]; ##because the ani_ids are numeric; maybe i should change to factor or do as.numeric. 
					p5 has the RECORD LIST of all ani_id in all parities
length(pa5)*5  # multiplied by 5 since there are 5 lactations per animal
#######################
sire<-siredota[duplicated(siredota[,c("ANIMAL_NUMBER", "SIRE_ID")])==FALSE,]
nrow(sire)
anim_per_Participant<-as.data.frame(table(zim_herd50$HERD_ID))
head(anim_per_Participant)
number of herds with more than/equal to 50 jersey cows
nrow(anim_per_Participant[which(anim_per_Participant$Freq>=50),])

zim_herd50<-dear[duplicated(dear[,c("ANIMAL_NUMBER", "HERD_ID")])==FALSE,]
nrow(zim_herd50)
anim_per_Participant<-as.data.frame(table(zim_herd50$HERD_ID))
head(anim_per_Participant)
number of herds with more than/equal to 50 jersey cows
nrow(anim_per_Participant[which(anim_per_Participant$Freq>=50),])

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
### THIS WAS DONE TO DTEREMINE CALVING INTERVAL FROM SIRES WITH 5 OR MORE DAUGHTERS (771 jersey SIRES HAS 5>=)
k<-1
finaldata$ANI_ID<-as.factor(finaldata$ANI_ID)
ANI<-levels(finaldata$ANI_ID)
for(i in ANI){
  extract<-subset(finaldata,finaldata$ANI_ID==i)
  date1<-extract$CALVING_DATE[2]-extract$CALVING_DATE[1]
  parity1<-extract$PARITY[2]
  date2<-extract$CALVING_DATE[3]-extract$CALVING_DATE[2]
  parity2<-extract$PARITY[3]
  date3<-extract$CALVING_DATE[4]-extract$CALVING_DATE[3]
  parity3<-extract$PARITY[4]
  date4<-extract$CALVING_DATE[5]-extract$CALVING_DATE[4]
  parity4<-extract$PARITY[5]
  write(paste(ANI[k]," ",as.numeric(date1)," ",as.numeric(parity1)),file="ci4",append=T)
  write(paste(ANI[k]," ",as.numeric(date2)," ",as.numeric(parity2)),file="ci4",append=T)
  write(paste(ANI[k]," ",as.numeric(date3)," ",as.numeric(parity3)),file="ci4",append=T)
  write(paste(ANI[k]," ",as.numeric(date4)," ",as.numeric(parity4)),file="ci4",append=T)
  k<-k+1
}

D$BIRTH_MTH <- factor(as.character(D$BIRTH_MTH), levels=1:12)
levels(D$BIRTH_MTH) <- c(2,2,2,1,1,1,1,1,1,2,2,2)
season 1 (winter): Apr, May, Jun, Jul, Aug, Sep
season 2 (dry): Oct, Nov, Dec, Jan, Feb, March


y90<-subset(J1,J1$CALVING_YR=="2007")
mean(y90$MILK_YLD_305D)
sd(y90$MILK_YLD_305D)
mean(y90$AGE,na.rm=T)/28
sd(y90$AGE,na.rm=T)/28
mean(y90$CI1,na.rm=T)
sd(y90$CI1,na.rm=T)
