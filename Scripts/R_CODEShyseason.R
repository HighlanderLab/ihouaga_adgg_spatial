total4 <- read.csv("D:\\Users\\OAbejide\\305D\\total4.csv")

total4$Month <- total4$LSDseason

## ##################################

total4$Herd_ID <- factor(as.character(total4$Herd_ID))

ncows <- tapply(total4$CowID, total4$Herd_ID, function(x){length(unique(x))})

herdkeep <- names(ncows)[ncows > =50]

recordkeep <- ! is.na(match(total4$Herd_ID, herdkeep))

sum(recordkeep)
sum(! recordkeep)

total4 <- total4[recordkeep,]

dim(total4)

## ##################################

wetmonths <- c(4,5,6,11)

total4$Season <- 2 - (! is.na(match(total4$Month, wetmonths)))

write.csv(total4, file="D:\\Users\\OAbejide\\305D\\total4-modified.csv")

## ##################################

total4$Combi <- factor(paste(total4$Herd_ID, total4$LSDyear, total4$Season))
 
ncow.combi <- tapply(total4$CowID, total4$Combi, function(x){length(unique(x))})

ncow.dat <- data.frame(Combi = names(ncow.combi), Ncow = as.numeric(ncow.combi))

write.csv(ncow.dat, file="D:\\Users\\OAbejide\\305D\\ncow.combi.csv")

ncow.combi <- tapply(total4$CowID, total4$Combi, function(x){length(unique(x))})
ncow.dat <- data.frame(Combi = names(ncow.combi), Ncow = as.numeric(ncow.combi))
combikeep <- names(ncow.combi)[ncow.combi > 3]

recordkeep.combi <- ! is.na(match(new_dates$Combi, combikeep))

sum(recordkeep.combi)
sum(! recordkeep.combi)

total4 <- total4[recordkeep.combi,]

write.csv(total4, file="D:\\Users\\OAbejide\\305D\\total4-modified2.csv")

## ##################################

newest <- data.frame(CowID = tapply(total4$CowID, total4$CowID, function(x){x[1]}))

newest$Herd_ID <- tapply(total4$Herd_ID, total4$CowID, function(x){x[1]})

newest$Season <- tapply(total4$Season, total4$CowID, function(x){x[1]})

newest$LSDyear <- tapply(total4$LSDyear, total4$CowID, function(x){x[1]})

newest$Combi <- tapply(total4$Combi, total4$CowID, function(x){x[1]})

newest <- newest[order(newest$Combi),]

newest <- newest[,c(2,3,4,1)]

row.names(newest) <- NULL

newest <- newest[,c(4,1:3)]

write.csv(newest, file="D:\\Users\\OAbejide\\305D\\cow-id-combination.csv")

## ##################################

#####################NUMBER OF DAUGHTERS PER SIRES
join <- join[order(join$CowID,join$LactNo),]
join$Sire<-as.factor(join$Sire)
listlevels<-levels(join$Sire)
listlevels[1:5]
#Kenyafinaldata<-matrix(,nrow=0, ncol=length(join[1,]),colnames(join))                             #Convert SIRES as a factor
singlecounts<-subset(join, !duplicated(CowID))          #Removes duplicates ANIMAL_IDs (creates a file containing one row per ANIMAL_ID)
count<-as.data.frame(table(singlecounts$CowID))                      #Counts in how many rows every SIRE appears (thus, equal to the number of daugthers each sire has)
extract5more<-subset(count,count$Freq>=5)                                  #Extracts a dataset with only those sires with 5 or more daughters
Kenyafinaldata<-matrix(, nrow = 0, ncol = length(join[1,]))             #Create void array to add results
colnames(Kenyafinaldata)=colnames(join)                                             #Give the void array the same column names as your original array
for (i in extract5more$Var1){                                                                       #Loop for all sires with 5 or more daughters
  Kenyaextractedinfo<-subset(join,join$CowID==i)  #extract subset of daugther information corresponding to the sire i
  Kenyafinaldata<-rbind(Kenyafinaldata,Kenyaextractedinfo)                            #Append the extracted subset to the array containing all the information
}

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