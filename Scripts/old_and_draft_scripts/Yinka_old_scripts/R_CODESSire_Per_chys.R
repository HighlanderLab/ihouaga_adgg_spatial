A = matrix( 
+   c(2, 4, 3, 1, 5, 7), # the data elements 
+   nrow=c,              # number of rows 
+   ncol=3,              # number of columns d
+   byrow = TRUE)

sire_count<-mp[duplicated(mp[,c("SIRE_ID", "ANIMAL_NUMBER")])==FALSE,]  # check sires and daughters common to SA and Kenya
> nrow(sire_count)
[1] 2468
> dota_persire<-as.data.frame(table(sire_count$SIRE_ID))  # sire and frequency of daughters to sires
> head(dota_persire)
> dota_persire<-as.data.frame(table(sire_count$SIRE_ID))
>  head(dota_persire)
##############################
dat_herd_freq<-as.data.frame(table(data_herd$HerdID,data_herd$Sire))
#####################
sire_locationt<-data2[duplicated(data2[,c("Sire_recode", "Address")])==FALSE,]
finaldata<-data.frame()
for (i in 1:72){                                    #Loop for all herds 
  ADD1<-subset(data2,data2$Address.1==i)               #Get only herd 1 from main data
  HERD1<-cbind(ADD1$Address.1,ADD1$Sire_recode)             #Only needed columns for sire and HYS
  HERD1<-as.data.frame(HERD1)                      #Convert to dataframe
  HERD1$V3<-paste(HERD1$V1,HERD1$V2,sep=" ")       #Combine columns V1 and V2 into V3
  FINALDATSETHERD1<-HERD1[!duplicated(HERD1$V3),]  #Remove duplicates from V3 (combinations of sire and HYR)
  finaldata<-rbind(finaldata,FINALDATSETHERD1)
}
table(finaldata$V2)                       #Table for 
################################################
**pa<-tapply(KHF$Sire, KHF$HerdID, function(x) length(table(x))) #number of sires per herd
#########
KHF<- KHF[order(KHF$Sire, KHF$C_HYS),]

#################### number of CHYS counts per sire to detect confounding efects of HYS
KHF_NEW_3CHYSm2_Nrec$Sire <- factor(KHF_NEW_3CHYSm2_Nrec$Sire)
nrecs.per.sire <- tapply(KHF_NEW_3CHYSm2_Nrec$C_HYS, KHF_NEW_3CHYSm2_Nrec$Sire, function(x){length(unique(x))})
KHF_NEW_3CHYSm2_Nrec$CHYS_Nrec <- nrecs.per.sire[match(KHF_NEW_3CHYSm2_Nrec$Sire, levels(KHF_NEW_3CHYSm2_Nrec$Sire))]
KHF_NEW_3CHYSm2_Nrecs <- KHF_NEW_3CHYSm2_Nrec[KHF_NEW_3CHYSm2_Nrec$CHYS_Nrec >= 2,]

####################
cow file= 6td_3hys_5d
khf=21,111 records
data4 <-merge(khf,cow, by.x=c("CowID"),by.y=c("CowID"))
tmp<-cow[which(cow$CowID %in% khf$CowID),]

tmp2<-cow[which(khf$CowID %in% cow$CowID),]

C22<-cow[cow$CowID%in%khf$CowID]
####################################################### #all cows must have a 1st lact
#
data.frame()
 for (i in 1:length(afc$CowID_recode)){
   a<-subset(R4,R4$CowID_recode==afc$CowID_recode[i])
 if(length(Cow_Lact15) == 0) {Cow_Lact15 = a}
 else {Cow_Lact15<-rbind(Cow_Lact15,a)}
 }
###
KHF_NEW_3CHYSm2_Nrecs6TD_1st<-data.frame()
 for (i in 1:length(afc$CowID)){
  a<-subset(KHF_NEW_3CHYSm2_Nrecs6TD,KHF_NEW_3CHYSm2_Nrecs6TD$CowID==afc$CowID[i])
  KHF_NEW_3CHYSm2_Nrecs6TD_1st<-rbind(KHF_NEW_3CHYSm2_Nrecs6TD_1st,a)
}
 nrow(KHF_NEW_3CHYSm2_Nrecs6TD_1st)
[1] 13121
> table(sol$LactNo)
########################## #HYS combinations
KHF_NEW_3CHYSm3_Nrecs6TD_1st5D$C_HYS<-as.factor(KHF_NEW_3CHYSm3_Nrecs6TD_1st5D$C_HYS)                                    #Convert as factor
niv<-levels(KHF_NEW_3CHYSm3_Nrecs6TD_1st5D$C_HYS)                                           #List of levels of the factor
extracted<-as.data.frame(which(table(KHF_NEW_3CHYSm3_Nrecs6TD_1st5D$C_HYS)>=5))[,1]         #List of levels that appear three or more times
rep_HYS<-vector(length=length(extracted))                      #Create vector to store HYS levels that appear three or more times
for (i in 1:length(extracted)){                                 #Loop to go through the entire list of repeated levels                       
  rep_HYS[i]<-niv[extracted[i]]                                #The levels we want is the one which is: niv[listed level in extracted]
}

#rep_HYS                                                        #List of herds that appear three or more times
KHF_NEW_5CHYSm3_Nrecs6TD_1st5D<-merge(KHF_NEW_3CHYSm3_Nrecs6TD_1st5D,rep_HYS,by="C_HYS",by.y=TRUE)	

##############
p4<- p4[order(p4$CowID_recode, p$TestDayDate),]
####################
p4$CowID_recode <- factor(p4$CowID_recode)

nrecs.per.cow <- tapply(p4$TestDayDate, p4$CowID_recode, function(x){length(unique(x))})

p4$Nrecs <- nrecs.per.cow[match(p4$CowID_recode, levels(p4$CowID_recode))]

p4_new <- p4[p4$Nrecs >= 1,]

p4_new$Interval <- NA

for(i in 1:length(levels(p4$CowID_recode))){

  zz <- (p4_new$CowID_recode == levels(p4$CowID_recode)[i])

  if(sum(zz) > 0){

    tmp <- p4_new[zz,]

    tmp$Interval <- c(NA, as.numeric(diff(tmp$TestDayDate)))

    p4_new$Interval[zz] <- tmp$Interval
  }
}


######################
get.lmy <- function(p5_new, refd){

 cows <- unique(p5_new$CowID_recode)

 ncows <- length(cows)

 ## ## ##

 p5_new$lmy <- NA
 
 p5_new$dim <- NA

 for(i in 1:ncows){

   ## print.noquote(paste(i, date()))

   test <- (p5_new$CowID_recode == cows[i])

   tmp <-  p5_new[which(p5_new$CowID_recode == cows[i]),]
   
   alldates <- c(tmp$LactStartDate[1], tmp$TestDayDate)

   ivals <- as.numeric(as.character(diff(alldates)))

   ## ################
   ## New code added 28 May 2015 to remove obs whose time from lactation is greater than 'refd'
   
   keepvals <- (cumsum(ivals) < refd)

   tmp <- tmp[keepvals,]

   alldates <- c(tmp$LactStartDate[1], tmp$TestDayDate)

   ivals <- as.numeric(as.character(diff(alldates)))
   
   ## #################
   
   n <- length(ivals)

   p5_new$dim[test] <- sum(ivals)
   
   ivals <- c(ivals, refd - sum(ivals))

   mvals <- c(tmp$Yield[1], (tmp$Yield[-n] + tmp$Yield[-1])/2, tmp$Yield[n])
   
   p5_new$lmy[test] <- sum(ivals * mvals)
 }
 
 p5_new
}
p5_new <- get.lmy(p5_new, refd=305)

#######################
KenyaCow_Lact15_sire5$CowID_recode<-as.factor(KenyaCow_Lact15_sire5$CowID_recode)
k<-1
#KenyaCow_Lact15_sire5$CI <- NA
KenyaCow_Lact15_sire5$CI1 <- NA
KenyaCow_Lact15_sire5$CI2 <- NA
KenyaCow_Lact15_sire5$CI3 <- NA
KenyaCow_Lact15_sire5$CI4 <- NA
ANI<-levels(KenyaCow_Lact15_sire5$CowID_recode)
for(i in ANI){

  #print.noquote(i)

  ikeep <- (KenyaCow_Lact15_sire5$CowID_recode==i)

  ikeep1 <- ikeep & (KenyaCow_Lact15_sire5$LactNo == 1)
  ikeep2 <- ikeep & (KenyaCow_Lact15_sire5$LactNo == 2)
  ikeep3 <- ikeep & (KenyaCow_Lact15_sire5$LactNo == 3)
  ikeep4 <- ikeep & (KenyaCow_Lact15_sire5$LactNo == 4)
  ikeep5 <- ikeep & (KenyaCow_Lact15_sire5$LactNo == 5)

  if(sum(ikeep1) > 0){

    if(sum(ikeep2) > 0){

      KenyaCow_Lact15_sire5$CI1[ikeep1] <- as.numeric(KenyaCow_Lact15_sire5$LactStartDate[ikeep2] - KenyaCow_Lact15_sire5$LactStartDate[ikeep1])
    }
  }

  if(sum(ikeep2) > 0){

    if(sum(ikeep3) > 0){

      KenyaCow_Lact15_sire5$CI2[ikeep2] <- as.numeric(KenyaCow_Lact15_sire5$LactStartDate[ikeep3] - KenyaCow_Lact15_sire5$LactStartDate[ikeep2])
    }
  }

  if(sum(ikeep3) > 0){

    if(sum(ikeep4) > 0){

      KenyaCow_Lact15_sire5$CI3[ikeep3] <- as.numeric(KenyaCow_Lact15_sire5$LactStartDate[ikeep4] - KenyaCow_Lact15_sire5$LactStartDate[ikeep3])
    }
  }

  if(sum(ikeep4) > 0){

    if(sum(ikeep5) > 0){

      KenyaCow_Lact15_sire5$CI4[ikeep4] <- as.numeric(KenyaCow_Lact15_sire5$LactStartDate[ikeep5] - KenyaCow_Lact15_sire5$LactStartDate[ikeep4])
    }
  }
}

###########
sire_count<-KenyaCow_Lact15[duplicated(KenyaCow_Lact15[,c("Sire_recode", "CowID_recode")])==FALSE,]
nrow(sire_count)
dota_persire<-as.data.frame(table(sire_count$Sire_recode))  # sire and frequency of daughters to sires
head(dota_persire)
dota_persire<-as.data.frame(table(sire_count$Sire_recode))
head(dota_persire)
###################
KenyaCow_Lact15_SIRE5$C_HYS<-as.factor(KenyaCow_Lact15_SIRE5$C_HYS)                                    #Convert as factor
niv<-levels(KenyaCow_Lact15_SIRE5$C_HYS)                                           #List of levels of the factor
extracted<-as.data.frame(which(table(KenyaCow_Lact15_SIRE5$C_HYS)>=3))[,1]         #List of levels that appear three or more times
rep_HYS<-vector(length=length(extracted))                      #Create vector to store HYS levels that appear three or more times
for (i in 1:length(extracted)){                                 #Loop to go through the entire list of repeated levels                       
  rep_HYS[i]<-niv[extracted[i]]                                #The levels we want is the one which is: niv[listed level in extracted]
}

#rep_HYS                                                        #List of herds that appear three or more times

KenyaCow_Lact15_SIRE5_3HYS<-merge(KenyaCow_Lact15_SIRE5,rep_HYS,by="C_HYS",by.y=TRUE)

#########################
hys_count<-KenyaCow_lmy_2TD_3HYS[duplicated(KenyaCow_lmy_2TD_3HYS[,c("C_HYS", "CowID_recode")])==FALSE,]
nrow(hys_count)
dota_hys<-as.data.frame(table(hys_count$C_HYS))  # sire and frequency of daughters to sires
head(dota_hys)
dota_hys<-as.data.frame(table(hys_count$C_HYS))
head(dota_hys)
##############
KHFsire_count<-KHF[duplicated(KHF[,c("HerdID", "Sire")])==FALSE,]
finaldata<-data.frame()
for (i in 1:203){                                    #Loop for all herds 
  herd1<-subset(KHF,KHF$HerdID==i)               #Get only herd 1 from main data
  HERD1<-cbind(herd1$Sire,herd1$C_HYS)             #Only needed columns for sire and HYS
  HERD1<-as.data.frame(HERD1)                      #Convert to dataframe
  HERD1$V3<-paste(HERD1$V1,HERD1$V2,sep=" ")       #Combine columns V1 and V2 into V3
  FINALDATSETHERD1<-HERD1[!duplicated(HERD1$V3),]  #Remove duplicates from V3 (combinations of sire and HYR)
  finaldata<-rbind(finaldata,FINALDATSETHERD1)
}
table(finaldata$V2)                       #Table for 
###############################
nrow(KHFsire_count)
sires_perherd<-as.data.frame(table(KHFsire_count$HerdID))  # sire and frequency of HERDS to sires
head(sires_perherd)
sires_perherd<-as.data.frame(table(KHFsire_count$HerdID))
head(sires_perherd)
#########
out <- tapply(new$C_HYS, list(new$Sire, new$HerdID), function(x){length(unique(x))})
###########
khf_Sherd<-data2[duplicated(data2[,c("Address", "Sire_recode")])==FALSE,]
nrow(khf_Sherd)
sire_per<-as.data.frame(table(data2$Address))
head(sire_per)
#number of herds with more than/equal to 2 jersey cows
nrow(anim_per_Participant[which(anim_per_Participant$Freq>=2),])
######################
khf_herd10<-khf[duplicated(khf[,c("ANIMAL_NUMBER", "PARTICIPANT")])==FALSE,]
nrow(khf_herd2)
anim_per_Participant<-as.data.frame(table(khf_herd2$PARTICIPANT))
head(anim_per_Participant)
#number of herds with more than/equal to 2 jersey cows
nrow(anim_per_Participant[which(anim_per_Participant$Freq>=2),])
#####
AFC$AGE[AFC$AGE > 1783.11 & ! is.na(AFC$AGE)] <- NA
AFC$AGE[AFC$AGE < 195.69 & ! is.na(AFC$AGE)] <- NA
AFC_AGE<-na.omit(AFC)
#################
DETERMINE SIRES THAT OCCUR IN MORE THAN ONE HERD
herd_vs_sire<-table(Sire=data2$Sire_recode,herd=data2$HerdID)
herd_vs_sire_logic<-herd_vs_sire>0 #sires appear at least once in a herd
nherds_per_sire<-rowSums(herd_vs_sire_logic)
nherds_per_sire[nherds_per_sire>1]
names(nherds_per_sire)[nherds_per_sire>1]
sires_nonunique<-names(nherds_per_sire)[nherds_per_sire>1]
data3<-data2[data2$Sire_recode%in%sires_nonunique,]
data3
head(data3)

