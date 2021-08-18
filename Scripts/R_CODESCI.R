get.lmy <- function(Lact1, refd){

 cows <- unique(Lact1$ANIMAL_NUMBER)

ncows <- length(cows)

Lact1$lmy <- NA
  Lact1$ivals <- NA
 for(i in 1:ncows){
 
  

    test <- (Lact1$ANIMAL_NUMBER == cows[i])
 
    tmp <-  Lact1[test,]
   
   alldates <- c(tmp$CALVING_DATE[1], tmp$TD)

   ivals <- as.numeric(as.character(diff(alldates)))
 ivals <- as.numeric(as.character(diff(c(tmp$CALVING_DATE[1], tmp$TD))))
  n <- length(ivals)
 
   ivals <- c(ivals, refd - sum(ivals))

    mvals <- c(tmp$MILK_YIELD[1], (tmp$MILK_YIELD[-n] + tmp$MILK_YIELD[-1])/2, tmp$MILK_YIELD[n])
 
   Lact1$lmy[test] <- sum(ivals * mvals)

   Lact1$ivals[test] <- cumsum(ivals[-length(ivals)])
  }
  
  Lact1
 }
 Lact1<- get.lmy(Lact1,refd=305)

####################################################
k<-1

DAT1$ANIMAL_NUMBER<-as.factor(DAT1$ANIMAL_NUMBER)

DAT1$CI <- NA

ANI<-levels(DAT1$ANIMAL_NUMBER)

for(i in ANI){

  #print.noquote(i)

  ikeep <- (DAT1$ANIMAL_NUMBER==i)

  ikeep1 <- ikeep & (DAT1$PARITY == 1)
  ikeep2 <- ikeep & (DAT1$PARITY == 2)
  ikeep3 <- ikeep & (DAT1$PARITY == 3)
  ikeep4 <- ikeep & (DAT1$PARITY == 4)
  ikeep5 <- ikeep & (DAT1$PARITY == 5)

  if(sum(ikeep1) > 0){

    if(sum(ikeep2) > 0){

      DAT1$CI[ikeep1] <- as.numeric(DAT1$CALVING_DTM[ikeep2] - DAT1$CALVING_DTM[ikeep1])
    }
  }

  if(sum(ikeep2) > 0){

    if(sum(ikeep3) > 0){

      DAT1$CI[ikeep2] <- as.numeric(DAT1$CALVING_DTM[ikeep3] - DAT1$CALVING_DTM[ikeep2])
    }
  }

  if(sum(ikeep3) > 0){

    if(sum(ikeep4) > 0){

      DAT1$CI[ikeep3] <- as.numeric(DAT1$CALVING_DTM[ikeep4] - DAT1$CALVING_DTM[ikeep3])
    }
  }

  if(sum(ikeep4) > 0){

    if(sum(ikeep5) > 0){

      DAT1$CI[ikeep4] <- as.numeric(DAT1$CALVING_DTM[ikeep5] - DAT1$CALVING_DTM[ikeep4])
    }
  }
}

#########################################################
to get EBVs from .sln file than manually deleting in R
SLN<-read.table("D:\\Users\\OAbejide\\JERSEY_SA290915\\JERSEY_HYS_MY.sln", head=FALSE)
> head(SLN)
SLN2<-SLN[which(SLN$V1=="ANIMAL_NUMBER"),]
> nrow(SLN2)
[1] 509492
> head(SLN2)
########################
EVERY COW MUST HAVE A FIRST LACTATION R-SCRIPT
cow file= 6td_3hys_5d
khf=21,111 records
data4 <-merge(khf,cow, by.x=c("CowID"),by.y=c("CowID"))
tmp<-cow[which(cow$CowID %in% khf$CowID),]

tmp2<-cow[which(khf$CowID %in% cow$CowID),]

C22<-cow[cow$CowID%in%khf$CowID]

sol<-data.frame()
for (i in 1:length(afc$CowID)){
  a<-subset(KHF3,KHF3$CowID==afc$CowID[i])
  sol<-rbind(sol,a)
}

sol$C_HYS<-as.factor(sol$C_HYS)                                    #Convert as factor
niv<-levels(sol$C_HYS)                                           #List of levels of the factor
extracted<-as.data.frame(which(table(sol$C_HYS)>=3))[,1]         #List of levels that appear three or more times
rep_HYS<-vector(length=length(extracted))                      #Create vector to store HYS levels that appear three or more times
for (i in 1:length(extracted)){                                 #Loop to go through the entire list of repeated levels                       
  rep_HYS[i]<-niv[extracted[i]]                                #The levels we want is the one which is: niv[listed level in extracted]
}

#rep_HYS                                                        #List of herds that appear three or more times
sol_3HYS_3TD<-merge(sol,rep_HYS,by="C_HYS",by.y=TRUE)	

##############################################

     as.numeric(extract$CALVING_DATE[extract$PARITY == 2] - extract$CALVING_DATE[extract$PARITY == 1])

  diffs <- c(as.numeric(extract$CALVING_DATE[extract$PARITY == 2] - extract$CALVING_DATE[extract$PARITY == 1]),
             as.numeric(extract$CALVING_DATE[extract$PARITY == 3] - extract$CALVING_DATE[extract$PARITY == 2]),
             as.numeric(extract$CALVING_DATE[extract$PARITY == 4] - extract$CALVING_DATE[extract$PARITY == 3]),
             as.numeric(extract$CALVING_DATE[extract$PARITY == 5] - extract$CALVING_DATE[extract$PARITY == 4]))

  


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