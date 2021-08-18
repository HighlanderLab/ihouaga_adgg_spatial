####sln2<-sln[which(sln$V1=="ANIMAL_NUMBER"),]#### selecting a column using "which" in R

ZIM$AGE[ZIM$AGE > 1465 & ! is.na(ZIM$AGE)] <- NA
ZIM$AGE[ZIM$AGE < 1099 & ! is.na(ZIM$AGE)] <- NA
summary(ZIM$AGE)
D1$CALVING_DATE<-as.Date(D1$CALVING.DATE,"%d-%b-%y")

D1<- D1[order(D1$ANIMAL_NUMBER, D1$TD),]
zimdata <- zimdata[orde$ANIMAL_NUMBER, zimdata$CALVING_DATE),]

### checking for negative ages at calving###  ag<-zZIM[which(ZIM$AGE  =< 0),]###
### age<-ZIM[drop(ZIM$AGE > 0),] #  drop 0 values of ages at calving
mapping_310316.csv has NAs removed and to be used for TIM for milk_yield.  there are cows with 2 or more calving dates  which will be seperated into lactations
#######################
 ZIM<-ZIM_AGE
> range(ZIM$AGE,na.rm=T)
[1]  491 3818
 ZIM$TD <- as.Date(as.character(ZIM$TD), "%d/%m/%Y")
ZIM$CALVING_DATE <- as.Date(as.character(ZIM$CALVING.DATE), "%d/%m/%Y")
 ZIM<-ZIM[order(ZIM$ANIMAL_NUMBER,ZIM$CALVING_DATE,ZIM$TD),]
ZIM$AGE[ZIM$AGE > 1098 & ! is.na(ZIM$AGE)] <- NA
>range(ZIM$AGE,na.rm=T)
[1]  491 1098
> L1<-ZIM
> range(L1$AGE,na.rm=T)
[1]  491 1098
> ZIM<-ZIM_AGE
> ZIM$TD <- as.Date(as.character(ZIM$TD), "%d/%m/%Y")
> ZIM$CALVING_DATE <- as.Date(as.character(ZIM$CALVING.DATE), "%d/%m/%Y")
> ZIM<-ZIM[order(ZIM$ANIMAL_NUMBER,ZIM$CALVING_DATE,ZIM$TD),]
> ZIM$AGE[ZIM$AGE > 1465 & ! is.na(ZIM$AGE)] <- NA
> ZIM$AGE[ZIM$AGE < 1099 & ! is.na(ZIM$AGE)] <- NA
> summary(ZIM$AGE)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   1099    1164    1233    1251    1328    1465  152213 
> summary(L1$AGE)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  491.0   740.0   796.0   824.7   886.0  1098.0   96193 
> L2<-ZIM
> summary(L2$AGE)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   1099    1164    1233    1251    1328    1465  152213 
ZIM$AGE[ZIM$AGE > 1830 & ! is.na(ZIM$AGE)] <- NA
ZIM$AGE[ZIM$AGE < 1465 & ! is.na(ZIM$AGE)] <- NA

########## This ORDERS cows by calving dates, test dates

ZIM<-ZIM[order(ZIM$ANIMAL_NUMBER,ZIM$CALVING_DATE,ZIM$TD),] # ORDERING CD, TD BY ANIMAL NUMBER
Q<-sum(! is.na(match(dope$ANIMAL_NUMBER, names(y[y == 1]))))
sum(! is.na(match(Lact1$Cow, names(y[y == 1]))))
#####################
y<- tapply(dope$CALVING.DATE, dope$ANIMAL_NUMBER, findunique) # order cows by calving dates
summary(factor(y))
Lact1 <- dat[! is.na(match(dat$ANIMAL_NUMBER, names(y[y == 1]))),] # match dataD with all groups of y=1 (lactation 1) with data frame
head(Lact1)
Lact2 <- dat[! is.na(match(dat$ANIMAL_NUMBER, names(y[y == 2]))),]

######################################################
dealing with dates
ZIM$CALVING_DATE <- as.Date(as.character(ZIM$CALVING.DATE), "%d/%m/%Y")
ZIM$TD <- as.Date(as.character(ZIM$TD), "%d/%m/%Y")
#############################################################
ANI$DOB<-as.Date(as.character(ANI$BIRTH.DATE),format="%d-%b-%y") #formating abbreviated month to number in R
#####################################
calculate 305D yield for mapping_310316 ~ zimdata ~ Lact1,Lact2,Lact3, Lact4 and Lact5(D1,D2,D3, and D4)
get.lmy <- function(D4, refd){

 cows <- unique(D4$ANIMAL_NUMBER)

 ncows <- length(cows)

 ## ## ##

D4$lmy <- NA
 
 D4$dim <- NA

 for(i in 1:ncows){

   ## print.noquote(paste(i, date()))

   test <- (D4$ANIMAL_NUMBER == cows[i])

   tmp <-  D4[which(D4$ANIMAL_NUMBER == cows[i]),]
   
   alldates <- c(tmp$CALVING_DATE[1], tmp$TD)

   ivals <- as.numeric(as.character(diff(alldates)))

   ## ################
   ## New code added 28 May 2015 to remove obs whose time from lactation is greater than 'refd'
   
   keepvals <- (cumsum(ivals) < refd)

   tmp <- tmp[keepvals,]

   alldates <- c(tmp$CALVING_DATE[1], tmp$TD)

   ivals <- as.numeric(as.character(diff(alldates)))
   
   ## #################
   
   n <- length(ivals)

   D4$dim[test] <- sum(ivals)
   
   ivals <- c(ivals, refd - sum(ivals))

   mvals <- c(tmp$MILK.YIELD[1], (tmp$MILK.YIELD[-n] + tmp$MILK.YIELD[-1])/2, tmp$MILK.YIELD[n])
   
   D4$lmy[test] <- sum(ivals * mvals)
 }
 
D4
}
 D4 <- get.lmy(D4, refd=305)
##############################
C2<-D2[!duplicated(D2$ANIMAL_NUMBER),] #  removing duplicate 305d milk yield and records
#################################
DF<-Lact1[!duplicated(Lact1$ANIMAL_NUMBER),] # remove duplicate records and animals to as to combine from parity 1 to 5

tmp$SEASON <- factor(as.character(tmp$CM), levels=1:12)
 summary(tmp$SEASON)
levels(tmp$SEASON ) <- c(2,2,2,1,1,1,1,1,1,1,2,2)

summary(B3$SEASON)
#######################################
k<-1

B3$ANIMAL_NUMBER<-as.factor(B3$ANIMAL_NUMBER)
B3$CI <- NA


ANI<-levels(B3$ANIMAL_NUMBER)

for(i in ANI){

  #print.noquote(i)
browse()
  ikeep <- (B3$ANIMAL_NUMBER==i)

  ikeep1 <- ikeep & (B3$PARITY == 1)
  ikeep2 <- ikeep & (B3$PARITY == 2)
  ikeep3 <- ikeep & (B3$PARITY == 3)
  ikeep4 <- ikeep & (B3$PARITY == 4)
 

  if(sum(ikeep1) > 0){

    if(sum(ikeep2) > 0){

      B3$CI[ikeep1] <- as.numeric(B3$CALVING_DATE[ikeep2] - B3$CALVING_DATE[ikeep1])
    }
  }

  if(sum(ikeep2) > 0){

    if(sum(ikeep3) > 0){

      B3$CI[ikeep2] <- as.numeric(B3$CALVING_DATE[ikeep3] - B3$CALVING_DATE[ikeep2])
    }
  }

  if(sum(ikeep3) > 0){

    if(sum(ikeep4) > 0){

      B3$CI[ikeep3] <- as.numeric(B3$CALVING_DATE[ikeep4] - B3$CALVING_DATE[ikeep3])
    }
  }
}
###########
animal appearing in lact 1 and lact 2
pa<-tapply(dat$PARITY, dat$ANIMAL_NUMBER, function(x) length(table(x))) ##hos is HF_SA HYSlactationfile
table(pa)# list of ani_id with 1, 2, 3, or more lact per animal
table(pa==5)
pa5<-names(pa[pa==5]) ##to pull out names of  cows from pa to pa5; NAMES of all ani_id in lact1-5
p5<-hos[hos$ANIMAL_NUMBER%in%pa5, ]
p5<-hos[hos$ANIMAL_NUMBER%in%as.numeric(pa5), ]; ##because the ani_ids are numeric; maybe i should change to factor or do as.numeric. 
					p5 has the RECORD LIST of all ani_id in all parities
length(pa5)*5  # multiplied by 5 since there are 5 lactations per animal
#######################################################################################################################################################################
 Y88<-subset(JOINT,JOINT$YOB=="2010")
> mean(Y88$EBV)
[1] 17.22512
> sd(Y88$EBV)
[1] 95.12092
> sd(Y88$EBV,na.rm=T)/sqrt(nrow(Y88))
[1] 27.45905
> Y88<-subset(JOINT,JOINT$YOB=="2011")
> sd(Y88$EBV,na.rm=T)/sqrt(nrow(Y88))
[1] 44.16427
> mean(Y88$EBV)
[1] -20.45667

 z1$ANIMAL_NUMBER <- factor(z1$ANIMAL_NUMBER, levels=unique(z1$ANIMAL_NUMBER))               #Prevent reordering of character factors
 
Freqs<-count(z1$ANIMAL_NUMBER)                                                    #Count frequencies of appearance of animals
 Uniques<-unique(z1$ANIMAL_NUMBER)                                                 #Get unique animal IDs
 
 Day<-rep(NA,length=12*length(unique(z1$ANIMAL_NUMBER)))                           #Create vectors
 Record<-rep(NA,length=12*length(unique(z1$MILK_YIELD)))
 Cows<-rep(NA,length=12*length(unique(z1$ANIMAL_NUMBER)))
 
 k<-0                                                                         #Counter for records in new table
 t<-0                                                                         #Counter for records in original table
 for(i in 1:length(Uniques)){
     l<-0
     for(j in 1:Freqs$freq[i]){
       k<-k+1
       t<-t+1
       l<-l+1
       Cows[k]<-as.character(z1$ANIMAL_NUMBER[t])                                                    #Create new records
       Record[k]<-z1$MILK_YIELD[t]
       Day[k]<-l
     }
     if((12-Freqs$freq[i])>0){                                                #Create NA records
       for(j in 1:(12-Freqs$freq[i])){
         k<-k+1
         l<-l+1
         Cows[k]<-as.character(z1$ANIMAL_NUMBER[t])
         Record[k]<-NA
         Day[k]<-NA
       }
     }
 }
 check<-cbind.data.frame(Cows,Record,Day)  

d1<-merge(check,Z1,by.x=c("Cows"),by.y=c("ANIMAL_NUMBER"))


ymp <-merge(P, by.x=c("SIRE_ID"),by.y=c("SIRE_ID"))

ymp <-merge(G,P, by.x=c("SIRE_ID"),by.y=c("SIRE_ID"))
#####################################################

PEDIGREE ERRORS FOR ZIMBABWE
1790850 appearing as both sire and dam.
207184648 appearing more than once in the progeny list.
3831000 appearing more than once in the progeny list.
9098197 appearing as both sire and dam.
5470562 appearing more than once in the progeny list.
5470579 appearing more than once in the progeny list.
5470562 appearing as both sire and dam.
5470562 appearing as both sire and dam.
60540099 appearing more than once in the progeny list.
1884404 appearing as both sire and dam.
8953097 appearing more than once in the progeny list.
KNIGHT appearing as both sire and dam.
S appearing as both sire and dam.



 summary(B3$CI)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    0.0   378.0   429.0   457.6   506.0  1655.0  101916 
> summary(B3$CI2)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    0.0   370.0   412.0   428.6   478.0   718.0  165474 
> summary(B3$CI3)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  132.0   378.0   435.0   498.6   541.0  1655.0  174930 
> summary(B3$CI1)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  113.0   382.0   436.0   460.1   513.0   910.0  136506 
