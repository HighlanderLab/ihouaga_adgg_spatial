convert dates month from letters to numbers
CALVING_DATE<-as.Date(Zim$CALVING_DATE,"%d-%m-%y")
Zim["CALVING_DTM"] <- NA
Zim$CALVING_DTM<-as.Date(Zim$CALVING_DATE,"%d-%b-%y")
Zim$CALVING_DATE<-format(Zim$CALVING_DTM,"%d/%m/%Y")
head(Zim)


TESTDAYDATE<-as.Date(Zim$TESTDAYDATE,"%d-%m-%y")
ZIM_LACTATION["TESTDAY_DTM"] <- NA
ZIM_LACTATION$TESTDAY_DTM<-as.Date(ZIM_LACTATION$TESTDAYDATE,"%d-%b-%y")
ZIM_LACTATION$TESTDAYDATE<-format(ZIM_LACTATION$TESTDAY_DTM,"%d/%m/%Y")
head(ZIM_LACTATION)

Zim$MLK_YLD[Zim$MLK_YLD > 19.3 & ! is.na(Zim$MLK_YLD)] <- NA
Zim$MLK_YLD[Zim$MLK_YLD < 1.89 & ! is.na(Zim$MLK_YLD)] <- NA
ZIM_LACTATION<-na.omit(Zim)

#################################
# order   dates
 ZIM_LACTATION<- ZIM_LACTATION[order(ZIM_LACTATION$ANIMAL_NUMBER, ZIM_LACTATION$TESTDAYDATE),]
ZIM_LACTATION <- ZIM_LACTATION[order(ZIM_LACTATION$ANIMAL_NUMBER, ZIM_LACTATION$CALVING_DATE),]
ZIM_LACTATION$TESTDAYDATE<- as.Date(as.character(ZIM_LACTATION$TESTDAYDATE), "%d/%m/%Y")

#####
ZIM_LACTATION$DIM <- as.numeric(as.character(diff(ZIM_LACTATION$TESTDAYDATE,ZIM_LACTATION$CALVING_DATE)))
as.numeric(ZIM_LACTATION$TESTDAYDATE-ZIM_LACTATION$CALVING_DATE)

ZIM_LACTATION$DIM = difftime(ZIM_LACTATION$TESTDAYDATE[2], ZIM_LACTATION$CALVING_DATE[1], units = "days") # days
> diff_in_days
Time difference of 435.9583 days

############
Data<- Data[order(Data$ANIMAL_NUMBER, Data$TESTDAY_DTM, Data$CALVING_D),]
#################################
MILK_305D<-Lactation$Milk_305 
PROTEIN_305D<-Lactation$Protein_305 
FAT_305D<-Lactation$Fat_305
Milk_yld<-Lactation$MILK_YLD                             #Save a vector with ages
Protein_yld<-Lactation$PROTEIN_CONC
Fat_yld<-Lactation$FAT_CONC
IDs<-Lactation$ANIMAL_NUMBER 
Herd_ID<-Lactation$HERD_ID
DIM<-Lactation$DIM
CALVING_YEAR<-Lactation$CALVING_YR                 
MILK_305D<-vector(length=length(IDs))		#Create  a vector for all animals  that will store the mean of MILK_YLD_305D
PROTEIN_305D<-vector(length=length(IDs))	#Create  a vector for all animals  that will store the mean of PROTEIN_YLD_305D
FAT_305D<-vector(length=length(IDs))	#Create  a vector for all animals  that will store the mean of FAT_YLD_305D
Milk_yld<-vector(length=length(IDs))
Protein_yld<-vector(length=length(IDs))
Fat_yld<-vector(length=length(IDs))
CALVING_YEAR<-vector(length=length(IDs))
DIM<-vector(length=length(IDs))
for(i in 1:length(IDs)){                            #Loop for all animals in ID vector
     
MILK_305D[i]<-mean(Lactation$Milk_305[which(Lactation$Milk_305==IDs[i])]) 
PROTEIN_305D[i]<-mean(Lactation$Protein_305[which(Lactation$ANIMAL_NUMBER ==IDs[i])])
FAT_305D[i]<-mean(Lactation$Fat_305[which(Lactation$ANIMAL_NUMBER==IDs[i])])       #Compute total as the sum of NEWYIELD and removing NAs
Milk_yld[i]<-mean(Lactation$MILK_YLD[which(Lactation$ANIMAL_NUMBER==IDs[i])])           #Compute total as the sum of NEWYIELD and removing NAs
Protein_yld[i]<-mean(Lactation$PROTEIN_CONC[which(Lactation$ANIMAL_NUMBER==IDs[i])])
Fat_yld[i]<-mean(Lactation$FAT_CONC[which(Lactation$ANIMAL_NUMBER==IDs[i])])
DIM[i]<-sum(Lactation$DIM[which(Lactation$ANIMAL_NUMBER==IDs[i])])
Herd_ID<-!duplicated(Lactation$HERD_ID==IDs[i])
CALVING_YEAR<-!duplicated(Lactation$CALVING_YR==IDs[i])
}

Lactation_305<-array(0,dim=c(length(IDs),10))
Lactation_305[,1]<-IDs
Lactation_305[,2]<-MILK_305D
Lactation_305[,3]<-Milk_yld
Lactation_305[,4]<-PROTEIN_305D
Lactation_305[,5]<-Protein_yld
Lactation_305[,6]<-FAT_305D
Lactation_305[,7]<-Fat_yld
Lactation_305[,8]<-DIM
Lactation_305[,9]<-CALVING_YR
Lactation_305[,10]<-Herd_ID
head(Lactation_305)
##########
dota_sire<-jer[duplicated(jer[,c("ANIMAL_NUMBER", "SIRE_ID")])==FALSE,]
nrow(dota_sire)
dotaper_sire<-as.data.frame(table(dota_sire$SIRE_ID))
head(dotaper_sire)
##########################
zim_herd50<-dear[duplicated(dear[,c("ANIMAL_NUMBER", "HERD_ID")])==FALSE,]
nrow(zim_herd50)
anim_per_Participant<-as.data.frame(table(zim_herd50$HERD_ID))
head(anim_per_Participant)
number of herds with more than/equal to 50 jersey cows
nrow(anim_per_Participant[which(anim_per_Participant$Freq>=50),])
####################################
dear$HERD_ID <- factor(as.character(dear$HERD_ID))
ncows <- tapply(dear$ANIMAL_NUMBER, dear$HERD_ID, function(x){length(unique(x))})
herdkeep <- names(ncows)[ncows >= 50]
cumsum(herdkeep) ## it gave me the herd count as 25 herds with 50>= animals
***   new<-dear[which(dear$HERD_ID %in% herdkeep),] *****
new<-dat[which(dat$PARTICIPANT %in% herdkeep),]
#########################################
ZIM$SEASON <- factor(as.character(ZIM$BIRTH_MTH), levels=1:12)
summary(ZIM$SEASON)
levels(ZIM$SEASON) <- c(1,1,1,1,2,2,2,2,1,2,1,1)

##### calving season
ZIM_LACT030816$C_SEASON <- factor(as.character(ZIM_LACT030816$CM), levels=1:12)
summary(ZIM_LACT030816$C_SEASON)
levels(ZIM_LACT030816$C_SEASON) <- c(2,2,2,2,1,1,1,1,2,1,2,2)
#######################################
** intervals between test day records
NEW$ANIMAL_NUMBER <- factor(NEW$ANIMAL_NUMBER)

nrecs.per.cow <- tapply(NEW$TD, NEW$ANIMAL_NUMBER, function(x){length(unique(x))})

NEW$Nrecs <- nrecs.per.cow[match(NEW$ANIMAL_NUMBER, levels(NEW$ANIMAL_NUMBER))]

ZIM_NEW <- NEW[NEW$Nrecs >= 6,]

ZIM_NEW$Interval <- NA

for(i in 1:length(levels(NEW$ANIMAL_NUMBER))){

  zz <- (ZIM_NEW$ANIMAL_NUMBER == levels(NEW$ANIMAL_NUMBER)[i])

  if(sum(zz) > 0){

    tmp <- ZIM_NEW[zz,]

    tmp$Interval <- c(NA, as.numeric(diff(tmp$TD)))

    ZIM_NEW$Interval[zz] <- tmp$Interval
  }
}
##
#######################################
3 or more HYS combinations
new$HHYS<-as.factor(new$HHYS)                                    #Convert as factor
niv<-levels(new$HHYS)                                           #List of levels of the factor
extracted<-as.data.frame(which(table(new$HHYS)>=3))[,1]         #List of levels that appear three or more times
rep_HYS<-vector(length=length(extracted))                      #Create vector to store HYS levels that appear three or more times
for (i in 1:length(extracted)){                                 #Loop to go through the entire list of repeated levels                       
  rep_HYS[i]<-niv[extracted[i]]                                #The levels we want is the one which is: niv[listed level in extracted]
}
head(rep_HYS)
############  3 or more HYS combinations. the below (zim_hys) seems more correct than rep_HYS
ncow.combi <- tapply(new$ANIMAL_NUMBER, new$HHYS, function(x){length(unique(x))})

ncow.dat <- data.frame(HHYS = names(ncow.combi), Ncow = as.numeric(ncow.combi))

write.csv(ncow.dat, file="D:\\Users\\OAbejide\\Zimbabwe Data\\ncow.combi.csv")

combikeep <- names(ncow.combi)[ncow.combi >= 3]

recordkeep.combi <- ! is.na(match(new$HHYS, combikeep))

sum(recordkeep.combi)
sum(! recordkeep.combi)

zim_hys<- new[recordkeep.combi,]

 #######
 IN ORDER TO "SORT" A PEDIGREE BASE PARENT TO PROGENY, THE FOLLOWING HAVE TO BE DONE!!!
 INSTALL PACKAGES#########
 install.packages("pedigree", lib="D:/Users/OAbejide/R_packages")
 install.packages("lme4", lib="D:/Users/OAbejide/R_packages")
 
 install.packages("reshape", lib="D:/Users/OAbejide/R_packages")
 install.packages("Rcpp", lib="D:/Users/OAbejide/R_packages")
 install.packages("dplyr", lib="D:/Users/OAbejide/R_packages")
 library(Rcpp, lib.loc="D:/Users/OAbejide/R_packages")
 library(pedigree, lib.loc="D:/Users/OAbejide/R_packages")
 library(Hmisc, lib.loc="D:/Users/OAbejide/R_packages")
 library(dplyr, lib.loc="D:/Users/OAbejide/R_packages")
 library(lme4, lib.loc="D:/Users/OAbejide/R_packages")
 install.packages("plyr", lib="D:/Users/OAbejide/R_packages")
 install.packages("Hmisc", lib="D:/Users/OAbejide/R_packages")
 install.packages("nlme", lib="D:/Users/OAbejide/R_packages")
  library(nlme, lib.loc="D:/Users/OAbejide/R_packages")
 #########
 install.packages("plyr", lib="D:/Users/OAbejide/R_packages") # tool for Splitting, Applying and Combining data
 library(plyr, lib.loc="D:/Users/OAbejide/R_packages")
 * AIC and BIC methods
 install.packages("AICcmodavg", lib="D:/Users/OAbejide/R_packages")
  library(AICcmodavg, lib.loc="D:/Users/OAbejide/R_packages")
 
 for correlation tables
 install.packages("picante", lib="D:/Users/OAbejide/R_packages")
 library(picante, lib.loc="D:/Users/OAbejide/R_packages")

 ################################
ord <- order(orderPed(ped))
ped1<- ped[ord, ]
Ped$gen<-countGen(ped) # count genetations
table(Ped$gen) # individuals per generation
########
match sire and daughetrs in the file
df1$match <- match(df1$location, df2$location, nomatch=0)
for instance df$id(sire) with df$cowid, nomatch=0)
sire (df1)$match <- match(df1$id(sire), df2$id(cow), nomatch=0)
##################
zim_herd<-D[duplicated(D[,c("ANIMAL_NUMBER", "HERD_ID")])==FALSE,]
nrow(zim_herd)
anim_perherd<-as.data.frame(table(zim_herd$HERD_ID))
head(anim_perherd)
number of herds with more than 50 jersey cows
nrow(anim_perherd[which(anim_perherd$Freq>=50),])
########################
anim_perherd
   Var1 Freq
1    14   92
2    77   24
3   125   23
4   179    1
5   209    7
6   210   23
7   287   32
8   288    5
9   353    2
10  469    4
11  557   30
12  627    5
13  641   10
14 1785    2
###############################
ord <- order(orderPed(PED))
ped1<- ped[ord, ]
Ped$gen<-countGen(ped) # count genetations
table(Ped$gen) # individua



ped<-order(orderPed(ped))
ped1<-ped[ord,]
Ped$gen<-countGen(ped)
table(Ped$gen)

############################
CALCULATING 305D YIELD WITH DAYS IN MILK INCLUDED AS "INEW")
get.lmy <- function(Lact1, refd){
 cows <- unique(Lact1$CALF)
ncows <- length(cows)
Lact1$lmy <- NA
   Lact1$ivals <- NA
   for(i in 1:ncows){
    test <- (Lact1$CALF == cows[i])
    tmp <-  Lact1[test,]
   alldates <- c(tmp$CALVING_DATE[1], tmp$TD)
   ivals <- as.numeric(as.character(diff(alldates)))
 ivals <- as.numeric(as.character(diff(c(tmp$CALVING_DATE[1], tmp$TD))))
  n <- length(ivals)
   ivals <- c(ivals, refd - sum(ivals))
    mvals <- c(tmp$MILK_YIELD[1], (tmp$MILK_YIELD[-n] + tmp$MILK_YIELD[-1])/2, tmp$MILK_YIELD[n])
   Lact1$lmy[test] <- sum(ivals * mvals)
   inew <- cumsum(ivals[-length(ivals)])
   Lact1$ivals[test] <- inew[length(inew)]
  }
  Lact1
 }
Lact2$CALF <- factor(as.numeric(factor(paste(Lact2$ANIMAL_NUMBER, Lact2$CALVING_DATE))))
Lact2 <- get.lmy(Lact2,refd=305)
#################################################################################################

ord<-orderPed(test)
> head(ord)
[1] 455 456 457 458 459 460
> pedigree <- test[order(ord),]
> head(pedigree)
         Progeny       Sire        Dam
18560 2011412075          0          0
14097 2011338601          0          0
15013 2011352750 2011412075 2011338601
14106 2011338718          0          0
14094 2011338551          0          0
11312   63300099 2011338718 2011338551

> pedigree$gen<-countGen(pedigree)
> table(pedigree$gen)

     0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20     21     22     23     24     25     26     27 
198887  39906  20919  11967   4661   2142   1898   8704  50134  74539  98583  77572  65700  59969  60881  66941  59316  49083  44805  27724  14321   6758   2307    658    152     31     11      4 
    28 
     2 
> TRIM/ REDUCE PEDIGREE

ped<-trimPed(pedigree, gmp,ngenback = NULL)

STANDARD ERROR IN R
y11<-subset(jer,jer$CALVING_YR=="2011")
> nrow(y11)
[1] 3
> sd(y11$MILK_YLD_305D,na.rm=T)/sqrt(3)
[1] 235.8496


Correlations between testday yield and 305day yield
converted TD to as.date and then to as.matrix as it said "x" must be numeric while X305D was numeric
cor(as.matrix(Z1$TD,Z1$X305D))
#####
>  cor(Z1$MILK_YIELD,Z1$X305D,use="complete.obs",method="pearson")
[1] 0.685208
>  cor(Z2$MILK_YIELD,Z2$X305D,use="complete.obs",method="pearson")
[1] 0.6519077
cor(Z1$X305D,Z1$MILK_YIELD,use="complete.obs",method="kendall")
>[1] 0.5280755
cor(Z2$MILK_YIELD,Z2$X305D,use="complete.obs",method="kendall")
[1] 0.3898499
cor(Z1$X305D,Z1$MILK_YIELD,use="complete.obs",method="spearman")
[1] 0.6942949
>  cor(Z2$X305D,Z2$MILK_YIELD,use="complete.obs",method="spearman")
[1] 0.5005425

remove default row names in R:  
write.table(newdata,"newdata.txt",row.names=F,col.names=F)
### PLYR was installed in order to do a correlation between test days and 305D yield for milk records in Zim data specifically (Z1,Z2 and mapping051115.csv)
install.packages("plyr", lib="D:/Users/OAbejide/R_packages") # tool for Splitting, Applying and Combining data
 library(plyr, lib.loc="D:/Users/OAbejide/R_packages")
 
 library(pedigree, lib.loc="D:/Users/OAbejide/R_packages")
library(plyr)                                                                #Library needed to count

Z1$ANIMAL_NUMBER <- factor(Z1$ANIMAL_NUMBER, levels=unique(Z1$ANIMAL_NUMBER))               #Prevent reordering of character factors

Freqs<-count(Z1$ANIMAL_NUMBER)                                                    #Count frequencies of appearance of animals
Uniques<-unique(Z1$ANIMAL_NUMBER)                                                 #Get unique animal IDs

Day<-rep(NA,length=11*length(unique(Z1$ANIMAL_NUMBER)))                           #Create vectors
Record<-rep(NA,length=11*length(unique(Z1$MILK_YIELD)))
Cows<-rep(NA,length=11*length(unique(Z1$ANIMAL_NUMBER)))

k<-0                                                                         #Counter for records in new table
t<-0                                                                         #Counter for records in original table
for(i in 1:length(Uniques)){
    l<-0
    for(j in 1:Freqs$freq[i]){
      k<-k+1
      t<-t+1
      l<-l+1
      Cows[k]<-as.character(Z1$ANIMAL_NUMBER[t])                                                    #Create new records
      Record[k]<-Z1$MILK_YIELD[t]
      Day[k]<-l
    }
    if((11-Freqs$freq[i])>0){                                                #Create NA records
      for(j in 1:(11-Freqs$freq[i])){
        k<-k+1
        l<-l+1
        Cows[k]<-as.character(Z1$ANIMAL_NUMBER[t])
        Record[k]<-NA
        Day[k]<-NA
      }
    }
}
Newdatatable<-cbind.data.frame(Cows,Record,Day)                              #Create new data table

To get information of Z1 to Newdatatable, Merge in R was used###
D1<-merge(Newdatatable,Z1,by.x=c("Cows"),by.y=c("ANIMAL_NUMBER"))

 cor(D1$Day,D1$X305D,use="complete.obs",method="pearson")