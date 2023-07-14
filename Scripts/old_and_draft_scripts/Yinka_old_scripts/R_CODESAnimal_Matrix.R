mat$Country <- 2
mat$Country[1:98399] <- 1

mat$ida <- factor(paste(mat$animal_J, mat$Rship))
mat$idb <- factor(paste(mat$animal_J, mat$Rship, mat$Country))

lida <- summary(mat$ida, maxsum=1000000)
lidb <- summary(mat$idb, maxsum=1000000)

mat$n.all <- as.numeric(lida)[match(mat$ida, names(lida))]
mat$n.this <- as.numeric(lidb)[match(mat$idb, names(lidb))]

mat$ida <- mat$idb <- NULL

mat$n1 <- mat$n.this * (mat$Country == 1) + (mat$n.all - mat$n.this) * (mat$Country == 2)
mat$n2 <- mat$n.this * (mat$Country == 2) + (mat$n.all - mat$n.this) * (mat$Country == 1)

to read a text file from cfc 
mat<-read.table("D:\\Users\\OAbejide\\KENYA&SA_HF\\SAKenya_Hols_Amatrix.txt", header=FALSE, fill=TRUE)



mean(y91$AGE)/28
sd(y91$AGE)/28
sd(y91$AGE)/28/sqrt(nrow(y91))
length(unique(y91$ANIMAL_NUMBER))

#####
MAT$idb <- factor(paste(MAT$animal_J, MAT$Rship, MAT$Country))
lidb <- summary(MAT$idb, maxsum=1000000)
MAT$n.all <- as.numeric(lidb)[match(MAT$idb, names(lidb))]

####################### Create output matrix file in chunks so as to handle easily

f <- file('amatrix.out', 'r')
header <- readLines(f, 1)
icon=1

while (TRUE) {
  out <- read.table(f, nrows=30000000, header=FALSE, colClasses=list('integer','integer','numeric'),col.names=c('ID1','ID2','rel'))
  name<-paste("data",icon)
  write.table(out,file=name) 
  icon=icon+1
  if (nrow(out) == 0) break
}
  
  m <- merged(out, mat_out, all.x=TRUE, all.y=FALSE)
  write.table(m, file='merged.out', col.names=FALSE, append=TRUE)

 
}
close(f)

###############################
size of matrix is too big so working with a smaller size to have maybe 4 generations
 read recoded file into R
  colnames(recoded)<-c("ANIMAL_NUMBER", "SIRE_ID","DAM_ID","newid", "newSire", "newDam")
  read common sire file
  head(sire)
         SIRE_ID
1  GBRM0000392457
2 USAM00001583197
3 USAM00001620273
4 USAM00001626813
5 USAM00001667366
6 USAM00001723741

 ped_gen0<-recode[which(sire$SIRE_ID%in%recode$oldid),]
> nrow(ped_gen0)
[1] 53
 ped_gen1<-recode[which(recode$SIRE_ID%in%ped_gen0$oldid),]
> tail(ped_gen1)
           oldid         SIRE_ID     DAM_ID newid newSire newDam
75972 7020646280  GBRM0000392457    5550566 75314   16152  22281
75973 1020200850 CANM00005457798    5550596 57106    1664  22286
75976 2020683301 USAM00002203706    5231177 73701    1090  20939
75980 6020820690 USAM00002124357 4020809190 75213    1066  33336
75985 8020811830 USAM00006215489 9020279480 75725    1946  75724
75992 3020811738 USAM00002271271   20683493 74238   12418  32947
> nrow(ped_gen1)
[1] 4574
ped_gen2<-recode[which(recode$SIRE_ID%in%ped_gen1$oldid | recode$DAM_ID%in%ped_gen1$oldid),]
> nrow(ped_gen2)
[1] 11593
> 11593+4574+53 = 16220 animals in generation 0 to generation 2

Joint Jersey  pedigree###################################
Jer_gen0<-recoded[which(recoded$ANIMAL_NUMBER%in%S$SIRE_ID),]
> nrow(Jer_gen0)
[1] 24
head(Jer_gen0)
       ANIMAL_NUMBER SIRE_ID DAM_ID newid newSire newDam
127 DNKM000000302265       *      *   543       0      0
200 CANM000000137656       *      *   839       0      0
218 USAM000000000388       *      *   865       0      0
 Jer_gen1<-recoded[which(recoded$SIRE_ID%in%Jer_gen0$ANIMAL_NUMBER),]
Jer_gen2<-recoded[which(recoded$SIRE_ID%in%Jer_gen1$ANIMAL_NUMBER | recoded$DAM_ID%in%Jer_gen1$ANIMAL_NUMBER),]

#####################
colnames(recoded)<-c("ANIMAL_NUMBER", "SIRE_ID","DAM_ID","newid", "newSire", "newDam")
>  colnames(mat)<-c("newid", "j", "rship")
> data1<-merge(mat,recode,by="newid",all.x=T)

24th februrary 2017
######
HFdota_sire<-HF[duplicated(HF[,c("ANIMAL_NUMBER", "ID_SIRE")])==FALSE,]
nrow(HFdota_sire)
HFper_sire<-as.data.frame(table(HFdota_sire$ID_SIRE))
head(HFper_sire)
SIRE5_dota<-HFper_sire[which(HFper_sire$Freq>=5),]
##############
Zim_sire<-data[duplicated(data[,c("ANIMAL_NUMBER", "ID_SIRE")])==FALSE,]
nrow(Zim_sire)
ZIMper_sire<-as.data.frame(table(Zim_sire$ID_SIRE))
head(ZIMper_sire)
ZIM5_dota<-ZIMper_sire[which(ZIMper_sire$Freq>=5),]
#######################################
Sire5_ped$C_HYS<-as.factor(Sire5_ped$C_AYS)                                    #Convert as factor
niv<-levels(Sire5_ped$C_AYS)                                           #List of levels of the factor
extracted<-as.data.frame(which(table(Sire5_ped$C_AYS)>=3))[,1]         #List of levels that appear three or more times
rep_AYS<-vector(length=length(extracted))                      #Create vector to store HYS levels that appear three or more times
for (i in 1:length(extracted)){                                 #Loop to go through the entire list of repeated levels                       
  rep_AYS[i]<-niv[extracted[i]]                                #The levels we want is the one which is: niv[listed level in extracted]
}

#rep_AYS                                                        #List of herds that appear three or more times

Sire5ped$_3AYS<-merge(Sire5_ped,rep_AYS,by="C_AYS",by.y=TRUE)