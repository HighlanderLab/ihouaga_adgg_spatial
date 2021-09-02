after reading cow data to do cow EBvs, there are duplicated rows
 c<-KJCow_EBV[!duplicated(KJCow_EBV$CowID_recode),] # remove duplicated rows
 for MY rel in Kenya
 c$MYRel <- 1 - (c$V4^2 / (34800.7))  # GV is  34800.7
 c$MYRel[ c$MYRel <= 0] <- 0.01 #where min Rel is less than 0.01
  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.008019 0.010000 0.017600 0.052850 0.097980 0.145900
#######
y76<-subset(Joint,Joint$YOB=="2000")
mean(y76$JCI1,na.rm=T)
sd(y76$JCI1,na.rm=T)/sqrt(nrow(y76))
mean(y76$JCI1Rel,na.rm=T)
sd(y76$JCI1Rel,na.rm=T)/sqrt(nrow(y76))
###
S771$DOB<-format(as.Date(as.character(S771$BIRTH_DATE), format="%Y%m%d"),"%d/%m/%Y") # formatting dates from 19811116 to 16/11/1981
 for AFC (of 1342days) rel in Kenya
AFC_EBV$AFCRel <- 1 - (AFC_EBV$V4^2 / (0.01)) # GV is 0.01 (using normalised AFC)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.03489 0.22790 0.22790 0.22730 0.22790 0.22790 

 for CI1 (of 1342days) rel in Kenya
AFC_EBV$CI1Rel <- 1 - (AFC_EBV$V4^2 / (0.01)) # GV for CI1 is 0.01 (using normalised AFC)
AFC_EBV$CI1Rel[ AFC_EBV$CI1Rel <= 0] <- 0.01 
  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.01000 0.07878 0.07878 0.07858 0.07878 0.07878 
############################# SOUTH AFRICA JERSEY using new data.
SA_JERMYEBV
SA_JERMYEBV$MYRel <- 1 - (SA_JERMYEBV$V4^2 / (159725)) #GV for MY (159725)
  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.1275  0.3736  0.4233  0.4096  0.4585  0.5567 

P1_EBV$AFCRel <- 1 - (P1_EBV$V4^2 / (4483.63)) #GV for AFC (4483.63)
P1_EBV$AFCRel[P1_EBV$AFCRel <= 0] <- 0.01
Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0002961 0.4446000 0.4968000 0.4759000 0.5289000 0.5947000 
###
P1_EBV$CI1Rel <- 1 - (P1_EBV$V4^2 / (239.408)) #GV for CI1 (239.408) 
P1_EBV$CI1Rel[P1_EBV$CI1Rel <= 0] <- 0.01
Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0003638 0.0601800 0.0974000 0.1095000 0.1459000 0.2919000 
##############JOINT JERSEY
JER_EBV$MYRel <- 1 - (JER_EBV$V4^2 / (167492)) #GV for Joint MY (167492)
  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.02746 0.37980 0.43030 0.41670 0.46700 0.56830 
#####
J1342_JAFC$AFCRel <- 1 - (J1342_JAFC$V4^2 / (92.4193)) #GV for Joint AFC (92.4193)
J1342_JAFC$AFCRel[J1342_JAFC$AFCRel <= 0] <- 0.01
Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0003111 0.0312700 0.0554800 0.0707400 0.0929100 0.2278000 
#####
J1342_JAFC$CI1Rel <- 1 - (J1342_JAFC$V4^2 / (318.688)) #GV for Joint CI1 (318.688)
J1342_JAFC$CI1Rel[J1342_JAFC$CI1Rel <= 0] <- 0.01
Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0002055 0.0760100 0.1196000 0.1299000 0.1755000 0.3219000 
#####################################################################################################################
SIRE RELIABILITY
Kenya Jersey