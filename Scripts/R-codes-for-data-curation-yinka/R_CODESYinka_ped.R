genoid<-read.table("pd_170216.txt",h=T)
dim(genoid)
head(genoid)
colnames(genoid)<-c("newid","newsire","newdam","id","Sire","Dam")
colnames(pd)<-c("newid","newsire","newdam","id","Sire","Dam")
pd<-genoid[,-c(4,5,6)]

pd<-genoid[,-c(4,5,6)]
dim(pd)
head(pd)

write.table(pd,"newped.txt",quote=F,row.names=F,col.names=F)

phen<-read.table("data.txt",h=T)
dim(phen)
head(phen)
######
removing default row names and column names in R (issue with joint holstein friesian pedigree)
write.table(newdata,"newdata.txt",row.names=F,col.names=F)
#############################
data2<-merge(phen,genoid, by.x=c("id"),by.y=c("oldid")) 
head(data2)
dim(data2)
#########
PD<-read.table("joint_hf3_Rtools_PD.txt",sep="",h=F)
names(PD)<-c("ANIMAL_NUMBER","SIRE_ID","DAM_ID","new_id","new_sire","new_dam")
data2<-merge(HF3_040516,PD, by.x=c("ANIMAL_NUMBER"),by.y=c("new_id")) 
##########
      
new<-merge(PED3,pedred, by.x=c("id"),by.y=c("id"), all.y=T) 


write.table(data2,"newpheno.txt",quote=F,row.names=F,col.names=T)
write.csv(data2,"newpheno.csv",quote=F,row.names=F)

### WRITE .txt, .csv TO FILE WITHOUT QUOTES, AND ROW NAMES
write.table(JOINT_HF_3COUNTRIES,"JOINT_HF_3COUNTRIES.txt",quote=FALSE,row.names=FALSE,col.names=T)

colnames(pede)<-c("ANIMAL_NUMBER", "SIRE_ID", "DAM_ID")
###################################
library(reshape2)
Acols<-read.table("amatrix.out", header=F) #can’t remember if the amatrix has  a header, if it does you need to change from F to T
colnames(Acols)<-c("i","j","relationship")

Amat<-acast(Acols, i ~ j, value.var = "relationship") #acast is from the reshape package and helps change the shape of your data
Amat2<-Amat
Amat2[upper.tri(Amat2)] <- (Amat2)[upper.tri(Amat2)]

View(Amat2) #its worth doing this to check that your matrix is symmetrical.

