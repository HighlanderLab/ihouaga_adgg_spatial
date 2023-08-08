
#Loading dplyr
library(dplyr)
# Import data
pheno<- read.table(file = "data5.dat", header = FALSE)
colnames(pheno) <- c("cow","milk","ward","herd","cyrsn","tyrmn","dgrp","lac","lacgr","age","long","lat", "intercept", "leg1","leg2", "mean", "ward_code", "Region_Cod")
head(pheno)
dim(pheno)
cow

# Import SNP data
geno<- read.table(file = "snpref-isi.txt", header = FALSE)
names(geno)[1] <- "cow"

#In GBLUP we can only work with genotyped animals - non-genotyped animals will be removed (see later),

# hence we will encode genotyped animals as 1:n
dim(geno)
geno$IId <- 1:nrow(geno)
names(geno)[1]
# Remove phenotypes for animals that are not genotyped
# (I am assuming that geno and pheno data.frames have Id columns
sel <- pheno$cow %in% geno$cow
pheno <- pheno[sel, ]

# Recode phenotyped animals into 1:n according to geno 1:n codes
sel <- match(pheno$cow, table = geno$cow)
pheno$IId <- geno$IId[sel]

# Save geno$cow and geno$IId as backup
genocow_IId <- subset(geno,select= c(IId,cow))
write.table(genocow_IId,"genocow_IId.txt")
# Remove old genoid and make

 
summary(k)
if(pheno$cow==pheno$IId){
  print("YES")
} else {
  print("No")
}
geno<- subset(geno,select= -c(cow))

#Make IId (last column as first)
geno <- geno[,c(ncol(geno),1:(ncol(geno)-1))]

#Remove all column  names

names(geno)<- NULL

head(geno, n=2)
head(pheno, n=2)

tail(geno, n=2)
tail(pheno,n=2)

#save pheno and geno
write.table(pheno, "pheno.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
write.table(geno, "geno.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")


