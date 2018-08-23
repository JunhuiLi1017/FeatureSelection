#####################################################
## library
#####################################################
library(magrittr)
#install.packages("InformationValue")
library(InformationValue)
#install.packages("gmodels")
library(gmodels)
library(ROCR)
#####################################################
## read data
#####################################################
setwd("E:/work/jiangjian/data/lung/rawData")
#E:\work\jiangjian\data\lung/rawData

#dat <- read.csv("lung.All.log10AFP.YvsC.csv")
dat <- read.csv("lung.lowAFP.YvsC.csv")
head(dat)
dim(dat)
posC <- which(dat[,2] == 0)
posH <- which(dat[,2] == 1)
#?wilcox.test
StatRes <- matrix(NA,7,3)
rownames(StatRes) <- c("Gender/Male:Female","Age",colnames(dat)[c(3,4,5,6,9)])
colnames(StatRes) <- c("Cirrhosis","HCC","p-value")
M.C <- sum(dat[posC,8] %in% 1)
M.H <- sum(dat[posH,8] %in% 1)
StatRes[1,1:2] <- c(paste(M.C,":",length(posC)-M.C,sep=""),paste(M.H,":",length(posH)-M.H,sep=""))



sexhcc <-matrix(c(10,18,5,6),ncol=2,byrow=T)
rownames(sexhcc)<-c("male","female")
colnames(sexhcc)<-c("HCC","Cirrhosis")

StatRes[1,3] <- prop.test(sexhcc)$p.value

StatRes[2,1] <-  paste(round(mean(dat[posC,7]),3)," ¡À ",round(sd(dat[posC,7]),3),sep="")
StatRes[2,2] <-  paste(round(mean(dat[posH,7]),3)," ¡À ",round(sd(dat[posH,7]),3),sep="")
StatRes[2,3] <-  t.test(dat[posC,7],dat[posH,7])$p.value

#i=3
ni = 2
for(i in c(3,4,5,6,9)){
  ni = ni + 1
  Wlc.test <- wilcox.test(dat[posC,i],dat[posH,i],alternative = c("two.sided"))
  #ci(dat[1:24,i])
  #sd(dat[1:24,i])/24^(1/2)
  StatRes[ni,1] <- paste(round(mean(dat[posC,i]),3)," ¡À ",round(sd(dat[posC,i]),3),sep="")
  StatRes[ni,2] <- paste(round(mean(dat[posH,i]),3)," ¡À ",round(sd(dat[posH,i]),3),sep="")
  StatRes[ni,3] <- Wlc.test$p.value
}

#write.csv(StatRes,"StatRes.lowAFP.csv")

#####################################################
## read data  -----------------------2018/07/13
#####################################################
setwd("E:/work/jiangjian/20180713/Cirrhosis.VS.HCC")

dat <- read.csv("data_YHvsHC.csv")
head(dat)
dim(dat)
log10Fe <- data.frame(log10(dat[,5]))
colnames(log10Fe) <- "log10Fer"
dat <- cbind(dat,log10Fe)
posC <- which(dat[,2] == 0)
posH <- which(dat[,2] == 1)
#?wilcox.test
##Build summary information for test
#get information besides age and gender
others <- c(3,4,5,8)
StatRes <- matrix(NA,ncol(dat)-2,3)
rownames(StatRes) <- c("Gender/Male:Female","Age",colnames(dat)[others])
colnames(StatRes) <- c("Cirrhosis","HCC","p-value")

#--* extract gender information
Ngeder <- 7
M.C <- sum(dat[posC,Ngeder] %in% 1)
M.H <- sum(dat[posH,Ngeder] %in% 1)
StatRes[1,1:2] <- c(paste(M.C,":",length(posC)-M.C,sep=""),paste(M.H,":",length(posH)-M.H,sep=""))
#--* proportion test for number of gender/male:female with Cirrhosis and HCC
sexhcc <-matrix(c(M.C,M.H,length(posC)-M.C,length(posH)-M.H),ncol=2,byrow=T)
rownames(sexhcc)<-c("male","female")
colnames(sexhcc)<-c("Cirrhosis","HCC")
StatRes[1,3] <- round(prop.test(sexhcc)$p.value,3)

#--* extract age information and test
Nage <- 6
StatRes[2,1] <-  paste(round(mean(dat[posC,Nage]),1)," ¡À ",round(sd(dat[posC,Nage]),1),sep="")
StatRes[2,2] <-  paste(round(mean(dat[posH,Nage]),1)," ¡À ",round(sd(dat[posH,Nage]),1),sep="")
StatRes[2,3] <-  round(t.test(dat[posC,Nage],dat[posH,Nage])$p.value,3)

#i=3
ni = 2
for(i in others){
  ni = ni + 1
  Wlc.test <- wilcox.test(dat[posC,i],dat[posH,i],alternative = c("two.sided"))
  #Wlc.test <- t.test(dat[posC,i],dat[posH,i],alternative = c("two.sided"))
  #ci(dat[1:24,i])
  #sd(dat[1:24,i])/24^(1/2)
  StatRes[ni,1] <- paste(round(mean(dat[posC,i]),3)," ¡À ",round(sd(dat[posC,i]),3),sep="")
  StatRes[ni,2] <- paste(round(mean(dat[posH,i]),3)," ¡À ",round(sd(dat[posH,i]),3),sep="")
  StatRes[ni,3] <- Wlc.test$p.value
}
write.table(StatRes,"BasicTest.txt",sep="\t",quote=F)



#####################################################
## read data  -----------------------2018/08/03
#####################################################
setwd("E:/work/jiangjian/20180803/HCvsHCC")

#dat <- read.csv("lung.All.log10AFP.YvsC.csv")
#dat <- read.csv("lung.lowAFP.YvsC.csv")

dat <- read.csv("data_HCvsHCC.csv")
head(dat)
dim(dat)
log10Fe <- data.frame(log10(dat[,6]))
colnames(log10Fe) <- "log10Fer"
dat[,7] <- log10Fe

posC <- which(dat[,2] == 0)
posH <- which(dat[,2] == 1)
others <- c(3,4,5,6,7)
StatRes <- matrix(NA,ncol(dat)-2,3)
rownames(StatRes) <- c("Gender/Male:Female","Age",colnames(dat)[others])
colnames(StatRes) <- c("Cirrhosis","HCC","p-value")

#--* extract gender information
Ngeder <- 9
M.C <- sum(dat[posC,Ngeder] %in% 1)
M.H <- sum(dat[posH,Ngeder] %in% 1)
StatRes[1,1:2] <- c(paste(M.C,":",length(posC)-M.C,sep=""),paste(M.H,":",length(posH)-M.H,sep=""))
#--* proportion test for number of gender/male:female with Cirrhosis and HCC
sexhcc <-matrix(c(M.C,M.H,length(posC)-M.C,length(posH)-M.H),ncol=2,byrow=T)
rownames(sexhcc)<-c("male","female")
colnames(sexhcc)<-c("Cirrhosis","HCC")
StatRes[1,3] <- round(prop.test(sexhcc)$p.value,3)

#--* extract age information and test
Nage <- 8
StatRes[2,1] <-  paste(round(mean(dat[posC,Nage]),1)," ¡À ",round(sd(dat[posC,Nage]),1),sep="")
StatRes[2,2] <-  paste(round(mean(dat[posH,Nage]),1)," ¡À ",round(sd(dat[posH,Nage]),1),sep="")
StatRes[2,3] <-  round(t.test(dat[posC,Nage],dat[posH,Nage])$p.value,3)

#i=3
ni = 2
for(i in others){
  ni = ni + 1
  Wlc.test <- wilcox.test(dat[posC,i],dat[posH,i],alternative = c("two.sided"))
  #Wlc.test <- t.test(dat[posC,i],dat[posH,i],alternative = c("two.sided"))
  #ci(dat[1:24,i])
  #sd(dat[1:24,i])/24^(1/2)
  StatRes[ni,1] <- paste(round(mean(dat[posC,i]),3)," ¡À ",round(sd(dat[posC,i]),3),sep="")
  StatRes[ni,2] <- paste(round(mean(dat[posH,i]),3)," ¡À ",round(sd(dat[posH,i]),3),sep="")
  StatRes[ni,3] <- Wlc.test$p.value
}
write.table(StatRes,"BasicTest.txt",sep="\t",quote=F)