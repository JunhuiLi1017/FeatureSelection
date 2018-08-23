##---------------------------------------------------
## 2018.08.03 tutorial for Basic Test and evaluate model and external validation
##---------------------------------------------------

##**************
## Basic Test
##**************


##**************
## EvaluationModelLogisticRegression
##**************

#setwd("E:/work/jiangjian/20180713/Cirrhosis.VS.HCC")
source("E:/work/jiangjian/20180713/CorrectPrograme/EvaluationModelLogisticRegression_V1.1.R")
wd = "E:/work/jiangjian/20180803/CvsHCC"
input = "data_CvsHCC.csv"
#mtd = "fold_cv:3"
mtd = "bootstrap"
times = 1000
fm = "factor" #fm = "formula"  
fmfile = "factor.txt"
CI =0.95
fpr.fixed <- c(0.05,0.1,0.15)
buildLogisticReg(wd,input,mtd,times,fm,fmfile,CI,fpr.fixed)


##*****************
##select top important formula
##*****************
setwd(wd)
#aicauc <- read.csv("AUC.summary.3fold_cv.1000.csv")
aicauc <- read.csv("AUC.summary.bootstrap.1000.csv")
aicauc[1:5,1:5]
rownames(aicauc) <- aicauc[,1]
aicauc <- aicauc[,-1]
alphaAIC <- 0.2
alphaAUC <- 0.75

w1 <- 0
w <- c(1:5)
for(i in w){
  w1 <- w1 + ncol(combn(w,i))
}

comX <- aicauc[aicauc[,1] < sort(aicauc[,1])[round(w1*alphaAIC,0)] & aicauc[,2] > sort(aicauc[,2])[round(w1*alphaAUC,0)],]
sum(aicauc[,1] < sort(aicauc[,1])[round(w1*alphaAIC,0)] & aicauc[,2] > sort(aicauc[,2])[round(w1*alphaAUC,0)])
a <- rownames(comX)
#a.boot <- a
a.fold <- a


##ROCplot
input = "data_CvsHCC.csv"
fmfile = "factor2.1.txt"
output = "ROC"
ROCplot(wd,input,fmfile,output)

##**************
## External validation
##**************

source("E:/work/jiangjian/20180713/CorrectPrograme/ExtValidation_v1.2.R")
#ExtValidation_v1.1.R error in test <- read.csv(train)
getwd()
dat <- read.csv("data_CvsHCC.csv")
dat.lowAFP <- dat[dat[,3] < 1000,]
#write.csv(dat.lowAFP,"data_CvsHCC.lowAFP.csv",row.names=F)
getwd()
fmfile = "factor2.2.txt"
foldfile = "CutOff.summary.3fold_cv.1000.csv"
bootfile = "CutOff.summary.bootstrap.1000.csv"
outfile = "factor2.3.txt"
extractCutoff(wd,fmfile,foldfile,bootfile,outfile)

wd = "E:/work/jiangjian/20180803/CvsHCC/Ext"
train = "data_CvsHCC.csv"
test = "data_CvsHCC.lowAFP.csv"
#mtd = "fold_cv:3"
mtd = "bootstrap"
times = 1000
fmfile = "factor2.3.txt"
CI =0.95
fpr.fixed <- c(0.05,0.1,0.15)
ExtValidation(wd,train,test,mtd,times,fmfile,CI,fpr.fixed)

##ROCplot
input ="data_CvsHCC.lowAFP.csv"
fmfile = "factor2.1.txt"
output = "ROC-lowAFP"
ROCplot(wd,input,fmfile,output)