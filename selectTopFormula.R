source("E:/work/jiangjian/data/lung/rawData/EvaluationModelLogisticRegression.R")

wd = "E:/work/jiangjian/data/lung/rawData"
input = "lung.lowAFP.YvsC.csv"
#input = "lung.All.log10AFP.YvsC.csv"
#mtd = "bootstrap"
mtd = "fold-cv:3"
times = 1000
fm = "formula"
#fm = "factor"
fmfile = "factor2.2.txt"
CI =0.95
fpr.fixed <- c(0.05,0.1,0.15)
#buildLogisticReg(wd,input,mtd,times,fm,fmfile,CI,fpr.fixed)

##ROCplot
wd = "E:/work/jiangjian/data/lung/rawData"
#input = "lung.All.log10AFP.YvsC.csv"
input = "lung.lowAFP.YvsC.csv"
fmfile = "factor2.1.txt"
output = "ROC_lowAFP"
#ROCplot(wd,input,fmfile,output)

##SignTest
wd = "E:/work/jiangjian/data/lung/rawData"
input = "AUC.raw.3fold-cv.1000.csv"
fmfile = "factor2.1.txt"
#SignTest(wd,input,fmfile)







#aicauc <- read.csv("AUC.summary.bootstrap.1000.csv")
aicauc <- read.csv("AUC.summary.3fold-cv.1000.csv")
aicauc[1:5,1:5]
rownames(aicauc) <- aicauc[,1]
aicauc <- aicauc[,-1]
alphaAIC <- 0.1
alphaAUC <- 0.9
comX <- aicauc[aicauc[,1] < sort(aicauc[,1])[round(63*alphaAIC,0)] & aicauc[,2] > sort(aicauc[,2])[round(63*alphaAUC,0)],]
sum(aicauc[,1] < sort(aicauc[,1])[round(63*alphaAIC,0)] & aicauc[,2] > sort(aicauc[,2])[round(63*alphaAUC,0)])
a <- rownames(comX)


?t.test






setwd("E:/work/hanyuqing")
dat1 <- read.csv("ROC.csv")

glm.AFP <- glm(Response~X,data=dat1,family = binomial(link='logit'))
p.AFP <- predict(glm.AFP, newdata=dat1, type="response")
pr.AFP <- prediction(p.AFP,dat1$Response)
prf.AFP <- performance(pr.AFP, measure = "tpr", x.measure = "fpr")
auc.AFP <- round(performance(pr.AFP, measure = "auc")@y.values[[1]],3)

plot(prf.AFP,main="ROC")









