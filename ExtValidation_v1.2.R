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
## Parameter description
#####################################################
#wd: work directory
#train:Train data
#test:External data
#mtd: resample method, only for bootstrap and fold-cv:k
#times: the number of resample
#fmfile: formula files
#CI:confidence interval

ExtValidation <- function(wd,train,test,mtd,times,fmfile,CI,fpr.fixed){
  #E:/work/jiangjian/data/lung/rawData
  setwd(wd)
  #read data
  dat_train <- read.csv(train)
  dat_test <- read.csv(test) #dat_test <- read.csv(train) motified by junhuili @ 20180803
  ##------------
  ## internal method
  ##------------
  mtd <- unlist(strsplit(mtd,":"))
  #if(mtd[1] == "bootstrap"){
  #  nM <- 2
  #}else{
  nM <- 1
  #}
  
  ##replace number
  nSam <- times
  ## formula for reading
  X <- read.table(fmfile,sep="\t",stringsAsFactors=F,header=T)
  X1 <- unlist(strsplit(X[1,1],"~"))[1]
  Xformula <- as.matrix(X[,1])
  
  auc <- matrix(NA,nrow(Xformula),nM*nSam)
  rownames(auc) <- Xformula[,1]
  
  #aic <- matrix(NA,nrow(Xformula),1)
  #rownames(aic) <- Xformula[,1]
  
  cutoff_mtr <- matrix(NA,nrow(Xformula),nM*nSam)
  rownames(cutoff_mtr) <- Xformula[,1]
  
  sen_mtr <- matrix(NA,nrow(Xformula),nM*nSam)
  rownames(sen_mtr) <- Xformula[,1]
  
  spe_mtr <- matrix(NA,nrow(Xformula),nM*nSam)
  rownames(spe_mtr) <- Xformula[,1]
  
  tpr_series <- matrix(NA,nrow(Xformula)*length(fpr.fixed),nM*nSam)
  tprrn <- NULL
  for(w in 1:length(fpr.fixed)){
    tprrn <- append(tprrn,paste(Xformula[,1],":",fpr.fixed[w],sep=""))
  }
  rownames(tpr_series) <- tprrn
  
  ik <- 0
  #i=1
  for(i in 1:nrow(Xformula)){
    ik <- ik+1
    asfl <- as.formula(Xformula[i])
    #print(asfl)
    glm.lgrg <- glm(asfl,data=dat_train,family = binomial(link='logit'))
    
    #n=1
    for(n in 1:nSam){
      #set.seed(n)
      dat1 <- dat_test
      ##------------
      ## spliting samples
      ##-----------    
      if(mtd[1] == "bootstrap"){
        ind <- sample(1:dim(dat1)[1],dim(dat1)[1],replace = TRUE)
        boot_data <- dat1[ind,]
        testdat <- boot_data
        
      }else if(mtd[1] == "fold_cv"){
        test1 <- sort(sample(1:dim(dat1)[1],round(dim(dat1)[1]/as.numeric(mtd[2]),0)))
        test0 <- c(1:dim(dat1)[1])[!1:dim(dat1)[1] %in% test1]
        
        traindat <- dat1[test0,]
        testdat <- dat1[test1,]
        
      }
      
      if(nlevels(as.factor(testdat[,X1]))==2){
        #glm.fit <- glm(Respons~AFP+KIN+ALT+ALP+Age+Gender,data=traindat,family = binomial(link='logit'))
        p <- predict(glm.lgrg, newdata=testdat, type="response")
        pr <- prediction(p, testdat[,X1])
        prf <- performance(pr, measure = "tpr", x.measure = "fpr")
        
        ##------------------------------
        ## TPR with a series of fiexd FPR
        ##------------------------------
        tpr <- unlist(prf@y.values)
        fpr <- unlist(prf@x.values)
        #ii=1
        for(ii in 1:length(fpr.fixed)){
          fpr.sub <- fpr.fixed[ii]
          if(fpr.sub %in% fpr){
            pos <- sum(fpr <= fpr.sub)
          }else{
            pos1 <- sum(fpr <= fpr.sub)
            pos <- c(pos1,pos1+1)
          }
          #tpr_series[ik+(ii-1)*length(fpr.fixed),n] <- mean(tpr[pos])
          tpr_series[ik+(ii-1)*nrow(Xformula),n] <- mean(tpr[pos])
        }
        
        #plot(prf)
        ##------------------------------
        ## AUC statistic
        ##------------------------------
        auc_list <- performance(pr, measure = "auc")
        auc_value <- auc_list@y.values[[1]]
        auc[ik,n] <- round(auc_value,3)
        
        ##------------------------------
        ## Specitivity and Sensetivity
        ##------------------------------
        prob = plogis(predict(glm.lgrg, testdat))
        #optCutOff <- optimalCutoff(testdat[,X1], p)[1]
        optCutOff <- X[i,mtd[1]]
        cutoff_mtr[ik,n] <- optCutOff
        misCEr <- misClassError(testdat[,X1], prob, threshold = optCutOff)
        confMtr <- confusionMatrix(testdat[,X1], prob, threshold = optCutOff)
        sen <- sensitivity(testdat[,X1], prob, threshold = optCutOff)
        spe <- specificity(testdat[,X1], prob, threshold = optCutOff)
        sen_mtr[ik,n] <- sen
        spe_mtr[ik,n] <- spe
      }
      if(n %% (nSam/10) == 0){
        cat(n,"\n")
      }
    }
    cat(Xformula[i],"\n")
  }
  
  #setwd("E:/work/jiangjian/data/lung/rawData/test1")
  colnames(auc) <- paste("AUC",1:nSam,sep="")
  colnames(tpr_series) <- paste("TPR",1:nSam,sep="")
  colnames(cutoff_mtr) <- paste("CutOff",1:nSam,sep="")
  colnames(sen_mtr) <- paste("Sen",1:nSam,sep="")
  colnames(spe_mtr) <- paste("Spe",1:nSam,sep="")
  
  write.csv(auc,file=paste("AUC.raw.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
  write.csv(tpr_series,file=paste("TPR.raw.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
  write.csv(cutoff_mtr,file=paste("CutOff.raw.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
  write.csv(sen_mtr,file=paste("Sensitivity.raw.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
  write.csv(spe_mtr,file=paste("Specitivity.raw.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
  
  auc_CI <- matrix(NA,nrow(Xformula),4)
  colnames(auc_CI) <- c("MeanAUC","CI lower","CI upper","Std. Error")
  rownames(auc_CI) <- Xformula[,1]

  tpr_CI <- matrix(NA,nrow(Xformula)*length(fpr.fixed),4)
  colnames(tpr_CI) <- c("MeanTPR","CI lower","CI upper","Std. Error")
  rownames(tpr_CI) <- rownames(tpr_series)
  
  cutoff_mtr_CI <- matrix(NA,nrow(Xformula),4)
  colnames(cutoff_mtr_CI) <- c("MeanCutOff","CI lower","CI upper","Std. Error")
  rownames(cutoff_mtr_CI) <- Xformula
  
  sen_mtr_CI <- matrix(NA,nrow(Xformula),4)
  colnames(sen_mtr_CI) <- c("MeanSen","CI lower","CI upper","Std. Error")
  rownames(sen_mtr_CI) <- Xformula
  
  spe_mtr_CI <- matrix(NA,nrow(Xformula),4)
  colnames(spe_mtr_CI) <- c("MeanSpe","CI lower","CI upper","Std. Error")
  rownames(spe_mtr_CI) <- Xformula
  
  #k=1
  for(k in 1:nrow(Xformula)){
    auc_CI[k,] <- ci(auc[k,], confidence=CI,na.rm=T)
    for(w in 1:length(fpr.fixed)){
      #tpr_CI[k+(w-1)*nrow(Xformula),] <- ci(tpr_series[k+(w-1)*length(fpr.fixed),1:nSam], confidence=CI,na.rm=T) ## motified by junhuili @ 20180717
      tpr_CI[k+(w-1)*nrow(Xformula),] <- ci(tpr_series[k+(w-1)*nrow(Xformula),1:nSam], confidence=CI,na.rm=T)
    }
    cutoff_mtr_CI[k,] <- ci(cutoff_mtr[k,], confidence=CI,na.rm=T)
    sen_mtr_CI[k,] <- ci(sen_mtr[k,], confidence=CI,na.rm=T)
    spe_mtr_CI[k,] <- ci(spe_mtr[k,], confidence=CI,na.rm=T)
  }
  #aicauc <- cbind(aic,auc_CI)
  #rownames(aicauc) <- Xformula
  #colnames(aicauc)[1] <- c("AIC")
  write.csv(auc_CI,file=paste("AUC.summary.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
  write.csv(tpr_CI,file=paste("TPR.summary.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
  #write.csv(cutoff_mtr_CI,file=paste("CutOff.summary.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
  write.csv(sen_mtr_CI,file=paste("Sen.summary.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
  write.csv(spe_mtr_CI,file=paste("Spe.summary.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
}

extractCutoff <- function(wd,fmfile,foldfile,bootfile,outfile){
  setwd(wd)
  a1 <- read.table(fmfile,sep="\t")
  a2 <- read.table(foldfile,sep=",",header=T)
  a3 <- read.table(bootfile,sep=",",header=T)
  
  a4 <- matrix(NA,dim(a1)[1],3)
  colnames(a4) <- c("Formula","fold_cv","bootstrap")
  a4[,1] <- as.character(a1[,1])
  a4[,2] <- a2[as.character(a2[,1]) %in% as.character(a1[,1]),2]
  a4[,3] <- a3[as.character(a3[,1]) %in% as.character(a1[,1]),2]
  
  write.table(a4,file=outfile,sep="\t",quote=F,row.names = F)
}

##External validation
#sorce("file/to/your/directory/ExtValidation.R")
#wd = "E:/work/jiangjian/data/lung/rawData"
#train = "lung.All.log10AFP.YvsC.csv"
#test = "lung.All.log10AFP.YvsC.csv"
#mtd = "bootstrap"
#times = 20
#fmfile = "factor2.txt"
#CI =0.95
#fpr.fixed <- c(0.05,0.1,0.15)
#ExtValidation(wd,train,test,mtd,times,fmfile,CI,fpr.fixed)


##extractCutoff
#source(extractCutoff.R)
#wd = "E:/work/jiangjian/20180713/Cirrhosis.VS.HCC"
#fmfile = "factor2.2.txt"
#foldfile = "CutOff.summary.3fold-cv.1000.csv"
#bootfile = "CutOff.summary.bootstrap.1000.csv"
#outfile = "factor2.2.txt"
