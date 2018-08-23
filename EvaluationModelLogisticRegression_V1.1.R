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
#input: input data file
#mtd: resample method, only for bootstrap and fold_cv:k
#times: the number of resample
#fm: 




buildLogisticReg <- function(wd,input,mtd,times,fm,fmfile,CI,fpr.fixed){
  #E:\work\jiangjian\data\lung/rawData
  setwd(wd)
  #read data
  dat <- read.csv(input)
  ##------------
  ## internal method
  ##------------
  mtd <- unlist(strsplit(mtd,":"))
  if(mtd[1] == "bootstrap"){
    nM <- 2
  }else{
    nM <- 1
  }
  
  ##replace number
  nSam <- times
  ## formula for reading
  X <- read.table(fmfile,sep="\t",stringsAsFactors=F)
  
  if(fm == "factor"){
    TolCom <- NULL
    #X <- c("Response","log10AFP","KIN","ALT","ALP","Age","Gender")
    X1 <- X[1,1]
    X2 <- X[1,-1]
    
    #i=2
    for(i in 1:length(X2)){
      combn6 <- combn(1:ncol(X2), i)
      #j=1
      for(j in 1:ncol(combn6)){
        TolCom <- rbind(TolCom,paste(X1,"~",paste(as.character(X2[1,combn6[,j]]),collapse="+"),sep=""))
      }
    }
    Xformula <- TolCom
  }else if(fm == "formula"){
    X1 <- unlist(strsplit(X[1,],"~"))[1]
    Xformula <- as.matrix(X)
  }
  auc <- matrix(NA,nrow(Xformula),nM*nSam)
  rownames(auc) <- Xformula[,1]
  
  aic <- matrix(NA,nrow(Xformula),1)
  rownames(aic) <- Xformula[,1]
  
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
  #n=1
  for(n in 1:nSam){
    #set.seed(n)
    dat1 <- dat
    ##------------
    ## spliting samples
    ##-----------    
    if(mtd[1] == "bootstrap"){
      ind <- sample(1:dim(dat)[1],dim(dat)[1],replace = TRUE)
      boot_data <- dat[ind,]
      
      traindat <- boot_data
      ##AUC of test data replaced by full data
      testdat <- dat
      
      ##AUC of test data replaced by boot data
      testdat1 <- boot_data
    }else if(mtd[1] == "fold_cv"){
      test1 <- sort(sample(1:dim(dat)[1],round(dim(dat)[1]/as.numeric(mtd[2]),0)))
      test0 <- c(1:dim(dat)[1])[!1:dim(dat)[1] %in% test1]
      
      traindat <- dat1[test0,]
      testdat <- dat1[test1,]
    }
    
    if(nlevels(as.factor(testdat[,X1]))==2){
      #glm.fit <- glm(Respons~AFP+KIN+ALT+ALP+Age+Gender,data=traindat,family = binomial(link='logit'))
      ik <- 0
      #i=1
      for(i in 1:nrow(Xformula)){
        ik <- ik+1
        asfl <- as.formula(Xformula[i])
        #print(asfl)
        glm.lgrg <- glm(asfl,data=traindat,family = binomial(link='logit'))
        if(n == nSam){
          glm.fit.full <- glm(asfl,data=dat,family = binomial(link='logit'))
          smyglm <- summary(glm.fit.full)
          aic[ik,] <- smyglm$aic
        }
        
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
        optCutOff <- optimalCutoff(testdat[,X1], p)[1]
        cutoff_mtr[ik,n] <- optCutOff
        misCEr <- misClassError(testdat[,X1], prob, threshold = optCutOff)
        confMtr <- confusionMatrix(testdat[,X1], prob, threshold = optCutOff)
        sen <- sensitivity(testdat[,X1], prob, threshold = optCutOff)
        spe <- specificity(testdat[,X1], prob, threshold = optCutOff)
        sen_mtr[ik,n] <- sen
        spe_mtr[ik,n] <- spe
        
        
        if(mtd[1] == "bootstrap"){
          
          p <- predict(glm.lgrg, newdata=testdat1, type="response")
          pr <- prediction(p, testdat1[,X1])
          prf <- performance(pr, measure = "tpr", x.measure = "fpr")
          tpr <- unlist(prf@y.values)
          fpr <- unlist(prf@x.values)
          ##------------------------------
          ## TPR with a series of fiexd FPR
          ##------------------------------
          #ii=1
          for(ii in 1:length(fpr.fixed)){
            fpr.sub <- fpr.fixed[ii]
            if(fpr.sub %in% fpr){
              pos <- sum(fpr <= fpr.sub)
            }else{
              pos1 <- sum(fpr <= fpr.sub)
              pos <- c(pos1,pos1+1)
            }
            #tpr_series[ik+(ii-1)*length(fpr.fixed),n+nSam] <- mean(tpr[pos])
            tpr_series[ik+(ii-1)*nrow(Xformula),n+nSam] <- mean(tpr[pos])
          }
          #plot(prf)
          
          ##------------------------------
          ## AUC statistic
          ##------------------------------
          auc_list <- performance(pr, measure = "auc")
          auc_value <- auc_list@y.values[[1]]
          auc[ik,n+nSam] <- round(auc_value,3)
          
          ##------------------------------
          ## Specitivity and Sensetivity
          ##------------------------------
          prob = plogis(predict(glm.lgrg, testdat))
          optCutOff <- optimalCutoff(testdat1[,X1], p)[1]
          cutoff_mtr[ik,n+nSam] <- optCutOff
          misCEr <- misClassError(testdat1[,X1], prob, threshold = optCutOff)
          confMtr <- confusionMatrix(testdat1[,X1], prob, threshold = optCutOff)
          sen <- sensitivity(testdat1[,X1], prob, threshold = optCutOff)
          spe <- specificity(testdat1[,X1], prob, threshold = optCutOff)
          sen_mtr[ik,n+nSam] <- sen
          spe_mtr[ik,n+nSam] <- spe
        }
      }
    }
    if(n %% (nSam/10) == 0){
      cat(n,"\n")
    }
  }
  
  if(mtd[1] == "bootstrap"){
    colnames(auc) <- c(paste("FullAUC",1:nSam,sep=""),paste("BootAUC",1:nSam,sep=""))
    colnames(tpr_series) <- c(paste("FullTPR",1:nSam,sep=""),paste("BootTPR",1:nSam,sep=""))
    colnames(cutoff_mtr) <- c(paste("FullCutOff",1:nSam,sep=""),paste("BootCutOff",1:nSam,sep=""))
    colnames(sen_mtr) <- c(paste("FullSen",1:nSam,sep=""),paste("BootSen",1:nSam,sep=""))
    colnames(spe_mtr) <- c(paste("FullSpe",1:nSam,sep=""),paste("BootSpe",1:nSam,sep=""))
    
    write.csv(auc,file=paste("AUC.raw.",mtd[1],".",nSam,".csv",sep=""))
    write.csv(tpr_series,file=paste("TPR.raw.",mtd[1],".",nSam,".csv",sep=""))
    write.csv(cutoff_mtr,file=paste("CutOff.raw.",mtd[1],".",nSam,".csv",sep=""))
    write.csv(sen_mtr,file=paste("Sensitivity.raw.",mtd[1],".",nSam,".csv",sep=""))
    write.csv(spe_mtr,file=paste("Specitivity.raw.",mtd[1],".",nSam,".csv",sep=""))
    
    auc1 <- auc[,1:nSam]
    auc2 <- auc[,1:nSam+nSam]
    auc_CI <- matrix(NA,nrow(Xformula),8)
    colnames(auc_CI) <- c(paste("Full.",c("MeanAUC","CI lower","CI upper","Std. Error"),sep=""),paste("Boot.",c("MeanAUC","CI lower","CI upper","Std. Error"),sep=""))
    
    tpr_series1 <- tpr_series[,1:nSam]
    tpr_series2 <- tpr_series[,1:nSam+nSam]
    tpr_CI <- matrix(NA,nrow(Xformula)*length(fpr.fixed),8)
    colnames(tpr_CI) <- c(paste("Full.",c("MeanTPR","CI lower","CI upper","Std. Error"),sep=""),paste("Boot.",c("MeanTPR","CI lower","CI upper","Std. Error"),sep=""))
    rownames(tpr_CI) <- rownames(tpr_series)
    
    cutoff_mtr1 <- cutoff_mtr[,1:nSam]
    cutoff_mtr2 <- cutoff_mtr[,1:nSam+nSam]
    cutoff_mtr_CI <- matrix(NA,nrow(Xformula),8)
    colnames(cutoff_mtr_CI) <- c(paste("Full.",c("Meancutoff","CI lower","CI upper","Std. Error"),sep=""),paste("Boot.",c("Meancutoff","CI lower","CI upper","Std. Error"),sep=""))    
    rownames(cutoff_mtr_CI) <- Xformula
    
    sen_mtr1 <- sen_mtr[,1:nSam]
    sen_mtr2 <- sen_mtr[,1:nSam+nSam]
    sen_mtr_CI <- matrix(NA,nrow(Xformula),8)
    colnames(sen_mtr_CI) <- c(paste("Full.",c("Meansen","CI lower","CI upper","Std. Error"),sep=""),paste("Boot.",c("Meansen","CI lower","CI upper","Std. Error"),sep=""))
    rownames(sen_mtr_CI) <- Xformula
    
    spe_mtr1 <- spe_mtr[,1:nSam]
    spe_mtr2 <- spe_mtr[,1:nSam+nSam]
    spe_mtr_CI <- matrix(NA,nrow(Xformula),8)
    colnames(spe_mtr_CI) <- c(paste("Full.",c("Meanspe","CI lower","CI upper","Std. Error"),sep=""),paste("Boot.",c("Meanspe","CI lower","CI upper","Std. Error"),sep=""))
    rownames(spe_mtr_CI) <- Xformula
    
  }else if(mtd[1] == "fold_cv"){
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
  }
  
  #k=1
  for(k in 1:nrow(Xformula)){
    if(mtd[1] == "bootstrap"){
      auc_CI[k,1:4] <- ci(auc1[k,], confidence=CI,na.rm=T)
      auc_CI[k,5:8] <- ci(auc2[k,], confidence=CI,na.rm=T)
      for(w in 1:length(fpr.fixed)){
        #w=3
        #tpr_CI[k+(w-1)*nrow(Xformula),1:4] <- ci(tpr_series1[k+(w-1)*length(fpr.fixed),1:nSam], confidence=CI,na.rm=T)   ## motified by junhuili @ 20180717
        #tpr_CI[k+(w-1)*nrow(Xformula),5:8] <- ci(tpr_series2[k+(w-1)*length(fpr.fixed),1:nSam], confidence=CI,na.rm=T)   ## motified by junhuili @ 20180717
        tpr_CI[k+(w-1)*nrow(Xformula),1:4] <- ci(tpr_series1[k+(w-1)*nrow(Xformula),1:nSam], confidence=CI,na.rm=T)
        tpr_CI[k+(w-1)*nrow(Xformula),5:8] <- ci(tpr_series2[k+(w-1)*nrow(Xformula),1:nSam], confidence=CI,na.rm=T)
        
      }
      
      cutoff_mtr_CI[k,1:4] <- ci(cutoff_mtr1[k,], confidence=CI,na.rm=T)
      cutoff_mtr_CI[k,5:8] <- ci(cutoff_mtr2[k,], confidence=CI,na.rm=T)
      
      sen_mtr_CI[k,1:4] <- ci(sen_mtr1[k,], confidence=CI,na.rm=T)
      sen_mtr_CI[k,5:8] <- ci(sen_mtr2[k,], confidence=CI,na.rm=T)
      
      spe_mtr_CI[k,1:4] <- ci(spe_mtr1[k,], confidence=CI,na.rm=T)
      spe_mtr_CI[k,5:8] <- ci(spe_mtr2[k,], confidence=CI,na.rm=T)
      
    }else if(mtd[1] == "fold_cv"){
      auc_CI[k,] <- ci(auc[k,], confidence=CI,na.rm=T)
      for(w in 1:length(fpr.fixed)){
        #tpr_CI[k+(w-1)*nrow(Xformula),] <- ci(tpr_series[k+(w-1)*length(fpr.fixed),1:nSam], confidence=CI,na.rm=T)   ## motified by junhuili @ 20180717
        tpr_CI[k+(w-1)*nrow(Xformula),] <- ci(tpr_series[k+(w-1)*nrow(Xformula),1:nSam], confidence=CI,na.rm=T)
      }
      cutoff_mtr_CI[k,] <- ci(cutoff_mtr[k,], confidence=CI,na.rm=T)
      sen_mtr_CI[k,] <- ci(sen_mtr[k,], confidence=CI,na.rm=T)
      spe_mtr_CI[k,] <- ci(spe_mtr[k,], confidence=CI,na.rm=T)
      
    }
  }
  aicauc <- cbind(aic,auc_CI)
  rownames(aicauc) <- Xformula
  colnames(aicauc)[1] <- c("AIC")
  
  #write.csv(aicauc,file=paste("AUC.summary.",mtd,".Full.",nSam,".lung.All.YC.csv",sep=""))
  #write.csv(aicauc,file=paste("AUC.summary.",mtd,".Boot.",nSam,".lung.All.YC.csv",sep=""))
  if(mtd[1] == "fold_cv"){
    write.csv(aicauc,file=paste("AUC.summary.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
    write.csv(tpr_CI,file=paste("TPR.summary.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
    write.csv(cutoff_mtr_CI,file=paste("CutOff.summary.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
    write.csv(sen_mtr_CI,file=paste("Sen.summary.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
    write.csv(spe_mtr_CI,file=paste("Spe.summary.",mtd[2],mtd[1],".",nSam,".csv",sep=""))
    
  }else if(mtd[1] == "bootstrap"){
    write.csv(aicauc,file=paste("AUC.summary.",mtd[1],".",nSam,".csv",sep=""))
    write.csv(tpr_CI,file=paste("TPR.summary.",mtd[1],".",nSam,".csv",sep=""))
    write.csv(cutoff_mtr_CI,file=paste("CutOff.summary.",mtd[1],".",nSam,".csv",sep=""))
    write.csv(sen_mtr_CI,file=paste("Sen.summary.",mtd[1],".",nSam,".csv",sep=""))
    write.csv(spe_mtr_CI,file=paste("Spe.summary.",mtd[1],".",nSam,".csv",sep=""))
  }
}


##----------------------
##plot t.test 
##----------------------
##plot
ROCplot <- function(wd,input,fmfile,output){
  setwd(wd)
  dat <- read.csv(input)
  X <- read.table(fmfile,sep="\t",stringsAsFactors=F)
  if(nrow(X)==1){
    glm.AFP <- glm(as.formula(X[1,1]),data=dat,family = binomial(link='logit'))
    p.AFP <- predict(glm.AFP, newdata=dat, type="response")
    pr.AFP <- prediction(p.AFP,dat$Response)
    prf.AFP <- performance(pr.AFP, measure = "tpr", x.measure = "fpr")
    auc.AFP <- round(performance(pr.AFP, measure = "auc")@y.values[[1]],3)
    png(paste(wd,"/",output,".png",sep=""),width = 480*2, height = 480*2,units = "px", pointsize = 12*2)
    plot(prf.AFP,main=X[1,1])
    legend(0.5,0.2,c(paste(X[1,2]," & AUC:",auc.AFP,sep="")),col = 1,lty=1)
    dev.off()
  }else{
    lgd <- NULL
    png(paste(wd,"/",output,".png",sep=""),width = 480*2, height = 480*2,units = "px", pointsize = 12*2)
    for(i in 1:nrow(X)){
      glm.AFP <- glm(X[i,1],data=dat,family = binomial(link='logit'))
      p.AFP <- predict(glm.AFP, newdata=dat, type="response")
      pr.AFP <- prediction(p.AFP,dat$Response)
      prf.AFP <- performance(pr.AFP, measure = "tpr", x.measure = "fpr")
      auc.AFP <- round(performance(pr.AFP, measure = "auc")@y.values[[1]],3)
      
      if(i==1){
        plot(prf.AFP,main="ROC plot")
      }else{
        lines(prf.AFP@x.values[[1]], prf.AFP@y.values[[1]], col = i,lty=i)
      }
      lgd <- append(lgd,c(paste(X[i,2]," & AUC:",auc.AFP,sep="")))
    }
    abline(v = 0.05, col="red", lty=3)
    abline(v = 0.1, col="red", lty=3)
    abline(v = 0.15, col="red", lty=3)
    legend(0.4,0.2,lgd,col = c(1:nrow(X)),lty=c(1:nrow(X)))
    dev.off()
  }
}
##t.test and Wilcoxon Signed-Rank Test
SignTest <- function(wd,input,fmfile){
  setwd(wd)
  dat <- read.csv(input,stringsAsFactors=F)
  X <- read.table(fmfile,sep="\t",stringsAsFactors=F)
  
  if(nrow(X)<2){
    stop("No. of formula should be more than 1! Good luck!")
  }else{
    test.ind <- combn(1:nrow(X),2)
    test.Res <- combn(X[1:nrow(X),2],2)
    test.Res <- rbind(test.Res,matrix(NA,2,ncol(test.Res)))
    colnames(test.Res) <- paste("Test",1:ncol(test.Res),sep="")
    rownames(test.Res) <- c("formula1","formula2","t.test","wilcox.test")
    for(i in 1:ncol(test.Res)){
      Xdat <- as.numeric(dat[which(dat[,1] %in% X[test.ind[1,i],1]),-1])
      Ydat <- as.numeric(dat[which(dat[,1] %in% X[test.ind[2,i],1]),-1])
      
      test.Res[3,i] <- t.test(Xdat,Ydat,var.equal = T,paired = T)$p.value
      test.Res[4,i] <- wilcox.test(Xdat, Ydat, paired=TRUE)$p.value
    }
  }
  write.table(test.Res,"SignificantTest.csv",sep=",")
}


##buildLogisticReg
#wd = "E:/work/jiangjian/data/lung/rawData"
#input = "lung.All.log10AFP.YvsC.csv"
#mtd = "bootstrap"
#times = 20
#fm = "formula"  #fm = "factor"
#fmfile = "factor2.txt"
#CI =0.95
#fpr.fixed <- c(0.05,0.1,0.15)
#buildLogisticReg(wd,input,mtd,times,fm,fmfile,CI,fpr.fixed)

##ROCplot
#wd = "E:/work/jiangjian/data/lung/rawData"
#input = "lung.All.log10AFP.YvsC.csv"
#fmfile = "factor2.1.txt"
#output = "ROC"
#ROCplot(wd,input,fmfile,output)

##SignTest
#wd = "E:/work/jiangjian/data/lung/rawData"
#input = "AUC.raw.3fold_cv.1000.csv"
#fmfile = "factor2.1.txt"
#SignTest(wd,input,fmfile)