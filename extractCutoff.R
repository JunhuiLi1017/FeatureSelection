extractCutoff <- function(wd,train,fmfile,foldfile,bootfile,outfile){
  setwd(wd)
  a1 <- read.table(fmfile,sep="\t",stringsAsFactors = F)
  a2 <- read.table(foldfile,sep=",",header=T)
  a3 <- read.table(bootfile,sep=",",header=T)
  dat <- read.csv(train)
  
  a4 <- matrix(NA,dim(a1)[1],4)
  colnames(a4) <- c("Formula","Coef","fold_cv","bootstrap")
  a4[,1] <- as.character(a1[,1])
  a4[,3] <- a2[as.character(a2[,1]) %in% as.character(a1[,1]),2]
  a4[,4] <- a3[as.character(a3[,1]) %in% as.character(a1[,1]),2]
  
  FormValue <- NULL
  Zvalue <- NULL
  #i=2
  for(i in 1:nrow(a1)){
    asfl <- as.formula(a1[i,1])
    glm.lgrg <- glm(asfl,data=dat,family = binomial(link='logit'))
    if(i == 1){
      Zvalue <- summary(glm.lgrg)$coefficients
    }else{
      Zvalue <- rbind(Zvalue,summary(glm.lgrg)$coefficients)
    }
    wei <- glm.lgrg$coefficients
    subFormValue <- append(a1[i,1],paste(wei,collapse = "+"))
    FormValue <- rbind(FormValue,subFormValue)
  }
  a4[,2] <- FormValue[,2]
  write.table(Zvalue,file="coef_pvalue.txt",sep="\t",quote=F)
  write.table(a4,file=outfile,sep="\t",quote=F,row.names = F)
}

##extractCutoff
#source(extractCutoff.R)
#wd = "E:/work/jiangjian/20180803/newExt"
#train = "data_CvsHCC.csv"
#fmfile = "factor2.2.txt"
#foldfile = "CutOff.summary.3fold_cv.1000.csv"
#bootfile = "CutOff.summary.bootstrap.1000.csv"
#outfile = "factor2.3.txt"
#extractCutoff(wd,train,fmfile,foldfile,bootfile,outfile)



library(magrittr)
#install.packages("InformationValue")
library(InformationValue)
#install.packages("gmodels")
library(gmodels)
library(ROCR)

setCoeffs <- function(frml, weights, len){
  el <- paste0("offset(", weights[-1], "*", 
               unlist(strsplit(as.character(frml)[-(1:2)], " +\\+ +")), ")")
  el <- c(paste0("offset(rep(", weights[1], ",", len, "))"), el)                                 
  as.formula(paste(as.character(frml)[2], "~", 
                   paste(el, collapse = " + "), " + -1"))
}
ROC_ext <- function(wd,input,fmfile,output){
  setwd(wd)
  dat <- read.csv(input)
  X <- read.table(fmfile,sep="\t",stringsAsFactors=F,header=T)
    lgd <- NULL
    png(paste(wd,"/",output,".png",sep=""),width = 480*2, height = 480*2,units = "px", pointsize = 12*2)
    for(i in 1:nrow(X)){
      mod1 <- glm(X[i,1],data=dat,family = binomial(link='logit'))
      givenWei <- as.numeric(unlist(strsplit(X[i,2], "[+]", perl = T)))
      mod2 <- update(mod1, setCoeffs(as.formula(X[i,1]), givenWei, nrow(dat)),family = binomial(link='logit'))
      #p.AFP <- predict(glm.AFP, newdata=dat, type="response")
      p.AFP <- predict(mod2,dat,type="response")
      pr.AFP <- prediction(p.AFP,dat$Response)
      prf.AFP <- performance(pr.AFP, measure = "tpr", x.measure = "fpr")
      auc.AFP <- round(performance(pr.AFP, measure = "auc")@y.values[[1]],3)
      
      if(i==1){
        plot(prf.AFP,main="ROC plot")
      }else{
        lines(prf.AFP@x.values[[1]], prf.AFP@y.values[[1]], col = i,lty=i)
      }
      lgd <- append(lgd,c(paste(unlist(strsplit(X[i,1],"~"))[2]," & AUC:",auc.AFP,sep="")))
    }
    abline(v = 0.05, col="red", lty=3)
    abline(v = 0.1, col="red", lty=3)
    abline(v = 0.15, col="red", lty=3)
    legend(0.6,0.15,lgd,col = c(1:nrow(X)),lty=c(1:nrow(X)),cex=0.5)
    dev.off()
}

#wd = "E:/work/jiangjian/20180803/newExt"
#input = "data_CvsHCC.lowAFP.csv"
#fmfile = "factor2.3.txt"
#output="ROC1"
#ROC_ext(wd,input,fmfile,output)





