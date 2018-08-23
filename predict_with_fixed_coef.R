setwd("E:/work/jiangjian/20180803/newExt")

FormValue <- NULL
for(i in 1:nrow(Xformula)){
  asfl <- as.formula(Xformula[i])
  glm.lgrg <- glm(asfl,data=dat_train,family = binomial(link='logit'))
  wei <- glm.lgrg$coefficients
  subFormValue <- append(Xformula[i],paste(wei,collapse = "+"))
  FormValue <- rbind(subFormValue, FormValue)
}
#write.table(FormValue,"coeff.txt",sep="\t",col.names=F,row.names = F,quote=F)

a <- read.table("coeff.txt",sep="\t",stringsAsFactors = F)
#n=1
for(n in 1:nrow(a)){
  fm1 <- as.formula(a[n,1])
  givenWei <- as.numeric(unlist(strsplit(a[n,2], "[+]", perl = T)))
  #mod11 <- glm(fm1,data=dat_train,family = binomial(link='logit'))
  mod1 <- glm(fm1,data=testdat,family = binomial(link='logit'))
  weights1 <- mod1$coefficients
  weights1 <- givenWei
  #for(i in 1:length(weights1)){
  #  weights1[i] <- givenWei[i]
  #}
  mod2 <- update(mod1, setCoeffs(fm1, weights1, nrow(testdat)),family = binomial(link='logit'))
  predict(mod2,testdat,type="response")
  #predict(mod11,testdat,type="response")
}

setCoeffs <- function(frml, weights, len){
  el <- paste0("offset(", weights[-1], "*", 
               unlist(strsplit(as.character(frml)[-(1:2)], " +\\+ +")), ")")
  el <- c(paste0("offset(rep(", weights[1], ",", len, "))"), el)                                 
  as.formula(paste(as.character(frml)[2], "~", 
                   paste(el, collapse = " + "), " + -1"))
}

