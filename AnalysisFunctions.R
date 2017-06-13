#used librarys for calculation or plotting

library("glmnet")
library("survival")
library("ggplot2")
library("survminer")
library("hdnom")
library("papeR")
library("risksetROC")
library("ggsci")
library("pROC")
library("grid")
library("ggthemes")

#preloading function------
library(papeR)

#Batch unicox regression analaysis with data matrix and a surival object from "surival" packages
#return a data.frame with several statistic figures 
batchUnivarCOXfun<-function(surv,data){
  df<-data.frame(name=c("coef","Hazard Ratio","CI (lower)","CI (upper)","se(coef)","z","Pr(>|Z|)","C-index","C-index-se","testCph.p"))
  pro=ncol(data)/10000
  j=1;
  k=1;
  for (i in 1:ncol(data)) {
    fml=as.formula(paste("surv ~ ",names(data)[i]))
    fit<-coxph(surv ~ data[,i])
    p<-data.frame(summary(fit)$coefficients)[1,5]
    test=cox.zph(coxph(surv ~ data[,i]))
    testp= test$table[3]
    #     x<-c(x,y)
    #     coefficients=summary(a)$coefficients
    y=prettify(summary(fit))
    y=y[,c(-1,-8,-9)]
    x=c(t(y)[,1],p,summary(fit)$concordance,testp)
    df=data.frame(df,x)
    if(i==pro*k){
      print(paste(j*k,"items finished"))
      k=k+1
    }
  }
  row.names(df)<-df[,1]
  df<-t(df[,-1])
  row.names(df)=names(data)
  return(df)
}

##Batch Wilcox test  analaysis with data matrix and category vector "class"--------
#return a data.frame with feature names and its match statistical p vaule based on the class 

batchWilcox=function(df,class,rownumber=2){
  
  singleWilcox=function(value,class){
    classV=unique(class)
    x=value[class==classV[1]]
    y=value[class==classV[2]]
    
    return(wilcox.test(x,y)$p.value)
  }
  
  valuelist=apply(df,rownumber,singleWilcox,class=class)
  nm=c()
  if(rownumber==2){
    nm=colnames(df)
  }else{
    nm=row.names(df)
  }
  return(data.frame(ID=nm,Pvalue=valuelist))
}


#cut continuous variable into binary variable with median value----

#1. return a factor of single vector
cutbymedian=function(vec){
  md=median(as.numeric(vec))
  vec2=rep("High",length(vec))
  vec2[vec<md]="Low"
  return(as.factor(vec2))
}

#2. dealing with a dataframe
cutbymedian_df=function(df){
  df2=data.frame(row.names(df))
  for(i in 1:ncol(df)){
    df2= cbind(df2,cutbymedian(df[,i]))
  }
  df2=df2[,-1]
  colnames(df2)=colnames(df)
  return(df2)
  
}
#3. return a binary output of 1/0
cugbymedian=function(vec){
  tempx=rep(1,length(vec))
  tempx[vec<median(vec)]=0
  vec=tempx
  return(vec)
  
}


#save ggplot object into pdf file with filenames adopted frome variable name---
savePlot <- function(myPlot) {
pdf(paste(deparse(substitute(myPlot)),".pdf",sep=""))
print(myPlot)
dev.off()
}
