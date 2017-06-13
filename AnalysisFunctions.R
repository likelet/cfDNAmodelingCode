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


# Ploting riscscore distribution across samples ----
# return a ggplot object 
# color have already been defined , one can change the color through modifying script directly
plotriskscoreInpatient=function(riskscorevec,status){
  df=data.frame(Riskscore=riskscorevec,Isalive=as.factor(status))
  df=df[with(df, order(-riskscorevec)), ]
  df=cbind(patientsnumber=seq(1:nrow(df)),df)
  df$patientsnumber=factor(df$patientsnumber,levels=seq(1:nrow(df)))
  p=ggplot(df,aes(x=patientsnumber,y=Riskscore,fill=as.factor(Isalive)))+geom_bar(stat="identity")+ggtitle(deparse(substitute(riskscorevec)))
  p=p+theme_base()+theme(axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank(),
                         legend.position=c(.9,.8))+
    scale_fill_manual(values=c("#56B4E9", "#E69F00"), 
                      name="")
  return(p)
}


# unicox analysis which always keep the safe factor as a reference when performing surival comparison
#Only categorize variables with two level were taken into calculations, otherwise ingnored.
printunicox2levelCategory=function(surv,df){
 datalist=list()
   for (i in 1:ncol(df)) {
    fml=as.formula(paste("surv ~ ",aquote(names(df)[i])))
    vec=unique(df[,i])
    vec=vec[!is.na(vec)]
    fit<-coxph(fml,df)
    tempsum=summary(fit)
    cof=tempsum$coefficients[1]
    # print(cof)
    if(length(vec)==2&&cof<0){
      vec=unique(df[,i])
      df[,i]=as.factor(df[,i])
      df[,i]= relevel(df[,i], ref =as.character(vec[1]))
      
      fit1<-suppressMessages(coxph(fml,data=df))
      sum1=summary(fit1)
      df[,i]= relevel(df[,i], ref =as.character(vec[2]))
      
      fit2<-suppressMessages(coxph(fml,data=df))
      sum2=summary(fit2)
      if(sum2$coefficients[1]>sum1$coefficients[1]){
        fit=fit2
      }else{
        fit=fit1
      }
      tempsum=summary(fit)
      cof=tempsum$coefficients[1]
     
    }
  
    y=data.frame(prettify(tempsum))
    colnames(y)=c("name","coef","Hazard Ratio","CI (lower)","CI (upper)","se(coef)","z","Pr(>|Z|)","na")
    datalist[[i]]=y
   }
 big_data = do.call(rbind, datalist)
 return(big_data)
}


#get survival information of a variable at a certain time ----
#status,binary vector
#df, data matrix
#cutpoint, surival time cut point when calculating survival rate
getsurvivalRatefun=function(surv,status,df,cutpoint){
  tsum=summary(survfit(surv~ status,data=df),times=cutpoint)
  tsumtable=cbind(tsum$table,tsum$time,tsum$n.risk,tsum$n.event,tsum$n.censor,tsum$surv,tsum$std.err,tsum$lower,tsum$upper)
  colnames(tsumtable)=c(colnames(tsum$table),names(tsum)[c(2,3,4,5,6,9,10,11)])
  rownames(tsumtable)=paste(deparse(substitute(status)),row.names(tsumtable),sep="_")
  return(tsumtable)
}

##filter ROC site wo avoid concave ROC plot
ROCfilter<-function(x){
  #sort data by column 1
  x <-x[order(x[,2]),]
  row.names(x)<-c(1:nrow(x))
  
  #set temp maxB
  maxB <- -1
  list=c()
  #remove decreased B
  for (i in 1:nrow(x)){
    if((x[i,2] > maxB)){
      maxB<-x[i,2]
      
      #print(maxB)
    }else {
      list<-c(list,i*-1)
      #    print(i)
    }
  }
  #in case that list equals null
  if(is.null(list)==FALSE){
    x<-unique(x[list,])
    row.names(x)<-c(1:nrow(x)) 
  }
  #add 0,0 and 1,1 for ROC
  x<-unique(rbind(c(0,0),x,c(1,1)))
  row.names(x)<-c(1:nrow(x))
  
  return(x)
}

#########multiple plot function from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#save ggplot object into pdf file with filenames adopted frome variable name---
savePlot <- function(myPlot) {
pdf(paste(deparse(substitute(myPlot)),".pdf",sep=""))
print(myPlot)
dev.off()
}
