#dependencies 
#data processing functions
suppressMessages(library("sva"))
suppressMessages(library("impute"))
#lasso functions
suppressMessages(library("glmnet"))
#surival analysis
suppressMessages(library("survival"))
#plot functions
suppressMessages(library("ggplot2"))
suppressMessages(library("survminer"))
suppressMessages(library("hdnom"))
suppressMessages(library("papeR"))
suppressMessages(library("ggsci"))
suppressMessages(library("grid"))
suppressMessages(library("ggthemes"))
suppressMessages(library("cowplot"))
suppressMessages(library("papeR"))
suppressMessages(library("pheatmap"))
#perfunmance evaluation functions
suppressMessages(library("risksetROC"))
suppressMessages(library("pROC"))
suppressMessages(library("ROCR"))

#functions 
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
batchTtest <- function(df,class,rownumber=2){
  
  singleTtest=function(value,class){
    classV=unique(class)
    x=value[class==classV[1]]
    y=value[class==classV[2]]
    
    return(wilcox.test(x,y)$p.value)
  }
  singlefoldchange=function(value,class){
    
    classV=unique(class)
    x=value[class==classV[1]]
    y=value[class==classV[2]]
    
    return(mean(x)-mean(y))
  }
  valuelist=apply(df,rownumber,singleTtest,class=class)
  fclist=apply(df,rownumber,singlefoldchange,class=class)
  fdrlist=p.adjust(valuelist)
  
  nm=c()
  if(rownumber==2){
    nm=colnames(df)
  }else{
    nm=row.names(df)
  }
  return(data.frame(ID=nm,fc=fclist,Pvalue=valuelist,fdr=fdrlist))
}
plotriskscoreInpatient=function(riskscorevec,status){
  df=data.frame(Riskscore=riskscorevec,Isalive=as.factor(status))
  df=df[with(df, order(-riskscorevec)), ]
  df=cbind(patientsnumber=seq(1:nrow(df)),df)
  df$patientsnumber=factor(df$patientsnumber,levels=seq(1:nrow(df)))
  p=ggplot(df,aes(x=patientsnumber,y=Riskscore,fill=as.factor(Isalive)))+geom_bar(stat="identity")+ggtitle(deparse(substitute(riskscorevec)))
  p=p+theme_classic()+theme(axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(),
                            legend.position=c(.9,.8))+
    scale_fill_manual(values=c("#56B4E9", "#E69F00"), 
                      name="")
  return(p)
}
#cut vector by median
cutbymedian=function(vec){
  md=median(as.numeric(vec))
  vec2=rep("High",length(vec))
  vec2[vec<md]="Low"
  return(as.factor(vec2))
}
cutbymedian_df=function(df){
  df2=data.frame(row.names(df))
  for(i in 1:ncol(df)){
    df2= cbind(df2,cutbymedian(df[,i]))
  }
  df2=df2[,-1]
  colnames(df2)=colnames(df)
  return(df2)
  
}
cugbymedian=function(vec){
  tempx=rep(1,length(vec))
  tempx[vec<median(vec)]=0
  vec=tempx
  return(vec)
  
}
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
batchfoldchange=function(df,class,rownumber=2){
  singleWilcox=function(value,class){
    classV=unique(class)
    x=value[class==classV[1]]
    y=value[class==classV[2]]
    
    return(mean(y)/mean(x))
  }
  valuelist=apply(df,rownumber,singleWilcox,class=class)
  print(paste(unique(class)[2],"/",unique(class)[1]))
  nm=c()
  if(rownumber==2){
    nm=colnames(df)
  }else{
    nm=row.names(df)
  }
  return(data.frame(ID=nm,foldchange=valuelist))
}

#plot function
savePlot <- function(myPlot) {
  pdf(paste(deparse(substitute(myPlot)),".pdf",sep=""))
  print(myPlot)
  dev.off()
}
#write out table ----
zhaoqi.write.csv.survival <- function(analysisDF){
  analysisDF$HRci=paste(analysisDF$`Hazard Ratio`," (",analysisDF$`CI (lower)`,"-",analysisDF$`CI (upper)`,")",sep="")
  df<-analysisDF[,c(1,2,10,7,8)]
  
  return(df)
}
#confusing matrix 
#preloading functions 
#from https://stackoverflow.com/questions/24716337/test-quality-of-logistic-regression-model-using-confusion-matrix-and-roc-curve
confusion_matrix <- function(score,status, cutoff = median(score), plot.it = TRUE,title = NULL) {
  
  
  
  cdscore <- ifelse(score<= cutoff, "Low", "High")
  confusion <- table(newscore,status)
  if (plot.it) fourfoldplot(confusion, color = c("#CC6666", "#99CC99"),
                            conf.level = 0, margin = 1, main = title)
  confusion
  
}

#time-dependent ROC
#get time dependent AUC dataframe of certain marker 
getTimeDentAUCdf_multi<-function(time,status,marker,marekrstr=" ",timeline=360){
  #1 remove NA 
  df=data.frame(time,status,marker)
  df<-df[complete.cases(df),]
  a=c()
  for(i in 3:ncol(df)){
    a=c(a,paste("df[,",i,"]"))
  }
  fm<- as.formula(paste("Surv(df[,1], df[,2]) ~ ", paste(a, collapse= "+")))
  fit=coxph(fm)
  
  res=risksetROC(Stime=df[,1], status=df[,2],
                 marker=fit$linear.predictors, predict.time=timeline, method="Cox",plot = FALSE) 
  outdf<-data.frame(TP=res$TP,FP=res$FP,Maker=rep(paste(marekrstr," AUC:",round(res$AUC,4)),length(res$FP)))
  return(outdf)
}
getTimeDentAUCdf<-function(time,status,marker,marekrstr=" ",timeline=360){
  #1 remove NA 
  df=data.frame(time,status,marker)
  df<-df[complete.cases(df),]
  
  fit=coxph(Surv(df[,1], df[,2]) ~ df[,3])
  res=risksetROC(Stime=df[,1], status=df[,2],
                 marker=fit$linear.predictors, predict.time=timeline, method="Cox",plot = FALSE) 
  outdf<-data.frame(TP=res$TP,FP=res$FP,Maker=rep(paste(marekrstr," AUC:",round(res$AUC,4)),length(res$FP)))
  return(outdf)
}
#plot time dependent AUC 
plot_time_depent_surival_AUC<-function(df,cuttime,datatype){
  p=ggplot(data=df,aes(x=FP,y=TP,color=Maker))+geom_line(size=1)+
    labs(title=paste("ROC for ",cuttime/30," months survival predicting in ",datatype,"dataset",sep=" "),x="False positive rate", y = "True negative rate")+
    scale_color_npg()+theme_base()+theme(legend.position="top")+geom_abline(intercept = 0,linetype="dotted")+
    scale_y_continuous(limits=c(0,1),expand = c(0,0),breaks=c(0.2,0.4,0.6,0.8,1))+ 
    scale_x_continuous(limits=c(0,1),expand = c(0,0),breaks=c(0,0.2,0.4,0.6,0.8,1))+
    theme(
      legend.position = c(0.6, 0.25), # c(0,0) bottom left, c(1,1) top-right
      legend.background = element_rect(colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 10, face = 'bold')
    )
  return(p)
}


#plot ROC with multiple predictors
plotmulti_ROC_with_test<- function(predictor.matrix,status,predictor.name="predictor"){
  #loading packages 
  suppressMessages(library("ROCR"))# get colors 
  suppressMessages(library("ggsci"))# get colors
  #build data list
  pred.number<-ncol(predictor.matrix)
  col.list<-pal_npg("nrc", alpha = 1)(pred.number)
  perf.list<-list()
  legend.list<-c()
  name.list<-names(predictor.matrix)
  rocobj.list<-list()
  for(i in 1:pred.number){
    pred<-prediction(predictor.matrix[,i],status)
    perf <- performance(pred, "tpr", "fpr" )
    auc <- round(performance(pred,"auc")@y.values[[1]][1],digits = 2)
    #make sure the relative value 
    if(auc<0.5){
      auc=1-auc
#       status2<-factor(status,levels = rev(levels(as.factor(status))))
#       pred<-prediction(predictor.matrix[,i],status2)
#       perf <- performance(pred, "tpr", "fpr" )
#       auc <- round(performance(pred,"auc")@y.values[[1]][1],digits = 2)
    }
    perf.list[[length(perf.list)+1]]<-perf
    #get ci 
    rocobj <- roc(status, predictor.matrix[,i])
    rocobj.list[[length(rocobj.list)+1]]<-rocobj
    legend.list=c(legend.list,
                  paste(name.list[i],
                        "  ",auc,"  ",
                        paste(round(as.numeric(ci.auc(rocobj)),digits = 2)[c(1,3)],collapse="-")))
  }
  #ROC of each performance 
  ROCR::plot( perf.list[[1]], box.lty=1, box.lwd=2,col=col.list[1],
              box.col="black", lwd=2, colorkey.relwidth=0.5, xaxis.cex.axis=1,
              xaxis.col='black', xaxis.col.axis="black", yaxis.col='black', yaxis.cex.axis=1,
              yaxis.at=c(0,0.5,0.8,0.85,0.9,1), yaxis.las=1, xaxis.lwd=1, yaxis.lwd=1,
              ylim=c(0, 1), xlim=c(0,1),
              yaxis.col.axis="black", cex.lab=1.5, cex.main=1,title=predictor.name)
  for(i in 2:pred.number)
  {
    ROCR::plot(perf.list[[i]], add = TRUE,colorize = F,col=col.list[i],lwd=2)
  }
  print(length(legend.list))
  legend(x=0.3,y=0.3, legend = legend.list , col=col.list, lty=1, cex=1)
  abline(0,1,col="red") 
  
  #get test pvalue 
  for(i in 2:pred.number)
  {
    p<-roc.test(rocobj.list[[1]], rocobj.list[[i]], paired=T, method="delong")
  print(paste(name.list[i]," vs ",name.list[i]," : ",round(p$p.value,digits = 4)))
  }
}
#get ROC plot of individual predictor
plotSingle_ROC_with_test<- function(predictor,status,predictor.name="predictor"){
  suppressMessages(library("ROCR"))
  pred <- prediction(predictor,status)
  pref <- performance(pred, "tpr", "fpr" )
  auc <- round(performance(pred,"auc")@y.values[[1]][1],digits = 4)
  ROCR::plot( pref, box.lty=1, box.lwd=2,col="#57719D",
              box.col="black", lwd=2, colorkey.relwidth=0.5, xaxis.cex.axis=1,
              xaxis.col='black', xaxis.col.axis="black", yaxis.col='black', yaxis.cex.axis=1,
              yaxis.at=c(0,0.5,0.8,0.85,0.9,1), yaxis.las=1, xaxis.lwd=1, yaxis.lwd=1,
              ylim=c(0, 1), xlim=c(0,1),
              yaxis.col.axis="black", cex.lab=1.5, cex.main=1,title=predictor.name)
  legend(x=0.5,y=0.2, legend =paste(predictor.name," AUC =",auc) , col='#57719D', lty=1, cex=1)
  abline(0,1,col="red") 
}

#get print version of HR 
get_print_HR_cox<-function(surv,variable){
  df<-prettify(summary(coxph(surv~variable)))
  
  a=paste(round(df$`Hazard Ratio`,digits = 2)," (",round(df$`CI (lower)`,digits = 2),"-",round(df$`CI (upper)`,digits = 2),")",sep="")
  return(a)
  
}
paste(analysisDF$`Hazard Ratio`," (",analysisDF$`CI (lower)`,"-",analysisDF$`CI (upper)`,")",sep="")

#customized survival plot 

cust_surp<-function(fit,df){
  ggsurvplot(fit, data = df,
             surv.median.line = "hv", # Add medians survival
             
             # Change legends: title & labels
             legend.title = "Risk",
             legend.labs = c("Low", "High"),
             # Add p-value and confidence intervals
             pval = TRUE,
             
             conf.int = TRUE,
             # Add risk table
             risk.table = TRUE,
             tables.height = 0.2,
             tables.theme = theme_cleantable(),
             
             # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
             # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
             palette = c("#E7B800", "#2E9FDF"),
             ggtheme = theme_bw() # Change ggplot2 theme
  )
}



