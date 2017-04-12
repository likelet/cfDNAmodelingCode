#for ROC plot
library(ROCR)
#for data imputation
library(mice) 
#for best cutoff 
library(pROC)
#for LASSO 
library(glmnet)
#for confusion table 
library(caret)
#for pretty output 
library(papeR)
#theme for plot 
library(ggsci)
library(ggthemes)
library(ggplot2)
library(survival)
library(risksetROC)
library(rms)

setwd("~/data/Survival/")
#read new survival data
survivalData=read.table("surivalInfo.txt",header=T,sep="\t")




#training data and validation dataset provided by binary analysis section which retained available markers 
trainingData3=cbind(trainingData2,ID=row.names(trainingData2))

ValidationData3=cbind(ValidationData2,ID=row.names(ValidationData2))


#merge two dataset seperately
trainingData2M=merge(trainingData3,survivalData,by="ID")
ValidationData2M=merge(ValidationData3,survivalData,by="ID")
names(trainingData2M)[414]="Time"
names(ValidationData2M)[414]="Time"
names(trainingData2M)[413]="Status"
names(ValidationData2M)[413]="Status"
trainingData2M=trainingData2M[!is.na(trainingData2M$Time),]
ValidationData2M=ValidationData2M[!is.na(ValidationData2M$Time),]



#univariable cox analysis for screening variables 
batchUnivarCOXfun<-function(surv,data){
  df<-data.frame(name=c("coef","Hazard Ratio","CI (lower)","CI (upper)","se(coef)","z","Pr(>|Z|)","C-index","C-index-se","testCph.p"))
  for (i in 1:ncol(data)) {
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
  }
  row.names(df)<-df[,1]
  df<-t(df[,-1])
  row.names(df)=names(data)
  return(df)
}
uniCoxDf=batchUnivarCOXfun(Surv(trainingData2M$Time, trainingData2M$Status),trainingData2M[,3:403])

# write.csv(uniCoxDf,"Training_uniCoxDf.csv")

markers_first=row.names(uniCoxDf[uniCoxDf[,7]<0.05,])
row.names(trainingData2M)=trainingData2M$ID
trainingData2MforSurvival=trainingData2M[,names(trainingData2M) %in% c(markers_first,"Status","Time")]

write.csv(trainingData2MforSurvival,"trainingData2MforSurvival.csv")
trainingData2MforSurvival$Status=as.factor(trainingData2MforSurvival$Status)


# markersSurvival=intersect(markers,markers_first)

trainingData2M$Status=as.numeric(trainingData2M$Status)

#lasso cox for select variables
#Lasso Cox with bootstrap model 
selecVlist=c()
for(i in 1:500){
  sampleindex2=sample(1:nrow(trainingData2MforSurvival),0.75*nrow(trainingData2MforSurvival),rep=F)
  effectdata=trainingData2MforSurvival[sampleindex2,]
  table(effectdata$Status)
  fit <- glmnet(as.matrix(effectdata[,c(1:243)]), Surv(effectdata$Time, effectdata$Status), family = "cox", maxit = 200)
  cv.fit <- cv.glmnet(as.matrix(effectdata[,c(1:243)]), Surv(effectdata$Time, effectdata$Status), family = "cox", maxit = 200)

  best_lambda <- cv.fit$lambda.1se
  result<-coef(fit, s = best_lambda)
  selecVlist=c(selecVlist,result@Dimnames[[1]][which(result != 0)])
  print(i)
}
tablecount=table(selecVlist)

markersSurvival=names(tablecount[tablecount>50])
fmla <- as.formula(paste("Surv(trainingData2M$Time, trainingData2M$Status) ~ ", paste(markersSurvival, collapse= "+")))
cox.fit <- coxph(fmla, data=trainingData2M)
#get data1 combine score
combinescoreData1=predict(cox.fit,trainingData2M)
trainingData2M=cbind(combinescoreData1,trainingData2M)
write.csv(trainingData2M,"Train_survival_combinescore.csv")
#get data2 combine score 
combinescoreData2=predict(cox.fit,ValidationData2M)
ValidationData2M=cbind(combinescoreData2,ValidationData2M)
write.csv(ValidationData2M,"Vali_survival_combinescore.csv")

#multiple cox in box data1 and data2
#refine names 
names(trainingData2M)[412]=c("tumor_burdeng")
names(ValidationData2M)[412]=c("tumor_burdeng")
#build formula 
adjustFactor=c("combinescoreData1","AFP","stage","sex","age")
adjustFactor2=c("combinescoreData2","AFP","stage","sex","age")
fmla1 <- as.formula(paste("Surv(trainingData2M$Time, trainingData2M$Status) ~ ", paste(adjustFactor, collapse= "+")))
fmla2 <- as.formula(paste("Surv(ValidationData2M$Time, ValidationData2M$Status) ~ ", paste(adjustFactor2, collapse= "+")))

mcox.fit1 =coxph(fmla1, data=trainingData2M,na.action=na.omit)
mcox.fit2 =coxph(fmla2, data=ValidationData2M,na.action=na.omit)

#AUC plot for compare combinescoreData1 and AFP
#Training set
train_cscore.fit=coxph(Surv(trainingData2M$Time, trainingData2M$Status) ~ combinescoreData1, data=trainingData2M,na.action=na.omit)
train_AFP.fit=coxph(Surv(trainingData2M$Time, trainingData2M$Status) ~ AFP, data=trainingData2M,na.action=na.omit)
train_stage.fit=coxph(Surv(trainingData2M$Time, trainingData2M$Status) ~ stage, data=trainingData2M,na.action=na.omit)

train_union.fit=coxph(Surv(trainingData2M$Time, trainingData2M$Status) ~ combinescoreData1+stage, data=trainingData2M,na.action=na.omit)
res1=risksetROC(Stime=trainingData2M$Time, status=trainingData2M$Status,
                marker=train_cscore.fit$linear.predictors, predict.time=180, method="Cox",plot = FALSE) 
res2=risksetROC(Stime=trainingData2M[!is.na(trainingData2M$AFP),]$Time, status=trainingData2M[!is.na(trainingData2M$AFP),]$Status,
                marker=train_AFP.fit$linear.predictors, predict.time=180, method="Cox",plot = FALSE)
resstage1=risksetROC(Stime=trainingData2M[!is.na(trainingData2M$stage),]$Time, status=trainingData2M[!is.na(trainingData2M$stage),]$Status,
                    marker=train_stage.fit$linear.predictors, predict.time=180, method="Cox",plot = FALSE)
resUnion=risksetROC(Stime=trainingData2M[!is.na(trainingData2M$stage),]$Time, status=trainingData2M[!is.na(trainingData2M$stage),]$Status,
                marker=train_union.fit$linear.predictors, predict.time=180, method="Cox",plot = FALSE)
ROCplotDF=data.frame(TP=res1$TP,FP=res1$FP,Maker=rep(paste("Cscore AUC:",round(res1$AUC,4)),length(res1$FP)))
ROCplotDF=rbind(ROCplotDF,
                data.frame(TP=res2$TP,FP=res2$FP,Maker=rep(paste("AFP AUC:",round(res2$AUC,4)),length(res2$FP))),
                data.frame(TP=resstage1$TP,FP=resstage1$FP,Maker=rep(paste("Stage AUC:",round(resstage1$AUC,4)),length(resstage1$FP)))
#                 data.frame(TP=resUnion$TP,FP=resUnion$FP,Maker=rep(paste("Cscore+stage:",round(resUnion$AUC,4)),length(resUnion$FP)))
#                 
)
ROCplotDF$Maker=factor(ROCplotDF$Maker,
                       levels=c(
                         # paste("Cscore+stage:",round(resUnion$AUC,4)),
                                paste("Cscore AUC:",round(res1$AUC,4)),
                                paste("Stage AUC:",round(resstage1$AUC,4)),
                                paste("AFP AUC:",round(res2$AUC,4))),
                       ordered = T)
ROCplot_trainp=ggplot(data=ROCplotDF,aes(x=FP,y=TP,color=Maker))+geom_line(size=1)+
  labs(title="ROC for 6 months survival predicting in training dataset",x="False positive rate", y = "True negative rate")+
  scale_color_npg()+theme_base()+theme(legend.position="top")+geom_abline(intercept = 0,linetype="dotted")+
  scale_y_continuous(limits=c(0,1),expand = c(0,0),breaks=c(0.2,0.4,0.6,0.8,1))+ 
  scale_x_continuous(limits=c(0,1),expand = c(0,0),breaks=c(0,0.2,0.4,0.6,0.8,1))+
  theme(
    legend.position = c(0.8, 0.25), # c(0,0) bottom left, c(1,1) top-right
    legend.background = element_rect(colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 15, face = 'bold'),
  )

#validation Set
valid_cscore.fit=coxph(Surv(ValidationData2M$Time, ValidationData2M$Status) ~ combinescoreData2, data=ValidationData2M,na.action=na.omit)
valid_AFP.fit=coxph(Surv(ValidationData2M$Time, ValidationData2M$Status) ~ AFP, data=ValidationData2M,na.action=na.omit)
validation_union.fit=coxph(Surv(ValidationData2M$Time, ValidationData2M$Status) ~ combinescoreData2+stage, data=ValidationData2M,na.action=na.omit)
validation_stage.fit=coxph(Surv(ValidationData2M$Time, ValidationData2M$Status) ~ stage, data=ValidationData2M,na.action=na.omit)

res3=risksetROC(Stime=ValidationData2M$Time, status=ValidationData2M$Status,
                marker=valid_cscore.fit$linear.predictors, predict.time=180, method="Cox",plot = FALSE) 

res4=risksetROC(Stime=ValidationData2M[!is.na(ValidationData2M$AFP),]$Time, status=ValidationData2M[!is.na(ValidationData2M$AFP),]$Status,
                marker=valid_AFP.fit$linear.predictors, predict.time=180, method="Cox",plot = FALSE)
resstage2=risksetROC(Stime=ValidationData2M[!is.na(ValidationData2M$stage),]$Time, status=ValidationData2M[!is.na(ValidationData2M$stage),]$Status,
                     marker=validation_stage.fit$linear.predictors, predict.time=180, method="Cox",plot = FALSE)

resUnion2=risksetROC(Stime=ValidationData2M[!is.na(ValidationData2M$stage),]$Time, status=ValidationData2M[!is.na(ValidationData2M$stage),]$Status,
                    marker=validation_union.fit$linear.predictors, predict.time=180, method="Cox",plot = FALSE)

ROCplotDF2=data.frame(TP=res3$TP,FP=res3$FP,Maker=rep(paste("Cscore AUC:",round(res3$AUC,4)),length(res3$FP)))
ROCplotDF2=rbind(ROCplotDF2,
                 data.frame(TP=res4$TP,FP=res4$FP,Maker=rep(paste("AFP AUC:",round(res4$AUC,4)),length(res4$FP))),
                 data.frame(TP=resstage2$TP,FP=resstage2$FP,Maker=rep(paste("Stage AUC:",round(resstage2$AUC,4)),length(resstage2$FP)))
#                  ,data.frame(TP=resUnion2$TP,FP=resUnion2$FP,Maker=rep(paste("Cscore+stage:",round(resUnion2$AUC,4)),length(resUnion2$FP)))
)
ROCplotDF2$Maker=factor(ROCplotDF2$Maker,
                        levels=c(
                          # paste("Cscore+stage:",round(resUnion2$AUC,4)),
                                 paste("Cscore AUC:",round(res3$AUC,4)),
                                 paste("Stage AUC:",round(resstage2$AUC,4)),
                                 paste("AFP AUC:",round(res4$AUC,4))),
                        ordered = T)
ROCplot_validationp=ggplot(data=ROCplotDF2,aes(x=FP,y=TP,color=Maker))+geom_line(size=1)+
  labs(title="ROC for 6 months survival predicting in validation dataset",x="False positive rate", y = "True negative rate")+
  scale_color_npg()+theme_base()+theme(legend.position="top")+geom_abline(intercept = 0,linetype="dotted")+
  scale_y_continuous(limits=c(0,1),expand = c(0,0),breaks=c(0.2,0.4,0.6,0.8,1))+ 
  scale_x_continuous(limits=c(0,1),expand = c(0,0),breaks=c(0,0.2,0.4,0.6,0.8,1))+
  theme(
    legend.position = c(0.8, 0.25), # c(0,0) bottom left, c(1,1) top-right
    legend.background = element_rect(colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 15, face = 'bold'),
  )


multiplot(ROCplot_trainp, ROCplot_validationp, cols=2)
#print dataset 

cpscore1=predict(train_union.fit,trainingData2M)
trainingData2MF=cbind(trainingData2M,cpscore1)
trainingData2MF[is.na(trainingData2MF$stage),"cpscore1"]=NA
write.csv(trainingData2MF,"trainingData2MF.csv")
rm(trainingData2MF)
cpscore2=predict(validation_union.fit,ValidationData2M)
ValidationData2MF=cbind(ValidationData2M,cpscore2)
write.csv(ValidationData2MF,"ValidationData2MF.csv")
ValidationData2MF[is.na(ValidationData2MF$stage),"cpscore2"]=NA
rm(ValidationData2MF)
