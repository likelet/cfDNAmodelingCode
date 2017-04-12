
cDep <- c("ROCR","mice","pROC","glmnet","caret","papeR","ggsci","ggthemes")

###INSTALLED PACKAGES
#get installed list
inst <- packageStatus()$inst

#check and install DEPENDENCIES from CRAN
for(i in 1:length(cDep)){
  tag = which(inst$Package == cDep[i])
  if(length(tag)){
    remove.packages(cDep[i])
  }
  install.packages(cDep[i])
}

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
#
library(ggthemes)

################preloading functions ################

#get wilox test among nearby groups 
wiloxtestPair<-function(values,groupV){
  groupV=as.factor(groupV)
  gnames=names(table(groupV))
  for(i in 1:(length(gnames)-1)){
    print(paste(" ",gnames[i]," vs ",gnames[i+1],
                wilcox.test(
                  values[groupV %in% gnames[i]], 
                  values[groupV %in% gnames[i+1]]
                  , conf.int = TRUE)$p.value,sep=" "))
  }
}

##########multiple plot function from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
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

########   Start analysis ###########################

setwd("~/data/luo")



############ Paired Sample analysis T & Cancer Blood##################
#Get detectable markers ----
expfra=0.1 #minimum expression threshold 
LBpairedDatamatrix=read.table("hisq13_imputed_blocks_liver_blood_and_tissue.txt",header=T, row.names=1,check.names = F)
LBpairdSampleListDF=read.table("tissue-serrum match sample list.txt",header=T, check.names = F,sep="\t")

LBpairedDatamatrix=t(LBpairedDatamatrix)
LBpairedDatamatrix=data.frame(LBpairedDatamatrix,ID=row.names(LBpairedDatamatrix))
LBpairdSampleListDF=cbind(LBpairdSampleListDF,index=row.names(LBpairdSampleListDF))

LBpairdSampleCancerDF=merge(LBpairedDatamatrix,LBpairdSampleListDF,by.x="ID",by.y="Cancer")
LBpairdSampleNormalDF=merge(LBpairedDatamatrix,LBpairdSampleListDF,by.x="ID",by.y="Normal")
corlist=c()
meanlistCamcer=c()
meanlistNormal=c()
for (i in 2:653) {
  corlist=c(corlist,cor(
    LBpairdSampleCancerDF[,i],LBpairdSampleNormalDF[,i],method ='spearman'))
  meanlistCamcer=c(meanlistCamcer,mean(LBpairdSampleCancerDF[,i]))
  meanlistNormal=c(meanlistNormal,mean(LBpairdSampleNormalDF[,i]))
}
cordf=data.frame(marker=names(LBpairdSampleCancerDF)[2:653],Cor=corlist,mc=meanlistCamcer,mn=meanlistNormal)
cordf$marker=as.character(cordf$marker)
detectableMarkers=intersect(cordf[(cordf$mc>expfra),"marker"],cordf[(cordf$mn>expfra),"marker"])

#plot paired sample dataset -----(optional)
names(LBpairdSampleNormalDF)[654]="Type"
names(LBpairdSampleCancerDF)[654]="Type"                            
heatmapdf=rbind(LBpairdSampleCancerDF,LBpairdSampleNormalDF)

row.names(heatmapdf)=heatmapdf$ID
markersannotation=data.frame(names(heatmapdf)[2:653],Ispassed=(names(heatmapdf)[2:653]%in% detectableMarkers))

heatmapplotDF=t(heatmapdf[,2:653])
heatmapplotDF=data.frame(heatmapplotDF,ispass=as.character(markersannotation$Ispassed))
heatmapplotDF=heatmapplotDF[order(heatmapplotDF$ispass),]
markersannotation=data.frame(Ispassed=heatmapplotDF$ispass)
row.names(markersannotation)=row.names(heatmapplotDF)

pheatmap(heatmapplotDF[,-c(57,58)],annotation_row=markersannotation,cluster_rows = T,cluster_cols=F,scale="none",show_rownames = F)


# start analysis with unpaired data----
#preprocess raw data
data<-read.table("hisq13_imputed_liver_normal_allMeth.tsv",header=T,sep="\t",row.names=1)
#remove bad samples
retainedSamplelist=read.table("retained_sample.txt")
data=data[row.names(data) %in% retainedSamplelist$V1,]




#seperate matrix into several predicting value-status subset(optional)----
invisible(readline(prompt="Press [enter] to continue"))
print("Alldata discrimination result for each marker")
for (i in 2:589){
  newdf=data[,c(i,1)]
  #omit NA samples
  newdf=newdf[!is.na(newdf[,1]),]
  pred=prediction(newdf[,1],newdf[,2])
  perf <- performance(pred, "tpr", "fpr")
  auc=performance(pred,"auc")@y.values[[1]][1]
  
  print(paste(colnames(data)[i],"    ",auc,"    ",nrow(newdf),"    ",table(newdf[,2])[1],"    ",table(newdf[,2])[2]))
  # plot(perf)
}
invisible(readline(prompt="Press [enter] to continue"))
#



#random devided dataset into two dataset (optional when reproducted analysis) ----
smp_size <- floor(0.66 * nrow(data))
train_idx=sample(1:nrow(data),smp_size,rep=F)
trainingData=data[train_idx,]
ValidationData=data[-train_idx,]
#save and reloading 2 dataset (code was removed)


#reread splited data set to replicate analysis---------- 
trainingData2=read.csv("traningData2.csv",row.names=1,header=T,sep=" ")
ValidationData2=read.csv("validataionData2.csv",row.names=1,header=T,sep=" ")
trainingData2=trainingData2[,-c(590:593)]

ValidationData2=ValidationData2[,-c(590:593)]
# retain high detectable makers----------- 
availablemarkers=intersect(names(trainingData2),detectableMarkers)
trainingData2=trainingData2[,which (names(trainingData2) %in%　c("Status",availablemarkers))]
ValidationData2=ValidationData2[,which (names(trainingData2) %in%　c("Status",availablemarkers))]

## get 401 markers 
#Lasso analysis with current markers----------
#get lasso selected variables using subsampling and bagging frame work
selecVlist2=c()
for(i in 1:500){
  sampleindex2=sample(1:nrow(trainingData2),0.75*nrow(trainingData2),rep=F)
  effectdata=trainingData[sampleindex2,]
  
  glmmod<-glmnet(as.matrix(effectdata[2:402]),y=as.factor(effectdata$Status),alpha=1,family='binomial')
  cv.glmmod<-cv.glmnet(as.matrix(effectdata[2:402]),y=as.factor(effectdata$Status),alpha=1,family='binomial')
  best_lambda <- cv.glmmod$lambda.1se
  result<-coef(glmmod, s = best_lambda)
  selecVlist2=c(selecVlist2,result@Dimnames[[1]][which(result != 0)])
  print(i)
}
tablecount1=table(selecVlist2)
markerslasso=names(tablecount1[tablecount1>450])[-1]

#Random forest with markers paralell with lasso analysis----
#random forest and varSelRF to select variables 
#see citation https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-328
library(varSelRF)
trainingData$Status=as.factor(trainingData2$Status)
training.rf=randomForest(Status ~ ., data=trainingData2, importance=TRUE)
# rf.rvi <-randomVarImpsRF(trainingData[,-1], trainingData$Status, 
#                          training.rf, 
#                          numrandom = 20, 
#                          usingCluster = FALSE)
rf.vs2<-varSelRF(trainingData2[,-1], trainingData2$Status,  ntree = 3000, ntreeIterat = 2000,
                 vars.drop.frac = 0.3, whole.range = FALSE,
                 keep.forest = TRUE)

markers.RF=rf.vs2$selected.vars

#overlap markers selected by two methods above ----
markers=intersect(markerslasso,markers.RF)



trainingData$Status=as.factor(trainingData$Status)
#Multivariable regression Ananlysis using the selected markers---
fmla <- as.formula(paste("Status ~ ", paste(markers, collapse= "+")))
trainingFit=glm(fmla, data = trainingData, family = "binomial")
summary(trainingFit)

#get roc of combine score from 13 markers--------
training_score=predict(trainingFit,trainingData)
Validation_score=predict(trainingFit,ValidationData)

#print performance of training dataset
# invisible(readline(prompt="Press [enter] to continue visualizing ROC step"))
pred1=prediction(training_score,trainingData$Status)
perf1 <- performance(pred1, "tpr", "fpr")

#plot ROC using combineRisk Score in training dataset
plot(perf1, main="Performance of combinescore in traning dataset",
     box.lty=2, box.lwd=2,
     box.col="black", lwd=2, colorkey.relwidth=0.5, xaxis.cex.axis=1,
     xaxis.col='blue', xaxis.col.axis="blue", yaxis.col='green', yaxis.cex.axis=1,
     yaxis.at=c(0,0.5,0.8,0.85,0.9,1), yaxis.las=1, xaxis.lwd=1, yaxis.lwd=1,
     yaxis.col.axis="orange", cex.lab=1.5, cex.main=1)
auc1=performance(pred1,"auc")@y.values[[1]][1]
print(paste("traning    ",auc1,"    "))
print(ci.auc(roc(trainingData$Status, training_score)))

#plot ROC using combineRisk Score in validation dataset
invisible(readline(prompt="Press [enter] to continue visualizing ROC in validataion data "))
pred2=prediction(Validation_score,ValidationData$Status)
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2, main="Performance of combinescore in Validation dataset",
     box.lty=2, box.lwd=2,
     box.col="black", lwd=2, colorkey.relwidth=0.5, xaxis.cex.axis=1,
     xaxis.col='blue', xaxis.col.axis="blue", yaxis.col='green', yaxis.cex.axis=1,
     yaxis.at=c(0,0.5,0.8,0.85,0.9,1), yaxis.las=1, xaxis.lwd=1, yaxis.lwd=1,
     yaxis.col.axis="orange", cex.lab=1.5, cex.main=1)
auc2=performance(pred2,"auc")@y.values[[1]][1]
print(paste("Validation    ",auc2,"    "))
print(ci.auc(roc(ValidationData$Status, Validation_score)))
invisible(readline(prompt="Press [enter] to continue"))


#get confusion Table-------
trainingData$Status=as.character(trainingData$Status)
my_roc <- roc(trainingData$Status, training_score)
cutoff=coords(my_roc, "best", ret = "threshold")
Pred_status_train=training_score
Pred_status_train[training_score>cutoff]="Pre_True"
Pred_status_train[training_score<=cutoff]="Pre_False"
status_train=trainingData$Status
status_train[trainingData$Status=="1"]="True"
status_train[trainingData$Status=="0"]="False"

Pred_status_validataion=Validation_score
Pred_status_validataion[Validation_score>cutoff]="Pre_True"
Pred_status_validataion[Validation_score<=cutoff]="Pre_False"
status_validation=ValidationData$Status
status_validation[ValidationData$Status=="1"]="True"
status_validation[ValidationData$Status=="0"]="False"


trainingData1=cbind(trainingData,Pred_status_train,status_train)
ValidationData1=cbind(ValidationData,Pred_status_validataion,status_validation)

###generate confusion Table 
table_training=table(Pred_status_train,status_train)
table_validataion=table(Pred_status_validataion,status_validation)


#heatmap
library(pheatmap)
#dataset Traingdataset
annotation_row=data.frame(Status=as.factor(status_train),Predict=as.factor(Pred_status_train))

rownames(annotation_row)=row.names(trainingData1)
trainingData=trainingData1[with(trainingData1, order(Status)),]
pheatmap(trainingData[,markers],
         color = colorRampPalette(c("green","black","red"))(50),
         annotation_row=annotation_row,
         
         # cluster_rows=T,
         show_rownames =F,scale="none")

annotation_row=data.frame(Status=as.factor(status_validation),Predict=as.factor(Pred_status_validataion))

rownames(annotation_row)=row.names(ValidationData1)
ValidationData=ValidationData[with(ValidationData1, order(Status)),]
pheatmap(ValidationData[,markers],
         color = colorRampPalette(c("green","black","red"))(50),
         annotation_row=annotation_row,show_rownames =F)





###########Exploring the clinical Relavence of combine scores in LIHC 
clinicalData=read.table("clinicalData.txt",header=T,sep="\t")
names(clinicalData)[1]="ID"
data1=cbind(data,ID=row.names(data))
datawithCli=merge(data1,clinicalData,by="ID")
#get predict combinescore using trained model
allcombinescore=predict(trainingFit,datawithCli)
datawithCli=cbind(datawithCli,cscore=allcombinescore)
#get normal data for integrating into plot
normaldata=data[data1$Status=='0',]
normaldataCscore=predict(trainingFit,normaldata)
normaldata=cbind(normaldata,ID=row.names(normaldata))
normaldata=cbind(normaldata,cscore=normaldataCscore)
normaldata$Status=rep("Normal",nrow(normaldata))


#violin plot by ggplot2

#order levels Days After surgery
responsedf1=datawithCli[datawithCli$response %in% c('No treatment','PR','SD','PD'),]
#remove CR case time less than 1M
responsedf1=responsedf1[!responsedf1$DayAS %in% '1M',]
responsedf1=responsedf1[,c("response","cscore")]
responsedf1[responsedf1$response=='SD',"response"]='PR'
names(responsedf1)[1]="Status"
respon1plotdf=rbind(responsedf1,normaldata[,c("Status","cscore")])
respon1plotdf$Status=factor(respon1plotdf$Status,
                            levels = c("Normal",'No treatment','PR','PD'),ordered = TRUE)
response1P=ggplot(respon1plotdf, aes(x=Status, y=cscore,fill=Status)) + 
  # geom_violin(trim=F)+
  labs(title="",x="Response", y = "Combine Score")+
  geom_boxplot(width=0.5)+
  scale_fill_npg()+theme_base()+ylim(-10,25)+ theme(legend.position="none")




#order levels response,add normaldata
responsedf2=datawithCli[datawithCli$response %in% c('No treatment','CR','RC'),]
#remove CR case time less than 1M
responsedf2=responsedf2[!responsedf2$DayAS %in% '1M',]
responsedf2=responsedf2[,c("response","cscore")]


names(responsedf2)[1]="Status"
responsedf2$Status=as.character(responsedf2$Status)
respon2plotdf=rbind(responsedf2,normaldata[,c("Status","cscore")])
respon2plotdf$Status=factor(respon2plotdf$Status,
                            levels = c("Normal",'No treatment','CR','RC'),ordered = TRUE)
response2P=ggplot(respon2plotdf, aes(x=Status, y=cscore,fill=Status)) + 
  # geom_violin(trim=F)+
  labs(title="",x="Response", y = "Combine Score")+
  geom_boxplot(width=0.5)+
  scale_fill_npg()+theme_base()+ylim(-10,25)+ theme(legend.position="none")


#order levels by tumor burden 
tumorburderndf=datawithCli[,c("tumor.burdeng","cscore")]
names(tumorburderndf)[1]="Status"
tumorburdernplotdf=rbind(tumorburderndf,normaldata[,c("Status","cscore")])
tumorburdernplotdf$Status=factor(tumorburdernplotdf$Status,
                                 levels = c('Normal','N','Y'),ordered = TRUE)
tumorburdernplotdf=tumorburdernplotdf[!is.na(tumorburdernplotdf$Status),]
tumorburdernP=ggplot(tumorburdernplotdf, aes(x=Status, y=cscore,fill=Status)) + 
  # geom_violin(trim=F)+
  labs(title="",x="Tumor Burden", y = "Combine Score")+
  geom_boxplot(width=0.5)+
  scale_fill_npg()+theme_base()+ylim(-10,25)+ theme(legend.position="none")



#order levels by time after surgery
dayasdf=datawithCli[,c("DayAS","cscore")]
dayasdf$DayAS=as.character(dayasdf$DayAS)
dayasdf[dayasdf$DayAS %in% c('1-3M','3-6M'),"DayAS"]='1-6M'
dayasdf[dayasdf$DayAS %in% c('6M'),"DayAS"]='>6M'
names(dayasdf)[1]="Status"
dayasplotdf=rbind(dayasdf,normaldata[,c("Status","cscore")])
dayasplotdf$Status=factor(dayasplotdf$Status,
                          levels = c('Normal','1M','1-6M','>6M'),ordered = TRUE)
dayasplotdf=dayasplotdf[!is.na(dayasplotdf$Status),]
dayasP=ggplot(dayasplotdf, aes(x=Status, y=cscore,fill=Status)) + 
  # geom_violin(trim=F)+
  labs(title="",x="Tumor Burden", y = "Combine Score")+
  geom_boxplot(width=0.5)+
  scale_fill_npg()+theme_base()+ylim(-10,25)+ theme(legend.position="none")




#order levels by Stage  
# datawithCli$stage=factor(datawithCli$stage,
# levels = c('1','2','3','4'),ordered = TRUE)
stagedf=datawithCli[,c("stage","cscore")]
stagedf$stage=as.character(stagedf$stage)
stagedf[stagedf$DayAS %in% c('1','2'),"stage"]='1'
stagedf[stagedf$DayAS %in% c('3','4'),"stage"]='2'
names(stagedf)[1]="Status"
stageplotdf=rbind(stagedf,normaldata[,c("Status","cscore")])
stageplotdf=stageplotdf[!is.na(stageplotdf$Status),]
stageplotdf$Status=factor(stageplotdf$Status,
                          levels = c('Normal','1','2'),ordered = TRUE)
stageP=ggplot(stageplotdf, aes(x=Status, y=cscore,fill=Status)) + 
  # geom_violin(trim=T)+
  labs(title="",x="Tumor Stage", y = "Combine Score")+
  geom_boxplot(width=0.5)+
  scale_fill_npg()+theme_base()+ylim(-10,25)+ theme(legend.position="none")

#stageSep
#remove CR
stageSepdf=datawithCli[!datawithCli$response %in% 'CR',c("stage","cscore")]
stageSepdf$stage=as.character(stageSepdf$stage)

names(stageSepdf)[1]="Status"
stageSepdf=rbind(stageSepdf,normaldata[,c("Status","cscore")])
stageSepplotdf=stageSepdf[!is.na(stageSepdf$Status),]
stageSepplotdf$Status=factor(stageSepplotdf$Status,
                             levels = c('Normal','1','2','3','4'),ordered = TRUE)
stageSepP=ggplot(stageSepplotdf, aes(x=Status, y=cscore,fill=Status)) + 
  # geom_violin(trim=T)+
  labs(title="",x="Tumor Stage", y = "Combine Score")+
  geom_boxplot(width=0.5)+
  scale_fill_npg()+theme_base()+ylim(-10,25)+ theme(legend.position="none")


#plot 4 figures into one graph
multiplot(response1P, response2P, tumorburdernP, stageSepP, cols=2)
#test difference 

#Test Respons1 difference using wilcox test
wiloxtestPair(respon1plotdf$cscore,respon1plotdf$Status)

#Test Respons2 difference using wilcox test
wiloxtestPair(respon2plotdf$cscore,respon2plotdf$Status)
#Test TumorBurden difference using wilcox test
wiloxtestPair(tumorburdernplotdf$cscore,tumorburdernplotdf$Status)
#Test DaysAS difference using wilcox test
wiloxtestPair(dayasplotdf$cscore,dayasplotdf$Status)
#Test Stage difference using wilcox test
# wiloxtestPair(stageplotdf$cscore,stageplotdf$Status)
#seperated data 
wiloxtestPair(stageSepplotdf$cscore,stageSepplotdf$Status)

####################  Associated with AFP && hepatasis   ############################################
afpdata=read.table("newSample/AFPinfo.txt",header=T,check.names = F)
alldataCbineScore=rbind(datawithCli[,c("ID","cscore")],normaldata[,c("ID","cscore")])
datawithCliAFP=merge(alldataCbineScore,afpdata,by="ID")
datawithCliAFP$Status=as.character(datawithCliAFP$hepatitis)
datawithCliAFP[!datawithCliAFP$Status %in% "Tumor","Status"]="0"
datawithCliAFP[datawithCliAFP$Status %in% "Tumor","Status"]="1"
datawithCliAFP$Status=as.factor(datawithCliAFP$Status)
#plot 2 ROC on the same graph

predCscore=prediction(datawithCliAFP$cscore,datawithCliAFP$Status)
predAFP=prediction(datawithCliAFP$AFP,datawithCliAFP$Status)
# get auc CI
performance(predCscore,"auc")@y.values[[1]][1]
print(ci.auc(datawithCliAFP$Status,datawithCliAFP$cscore))
performance(predAFP,"auc")@y.values[[1]][1]
print(ci.auc(datawithCliAFP$Status,datawithCliAFP$AFP))

#plot ROC
perfCscore <- performance( predCscore, "tpr", "fpr" )
perfAFP <- performance(predAFP, "tpr", "fpr")
plot( perfCscore, box.lty=2, box.lwd=2,
      box.col="black", lwd=2, colorkey.relwidth=0.5, xaxis.cex.axis=1,
      xaxis.col='blue', xaxis.col.axis="blue", yaxis.col='green', yaxis.cex.axis=1,
      yaxis.at=c(0,0.5,0.8,0.85,0.9,1), yaxis.las=1, xaxis.lwd=1, yaxis.lwd=1,
      yaxis.col.axis="orange", cex.lab=1.5, cex.main=1)
plot(perfAFP, add = TRUE,colorize = TRUE)
# plot boxplot for compare score distribution among case with or without hepastasis
datawithCliAFP$hepatitis=factor(datawithCliAFP$hepatitis,
                                levels = c('N','Y','Tumor'),ordered = TRUE)
hepP=ggplot(datawithCliAFP, aes(x=hepatitis, y=cscore,fill=hepatitis)) + 
  # geom_violin(trim=T)+
  labs(title="",x="Hepatitis", y = "Combine Score")+
  geom_boxplot(width=0.5)+
  scale_fill_npg()+theme_base()+ylim(-10,25)+ theme(legend.position="none")
# do wilox test  between groups 

wiloxtestPair(datawithCliAFP$cscore,datawithCliAFP$hepatitis)

# plot boxplot for AFP distribution among case with or without hepastasis

afpP=ggplot(datawithCliAFP, aes(x=hepatitis, y=AFP,fill=hepatitis)) + 
  # geom_violin(trim=T)+
  labs(title="",x="Hepatitis", y = "AFP")+
  geom_boxplot(width=0.5)+
  scale_y_log10()+
  scale_fill_npg()+theme_base()+ theme(legend.position="none")
# do wilox test  between groups 
wiloxtestPair(datawithCliAFP$AFP,datawithCliAFP$hepatitis)

#combine plot 
multiplot(hepP, afpP, cols=2)


################### samples with all clinical data ####################
#merge data
datawithCliAFP1=datawithCliAFP[,c(1,3,4)]
DatawithAllCli=merge(datawithCliAFP1,datawithCli,by="ID")

#subset the plot data
plotcomparisonData=DatawithAllCli[,c("ID","cscore","AFP","stage","response")]
#define PR SD as after treatment 
plotcomparisonData$response=as.character(plotcomparisonData$response)

plotcomparisonData[plotcomparisonData$response %in% c("PR","SD"),"response"]="AF"

# plotcomparisonData$response=factor(plotcomparisonData$response,
#                           levels = c('No treatment','AF','CR','RC',"PD"),ordered = TRUE)

#plot correlation boxplot

#score Res
scoreResdf=plotcomparisonData[!is.na(plotcomparisonData$response),]
scoreResdf$response=as.character(scoreResdf$response)
# scoreResdf=scoreResdf[!scoreResdf$response %in% "SD",]
scoreResdf$response=factor(scoreResdf$response,
                           levels = c('No treatment','AF','CR','RC',"PD"),ordered = TRUE)
scoreResp=ggplot(scoreResdf, aes(x=response, y=cscore,fill=response)) + 
  # geom_violin(trim=F)+
  labs(title="",x="Treatmant", y = "Combine Score")+
  geom_boxplot(width=0.5)+
  scale_fill_npg()+theme_base()+ylim(-10,25)+ theme(legend.position="none")
#AFP Res
afpRespdf=plotcomparisonData[!is.na(plotcomparisonData$response),]
afpRespdf$response=as.character(afpRespdf$response)
# afpRespdf=afpRespdf[!afpRespdf$response %in% "SD",]
afpRespdf$response=factor(afpRespdf$response,
                          levels = c('No treatment','AF','CR','RC',"PD"),ordered = TRUE)
afpResp=ggplot(afpRespdf, aes(x=response, y=AFP,fill=response)) + 
  # geom_violin(trim=F)+
  labs(title="",x="Treatmant", y = "Log AFP")+
  geom_boxplot(width=0.5)+
  scale_y_log10()+
  scale_fill_npg()+theme_base()+ theme(legend.position="none")

#score stage
#score Res
scoreResdf=plotcomparisonData[!is.na(plotcomparisonData$response),]
scoreResdf$response=as.character(scoreResdf$response)
scoreResdf[scoreResdf$response %in% c("AF","CR"),"response"]="CR+PR+SD"
# scoreResdf=scoreResdf[!scoreResdf$response %in% "SD",]
scoreResdf$response=factor(scoreResdf$response,
                           levels = c('No treatment','PD','CR+PR+SD',"RC"),ordered = TRUE)
scoreResp=ggplot(scoreResdf, aes(x=response, y=cscore,fill=response)) + 
  # geom_violin(trim=F)+
  labs(title="",x="Treatmant", y = "Combine Score")+
  geom_boxplot(width=0.5)+
  scale_fill_npg()+theme_base()+ylim(-10,25)+ theme(legend.position="none")
#AFP Res
afpRespdf=plotcomparisonData[!is.na(plotcomparisonData$response),]
afpRespdf$response=as.character(afpRespdf$response)
# afpRespdf=afpRespdf[!afpRespdf$response %in% "SD",]
afpRespdf[afpRespdf$response %in% c("AF","CR"),"response"]="CR+PR+SD"
# scoreResdf=scoreResdf[!scoreResdf$response %in% "SD",]
afpRespdf$response=factor(afpRespdf$response,
                          levels = c('No treatment','PD','CR+PR+SD',"RC"),ordered = TRUE)

afpResp=ggplot(afpRespdf, aes(x=response, y=AFP,fill=response)) + 
  # geom_violin(trim=F)+
  labs(title="",x="Treatmant", y = "Log AFP")+
  geom_boxplot(width=0.5)+
  scale_y_log10()+
  scale_fill_npg()+theme_base()+ theme(legend.position="none")

#test difference 
wiloxtestPair(scoreResdf$cscore,scoreResdf$response)
wiloxtestPair(afpRespdf$AFP,afpRespdf$response)
wiloxtestPair(scoreStageplotDF$cscore,scoreStageplotDF$stage)
wiloxtestPair(scoreStageplotDF$AFP,scoreStageplotDF$stage)

multiplot(scoreResp,scoreStagep,afpResp,afpStageP, cols=2)


#write pairedData as well as combinescore with 10markers(optional)----
cscorePairedTorTB=predict(trainingFit,heatmapdf)
write.csv(data.frame(heatmapdf,cscorePairedTorTB),"PairedDataTandTB.csv")


#Modified heatmap
#rearranged sample names like "T N T N T N"  and add normal samples into matrix
cooccuranceMarker=intersect(cordf$marker,names(trainingData))

forPairNormaldata=trainingData[trainingData$Status==0,which (names(trainingData) %in% cooccuranceMarker)]
row.names(forPairNormaldata)=row.names(trainingData[trainingData$Status==0,])
#get 100 normalsample 
forPairNormaldata=forPairNormaldata[sample(1:560,30,rep=F),]
forPairNormaldata=data.frame(t(forPairNormaldata))
forPairNormaldata=cbind(forPairNormaldata,ID=row.names(forPairNormaldata))
forPairNormaldata$ID=as.character(forPairNormaldata$ID)

# heatmapplotDF=cbind(heatmapplotDF,ID=row.names(heatmapplotDF))

# subheatmapplotDF=heatmapplotDF[which (as.character(heatmapplotDF$ID) %in% cooccuranceMarker),]
finalHeatmapdf=merge(heatmapplotDF,forPairNormaldata,by="ID")
row.names(finalHeatmapdf)=as.character(finalHeatmapdf$ID)

finalHeatmapdf=finalHeatmapdf[,-58]
finalHeatmapdf$mt=apply(finalHeatmapdf[,2:39],1,mean)
finalHeatmapdf$mn=apply(finalHeatmapdf[,30:57],1,mean)
finalHeatmapdf$ispass=finalHeatmapdf$mt>0.1 & finalHeatmapdf$mt >0.1
finalHeatmapdf=finalHeatmapdf[order(finalHeatmapdf$ispass),]
markersannotation2=data.frame(Ispassed=as.character(finalHeatmapdf$ispass))
row.names(markersannotation2)=finalHeatmapdf$ID
markersannotationCol=data.frame(Sampletype=c(rep("Tumor tissue",28),rep("Tumor blood",28),rep("Normal tissue",30)))
row.names(markersannotationCol)=names(finalHeatmapdf)[2:87]
pheatmap(finalHeatmapdf[,2:87],annotation_row=markersannotation2,annotation_col=markersannotationCol,cluster_rows = F,cluster_cols=F,scale="none",show_rownames = F)

