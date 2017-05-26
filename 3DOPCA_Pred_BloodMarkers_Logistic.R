#####
#Load Relevant Packages
install.packages("ROCR")
require("ROCR")


install.packages("AUC")
require("AUC")
#####
#Load Data
#Different Possible Paths based on Computer


for (sexLoop in 1:2){

setwd("W:\\3D Optical_Data\\Fit 3D\\Michelle Shape Up Optical Models\\0515 Updated Shape Models\\")
  
if (sexLoop==1){  
allData <- read.csv('combined_measurements_Male.csv', header = TRUE)
sexType <- "Male"
}else if (sexLoop ==2){
allData <- read.csv('combined_measurements_Female.csv', header = TRUE)
sexType <- "Female"
}
  
######
  #Add Thresholds to each thing   
        # for ("GLU", "CHOL", "TRIG", "HDL", "LDL","INS","AIC")
bMarkThresh <- c(100, 100, 100, 100, 100, 100, 6.5)
# Glucose
allData$GLUThresh[allData$GLU<=100]<-0  #Healthy
allData$GLUThresh[allData$GLU>100]<-1  #Prediabetic
# allData$GLUThresh[allData$GLU>126]<-2  #Diabetes
# Cholesterol
allData$CHOLThresh[allData$CHOL<=200]<-0
allData$CHOLThresh[allData$CHOL>200]<-1
allData$CHOLThresh[allData$CHOL>240]<-2
# Triglycerides
allData$TRIGThresh[allData$TRIG<=150]<-0  #Healthy
allData$TRIGThresh[allData$TRIG>150]<-1   #Moderate Risk
allData$TRIGThresh[allData$TRIG>200]<-2   #High Heart Disease Risk
# HDL
allData$HDLThresh[allData$HDL<40]<-1   #High risk Heart Disease
allData$HDLThresh[allData$HDL>=40]<-0   #Healthy
# LDL
allData$LDLThresh[allData$LDL>130]<-0
allData$LDLThresh[allData$LDL<=130]<-1
allData$LDLThresh[allData$LDL<100]<-2
# INS
allData$INSThresh[allData$INS>27]<-1  #Insulin Resistance
allData$INSThresh[allData$INS<=27]<-0  #Healthy
allData$INSThresh[allData$INS<7]<-0  #Diabetes
# AIC
allData$AICThresh[allData$AIC<=5.7]<-0  #Healthy
allData$AICThresh[allData$AIC>5.7]<-1  #Prediabetes
allData$AICThresh[allData$AIC>6.3]<-2  #Diabetes
#####
#Description of Data
print(colnames(allData))

keep <- c("SubjectID", "Ethnicity", "BMI","WaistCircum1", "HipCircum","UpperArmCircum",
          "GLUThresh", "CHOLThresh", "TRIGThresh", "HDLThresh", "LDLThresh","INSThresh","AICThresh",
          "WB_TOT_BMC", "WB_TOT_BMD","WB_TOT_FAT","WB_TOT_Lean","WB_TOT_PFAT","WB_VFAT_AREA","WB_SAT_AREA",
          "age","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","HEIGHT","WEIGHT" )
relData <- allData[keep]

#####
#Set Save Directory
setwd("X:\\bhinton\\ShapeUp\\PCA_To_BloodMarkersGender\\Analysis\\")
savePathBase <- ("X:\\bhinton\\ShapeUp\\PCA_To_BloodMarkersGender\\Analysis\\")
#####

#Perform analysis with Glucose
bMarkerVec <- c("GLU", "CHOL", "TRIG", "HDL", "LDL","INS","AIC")
bMarkerThresh <- c("GLUThresh", "CHOLThresh", "TRIGThresh", "HDLThresh", "LDLThresh","INSThresh","AICThresh")

#####
nFolds <-3
for (bMarkLoop in 1:7){
bMarker <- bMarkerVec[bMarkLoop]

savePath <- paste(savePathBase,bMarker, "\\", sep="")
setwd(savePath)
keep <- c(bMarkerThresh,'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7','PC8', 'PC9' )
pcsKeep <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7','PC8', 'PC9' )
predData <- na.omit(relData[keep])
#Prep for K-Folds
#Randomly shuffle the data (CC)
predData<-predData[sample(nrow(predData)),]
#Create 10 equally size folds (CC)
predDataFolds <- cut(seq(1,nrow(predData)),breaks=nFolds,labels=FALSE)

regResults <- data.frame(matrix(ncol = 6, nrow = 11))
colnames(regResults) <- c('train_r','train_r2', 'train_rmse', 'test_r', 'test_r2', 'test_rmse')

allPredicted <- NULL
allActual <- NULL
txtRegFile <-paste("Male_regression",bMarker,".txt", sep="")
fileConn<-file(txtRegFile)
for(j in 1:nFolds){
  testSetIndexes <- which(predDataFolds==j,arr.ind=TRUE)
  testSet <- rbind(predData[testSetIndexes, ])
  trainSet <- rbind(predData[-testSetIndexes, ])
  
  model <- lm(trainSet[,1]~., family=binomial(logit), data=trainSet[pcsKeep])
  modelStep <- step(model)
  trainRegSumm <- summary(modelStep)
  out <- capture.output(summary(modelStep))
  cat(out,file=txtRegFile,sep="\n",append=TRUE)
  
  #PLOT TRAINING SET AUC
  png(filename=paste(savePath,bMarker,"TrainAUCPlot",j,".png", sep=""))
  plot.new()
  p <- predict(model, newdata=trainSet[pcsKeep], type="response")
  pr <- prediction(p, trainSet[,1])
  prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  plot(prf, lty=1,lwd=5,cex = 1,cex.axis=2,cex.lab = 2,col="red",
            xlim=c(min(predData[,1]),max(predData[,1])), ylim=c(min(predData[,1]),max(predData[,1])),
            xlab = paste("FPR of ", bMarker, sep=""), ylab = paste("TPR of ", bMarker, sep=""))
  title(main = paste("Regression of ",sexType ," Train Set Version ", j, sep= ""))
  dev.off()
  
  fitted.results <- predict(modelStep,newdata=testSet,type='response')
  forced.results <- ifelse(fitted.results > 0.5,1,0)
  misClasificError <- mean(forced.results != testSet[,1])
  Accuracy<- 1-misClasificError
  
  png(filename=paste(savePath,bMarker,"TestAUCPlot",j,".png", sep=""))
  plot.new()
  p <- predict(model, newdata=testSet[pcsKeep], type="response")
  pr <- prediction(p, testSet[,1])
  prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  plot(prf, lty=1,lwd=5,cex = 1,cex.axis=2,cex.lab = 2,col="red",
       xlim=c(min(predData[,1]),max(predData[,1])), ylim=c(min(predData[,1]),max(predData[,1])),
       xlab = paste("FPR of ", bMarker, sep=""), ylab = paste("TPR of ", bMarker, sep=""))
  title(main = paste("Regression of ",sexType ," Test Set Version ", j, sep= ""))
  dev.off()
  
  
  allPredicted <- c(allPredicted, p)
  allActual <- c(allActual, testSet[,1])
  
  
} 
#####SOMEHOWCOMBINE ALL THE PREDICTIONS (of test set) INTO AN AUC TYPE CURVE THING#####*******(**&)

#Need to figure out how to turn the allPredicted and allActual into an AUC
# https://cran.r-project.org/web/packages/AUC/AUC.pdf
auc(sensitivity(allPredicted,allActual))

roc(allPredicted, allActual)
plot(roc(allPredicted, allActual))
auc(roc(allPredicted,allActual))

data(churn)
auc(sensitivity(churn$predictions,churn$labels))
plot(roc(sensitivity(churn$predictions,churn$labels)))
auc(specificity(churn$predictions,churn$labels))
auc(accuracy(churn$predictions,churn$labels))


data(churn)
plot(sensitivity(churn$predictions,churn$labels))
plot(specificity(churn$predictions,churn$labels))
plot(accuracy(churn$predictions,churn$labels))
plot(roc(churn$predictions,churn$labels))


png(filename=paste(savePath,"Combined",bMarker,sexType,"TestPlot.png", sep=""))
plot.new()
plot(allActual, allPredicted,lty=9,lwd=9,cex = 1,cex.axis=2,cex.lab = 2,col="red",
     xlim=c(min(predData[,1]),max(predData[,1])), ylim=c(min(predData[,1]),max(predData[,1])),
     xlab = paste("Actual Value of ", bMarker, sep=""), ylab = paste("Predicted Value of ", bMarker, sep=""))
title(main = paste("Combined Predictions from All KFolds for ",sexType," " , bMarker, sep=""))
dev.off()

rVal <- cor(allActual,allPredicted)
regResults[11,4] <- rVal
regResults[11,5] <- rVal*rVal
regResults[11,6] <- sqrt(mean((fitted.results - testSet[,1])^2, na.rm = TRUE) )

write.csv(regResults, file = paste(savePath,sexType,"_",bMarker,"_kFoldRegression.csv", sep=""),row.names=TRUE)

close(fileConn)


}

}