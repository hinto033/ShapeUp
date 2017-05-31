#Load Necessary Packages#####

install.packages("ROCR")
require("ROCR")
install.packages("AUC")
require("AUC")

#REQUIRED INPUT FOR EACH MACHINE#####

#All relevant variable that need changing are done in this section
dataDir <- "W:\\3D Optical_Data\\Fit 3D\\Michelle Shape Up Optical Models\\0515 Updated Shape Models\\7k template registration\\"
saveDir <- "X:\\bhinton\\ShapeUp\\PCA_To_BloodMarkersGender\\Analysis\\Logistic\\"
nFolds <-3 #Number of K Folds cross validations to do

#LOOP for entire analysis#####
for (sexLoop in 1:2){
  #Does loop for males and females#####
  setwd(dataDir)
  
  #Loads appropriate document
  if (sexLoop==1){  
    allData <- read.csv('combined_measurements_Male.csv', header = TRUE)
    sexType <- "Male"
  }else if (sexLoop ==2){
    allData <- read.csv('combined_measurements_Female.csv', header = TRUE)
    sexType <- "Female"
  }
    
  #Add Safe/Unhealthy thresholds for each Marker#####
  # 0 is good, 1 is bad
  
  # Glucose
  allData$GLUThresh[allData$GLU<=100]<-0  #Healthy
  allData$GLUThresh[allData$GLU>100]<-1  #Prediabetic
  # allData$GLUThresh[allData$GLU>126]<-2  #Diabetes
  
  # Cholesterol
  allData$CHOLThresh[allData$CHOL<=200]<-0 #Good
  allData$CHOLThresh[allData$CHOL>200]<-1 #Bad
  # allData$CHOLThresh[allData$CHOL>240]<-2
  
  # Triglycerides
  allData$TRIGThresh[allData$TRIG<=150]<-0  #Healthy
  allData$TRIGThresh[allData$TRIG>150]<-1   #Moderate Risk
  # allData$TRIGThresh[allData$TRIG>200]<-2   #High Heart Disease Risk
  
  # HDL
  allData$HDLThresh[allData$HDL<40]<-1   #High risk Heart Disease
  allData$HDLThresh[allData$HDL>=40]<-0   #Healthy
  
  # LDL
  allData$LDLThresh[allData$LDL>130]<-0 #Good
  allData$LDLThresh[allData$LDL<=130]<-1 #Bad
  # allData$LDLThresh[allData$LDL<100]<-2
  
  # INS
  allData$INSThresh[allData$INS>27]<-1  #Insulin Resistance
  allData$INSThresh[allData$INS<=27]<-0  #Healthy
  # allData$INSThresh[allData$INS<7]<-0  #Diabetes
  
  # AIC
  allData$AICThresh[allData$AIC<=5.7]<-0  #Healthy
  allData$AICThresh[allData$AIC>5.7]<-1  #Prediabetes
  # allData$AICThresh[allData$AIC>6.3]<-2  #Diabetes
  
  #Isolate variables of interest#####
  print(colnames(allData))
  
  #Saving relevant data for logistic regression
  keep <- c("SubjectID", "Ethnicity", "BMI","WaistCircum1", "HipCircum","UpperArmCircum",
            "GLUThresh", "CHOLThresh", "TRIGThresh", "HDLThresh", "LDLThresh","INSThresh","AICThresh",
            "WB_TOT_BMC", "WB_TOT_BMD","WB_TOT_FAT","WB_TOT_Lean","WB_TOT_PFAT","WB_VFAT_AREA","WB_SAT_AREA",
            "age","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","HEIGHT","WEIGHT" )
  relData <- allData[keep]
  
  #Set up vectors to run loops for analysis of each blood marker
  bMarkerVec <- c("GLU", "CHOL", "TRIG", "HDL", "LDL","INS","AIC")
  bMarkerThresh <- c("GLUThresh", "CHOLThresh", "TRIGThresh", "HDLThresh", "LDLThresh","INSThresh","AICThresh")

  #Set Save Directory for once analysis begins
  setwd(saveDir)
  savePathBase <- (saveDir)
  
  #LOOP for Blood Marker Selection#####
  for (bMarkLoop in 1:7){
    #Start Analysis of Each Blood Marker #####  
    
    #Set savepath to right folder
    bMarker <- bMarkerVec[bMarkLoop]
    savePath <- paste(savePathBase,bMarker, "\\", sep="")
    setwd(savePath)
    
    #Save just the PCs and the one blood marker threshold
    keep <- c(bMarkerThresh[bMarkLoop],'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7','PC8', 'PC9' )
    pcsKeep <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7','PC8', 'PC9' )
    predData <- na.omit(relData[keep])
    #Check if enough data for logistic regression#####
    
    #Searches if possible to have at least 1 extreme case in each fold
    if (sum(predData[,1]) < nFolds){
      txtRegFile <-paste ("ErrorNotEnough",sexType,"Cases.txt", sep="")
      fileConn<-file(txtRegFile)
      cat(paste ("NOT ENOUGH CASES FOR ", sexType, sep=""),file=txtRegFile,sep="\n",append=TRUE)
      close(fileConn)
      next
    }
    
    #Reshuffles and reorders kfolds until 1 extreme marker in each fold
    goodSet <- 0
    counter <-0
    while(goodSet < 1){
      goodSet <-1
      
      #Randomly shuffle the data
      predData<-predData[sample(nrow(predData)),]
      #Create nFolds equally size folds
      predDataFolds <- cut(seq(1,nrow(predData)),breaks=nFolds,labels=FALSE)
      dataCheck <- data.frame(cbind(predData[,1],predDataFolds))
      for (p in 1:3){ #If any fold has no extreme cases, then reshuffles
        if (sum(dataCheck[which(dataCheck$predDataFolds == p),1]) == 0){
          goodSet <- 0
        }
      }
      
      #If shuffled over 50 times with no good result, then produce error document
      counter<- counter+1
      if (counter > 50){ 
        goodSet <- 0
        bMarkLoop <- bMarkLoop+1
        txtRegFile <-paste ("ErrorNotEnough",sexType,"Cases.txt", sep="")
        fileConn<-file(txtRegFile)
        cat(paste ("NOT ENOUGH CASES FOR ", sexType, sep=""),file=txtRegFile,sep="\n",append=TRUE)
        close(fileConn)
        
        #Then moves on to the next blood marker
        bMarker <- bMarkerVec[bMarkLoop]
        savePath <- paste(savePathBase,bMarker, "\\", sep="")
        setwd(savePath)
        keep <- c(bMarkerThresh[bMarkLoop],'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7','PC8', 'PC9' )
        pcsKeep <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7','PC8', 'PC9' )
        predData <- na.omit(relData[keep])
      }
    }
    #Initialize Regression Tables#####
    
    #Create, name columns, etc.
    regResults <- data.frame(matrix(ncol = 4, nrow = nFolds+1))
    colnames(regResults) <- c('train_pctCorrect','train_AUC', 'test_pctCorrect', 'test_auc')
    allPredicted <- NULL
    allActual <- NULL
    txtRegFile <-paste(sexType,"_logisticRegression_",bMarker,".txt", sep="")
    fileConn<-file(txtRegFile)
    
    #LOOP for k-folds analysis of that blood marker#####
      for(j in 1:nFolds){
        #Perform KFolds Regression#####
        
        #Set test and train set for each fold
        testSetIndexes <- which(predDataFolds==j,arr.ind=TRUE)
        testSet <- rbind(predData[testSetIndexes, ])
        trainSet <- rbind(predData[-testSetIndexes, ])
        
        #Perform the model & stepwise regerssion
        model <- lm(trainSet[,1]~., family=binomial(logit), data=trainSet[pcsKeep])
        modelStep <- step(model)
        trainRegSumm <- summary(modelStep)
        out <- capture.output(summary(modelStep)) #Saves to .txt.
        cat(out,file=txtRegFile,sep="\n",append=TRUE)
        
        #Plots AUC of training set and saves it.
        png(filename=paste(savePath,bMarker,sexType,"TrainAUCPlot",j,".png", sep=""))
        plot.new()
        p <- predict(model, newdata=trainSet[pcsKeep], type="response")
        pr <- prediction(p, trainSet[,1])
        prf <- performance(pr, measure = "tpr", x.measure = "fpr")
        plot(prf, lty=1,lwd=5,cex = 1,cex.axis=2,cex.lab = 2,col="Blue")
        title(main = paste("ROC AUC of ",sexType ," Train Set Version ", j, sep= ""))
        dev.off()
        
        #Saves training set statistics (Accuracy, ROC)
        forced.results <- ifelse(p > 0.5,1,0)
        misClasificError <- mean(forced.results != trainSet[,1])
        Accuracy<- 1-misClasificError
        auc <- performance(pr, measure = "auc")
        aucMarkers <- auc@y.values[1]
        regResults[j,1] <- Accuracy 
        regResults[j,2] <- aucMarkers
      
        #Applies the model to the test set
        fitted.results <- predict(modelStep,newdata=testSet,type='response')
        forced.results <- ifelse(fitted.results > 0.5,1,0)
        misClasificError <- mean(forced.results != testSet[,1])
        Accuracy<- 1-misClasificError
        
        #Saves a plot of the test set
        png(filename=paste(savePath,bMarker,sexType,"TestAUCPlot",j,".png", sep=""))
        plot.new()
        p <- predict(model, newdata=testSet[pcsKeep], type="response")
        pr <- prediction(p, testSet[,1])
        prf <- performance(pr, measure = "tpr", x.measure = "fpr")
        plot(prf, lty=1,lwd=5,cex = 1,cex.axis=2,cex.lab = 2,col="Blue")
        title(main = paste("ROC AUC of ",sexType ," Test Set Version ", j, sep= ""))
        dev.off()
        
        #Adding to table for saving later
        auc <- performance(pr, measure = "auc")
        aucMarkers <- auc@y.values[1]
        regResults[j,3] <- Accuracy
        regResults[j,4] <- aucMarkers
      
        #Creates running list of all predictions/actual results for combined k-folds later
        allPredicted <- c(allPredicted, p)
        allActual <- c(allActual, testSet[,1])
        
      } 
    #Combining results from all kFolds#####
    
    allForcedChoice <- ifelse(allPredicted > 0.5,1,0)
    misClasificError <- mean(allForcedChoice != allActual)
    print(paste('Accuracy',1-misClasificError))
    Accuracy<- 1-misClasificError
    
    #Predicts with the combined data
    predictDF <- data.frame(cbind(allPredicted,allForcedChoice,allActual))
    pred <- prediction(predictDF$allPredicted,predictDF$allActual)
    prf <- performance(pred, measure = "tpr", x.measure = "fpr")
    auc <- performance(pred, measure = "auc")
    aucMarkers <- auc@y.values[1]
    
    #Plots and saves combined prediction
    png(filename=paste(savePath,"AUCCombinedKfold_",bMarker,sexType,".png", sep=""))
    plot.new()
    plot(prf,lty=1,lwd=9,cex = 1,cex.axis=2,cex.lab = 2,col="Blue")#
    title(main = paste("AUC Curve from All KFolds for ",sexType," " , bMarker, sep=""))
    dev.off()
    
    #Averages training set statistics
    regResults[nFolds+1,1] <- mean(na.omit(regResults[,1]))
    regResults[nFolds+1,2] <- mean(na.omit(regResults[,2]))
    
    #adds the combined k-folds accuracy/AUC as test set avg
    regResults[nFolds+1,3] <- Accuracy
    regResults[nFolds+1,4] <- aucMarkers
    
    write.csv(regResults, file = paste(savePath,sexType,"_",bMarker,"_kFoldRegression.csv", sep=""),row.names=TRUE)
    close(fileConn)
  }
}