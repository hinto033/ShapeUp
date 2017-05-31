#Load Relevant Packages#####

#REQUIRED INPUT FOR EACH MACHINE#####

#All relevant variable that need changing are done in this section
dataDir <- "W:\\3D Optical_Data\\Fit 3D\\Michelle Shape Up Optical Models\\0515 Updated Shape Models\\7k template registration\\"
saveDir <- "X:\\bhinton\\ShapeUp\\PCA_To_BloodMarkersGender\\Analysis\\LinRegression\\"
nFolds <-10 #Number of K Folds cross validations to do

#LOOP for entire analysis#####
for (sexLoop in 1:2){
  #Loads Data for Males and Females#####
  setwd(dataDir)
  if (sexLoop==1){  
    allData <- read.csv('combined_measurements_Male.csv', header = TRUE)
    sexType <- "Male"
  }else if (sexLoop ==2){
    allData <- read.csv('combined_measurements_Female.csv', header = TRUE)
    sexType <- "Female"
  }
  #Isolate variables of interest#####
  
  print(colnames(allData))
  keep <- c("SubjectID", "Ethnicity", "BMI","WaistCircum1", "HipCircum","UpperArmCircum",
            "GLU", "CHOL", "TRIG", "HDL", "LDL","INS","AIC",
            "WB_TOT_BMC", "WB_TOT_BMD","WB_TOT_FAT","WB_TOT_Lean","WB_TOT_PFAT","WB_VFAT_AREA","WB_SAT_AREA",
            "age","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","HEIGHT","WEIGHT" )
  relData <- allData[keep]
  
  #Set Save Directory
  setwd(saveDir)
  savePathBase <- (saveDir)
  
  #LOOP for Blood Marker Selection#####
  #Perform analysis with Each Blood Marker
  bMarkerVec <- c("GLU", "CHOL", "TRIG", "HDL", "LDL","INS","AIC")
  for (bMarkLoop in 1:7){
    #Start Analysis & kFold Prep of Each Blood Marker ##### 
    
    bMarker <- bMarkerVec[bMarkLoop]
    savePath <- paste(savePathBase,bMarker, "\\", sep="")
    setwd(savePath)
    keep <- c(bMarker,'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7','PC8', 'PC9' )
    pcsKeep <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7','PC8', 'PC9' )
    predData <- na.omit(relData[keep])
    
    #Prep for K-Folds
    #Randomly shuffle the data
    predData<-predData[sample(nrow(predData)),]
    #Create nFolds equally size folds
    predDataFolds <- cut(seq(1,nrow(predData)),breaks=nFolds,labels=FALSE)
    
    #Initialize Regression Tables#####
    
    regResults <- data.frame(matrix(ncol = 6, nrow = nFolds+1))
    colnames(regResults) <- c('train_r','train_r2', 'train_rmse', 'test_r', 'test_r2', 'test_rmse')
    allPredicted <- NULL
    allActual <- NULL
    txtRegFile <-paste("Male_regression",bMarker,".txt", sep="")
    fileConn<-file(txtRegFile)
    #LOOP For KFolds Regression#####
    for(j in 1:nFolds){
      #Perform KFolds Regression#####
      
      #Set train and test set for each fold
      testSetIndexes <- which(predDataFolds==j,arr.ind=TRUE)
      testSet <- rbind(predData[testSetIndexes, ])
      trainSet <- rbind(predData[-testSetIndexes, ])
      
      #Define model
      model <- lm(trainSet[,1]~., data=trainSet[pcsKeep])
      modelStep <- step(model)
      trainRegSumm <- summary(modelStep)
      out <- capture.output(summary(modelStep))
      cat(out,file=txtRegFile,sep="\n",append=TRUE)
      
      #Plot training regression and save
      png(filename=paste(savePath,bMarker,"TrainPlot",j,".png", sep=""))
      plot.new()
      plot(trainSet[,1], modelStep$fitted.values,lty=9,lwd=9,cex = 1,cex.axis=2,cex.lab = 2,col="red",
           xlim=c(min(predData[,1]),max(predData[,1])), ylim=c(min(predData[,1]),max(predData[,1])),
           xlab = paste("Actual Value of ", bMarker, sep=""), ylab = paste("Predicted Value of ", bMarker, sep=""))
      title(main = paste("Regression of ",sexType ," Train Set Version ", j, sep= ""))
      dev.off()
      
      #Save statistics of regression
      rValTrain <- cor(modelStep$fitted.values,trainSet[,1])
      regResults[j,1] <- rValTrain
      regResults[j,2] <- rValTrain*rValTrain
      regResults[j,3] <- sqrt(mean((modelStep$fitted.values - trainSet[,1])^2, na.rm = TRUE) )
      
      #Perform analysis on test set.
      fitted.results <- predict(modelStep,newdata=testSet,type='response')
      allPredicted <- c(allPredicted, fitted.results)
      allActual <- c(allActual, testSet[,1])
      
      #plot and save test set regression
      png(filename=paste(savePath,bMarker,sexType,"TestPlot",j,".png", sep=""))
      plot.new()
      plot(testSet[,1], fitted.results,lty=9,lwd=9,cex = 1,cex.axis=2,cex.lab = 2,col="red",
           xlim=c(min(predData[,1]),max(predData[,1])), ylim=c(min(predData[,1]),max(predData[,1])),
           xlab = paste("Actual Value of ", bMarker, sep=""), ylab = paste("Predicted Value of ", bMarker, sep=""))
      title(main = paste("Regression of ",sexType ," Test Set Version ", j, sep= ""))
      dev.off()
      
      #Save statistics of test set
      rVal <- cor(fitted.results,testSet[,1])
      regResults[j,4] <- rVal
      regResults[j,5] <- rVal*rVal
      regResults[j,6] <- sqrt(mean((fitted.results - testSet[,1])^2, na.rm = TRUE) )
    } 
    
    #GeneratePlot for Combining kFolds and save#####
    
    #plot and Save the combined k-folds regression
    png(filename=paste(savePath,"Combined",bMarker,sexType,"TestPlot.png", sep=""))
    plot.new()
    plot(allActual, allPredicted,lty=9,lwd=9,cex = 1,cex.axis=2,cex.lab = 2,col="red",
         xlim=c(min(predData[,1]),max(predData[,1])), ylim=c(min(predData[,1]),max(predData[,1])),
         xlab = paste("Actual Value of ", bMarker, sep=""), ylab = paste("Predicted Value of ", bMarker, sep=""))
    title(main = paste("Combined Predictions from All KFolds for ",sexType," " , bMarker, sep=""))
    dev.off()
    
    #Save R, r^2, RMSE values
    rVal <- cor(allActual,allPredicted)
    regResults[11,4] <- rVal
    regResults[11,5] <- rVal*rVal
    regResults[11,6] <- sqrt(mean((fitted.results - testSet[,1])^2, na.rm = TRUE) )
    
    #Save tables
    write.csv(regResults, file = paste(savePath,sexType,"_",bMarker,"_kFoldRegression.csv", sep=""),row.names=TRUE)
    close(fileConn)
  }
}