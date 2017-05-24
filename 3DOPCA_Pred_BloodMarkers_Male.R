#####
#Load Relevant Packages


#####
#Load Data
#Different Possible Paths based on Computer

setwd("W:\\3D Optical_Data\\Fit 3D\\Michelle Shape Up Optical Models\\0515 Updated Shape Models\\")
allData <- read.csv('combined_measurements_Male.csv', header = TRUE)


#####
#Description of Data
print(colnames(allData))

keep <- c("SubjectID", "Ethnicity", "BMI","WaistCircum1", "HipCircum","UpperArmCircum",
          "GLU", "CHOL", "TRIG", "HDL", "LDL","INS","AIC",
          "WB_TOT_BMC", "WB_TOT_BMD","WB_TOT_FAT","WB_TOT_Lean","WB_TOT_PFAT","WB_VFAT_AREA","WB_SAT_AREA",
          "age","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","HEIGHT","WEIGHT" )
relData <- allData[keep]

#####
#Regression w.r.t. Only Shape Markers

# setwd("D:\\LabData\\ShapeUp\\Shapeup-PCAToBloodMarkers\\")
setwd("X:\\bhinton\\ShapeUp\\PCA_To_BloodMarkersGender\\Code\\")
# fileConn<-file("3DO_regressionWRTShapePCsOnly.txt")


#####



#####

# "GLU", "CHOL", "TRIG", "HDL", "LDL","INS","AIC",
keep <- c('INS','PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7','PC8', 'PC9' )
predData <- na.omit(relData[keep])
#Prep for K-Folds
#Randomly shuffle the data (CC)
predData<-predData[sample(nrow(predData)),]
#Create 10 equally size folds (CC)
predDataFolds <- cut(seq(1,nrow(predData)),breaks=10,labels=FALSE)

regResults <- data.frame(matrix(ncol = 2, nrow = 10))
colnames(regResults) <- c('r','r2')
# j<-1

for(j in 1:10){
  testSetIndexes <- which(relDataFolds==j,arr.ind=TRUE)
  testSet <- rbind(predData[testSetIndexes, ])
  trainSet <- rbind(predData[-testSetIndexes, ])

  model <- lm(trainSet$INS~., data=trainSet)
  modelStep <- step(model, direction='forward')
  modelStep <- step(model)
  summary(modelStep)
  
  fitted.results <- predict(modelStep,newdata=testSet,type='response')
  
  # Add Plots of actual vs predicted? ################################################
  # Save the Actual Regression Equations ###########################################
  # Save the training R2 and the test R2 #############################################
  
  
  rVal <- cor(fitted.results,testSet$INS)
  regResults[j,1] <- rVal
  regResults[j,2] <- rVal*rVal
} 







#####









keep <- c('HDL','PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$HDL~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionWRTShapePCsOnly.txt",sep="\n",append=TRUE)

keep <- c('LDL','PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$LDL~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionWRTShapePCsOnly.txt",sep="\n",append=TRUE)

keep <- c('TRIG','PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$TRIG~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionWRTShapePCsOnly.txt",sep="\n",append=TRUE)

keep <- c('CHOL','PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$CHOL~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionWRTShapePCsOnly.txt",sep="\n",append=TRUE)

keep <- c('GLU','PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$GLU~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionWRTShapePCsOnly.txt",sep="\n",append=TRUE)

keep <- c('INS','PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$INS~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionWRTShapePCsOnly.txt",sep="\n",append=TRUE)

keep <- c('URIC','PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$URIC~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionWRTShapePCsOnly.txt",sep="\n",append=TRUE)

close(fileConn)
#####


#Regression w.r.t. Shape After Controls

# setwd("X:\\bhinton\\Shapeup-PCAToBloodMarkers\\")
setwd("D:\\LabData\\ShapeUp\\Shapeup-PCAToBloodMarkers\\")
fileConn<-file("3DO_regressionShapeAfterDemogControl.txt")

keep <- c('HBAIC', 'age' ,'Gender', 'Ethnicity', 'BMI', 'HeightCMAvg', 'WeightKGAvg',
          'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7' )
predData <- na.omit(relData[keep])
model <- lm(predData$HBAIC~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionShapeAfterDemogControl.txt",sep="\n",append=TRUE)

keep <- c('HDL','age' ,'Gender', 'Ethnicity', 'BMI', 'HeightCMAvg', 'WeightKGAvg',
          'Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$HDL~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionShapeAfterDemogControl.txt",sep="\n",append=TRUE)

keep <- c('LDL','age' ,'Gender', 'Ethnicity', 'BMI', 'HeightCMAvg', 'WeightKGAvg',
          'Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$LDL~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionShapeAfterDemogControl.txt",sep="\n",append=TRUE)

keep <- c('TRIG','age' ,'Gender', 'Ethnicity', 'BMI', 'HeightCMAvg', 'WeightKGAvg',
          'Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$TRIG~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionShapeAfterDemogControl.txt",sep="\n",append=TRUE)

keep <- c('CHOL','age' ,'Gender', 'Ethnicity', 'BMI', 'HeightCMAvg', 'WeightKGAvg',
          'Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$CHOL~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionShapeAfterDemogControl.txt",sep="\n",append=TRUE)

keep <- c('GLU','age' ,'Gender', 'Ethnicity', 'BMI', 'HeightCMAvg', 'WeightKGAvg',
          'Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$GLU~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionShapeAfterDemogControl.txt",sep="\n",append=TRUE)

keep <- c('INS', 'age' ,'Gender', 'Ethnicity', 'BMI', 'HeightCMAvg', 'WeightKGAvg',
          'Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$INS~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionShapeAfterDemogControl.txt",sep="\n",append=TRUE)

keep <- c('URIC', 'age' ,'Gender', 'Ethnicity', 'BMI', 'HeightCMAvg', 'WeightKGAvg',
          'Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$URIC~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionShapeAfterDemogControl.txt",sep="\n",append=TRUE)


close(fileConn)


