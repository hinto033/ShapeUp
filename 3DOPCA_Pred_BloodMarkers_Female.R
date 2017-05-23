#####
#Load Relevant Packages
install.packages("sas7bdat")
require("sas7bdat")

#####
#Load Data

#Load Demog/DXA/Blood
# setwd("X:\\bhinton\\Shapeup-PCAToBloodMarkers\\")
setwd("D:\\LabData\\ShapeUp\\Shapeup-PCAToBloodMarkers\\")

allData <- read.sas7bdat('combine_all.sas7bdat', debug=FALSE)


#####
#Description of Data
print(colnames(allData))

keep <- c('HBAIC','HDL','LDL','TRIG','CHOL', 'GLU', 'INS', 'URIC', 
          'Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7',
          'Gender', 'Ethnicity', 'HeightCMAvg', 'WeightKGAvg', 'BMI', 'age',
          'TOTAL_PFAT', 'TOTAL_MASS', 'RestingHR60',
          'WB_TOT_PFAT', 'WB_TOT_FAT', 'WB_TOT_Lean')
relData <- allData[keep]

#####
#Regression w.r.t. Only Shape Markers

# setwd("X:\\bhinton\\Shapeup-PCAToBloodMarkers\\")
setwd("D:\\LabData\\ShapeUp\\Shapeup-PCAToBloodMarkers\\")
fileConn<-file("3DO_regressionWRTShapePCsOnly.txt")

keep <- c('HBAIC','Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$HBAIC~., data=predData)
modelStep <- step(model, direction='forward')
modelStep <- step(model)
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionWRTShapePCsOnly.txt",sep="\n",append=TRUE)

keep <- c('HDL','Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$HDL~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionWRTShapePCsOnly.txt",sep="\n",append=TRUE)

keep <- c('LDL','Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$LDL~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionWRTShapePCsOnly.txt",sep="\n",append=TRUE)

keep <- c('TRIG','Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$TRIG~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionWRTShapePCsOnly.txt",sep="\n",append=TRUE)

keep <- c('CHOL','Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$CHOL~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionWRTShapePCsOnly.txt",sep="\n",append=TRUE)

keep <- c('GLU','Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$GLU~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionWRTShapePCsOnly.txt",sep="\n",append=TRUE)

keep <- c('INS','Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
predData <- na.omit(relData[keep])
model <- lm(predData$INS~., data=predData)
modelStep <- step(model)
rmseModel <- sqrt(mean((modelStep$fitted.values - modelStep$y)^2, na.rm = TRUE) )
r2Model <- summary(modelStep)$r.squared
out <- capture.output(summary(modelStep))
cat(out,file="3DO_regressionWRTShapePCsOnly.txt",sep="\n",append=TRUE)

keep <- c('URIC','Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7')
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
          'Op_PC1', 'Op_PC2', 'Op_PC3', 'Op_PC4', 'Op_PC5', 'Op_PC6','Op_PC7' )
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


