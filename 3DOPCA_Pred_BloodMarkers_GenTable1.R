#Load Necessary Packages#####

install.packages("ROCR")
require("ROCR")
install.packages("AUC")
require("AUC")

require('fBasics')
library(matrixStats)

library("psych")    

#REQUIRED INPUT FOR EACH MACHINE#####

#All relevant variable that need changing are done in this section
dataDir <- "W:\\3D Optical_Data\\Fit 3D\\Michelle Shape Up Optical Models\\0515 Combined Updated Shape Models\\7k template registration\\"
saveDir <- "X:\\bhinton\\ShapeUp\\PCA_To_BloodMarkersGender\\Analysis\\"
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
  
  
  
  #Generate Table 1: Demographics #####
  keep <- c('BMI',  'age', 'HEIGHT', 'WEIGHT', 'GLU', 'CHOL','TRIG','HDL', 'LDL', 'INS', 'AIC',
            'GLUThresh', 'CHOLThresh', 'TRIGThresh', 'HDLThresh', 'LDLThresh', 'INSThresh', 'AICThresh')
   relData<-(allData[keep])

  demogMedians<-colMeans(relData, na.rm=TRUE)
  demogsds<-colStdevs(relData, na.rm=TRUE)
  demogSums<- colSums(relData, na.rm=TRUE)  #Sum of the threshold ones is number of positive cases.

  demogTable<- rbind(demogMedians, demogsds, demogSums)  

  nEth <- length(list(table(allData$Ethnicity))[[1]])
  ethTable <- data.frame(matrix(ncol = 2, nrow = nEth))
  
  #Generate Table 2: Correlations #####
  
  pcsKeep<-c('PC1',  'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7','PC8', 'PC9')
  corDataKeep<-c('BMI',  'age', 'HEIGHT', 'WEIGHT', 'GLU', 'CHOL','TRIG','HDL', 'LDL', 'INS', 'AIC',
                 'GLUThresh', 'CHOLThresh', 'TRIGThresh', 'HDLThresh', 'LDLThresh', 'INSThresh', 'AICThresh',
                 "WaistCircum1", "HipCircum","UpperArmCircum","WB_TOT_BMC", "WB_TOT_BMD","WB_TOT_FAT",
                 "WB_TOT_Lean","WB_TOT_PFAT","WB_VFAT_AREA","WB_SAT_AREA")
  
  pcsTable <- allData[pcsKeep]
  corDataTable<- allData[corDataKeep]
  corData<-cor(pcsTable,corDataTable, use = "complete.obs")
  
  
  relDataKeep<-c('PC1',  'PC2', 'PC3', 'PC4', 'PC5', 'PC6','PC7','PC8', 'PC9',
                 'BMI',  'age', 'HEIGHT', 'WEIGHT', 'GLU', 'CHOL','TRIG','HDL', 'LDL', 'INS', 'AIC',
                 'GLUThresh', 'CHOLThresh', 'TRIGThresh', 'HDLThresh', 'LDLThresh', 'INSThresh', 'AICThresh',
                 "WaistCircum1", "HipCircum","UpperArmCircum","WB_TOT_BMC", "WB_TOT_BMD","WB_TOT_FAT",
                 "WB_TOT_Lean","WB_TOT_PFAT","WB_VFAT_AREA","WB_SAT_AREA")
  relDataTable<-allData[relDataKeep]
  
  corpVals<- corr.test(relDataTable, use = "complete.obs", adjust = "none")
  rVals<- corpVals$r
  pVals<- corpVals$p
  
  # Save ethTable and demogTable and rVal Table and pVal Table
  write.csv(demogTable, file = paste(saveDir,sexType,"_DemographicSummaryTable.csv", sep=""),row.names=TRUE)
  write.csv(ethTable, file = paste(saveDir,sexType,"_EthnicitySummaryTable.csv", sep=""),row.names=TRUE)
  write.csv(rVals, file = paste(saveDir,sexType,"_CorrelationRTable.csv", sep=""),row.names=TRUE)
  write.csv(pVals, file = paste(saveDir,sexType,"_CorrelationPValues.csv", sep=""),row.names=TRUE)
}