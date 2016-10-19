library(REIDS)

DataProcessing(chipType="HuEx-1_0-st-v2",tags="coreR3,A20071112,EP",Name="TissueData",ExonSummarization=TRUE,GeneSummarization=TRUE,FIRMA=TRUE,location="TissueData",verbose=TRUE)

load("TissueData.RData")

PivotTransformData(Data=TissueData,GeneID=TissueData$GeneID,ExonID=TissueData$ExonID,savecsv=FALSE,Name="TissueData_Pivot",location="TissueData")