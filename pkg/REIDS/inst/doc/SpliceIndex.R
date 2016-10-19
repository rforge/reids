library(REIDS)

load("TissueData/TissueData_ExonLevelSummarized.RData")
load("TissueData/TissueData_GeneLevelSummarized.RData")


groupHMPT=c(7,8,9,16,17,18,22,23,24,31,32,33)
groupOthers=c(1,2,3,4,5,6,10,11,12,13,14,15,19,20,21,25,26,27,28,29,30)
groups=list(group1=groupHMPT,group2=groupOthers)

SI_Output=SpliceIndex(GeneData=TissueData_ExonLevelSummarized,ExonData=TissueData_GeneLevelSummarized,InformativeExons=NULL,groups=groups,paired=FALSE,significancelevel=NULL)
save(SI_Output,file="TissueData/SI_Output.RData")