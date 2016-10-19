# Packages:
library(REIDS)

#Read in the .csv files of the Gene IDs
geneID=read.csv(file="GeneID.csv",header=TRUE)  

# Bind all output files of the REIDS function together
CreateOutput(geneID)