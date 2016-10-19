## THE REIDS PACKAGE ##


## IMPORTS ##

#' @import aroma.affymetrix
#' @import aroma.core
#' @importFrom MCMCpack riwish
#' @importFrom data.table rbindlist 
#' @importFrom RColorBrewer brewer.pal


##### DATA PROCESSING FUNCTIONS #####


## DataProcessing ##
## A wrapper function of some of the functions in the aroma.affymetrix package to make sure the data is processed correctly



#' @title DataProcessing
#' 
#' @description The DataProcessing function processes raw .CEL files to probe intensities values with the help of functions of the aroma.affymetrix package. It returns a data frame and saves it as an .RData file.
#' 
#' @export
#' @param chipType The name of the chip type of the array data.
#' @param tags Tags that is added to the chipType.
#' @param Name The name of the data.
#' @param ExonSummarization Logical. Should the data be summarized at the exon level?
#' @param GeneSummarization Logical. Should the data be summarized at the gene level?
#' @param FIRMA Logical. Should the FIRMA model be performed on the data?
#' @param location The location where the .rda file is to be stored. Defaults to the current working directory.
#' @param verbose Logical. If TRUE, messages are printed during the data processing.
#' @return An .rda file that is saved at the specified location.
#' @details The DataProcessing function is a wrapper of several functions of the aroma.affymetrix package. To obtain the data to perform the REIDS model on the raw .CEL files are 
#' background corrected with the rma background correction and normalization is performed with the quantile normalization. In order for the function to run properly, a chipType and
#' its possible tags need to be specified. It is also important to have the same folder structure as required by the aroma.affymetrix package. This implies the following:
#' a rawData folder with therein a folder with the "Name" parameter. This "Name" folder should contain a folder with the chipType name and herein the .CEL files should be placed.
#' Also a folder annotationData should be present. Herein a folder chipTypes should be make which contains folders for type of chips with the respective names. In the folder of
#' each chiptype the corresponding .cdf file should be saved. The processed data will be at the specified location save as data frame with the first colum the gene IDs and the second column the exon IDs. All
#' other columns contain the sample values. Further the object also contains a vector of the unique gene ID and a vector of the  unique exon IDs. If requested, exon and gene level summarization are performed and saved as data frames at the specified location. Further,the option is provided to perform the FIRMA model on the data as well.
#' @examples
#' \dontrun{
#' DataProcessing(chipType="HuEx-1_0-st-v2",tags="coreR3,A20071112,EP",
#' Name="ColonCancer",ExonSummarization=TRUE,GeneSummarization=TRUE,
#' FIRMA=TRUE,location="ColonCancer",verbose=TRUE)
#' }
DataProcessing<-function(chipType="HuEx-1_0-st-v2",tags="coreR3,A20071112,EP",Name="ColonCancer",ExonSummarization=TRUE,GeneSummarization=TRUE,FIRMA=TRUE,location="test",verbose=TRUE){
	
	if(verbose==TRUE){
		verbose <- Arguments$getVerbose(-8, timestamp=TRUE)
	}
	
	#Setting up the CDF file
	if(is.null(tags)){
		cdf <- AffymetrixCdfFile$byChipType(chipType)
	}
	else{	
		cdf <- AffymetrixCdfFile$byChipType(chipType, tags=tags)
	}
	print(cdf)
	
	#Setting up CELset
	cs <- AffymetrixCelSet$byName(Name, cdf=cdf)
	print(cs)
	
	setCdf(cs, cdf)
	
	
	#Backgroundcorrection
	bc <- RmaBackgroundCorrection(cs,tags="*,coreR3")
	csBC <-process(bc,verbose=verbose)
	
	
	#(Quantile) Normalization
	qn <- QuantileNormalization(csBC, typesToUpdate="pm")
	csN <-process(qn, verbose=verbose)
	
	
	flattenCellIndices <- function(cells, ...) {
		# Flatten cell data
		cells <- unlist(cells);
		
		# Do some tricks to clean up the names
		names(cells) <- gsub("[.](groups|indices)", "", names(cells));
		
		# Extract the vector of unit names
		unitNames <- gsub("[.].*", "", names(cells))
		
		# Extract the unit group names
		groupNames <- gsub(".*[.]", "", names(cells))
		
		# Merge data
		data.frame(unitNames=unitNames, groupNames=groupNames, cell=cells)
	} 
	
	
	cells1 <-getCellIndices(cdf, units=1:nbrOfUnits(cdf))
	cells2 <-flattenCellIndices(cells1)

	probeintensities <- getData(csN, indices=cells2$cell, fields=c("intensities"))$intensities
	probeintensities <- log2(probeintensities)
	colnames(probeintensities) <- cs$names
	
	probeintensities=as.data.frame(probeintensities)
	
	if(length(grep("hta",chipType))!=0){
		
		ExonID <- as.character(cells2$unitNames)
		
		UniqueExonID <- unique(ExonID)
		
		Data=cbind(ExonID,probeintensities)
		Data$ExonID=as.character(Data$ExonID)
				
	}
	
	else{ 
		
		GeneID <- as.character(cells2$unitNames)
		ExonID <- as.character(substring(cells2$groupNames,1,7))
		
		UniqueGeneID <- unique(GeneID)
		UniqueExonID <- unique(ExonID)
		
		Data=cbind(GeneID,ExonID,probeintensities)
		Data$GeneID=as.character(Data$GeneID)
		Data$ExonID=as.character(Data$ExonID)
		
	}
	
	assign(Name,Data,envir=environment())
	if(is.null(location)){
		location=getwd()
	}
	eval(parse(text=paste("save(",Name,",UniqueGeneID,UniqueExonID,file=\"",location,"/", Name, ".RData\")", sep=""))) 	
	
	
	if(ExonSummarization==TRUE){
		
		getCdf(csN)
		plmEx <- ExonRmaPlm(csN, mergeGroups=FALSE)
		fit(plmEx, verbose=verbose)
		
		cesEx <- getChipEffectSet(plmEx)
		exFit <- extractDataFrame(cesEx, units=NULL, addNames=TRUE)
		
		ExonLevelSummarized_rma=exFit[,-c(3,4,5)]
		colnames(ExonLevelSummarized_rma)[c(1,2)]=c("GeneID","ExonID")
		
		assign(paste(Name,"ExonLevelSummarized",sep="_"),ExonLevelSummarized_rma,envir=environment())
		eval(parse(text=paste("save(",Name, "_ExonLevelSummarized, file=\"",location,"/",Name,"", "_ExonLevelSummarized.RData\")", sep=""))) 
		
		message("Warning: Data at exon level is not log2 transformed.")
		
	}
	
	if(GeneSummarization==TRUE){
		
		plmTr <- ExonRmaPlm(csN, mergeGroups=TRUE)
		fit(plmTr, verbose=verbose)
		
		cesTr <- getChipEffectSet(plmTr)
		TrFit <- extractDataFrame(cesTr, units=NULL, addNames=TRUE)
		
		GeneLevelSummarized_rma=TrFit[,-c(2,3,4,5)]
		colnames(GeneLevelSummarized_rma)[1]=c("GeneID")
		
		
		assign(paste(Name,"GeneLevelSummarized",sep="_"),GeneLevelSummarized_rma,envir=environment())
		eval(parse(text=paste("save(",Name, "_GeneLevelSummarized, file=\"",location,"/",Name,"", "_GeneLevelSummarized.RData\")", sep=""))) 
		
		
		message("Warning: Data at gene level is not log2 transformed.")
		
	}
	
	if(FIRMA==TRUE){
		plmTr <- ExonRmaPlm(csN, mergeGroups=TRUE)
		fit(plmTr,verbose=verbose)

		firma <- FirmaModel(plmTr)
		fit(firma, verbose=verbose)
		
		fs <- getFirmaScores(firma)
		exFit <- extractDataFrame(fs, units=NULL, addNames=TRUE)
		exFit=exFit[,-c(3:5)]
		colnames(exFit)[c(1,2)]=c("GeneID","ExonID")

		FIRMA_Output=exFit
		
		assign(paste(Name,"FIRMA_Output",sep="_"),FIRMA_Output,envir=environment())
		eval(parse(text=paste("save(",Name, "_FIRMA_Output, file=\"",location,"/",Name,"", "_FIRMA_Output.RData\")", sep=""))) 

	}
}


## Pivot Transformation
## Data format should be an data.frame with an GeneID colum, ExonID column and a column per array ##

#' "PivotTransformation"
#' 
#' The PivotTransformation function converts a data frame with multiple rows per gene into a .csv file with one row per gene. This is the first step in data transformation to apply the REIDS function on a HPC Cluster.
#' @export
#' @param Data The data frame to be transformed.
#' @param GeneID A character vector of the the gene IDs that correspond to the rows of the data frame. Necessary if no GeneID column is present in the data frame
#' @param ExonID A character vector of the the gene IDs that correspond to the rows of the data frame. Necessary if no ExonID column is present in the data frame
#' @param savecsv Logical. Should the file be saved as a .csv?
#' @param Name The name of the file if it is saved a a .csv file.
#' @param location The location where the file should be saved. The value NULL default to the current working directory.
#' @return A data.frame with one row per gene. This row contains the values for each exon per sample and is convenient for processing on a HPC cluster.
#' @details All information concerning one gene is gathered. The first column of the returned data frame is the gene ID, the second column contains the exon IDs of all exons of that gene. The third colum indicates the number of probes per exon, the fourth contains the values of thos probes per sample and the last column contains the sample names.This way a .csv file is created for processing on a HPC cluster.
#' @examples
#' \dontrun{
#' data=data.frame(GeneID=c("1","1","1","2","2","3","3","3","3"),
#' ExonID=c("a","b","b","c","c","d","d","e","e"),
#' Sample1=c(1.2,3.5,6.2,8.5,9.6,7.8,2.2,2.1,6.3),
#' Sample2=c(1.5,1.2,1.3,2.9,2.5,2.4,7.6,7.2,7.3))
#' 
#' PivotTransformData(Data=data,GeneID=NULL,ExonID=NULL,savecsv=TRUE,
#' Name="test",location=NULL)
#' }
PivotTransformData<-function(Data, GeneID=NULL,ExonID=NULL,savecsv=FALSE,Name=NULL,location=NULL){
	Data=as.data.frame(Data)
	
	if(is.null(Data$GeneID)&is.null(Data$ExonID)){
		DataTemp=data.frame(GeneID=GeneID,ExonID=ExonID)
		DataTemp=cbind(DataTemp,Data)
		Data=DataTemp
	}else if(is.null(GeneID)){
		GeneID=Data$GeneID
	}else if(is.null(ExonID)){
		ExonID=Data$ExonID
	}
	
	## From this point we assume that the first column of the data is a Gene ID and the second column is the Exon ID
	
	Transformation<-function(gID,Data){
		Subset=Data[which(Data$GeneID==gID),]
		
		generow=c()
		gID=as.character(gID)
		
		eID=unique(Subset$ExonID)
		lengthe=sapply(eID,function(x) return(length(which(Subset$ExonID==x))))
		
		eID=paste(as.character(eID),sep="",collapse=",")
		lengthe=paste(as.character(lengthe),sep="",collapse=",")
		
		#get all samples at once: apply on columns of the Subset data and discard the geneID and exonID columns
		samples=apply(Subset[,-c(1,2)],2,function(x) paste(as.character(x),sep="",collapse=","))
		allsamples=paste(samples,sep="",collapse=",")
		samplenames=paste(colnames(Data)[-c(1,2)],sep="",collapse=",")
		#put everything intro c()
		generow=c(gID,eID,lengthe,allsamples,samplenames)
		
	}
	#use rbindlist to get a full data file
	DataPivot=lapply(unique(GeneID),Transformation,Data)
	DataBind=t(rbindlist(list(DataPivot)))
	colnames(DataBind)=c("GeneID","ExonID","lengthexons","allsamples","samplenames")
	rownames(DataBind)=unique(GeneID)
	
	DataBind=as.data.frame(DataBind)
	DataBind$GeneID=as.character(DataBind$GeneID)
	DataBind$ExonID=as.character(DataBind$ExonID)
	DataBind$lengthexons=as.character(DataBind$lengthexons)
	DataBind$allsamples=as.character(DataBind$allsamples)
	DataBind$samplenames=as.character(DataBind$samplenames)
	
	if(is.null(location)){
			location=getwd()
	}
	location2=paste(location,"/",Name,".csv",sep="")	
	
	assign(Name,DataBind,envir=environment())
	utils::write.table(get(Name),file=location2,row.names=FALSE,col.names=TRUE,sep=",",quote=TRUE,qmethod="double")
}


#Pivot=PivotTransformData(Data=ColonCancer, GeneID=ColonCancer$GeneID,ExonID=ColonCancer$ExonID,savecsv=FALSE,Name="ColonCancer_Pivot",location="test")


##### REIDS FUNCTIONS ######


# I/NI Calls model

#' "iniREIDS"
#' 
#' The function is part of the larger REIDS function and performs the I/NI calls filtering model on the data for a single gene to select only the informative probesets. The method is performed if Informative=TRUE in the REIDS function before
#' before applying the REIDS model itself. The function is written to fit in the flow of the REIDS model
#' @export
#' @param SubgeneData A subset of the data. Particularly, a subset of the data corresponding to one gene.
#' @param nsim The number of iterations to perform.
#' @return A list with an item per gene. Per gene it is mentioned wich probesets are informative and which are not.
iniREIDS<- function(SubgeneData, nsim=1000) {
	
	# Preparation 
	nprobes <- nrow(SubgeneData)
	nsamples<- ncol(SubgeneData)
	nexon <- G <- length(unique(rownames(SubgeneData)))
	totalObs <- nprobes*nsamples
	geneExp<- as.vector(as.matrix(SubgeneData))
	probes <- as.factor(as.vector(matrix(c(1:nprobes),nrow=nprobes,ncol=nsamples)))
	grp <- as.vector(matrix(rownames(SubgeneData),nrow=nprobes,ncol=nsamples))
	id <-   as.factor(as.vector(t(matrix(c(1:nsamples),nrow=nsamples,ncol=nprobes))))
	
	# Fit of model
	if(length(levels(probes))==1){
		geneRes=rep(0,length(probes))
	}
	else{
		ft <- stats::lm(geneExp~probes-1,x=TRUE) # Why is the model fitted? To get resiudal for estimation of random intercepts?
		
		beta<- as.vector(summary(ft)$coefficient[,1]) #not used?
		fixedDesignMatrix <- as.matrix(as.data.frame(ft$x)) ## create designsamples matrix: not used?
		geneRes <- as.vector(summary(ft)$residuals)  #get residuals ==> for estimation of the ramdom intercepts??
	}
	geneResMat <- matrix(geneRes,nprobes,nsamples) #filled in by column: ok, one column per sample, calculation happens per sample
	
	zid <- matrix(c(1:nsamples),nrow=nsamples,ncol=G)
	uz <- unique(grp)
	zcls2<- sapply(uz,function(x)ifelse(grp==x,1,0))
	
	idv <- as.numeric(id)
	Z<-matrix(0,totalObs ,G*nsamples)
	for (ii in 1:nsamples) Z[idv==ii,(G*(ii-1)+1):(G*ii)]<-zcls2[idv==ii,]
	
	###########
	# Priors  #
	###########
	d0<-g0<-.001		# alpha and beta	
	tau <- 1 #sigma^-2 of measurement error
	
	#################
	# Store Results #
	#################
	nused <- floor(nsim/2)
	nburnin <-  nsim-nused
	outputbik <- matrix(0,nused,nsamples*G)
	outputerror <- NULL
	outputsigma <- matrix(0,nused,G)	# Fixed Effects
	
	
	###############################
	# Fixed Posterior Hyperparameter #
	#    for tau and taub		#
	###############################
	d<-d0+totalObs/2  #alpah+0.5*N
	
	###################
	# GIBBS SAMPLER   #
	###################
	
	nu0 <- G
	c0<-diag(G)			# Scale matrix for Wishart Prior onsamples Sigma.b
	nu<-nu0+nsamples   #gamma+n 
	Taub<-diag(G)			# Ransamplesdom Effects Prec Matrix
	z <- grp
	Taub <- diag(1,G,G) #not needed
	
	for (i in 1:nsim){
		
		set.seed(123*i+i)
		
		# Update b
		cZZ <- diag(colSums(Z))
		#vb <- diag(nsamples)%x%Taub+tau*cZZ #This is Var =D^-1+ sigma^(-2)*ni)
		vb <- solve(diag(nsamples)%x%Taub+tau*cZZ) # This is Var^(-1) =(D^-1+ sigma^(-2)*ni)^-1
		#vbInv=solve(vb)
		mb2<-(tau*crossprod(Z,geneRes ))  #Theta
		mb <- apply(vb,1,function(x) sum(x*mb2) ) #Var^-1*Theta
		vb1 <- diag(vb)  #only need variances here
		b<-sapply(c(1:(nsamples*G)),function(x) stats::rnorm(1,mean=mb[x],sd=sqrt(vb1[x])))  #aren't these the variances^-1 instead of the variances?
		bmat<-matrix(b,ncol=G,byrow=T)  		# Put insamples nsamples x G matrix form for updatinsamplesg taub
		
		
		# Update tau
		zb <- rowSums(sapply(c(1:G),function(x) t(matrix(bmat[,x],nsamples,nprobes))*matrix(zcls2[,x],nprobes,nsamples)))	
		geneReszb <-geneRes-zb	
		g<-g0+crossprod(geneReszb,geneReszb)/2
		tau<-stats::rgamma(1,d,g)
		
		
		# Update Taub 
		
		Sigma.b<-riwish(nu,c0+crossprod(bmat,bmat))
		Taub<-solve(Sigma.b)
		
		if(i > nburnin ){ 
			outputbik[(i-nburnin),] <- b
			outputerror[(i-nburnin)] <- 1/tau
			outputsigma[(i-nburnin),] <- as.vector(diag(Sigma.b))
			
		}
	}
	
	tp1 <- apply(outputsigma,2,function(x) mean(x))
	tp2 <- mean(outputerror)
	icc <- tp1/(tp1+tp2)
	
	Output <- icc
	return(Output)
}


# REIDS model

#' "REIDSmodel_intern"
#' 
#' The function is part of the larger REIDS function and performs the REIDS model on the data for a single gene.
#' @export
#' @param SubgeneData A subset of the data. Particularly, a subset of the data corresponding to one gene.
#' @param nsim The number of iterations to perform.
#' @return A list with 2 items per gene. The first item is the exon scores of the corresponding probesets and the second contains a data frame with the array scores of the exons across the samples.
#' If the iniREIDS model was performed. The items will be added to the previously made list.
REIDSmodel_intern<- function(SubgeneData, nsim=100) {
	
	nprobes <- nrow(SubgeneData)
	nsamples<- ncol(SubgeneData)
	nexon <- G <- length(unique(rownames(SubgeneData)))
	totalObservations <- nprobes*nsamples
	geneExp<- as.vector(as.matrix(SubgeneData))
	probes <- as.factor(as.vector(matrix(c(1:nprobes),nrow=nprobes,ncol=nsamples)))
	grp <- as.vector(matrix(rownames(SubgeneData),nrow=nprobes,ncol=nsamples))
	id <-   as.factor(as.vector(t(matrix(c(1:nsamples),nrow=nsamples,ncol=nprobes))))
	
	if(length(levels(probes))==1){
		geneRes=rep(0,length(probes))
	}
	else{
		ft <- stats::lm(geneExp~probes+id-1,x=TRUE)  #Difference from iniREIDS: id(samples) is involved in the calculation
		
		beta<- as.vector(summary(ft)$coefficient[,1])
		fixedDesignMatrix <- as.matrix(as.data.frame(ft$x)) ## create designsamples matrix
		geneRes <- as.vector(summary(ft)$residuals)
	}
	geneResMat <- matrix(geneRes,nprobes,nsamples)
	zid <- matrix(c(1:nsamples),nrow=nsamples,ncol=G)
	uz <- unique(grp)
	zcls2<- sapply(uz,function(x)ifelse(grp==x,1,0))
	idv <- as.numeric(id)
	Z<-matrix(0,totalObservations ,G*nsamples)
	for (ii in 1:nsamples) Z[idv==ii,(G*(ii-1)+1):(G*ii)]<-zcls2[idv==ii,]
	
	###########
	# Priors  #
	###########
	d0<-g0<-.001			
	tau <- 1
	
	#################
	# Store Results #
	#################
	nused <- floor(nsim/2)
	nburnin <-  nsim-nused
	outputbik <- matrix(0,nused,nsamples*G)
	outputerror <- NULL
	outputsigma <- matrix(0,nused,G)	# Fixed Effects
	###############################
	# Fixed Posterior Hyperparameter #
	#    for tau and taub		#
	###############################
	d<-d0+totalObservations/2
	
	###################
	# GIBBS SAMPLER   #
	###################
	
	nu0 <- G
	c0<-diag(G)			# Scale matrix for Wishart Prior onsamples Sigma.b
	nu<-nu0+nsamples
	Taub<-diag(G)			# Ransamplesdom Effects Prec Matrix
	z <- grp
	Taub <- diag(1,G,G)
	
	for (i in 1:nsim){
		
		set.seed(123*i+i)
		
		# Update b
		cZZ <- diag(colSums(Z))
		#vb <- diag(nsamples)%x%Taub+tau*cZZ #This is Var =D^-1+ sigma^(-2)*ni)
		vb <- solve(diag(nsamples)%x%Taub+tau*cZZ) # This is Var^(-1) =(D^-1+ sigma^(-2)*ni)^-1
		#vbInv=solve(vb)
		mb2<-(tau*crossprod(Z,geneRes))
		mb <- apply(vb,1,function(x) sum(x*mb2))
		vb1 <- diag(vb)
		b<-sapply(c(1:(nsamples*G)),function(x) stats::rnorm(1,mean=mb[x],sd=sqrt(vb1[x])))
		bmat<-matrix(b,ncol=G,byrow=T)  		# Put insamples nsamples x G matrix form for updatinsamplesg taub
		
		
		# Update tau
		zb <- rowSums(sapply(c(1:G),function(x) t(matrix(bmat[,x],nsamples,nprobes))*matrix(zcls2[,x],nprobes,nsamples)))	
		geneReszb <-geneRes-zb	
		g<-g0+crossprod(geneReszb,geneReszb)/2
		tau<-stats::rgamma(1,d,g)
		
		
		# Update Taub 
		
		Sigma.b<-riwish(nu,c0+crossprod(bmat,bmat))
		Taub<-solve(Sigma.b)
		
		if(i > nburnin ){ 
			outputbik[(i-nburnin),] <- b
			outputerror[(i-nburnin)] <- 1/tau
			outputsigma[(i-nburnin),] <- as.vector(diag(Sigma.b))
			
		}
	}
	
	tp1 <- apply(outputsigma,2,function(x) stats::quantile(x, probs = c(0.025, 0.5, 0.95)))
	tp2 <- stats::quantile(outputerror, probs = c(0.025, 0.5, 0.95))
	icc <- t(apply(tp1,2,function(x) x/(x+tp2)))
	icc <- data.frame(exon=uz,type="icc",icc )
	nik <- as.vector(t(sapply(uz,function(x) rep(x,nsamples))))
	ikEffect2  <- t(apply(outputbik,2,function(x) stats::quantile(x, probs = c(0.025, 0.5, 0.95))))
	ikEffect<- data.frame(exon=nik,type="bik",ikEffect2)
	exonscore <- icc[,c(1,4)]
	arrayscore <- ikEffect[,4] 
	
	dim(arrayscore)<- c(G,nsamples)
	rownames(arrayscore)<-ikEffect$exon[1:G] 
	colnames(arrayscore) <- colnames(SubgeneData)
	Output <- list(exonScores=exonscore,arrayScores=arrayscore)
	return(Output)
}


# REIDS Function - Cluster version
# This function is accompagnied with a .pbs file in the documentation folder. It is advised to run this model an a HPC cluster and not on a regular laptop as it will
# consume time and memory

#' "REIDS_ClusterVersion"
#' 
#' The REIDS_ClusterVersion performs the REIDS model and was adapted for use on a HPC cluster. This function should be used with the REIDS_ClusterVersion.R file and REIDS_ClusterVersion.pbs script in the documentation folder of the package.
#' After running this function on the cluster, the output files should be binded together with the CreateOutput function.
#' @export
#' @param geneData The data with as rows the probesets and as columns the samples. Note that the first column should contain the gene IDs and the second column the exon IDs
#' @param nsim The number of iterations to perform.
#' @param geneID A vector of the gene IDs.
#' @param informativeCalls Logical. Should the I/NI calls method be perform before applying the REIDS model?
#' @param alpha The threshold for filtering in the I/NI calls method. Probesets with scores higher than alpha are kept.
#' @return A .RData file will be saved for each gene with the elements returned by the iniREIDS and REIDS functions.
REIDSFunction_ClusterVersion<- function(geneData,nsim=1000,geneID,informativeCalls=TRUE,alpha=0.5){
	
	exonScore <- arrayScore <- informativeData<- NULL
	
	output=list()
	output[[1]]=list()
	names(output)[1]=geneID
	lcmmData <- geneData[,-c(1,2)]
	lcmmData=as.matrix(lcmmData)
	enames <- geneData$exonID
	names(enames)=NULL
	
	rownames(lcmmData)<- enames
	
	##informative calls
	i=1
	if(informativeCalls){			
		fit <- iniREIDS(SubgeneData=lcmmData, nsim) 
		fit2 <- data.frame(exonNames=unique(enames),Score=fit,informative=fit>alpha) # iniREIDS returns one value per exon: filtering on exon level, no replicates
		output[[1]][[i]]=fit2
		names(output[[1]])[i]="Informative"
		iniData <- lcmmData[which(rownames(lcmmData)%in%fit2$exonNames[fit2$informative]),] # Of those that pass filtering step, retrieve the replicates and samples
		lcmmData <-  iniData   
		i=i+1
	}
	
	
	if(!is.null(lcmmData)&length(unique(rownames(lcmmData)))>1){  
		fit <- REIDSmodel_intern(SubgeneData=lcmmData, nsim) 
		exonScore <-fit$exonScores
		arrayScore <- fit$arrayScores
		output[[1]][[i]]=exonScore
		names(output[[1]])[i]="exonScore"
		output[[1]][[i+1]]=arrayScore
		names(output[[1]])[i+1]="arrayScore"
	}
	
	return(output)
}


# CreateOutput

#' "CreateOutput"
#' 
#' The CreateOutput functions bind the .RData files returned by the REIDS_ClusterVersion together into one list with an element per gene. The function is advised to be used with the CreateOutput.R and CreateOutput.pbs file in the documentation folder.
#' @export
#' @param ID A list of gene IDs. These are read in via a csv file in the CreateOutput.R file in the documentation folder.
#' @param Name A name for the returned list.
#' @return A list with the output of the REIDS_ClusterVersion binded together for all genes.
CreateOutput<-function(ID,Name){
	Output=list()
	for(i in ID$geneID){
		Data=get(load(paste("REIDS_Gene_",as.character(i),".RData",sep="")))
		
		Output[length(Output)+1]=Data
		names(Output)[length(Output)]=i
	}
	
	assign(paste(Name,"REIDS_Output",sep="_"),Output,envir=environment())
	eval(parse(text=paste("save(",Name, "_REIDS_Output, file=\"/folder/",Name, "_REIDS_Output.RData\")", sep=""))) 
}


# REIDS Function - Regular version
# It is advised to run this model on a HPC cluster and not on a regular laptop as it will consume time and memory

#' "REIDSFunction"
#' 
#' The REIDSFunction performs the REIDS model on the output of the DataProcessing function and is a regular R function. It is advised to use this function in its cluster version and not on a regular laptop as it will consume time and memory.
#' @export
#' @param geneData The data with as rows the probesets and as columns the samples. Note that the first column should contain the gene IDs and the second column the exon IDs
#' @param nsim The number of iterations to perform.
#' @param geneID A vector of the gene IDs.
#' @param exonID A vector of the exon IDs
#' @param informativeCalls Logical. Should the I/NI calls method be perform before applying the REIDS model?
#' @param alpha The threshold for filtering in the I/NI calls method. Probesets with scores higher than alpha are kept.
#' @return A list with an element for each gene with per gene the values returned by the iniREIDS and REIDS functions.
REIDSFunction <- function(geneData,nsim=1000,geneID,exonID,informativeCalls=TRUE,alpha=0.5){
	
	exonScore <- arrayScore <- informativeData<- NULL
	
	if(is.null(rownames(geneData))|is.null(names(exonID))){
		message("No rownames were specified for geneData or exonID. The geneID's in the provided order will be taken as the names.")
		rownames(geneData)=geneID
		names(exonID)=geneID
	}
	else if(all(!(unique(rownames(geneData))%in%unique(geneID)))|all(!(unique(names(exonID))%in%unique(geneID)))){ 
		stop("The rownames of geneData do not match the geneID's. Please put the correctly oredered geneID's as rownames for geneData.")
	}
	
	output=list()
	for(e in 1:length(unique(geneID))){ 
		output[[e]]=list()
		names(output)[e]=unique(geneID)[e]
		lcmmData <- geneData[which(rownames(geneData)==unique(geneID)[e]),]  
		enames <- exonID[which(names(exonID)==unique(geneID)[e])]
		names(enames)=NULL
		rownames(lcmmData)<- enames
		
		##informative calls
		i=1
		if(informativeCalls){			
			fit <- iniREIDS(SubgeneData=lcmmData, nsim) 
			fit2 <- data.frame(exonNames=unique(enames),Score=fit,informative=fit>alpha) # iniREIDS returns one value per exon: filtering on exon level,no replicates
			output[[e]][[i]]=fit2
			names(output[[e]])[i]="Informative"
			iniData <- lcmmData[which(rownames(lcmmData)%in%fit2$exonNames[fit2$informative]),] # Of those that pass filtering step, retrieve the replicates and samples
			lcmmData <-  iniData   
			i=i+1
		}
		
		
		if(!is.null(lcmmData)&length(unique(rownames(lcmmData)))>1){  
			fit <- REIDSmodel_intern(SubgeneData=lcmmData, nsim) 
			exonScore <-fit$exonScores
			arrayScore <- fit$arrayScores
			output[[e]][[i]]=exonScore
			names(output[[e]])[i]="exonScore"
			output[[e]][[i+1]]=arrayScore
			names(output[[e]])[i+1]="arrayScore"
		}
		
	}
	return(output)
}



##### OUTPUT ANALYSIS FUNCTIONS #####

#FILTERING FUNCTION

#' "FilterInformativeGenes"
#'
#' The FilterInformativeGenes function works on the list returned by the REIDS functions. It returns a list that only retains the genes that have informative probesets.
#' @export 
#' @param Data The output of the REIDS function.
#' @return A subset of the list provided in Data.
FilterInformativeGenes <- function(Data){
	Out1=list()
	filter <- function(Subset){
		if("exonScore"%in% names(Subset)){
			return(Subset)
		}
	}
	Out1=lapply(Data,filter)
	
	Out2=Out1[vapply(Out1, Negate(is.null), NA)]
	
	return(Out2)
}

#Testing Search on 1) list of the REIDS Output with and without filtering
# on the REIDS Output
#load("TestData/ColonCancer_OutputREIDSModel.RData")
#test1=FilterInformativeGenes(ColonCancer_OutputREIDSModel)
#
#load("TestData/ColonCancer_OutputREIDSModel_NoFilter.RData")
#test2=FilterInformativeGenes(ColonCancer_OutputREIDSModel_NoFilter)


#SEARCH FUNCTION

#' "Search"
#' 
#' The Search function investigates whether specific genes and/or exons are present in the data. It returns the elements of the found elements and a list of those that were not found.
#' @export 
#' @param WhatToLookFor A data frame with a GeneID and ExonID column. Only the ExonID column is necessary, the corresponding gene ID will be sought in the data
#' @param Data The data in which the given gene and exon ID should be sought.
#' @param AggregateResults Logical. Should the results be aggregated on gene level? This results in a list with one item per gene
#' @param NotFound Not be specified by the user.
#' @return The returned value is a list with two elements. The first element is SearchResults which contains a list of the found intstances in the data. The list contains an element per found instance. If AggregateResults is TRUE, the list is reduced to an element per gene. The secod is a data frame calles NotFound which contain the gene and exon ID which were not found in the data. If only exon ID were specified, the gene ID are NA in this data frame.
Search <- function(WhatToLookFor=data.frame(GeneID=NULL,ExonID=NULL), Data, AggregateResults=FALSE,NotFound=NULL){
	if(!(is.null(WhatToLookFor$GeneID))){
		
		if(is.null(NotFound)){
			NotFound<-data.frame(GeneID=integer(),ExonID=integer())
		}
		
		if(class(Data[[1]])=="list"){
			retrieve<-function(GID,EID,Data){
				
				listnumber=which(names(Data)==as.character(GID))
				
				if(length(listnumber)==0){
					NotFound=rbind(NotFound,data.frame(geneID=GID,exonID=EID))
					#assign("NotFound", notfound, envir = .GlobalEnv)
					
					ExonInfo=NULL
					#return(ExonInfo)
				}
				else{
					SelectGeneData=Data[[listnumber]]
					

					if(class(SelectGeneData)=="list"){ #here we assume that the data is the direct output of the REIDSfunction
						
						if(names(SelectGeneData)[1]=="exonScore"){
							ExonScore=SelectGeneData$exonScore[which(SelectGeneData$exonScore$exon==EID),,drop=FALSE]
							ArrayScore=SelectGeneData$arrayScore[which(rownames(SelectGeneData$arrayScore)==EID),,drop=FALSE]	
							
							ExonInfo=list(exonScore=ExonScore,arrayScore=ArrayScore)
						}
						else{
							InformativeExon=SelectGeneData$Informative[which(SelectGeneData$Informative$exonNames==EID),,drop=FALSE]
						
							if(InformativeExon$informative==TRUE & length(which(SelectGeneData$Informative$informative==TRUE))>1){
								ExonScore=SelectGeneData$exonScore[which(SelectGeneData$exonScore$exon==EID),,drop=FALSE]
								ArrayScore=SelectGeneData$arrayScore[which(rownames(SelectGeneData$arrayScore)==EID),,drop=FALSE]
						
								ExonInfo=list(Informative=InformativeExon,exonScore=ExonScore,arrayScore=ArrayScore)
							}
							else{
								ExonInfo=list(Informative=InformativeExon)
							}
						}
					
					}
					else if(class(SelectGeneData)=="data.frame"){
						ExonInfo=SelectGeneData[which(SelectGeneData$ExonID==EID),,drop=FALSE]
					
						if(nrow(ExonInfo)==0){
							ExonInfo=NULL
						
							NotFound=rbind(NotFound,data.frame(geneID=GID,exonID=EID))
							#assign("NotFound", notfound, envir = .GlobalEnv)
						}
					
					}
				}	
				
				return(ExonInfo)
				
			}
			
		}
		
		else if(class(Data)=="data.frame"){
			retrieve<-function(GID,EID,Data){
				
				Rownumbers=which(Data$GeneID==GID)
				if(length(Rownumbers)==0){
					NotFound=rbind(NotFound,data.frame(geneID=GID,exonID=EID))
					
					ExonInfo=NULL

				}
				else{
					
					SelectGeneData=Data[Rownumbers,]
					
					ExonInfo=SelectGeneData[which(as.character(SelectGeneData$ExonID)==as.character(EID)),,drop=FALSE]
					
					if(nrow(ExonInfo)==0){
						ExonInfo=NULL
						
						NotFound=rbind(NotFound,data.frame(geneID=GID,exonID=EID))
						#assign("NotFound", notfound, envir = .GlobalEnv)
					}			
				}
			
				return(ExonInfo)		
			}
		}
		
		Out1=lapply(1:nrow(WhatToLookFor),function(x) retrieve(GID=WhatToLookFor[x,1],EID=WhatToLookFor[x,2],Data = Data))
		names(Out1)=WhatToLookFor$GeneID
		Out2=Out1[vapply(Out1, Negate(is.null), NA)]
		
		SearchResults=Out2
		
		if(AggregateResults==TRUE){
			message("Results are aggregated on the gene level...")
			
			genes=unique(names(Out2))
			
			aggregate<-function(gene,Data){
				if(class(Data[[1]])=="data.frame"){  #assume a data.frame per gene
					Out3=rbindlist(Data)
					Out3=as.data.frame(Out3)
				}
				else if(class(Data[[1]])=="list"){ #assume a list of data frames per gene
					
					ntables=max(sapply(Data, length))
					
					aggregatetables<-function(data){
						data=lapply(data,as.data.frame)
						Out4=rbindlist(data)
						Out4=as.data.frame(Out4)
						return(Out4)					 
					}
					
					Out3=lapply(1:ntables,function(i) aggregatetables(data=lapply(Data,function(x) ifelse(length(x)>=i,x[i],x[NA]))))	
					names(Out3)=names(Data[[ which(sapply(Data,length)==max(sapply(Data, length)))[1]]])
					
				}
				
				return(Out3)
				
			}
			
			Out5=lapply(1:length(genes),function(x) aggregate(gene = genes[x], Data=Out2[which(names(Out2)==genes[x])]))
			names(Out5)=genes
			
			SearchResults=Out5
			
		}
		
		
	}
	
	else if(is.null(WhatToLookFor$GeneID) & !(is.null(WhatToLookFor$ExonID))){
		
		if(class(Data[[1]])=="list"){
			
			if(names(Data[[1]])[1]=="exonScore"){
				ExonNames=lapply(1:length(Data),function(x) as.character(Data[[x]]$exonScore$exon))
			}
			else{
				ExonNames=lapply(1:length(Data),function(x) as.character(lapply(Data[x],"[[",1)[[1]]$exonNames))
			}
			
			names(ExonNames)=names(Data)		
			GeneNames=sapply(1:length(ExonNames),function(i) rep(names(ExonNames)[i],length(ExonNames[[i]])))			
			GenesAndExons=data.frame(GeneID=as.integer(unlist(GeneNames)),ExonID=as.integer(unlist(ExonNames)))
			
			exons=WhatToLookFor$ExonID
			
			LookFor=lapply(exons, function(x) GenesAndExons[which(GenesAndExons$ExonID==x),])			
			LookFor=as.data.frame(rbindlist(LookFor))
			
			exons_notfound=exons[which(!exons%in%LookFor$ExonID)]
			
			if(length(exons_notfound)==0){
				NF=NULL
			}
			else{
				NF=data.frame(GeneID=rep(NA,length(exons_notfound)),ExonID=exons_notfound)
			}
			
		}
		else if(class(Data)=="data.frame"){
			
			#ExonNames=lapply(1:length(Data),function(x) Data[x][[1]]$ExonID)
			ExonNames=Data$ExonID
			GeneNames=Data$GeneID
		
			
			#GeneNames=sapply(1:length(ExonNames),function(i) rep(names(ExonNames)[i],length(ExonNames[[i]])))			
			#GenesAndExons=data.frame(geneID=as.integer(unlist(GeneNames)),exonID=as.integer(unlist(ExonNames)))
			GenesAndExons=data.frame(GeneID=GeneNames,ExonID=ExonNames)
			exons=WhatToLookFor$ExonID
			
			LookFor=lapply(exons, function(x) ifelse(nrow(GenesAndExons[which(GenesAndExons$ExonID==x),])!=0,return(GenesAndExons[which(GenesAndExons$ExonID==x),]),return(data.frame(GeneID="NA",ExonID=x))))
			LookFor=as.data.frame(rbindlist(LookFor))
			LookFor$GeneID=suppressWarnings(as.integer(as.character(LookFor$GeneID)))
			exons_notfound=LookFor[which(is.na(LookFor$GeneID)),]
			LookFor=LookFor[-which(is.na(LookFor$GeneID)),]
			
			if(nrow(exons_notfound)==0){
				NF=NULL
			}
			else{
				NF=exons_notfound
			}
			
		}
		
		return(Search(WhatToLookFor=LookFor, Data=Data, AggregateResults=AggregateResults,NotFound=NF))
	}
	
	NMatched=length(Out2)
	message(paste(NMatched," instances found in the Data..."))
	return(list("SearchResults"=SearchResults,"NotFound"=NotFound))
}

##Testing Search on 1) list of the REIDS Output and 2) a data frame of processed REIDS Output
#exonID <- c(3252129,3597384,3333718,3735208,2598321,3338589,2605391,2605390,
#		2605386,3025632,2375766,3569827,3569830,2334499,3972987,2516011,
#		2989068,3422189)
#
## on the REIDS Output
#load("TestData/ColonCancer_OutputREIDSModel.RData")
#test1=Search(WhatToLookFor=data.frame(ExonID=exonID), Data=ColonCancer_OutputREIDSModel, AggregateResults=FALSE,NotFound=NULL)
#
#load("TestData/ColonCancer_OutputREIDSModel_NoFilter.RData")
#test2=Search(WhatToLookFor=data.frame(ExonID=exonID), Data=ColonCancer_OutputREIDSModel_NoFilter, AggregateResults=FALSE,NotFound=NULL)
#
#
## on a data.frame
#load("TestData/MeanPairedDiff_REIDSOutput.RData")
#colnames(MeanPairedDiff_REIDSOutput)[1]="GeneID"
#test3=Search(WhatToLookFor=data.frame(ExonID=exonID), Data=MeanPairedDiff_REIDSOutput, AggregateResults=FALSE,NotFound=NULL)


#TESTING FUNCTION
# input : output of the REIDS model ( a list per gene)
# procedure : filter on ICC, testing on arrayscores, adjusted p-values, filter on adjuster p-values
# output : a data frame

#' "ExonTesting"
#' 
#' The ExonTesting function performs a t-test on the array score of predefined groups. If specified, probesets are filtered out on exon scores and test significance.
#' @export
#' @param Data The Data on which testing of the array scores should be conducted. This is preferably output of the REIDS function
#' @param Exonthreshold The exon score threshold to be maintained. If not NULL, probesets with an exon score lower than this value are not considered further and the p-values will be adjusted for multiplicity after testing. If NULL, all probesets are considered and a multiplicity correction is not performed.
#' @param groups A list with two elements speficifing the columns of the data of group 1 in group1 and those of group 2 in group2.
#' @param paired Logical. Are the groups paired? If TRUE the mean paired differences are calculated and tested whether these are significantly different from zero or not.
#' @param significancelevel The significance level to be maintained on the p-values. The filtering on the significance is conducted only if an Exonthreshold is specified and the p-value are adjusted for multiplicity.
#' @return A data frame with one line per exon. The columns contain the gene ID, the exon ID, the test statistic, a p-value and an adjusted p-value. If the groups are paired also the mean paired difference is given. The p-values are adjusted for multiplicity and filtered on significance if significancelevel is not NULL.
ExonTesting <- function(Data, Exonthreshold=NULL,groups=list(group1=NULL,group2=NULL),paired=FALSE,significancelevel=NULL){
	if(is.null(Exonthreshold)){
		Exonthreshold=0
	}
	
	if(is.null(groups$group1)){
		message("A test on the array scores will NOT be performed")
	}
	
	message("Data is filtered to only contain genes with informative exons")
	Data_Filtered1=FilterInformativeGenes(Data)  #only those genes with informative exons remain
	
	filterASExons<-function(i,geneID,Subset,Exon,groupings,pairing){
		print(i)
		
		## Step 1 : filtering on the Exon value. If 0, no filtering occurs
		ExonsExon = Subset$exonScore
		ArrayScores = Subset$arrayScore
		
		SelectExonsExon = ExonsExon[which(ExonsExon$X50>Exon),]
		ExonsPassedExon = SelectExonsExon$exon
		
		SelectRowsArray = which(rownames(ArrayScores)%in%ExonsPassedExon)
		
		## Step 2 : Testing of the Array Scores -- paired or not paired
		if(pairing==FALSE){  # Test between two groups of Array Scores
			
			if(!(is.null(groupings$group1)) & length(SelectRowsArray)!=0){
				
				ArrayScore_group1=ArrayScores[SelectRowsArray,groupings$group1,drop=FALSE]
				ArrayScore_group2=ArrayScores[SelectRowsArray,groupings$group2,drop=FALSE]
				
				ttest<-function(i,g1,g2,pairs){
					out1=stats::t.test(x=g1,y=g2)
					out2=cbind(out1$statistic,out1$p.value)
					
					return(out2)
				}
				
				ArrayScoreTest = t(sapply(c(1:length(ExonsPassedExon)),function(i) ttest(g1=ArrayScore_group1[i,],g2=ArrayScore_group2[i,], pairs = pairing) ))
				ArrayScoreTest=as.data.frame(ArrayScoreTest)
				colnames(ArrayScoreTest)=c("t.statistic","p.value")
				rownames(ArrayScoreTest)=rownames(ArrayScore_group1)
				
				
			}
			else{
				ArrayScoreTest = NULL
			}
		}

		
		if(pairing==TRUE){ # Test the mean paired difference against zero
			
			
			if(!(is.null(groupings$group1)) & length(SelectRowsArray)!=0){
				
				ArrayScore_group1=ArrayScores[SelectRowsArray,groupings$group1,drop=FALSE]
				ArrayScore_group2=ArrayScores[SelectRowsArray,groupings$group2,drop=FALSE]
				
				mean_paired_diff<-function(g1,g2){
					Paired_Diff=g1-g2
					Mean_Diff=mean(Paired_Diff)
					
					out1=stats::t.test(x=Paired_Diff)
					out2=cbind(out1$statistic,out1$p.value)
					out3=cbind(Mean_Diff,out2)
					
					return(out3)
				}
				
				ArrayScoreTest = t(sapply(c(1:length(ExonsPassedExon)),function(i) mean_paired_diff(g1=ArrayScore_group1[i,],g2=ArrayScore_group2[i,]) ))
				ArrayScoreTest=data.frame("Mean_Diff"=ArrayScoreTest[,1],"t-statistic"=ArrayScoreTest[,2],p.value=ArrayScoreTest[,3])
				rownames(ArrayScoreTest)=rownames(ArrayScore_group1)
			
			
			}
			else{
				ArrayScoreTest = NULL
			}
		}	
		
		if(!is.null(ArrayScoreTest)){
			Output=cbind("geneID"=rep(geneID,length(ExonsPassedExon)), SelectExonsExon,ArrayScoreTest)			
		}
		else{
			Output=cbind("geneID"=rep(geneID,length(ExonsPassedExon)),SelectExonsExon)
		}
		
		colnames(Output)[2]="ExonID"
		
		if(nrow(Output)==0){
			Output=NULL
		}
		
		return(Output)
		
	}
	
	Out1=lapply(1:length(Data_Filtered1), function(i) filterASExons(i,geneID = names(Data_Filtered1[i]), Subset = Data_Filtered1[[i]], Exon = Exonthreshold , groupings=groups, pairing = paired))
	names(Out1)=names(Data_Filtered1)
	
	Out2=Out1[vapply(Out1, Negate(is.null), NA)]
	
	Out3<-do.call(rbind.data.frame, Out2)
	
	
	## Step 3 : if Exon threshold specified: we adjust for multiplicity ;  filter on significance level 
	if(nrow(Out3)!=0){
		message("The p-value are adjusted")
		if(!is.null(Out3$p.value)){
			Out3$adj.p.value=stats::p.adjust(Out3$p.value,"fdr")
		}	
		if(!is.null(significancelevel)){
			Out3=Out3[which(Out3$adj.p.value<=significancelevel),]
		}
		rownames(Out3)=seq(1:nrow(Out3))
	}
	return(Out3)
	
}

## Identification of the AS exons

#' "ASExons"
#' 
#' The ASExons functions can be performed either on the output of the REIDS model or of the ExonTesting model and identifies the alternatively spliced exons. It filters probesets on their exon scores, adjusts p-values for multiplicity and only keeps the significant probesets.
#' @export
#' @param Data The Data on which testing of the array scores should be conducted. This can be either the output of the REIDS model or the ExonTesting function.
#' @param Exonthreshold The exon score threshold to be maintained. If not NULL, probesets with an exon score lower than this value are not considered further and the p-values will be adjusted for multiplicity after testing. If NULL, all probesets are considered and a multiplicity correction is not performed.
#' @param groups A list with two elements speficifing the columns of the data of group 1 in group1 and those of group 2 in group2.
#' @param paired Logical. Are the groups paired? If TRUE the mean paired differences are calculated and tested whether these are significantly different from zero or not.
#' @param significancelevel The significance level to be maintained on the p-values. The filtering on the significance is conducted only if an Exonthreshold is specified and the p-value are adjusted for multiplicity.
#' @return A data frame with one line per exon. The columns contain the gene ID, the exon ID, the test statistic, a p-value and an adjusted p-value. If the groups are paired also the mean paired difference is given. Only the probesets with high enough exon scores and a significant test are kept in the data frame.
ASExons<-function(Data,Exonthreshold=0.5,groups=list(group1=NULL,group2=NULL),paired=FALSE,significancelevel=0.05){
	
	if(class(Data)=="list"){
		message("The data is assumed to be output of the REIDS model. Filtering of the probesets and testing of the array scores will be performed")
		message("The used threshold for the exon scores is 0.5")
		message("The used significance level for the p-values is 0.05")
		
		TestedData=ExonTesting(Data=Data,Exonthreshold=Exonthreshold,groups=groups,paired=paired,significancelevel=significancelevel)
		
		if(nrow(TestedData)!=0){
#			message(paste("Keep probesets with exon score greater than",Exonthreshold,sep=" "))
#			Data_Filt1=TestedData[which(TestedData$X50.>Exonthreshold),]
#		
#			message("Adjusting p-values for multiplicity")	
#			Data_Filt1$adj.p.value=p.adjust(Data_Filt1$p.value,"fdr")
#		
#			message(paste("Keep probesets with a p-value lower than",significancelevel,sep=" "))
#			Data_Sign=Data_Filt1[which(Data_Filt1$adj.p.value<0.05),]
		
		
			message("Ordering data in from high tolow exon scores")
			Data_Sign_Ordered=TestedData[order(-TestedData$X50.),]
			rownames(Data_Sign_Ordered)=c(1:nrow(Data_Sign_Ordered))
		}
		else{
			Data_Sign_Ordered=NULL
		}
		
	}
	else if(class(Data)=="data.frame"){
		
		message("In using this function please make sure that the data has not been filtered yet and still has unadjusted p-values.")
		
		message(paste("Keep probesets with exon score greater than",Exonthreshold,sep=" "))
		Data_Filt1=Data[which(Data$X50.>Exonthreshold),]
		
		if(nrow(Data_Filt1)!=0){
			message("Adjusting p-values for multiplicity")	
			Data_Filt1$adj.p.value=stats::p.adjust(Data_Filt1$p.value,"fdr")
		
			message(paste("Keep probesets with a p-value lower than",significancelevel,sep=" "))
			Data_Sign=Data_Filt1[which(Data_Filt1$adj.p.value<significancelevel),]

			if(nrow(Data_Sign)!=0){
				message("Ordering data in from high tolow exon scores")
				Data_Sign_Ordered=Data_Sign[order(-Data_Sign$X50.),]
				rownames(Data_Sign_Ordered)=c(1:nrow(Data_Sign_Ordered))
			}
			else{
				Data_Sign_Ordered=NULL
			}

		}
		else{
			Data_Sign_Ordered=NULL
		}
	}
	return(Data_Sign_Ordered)
	
}

#testing function ExonTesting and ASExons
#load("TestData/ColonCancer_OutputREIDSModel.RData")
#testdata=ColonCancer_OutputREIDSModel[1:500]
#groupT=c(1,14,16,18,20,3,5,7,9,11)
#groupN=c(12,15,17,19,2,4,6,8,10,13)
#groups=list(group1=groupN,group2=groupT)
#
#
#test1=ExonTesting(Data=testdata,Exonthreshold=NULL,groups=groups,paired=TRUE,significancelevel=NULL)
#
#test2=ASExons(Data=testdata,Exonthreshold=NULL,groups=groups,paired=TRUE,significancelevel=NULL)
#
#test3=ASExons(Data=test2,Exonthreshold=0.5,groups=groups,paired=TRUE,significancelevel=0.05)

## FIRMA Model Analysis

#' "FIRMAScores"
#' 
#' The FIRMAScores function performs an analyais on the FIRMA scores of the sanples. If the groups are not paired, the Firma all sample score will be the minimum value of group 1. A test statistic is perforrmed on the grouping of interest to see if there is a significant difference between them. If the data is paired, the mean paired differences are obtained and tested versus zero.
#' @export
#' @param Data The output of the FIRMA model
#' @param InformativeExons A character vector of exon IDs. As for the REIDS model probesets are filtered out by I/NI calls model and later on exon score, the remaining exons can be specified here. Only these shall be considered in the FIRMA analysis to make the results between REIDS and FIRMA more comparable
#' @param groups A list with two elements speficifing the columns of the data of group 1 in group1 and those of group 2 in group2. Here is the possibility to specifify multiple subgroups in group 1 for the calculation of the all FIRMA sample score.
#' @param paired Logical. Are the groups paired? If TRUE the mean paired differences are calculated and tested whether these are significantly different from zero or not.
#' @param significancelevel If specified, filtering is conducted on the p-values.
#' @return A data frame with one line per exon. The columns contain the gene ID, the exon ID, the test statistic, a p-value and an adjusted p-value. If the groups are paired also the mean paired difference is given. Only the probesets with high enough exon scores and a significant test are kept in the data frame.
#' @details The input to this function is the output of the FIRMA model. On the FIRMA scores, a t-test is performed. This is either between groups if the data is not paired and on the mean paired differences if the data is paired. If no pairing is present, an all sample FIRMA score is computed as the maximum of the minimum values of the 
#' of the subgroups in group 1. The returned p-value is adjusted for multiplicity. If a vector is given in InformativeExons, only these exon IDS are considered in the calculations.
FIRMAScores <- function(Data,InformativeExons=NULL,groups=list(group1=list(group1a=NULL,group1b=NULL),group2=NULL),paired=FALSE,significancelevel=NULL){
#	if(is.null(groups$group1)){
#		message("A test on the array scores will NOT be performed")
#	}
	
	if(!is.null(InformativeExons)){
		Data=Data[which(Data$ExonID%in%as.character(InformativeExons)),]
	}
	
	data=Data[,-c(1,2)]
	Data=Data[-which(is.na(data)),]
	data=data[-which(is.na(data)),]
	
	if(paired==FALSE){
		
		# 1) Compute the FIRMA scores as Sj=max(min(Fij_group1a),min(Fij_group1b))
		subgroupsgroup1=groups$group1
		minofeachgroup=sapply(subgroupsgroup1,function(x) apply(data[,x],1,min))
	
		AllSample_FIRMAScore=apply(minofeachgroup,1,max)	
	
		Out1=data.frame('geneID'=Data[,1],'exonID'=Data[,2],'AllSample_FIRMAScore'=AllSample_FIRMAScore)
	
	
		# 2) Perform a t-test on the FIRMA scores
	
		data_group1=data[,unlist(groups$group1)]
		data_group2=data[,groups$group2]
	
		filterASExons1<-function(i,group1,group2){
			out1=stats::t.test(x=group1,y=group2)
			out2=cbind(out1$statistic,out1$p.value)
			return(out2)
		
		}
	
		Out_ttest=t(sapply(1:nrow(data),function(i) filterASExons1(i,group1=data_group1[i,],group2=data_group2[i,])))
	
		Out2=cbind(Out1,Out_ttest)	
		colnames(Out2)=c("GeneID","ExonID","AllSample_FIRMAScore","t.statistic","p.value")
		Out2$adj.p.value=stats::p.adjust(Out2$p.value,"fdr")
		Out2$GeneID=as.character(Out2$GeneID)
		Out2$ExonID=as.character(Out2$ExonID)
	}
	
	else if(paired==TRUE){
		
		filterASExons2<-function(Subset){
			out1=stats::t.test(x=Subset)
			out2=cbind(out1$statistic,out1$p.value)
			return(out2)
			
		}
		data_group1=data[,groups$group1]
		data_group2=data[,groups$group2]
		
		Paired_Diff=as.matrix(data_group1-data_group2)
		
		Mean_Diff=apply(Paired_Diff,1,mean)
		
		Out=data.frame('geneID'=Data[,1],'exonID'=Data[,2],'Mean_Diff'=Mean_Diff)
		Paired_Diff=Paired_Diff[-which(is.na(Out$Mean_Diff)),]
		Out=Out[-which(is.na(Out$Mean_Diff)),]
		
		Out_ttest=t(sapply(1:nrow(Paired_Diff),function(i) filterASExons2(Subset=Paired_Diff[i,])))
		
		Out2=cbind(Out,Out_ttest)
		
		colnames(Out2)=c("GeneID","ExonID","Mean_Diff","t.statistic","p.value")
		Out2$adj.p.value=stats::p.adjust(Out2$p.value,"fdr")
		Out2$GeneID=as.character(Out$GeneID)
		Out2$ExonID=as.character(Out$ExonID)

	}
	
	
	if(!is.null(significancelevel)){
		Out2$adj.p.value=stats::p.adjust(Out2$p.value,"fdr")
		Out2=Out2[which(Out2$adj.p.value<0.05),]
	}
	
	return(Out2)
	
}


#PLOTTING FUNCTION

#' "PlotFunction"
#' 
#' The PlotFunction produces three plots concerning a specific exon and its corresponding gene.
#' @export
#' @param GeneID The gene ID of the gene of interest.
#' @param ExonID The exon ID of the exon of interest.
#' @param Data The processed data as returned by DataProcessing. This is were the observed probe intensities will be retrieved.
#' @param REIDS_Output The output of the REIDS model. This is were the array scores will be retrieved.
#' @param GeneLevelData The gene level summarized data to retrieve the gene level values.
#' @param ExonLevelData The exon level summarized data to retrieve the exon level values.
#' @param FIRMA_Data The output of the FIRMA model to retrieve the FIRMA scores of the samples.
#' @param groups The groups of interest in the data.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a location should be provided in "location" and the figure is saved there. If "new" a new graphic device is opened and if "sweave", the figure is made compatible to appear in a sweave or knitr document, i.e. no new device is opened and the plot appears in the current device or document.
#' @param location If plottype is "pdf", a location should be provided in "location" and the figure is saved here. the location should include the name of the plot. 

PlotFunction<-function(GeneID=NULL,ExonID=NULL,Data,REIDS_Output,GeneLevelData=NULL,ExonLevelData=NULL,FIRMA_Data=NULL,groups,plottype="new",location){
	message("The gene and exon level data will be log2 transformed")
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	
	if(is.null(GeneID) | is.null(ExonID)){
		stop("no GeneID and/or ExonID specified")
	}
	
	# Array Scores
	ArrayScores=REIDS_Output[[which(names(REIDS_Output)==GeneID)]]$arrayScore # the arrayscores for the entire gene
	
	# Exon Scores
#	ExonScores=REIDS_Output[[which(names(REIDS_Output)==GeneID)]]$exonScore

	# Exon Level Data
#	ExonLevel=ExonLevelData[which(ExonLevelData$GeneID==GeneID),]
#	ExonLevel=ExonLevel[,-c(1,2)]
#	rownames(ExonLevel)=ExonLevel[,2]
#	
#	ExonLevel=ExonLevel[which(rownames(ExonLevel%in%rownames(ArrayScores)),]
#	message("The Exon Level Data is log2 transformed")
#	ExonLevel=log2(ExonLevel)
#	ExonLevel=ExonLevel[,unlist(Groups)]
	
	Exon_ExonID=ExonLevelData[which(ExonLevelData$ExonID==ExonID),]
	Exon_ExonID=Exon_ExonID[,-c(1,2)]
	Exon_ExonID=log2(Exon_ExonID)
	Exon_ExonID=Exon_ExonID[,unlist(groups)]
	
	
	# Gene Level Data
	Gene_GeneID=GeneLevelData[which(GeneLevelData$GeneID==GeneID),]
	Gene_GeneID[,-c(1)]=log2(Gene_GeneID[,-c(1)])
	Gene_GeneID=Gene_GeneID[,-c(1)]
	Gene_GeneID=Gene_GeneID[,unlist(groups)]
	
	group1=seq_along(groups$group1)
	group2=length(groups$group1)+seq_along(groups$group2)
	print(stats::t.test(x=Gene_GeneID[1,group1],y=Gene_GeneID[1,group2]))
	
	
	# Observed Probe intensities
	ProbeIntensities=Data[which(Data$ExonID%in%rownames(ArrayScores)),]
	ProbeIntensities_ExonID=ProbeIntensities[which(ProbeIntensities$ExonID==ExonID),]
	ProbeIntensities_ExonID=ProbeIntensities_ExonID[,-c(1,2)]
	ProbeIntensities_ExonID=ProbeIntensities_ExonID[,unlist(groups)]
	
	
	#plot 1 : gene level + exon level + probe intensities
	
	max_ylim=max(ProbeIntensities_ExonID,Gene_GeneID,Exon_ExonID)
	if(!is.null(location)){
		location1=paste(location,"_plot1",sep="")
	}
	else{
		location1=NULL
	}
	plottypein(plottype,location1)
	graphics::par(mar=c(6.5,4.5,2.5,4.5))
	graphics::plot(0,0,typ="n",xlab="",ylab=paste("Gene",GeneID," - Exon",ExonID,sep=" "),ylim=c(0,max_ylim+2),xlim=c(1,ncol(ProbeIntensities_ExonID)),xaxt="n")
	graphics::lines(x=c(1:ncol(ProbeIntensities_ExonID)),y=Gene_GeneID) #Gene level data...
	graphics::lines(x=c(1:ncol(ProbeIntensities_ExonID)),y=Exon_ExonID,col="blue") #Exon Level data
	for(i in 1:nrow(ProbeIntensities_ExonID)){
		graphics::points(x=c(1:ncol(ProbeIntensities_ExonID)),y=as.matrix(ProbeIntensities_ExonID[i,]),pch=19,col="blue")
	}
	graphics::axis(1,labels=colnames(Gene_GeneID),at=c(1:ncol(ProbeIntensities_ExonID)),las=2)
	plottypeout(plottype)

	#plot reserve : plot of cdf of exon scores

#	cdf=ecdf(ExonScores$X50.)
#	plot(sort(ExonScores$X50.),cdf(sort(ExonScores$X50.)),pch=19,xlim=c(0,1),ylim=c(0,1),xlab="Exon scores",ylab="Fn(x)")
#	points(sort(ExonScores$X50.[which(ExonScores$exon%in%ASExons_2736322$ExonID)]),cdf(sort(ExonScores$X50.[which(ExonScores$exon%in%ASExons_2736322$ExonID)])),pch=19,col="red")
	

	#plot 2 : density plot of the array scores
	if(!is.null(location)){
		location2=paste(location,"_plot2",sep="")
	}
	else{
		location2=NULL
	}
	plottypein(plottype,location2)
	graphics::plot(stats::density(as.matrix(ArrayScores[which(rownames(ArrayScores)==ExonID),unlist(groups)])),xlab=paste("Array Scores of Exon",ExonID,sep=" "),lty=6,,ylab=paste("Density of the Array Scores of Exon",ExonID,sep=" "),main="",lwd=2)
	graphics::points(x=as.matrix(ArrayScores[which(rownames(ArrayScores)==ExonID),groups$group1]),y=rep(0,length(groups$group1)),col="red",pch=8,cex=1.2)
	graphics::points(x=as.matrix(ArrayScores[which(rownames(ArrayScores)==ExonID),groups$group2]),y=rep(0,length(groups$group2)),col="blue",pch=8,cex=1.2)
	graphics::legend(x="topright",legend=c("Group 1","Group 2"),pch=19,col=c("red","blue"))
	plottypeout(plottype)
	
	#plot 3 : heatmaps of array and FIRMA scores
	FirmaScores=FIRMA_Data[which(FIRMA_Data$GeneID==GeneID),]
	FirmaScores=FirmaScores[which(FirmaScores$ExonID%in%rownames(ArrayScores)),]

	FirmaScores=FirmaScores[order(match(as.character(FirmaScores$ExonID),rownames(ArrayScores))),]
	
	
	#heatmap array scores
	if(!is.null(location)){
		location3a=paste(location,"_plot3a",sep="")
	}
	else{
		location3a=NULL
	}
	plottypein(plottype,location3a)
	mypalette<-brewer.pal(9,"PuBu")
	heatmapdata=ArrayScores
	colnames(heatmapdata)=colnames(heatmapdata)
	stats::heatmap(as.matrix(t(heatmapdata)),Rowv=NA,Colv=NA,col=mypalette,margins=c(6,6))
	plottypeout(plottype)
	
	
	#heatmap FIRMA scores
	if(!is.null(location)){
		location3b=paste(location,"_plot3b",sep="")
	}
	else{
		location3b=NULL
	}
	plottypein(plottype,location3b)
	mypalette<-brewer.pal(9,"PuBu")
	heatmapdata=FirmaScores[,-c(1,2)]
	rownames(heatmapdata)=FirmaScores[,2]
	colnames(heatmapdata)=colnames(heatmapdata)
	stats::heatmap(as.matrix(t(heatmapdata)),Rowv=NA,Colv=NA,col=mypalette,margins=c(6,6))
	plottypeout(plottype)
	
}

#Top gene 3762198  Top exon 3762266
#load("TestData/ColonCancerData.rda")
#load("TestData/ExonLevelSummarized_rma.RData")
#load("TestData/GeneLevelSummarized_rma.RData")
#load("TestData/ColonCancer_OutputREIDSModel.RData")
#load("TestData/Output_FIRMA.RData")
#
#groupT=c(1,14,16,18,20,3,5,7,9,11)
#groupN=c(12,15,17,19,2,4,6,8,10,13)
#groups=list(group1=groupN,group2=groupT)
#
#
#PlotFunction(GeneID="3762198",ExonID="3762266",Data=ColonCancerData,REIDS_Output=ColonCancer_OutputREIDSModel,GeneLevelData=GeneLevelSummarized_rma,ExonLevelData=ExonLevelSummarized_rma,FIRMA_Data=Output_FIRMA,groups=groups,plottype="pdf",location="Testdata/Gene3762198_Exon3762266")

# SI Index
#' "SpliceIndex"
#' 
#' The SpliceIndex function computes the ratio of the splice indices of defined groups. Further, it performs the SI algorithm if the length of the groups is two and the MiDAS algorithm is more groups are specified. Both algorithms are implemented as defined by Affymetrix.
#' @export
#' @param GeneData The microarray data summarized at gene level.
#' @param ExonData The microarray data summairzed at exon level.
#' @param InformativeExons A character vector of exon IDs. As for the REIDS model probesets are filtered out by I/NI calls model and later on exon score, the remaining exons can be specified here. Only these shall be considered in the FIRMA analysis to make the results between REIDS and FIRMA more comparable
#' @param groups The groups of interest in the data. Default two groups are specificied byut more can added as group3, group4,...
#' @param paired Logical. Are the groups paired? only used if two groups are present.
#' @param significancelevel If specified, filtering is conducted on the p-values.
#' @return A data frame wiith one line per exon. The columns conatin the gene ID, the exon ID, the ratio of the splice indices if two groups are present, a t- or F-statitic, a p-value and an adjusted p-value.
#' @details Given the gene level and exon level summarized data, the splice index method for the detection of alternative splicing is performed. The first step is to normalize the exon
#' data by taking the ratio with the gene level data. These values are referred to as the splice indices. If only two groups are specified, the ratio of their splice indices is taken as a measure for alternative splicing. The more the ratio deviates from zero, the more there is an indication of alternative splicing.
#' A t-test is conducted on the splice indices of the two groups to test their difference. If more than two groups are specified, an ANOVA model is fitted on the splice indices to discover with an F-test whether there is a difference between the groups somewhere. If a vector of informative exons is given 
#' to the function, only these are considered for the analysis. Finally, the p-values are adjusted for multiplicity and if a significance level is specified only the significant p-valuesare kept in the data frame.
SpliceIndex<-function(GeneData,ExonData,InformativeExons=NULL,groups=list(group1=NULL,group2=NULL),paired=FALSE,significancelevel=NULL){
	message("The gene and exon level data will be log2 transformed")
	
	Genedata[,-c(1)]=log2(Genedata[,-c(1)])
	Exondata[,-c(1,2)]=log2(Exondata[,-c(1,2)])
	
	if(!is.null(InformativeExons)){
		ExonData=ExonData[which(ExonData$ExonID%in%as.character(InformativeExons)),]
	}

	GeneID=intersect(GeneData$GeneID,ExonData$GeneID)  #special measure for the HTA Data : ExonLevel has Gene IDs from the mapping while GeneLevel has those of the ENSg cdf: different genes
	
	
	si<-function(ID,DG,DE,groups,paired){
		DataG=as.matrix(DG[,-c(1)])  #Gene Expression Data Level
		DataE=as.matrix(DE[,-c(1,2)]) #Exon Expression Data Level (exons coming from the same gene)
		
		DataE=DataE[!duplicated(rownames(DataE)), ,drop=FALSE]
		
		grouping=c(groups$group1,groups$group2)
		
		DataE=DataE[,grouping,drop=FALSE]
		DataG=DataG[,grouping,drop=FALSE]
		
		group1=seq_along(groups$group1)
		group2=length(group1)+seq_along(groups$group2)
		groups=list(group1=group1,group2=group2)
		
		NormData=sweep(DataE,2,DataG,"-")  #The Normalized Expression Data ( "-" because  log2 is already taken)
		
		SI_group1=apply(NormData[,group1,drop=FALSE],1,mean)
		
		SI_group2=apply(NormData[,group2,drop=FALSE],1,mean)
		
		SI=as.matrix(SI_group1-SI_group2,drop=FALSE)
		colnames(SI)="Splice Index"
		rownames(SI)=DE$ExonID
		
		stats<-function(i,Data,groups,paired){
			out1=stats::t.test(Data[i,groups$group1],y=Data[i,groups$group2],paired=paired)
			out2=cbind(out1$statistic,out1$p.value)
			
			return(out2)
		}
		
		tstats=t(sapply(c(1:nrow(NormData)),stats, Data=NormData, groups=groups,paired=paired)) 
		colnames(tstats)=c("t-statistic","p.value")
		
		out_temp=cbind(SI,tstats)
		
		return(out_temp)
	}
	

	midas<-function(ID,DG,DE,groups){
		DataG=as.matrix(DG[,-c(1)])  #Gene Expression Data Level
		DataE=as.matrix(DE[,-c(1,2)]) #Exon Expression Data Level (exons coming from the same gene)
		
		DataE=DataE[!duplicated(rownames(DataE)), ,drop=FALSE]
		
		grouping=as.vector(unlist(groups))
		
		DataE=DataE[,grouping,drop=FALSE]
		DataG=DataG[,grouping,drop=FALSE]
	
		
		groupstemp=list()
		for(i in 1:length(groups)){
			groupstemp[[i]]=rep(i,length(groups[[i]]))		
		}
		names(groupstemp)=names(groups) #in the assumptution that the names of groups are group1, group2, group3,...
		groups=as.factor(unlist(groupstemp))
		
		#normalized intensities : responses of the ANOVA model
		NormData=sweep(DataE,2,DataG,"-") 
	
		
		#anova model
		stats=t(sapply(1:nrow(DataE),function(i) summary(stats::aov(NormData[i,]~groups))[[1]][1,c(4,5)]))
		colnames(stats)=c("F-statistic","p.value")
		
		out_temp=stats
		return(out_temp)
		
	}
	
	
	if(length(groups)==2){
		#2 groups ; perform splice index analysis with t-test
		out=lapply(GeneID,function(x) si(ID=x,DG=GeneData[which(GeneData$GeneID==x),],DE=ExonData[which(ExonData$GeneID==x),],groups=groups,paired=paired)) #for every gene
		names(out)=GeneID
	}
	
	if(length(groups)>2){
		out=lapply(GeneID,function(x) midas(ID=x,DG=GeneData[which(GeneData$GeneID==x),],DE=ExonData[which(ExonData$GeneID==x),],groups=groups)) #for every gene
		names(out)=GeneID
	}
	
	out1<-do.call(rbind.data.frame, out)
	out1$adj.p.value=stats::p.adjust(out1$p.value,"fdr")
	
	fac=as.factor(sapply(c(1:nrow(out1)),function(x) strsplit(rownames(out1)[x],"[.]")[[1]][1]))
	out2<-split(out1,fac)
	
	replacerownames<-function(ID,names,Data){
		print(ID)
		rownames(Data)=names
		return(Data)
	}
	
	output=lapply(1:length(GeneID),function(x) replacerownames(ID=x,names=ExonData[which(ExonData$GeneID==GeneID[x]),2],Data=out2[[which(names(out2)==GeneID[x])]])) 
	names(output)=GeneID
	
	output2=do.call(rbind.data.frame, output)
	
	Names=strsplit(rownames(output2),"[.]")
	Names=do.call(rbind.data.frame, output)
	colnames(Names)=c("GeneID","ExonID")
	Names$GeneID=as.character(Names$GeneID)
	Names$ExonD=as.character(Names$ExonID)
	
	SI=cbind(GeneID=Names$GeneID,ExonID=Names$ExonID,output2)
	rownames(SI)=seq(1:nrow(SI))
	
	if(!is.null(significancelevel)){
		SI$adj.p.value=stats::p.adjust(SI$p.value,"fdr")
		SI=SI[which(SI$adj.p.value<0.05),]
	}

	return(SI)
	
}


if(getRversion() >= "2.15.1"){
	globalVariables(c("Arguments"))
}
