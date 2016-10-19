
# Necessary Pacackages
library(REIDS)

# Reading the data in form the line indexing file created with the cluster function
args <- commandArgs(TRUE)
file_name <- as.character(args[1])
file_pos <- as.numeric(args[2])
line_length <- as.numeric(args[3])

conn <- file(file_name, 'rb')
current.pos <- seek(conn, where = file_pos, origin = 'start')
data <- readBin(conn, 'raw', n = line_length)
s <- paste(sapply(data, rawToChar), collapse='')

t=strsplit(s,"\",\"")

d=unlist(t)

d[1]=substr(d[1],start=2,stop=nchar(d[1]))
d[5]=substr(d[5],start=1,stop=nchar(d[5])-1)


if(d[1]!="geneID"){
	
	geneID=as.character(d[1])
	exonID=as.character(unlist(strsplit(d[2],",")))
	lengthexons=as.integer(unlist(strsplit(d[3],",")))
	
	npersample=sum(lengthexons)
	
	allsamples=d[4]
	samples=as.numeric(unlist(strsplit(allsamples,",")))
	samplenames=as.character(unlist(strsplit(d[5],",")))
	nsamples=length(samplenames)
	
	splitsamples<-function(x,samples,npersample){
		start=1+npersample*(x-1)
		end=npersample*x
		values=samples[start:end]
		return(values)
	}
	
	samplevalues=lapply(c(1:nsamples),function(i) splitsamples(i,samples,npersample) )
	TempData=rbindlist(list(samplevalues))
	setnames(TempData,colnames(TempData),samplenames)
	
	geneData=data.frame(geneID=rep(geneID,npersample),exonID=rep(exonID,lengthexons))
	geneData=data.frame(lapply(geneData, as.character), stringsAsFactors=FALSE)
	geneData=cbind(geneData,TempData)
	
	REIDS_Gene=REIDS_ClusterVersion(geneData=geneData,nsim=5000,geneID=geneID,informativeCalls=TRUE,alpha=0.5)
	assign(paste("REIDS_Gene",geneID,sep="_"),REIDS_Gene)
	eval(parse(text=paste("save(REIDS_Gene", geneID, ", file=\"/.../REIDS_Gene", geneID, ".RData\")", sep=""))) ## FILL IN  A LOCATION TO SAVE THE DATA
	
}