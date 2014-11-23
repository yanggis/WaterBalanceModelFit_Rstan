dir_master <- "~/Desktop/Catchment Analysis 2011/AA_CA_Catchments_Master/GAGESII_CATCHMENTS_CA/GAGESII_CATCHMENTS_219"

# IMPORT LIST OF DIRECTORIES FOR STATION DATA
ImportStaDirs <- function() {
	sta_file <- 'FILTERED_CATCHMENTS_219.csv'
	StaDirs  <- read.csv(file.path(dir_master, sta_file), header=F)
	StaDirs  <- data.frame(lapply(StaDirs, as.character), stringsAsFactors=FALSE)
}

ImportData <- function(ixWS=52, ixPtype=1, plotTF=F){
	Ptypes <- c('VICs', 'PRISM', 'VIC', 'GHCN')
	ptype <- Ptypes[ixPtype]
	
	dirData <- paste(dir_master, StaDirs$V1[ixWS], 'PARAM_FIT', 
									 paste('PDATA', ptype, sep='_'), sep='/') 
	tprFile <- file.path(dirData, 'tpr.csv')
	
	# IMPORT P, R DATA 
	TPR 	<- read.csv(tprFile, header=F) 		# no header
	xchk 	<- which(TPR[,2] < 0) 			# check for invalid data
	ychk 	<- which(TPR[,3] < 0) 	
	chk0 	<- unique(c(xchk, ychk))
	if (length(chk0) > 0)
		TPR <- TPR[-chk0, ]
	Data <- as.data.frame(TPR)
	names(Data) <- c('yrobs', 'xobs', 'yobs')
	
	# IMPORT LONG-TERM P RECORD FOR CHI DISTRIB ESTIM
	Pall <- read.csv(file.path(dirData, 'p_all.csv'), header=F)
	Pall <- Pall$V2
	
	# IMPORT PET
	pet <- read.csv(file.path(dirData, '../pet_cimis_mean.csv'), header=F)
	pet <- round(pet[[1]])
	
	if (plotTF){
		plot(yobs~xobs, data=Data, pch=16, col='indianred3')
	}
	
	return(list(Data=Data, Pall=Pall, PET=pet, dirData=dirData))
}
