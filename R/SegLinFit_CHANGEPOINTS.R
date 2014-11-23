# SegLinFit_CHANGEPOINTS
# used infrequently as of Sep 2014

# CHANGEPOINT DETECTION SETTINGS
tParamSettings <- function(tType='All'){
	if (tType=='All'){
		cpParams <- c(12, 12) 			# min number of years in first/last intvl of cp detection
	}
	if (tType=='WindowTcp'){
		cpParams <- c(15, 15)
	}
	if (tType=='WindowYrs'){
		cpParams <- c(12, 12)
	}
	if (tType=='WindowIdx'){
		cpParams <- c(12, 12)
	}
	return(cpParams)
}

DataWindowTcp <- function(Data, TcpYr, cpParams) {
	cpBuff1 <- cpParams[1]
	cpBuff2 <- cpParams[2]
	yrs 		<- Data$yrobs
	ixTcp 	<- min(which(yrs >= TcpYr))
	if (ixTcp==Inf)
		return(F)
	
	ix1 <- ixTcp - cpBuff1 - 4 + 1 	# eg if Tcp = 1978, start data in 1960, TcpWindow extends back 4 years
	if (ix1 < 1)
		return(F)
	
	ixEnd <- ixTcp + cpBuff2 + 6
	if (ixEnd > length(yrs))
		return(F)
	
	ixYrs <- ix1:ixEnd 					# total interval = 40 years
	Data <- Data[ixYrs, ]
}

DataWindowYrs <- function(Data, tWindow) {
	yr1 <- tWindow[1]
	yr2 <- tWindow[2]
	yrs <- Data$yrobs
	ix1 <- min(which(yrs >= yr1))
	if (ix1==Inf)
		return(F)
	ix2 <- min(which(yrs >= yr2))
	if (ix2==Inf)
		return(F)
	if (ix2>length(yrs))
		return(F)
	ixYrs <- ix1:ix2
	Data <- Data[ixYrs, ]
}

DataWindowIdx <- function(Data,tWindow){
	Data <- Data[tWindow[1]:tWindow[2], ]
}

# **** GETDATA FUNCTIONS
# if (tType=='WindowTcp'){
# 	Data$Data <- DataWindowTcp(Data$Data,TcpYr,tParams)
# }
# if (tType=='WindowYrs'){
# 	Data$Data <- DataWindowYrs(Data$Data,tWindow)
# }
# if (tType=='WindowIdx'){
# 	Data$Data <- DataWindowIdx(Data$Data,tWindow)
# }