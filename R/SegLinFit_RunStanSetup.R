# RunStanSetup

# INITIALIZE STAN
out.StanInit <- InitStan(intvlType)
segfit <- out.StanInit$segfit

modelXY <- function(k1, k2, xmax, s1, s2){
	ModelVals.xy <- matrix(nrow=4,ncol=2,dimnames=list(c('0','K1','K2','xmax'), c('x','y'))) 
	ModelVals.xy[,'x'] <- c(0,k1,k2,xmax)
	ModelVals.xy[,'y'] <- modY(ModelVals.xy[,'x'],k1,k2,s1,s2)
	return(ModelVals.xy)
}

SaveOutput <- function(ixws,Out.1intvl,dirOut=dirOut){
# 	stadir <- StaDirs[ixws,1]
# 	ixStr <- gregexpr('Site',stadir); ixStr <- ixStr[[1]][1]
# 	staid <- substr(stadir,ixStr+4,ixStr+11)
	staid <- getSID2(ixws)
	fname <- paste('Data_',staid,'.Rdata',sep='')
	out.path <- file.path(dirOut, fname)
	
	# [year, P, R] data frame
	ypr <- Out.1intvl$Data$Data
	
	# Summary quantiles matrix
	S <- Out.1intvl$StanOut1$S
	s <- S[c('K1','K2','ETo','S1','S2'),c('mean','10%','50%','90%')]
	PDquants <- s
	
	# Model values for avg of posteriors
	avgType = '50%'
	k1 <- s['K1',avgType]; k2 <- s['K2',avgType] 
	eto<- s['ETo',avgType];
	s1 <- s['S1',avgType]; s2 <- s['S2',avgType]
	xmax <- max(ypr[,'xobs'])*1.2
	ModelVals.xy <- modelXY(k1,k2,xmax,s1,s2)
	# 	ModelVals.xy <- matrix(nrow=4,ncol=2,dimnames=list(c('0','K1','K2','xmax'), c('x','y'))) 
	# 	ModelVals.xy[,'x'] <- c(0,k1,k2,xmax)
	# 	ModelVals.xy[,'y'] <- modY(ModelVals.xy[,'x'],k1,k2,s1,s2)
	
	save(list=c('ypr','PDquants','ModelVals.xy'),file=out.path)
	return(list(staid=staid,ypr=ypr,PDquants=PDquants,ModelVals.xy=ModelVals.xy))
}

getWSix <- function(WStype){
	if (WStype=='nonalp_40'){  # Non-alp, low flow chk, 40 yrs+, record ends on/before 1990
		ixWS <- c(4, 19, 28, 42,44,  48,  52,  66,  79,  82,  89, 101, 103, 108, 115, 117, 119,
							124, 130, 146, 149, 151, 154, 155, 157, 166 ,169, 171, 173, 174, 182, 183, 194,  46, 112)
		m <- c(20, 30, 30, 20, 20, 30, 20, 20, 20, 40, 20, 20, 20, 20, 20, 30, 20, 20, 20, 20, 20, 20, 
					 20, 20, 20, 20, 20, 20, 20, 20, 40, 20, 20, 20, 30)
		m <- m + round(m/10)
	}
	
	if (WStype=='alp_40'){ 		# alpine, 40 yrs+, record ends on/before 1990
		ixWS <- c(9,38,63,78,95,110,122,123,139,143,165,178,186,187,189,196,199,204,211,152,215,59,148,216,51)
		m <- c(20,30,30,20,30,20,20,20,20,20,20,20,20,20,20,30,20,20,30,20,30,20,20,20,30)
		m <- m + round(m/10)
	}
	
	if (WStype=='all_20'){ 		# all ws 20 yrs+
		ixWS <- c(186,82,187,152,166,139,117,9,63,95,215,112,196,78,189,108,123,154,110,122,
							211,119,143,56,130,182,19,66,183,42,79,38,89,51,113,115,146,174,48,128,169,
							106,46,149,173,148,44,52,156,157,207,194,178,171,28,101,124,59,127,151,73,
							32,164,216,33,103,40,165,199,204,4,145,155,185,35,7,90,163,191,11,74,76,109,
							118,120,140,29,197,138,69,105,47,54,26,210,17,18,94,147,190,193,136,137,188,
							68,98,153,177,2,91,158,184,55,159,179,21,25,27,36,83,195,12,50,57,114,201,85,
							86,96,208,15,99,150,142,20,62,75,97,141,144,192,202,61,71,100,
							134,121,198,205,64,133,160,181,1,45,84,87,88,107,129,214,217,65,72,167,200)
	}
	
	if (WStype=='nonalp_20'){
		ixWS <- c(130,	66,	4	,68,	98,	96,	75,	129,	2,	69,	158,	25,	156,	91,	64,	83,	56,	20,	101,	18,	133,	45,	76,	85,
							191,	73,	90,	159,	57,	86,	71,	84,	42,	115,	169,	160,	108,	28,	21,	1,	19,	74,	184,	12,	114,	185,	
							100,	140,	99,	89,	44,	217,	109,	50,	17,	61,	149,	128,	171,	190,	52,	62,	33,	29,	88,	181,	173,	179,
							82,	163,	214,	134,	46,	32,	182,	154,	155,	141,	55,	11,	105,	150,	119,	103,	183,	94,	197,	117,	146,
							144,	54,	174,	7,	48,	79,	194,	87,	127,	97,	157)
	}
	return(ixWS)
}