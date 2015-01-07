# FUNCTIONS TO RUN STAN OVER SINGLE INTERVAL

# MODEL PARAMETER UNCERTAINTY BOUNDS
ParamUncertSettings <- function(PET) {
	# Uncertainty parameters for RStan priors
	sdX 		<- c(10, 0, 0.1) 			# [minSD,minCoef,maxCoef]
	sdY 		<- c(0.1, 0, 0.1)			# [minSD,minCoef,maxCoef] - minSD = 0 for yobs
	sdI			<- c(25, 500) 				# [min, max] intrinsic uncert stddev
	S2_lims <- c(0.6, 1.1) 				# (minS2, maxS2)
	PETlim 	<- 1.0*PET 						# PET estimate from CIMIS ETo product
	ParamUncert <- list(sdX=sdX, sdY=sdY, sdI=sdI, S2_lims=S2_lims, PETlim=PETlim)
}

StanInit <- function(stanFile) {
	standir  <- '~/Documents/CODE/R/WaterBalanceModelFit_Rstan/RStan/'
	stanPath <- paste(standir, stanFile, sep='')
	stanfit  <- stan(file=stanPath ,chains=0)
}

InitStan <- function(intvlType) {
	# INITIALIZE STAN MODELS
	if (intvlType=='1'){ 	# single interval 	
		# stanFile <- 'SegLinFit14_1e.stan'
		stanFile <- 'SegLinFit14_1e_k1_0.stan'
		if (etiTF)
			stanFile <- 'SegLinFit14_1e_k1_0_aPx.stan'
		if (!exists('segfit')) {					# only re-compile stan file if not already done
			segfit <- StanInit(stanFile) 
			cat('IGNORE THIS ERROR')
		}
		pars=c('S1','S2','K1','K2','ETo','sdI','ax','ay','A')
	}
	if (intvlType=='2'){ 	
		stanFile <- 'SegLinFit14_2b.stan'
		if (!exists('segfit2')){					# only re-compile stan file if not already done
			segfit <- StanInit(stanFile) 
			cat('IGNORE THIS ERROR')
		}
		pars <- c('Tcp','K1','dK1','S1','dS1','ETo','dETo','sdI','ax','ay','B')
	}
	return(list(segfit=segfit, pars=pars))
}

GetData <- function(ixWS, dataType, tType, tWindow, TcpYr, ixPtype=1) {
	if (dataType == 'obsv'){ 						# OBSERVED DATA
		Data <- ImportData(ixWS=ixWS, ixPtype=ixPtype, plotTF=F)
	} else if (dataType == 'synth'){
		Data <- SynthDataCalc(N=80) 	# SYNTHETIC DATA
	}
	return(Data)
}

StanData <- function(intvlType, PRdata, XdataMix, ParamUncert, cpParams=F, Tcp=F) {
	# OBSV DATA
	xobs <- PRdata$Pmm
	xobs[xobs<3] <- 3 		# need non-zero values for lognormal distribs in Stan
	yobs <- PRdata$Rmm
	yobs[yobs<3] <- 3
	N <- length(xobs)	
	
	# UNCERTAINTY PARAMS
	sdX 	 <- ParamUncert$sdX 			# [minSD,minCoef,maxCoef]
	sdY 	 <- ParamUncert$sdY				# [minSD,minCoef,maxCoef] - minSD = 0 for yobs
	sdI		 <- ParamUncert$sdI 			# [min, max] intrinsic uncert stddev
	s2lims <- ParamUncert$S2_lims		# (minS2, maxS2)
	PETlim <- ParamUncert$PETlim 		# PET is upper bound for ETo
	
	# STAN INPUTS
	if (intvlType == '1') {
		D <- list(N=N, xobs=xobs, yobs=yobs, ox=sdX, oy=sdY, s2lim=s2lims, oI=sdI, 
							chi_mu=XdataMix$MixMu, chi_sig=XdataMix$MixSD, 
							chi_theta=XdataMix$MixTheta,
							PETlim = PETlim)
	}
	if (etiTF)
		D$ETi <- PRdata$ETi
	
	# MUST UPDATE 2 INTVL REGRESSION TO MATCH 1 INTVL VERSION!!!
# 	if (intvlType == '2'){ 
# 		D <- list(N=N, xobs=xobs, yobs=yobs, ox=sdX, oy=sdY, s2lim=s2lims, oI=sdI, 
# 							chi_mu=XdataMix$MixMu, chi_sig=XdataMix$MixSD, 
# 							chi_theta=XdataMix$MixTheta,
# 							cpBuff = cpParams)
# 	}
# 	if (intvlType == '2Tcp'){
# 		D <- list(N=N, xobs=xobs, yobs=yobs, ox=sdX, oy=sdY, s2lim=s2lims, oI=sdI, 
# 							chi_mu=XdataMix$MixMu, chi_sig=XdataMix$MixSD, 
# 							chi_theta=XdataMix$MixTheta,
# 							cpBuff = cpParams, Tcp=Tcp)
# 	}
	return(D)
}

initList <- function(xobs, yobs, s2lim1, pet=3000, parInfo) {
	S2  <- 1
	k1 	<- runif(1, 1, pet/4)
	eto	<- runif(1, k1, pet*1/2)
	dk 	<- runif(1, (eto-k1), ((eto-k1)/(1-s2lim1))) # dk constrained by max s2
	initSD <- 10
	sdI <- rnorm(1, parInfo$sdI[2]/2, initSD)
	ax  <- parInfo$sdX[3]/2
	ay  <- parInfo$sdY[3]/2
	A   <- 5 		# inv logistic fcn shape param (as applicable)
	
	initList <- list(S2=S2, K1=k1, dK=dk, ETo=eto, nu=yobs, chi=xobs,
			 						 sdI=sdI, ax=ax, ay=ay, A=5 )
	
	if (etiTF) {
		initList$c <- 1	
	}
	return(initList)
}

initFcn <- function(Dstan, ParamUncert, intvlType, chains) {
	xobs <- Dstan$xobs
	yobs <- Dstan$yobs
	pet  <- Dstan$PETlim
	s2lim1 <- Dstan$s2lim[1]
	if (intvlType=='1') {
		init <- lapply(1:chains, 
									 function(id) initList(xobs, yobs, s2lim1, pet, ParamUncert)
									 )
	} else {
		init <- 'random'
	}
	return(init)
}
