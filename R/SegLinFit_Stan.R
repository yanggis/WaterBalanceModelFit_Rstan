# SegLinFitDATA fits a constrained tri-linear model to the data 'xobs' and 'yobs' 
# using RSTAN
# TC Moran, UC Berkeley 2014

# REQUIRED PACKAGES FOR THIS SCRIPT
library('rstan')
library('parallel')

StanInit <- function(stanFile) {
	standir  <- '~/Documents/CODE/R/RSTAN_SegLinFit/RSTAN/'
	stanPath <- paste(standir, stanFile, sep='')
	stanfit  <- stan(file=stanPath ,chains=0)
}

initList <- function(xobs, yobs, s2lim1, pet=3000, parInfo) {
	k1 	<- runif(1, 1, pet/4)
	eto	<- runif(1, k1, pet*1/2)
	dk 	<- runif(1, (eto-k1), ((eto-k1)/(1-s2lim1)))
	initSD <- 10
	
	list(S2=1, K1=k1, dK=dk, ETo=eto, nu=yobs, chi=xobs,
			 sdI=rnorm(1, parInfo$sdI[2]/2, initSD),
			 ax=parInfo$sdX[3]/2, ay=parInfo$sdY[3]/2, A=5 )
}

initFcn <- function(Dstan, ParamUncert, intvlType, chains) {
	xobs <- Dstan$xobs
	yobs <- Dstan$yobs
	pet  <- Dstan$PETlim
	s2lim1 <- Dstan$s2lim[1]
	if (intvlType=='1'){
		init <- lapply(1:chains, function(id) initList(xobs, yobs, s2lim1, pet, ParamUncert))
	} else {
		init <- 'random'
	}
	return(init)
}

StanData <- function(intvlType, Data, XdataMix, ParamUncert, cpParams=F, Tcp=F) {
	# OBSV DATA
	xobs <- Data$xobs
	xobs[xobs<3] <- 3 		# need non-zero values for lognormal distribs in Stan
	yobs <- Data$yobs
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
	if (intvlType == '2'){
		D <- list(N=N, xobs=xobs, yobs=yobs, ox=sdX, oy=sdY, s2lim=s2lims, oI=sdI, 
							chi_mu=XdataMix$MixMu, chi_sig=XdataMix$MixSD, 
							chi_theta=XdataMix$MixTheta,
							cpBuff = cpParams)
	}
	if (intvlType == '2Tcp'){
		D <- list(N=N, xobs=xobs, yobs=yobs, ox=sdX, oy=sdY, s2lim=s2lims, oI=sdI, 
							chi_mu=XdataMix$MixMu, chi_sig=XdataMix$MixSD, 
							chi_theta=XdataMix$MixTheta,
							cpBuff = cpParams, Tcp=Tcp)
	}
	return(D)
}

StanSegLinFit <- function(Dstan, segfit, Init, niter=1000, chains=4, parTF=F) {
	if(parTF) { # parallel RSTAN runs
		sflist <- mclapply(1:2, mc.cores = 2,      
											 function(i) stan(fit=segfit, data=Dstan, iter=niter, init=Init,
											 								  chains = 1, chain_id = i, refresh = -1)
											 )
		SegFit <- sflist2stanfit(sflist) 			# extract Stan output from parallel object
	} else {
		SegFit <- stan(fit=segfit, data=Dstan, iter=niter, chains = chains, init=Init)
	}
	return(SegFit)
}

GetData <- function(ixWS, dataType, tType, tWindow, TcpYr, ixPtype=1) {
	if (dataType == 'obsv'){ 						# OBSERVED DATA
		Data <- ImportData(ixWS=ixWS, ixPtype=ixPtype, plotTF=F)
	} else if (dataType == 'synth'){
		Data <- SynthDataCalc(N=80) 	# SYNTHETIC DATA
	}
	return(Data)
}

checkRhat <- function(S){
	rhat <- S[c('S1', 'K1', 'ETo', 'sdI'), 'Rhat']
	return(max(rhat)<1.1)
}

RunStan <- function(segfit, pars, init, chains, Dstan, parTF, niter) {	
	
	SegFit <- StanSegLinFit(Dstan, segfit, init, niter=niter, chains=chains, parTF=parTF)
	S <- summary(SegFit)
	rr <- 1
	while (!checkRhat(S$summary)) {
		print(paste('Poor convergence: retry number', rr))
		SegFit <- StanSegLinFit(Dstan, segfit, init, niter=niter, chains=chains, parTF=parTF)
		S <- summary(SegFit)
		rr <- rr + 1
		if (rr > 5) 
			break
	}	
	PD 	<- extract(SegFit, permuted=T)
	S 	<- summary(SegFit, probs=c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975))
	S 	<- S$summary

	# PRINT POSTERIOR SUMMARIES
	print(SegFit, digits_summary=3, pars=pars)
	
	return(list(SegFit=SegFit, PD=PD, S=S))
}