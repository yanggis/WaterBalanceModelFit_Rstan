# Relatively generic functions to execute Rstan for segmented linear regression

library('rstan')
library('parallel')

StanSegLinFit <- function(Dstan, segfit, Init, niter=1000, chains=4, parTF=F) {
	# Runs Stan model in either parallel (parTF=T) or serial mode
	if(parTF) { # parallel RSTAN runs
		sflist <- mclapply(1:2, mc.cores = 2,      
											 function(i) stan(fit=segfit, data=Dstan, iter=niter, init=Init,
											 								 chains=1, chain_id=i, refresh=-1)
		)
		SegFit <- sflist2stanfit(sflist) 			# extract Stan output from parallel object
	} else {
		SegFit <- stan(fit=segfit, data=Dstan, iter=niter, chains=chains, init=Init)
	}
	return(SegFit)
}

checkRhat <- function(S, rhat_max=1.1){
	# check rhat below specified threshold rhat_max
	rhat <- S[c('S1', 'K1', 'ETo', 'sdI'), 'Rhat']
	return(max(rhat)<rhat_max)
}

StanIterate <- function(segfit, init, chains, Dstan, parTF, niter) {	
	# Iterates Stan runs until satisfactory Rhat metrics or Nmax runs
	Nmax <- 5
	rr <- 1
	rhatChk <- F
	while (!rhatChk && rr<=Nmax) {
		cat('Running Stan: iteration', rr, '\n')
		SegFit <- StanSegLinFit(Dstan, segfit, init, niter=niter, chains=chains, parTF=parTF)
		rr <- rr + 1
		rhatChk <- checkRhat(summary(SegFit)$summary)
	}	
	return(SegFit)
}

StanRun <- function(segfit, PRdata, ParamUncert, PdataGMix, fit.pars, niter=1000, chains=2, parTF=F) {
	# PARAMETERS
	if (parTF)
		chains=1
	
	Dstan <- StanData(intvlType, PRdata, PdataGMix, ParamUncert)
	Init 	<- initFcn(Dstan, ParamUncert, intvlType, chains)
	
	# RUN STAN
	SegFit <- StanIterate(segfit=segfit, init=Init, chains=chains, 
													Dstan=Dstan, parTF=parTF, niter=niter)
	
	# PRINT POSTERIOR SUMMARIES
	print(SegFit, digits_summary=3, pars=fit.pars)
	
	return(list(SegFit=SegFit, Dstan=Dstan))
}

