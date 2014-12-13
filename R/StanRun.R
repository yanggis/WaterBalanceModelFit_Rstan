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
# 	PD 	<- extract(SegFit, permuted=T)
# 	S 	<- summary(SegFit, probs=c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975))$summary
	
	# PRINT POSTERIOR SUMMARIES
# 	print(SegFit, digits_summary=3, pars=pars)
	return(SegFit)
}

StanRun <- function(segfit, PRdata, ParamUncert, PdataGMix, fit.pars, niter=1000, chains=2, parTF=F) {
	# PARAMETERS
# 	tWindow 		<- c(1,22) 				# [first,last]  either idx or yrs
# 	tType  			<- 'All' 					# ['All','WindowTcp','WindowYrs','WindowIdx']
# 	TcpYr 			<- 1978 					# first year of changed behavior
# 	dataType 		<- 'obsv' 				# ['obsv','synth']
# 	ixPtype 		<- 1							# Ptypes <- c('VICs','PRISM','VIC','GHCN')
# 	plotTF  		<- T
# 	parTF 			<- T
	if (parTF)
		chains=1
	
	# CALC DATA AND PROVIDE AS INPUT
	# Data 				<- GetData(ixWS=ixWS,dataType=dataType,tType=tType,tWindow=tWindow,TcpYr=TcpYr,ixPtype=ixPtype) # input data
	# PROVIDE PARAM SETTINGS AS INPUT
	# ParamUncert <- ParamUncertSettings(PET=Data$PET)
	# NOT USED
	# tParams 		<- tParamSettings(tType)
	# PROVIDE XDATAMIX AS INPUT
	# XdataMix 		<- MixNorm_Estimate(K=2, d=Data$Pall, plotTF=F) 	# gauss mix estimate of source data
	Dstan <- StanData(intvlType, PRdata, PdataGMix, ParamUncert)
	Init 	<- initFcn(Dstan, ParamUncert, intvlType, chains)
	
	# RUN STAN
	SegFit <- StanIterate(segfit=segfit, init=Init, chains=chains, 
													Dstan=Dstan, parTF=parTF, niter=niter)
	
	# PRINT POSTERIOR SUMMARIES
	print(SegFit, digits_summary=3, pars=fit.pars)
	
	return(list(SegFit=SegFit, Dstan=Dstan))
# 	dir.save <- Data$dirData
# 	if (!is.null(saveSubDir)){
# 		dir.save <- file.path(dir.save,saveSubDir)
# 		dir.create(dir.save, showWarnings = FALSE)
# 	}
# 	
# 	if (plotTF){
# 		PlotRvP(intvlType,Dstan,StanOut1$PD,StanOut1$S,ixWS,dir.save)
# 		# traceplot(SegFit,pars=pars)
# 	}
# 	
# 	# SAVE DATA
# 	if (saveTF){
# 		# SUMMARY TABLE
# 		fname <- file.path(dir.save,'RSTAN_1IntvlFitSummary.RData')
# 		S <- StanOut1$S
# 		save(list=c('S'),file=fname)
# 		
# 		# MODEL, POSTERIORS, PARAMETERS
# 		fname <- file.path(dir.save,'RSTAN_1IntvlFitData.RData')
# 		save(list=c('StanOut1','segfit','ParamUncert','tParams','XdataMix','Dstan','Init'),file=fname)
# 	}
# 	return(list(StanOut1=StanOut1,Data=Data,Dstan=Dstan,XdataMix=XdataMix))
}

