# FUNCTIONS TO RUN STAN OVER SINGLE INTERVAL

if (!exists('StaDirs')){ 			# station directory listing
	StaDirs <- ImportStaDirs()
}

InitStan <- function(intvlType='1'){
	# INITIALIZE STAN MODELS
	if (intvlType == '1'){ 	
		stanFile <- 'SegLinFit14_1e.stan'
		if (!exists('segfit')){					# only re-compile stan file if not already done
			segfit <- StanInit(stanFile) 
		}
		pars=c('S1','S2','K1','K2','ETo','sdI','ax','ay','A')
	}
	if (intvlType == '2'){ 	
		stanFile <- 'SegLinFit14_2b.stan'
		if (!exists('segfit2')){					# only re-compile stan file if not already done
			segfit <- StanInit(stanFile) 
		}
		pars <- c('Tcp','K1','dK1','S1','dS1','ETo','dETo','sdI','ax','ay','B')
	}
	return(list(segfit=segfit,pars=pars))
}

GoStanGo <- function(segfit,pars,ixWS,niter=1000,saveTF=T,saveSubDir=NULL){
	# PARAMETERS
	tWindow 		<- c(1,22) 				# [first,last]  either idx or yrs
	tType  			<- 'All' 					# ['All','WindowTcp','WindowYrs','WindowIdx']
	TcpYr 			<- 1978 					# first year of changed behavior
	dataType 		<- 'obsv' 				# ['obsv','synth']
	ixPtype 		<- 1							# Ptypes <- c('VICs','PRISM','VIC','GHCN')
	plotTF  		<- T
	parTF 			<- T
	if (parTF){	chains<-1 } else{ chains<-2 }
	
	Data 				<- GetData(ixWS=ixWS,dataType=dataType,tType=tType,tWindow=tWindow,TcpYr=TcpYr,ixPtype=ixPtype) # input data
	ParamUncert <- ParamUncertSettings(PET=Data$PET)
	tParams 		<- tParamSettings(tType)
	XdataMix 		<- MixNorm_Estimate(K=2, d=Data$Pall, plotTF=F) 	# gauss mix estimate of source data
	Dstan 			<- StanData(intvlType, Data$Data, XdataMix, ParamUncert)
	Init 				<- initFcn(Dstan, ParamUncert, intvlType, chains)
	
	# RUN STAN
	StanOut1 <- RunStan(segfit=segfit, pars=pars, init=Init, chains=chains, 
											Dstan=Dstan, parTF=parTF, niter=niter)
	
	dir.save <- Data$dirData
	if (!is.null(saveSubDir)){
		dir.save <- file.path(dir.save,saveSubDir)
		dir.create(dir.save, showWarnings = FALSE)
	}
	
	if (plotTF){
		PlotRvP(intvlType,Dstan,StanOut1$PD,StanOut1$S,ixWS,dir.save)
		# traceplot(SegFit,pars=pars)
	}
	
	# SAVE DATA
	if (saveTF){
		# SUMMARY TABLE
		fname <- file.path(dir.save,'RSTAN_1IntvlFitSummary.RData')
		S <- StanOut1$S
		save(list=c('S'),file=fname)
		
		# MODEL, POSTERIORS, PARAMETERS
		fname <- file.path(dir.save,'RSTAN_1IntvlFitData.RData')
		save(list=c('StanOut1','segfit','ParamUncert','tParams','XdataMix','Dstan','Init'),file=fname)
	}
	return(list(StanOut1=StanOut1,Data=Data,Dstan=Dstan,XdataMix=XdataMix))
}
