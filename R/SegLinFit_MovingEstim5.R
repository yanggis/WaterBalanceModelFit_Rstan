# RUN PARAM ESTIM FOR MULTIPLE INTERVALS

DropData <- function(data,PD,dtype='lp',yrIgnore=1977){
	if (dtype=='lp'){
		# Ignore least likely data point
		lp_nu 		<- colMeans(PD$nu_lp) 	# log prob of 'nu' parameters
		ix_ignore <- which.min(lp_nu)
	}
	
	if (dtype=='yr'){
		ix_ignore <- which(data$yrobs==yrIgnore)
	}
	if (length(ix_ignore)>0){
		data <- data[-ix_ignore,]
	}
	return(data)
}

MovingFit <- function(Out.1intvl,pars,niter,m=22){ 	# m = regression interval span
	S1 	<- Out.1intvl$StanOut1$S
	PD1 <- Out.1intvl$StanOut1$PD
	D1  <- Out.1intvl$StanOut1
	Data<- Out.1intvl$Data
	XdataMix <- Out.1intvl$XdataMix
	
	dirData <- paste(Data$dirData,'ANALYSIS_MovingFit',sep='/')
	dir.create(dirData, showWarnings = FALSE)
	wd <- setwd(dirData)
	
	dy <- 1 		# time step
	IntvlParams <- list(m=m,dy=dy)
	ix <- 1:m
	DATA 		<- Data$Data
	N    		<- length(DATA$xobs)
	M 			<- floor((N-m)/dy)+1 		# number of intervals
	dim.summ<- c(M,dim(S1)[2]) 			# columns correspond to StanFit summary columns, set in RunStan
	col.summ<- colnames(S1)
	plotTF 	<- F
	parTF  	<- T
	if (parTF){chains=1} else{chains=2}
	yr <- 0; 		# initialize year vector
	
	ParamUncert <- ParamUncertSettings(Data$PET)
	
	# Posterior Distributions
	PDlist <- list(j1=list(),j2=list())
	
	# Posterior Summaries
	pnames <- list('ETo','K1','dK','K2','S2','S1','sdI','ax','ay','lp__','A')
	Slist <- lapply(1:length(pnames),function(n)list(j1=array(dim=dim.summ,dimnames=list(NULL,col.summ)),j2=array(dim=dim.summ,dimnames=list(NULL,col.summ))))
	names(Slist) <- pnames
	
	# Checks for informed intervals
	dataChk <- lapply(1:2,function(n)list(j1=vector(length=M),j2=vector(length=M)))
	names(dataChk) <- list('K1','K2')
	
	k <- 1  									# set to 2 to run 2nd iteration removing some value
	for (ii in 1:M){
		data <- DATA[ix,]
		yr[ii] <- data$yrobs[m]
		for (jj in 1:k){
			if (jj==2) data<-DropData(data,PD,dtype='lp')
			
			print(data$yrobs)
			
			Dstan <- StanData(intvlType='1', data, XdataMix, ParamUncert)
			if (Dstan$oI[1] > max(data$yobs)){ 	# Adjust sdI for very dry watersheds
				Dstan$oI[1] 	<-max(data$yobs)/2
			} 
			
			Init 				<- initFcn(Dstan,ParamUncert,intvlType='1',chains)
			ModelOutput <- RunStan(segfit,pars,Init,chains,Dstan,parTF,niter) 				# Run RStan
			
			S  <- ModelOutput$S
			PD <- ModelOutput$PD
			PDlist[[jj]][[ii]] <- PD
			
			# PLOT: overlay this inverval on all data and mean
			if (jj==1){
				PlotRvP_1Intvl_AddSubIntvl(DATA,S1,data,S,Title=paste(yr[ii],' v',jj,sep=''))
				plot_name <- paste('SegLinFitMoving_',yr[ii],'_m',m,'_dy',dy,'_v',jj,'.png',sep='')
				dev.copy(png,plot_name,width=8,height=6,units="in",res=200)
				dev.off()
			}
			
			# PARAMETER SUMMARIES
			for (n in 1:length(pnames)){
				pname=pnames[[n]] 
				Slist[[pname]][[jj]][ii,]=S[pname,]
			}
			
			# MCMC ESTS OF 'TRUE' VALUES
			xT <- colMeans(PD$chi) 				# mean value of 'true' estimated x
			yT <- colMeans(PD$nu) 				# mean value of 'true' estimated y
			
			# K2 FILTER: FLAG ESTS WITH NO/FEW CHI > K2
			k2lo <- Slist[['K2']][[jj]][ii,'25%']
			if (sum(xT > k2lo) > 2){TF<-T}	else{TF<-F}
			dataChk$K2[[jj]][ii] <- TF
			
			# K1 FILTER: FLAG ESTS WITH NO/FEW CHI < K2
			if (sum(xT < k2lo) > 2){TF<-T}	else{TF<-F}
			dataChk$K1[[jj]][ii] <- TF	
		}
		ix <- ix+dy
	}
	
	# PLOTS
	Qpd <- matrix(nrow=M,ncol=3)
	Mpd <- matrix(nrow=M,ncol=1)
	Vars <- c('K1','S1','ETo')
	plotProbs <- c(0.1,0.5,0.9)
	yr1Fit <- 1950 									# regression on data since 1950
	idxRegr <- yr %in% yr1Fit:2003
	LinFits <- list()
	
	for (ii in 1:length(Vars)){
		Var <- Vars[ii]
		for (n in 1:M){
			dPD <- PDlist[[k]][[n]][[Var]] - PD1[[Var]] 	# EACH INTERVAL MINUS LONG-TERM INTERVAL
			Qpd[n,] <- quantile(dPD,probs=plotProbs)
			Mpd[n,] <- median(dPD)
		}
		plot(yr,Mpd,ylim=range(Qpd),xlim=c(yr[1]-m-1,yr[M]),xlab='Water Year',ylab='Change vs Long-term Mean')
		segments(yr-m,Mpd,yr,Mpd,col='black',lty=3)
		points(yr-m,Mpd,pch=16,col='black',cex=0.5)
		
		if (Var=='K1' | Var=='S1')	{chkP<-dataChk$K1[[k]]}
		if (Var=='ETo')							{chkP<-dataChk$K2[[k]]}
		arrows(yr[chkP], Qpd[,1][chkP], yr[chkP], Qpd[,3][chkP], length=0.05, angle=90, code=3)
		arrows(yr[!chkP], Qpd[,1][!chkP], yr[!chkP], Qpd[,3][!chkP], length=0.05, angle=90, code=3,lty=2)
		title(paste(Vars[ii],' v',k,' m=',m,' dy=',dy, sep=''))
		abline(h=0)
		grid()
		
		# REGRESSION OF MEDIAN VALUES
		linfit <- lm(Mpd[idxRegr]~yr[idxRegr])
		abline(linfit,col='blue')
		legend('topleft',paste('slope=',format(linfit$coefficients[2],digits=2),sep=''),lty=1,col='blue')
		LinFits[[Var]] <- linfit
		
		# Save plot
		plot_name <- paste('MovingFitRegress_',Var,'.png',sep='')
		dev.copy(png,plot_name,width=8,height=6,units="in",res=200)
		dev.off()
	}
	
	# SAVE DATA
	fname <- paste('RSTAN_MOVING_FIT_m',m,'_dy',dy,'.RData',sep='')
	save(list=c('Slist','LinFits','dataChk','yr','IntvlParams'),file=fname)
	setwd(wd)
	return(list(PDlist=PDlist,Slist=Slist,LinFits=LinFits,dataChk=dataChk,yr=yr,IntvlParams=IntvlParams))	
}
