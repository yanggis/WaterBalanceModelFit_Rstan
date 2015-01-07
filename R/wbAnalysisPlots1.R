# wbAnalysisPlots1
library(ggplot2)

# SCATTER PLOT WITH X-AXIS RUG AND MODEL FIT
plotRvP <- function(data, sta.info, mod.vals, pd.quants, plotPrefs){
	clr <- 'deepskyblue3'
	# P, R, PET data
	pet <- data$PET
	ypr <- as.data.frame(data$Data)
	Pall<- as.data.frame(data$Pall); names(Pall)<-'P'
	pscale <- 1/pet
	if (plotPrefs$scaleTF){
		ypr[,c('xobs','yobs')] <- ypr[,c('xobs','yobs')]*pscale
		Pall <- Pall*pscale
	}
	y.max <- max(ypr$yobs) 		# used for plot aesthetics
	y.lim <- c(-y.max/15,1.2*max(ypr$yobs))
	x.txt <- diff(range(ypr$xobs))/10
	y.txt <- diff(range(ypr$yobs))/3
	
	# Site info
	sid <- sta.info$site_num
	sta.name <- sta.info$site_name
	sta.area <- sta.info$ws_area_km2
	
	# MODEL AND ERROR BAR VALUES
	dy <- pmax(10,y.max/50)
	y.lab <- 'R (mm)'
	x.lab <- 'P (mm)'
	if (plotPrefs$scaleTF){
		mod.vals <- mod.vals*pscale
		pd.quants[c('K1','K2','ETo'),] <- pd.quants[c('K1','K2','ETo'),]*pscale
		dy <- y.max/50
		y.lab <- 'R/PET'
		x.lab <- 'P/PET'
	}
	
	pd.quants$y <- c(0,0,-dy/2,0,0) 	# add y vals for plotting
	pd.quants$height <- array(dy,dim = c(5,1))
	names(pd.quants) <- c('mean','ten','fifty','ninety','y','height')
	s1.txt<- paste0('S1 = ',format(pd.quants['S1','fifty'],digits=2),'\n[',format(pd.quants['S1','ten'],digits=2),',',format(pd.quants['S1','ninety'],digits=2),']')
	
	# scatter plot with rug 
	p <- ggplot(data=ypr,aes(x=xobs,y=yobs)) + 
		ylab(y.lab) + xlab(x.lab) + ggtitle(paste0(sid,': ',sta.name,', ',sta.area,' km2, PET=',pet,'mm')) + 
		geom_point(size=2) + 
		geom_rug(data=Pall,aes(x=P,y=0),side='b') + 
		theme_bw(16) + 
		geom_line(data = mod.vals,aes(x=x,y=y),colour=clr) + 
		geom_point(data = mod.vals,aes(x=x,y=y),shape=18,size=4,colour=clr) + 
		geom_errorbarh(data=pd.quants[c(1,3),], aes(x=fifty, xmin=ten, xmax=ninety, y=y, height=height)) +
		geom_point(data=pd.quants['ETo',], aes(x=fifty, y=y), shape=18, size=4, colour=clr) + 
		annotate('text',label=s1.txt,x=x.txt,y=y.txt,hjust=0) + 
		coord_fixed(ratio=1,xlim=c(0,1.1*max(Pall)), ylim=y.lim)  
	
	if(plotPrefs$plotTF) plot(p)
	if(plotPrefs$saveTF){
		out.dir <- file.path(plotPrefs$output.dir,'plots')
		plot.name <- file.path(out.dir,paste0(sid,'_RvP','.png'))
		if (plotPrefs$scaleTF) plot.name <- file.path(out.dir,'_scaled',paste0(sid,'_RvP_scaled','.png'))
		ggsave(plot.name,plot=p,width=11,height=8)
	}
	return(p)
}

plotRvPspread <- function(sid, pd, Summ, P.all, pet, PR,  scatTF=F) {
	# HEAVISIDE FUNCTION
	H <- function(x){
		return((sign(x) + 1)/2)
	}
	
	# 3 segment linear model
	modY <- function(x,k1,k2,s1,s2){	
		modY <- s1*(x-k1)*H(x-k1) + (s2-s1)*(x-k2)*H(x-k2)	
	}
	
	# model xy values
	modelXY <- function(k1, k2, xmax, s1, s2){
		ModelVals.xy <- matrix(nrow=4, ncol=2, dimnames=list(c('0','K1','K2','xmax'), c('x','y'))) 
		ModelVals.xy[ , 'x'] <- c(0, k1, k2, xmax)
		ModelVals.xy[ , 'y'] <- modY(ModelVals.xy[,'x'], k1, k2, s1, s2)
		return(ModelVals.xy)
	}
	
	# GET NEEDED QUANTILES FROM SegFit SUMMARY 
	pd.quants <- as.data.frame(Summ[c('S1', 'S2', 'K1', 'K2', 'ETo'), c('mean', '10%', '50%', '90%')])
	
	sta.info <- SqlGetWsInfo(sid=sid)
	x <- PR[ , 'Pmm']
	y <- PR[ , 'Rmm']
	p.all <- as.data.frame(P.all)
	names(p.all) <- 'P'
	
	# plot aesthetic settings
	xmax <- max(x)*1.2
	y.max <- max(y) 		
	y.lim <- c(-y.max/15, 1.2*y.max)
	x.lim <- c(0,1.1*max(P.all))
	x.txt <- diff(range(x))/10
	y.txt <- diff(range(y))
	
	clrS <- 'cadetblue4'
	clr <- 'deepskyblue4'
	clrB <- 'black'
	clrP <- 'tomato3'
	
	# parameter plotting and text
	dy <- pmax(10, y.max/50)
	pd.quants$y <- c(0, 0, -dy/2, 0, 0) 	# add y vals for plotting
	pd.quants$height <- array(dy, dim=c(5,1))
	names(pd.quants) <- c('mean', 'q10', 'q50', 'q90', 'y', 'height')
	par.txt <- paste0('S1 = ', format(pd.quants['S1', 'q50'], digits=2), ' [', format(pd.quants['S1','q10'], digits=2), ',', format(pd.quants['S1','q90'], digits=2), ']\n',
										'S2 = ', format(pd.quants['S2', 'q50'], digits=2), ' [', format(pd.quants['S2','q10'], digits=2), ',', format(pd.quants['S2','q90'], digits=2), ']\n',
										'K2 = ', format(pd.quants['K2', 'q50'], digits=2), 'mm [', format(pd.quants['K2','q10'], digits=2), ',', format(pd.quants['K2','q90'], digits=2), ']\n'
	)
	
	# INITIALIZE WITH RUG PLOT 
	pg <- ggplot(data=PR, aes(x=Pmm, y=Rmm)) + 
		ylab('R (mm)') + xlab('P (mm)') + 
		ggtitle(paste0(sid,': ',sta.info$site_name,', ',round(sta.info$ws_area_km2), ' km2, PET=', pet,'mm')) + 
		geom_rug(data=p.all, aes(x=P, y=0),side='b') +
		theme_bw(14)
	
	# ADD MANY REALIZATIONS OF MODEL ESTIMS
	Ni <- length(pd$K1)
	Nm <- min(Ni, 250)
	cat('Adding', Nm, 'model realizations to RvP plot, be patient...\n')
	for (ii in sample(1:Ni, Nm)){
		mod.vals <- as.data.frame(modelXY(pd$K1[ii], pd$K2[ii], xmax, pd$S1[ii], pd$S2[ii]))
		pg <- pg + geom_line(data = mod.vals, aes(x=x, y=y), colour=clrS, alpha=0.2, size=0.5)
	}
	
	# ADD MEAN FIT
	avg = 'q50'
	mod.xy.avg <- as.data.frame(modelXY(pd.quants['K1', avg], pd.quants['K2', avg], xmax,
																			pd.quants['S1', avg], pd.quants['S2', avg]))
	pg <- pg + geom_line(data=mod.xy.avg, aes(x=x,y=y), colour=clrB, size=1.5) + 
		geom_point(data=mod.xy.avg, aes(x=x, y=y), shape=18, size=4, colour=clrP)
	
	# ADD ERROR BARS AND ANNOTATION
	pg <- pg + geom_errorbarh(data=pd.quants[c('K1', 'K2', 'ETo'), ], aes(x=q50, xmin=q10, xmax=q90, y=y, height=height), colour=clrP) +
		geom_point(data=pd.quants['ETo', ], aes(x=q50, y=y), shape=18, size=4, colour=clrP) + 
		annotate('text', label=par.txt, x=x.txt, y=y.txt, hjust=0) +
		geom_point(size=3) +
		coord_fixed(ratio=1, xlim=x.lim, ylim=y.lim) 
	
	if(scatTF) {
		# S1 vs K1 scatterplot, used to assess correlations
		plot(S1~K1, data=pd, pch=20, cex=0.6)
		abline(v=median(pd$K1), h=median(pd$S1), col='blue')
	}
	return(pg)
}

plotRvPspreadETi <- function(sid, pd, Summ, P.all, pet, PR,  scatTF=F) {
	# HEAVISIDE FUNCTION
	H <- function(x){
		return((sign(x) + 1)/2)
	}
	
	# 3 segment linear model
	modY <- function(x,k1,k2,s1,s2){	
		modY <- s1*(x-k1)*H(x-k1) + (s2-s1)*(x-k2)*H(x-k2)	
	}
	
	# model xy values
	modelXY <- function(k1, k2, xmax, s1, s2){
		ModelVals.xy <- matrix(nrow=4, ncol=2, dimnames=list(c('0','K1','K2','xmax'), c('x','y'))) 
		ModelVals.xy[ , 'x'] <- c(0, k1, k2, xmax)
		ModelVals.xy[ , 'y'] <- modY(ModelVals.xy[,'x'], k1, k2, s1, s2)
		return(ModelVals.xy)
	}
	
	# GET NEEDED QUANTILES FROM SegFit SUMMARY 
	pd.quants <- as.data.frame(Summ[c('S1', 'S2', 'K1', 'K2', 'ETo'), c('mean', '10%', '50%', '90%')])
	
	sta.info <- SqlGetWsInfo(sid=sid)
	x <- PR[ , 'Pmm']
	y <- PR[ , 'Rmm']
	
	# Calculate Precip Excess Px = P - c*ETi
	c.mean <- Summ['c', 'mean']
	PR$Px <- PR$P - c.mean * PR$ETi
	
	p.all <- as.data.frame(P.all)
	names(p.all) <- 'P'
	
	# plot aesthetic settings
	xmax <- max(x)*1.2
	y.max <- max(y) 		
	y.lim <- c(-y.max/15, 1.2*y.max)
	x.lim <- c(0,1.1*max(P.all))
	x.txt <- diff(range(x))/10
	y.txt <- diff(range(y))
	
	clrS <- 'cadetblue4'
	clr <- 'deepskyblue4'
	clrB <- 'black'
	clrP <- 'tomato3'
	
	# parameter plotting and text
	dy <- pmax(10, y.max/50)
	pd.quants$y <- c(0, 0, -dy/2, 0, 0) 	# add y vals for plotting
	pd.quants$height <- array(dy, dim=c(5,1))
	names(pd.quants) <- c('mean', 'q10', 'q50', 'q90', 'y', 'height')
	par.txt <- paste0('S1 = ', format(pd.quants['S1', 'q50'], digits=2), ' [', format(pd.quants['S1','q10'], digits=2), ',', format(pd.quants['S1','q90'], digits=2), ']\n',
										'K1 = ', format(pd.quants['K1', 'q50'], digits=2), 'mm [', format(pd.quants['K1','q10'], digits=2), ',', format(pd.quants['K1','q90'], digits=2), ']\n',
										'ETo= ', format(pd.quants['ETo', 'q50'], digits=2), 'mm [', format(pd.quants['ETo', 'q10'], digits=2), ',', format(pd.quants['ETo', 'q90'], digits=2), ']'
	)
	
	# INITIALIZE WITH RUG PLOT 
	pg <- ggplot(data=PR, aes(x=Px, y=Rmm)) + 
		ylab('R (mm)') + xlab('Px (mm)') + 
		ggtitle(paste0(sid,': ',sta.info$site_name,', ',round(sta.info$ws_area_km2), ' km2, PET=', pet,'mm')) + 
		geom_rug(data=p.all, aes(x=P, y=0),side='b') +
		theme_bw(14)
	
	# ADD MANY REALIZATIONS OF MODEL ESTIMS
	Ni <- length(pd$K1)
	Nm <- min(Ni, 250)
	cat('Adding', Nm, 'model realizations to RvP plot, be patient...\n')
	for (ii in sample(1:Ni, Nm)){
		mod.vals <- as.data.frame(modelXY(pd$K1[ii], pd$K2[ii], xmax, pd$S1[ii], pd$S2[ii]))
		pg <- pg + geom_line(data = mod.vals, aes(x=x, y=y), colour=clrS, alpha=0.2, size=0.5)
	}
	
	# ADD MEAN FIT
	avg = 'q50'
	mod.xy.avg <- as.data.frame(modelXY(pd.quants['K1', avg], pd.quants['K2', avg], xmax,
																			pd.quants['S1', avg], pd.quants['S2', avg]))
	pg <- pg + geom_line(data=mod.xy.avg, aes(x=x,y=y), colour=clrB, size=1.5) + 
		geom_point(data=mod.xy.avg, aes(x=x, y=y), shape=18, size=4, colour=clrP)
	
	# ADD ERROR BARS AND ANNOTATION
	pg <- pg + geom_errorbarh(data=pd.quants[c('K1', 'K2', 'ETo'), ], aes(x=q50, xmin=q10, xmax=q90, y=y, height=height), colour=clrP) +
		geom_point(data=pd.quants['ETo', ], aes(x=q50, y=y), shape=18, size=4, colour=clrP) + 
		annotate('text', label=par.txt, x=x.txt, y=y.txt, hjust=0) +
		geom_point(size=3) +
		coord_fixed(ratio=1, xlim=x.lim, ylim=y.lim) 
	
	if(scatTF) {
		# S1 vs K1 scatterplot, used to assess correlations
		plot(S1~K1, data=pd, pch=20, cex=0.6)
		abline(v=median(pd$K1), h=median(pd$S1), col='blue')
	}
	return(pg)
}

PlotPD_1Intvl <- function(PD, path.save=NULL, saveTF=T){
	par(mfrow=c(2,2))
	plot(density(PD$K1))
	plot(density(PD$S2))
	plot(density(PD$ETo))
	plot(density(PD$S1))
	par(mfrow=c(1,1))
	if (saveTF) {
		dev.copy(png,path.save, width=8, height=6, units="in", res=72)
		dev.off()
	}
}

PlotPD_1Intvl2 <- function(PD, path.save=NULL, saveTF=T){
	par(mfrow=c(2,2))
	plot(density(PD$K1))
	plot(density(PD$K2))
	plot(density(PD$S1))
	plot(density(PD$S2))
	par(mfrow=c(1,1))
	if (saveTF) {
		dev.copy(png, path.save, width=8, height=6, units="in", res=72)
		dev.off()
	}
}