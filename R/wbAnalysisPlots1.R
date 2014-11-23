# wbAnalysisPlots1
library(ggplot2)

# SCATTER PLOT WITH X-AXIS RUG AND MODEL FIT
plotRvP <- function(data,sta.info,mod.vals,pd.quants,plotPrefs){
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
		geom_errorbarh(data=pd.quants[c(1,3),],aes(x=fifty,xmin=ten,xmax=ninety,y=y,height=height)) +
		geom_point(data=pd.quants['ETo',], aes(x=fifty,y=y),shape=18,size=4,colour=clr) + 
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

plotRvPspread <- function(pd,wbdata,sta.info,pd.quants,mod.xy.avg,scatTF=F){
	# HEAVISIDE FUNCTION
	H <- function(x){
		return((sign(x) + 1)/2)
	}
	# 3 segment linear model
	modY <- function(x,k1,k2,s1,s2){	
		modY <- s1*(x-k1)*H(x-k1) + (s2-s1)*(x-k2)*H(x-k2)	
	}
	# model xy values
	modelXY <- function(k1,k2,xmax,s1,s2){
		ModelVals.xy <- matrix(nrow=4,ncol=2,dimnames=list(c('0','K1','K2','xmax'), c('x','y'))) 
		ModelVals.xy[,'x'] <- c(0,k1,k2,xmax)
		ModelVals.xy[,'y'] <- modY(ModelVals.xy[,'x'],k1,k2,s1,s2)
		return(ModelVals.xy)
	}
	
	Pall<- as.data.frame(wbdata$Pall); names(Pall)<-'P'
	PET <- wbdata$pet
	ypr <- wbdata$Data
	
	# used for plot aesthetics
	xmax <- max(ypr[,'xobs'])*1.2
	y.max <- max(ypr$yobs) 		
	y.lim <- c(-y.max/15,1.2*max(ypr$yobs))
	x.lim <- c(0,1.1*max(wbdata$Pall))
	x.txt <- diff(range(ypr$xobs))/10
	y.txt <- diff(range(ypr$yobs))
# 	mod.xy.avg <- as.data.frame(modelXY(pD$pd.summ['K1','50%'],pD$pd.summ['K2','50%'],xmax,
# 																				pD$pd.summ['S1','50%'],pD$pd.summ['S2','50%']))
	clrS <- 'cadetblue4'
	clr <- 'deepskyblue4'
	clrP <- 'tomato3'
	
	# parameter plotting and text
	dy <- pmax(10,y.max/50)
	pd.quants$y <- c(0,0,-dy/2,0,0) 	# add y vals for plotting
	pd.quants$height <- array(dy,dim = c(5,1))
	names(pd.quants) <- c('mean','ten','fifty','ninety','y','height')
	par.txt<- paste0('S1 = ',format(pd.quants['S1','fifty'],digits=2),' [',format(pd.quants['S1','ten'],digits=2),',',format(pd.quants['S1','ninety'],digits=2),']\n',
									 'K1 = ',format(pd.quants['K1','fifty'],digits=2),'mm [',format(pd.quants['K1','ten'],digits=2),',',format(pd.quants['K1','ninety'],digits=2),']\n',
									 'ETo= ',format(pd.quants['ETo','fifty'],digits=2),'mm [',format(pd.quants['ETo','ten'],digits=2),',',format(pd.quants['ETo','ninety'],digits=2),']')
	
	pg <- ggplot(data=ypr,aes(x=xobs,y=yobs)) + 
		ylab('R (mm)') + xlab('P (mm)') + 
		ggtitle(paste0(sid,': ',sta.info$site_name,', ',sta.info$ws_area_km2,' km2, PET=',PET,'mm')) + 
		geom_rug(data=Pall,aes(x=P,y=0),side='b') +
		theme_bw(16)
	
	Ni <- length(pd$K1)
	for (ii in sample(1:Ni,round(Ni/5))){
		mod.vals <- as.data.frame(modelXY(pd$K1[ii],pd$K2[ii],xmax,pd$S1[ii],pd$S2[ii]))
		pg <- pg + geom_line(data = mod.vals,aes(x=x,y=y),colour=clrS,alpha=0.2,size=0.5)
	}
	pg <- pg +  geom_line(data = mod.xy.avg,aes(x=x,y=y),colour=clr,size=1.5) + 
		geom_point(data = mod.xy.avg,aes(x=x,y=y),shape=18,size=4,colour=clrP) + 
		geom_errorbarh(data=pd.quants[c(1,3),],aes(x=fifty,xmin=ten,xmax=ninety,y=y,height=height),colour=clrP) +
		geom_point(data=pd.quants['ETo',], aes(x=fifty,y=y),shape=18,size=4,colour=clrP) + 
		annotate('text',label=par.txt,x=x.txt,y=y.txt,hjust=0) +
		geom_point(size=3) +
		coord_fixed(ratio=1,xlim=x.lim, ylim=y.lim) 
	
	if(scatTF){
		plot(S1~K1,data=pd,pch=20,cex=0.6)
		abline(v=median(pd$K1),h=median(pd$S1),col='blue')
	}
	return(pg)
}

PlotPD_1Intvl <- function(PD){
	plot_name <- 'PD_1intvl.png'
	par(mfrow=c(2,2))
	plot(density(PD$K1))
	plot(density(PD$S2))
	plot(density(PD$ETo))
	plot(density(PD$S1))
	par(mfrow=c(1,1))
	dev.copy(png,plot_name,width=8,height=6,units="in",res=200)
	dev.off()
}
