# SegLinFit_MultIntvl_Plots_# characterizes and plots parameter estimates for
# multiple intervals as output from SegLinFit_MultIntvl_#
#
# TC Moran, UC Berkeley 2014


# HEAVISIDE FUNCTION - also in SYNTHDATA
H <- function(x) {
	return((sign(x) + 1)/2)
}

# 3 segment linear model - also in SYNTHDATA
modY <- function(x, k1, k2, s1, s2) {	
	modY <- s1*(x-k1)*H(x-k1) + (s2-s1)*(x-k2)*H(x-k2)	
}

# PlotRvP <- function(intvlType, Dstan, PD, S, ixWS, dirData) {
# 	wd <- setwd(dirData)
# 	# PLOTS 
# 	if (intvlType == '2') {
# 		PlotPD_2Intvl(PD) 												# PD densities
# 		PlotRvP_2Intvl(Data$Data, S, Title=ixWS) 	# model
# 	}
# 	if (intvlType == '1'){
# 		PlotPD_1Intvl(PD)
# 		PlotRvP_1Intvl(D=Dstan, PD=PD, S=S, Title=ixWS)	
# 	}
# 	setwd(wd)
# }

PlotRvP_1Intvl <- function(D, PD, S, Title='', avg_type='50%') {
	plot_name <- 'RvP_SegFit_1intvl.png'
	clrs <- c('blue3','slategrey')
	x 	<- D$xobs
	y 	<- D$yobs
	xT 	<- colMeans(PD$chi)
	yT 	<- colMeans(PD$nu)
	# PARAMETER ESTIMATES
	k1 	<- S['K1', avg_type]
	k2 	<- S['K2', avg_type] 
	s1 	<- S['S1', avg_type] 
	eto	<- S['ETo',avg_type]
	s2 	<- S['S2', avg_type]
	# PLOT LIMITS
	Xmax <- 1.2 * max(x)
	Xlim <- c(0, Xmax)
	Ymax <- 1.2 * max(y)
	Ymin <- -40
	Ylim <- c(Ymin, Ymax)
	
	# PLOT DATA
	matplot(x,y,pch=19,col=clrs[1],xlab='P (mm)',ylab='R (mm)',xlim=Xlim,ylim=Ylim) 	# obsv data
	matplot(xT,yT,pch=1,col=clrs[1],add=T,cex=0.7) 																		# truth estimate
	
	# lines connecting obsv to truth estim
	for (n in 1:length(x)){
		segments(x[n],y[n],xT[n],yT[n],col=clrs[1],lwd=1,lty=3)
	}
	
	# Model Segments
	lx <- c(0,k1,k2,Xmax) 						
	ly <- modY(lx,k1,k2,s1,s2)
	segments(lx[1:3],ly[1:3],lx[2:4],ly[2:4],col=clrs[2],lwd=1)
	matplot(lx,ly,pch=18,col=clrs[2],add=T,cex=1.2)
	
	# Parameter estimates
	# S1
	dx <- diff(Xlim)/50
	dy <- diff(Ylim)/4
	y1 <- diff(Ylim)/4
	y2 <- y1+dy
	# scaled subplot
	segments(0,y1,dx,y1)
	segments(0,y2,dx,y2)
	text(1.2*dx,y1,'0',cex=0.8)
	text(1.2*dx,y2,'1',cex=0.8)
	text(0.5*dx,y2+dy*0.1,'S1',cex=0.9)
	
	matplot(0.5*dx,y1+dy*s1,pch=20,add=T,col=clrs[1])
	arrows(dx*0.5, y1+dy*S['S1','10%'], dx*0.5, y1+dy*S['S1','90%'], length=0.05,angle=90,code=3,lty=1)
	
	# K1
	arrows(S['K1','10%'], -20, S['K1','90%'], -20, length=0.05,angle=90,code=3,lty=1) # error bars
	text(k1,-50,'K1',cex=0.9)
	# ETo
	matplot(eto,-20,col=clrs[2],add=T,pch=5)
	text(eto,-50,'ETo',cex=0.9)
	arrows(S['ETo','10%'], -20, S['ETo','90%'], -20, length=0.05,angle=90,code=3,lty=1)
	
	# Legend
	leg_txt <- c('Obsv','Est','Model')
	legend('topleft',legend=leg_txt,col=c(clrs[1],clrs[1],clrs[2]),pch=c(19,1,18) )
	title(Title)
	grid()
	
	# Save
	dev.copy(png,plot_name,width=8,height=6,units="in",res=200)
	dev.off()
}

PlotRvP_1Intvl_AddSubIntvl <- function(D1,S1,D,S,Title='',avg_type='50%'){ 	# D1, S1 are from entire record, D, S from this sub-intvl
	clrs <- c('blue3','firebrick1','slategrey')	
	plotRvP <- function(D,clr,pch,addTF){
		x <- D$xobs; y <- D$yobs 				# x, y data
		# PLOT DATA
		matplot(x,y,pch=pch,col=clr,xlab='P (mm)',ylab='R (mm)',xlim=Xlim,ylim=Ylim,add=addTF) 	# obsv data
	}
	plotSegs <- function(S,clr){
		k1 <- S['K1',avg_type]; k2 <- S['K2',avg_type]  			# PARAMETER ESTIMATES
		s1 <- S['S1',avg_type]; eto<- S['ETo',avg_type]; s2 	 <- S['S2',avg_type]
		# Model Segments
		lx <- c(0,k1,k2,Xmax) 						
		ly <- modY(lx,k1,k2,s1,s2)
		segments(lx[1:3],ly[1:3],lx[2:4],ly[2:4],col=clr,lwd=1)
		matplot(lx,ly,pch=18,col=clr,add=T,cex=1.2)
	}
	
	# PLOT PARAMETER RANGES
	plotParams <- function(S,clr){
		k1 <- S['K1',avg_type]; k2 <- S['K2',avg_type]  			# PARAMETER ESTIMATES
		s1 <- S['S1',avg_type]; eto<- S['ETo',avg_type]; s2 	 <- S['S2',avg_type]
		# S1
		matplot(0.5*dx,y1+dy*s1,pch=20,add=T,col=clr)
		arrows(dx*0.5, y1+dy*S['S1','10%'], dx*0.5, y1+dy*S['S1','90%'], length=0.05,angle=90,code=3,lty=1,col=clr)
		# K1
		arrows(S['K1','10%'], -20, S['K1','90%'], -20, length=0.05,angle=90,code=3,lty=1,col=clr) # error bars
		text(k1,-50,'K1',cex=0.9,col=clr)
		# ETo
		matplot(eto,-20,col=clr,add=T,pch=5)
		text(eto,-50,'ETo',cex=0.9,col=clr)
		arrows(S['ETo','10%'], -20, S['ETo','90%'], -20, length=0.05,angle=90,code=3,lty=1,col=clr)
	}
	
	# plot limits
	Xmax <- 1.2*max(D1$xobs); 
	Xlim <- c(0,Xmax)
	Ymax <- 1.2*max(D1$yobs); Ymin <- -40
	Ylim <- c(Ymin,Ymax)
	
	plotRvP(D1,clrs[1],pch=1,addTF=F)
	plotSegs(S1,clrs[1])
	
	plotRvP(D,clrs[2],pch=19,addTF=T)
	plotSegs(S,clrs[2])
	
	# values for scaled S1 subplot
	dx <- diff(Xlim)/40; dy <- diff(Ylim)/4
	y1 <- diff(Ylim)/4;  y2 <- y1+dy
	# scaled subplot axis
	segments(0,y1,dx,y1)
	segments(0,y2,dx,y2)
	text(1.2*dx,y1,'0',cex=0.8)
	text(1.2*dx,y2,'1',cex=0.8)
	text(0.5*dx,y2+dy*0.1,'S1',cex=0.9)
	
	plotParams(S1,clrs[1])
	plotParams(S,clrs[2])
	
	# Legend
	leg_txt <- c('All','Intvl','Model')
	legend('topleft',legend=leg_txt,col=clrs,pch=c(1,19,18) )
	plot_name <- paste('PvR_1Intvl_AddSubIntvl.png',sep='')
	title(Title)
	grid()
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

PlotRvP_2Intvl <- function(D,S,Title='',avg_type='50%'){ 		# (data, params, changepoint index)
	iTend <- dim(D)[1]
	iTcp <- S['Tcp','mean']
	clrs <- c('blue3','indianred3')
	
	x <- D$xobs
	y <- D$yobs
	
	# MODEL ESTIMATES
	k1.1 <- S['K1',avg_type]; k1.2 <- k1.1+S['dK1',avg_type]
	k2.1 <- S['K2',avg_type]; k2.2 <- k2.1+S['dK2',avg_type]
	s1.1 <- S['S1',avg_type]; s1.2 <- s1.1+S['dS1',avg_type]
	eto.1<- S['ETo',avg_type];eto.2<- eto.1+S['dETo',avg_type]
	s2   <- S['S2',avg_type] 
	
	Xmax <- 1.2*max(x)
	Xlim <- c(0,Xmax)
	Ymax <- 1.2*max(y)
	Ymin <- -50
	Ylim <- c(Ymin,Ymax)
	
	# PLOT OBSV DATA
	x1 <- x[1:iTcp-1]
	y1 <- y[1:iTcp-1]
	matplot(x1,y1,pch=1,col=clrs[1],xlab='P (mm)',ylab='R (mm)',xlim=Xlim,ylim=Ylim)
	x2 <- x[iTcp:iTend]
	y2 <- y[iTcp:iTend]
	matplot(x2,y2,pch=1,col=clrs[2],add=T)
	
	# INTERVAL 1
	lx.1 <- c(0,k1.1,k2.1,Xmax) 				# x and y values for linear segments
	ly.1 <- modY(lx.1,k1.1,k2.1,s1.1,s2)
	segments(lx.1[1:3],ly.1[1:3],lx.1[2:4],ly.1[2:4],col=clrs[1],lwd=1)
	matplot(lx.1,ly.1,pch=18,col=clrs[1],add=T,cex=1.2)
	
	# INTERVAL 2
	lx.2 <- c(0,k1.2,k2.2,Xmax) 				# x and y values for linear segments
	ly.2 <- modY(lx.2,k1.2,k2.2,s1.2,s2)
	segments(lx.2[1:3],ly.2[1:3],lx.2[2:4],ly.2[2:4],col=clrs[2],lwd=1)
	matplot(lx.2,ly.2,pch=18,col=clrs[2],add=T,cex=1.2)
	
	# K1 UNCERTAINTY BARS	
	arrows(k1.1+S['dK1','10%'], -20, k1.1+S['dK1','90%'], -20,
				 length=0.05,angle=90,code=3,lty=1) # error bars
	
	# ETo + error bars
	matplot(eto.1,Ymin,col=clrs[1],add=T,pch=5)
	matplot(eto.2,Ymin,col=clrs[2],add=T,pch=5)
	arrows(eto.1+S['dETo','10%'], Ymin, eto.1+S['dETo','90%'], Ymin,
				 length=0.05,angle=90,code=3,lty=1)
	
	# S1
	dx <- diff(Xlim)/20
	dy <- diff(Ylim)/6
	y1 <- diff(Ylim)/3
	y2 <- y1+dy
	segments(0,y1,dx,y1)
	matplot(dx+0.1*dx,y1,add=T,cex=0.8,pch='0')
	segments(0,y2,dx,y2)
	matplot(0.5*dx,y2+dy*0.2,add=T,pch='S',cex=0.8)
	matplot(dx+0.1*dx,y2,add=T,cex=0.8,pch='1')
	matplot(0.2*dx,y1+dy*s1.1,pch=18,add=T,col=clrs[1])
	matplot(0.8*dx,y1+dy*s1.2,pch=18,add=T,col=clrs[2])
	arrows(dx*0.8, y1+dy*(s1.1+S['dS1','10%']), dx*0.8, y1+dy*(s1.1+S['dS1','90%']),
				 length=0.05,angle=90,code=3,lty=1)
	# Legend
	leg_txt <- c('pre','post','ETo')
	legend('topleft',legend=leg_txt,col=c(clrs,'blue'),pch=c(1,1,5) )
	grid()
	plot_name <- paste('PvR_2segs.png',sep='')
	title(Title)
	# 	dev.copy(png,plot_name,width=8,height=6,units="in",res=200)
	# 	dev.off()
}

PlotPD_2Intvl <- function(PD){
	par(mfrow=c(4,2))
	plot(density(PD$K1))
	plot(density(PD$K2))
	plot(density(PD$dK1))
	plot(density(PD$dK2))
	plot(density(PD$ETo))
	plot(density(PD$S1))
	plot(density(PD$dETo))
	plot(density(PD$dS1))
	par(mfrow=c(1,1))
}