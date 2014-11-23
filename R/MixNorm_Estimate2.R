# ESTIMATE MIXED GAUSSIAN FIT TO INPUT DATA
# estimates the components of mixed Gaussians fit to input data

library('mixtools')

MixNorm_Estimate <- function(K, d, plotTF=T) {
	
	# ESTIMATES ARE RANDOM, SO CHOOSE MAX LIKELIHOOD ESTIMATE OF 10 ITERATIONS
	mixmdl <- normalmixEM(d, k=K, arbvar=F, fast=F) 				# initial estimate
	for(ii in 1:10) {
		mixmdl.ii <- normalmixEM(d, k=K, arbvar=F, fast=F); 	# FIXED VARIANCE SOLVES PROB OF OVERLY-PEAKED MIXTURES
		# print(mixmdl$loglik) 		# display log likelihood
		if(mixmdl.ii$loglik > mixmdl$loglik) 	# replace model if log-lik is greater 
			mixmdl <- mixmdl.ii
	}
	
	MixTheta <- mixmdl$lambda 	# mixture fraction
	MixMu    <- mixmdl$mu 			# distribution means
	MixSD    <- mixmdl$sigma 		# distribution SDs
	
	if (plotTF) {
		# PLOT DATA HISTOGRAM, DENSITY PLOT, INDIVIDUAL GAUSSIANS, AND SUM OF MIXED GAUSSIANS
		plot(mixmdl, which=2, main2='Data Distrib & Gauss Fits')
		lines(density(d), lty=2)
		mx <- seq(min(d) - 0.2*abs(min(d)), max(d) + 0.2*abs(max(d)))
		my <- matrix(0, length(mx))
		for(k in 1:K)
			my <- my + MixTheta[k]*dnorm(x=mx,MixMu[k],MixSD[k])

		lines(mx, my, lty=1, lwd=2, col='blue')
		legend('topright',c('Sum ','Gauss1 ','Gauss2 ','Density ','Histogram'),
					 col=c('blue','red','green','black','black'), 
					 lty=c(1,1,1,2,1),lwd=c(2,1,1,1,1))
		dev.copy(png, 'PRECIP_DISTRIB_2MIXED_GAUSS.png', width=8, height=6,
						 units="in",res=200)
		dev.off()
	}
	
	XdataMix <- as.data.frame(cbind(MixTheta, MixMu, MixSD))
}
# XdataMix <- MixNorm_Estimate(2, DataSynth$xobs)