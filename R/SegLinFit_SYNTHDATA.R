# SYNTHETIC DATA FOR SEGMENTED LINEAR .STAN SCRIPT
# 3 segment constrained linear with intrinsic and y & x errors
#		 specify sig_x and sig_y as multipliers, not single value

# REQUIRED PACKAGES
require('truncnorm')

SynthDataCalc <- function(N,plotTF=T){
	# Produces synthetic data of similar form to observed data for annual water 
	# balances in California. Allows for parameter shifts following some change
	# point t_cp
	# Args:
	# 	N: number of years of data
	# Returns:
	# 	$DataObs: synthetic data observations [yrobs, xobs, yobs, Tx, Ty]
	# 	$Pall: 	  superset of 'P' observations, actually same as xobs for this data
	
	# TRUE MODEL PARAMS
	NxAdd <- 50 					# additional number of x data for distrib calc but not obsv
	yrobs <- 1950:(1950+N-1)
	# basic model params in terms of the params directly estimated by RStan
	mod 			<- list()
	mod$K1 		<- 500					# starting value of K1
	mod$DK 		<- 1500					# starting value of K2-K1
	mod$ETo   <- 1000					# starting value of ETo
	mod$S2 		<- 1 						# slope of 3rd segment
	
	mod$t_cp 	<- round(N*1/2) # time index for change in parameters - first year of new behavior
	mod$dK1 	<- -300					# change in param values for t > t_cp
	mod$dDK 	<- 0
	mod$dETo 	<- 200 					
	
	# derived params - done this way for consistency with RStan model in/out
	mod$K2  <- with(mod, K1 + DK)
	mod$S1  <- with(mod, S2*(K2-ETo)/DK)
	mod$dK2 <- with(mod, dDK + dK1)
	mod$dS1 <- with(mod, S2*((K2+dK2)-(ETo+dETo))/(DK+dDK) - S1)
	
	# UNCERTAINTY PARAMS
	mod$sd_I 			<- 30				# std dev of intrinsic model uncertainty
	sig_x_min <- 5				# min sd of x obsv
	sig_x_coef<- 0.05 		# x obsv multiplicative uncert coef
	sig_y_min <- 0 				# min sd of y obsv
	sig_y_coef<- 0.05 		# y obsv multiplicative uncert coef
	
	# HEAVISIDE FUNCTION
	H <- function(x) {
		# Heaviside (step) function: 0 for x < 0, 1/2 for x==0, 1 for x > 0
		return((sign(x) + 1)/2)
	}
	
	# 3-SEGMENT LINEAR MODEL
	modY <- function(x, k1, k2, s1, s2) {	
		# Tri-linear constrained water balance model
		# Returns modeled Y values given x values and params [k1, k2, s1, s2]
		modY <- s1*(x-k1)*H(x-k1) + (s2-s1)*(x-k2)*H(x-k2)	
	}
	
	# True X values - CALCULATE MORE THAN OBSERVED, AS WITH REAL OBSV WHERE NP > NR
	TxCalc <- function(N, K1, K2) {
		# Models x values as drawn from mixed distribution of two Gaussians
		# Args:
		# 	N: 	number of samples
		# 	K1:	true value of model param K1, used to scale distribution means
		# 	K2: true value of model param K1, used to scale distribution means
		# Returns: 
		# 	N samples of x drawn from mixed Gaussian model
		
		meanK 	<- mean(c(K1, K2)) 		
		mu_Tx  	<- meanK*c(1, 2) 			# means of the two mixed Gaussians x is drawn from
		sig_Tx 	<- meanK*c(0.4, 0.5) 	# sigma of mixed Gauss
		theta_Tx<- c(1/2, 1/2) 				# mixture of Gaussians, sum = 1
		nTx 		<- floor(theta_Tx[1]*N)# fraction of occurrences in MG1 
		nTx 		<- c(nTx, N-nTx)
		Tx 			<- array(c(rnorm(nTx[1], mu_Tx[1], sig_Tx[1]), 
										 rnorm(nTx[2], mu_Tx[2], sig_Tx[2]))) # mixture of 2 Gaussians
		Tx  		<- sample(Tx) 				# random permutation of Tx
		Tx[Tx < 0]<- -Tx[Tx < 0] 			# wrap data around 0, a little wonky, but works
		return(Tx)
	}
	
	# True Y Values: piecewise linear model
	TyCalc <- function(Tx, mod) {
		# Calculates values of Y given X and adds Gaussian noise
		# Args
		# 	Tx:	x values
		# 	mod: model parameter list
		# Returns
		# 	Y values with intrinsic noise added
		
		Ty <- vector()
		for (t in 1:length(Tx)){ 	# calculate Ty for each Tx given changepoint t_cp
			dt <- t - mod$t_cp + 0.5
			s1 <- mod$S1 + mod$dS1 * H(dt)
			k1 <- mod$K1 + mod$dK1 * H(dt)
			k2 <- mod$K2 + mod$dK2 * H(dt)
			Ty[t] <- modY(Tx[t], k1, k2, s1, mod$S2)
		}
		# add intrinsic noise
		Ty <- rtruncnorm(N, a=0, mean=Ty, sd=mod$sd_I) 
	}
	
	# generic fcn to add observation noise
	ObsvNoisy <- function(Tz, sdmin, az, lowerBound=-Inf) {
		# Add multiplicative noise to inputs
		# Args
		# 	Tz: input data
		# 	sdmin: minimum standard deviation
		# 	az: uncertainty coefficient
		# 	lowerBound: minimum allowed observaton value (if any)
		# Returns
		# 	input data with multiplicative noise added
		
		sdz  <- sdmin + az*abs(Tz) 				# multiplicative error
		oz   <- rnorm(length(Tz),0,sdz)  	# random draws
		zobs <- Tz + oz
		zobs[zobs<lowerBound] <- lowerBound
		return(zobs)
	}
	
	# X TRUE, OBSERVED
	Nx 		<- N + NxAdd
	TOT_Tx<- TxCalc(Nx, mod$K1, mod$K2) 	
	TOT_xobs <- ObsvNoisy(TOT_Tx, sig_x_min, sig_x_coef)
	Tx 		<- TOT_Tx[1:N]
	xobs 	<- TOT_xobs[1:N]
	
	# X, Y OBSERVED
	Ty 		<- TyCalc(Tx, mod)
	yobs 	<- ObsvNoisy(Ty, sig_y_min, sig_y_coef, lowerBound=0)
	
	if (plotTF) {
		# PLOT
		matplot(xobs, yobs, pch=16, col='indianred3')
		matplot(Tx, Ty, pch=3, col='gray48', add=T)
		
		# PLOT MODEL SEGMENTS
		k1 <- with(mod, c(K1, K1+dK1))
		k2 <- with(mod, c(K2, K2+dK2))
		s1 <- with(mod, c(S1, S1+dS1))
		s2 <- with(mod, c(S2, S2))
		for (jj in 1:length(k1)) {
			lx <- c(0, k1[jj], k2[jj], max(xobs)) 	# x and y values for linear segments
			ly <- modY(lx, k1[jj], k2[jj], s1[jj], s2)
			segments(lx[1:3], ly[1:3], lx[2:4], ly[2:4], col='gray48', lwd=1)
			matplot(lx, ly, pch=18, col='royalblue1', add=T, cex=1.2)
		}
	}
	DataSynth <- as.data.frame(cbind(yrobs, xobs, yobs, Tx, Ty))
	return(list(DataObs=DataSynth, Pall=TOT_xobs))
}

# DataSynth <- SynthDataCalc(80) 			# FOR TESTING ONLY!!! 
