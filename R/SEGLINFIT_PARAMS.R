# SEGLINFIT MODEL PARAMS

# MODEL PARAMETER UNCERTAINTY BOUNDS
ParamUncertSettings <- function(PET) {
	# Uncertainty parameters for RStan priors
	sdX 		<- c(10, 0, 0.1) 			# [minSD,minCoef,maxCoef]
	sdY 		<- c(0.1, 0, 0.1)			# [minSD,minCoef,maxCoef] - minSD = 0 for yobs
	sdI			<- c(25, 500) 				# [min, max] intrinsic uncert stddev
	S2_lims <- c(0.9, 1.1) 				# (minS2, maxS2)
	PETlim 	<- 1.0*PET 						# PET estimate from CIMIS ETo product
	ParamUncert <- list(sdX=sdX, sdY=sdY, sdI=sdI, S2_lims=S2_lims, PETlim=PETlim)
}
