# SegLinFit_RunSuite2
source('R/SEGLINFIT_PARAMS.R') 				# MODEL PARAMETERS
# source('R/SegLinFit_SYNTHDATA.R')			# SYNTHETIC DATA
source('R/sqlFuncs.R') 								# GET WS DATA FROM SQL DB
# source('R/SegLinFit_CHANGEPOINTS.R') 	# CHANGEPOINT FUNCTIONS
source('R/MixNorm_Estimate2.R') 			# MIXED GAUSS EST FOR SOURCE DATA
source('R/SegLinFit_Stan.R') 					# RSTAN INITIALIZE & RUN
source('R/wbAnalysisPlots1.R') 				# RvP Plots, including uncert spread
# source('R/SegLinFit_PLOTS_RvP.R') 		# RvP Plots, including multi-interval
source('R/wbAnalysisHelpers.R')				# MISC HELPERS
source('R/SegLinFit_RUN_STAN_1intvl.R')# RUN STAN FUNCTIONS: 1 INTVL
# source('R/SegLinFit_MovingEstim5.R') 	# RUN STAN: MOVING INTVL
dirOut <- file.path(dir_master,'_SegLinFit_DataOut')

# Stan setup params
niter 			<- 2000
intvlType  	<- '1'
source('R/SegLinFit_RunStanSetup.R') 	# STAN SETUP FCNS

WStype 	<- 'all_20'
ixWS 		<- getWSix(WStype)
ixWS 		<- getIX(key='Ukiah')
saveSubDir <- NULL
M 			<- length(ixWS)
OUT.DATA<- list()
for (ii in 1:M) {
	Out.1intvl<- GoStanGo(segfit=segfit, pars=out.StanInit$pars, ixWS=ixWS[ii], 
												niter=niter, saveSubDir=saveSubDir)
	# Out.Moving<-MovingFit(Out.1intvl,pars=out.StanInit$pars,niter=niter,m=m[ii])
	if (ii>1) 
		graphics.off()
	out.data <- SaveOutput(ixWS[ii], Out.1intvl, dirOut)
	OUT.DATA[paste('id', out.data$staid, sep='')] <- 
		list(list(ypr=out.data$ypr, PDquants=out.data$PDquants, 
							ModelVals.xy=out.data$ModelVals.xy)
				 )
}
# save.OUT <- file.path(dirOut, paste0('Data_', WStype,  '_',Sys.Date(), '.Rdata'))
# save('OUT.DATA', file=save.OUT)
