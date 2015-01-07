# SegLinFit_RunSuite3
# source('R/SEGLINFIT_PARAMS.R') 				# MODEL PARAMETERS
# source('R/SegLinFit_SYNTHDATA.R')			# SYNTHETIC DATA
# source('R/SegLinFit_CHANGEPOINTS.R') 	# CHANGEPOINT FUNCTIONS
# source('R/SegLinFit_Stan2.R') 					# RSTAN INITIALIZE & RUN

# source('R/SegLinFit_PLOTS_RvP.R') 		# RvP Plots, including multi-interval
# source('R/wbAnalysisHelpers.R')				# MISC HELPERS
# source('R/SegLinFit_RUN_STAN_1intvl.R')# RUN STAN FUNCTIONS: 1 INTVL
# source('R/SegLinFit_MovingEstim5.R') 	# RUN STAN: MOVING INTVL

source('R/MixNorm_Estimate2.R') 			# MIXED GAUSS EST FOR SOURCE DATA
source('R/StanSetup.R')
source('R/StanRun.R')
source('~/Documents/CODE/R/WB_shared/wsList.R')
source('~/Documents/CODE/R/WB_shared/sqlGetFuncs2.R')
source('~/Documents/CODE/R/WB_shared/SqlWriteFcns.R')
source('~/Documents/CODE/R/WB_shared/filterPRdata.R')
source('R/wbAnalysisPlots1.R') 				# RvP Plots, including uncert spread

dir.master = '/Users/tcmoran/Desktop/Catchment Analysis 2011/AA_CA_Catchments_Master/GAGESII_CATCHMENTS_CA/GAGESII_CATCHMENTS_219'
dir.save <- file.path(dir.master, 'OUTPUT_JAN2015')

# STAN PARAMS
niter 			<- 2000
chains 			<- 2
intvlType  	<- '1'
etiTF <- F
# INITIALIZE STAN
out.StanInit <- InitStan(intvlType)
segfit <- out.StanInit$segfit
fit.pars <- out.StanInit$pars

# WATERSHED IDs
SIDs <- WsList(Type='Custom')
# SIDs <- 10343500
# sid <- 11378800
M <- length(SIDs)

plotTF = T
saveTF = T
i1 = 1
# DATA PARAMS
p.type = 'PEXC'
p.source='VICs_TVICs'
file.note = '_k1_0'
for (sid in SIDs[i1:M]) {
	PRdata <- SqlGetWB2(sid=sid, p.type=p.type, p.source=p.source)
	PRdata <- filterPRref(sid=sid, wyPR=PRdata, p.source='VICs')
	if (etiTF) {
		eti <- SqlGetETi(sid=sid, p.source=p.source, wys=PRdata$wy)
		PRdata <- merge(PRdata, eti, by = 'wy')
	}
	P.all <- SqlGetPall2(sid, p.type=p.type, p.source=p.source)[ , -1] # exclude WY, first col
	PdataGMix <- MixNorm_Estimate(K=2, d=P.all, plotTF=F)
	pet.mean <- SqlGetPetMean(sid=sid)$MeanPetCimis
	ParamUncert <- ParamUncertSettings(pet.mean)
	
	StanOut <- StanRun(segfit=segfit, PRdata=PRdata, ParamUncert=ParamUncert, fit.pars=fit.pars,
										 PdataGMix=PdataGMix, niter=niter, chains=chains, parTF=T)
	PD 	<- extract(StanOut$SegFit, permuted=T)
	Dstan <- StanOut$Dstan
	pdSumm <- summary(StanOut$SegFit, probs=c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975))$summary
	
	if (plotTF) {
		cat('Saving Plots \n')
		# Posterior Densities
		path.dens <- file.path(dir.save, 'plots', paste0(sid, '_', p.type, '_', p.source, file.note, '_PD_density.png'))
		PlotPD_1Intvl2(PD, path.dens, saveTF=T)

		# RvP
		pRvP <- plotRvPspread(sid=sid, pd=PD, Summ=pdSumm, P.all=P.all, pet=pet.mean, PR=PRdata, scatTF=F)
		# pRvP <- plotRvPspreadETi(sid=sid, pd=PD, Summ=pdSumm, P.all=P.all, pet=pet.mean, PR=PRdata, scatTF=F)
		path.rvp <- file.path(dir.save, 'plots', paste0(sid,'_R_v_', p.type, '_', p.source, file.note, '_spread.png'))
		ggsave(path.rvp, plot = pRvP, width = 11, height = 8, units = 'in', dpi=72)
		dev.off()
	}
	
	# SAVE DATA
	if (saveTF){
		cat('Saving PD data \n')
		# SUMMARY TABLE
		path.summ <- file.path(dir.save, 'data', paste0(sid, '_', p.type, '_', p.source, file.note, '_RSTAN_1IntvlFitSummary.RData'))
		save('pdSumm', file=path.summ)
		
		# MODEL, POSTERIORS, PARAMETERS
		path.pd <- file.path(dir.save, 'data', paste0(sid, '_', p.type, '_', p.source, file.note, '_RSTAN_1IntvlFitData.RData'))
		save(list=c('PD', 'Dstan', 'ParamUncert', 'PdataGMix'), file=path.pd)
	}
}
