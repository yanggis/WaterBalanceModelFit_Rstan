# SegLinFit_ReAnalysis
# A GENERAL FRAMEWORK FOR RE-ANALYZING WATER BALANCE PARAMETER ESTIMATES
# MODIFY AS NEEDED

# Helper functions
source('~/Documents/CODE/R/WB_shared/wsList.R')
source('~/Documents/CODE/R/WaterBalanceModelFit_Rstan/R/wbAnalysisPlots1.R')
source('~/Documents/CODE/R/WB_shared/sqlGetFuncs2.R')
source('~/Documents/CODE/R/WB_shared/filterPRdata.R')

# specify data to load
dir.data.load <- '/Users/tcmoran/Desktop/Catchment Analysis 2011/AA_CA_Catchments_Master/GAGESII_CATCHMENTS_CA/GAGESII_CATCHMENTS_219/OUTPUT_JAN2015/data'
data.load <- '_PEXC_VICs_TVICs_k1_0_RSTAN_1IntvlFitData.Rdata'
# where to save new/updated data/plots
dir.data.save <- '/Users/tcmoran/Desktop/Catchment Analysis 2011/AA_CA_Catchments_Master/GAGESII_CATCHMENTS_CA/GAGESII_CATCHMENTS_219/OUTPUT_JAN2015/plots'
data.save <- '_PEXC_VICs_TVICs_k1_0_PD_density.png'

# which SIDs to process
wsSids <- WsSidList(Type='All')

# misc settings
p.type = 'PEXC'
p.source='VICs_TVICs'

for (sid in wsSids) {
	# load saved data
	fname.load1 <- paste0(sid, data.load)
	path.load1  <- file.path(dir.data.load, fname.load1)
	load(path.load1)
	fname.load2 <- paste0(sid, '_PEXC_VICs_TVICs_k1_0_RSTAN_1IntvlFitSummary.RData')
	path.load2  <- file.path(dir.data.load, fname.load2)
	load(path.load2)
	P.all <- SqlGetPall2(sid, p.type=p.type, p.source=p.source)[ , -1] # exclude WY, first col
	pet.mean <- SqlGetPetMean(sid=sid)$MeanPetCimis
	PRdata <- SqlGetWB2(sid=sid, p.type=p.type, p.source=p.source)
	PRdata <- filterPRref(sid=sid, wyPR=PRdata, p.source='VICs')
	
	# process
	
	# save new/updated data/plots
	# 	fname.save <- paste0(sid, data.save)
	# 	path.save  <- file.path(dir.data.save, fname.save)
	# RvP
	pRvP <- plotRvPspread(sid=sid, pd=PD, Summ=pdSumm, P.all=P.all, pet=pet.mean, PR=PRdata, scatTF=F)
	# pRvP <- plotRvPspreadETi(sid=sid, pd=PD, Summ=pdSumm, P.all=P.all, pet=pet.mean, PR=PRdata, scatTF=F)
	path.rvp <- file.path(dir.data.save, paste0(sid, '_R_v_PEXC_VICs_TVICs_k1_0_spread.png'))
	ggsave(path.rvp, plot = pRvP, width = 11, height = 8, units = 'in', dpi=72)
}