# SegLinFit_ReAnalysis
# A GENERAL FRAMEWORK FOR RE-ANALYZING WATER BALANCE PARAMETER ESTIMATES
# MODIFY AS NEEDED

# Helper functions
source('~/Documents/CODE/R/WB_shared/wsList.R')


# specify data to load
dir.data.load <- '/Users/tcmoran/Desktop/Catchment Analysis 2011/AA_CA_Catchments_Master/GAGESII_CATCHMENTS_CA/GAGESII_CATCHMENTS_219/OUTPUT_NOV2014/data'
data.load <- '_RSTAN_1IntvlFitData.Rdata'
# where to save new/updated data/plots
dir.data.save <- '/Users/tcmoran/Desktop/Catchment Analysis 2011/AA_CA_Catchments_Master/GAGESII_CATCHMENTS_CA/GAGESII_CATCHMENTS_219/OUTPUT_NOV2014/plots'
data.save <- '_PD_density.png'

# which SIDs to process
# wsSids <- WsSidList(Type='All')
wsSids <- WsSidModList()

for (sid in wsSids) {
	# load saved data
	fname.load <- paste0(sid, data.load)
	path.load  <- file.path(dir.data.load, fname.load)
	load(path.load)
	
	# process
	
	# save new/updated data/plots
	fname.save <- paste0(sid, data.save)
	path.save  <- file.path(dir.data.save, fname.save)
	PlotPD_1Intvl(PD, path.save)	
	dev.off()
	
}