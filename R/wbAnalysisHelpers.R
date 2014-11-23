
source('~/Documents/CODE/R/RSTAN_SegLinFit/SegLinFit_IMPORT_DATA.R')		# IMPORT OBSV DATA
dir_master <- "~/Desktop/Catchment Analysis 2011/AA_CA_Catchments_Master/GAGESII_CATCHMENTS_CA/GAGESII_CATCHMENTS_219"

if (!exists('StaDirs')){ 			# station directory listing
	StaDirs <- ImportStaDirs()
}
getIX <- function(key){
	grep(paste0(key),as.list(StaDirs$V1),ignore.case=T)
}
getSval <- function(S,val){
	ix <- which(regexpr(val,rownames(S))>0)
	if (length(ix)>1) ix <- which(regexpr(paste0(val,'\\['),rownames(S))>0)
	S[ix,]
}
getStaInfo <- function(data){
	n <- nchar(data$dirData)
	ddir <- substr(data$dirData,1,n-21)
	sta.info <- read.csv(file.path(ddir,'site_info.csv'))
}
getStaInfo2<- function(ddir){
	sta.info <- read.csv(file.path(ddir,'site_info.csv'),stringsAsFactors=F)
}

loadModelData <- function(data.dir,sid){
	load(file.path(data.dir,'RSTAN_1IntvlFitSummary.Rdata'))
	load(file.path(dir_master,'_SegLinFit_DataOut',paste0('Data_',sid,'.Rdata')))
	return(list(S=S,mod.vals=ModelVals.xy,PDquants=PDquants))
}

ixType <- function(Type='All'){
	# Alp with >= 30% precip as snow
	ixNonAlp <- c(130,	66,	4	,68,	98,	96,	75,	129,	2,	69,	158,	25,	156,	91,	64,	83,	56,	20,	101,	18,	133,	45,	76,	85,
								191,	73,	90,	159,	57,	86,	71,	84,	42,	115,	169,	160,	108,	28,	21,	1,	19,	74,	184,	12,	114,	185,	
								100,	140,	99,	89,	44,	217,	109,	50,	17,	61,	149,	128,	171,	190,	52,	62,	33,	29,	88,	181,	173,	179,
								82,	163,	214,	134,	46,	32,	182,	154,	155,	141,	55,	11,	105,	150,	119,	103,	183,	94,	197,	117,	146,
								144,	54,	174,	7,	48,	79,	194,	87,	127,	97,	157, 112,	118,	113,	208,	151,	124,	107,	145,	147,	215,
								51,	36,	65,	193,	95,	165,	38,	192,	152,	189,	120,	196)
	ixAlp <- c(188,	166,	142,	27,	137,	35,	47,	72,	164,	195,	205,	15,	59,	136,	202,	148,	187,	78,	177,	198,	110,
						121,	167,	178,	26,	204,	199,	200,	211,	186,	216,	207,	201,	123,	106,	122,	210,	63,	40,	153,	143,	
						139,	138)
	if (Type=='NonAlp')	ixWS <- ixNonAlp
	if (Type=='Alp') ixWS <- ixAlp
	if (Type=='All') ixWS <- c(ixNonAlp,ixAlp)
	return(ixWS)
}

sidType <- function(Type='All'){
	sidNonAlp <- c(11132500,	11120500,	11075720,	11086500,	11134800,	11046300, 11154100,	11120550,	
		11160300,	11117600,	11141500,	11110500,	11337500,	11467500,	11162570,	11084500,
		11180500,	11147070,	11153900,	11463900,	11182100,	11180960,	11464500,	11096500,
		11309000,	11147000,	11271320,	11197250, 11451720,	11180825,	11147500,	11449500,
		11468000,	11152900,	11274500,	11458000,	11458500,	11407300,	11456000,	11149900,
		11139000,	11454000, 11390672,	11143500,	11169800,	11253310,	11160020,	11224500,
		11453500,	11299600,	11274630,	11172100,	11452000,	11260480,	11468500,	11116000,
		11481200,	11451100,	11461000,	11378800,	11451500,	11407500,	11220500,	11114500,
		10255810,	11033000,	11098000,	11176400,	11269300,	11173200,	11176000,	11259000,
		11124500,	10258500,	10259200,	11258900,	11306000,	11472200,	11211300,	11475500,
		11476500,	11199500,	11138500,	11475800,	10257600,	11470500,	11379500,	11334300,
		11284400,	10258000,	11336000,	11374000,	11015000,	11482500,	11220000,	11080500,
		11094000,	11111500,	11382000,	11473900,	11384000,	11474500,	11521500,	10264000,
		11257100,	11379000,	11409300,	11281000, 11204500,	11525600,	11433260,	11481500,
		10263500,	11316800,	11342000,	11480390,	11282000,	11381500,	11528700,	11522500,
		11118000, 11465200)
	sidAlp <- c(11283500,	11383500,	11284700,	11341400, 11210100,	11404500,	11394500,	11401200,
		11278000,	11489500,	10360900,	11251000,	10356500,	11209900,	10336645,	11414000,
		11266500,	11237500,	10336676,	10292000,	11296500,	10289000,	11426150,	10343500,
		11294000,	10336780,	10308200,	10308201,	10310000,	11264500,	11292500,	11276500,
		10336600,	10296500,	11213500,	10296000,	10291500,	11436000,	10287400, 10290500,	
		11292000,	11230500,	11214000)
	if (Type=='NonAlp')	sids <- sidNonAlp
	if (Type=='Alp') sids <- sidAlp
	if (Type=='All') sids <- c(sidNonAlp,sidAlp)
	return(sids)
}

outDir <- function(Type){
	if (Type=='NonAlp')	output.dir <- file.path(dir_master,'OUTPUT','NonAlp_20yrs')
	if (Type=='Alp')	output.dir <- file.path(dir_master,'OUTPUT','Alp_20yrs')
	if (Type=='All')	output.dir <- file.path(dir_master,'OUTPUT','All_20yrs')
	return(output.dir)
}

getSID <- function(data){
	idx <- regexpr('_Site',data$dirData)
	as.numeric(substr(data$dirData,idx[1]+5,idx[1]+12))
}
getSID2 <- function(ixWS) getSID(ImportData(ixWS,plotTF=F))

getSegFitData <- function(sid,pdTF=F){
	dir.data <- file.path(dir_master,StaDirs[getIX(sid),],'PARAM_FIT','PDATA_VICs')
	load(file.path(dir.data,'RSTAN_1IntvlFitSummary.Rdata'))
	if (!pdTF) return(list(pd.summ=S))
	load(file.path(dir.data,'RSTAN_1IntvlFitData.Rdata'))
	return(list(pd.summ=S,pd.data=StanOut1))
}



