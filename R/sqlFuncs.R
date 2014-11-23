# sqlFuncs

# WS INFO
SqlGetWsInfo <- function(sid) {
  require(RSQLite)
  path.sql <- file.path('/Users/tcmoran/Documents/CODE/R/WB_dB', 'gagesIIca.sqlite')
  db <- dbConnect(SQLite(), dbname=path.sql)
  sql.query <- paste('SELECT * FROM catchment_info WHERE site_num =', sid, ';')
  ws.info <- dbGetQuery(db, sql.query)
  dbDisconnect(db)
  return(ws.info)
}

SqlGetPall <- function(sid, Ptype = 'VICs') {
	# Args:
	# 	sid = USGS ID, must be one of 213 GagesII watersheds
	# 	Ptype = precip data type ('PRISM', 'VIC' , 'VICs', 'GHCN')
	# Returns:
	# 	[wy, P(mm)]
	require(RSQLite)
	path.sql <- file.path('/Users/tcmoran/Documents/CODE/R/WB_dB', 'gagesIIca213.sqlite')
	db <- dbConnect(SQLite(), dbname=path.sql)
	sql.query <- paste0('SELECT wy, Pmm FROM PrecipAll_', Ptype, ' WHERE STAID = ', sid, ';')
	wy.Pmm <- dbGetQuery(db, sql.query)
	dbDisconnect(db)
	if (nrow(wy.Pmm)==0) {
		cat('No precip data for', sid, '\n')
		wy.Pmm = NULL
	}
	return(wy.Pmm)
}

SqlGetWB <- function(sid, Ptype = 'VICs') {
	# Args:
	# 	sid = USGS ID, must be one of 213 GagesII watersheds
	# 	Ptype = precip data type ('PRISM', 'VIC' , 'VICs', 'GHCN')
	# Returns:
	# 	[wy, P(mm), R(mm)]
	require(RSQLite)
	path.sql <- file.path('/Users/tcmoran/Documents/CODE/R/WB_dB', 'gagesIIca213.sqlite')
	db <- dbConnect(SQLite(), dbname=path.sql)
	sql.query <- paste0('SELECT wy, Pmm, Rmm FROM WaterBalance_', Ptype, ' WHERE STAID = ', sid, ';')
	wy.P.R.mm <- dbGetQuery(db, sql.query)
	dbDisconnect(db)
	if (nrow(wy.P.R.mm)==0) {
		cat('No water balance data for', sid, '\n')
		wy.P.R.mm = NULL
	}
	return(wy.P.R.mm)
}

SqlGetPetMean <- function(sid) {
	# Args:
	# 	sid = USGS ID, must be one of 213 GagesII watersheds
	# Returns:
	# 	Mean PET as calculated using CIMIS geospatial data [PETmean(MM)]
	require(RSQLite)
	path.sql <- file.path('/Users/tcmoran/Documents/CODE/R/WB_dB', 'gagesIIca213.sqlite')
	db <- dbConnect(SQLite(), dbname=path.sql)
	sql.query <- paste0('SELECT MeanPetCimis FROM PETmean_CIMIS WHERE STAID = ', sid, ';')
	pet.mm <- dbGetQuery(db, sql.query)
	dbDisconnect(db)
	if (nrow(pet.mm)==0) {
		cat('No CIMIS PET data for', sid, '\n')
		pet.mm = NULL
	}
	return(pet.mm)
}

# SEGMENTED LINEAR MODEL PARAMETERS AND UNCERTAINTIES
SqlGetParams <- function(sid) {
  require(RSQLite)
  # Get Model Parameters from SQL
  path.sql <- file.path('/Users/tcmoran/Documents/CODE/R/WB_dB', 'gagesIIca213.sqlite')
  db <- dbConnect(SQLite(), dbname=path.sql)
  sql.query <- paste('SELECT * FROM ModelParamFit WHERE STAID =', sid, ';')
  model.fit <- dbGetQuery(db, sql.query)
  dbDisconnect(db)
  if (nrow(model.fit)==0) {
    cat('No model params for', sid)
    return(NULL)
  }
  rownames(model.fit) <- model.fit[['Param']] # make param types rownames
  return(model.fit)
}

