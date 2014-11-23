# importData.R
# Functions for importing various data
library(raster)


LoadBdryCA <- function() {  # California boundary
  # Load CA border file from shapefile, return as SpatialPolygonsDataFrame
  dir.ca <- file.path('data','CA_boundary')
  if (!file.exists(dir.ca))
    stop('Need valid directory for CA boundary shapefile, importData:LoadBdryCA()')
  bdry.CA <- readOGR(dir.ca,'CA_Bdry_Simpler')
}

BuildRasterCA <- function(path.Pbrick) {
  # Makes a raster brick of monthly PRISM data from a directory of .bil files
  stop('Build the PRISM raster with CA_Drought')
}

GetRasterCA <- function() {
  # Load PRISM raster brick if it has been saved, otherwise build it
  fname.Pbrick<- 'P_brick_ca.grd'
  path.Pbrick <- file.path('/Users/tcmoran/Documents/CODE/R/CA_drought/data/PRISM',
                           fname.Pbrick)
  if (file.exists(path.Pbrick)) {
    brick(path.Pbrick)
  } else {
    BuildRasterCA(path.Pbrick)
  }
}

LoadBoundary <- function(sid, polyTF=F) {  
  # Load WS border file from shapefile, return as SpatialPolygonsDataFrame
  dir.bdry <- '/Users/tcmoran/Documents/ENV_DATA/Boundaries/gagesII_CA'
  bdry.yx <- read.csv(file.path(dir.bdry, paste0(sid, '.txt')))
  bdry.xy <- subset(bdry.yx, select=c(2,1))
  colnames(bdry.xy) <- c('x', 'y')
  bdry.sp <- Polygon(bdry.xy)
  if (polyTF)
    return(bdry.sp)
  bdry.sp <- Polygons(list(bdry.sp), ID='boundary')
  bdry.sp <- SpatialPolygons(list(bdry.sp), proj4string=CRS("+init=epsg:4326"))
  bdry.sp <- spTransform(bdry.sp, CRS('+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0'))
}