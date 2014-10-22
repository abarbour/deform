#' Compute the length of 1-degree of latitude and longitude
#' @export
#' @param lat numeric; degrees
#' @param model (optional) character; which ellipsoidal correction should be used?
#' @author A.J.Barbour, adapted from \code{degreelen.m} by B.W.Crowell
#' @references \url{http://en.wikipedia.org/wiki/World_Geodetic_System}
#' @examples
#' lats <- 30:32
#' # Lon should be about 111.4 km, and Lat 96 - 94 km
#' degreelen(lats)
#' 
#' # In most cases there is no discerible effect
#' # caused by varying the semi-major axis.
#' #   Not much difference according to this machine:
#' all.equal(degreelen(lats,"IRTF"), degreelen(lats,"WGS84"), tolerance = .Machine$double.eps)
#' 
#' #   unless the precision of the input is in micrometers, 
#' #   which is essentially impossible to achieve:
#' options(digits=13)
#' lats <- 30.1111111123
#' rbind(degreelen(lats,"IRTF"),degreelen(lats,"WGS84"))
degreelen <- function(lat, model){
  #
  pic <- pi/180
  latrad <- as.numeric(lat) * pic
  #
  if (missing(model)) model <- "IRTF"
  Ell <- ellipsoid(model)
  Ea. <- Ell[["major"]] # major axis
  Fe. <- Ell[["flat"]]  # flattening
  EE. <- Ell[["ecc"]]   # eccentricity
  #
  #http://surveying.wb.psu.edu/sur162/control/surfaces.htm
  slp <- sin(latrad)
  slpsq <- slp * slp
  #
  degDist <- 1 - EE. * slp * slp
  degDistsq <- sqrt(degDist)
  #
  N. <- Ea. / degDistsq
  M. <- Ea. * cos(latrad) * (1 - EE.) / degDistsq / degDist ## equiv to ^1.5
  #
  return(cbind(Lat=lat, (cbind(Meters.lon=N., Meters.lat=M.) * pic)))
}

#' @rdname degreelen
#' @export
ellipsoid <- function(model=c("IRTF","WGS84")){
  # ITRF ellipsoid (GRS80) is default
  mdl <- match.arg(model)
  message(mdl)
  Ea <- 6378137.0
  Eb <- 6356752.0 + switch(mdl, IRTF=0.314140, WGS84=0.314245)
  Flat <- 1 - Eb/Ea       # flattening
  Ecc <- Flat*(2 - Flat)  # eccentricity
  return(data.frame(model=mdl, major=Ea, minor=Eb, flat=Flat, ecc=Ecc))
}
