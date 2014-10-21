#degreelen.m
#This takes in the latitude and computes the length of 1 degree of latitude
#and longitude
#lat is entered in degrees

# ITRF ellipsoid (GRS80), or WGS84 ??
degreelen <- function(lat, Ea=6378137.0, Eb=6356752.3142){
  #
  #http://surveying.wb.psu.edu/sur162/control/surfaces.htm
  Fe. <- 1 - Eb/Ea # flattening
  EE. <- 2*Fe. - Fe. * Fe. # eccentricity
  #
  pic = pi/180
  latrad = lat * pic
  slp <- sin(latrad)
  slpsq <- slp * slp
  #
  degDist <- 1 - EE. * slpsq
  degDistsq <- sqrt(degDist)
  #
  N. <- Ea / degDistsq
  M. <- Ea * cos(latrad) * (1 - EE.) / degDistsq / degDist ## equiv to ^1.5
  #
  toret <- zapsmall(cbind(N., M.) * pic)
  rownames(toret) <- lat
  return(toret)
}
##
##
