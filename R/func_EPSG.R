##
agency.root <- function(agency="epsg"){
  agency <- match.arg(agency)
  toret <- sprintf("http://spatialreference.org/ref/%s",agency)
  return(toret)
}
##
getEPSG <- function(ref, ref.type=c("proj"), verbose=TRUE){
  #http://spatialreference.org/ref/epsg/2805/proj4/
  #http://spatialreference.org/ref/epsg/4326/proj4/
  ref.type=switch(match.arg(ref.type), proj="proj4")
  ref.url <- sprintf("%s/%i/%s/", agency.root("epsg"), as.numeric(ref), ref.type)
  if (verbose) message(ref.url)
  require(RCurl)
  stopifnot(RCurl::url.exists(ref.url))
  prj <- RCurl::getURL(ref.url)
  return(prj)
}
##
latlon_project <- function(latlon, to.epsg="4326"){
  latlon <- as.data.frame(latlon)
  coordinates(latlon) <- names(latlon)[1:2]
  # proj4 string
  eprj <- getEPSG(to.epsg, ref.type="proj", verbose=F)
  stopifnot(exists("eprj"))
  ##
  proj4string(latlon) <- eprj
  toret <- spTransform(latlon, CRS=CRS(eprj))
  return(toret)
}
##