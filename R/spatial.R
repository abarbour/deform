#' Transform cartesian coordinates to polar
#' @param x numeric vector, or numeric matrix (see \code{y})
#' @param y numeric vector; if missing, then \code{x} must have at least
#' two columns: the first to represent \code{x} and the second to represent \code{y}
#' @param z numeric vector
#' @param dist.type character; the type of distance calculation
#' @param radians logical; should the output angle be in radians? (\code{FALSE} returns in degrees)
#' @author modified from function in beadarrayMSV
#' @export
#' @examples
#' 
#' cart2pol(1:10, 1:10)
#' all.equal(cart2pol(1:10, 1:10), cart2pol(cbind(1:10, 1:10)))
#' 
#' # Change radius computation
#' cart2pol(1:10, 1:10, 1:10) # assumes 'euclidean'
#' cart2pol(1:10, 1:10, 1:10, 'manhattan')
#' 
#' # use degrees
#' cart2pol(1:10, 1:10, radians=FALSE)
#' 
cart2pol <- function(x, y, z = NULL, dist.type = c('euclidean','manhattan'), radians=TRUE){
  pol <- list()
  dist.type <- match.arg(dist.type)
  if (missing(y)){
    X <- as.matrix(x)
    stopifnot(ncol(X) >= 2)
    y <- X[,2]
    x <- X[,1]
  }
  x <- as.vector(x)
  y <- as.vector(y)
  r. <- switch(dist.type,
         euclidean = sqrt(x^2 + y^2),
         manhattan = x + y
  )
  th. <- atan2(y, x)
  if (!radians) th. <- th. * 180 / pi
  rth <- cbind(Theta = th., Radius = r.)
  if (!is.null(z)){
    rth <- cbind(rth, Z = as.vector(z))
  }
  return(rth)
}

#' Functions to work with EPSG projection specifications
#' @name epsg
#' @aliases EPSG
#' @param agency character;
#' @param ref.id integer;
#' @param verbose logical;
#' @param Df an object to be coerced into a \code{\link{data.frame}}
#' @param coords character; the names of the coordinates in \code{Df}
#' @param p4,p4old; these are the proj4 strings for the new (\code{p4})
#'  and old (\code{p4old}) data.  More specifically, \code{p4} will be
#'  the new projection (set by \code{\link{CRS}}), and  
#'  \code{p4old} is the old projection (set by \code{\link{proj4string}}).
#' @param ... additional parameters to \code{\link{elide}}
#' @references 
#' A list of EPSG projections:
#'   \url{http://spatialreference.org/ref/epsg/}
#' @references 
#' A list of PROJ commands:
#'   \url{https://trac.osgeo.org/proj/wiki/GenParms}
#' @examples
#' \dontrun{
#' library(maps)
#' library(maptools)
#' library(datasets)
#' data(state)
#' states <- data.frame(state.x77, state.center)
#' # HI/AK are in the middle of the ocean, so take them out:
#' states <- states[states$x > -121,]
#' 
#' # uses getEPSG to project if p4/p4old are missing
#' d1 <- D_project(states, c("x","y"))
#' 
#' # set p4
#' p <- getEPSG()
#' d2 <- D_project(states, c("x","y"), p)
#' all.equal(d1,d2) # TRUE
#' 
#' # reproject into UTM-11 in kilometers
#' d2u <- spTransform(d2, CRS("+proj=utm +zone=11 +ellps=WGS84 +units=km"))
#' 
#' # change the projection
#' pold <- p
#' p <- getEPSG(2153)  # UTM-11 in meters
#' d3 <- D_project(states, c("x","y"), p, pold)
#' 
#' layout(matrix(1:6,ncol=3))
#' us <- map("usa",plot=FALSE)
#' plot(d1, axes=TRUE); title("d1"); lines(us)
#' plot(d2, axes=TRUE);title("d2 (== d1)"); lines(us)
#' plot(d2u, axes=TRUE);title("d2 - utm11/km")
#' plot(d3, axes=TRUE);title("d3 - utm11/m")
#' #
#' # rotate with maptools (can also do this using ...)
#' plot(elide(d2u, rotate=45));title("d2 at 45") # clockwise
#' plot(elide(d2u, rotate=-45));title("d2 at -45") # counter-clockwise
#' #
#' }
NULL

#' @rdname epsg
#' @export
agency.root <- function(agency=NULL){
  agency <- match.arg(agency,c("epsg","esri","iau2000","sr-org"))
  return(sprintf("http://spatialreference.org/ref/%s",agency))
}

#' @rdname epsg
#' @export
getEPSG <- function(ref.id=4326, verbose=TRUE, agency=NULL){
  ref.id <- as.integer(ref.id)
  #http://spatialreference.org/ref/epsg/2805/proj4/
  #http://spatialreference.org/ref/epsg/4326/proj4/
  ref.url <- sprintf("%s/%s/proj4/", agency.root(agency), ref.id)
  if (verbose) message(ref.url)
  stopifnot(url.exists(ref.url))
  prj <- getURL(ref.url)
  return(prj)
}

#' @rdname epsg
#' @export
D_project <- function(Df, coords=names(Df)[1:2], p4, p4old, verbose=TRUE, ...){
  Df <- data.frame(Df)
  stopifnot(length(coords)>=2)
  #
  # proj4 string
  if (missing(p4)) p4 <- getEPSG(verbose=verbose)
  if (missing(p4old)) p4old <- p4
  eprj <- CRS(p4)
  stopifnot(exists("eprj"))
  ##
  coordinates(Df) <- coords
  proj4string(Df) <- p4old
  #
  toret <- try(elide(spTransform(Df, CRSobj=eprj), ...))
  stopifnot(!inherits(toret,"try-error"))
  return(toret)
}
