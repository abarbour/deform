#' Error functions
#' 
#' @note These remain for posterity; we now use RcppFaddeeva's implementation (\code{\link[RcppFaddeeva]{erf}}, etc.).
#' 
#' @param x numeric
#' @export
erf_re <- function(x){
  2 * pnorm(x * sqrt(2)) - 1
}
#' @rdname erf_re
#' @export
erfc_re <- function(x){
  2 * pnorm(x * sqrt(2), lower.tail = FALSE)
}
#' @rdname erf_re
#' @export
ierf_re <- function (x){
  qnorm((1 + x) / 2) / sqrt(2)
}
#' @rdname erf_re
#' @export
ierfc_re <- function (x){
  qnorm(x/2, lower.tail = FALSE) / sqrt(2)
}
#' @rdname erf_re
#' @export
ierfc2_re <- function(x){
  exp(-1 * x^2) / sqrt(pi) - x * erfc(x)
}

#' Check whether a quantity is wihin an acceptable range
#' 
#' @note This function throws and error if range of \code{x} is outside [\code{xmin},\code{xmax}].
#' 
#' @author Andrew J. Barbour <andy.barbour@@gmail.com> 
#' @param x numeric, or coercible to numeric
#' @param xmin numeric; the minimum value of the acceptable range
#' @param xmax numeric; the maximum value of the acceptable range
#' 
#' @export
#' 
#' @return see \code{\link{stopifnot}}
check_range <- function(x, xmin=0, xmax=1){
  x <- na.omit(as.numeric(x))
  stopifnot(all(x >= xmin) & all(x <= xmax))
}

#' Simple numerical deformation estimates
#' 
#' Calculate tilts and extensions based on spatially varying displacements
#' 
#' @details
#' \code{\link{Uniaxial_extension}} calculates the component of
#' strain associated with deformation along a single axis.
#' For example, the change in the Eastward displacements (or rates)
#' in the East direction.
#' 
#' \code{\link{Tilt}} calculates the tilt field associates with
#' spatial variations in vertical positions (or rates of change);
#' the sign convention used is such that
#' a ball placed on the x-y plane would roll in the direction of the
#' tilt vector.  Or, in other words, the tilt vector is the direction a
#' plumb bob would move.
#' 
#' Calculations are done with \code{\link{diff}}.
#' 
#' @note \code{\link{Tilt}} has not been well-tested for two-dimensional
#' results!
#' @name Simple-deformation
#' @param x numeric; the spatial vector in real units (not an index!)
#' @param X numeric; the quantity varying in the direction of \code{x}
#' @param y numeric; the spatial vector perpendicular to \code{x} in real units
#' @param z numeric; values which vary in the direction perpendicular to the \code{x}-\code{y} plane
#' @param left logical; should padding be on the "left" side of the vector (i.e., index number 1)?
#' @export
#' @seealso \code{\link{segall85}}
#' @examples
#' # Some x values
#' xval. <- sort(unique(c(-7:-3, seq(-3.90,3.90,by=0.05), 3:7)))
#' set.seed(1221)
#' xanis <- rnorm(length(xval.), sd = 10)
#' 
#' #  surface displacements
#' su <- surface_displacement(xval.*1e3, C.=1e13, z_src=0.7e3)
#' 
#' # Vertical tilt -- assumes axial symmetry
#' sut <- with(su, Tilt(x, z=uz))
#' #               -- including anisotropic effects
#' sut.anis <- with(su, Tilt(x, x = sort(x + xanis), z=uz))
#' 
#' plot(ztilt ~ x, sut.anis, col='blue', pch=16, cex=0.5)
#' lines(ztilt ~ x, sut, lwd=2)
#' 
#' plot(ztilt ~ abs(x), sut.anis, col='blue', pch=ifelse(sign(x)==1,16,1), cex=0.5)
#' lines(ztilt ~ abs(x), sut, lwd=2)
#'  
#' # Uniaxial strain in the 'x' direction
#' sue <- with(su, Uniaxial_extension(x, X=ux))
#' sue.anis <- with(su, Uniaxial_extension(sort(x + xanis), X=ux))
#' 
#' plot(dXdx ~ x, sue.anis, col='blue', pch=16, cex=0.5)
#' lines(dXdx ~ x, sue, lwd=2)
#'  
#' plot(dXdx ~ abs(x), sue.anis, col='blue', pch=ifelse(sign(x)==1,16,1), cex=0.5)
#' lines(dXdx ~ abs(x), sue, lwd=2)
#' 
#' # see how the 'left' argument affects things:
#' .setleft(1,TRUE) # NA  1
#' .setleft(1,FALSE) # 1  NA
Uniaxial_extension <- function(x, X, left=TRUE){
  deriv <- diff(X)/diff(x)
  dXdx <- .setleft(deriv, left)
  data.frame(x, dXdx)
}
#' @rdname Simple-deformation
#' @export
Tilt <- function(x, y=NULL, z, left=TRUE){
  deriv1 <- diff(z)/diff(x)
  dzdx <- .setleft(deriv1, left)
  dzdy <- if (is.null(y)){
    dzdx
  } else {
    # this needs to be more intelligent:
    # what if z is a matrix!?
    deriv2 <- diff(z)/diff(x)
    .setleft(deriv2, left)
  }
  # reverse sign so that the sign convention
  # is such that a positive tilt corresponds to
  # the direction a plumb bob would head, or
  # the direction a ball on the surface would roll
  tlt <- -1*(dzdx + dzdy)
  ang <- atan2(dzdy, dzdx) * 180/pi 
  data.frame(x = x, ztilt = tlt, xy.direction = ang)
}
#' @rdname Simple-deformation
#' @export
.setleft <- function(x, left=TRUE){
  x <- as.vector(x)
  if (left){
    c(NA, x)
  } else {
    c(x, NA)
  }
}