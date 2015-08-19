#' Surface deformation associated with fluid withdrawl: an alternative formulation
#' @name vasco98
#' @aliases vasco1998
#' @export
#' @param help logical; load documentation for \code{\link{vasco98}}
#' @references Vasco, D., et al. (1998), Monitoring of Fluid Injection and 
#' Soil Consolidation Using Surface Tilt Measurements,
#' \url{http://ascelibrary.org/doi/abs/10.1061/(ASCE)1090-0241(1998)124:1(29)}
#' @seealso \code{\link{segall85}}, or \code{\link{segall89}}; \code{\link{Simple-deformation}}
#' @export
#' @examples
#' # Reproduce...
#' r <- seq(-10,10,by=0.2) # km
#' r.m <- r*1e3
#' xx <- surface_displacement_point(r.m, depth=2e3, delV.=1e7)
#' plot(uz ~ x, xx, type='l')
vasco98 <- function(){
  cat("\nThis function is simply a placeholder. See the documentation ( `?vasco98` ).\n")
  if (help) ?vasco98
}

#' @rdname vasco98
#' @export
.surface_g_pt <- function(x=0, x_src=0, z_src=0, nuu=1/3){
  # vasco98 11
  sc <- (1 + nuu)/(3*pi)
  xrel <- (x - x_src)
  g <- sc * xrel / (xrel^2 + z_src^2)**(3/2)
  data.frame(x, g, xz = x/z_src)
}

#' @rdname vasco98
#' @export
#' @inheritParams surface_displacement_reservoir
#' @param depth numeric; the depth below the surface to the source
#' @param delV. numeric; uniform change in pore-fluid volume (positive = loss)
#' @param tol numeric; the numerical tolerance in the integration; if 
#' supplied, this should be much smaller than the smallest difference between
#' any values in \code{x} divided by the depth
#' @param ... additional arguments passed to \code{\link{.surface_g_pt}}
surface_displacement_point <- function(x, depth, delV., B.=1, rho_f.=1000, tol, ...){
  x <- as.vector(x)
  depth <- depth[1]
  if (depth <= 0) stop('source must be below the surface (positive depth)')
  delV. <- delV.[1]
  # Setup some vectors/info for interpolation
  if (missing(tol)) tol <- 0.01
  xr <- range(pretty(x))
  mx <- ceiling(max(abs(x), na.rm = TRUE)/depth)
  mxsc <- 15 # this many times the x/z distance
  xcalc <- depth * seq.int(from = -1* mxsc * mx, to = mxsc * mx, by = tol)
  # Calculate point source response
  sgcalc <- deform::.surface_g_pt(xcalc, z_src = depth, ...)
  # Calculate response function
  ix <- sgcalc$x
  iy <- sgcalc$g
  n <- length(ix)
  stopifnot(length(n) < 2)
  ih <- (ix[2] - ix[1])
  ca <- (iy[2] - iy[1]) / ih
  cb <- (iy[n] - iy[n - 1]) / ih
  Gfun <- stats::approxfun(ix, pracma::cumtrapz(ix, iy) - (cb - ca) * ih^2/12)
  # Use interpolation function to resample integrated response function
  g <- Gfun(x)
  C. <- B. * delV. / rho_f.
  sg <- data.frame(x, g, xz = x/depth, uz = g * C.)
  return(sg)
}
