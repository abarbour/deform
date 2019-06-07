#' @title Elastostatic subsidence of a half space
#' 
segall92 <- function(){
  cat("\nThis function is simply a placeholder. See the documentation ( `?segall92` ).\n")
}

# Integral in the Green's function that must be integrated
# @rdname segall92
# @export
.green_segall92 <- function(k., r., a., z_prime., radial. = TRUE){
  # r radial observation point
  # a radius of disk
  # z_prime depth disk
  # k vector along direction of integration (k >= 0)
  ka <- k. * a.
  kr <- k. * r.
  kz <- k. * z_prime.
  J0 <- Bessel::BesselJ(kr, 0)
  J1 <- Bessel::BesselJ(kr, 1)
  J1a <- Bessel::BesselJ(ka, 1)
  E <- exp(-1 * kz)
  if (radial.){
    # radial displacements at the surface [ur(r,0)]
    # Wang (2000) equation 9.29
    J1 * J1a * E
  } else {
    # vertical displacements at the surface [uz(r,0)]
    # Wang (2000) equation 9.30
    J0 * J1 * E
  }
}

#' @rdname segall92
#' @param x 
#' @param radius,depth,thickness numeric; the radius, midpoint depth, and thickness of the circular-disk pressure source; the
#' midpoint-depth must be at least \code{thickness/2}
#' @param p0 numeric; the (uniform) pressure change inside the circular-disk source
#' @param nuu numeric; the undrained Poisson's ratio
#' @param ShearModulus numeric; the elastic shear modulus
#' @param BiotCoef numeric; Biot's coefficient (usually written as \eqn{\alpha})
#' @param verbose logical; should messages be given?
#' @param x.lim numeric; the limit of integration in the radial direction
#' @param tol_pow numeric; the exponent applied to \code{.Machine$double.eps} used as the integration tolerance (see \code{\link{integrate}})
#' @param nsub integer; the number of subdivisions in the integration (see \code{\link{integrate}})
#' 
#' @references Segall, P. (1992), Induced stresses due to fluid extraction from axisymmetric reservoirs,
#' PAGEOPH, 139: 535. doi: 10.1007/BF00879950
#' @export
#' 
#' @examples 
#' 
#' # Calculate surface displacements from a pressure-disk source
#' 
#' # Radial distances (e.g., observation points) and pressure source strength
#' # (in this case, keep distances in km to stabilize integration
#' #  and then rescale pressure to maintain unit consistency)
#' r <- seq(0, 15, length.out=151)
#' Pa <- 80e3 #Pa (N / m^2)
#' pressure <- 2.5 * Pa * (1000^2)
#' 
#' # dimensions of inflating (pressure > 0) or deflating (pressure < 0) reservoir
#' disk.radius <- 4
#' disk.thickness <- 2
#' disk.middepth <- 3
#' 
#' # calculations...
#' res <- surface_displacement_ringdisk(r, radius = disk.radius, thickness = disk.thickness, depth = disk.middepth, p0 = pressure)
#' head(res)
#' 
#' with(res,{
#'   uz_sc <- -uz
#'   ur_sc <- ur
#'   yl <- range(c(uz_sc, ur_sc, err), na.rm=TRUE)
#'   plot(x, uz_sc, col=NA, log='', ylim=yl, xaxs='i',
#'        ylab="Displacement, mm", xlab="Distance from source, km", main="Surface Displacements: Pressure-disk Source")
#'   abline(h=0, col='grey')
#'
#'   lines(x, uz_sc, lwd=1.5)
#'   text(10, 7, "Vertical")
#'
#'   lines(x, ur_sc, lty=2, lwd=1.5)
#'   text(5, 5, "Radial")
#'
#'   lines(x, err, lty=3, lwd=1.5)
#'   text(4.5, -0.75, "Radial strain\n(mm/km)", pos=2)
#'
#'   mtext(sprintf("Segall (1992) model with radius: %s, Thickness: %s, Depth: %s", disk.radius, disk.thickness, disk.middepth), font=3)
#' })
#' 
surface_displacement_ringdisk <- function(x, radius, depth, thickness, 
                                          p0=1, nuu=1/4, ShearModulus=15e9, BiotCoef=1, verbose=TRUE, 
                                          x.lim, tol_pow, nsub){

	if (depth < thickness/2) stop("must be from: <depth - thickness/2>  to  <depth + thickness/2>")
	
	xpos <- x
	which_neg <- xpos < 0
	xpos[xpos < 0] <- abs(xpos[xpos < 0])


	x.lim <- if (missing(x.lim)){
		ceiling(2*max(xpos, na.rm=TRUE))
	} else {
		as.numeric(x.lim)
	}
	
	# Geetsma's uniaxial poroelastic coefficient (Young's modulus?)
	# Wang (2000) 3.71
	c_m <- BiotCoef / ShearModulus
	
	scaling <- 2 * c_m * (1 - nuu) * p0 * thickness * radius

	# Integration parameters
	nsub <- ifelse(missing(nsub), 101L, as.integer(nsub))
	tol <- .Machine$double.eps^ifelse(missing(tol_pow), 0.35, tol_pow)
	
	.fun <- function(r.., ...){
		int <- stats::integrate(.green_segall92, lower=0, upper=x.lim, stop.on.error=FALSE, r. = r.., ..., subdivisions = nsub, rel.tol=tol)
		int[['value']]
	}
	   
   	if (verbose) message("Calculating displacements... r")
	ur_ <- sapply(xpos, .fun, radial.=TRUE, a. = radius, z_prime. = depth)
	ur <- scaling * ur_
	
	if (verbose) message("Calculating displacements... z")
	uz_ <- sapply(xpos, .fun, radial.=FALSE, a. = radius, z_prime. = depth)
	uz <- -1 * scaling * uz_
	
	if (verbose) message("Calculating strain...")
	err <- deform::Uniaxial_extension(x, ur)$dXdx
	
	data.frame(x, scaling, ur, uz, err)
	
} 
