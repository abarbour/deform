segall92 <- function(help=FALSE){
  cat("\nThis function is simply a placeholder. See the documentation ( `?segall92` ).\n")
  if (help) ?segall92
}

#' @rdname segall92
#' @export
.surface_J_ringdisk <- function(x, ring.depth, radius, thickness){
	
	z_prime <- ring.depth
	r <- x
	r_prime <- radius
	
}

	# Integral in the Green's function that must be integrated
	.green_integ <- function(k., r., a., z_prime., radial. = TRUE){
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
	
surface_displacement_ringdisk <- function(x, radius, depth, thickness, p0=1, nu=1/4, ShearModulus=15e9, BiotCoef=1, x.lim, verbose=TRUE, tol_pow){

	if (depth < thickness/2) stop("must be from: <depth - thickness/2>  to  <depth + thickness/2>")
	
	xpos <- x
	which_neg <- xpos < 0
	xpos[xpos < 0] <- abs(xpos[xpos < 0])


	x.lim <- if (missing(x.lim)){
		ceiling(2*max(xpos, na.rm=TRUE))
	} else {
		as.numeric(x.lim)
	}
	
	# Geetsma's uniaxial poroelastic coefficient
	# Wang (2000) 3.71
	c_m <- BiotCoef / ShearModulus
	
	scaling <- 2 * c_m * (1 - nu) * p0 * thickness * radius

	# Integration parameters
	nsub <- 101L
	tol <- .Machine$double.eps^ifelse(missing(tol_pow), 0.35, tol_pow)
	
	.fun <- function(r.., ...){
		int <- integrate(.green_integ, lower=0, upper=x.lim, stop.on.error=FALSE, r. = r.., ..., subdivisions = nsub, rel.tol=tol)
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
