# Functions from https://github.com/vlambert/ReservoirStresses
# converted with matconv::mat2r, and then hand-edited
#
# AJB NOTES: 
# * Verified from Segall 1989 that epsneg/epspos are indeed (x+-resHalfWidth)/resTopDepth, NOT x+-resHalfWidth/resTopDepth as written in LT2020 appendix
# * Added erfc function -- (real)
# ?  What is y1?  horizontal interval −resHalfWidth < y1 < resHalfWidth, as well as with depth resTopDepth < y2 < resTopDepth +resThickness, reflecting resHalfWidth producing layer of fixed length 2a and thickness resThickness

# y1, y2 are horizontal and vertical positions of center of dilatation source

erfc_ <- function(x){
	2 * pnorm(x * sqrt(2), lower.tail = FALSE)
}


# Solutions in Lambert & Tsai (2020)
LambertTsai2020_Diffusive <- function(x1, x2, y1, resTopDepth, resThickness, resDiffusivity, Time){
	
	#                     Free Surface
	#   <----------------------------------------------------> x1
	#                           |  ^
	#                           |  |          mu = shear modulus
	#                           |  resTopDepth
	#                           |  |
	#                           |  v
	#               <-----------|----------->  dm(Time,resDiffusivity)     Reservoir
	#                           |                        (thickness resThickness)
	#                           |
	#                           v
	#                           x2
	
	sfct <- sqrt(4 * resDiffusivity * Time)
	
	# Change in fluid mass distribution
	dm <- sqrt(Time / resDiffusivity) * (exp(-y1^2 / sfct^2) / sqrt(pi) - abs(y1) / sfct * erfc_(abs(y1) / sfct))
	dm[Time == 0] <- 0
	
	DT <- resTopDepth + resThickness
	dxy1 <- x1 - y1

	#Horizontal extensional surface strain, eps11 
	f_e11 <- (resTopDepth*DT - dxy1^2) / (resTopDepth^2 + dxy1^2) / (DT^2 + dxy1^2)
	eps11 <- 2 * resThickness * f_e11 * dm
	
	# Horizontal normal stress, Sigma_11 / mu 
	f1_s11 <- (resTopDepth - x2) / ((resTopDepth - x2)^2 + dxy1^2) + 
		(3 * resTopDepth + x2) / ((resTopDepth + x2)^2 + dxy1^2) + 
		(4 * x2 * dxy1^2) / ((resTopDepth + x2)^2 + dxy1^2)^2
	f2_s11 <- (DT - x2) / ((DT - x2)^2 + dxy1^2) + 
		(3 * DT + x2) / ((DT + x2)^2 + dxy1^2) + 
		(4 * x2 * dxy1^2) / ((DT + x2)^2 + dxy1^2)^2
	sig11 <- dm * (f1_s11 - f2_s11)
	
	# Shear stress, Sigma_12 / mu 
	f1_s12 <- 2 * resTopDepth^3 + x2 * resTopDepth^2 + x2 * (x2^2 + dxy1^2) + 2 * resTopDepth * dxy1^2
	f1p_s12 <- f1_s12 / (((resTopDepth - x2)^2 + dxy1^2) * ((resTopDepth + x2)^2 + dxy1^2)^2)
	f2_s12 <- 2 * DT^3 + x2 * DT^2 + x2 * (x2^2 + dxy1^2) + 2 * DT * dxy1^2
	f2p_s12 <- f2_s12 / (((DT - x2)^2 + dxy1^2) * ((DT + x2)^2 + dxy1^2)^2)
	sig12 <- 4 * dm * x2 * dxy1 * (f1p_s12 - f2p_s12)
	
	# Vertical normal stress, Sigma_22 / mu 
	f1_s22 <- resTopDepth * (x2^2 + 3 * dxy1^2) + x2 * (x2^2 + dxy1^2) - x2 * resTopDepth^2 - resTopDepth^3
	f1p_s22 <- f1_s22 / (((x2 - resTopDepth)^2 + dxy1^2) * ((x2 + resTopDepth)^2 + dxy1^2)^2)
	f2_s22 <- DT * (x2^2 + 3 * dxy1^2) + x2 * (x2^2 + dxy1^2) - x2 * DT^2 - DT^3
	f2p_s22 <- f2_s22 / (((x2 - resTopDepth - resThickness)^2 + dxy1^2) * ((x2 + DT)^2 + dxy1^2)^2)
	sig22 <- 4 * x2^2 * dm * (f1p_s22 - f2p_s22)
	
	# Horizontal surface displacements, u1
	f_u1 <- atan(DT / dxy1) - atan(resTopDepth / dxy1)
	u1 <- 2 * dm * f_u1
	
	#_u2 Vertical surface displacements, u2  
	f_u2  <- log(resTopDepth^2 + dxy1^2) - log(DT^2 + dxy1^2)
	u2 <- dm * f_u2
	
	Sig <- cbind(S11=sig11, S12=sig12, S22=sig22)
	
	ret <- list(
		diffusion.kernels = TRUE,
		pars = data.frame(x1=x1, x2=x2, y1=y1, Time=Time),
		res = data.frame(top.depth=resTopDepth, thickness=resThickness, half.width=NA, diffusivity=resDiffusivity),
		dM = dm,
		Def = cbind(U1=u1, U2=u2, Sig, E11=eps11)
	)
	class(ret) <- c('LT.diffusive', 'list')
	ret
}


LambertTsai2020_Uniform <- function(x1, x2, resTopDepth, resThickness, resHalfWidth){
	x1.o <- x1
	x1 <- abs(x1) # weird results if x1<0 are computed, but this is symmetric model anyway ??
	
	stopifnot(all(c(x2, resThickness, resHalfWidth) >= 0)) # depth reckoned positive
	stopifnot(all(c(resTopDepth) > 0)) # topdepth may not be zero
	
	#see also deform::segall89
	
	#                     Free Surface
	#   <----------------------------------------------------> x1
	#                           |  ^
	#                           |  |          mu = shear modulus
	#                           |  resTopDepth
	#                           |  |
	#                           |  v
	#             -resHalfWidth -----------|----------- resHalfWidth      Reservoir
	#                           |                 (thickness resThickness)
	#                           |
	#                           v
	#                           x2
	
	x2D <- x2 + resTopDepth
	x2Dm <- x2 - resTopDepth
	
	x1a <- x1 + resHalfWidth
	x1am <- x1 - resHalfWidth
	
	# Horizontal normal stress, Sigma_11 / mu 
	f1_s11 <- 4 * resHalfWidth * x2 * ( 
		-(resHalfWidth^2 - x1^2 + x2D^2) / ((x1am^2 + x2D^2) * (x1a^2 + x2D^2)) +
    	(resHalfWidth^2 - x1^2 + (x2D + resThickness)^2) / ((x1am^2 + (x2D + resThickness)^2) * (x1a^2 + (x2D + resThickness)^2))
    )
	f2_s11 <- atan((x2Dm - resThickness) / x1am) -
		atan((x2Dm - resThickness) / x1a) + 
		atan(x2Dm/ x1a) - 
		atan(x2Dm/ x1am)
	f3_s11 <- 3 * (
		atan(x2D/ x1am) - 
		atan(x2D/ x1a) - 
		atan((x2D + resThickness)/ x1am) + 
		atan((x2D + resThickness)/ x1a)
	)
	sig11 <- f1_s11 + f2_s11 + f3_s11
	
	
	# Shear stress, Sigma_12 / mu 
	f1_s12 <- -log(x1am^2 + x2Dm^2) + log(x1am^2 + (x2Dm - resThickness)^2) 
	#f1 = -log((x1-a).^2+(x2-D).^2) + log((x1-a).^2+(x2-D-T).^2);
	f2_s12 <- (16 * resHalfWidth * x1 * x2 * x2D + (x1am^2 + x2D^2) * (x1a^2 + x2D^2) * (log(x1a^2 + x2Dm^2) + log(x1am^2 + x2D^2) - log(x1a^2 + x2D^2))) / ((x1am^2 + x2D^2) * (x1a^2 + x2D^2))
	f3_s12 <- -(16 * resHalfWidth * x1 * x2 * (x2D + resThickness) + (x1am^2 + (x2D + resThickness)^2) * (x1a^2 + (x2D + resThickness)^2) * (log(x1a^2 + (x2Dm + resThickness)^2) + log(x1am^2 + (x2D + resThickness)^2) - log(x1a^2 + (x2D + resThickness)^2))) / ((x1am^2 + (x2D + resThickness)^2) * (x1a^2 + (x2D + resThickness)^2))
	sig12 <- (f1_s12 + f2_s12 + f3_s12) / 2
	
	
	# Vertical normal stress, Sigma_22 / mu [seems okay]
	f1_s22 <- 4 * resHalfWidth * x2 * (
		(resHalfWidth^2 - x1^2 + x2D^2) / ((x1am^2 + x2D^2) * (x1a^2 + x2D^2)) - 
		(resHalfWidth^2 - x1^2 + (x2D + resThickness)^2) / ((x1am^2 + (x2D + resThickness)^2) * (x1a^2 + (x2D + resThickness)^2))
	)
	f2_s22 <- atan((x2Dm - resThickness) / x1a) - 
		atan((x2Dm - resThickness) / x1am) - 
		atan(x2Dm / x1a) + 
		atan(x2Dm / x1am)    
	f3_s22 <- atan(x2D/ x1am) - 
		atan(x2D/ x1a) - 
		atan((x2D + resThickness)/ x1am) + 
		atan((x2D + resThickness)/ x1a)
	sig22 <- f1_s22 + f2_s22 + f3_s22
	
	
	epspos <- (x1.o + resHalfWidth) / resTopDepth
	epsneg <- (x1.o - resHalfWidth) / resTopDepth
	
	# Horizontal surface displacements, u1 
	u1 <- log((1 + epspos^2) / (1 + epsneg^2)) / 2

	# Vertical surface displacements, u2 
	u2 <- atan(epsneg) - atan(epspos)

	#  Horizontal extensional surface strain, epsilon_11 
	eps11 <- ((epspos / (1 + epspos^2) - epsneg / (1 + epsneg^2))) / resTopDepth

	Sig <- cbind(S11=sig11, S12=sig12, S22=sig22)
	#Tau <- Shear(Sig[,'S11'], Sig[,'S12'], Sig[,'S22'])
	
	ret <- list(
		diffusion.kernels = FALSE,
		pars = data.frame(x1=x1.o, x2=x2, y1=NA, Time=NA),
		res = data.frame(top.depth=resTopDepth, thickness=resThickness, half.width=resHalfWidth, diffusivity=NA),
		dM = NA,
		Def = cbind(U1=u1, U2=u2, Sig, E11=eps11)
	)
	class(ret) <- c('LT.uniform', 'list')
	ret
}

#--------------------------------------------------------------
# Plane-stress computations

shear <- function(x, ...) UseMethod('shear')
shear.LT.diffusive <- shear.LT.uniform <- function(x, angle=FALSE){
	Def <- x[['Def']]
	Def[,'S12']
}
ext1 <- function(x, ...) UseMethod('ext1')
ext1.LT.diffusive <- ext1.LT.uniform <- function(x, angle=FALSE){
	Def <- x[['Def']]
	Def[,'S11']
}
ext2 <- function(x, ...) UseMethod('ext2')
ext2.LT.diffusive <- ext2.LT.uniform <- function(x, angle=FALSE){
	Def <- x[['Def']]
	Def[,'S22']
}

max_shear_principal <- function(x, ...) UseMethod('max_shear_principal')
max_shear_principal.LT.diffusive <- max_shear_principal.LT.uniform <- function(x){
	max_shear_principal.default(S1=ext1(x), S3=ext2(x))
}
max_shear_principal.default <- function(S1, S3) (S1 - S3)/2

max_shear <- function(x, ...) UseMethod('max_shear')
max_shear.LT.diffusive <- max_shear.LT.uniform <- function(x, angle=FALSE){
	p <- PSPS(x)
	p[,ifelse(angle,'tau','tau_theta')]
}

mean_stress <- function(x, ...) UseMethod('mean_stress')
mean_stress.LT.diffusive <- mean_stress.LT.uniform <- function(x){
	Def <- x[['Def']]
	mean_stress.default(S1=ext1(x), S3=ext2(x))
}
mean_stress.default <- function(S1, S3) (S1 + S3)/2

PSPS <- function(x, ...) UseMethod('PSPS')

PSPS.LT.diffusive <- PSPS.LT.uniform <- function(x, ...){
	PSPS.default(s11=ext1(x), s12=shear(x), s22=ext2(x))
}

PSPS.default <- function(s11, s12, s22){ # plane-stress principal stresses (i.e., 2D Mohr's circle)
	ms <- mean_stress.default(s11, s22)
	maxs <- max_shear_principal.default(s11, s22)

	shearstress <- sqrt(maxs^2 + s12^2)
	shearangle <- atan(-maxs / s12) / 2

	S1 <- ms + sqrt(maxs^2 + s12^2)
	S2 <- ms - sqrt(maxs^2 + s12^2)

	cbind(S1=S1, S2=S2, tau = shearstress, tau_theta = shearangle)
}

#--------------------------------------------------------------
# Fault-plane projection computations

fault_normal_stress <- function(x, theta=0, ...) UseMethod('fault_normal_stress')

fault_normal_stress.LT.diffusive <- fault_normal_stress.LT.uniform <- function(x, theta=0, ...){
	fault_normal_stress.default(s11=ext1(x), s12=shear(x), s22=ext2(x), theta = theta, ...)
}

fault_normal_stress.default <- function(s11, s12, s22, theta, is.deg=TRUE){
	# normal stress acting on a plane oriented by theta radians from the 1-2 coord system
	ms <- mean_stress.default(s11, s22)
	maxs <- max_shear_principal.default(s11, s22)
	if (is.deg) theta <- theta * pi / 180
	ms + maxs * cos(2*theta) + s12 * sin(2*theta)
}

fault_shear_stress <- function(x, theta=0, ...) UseMethod('fault_shear_stress')

fault_shear_stress.LT.diffusive <- fault_shear_stress.LT.uniform <- function(x, theta=0, ...){
	fault_shear_stress.default(s11=ext1(x), s12=shear(x), s22=ext2(x), theta = theta, ...)
}

fault_shear_stress.default <- function(s11, s12, s22, theta, is.deg=TRUE){
	# shear stress acting on a plane oriented by theta radians from the 1-2 coord system
	maxs <- max_shear_principal.default(s11, s22)
	if (is.deg) theta <- theta * pi / 180
	# use double-angle identity: sin(2x) = 2*sin(x)cos(x)
	-maxs*sin(2*theta) + s12*(cos(theta)^2 - sin(theta)^2)
}

# (cfs_isoporo.defaultin beeler)
cfs_isoporo.LT.diffusive <- cfs_isoporo.LT.uniform <- function(x, theta=0, fric=0.65, Skemp=0.6, verbose=TRUE, ...){
	if (verbose) message('CFS {S1-to-fault angle = ', theta, '°, friction = ', fric, "} for isotropic undrained poroelastic response {B = ", Skemp, "}")
	tauf <- fault_shear_stress(x, theta)
	Snorm <- fault_normal_stress(x, theta)
	Smean <- mean_stress(x)
	cfsip <- cfs_isoporo.default(tau=tauf, mu=fric, Snormal = Snorm, Smean = Smean, B=Skemp)
	zapsmall(cfsip)
}

