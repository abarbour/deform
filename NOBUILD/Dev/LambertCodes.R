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
	f2_s12 <- (
		16 * resHalfWidth * x1 * x2 * x2D + 
		(x1am^2 + x2D^2) * (x1a^2 + x2D^2) * (log(x1a^2 + x2Dm^2) + 
		log(x1am^2 + x2D^2) - 
		log(x1a^2 + x2D^2))
	) / (
		(x1am^2 + x2D^2) * (x1a^2 + x2D^2)
	)
	f3_s12 <- -(
		16 * resHalfWidth * x1 * x2 * (x2D + resThickness) + 
		(x1am^2 + (x2D + resThickness)^2) * (x1a^2 + (x2D + resThickness)^2) * (log(x1a^2 + (x2Dm + resThickness)^2) + 
		log(x1am^2 + (x2D + resThickness)^2) - 
		log(x1a^2 + (x2D + resThickness)^2))
	) / (
		(x1am^2 + (x2D + resThickness)^2) * (x1a^2 + (x2D + resThickness)^2)
	)
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

ang_max_shear <- function(x, ...) UseMethod('ang_max_shear')
ang_max_shear.LT.diffusive <- ang_max_shear.LT.uniform <- function(x){
	s11 <- ext1(x)
	s22 <- ext2(x)
	s12 <- shear(x)
	atan((s22 - s11)/(2 * s12)) / 2
}

planestress_principal_stresses <- function(x, ...) UseMethod('planestress_principal_stresses')
planestress_principal_stresses.LT.diffusive <- planestress_principal_stresses.LT.uniform <- function(x){
	s1 <- ext1(x)
	s2 <- ext2(x)
	s12 <- shear(x)
	zer <- 0*s1
	S <- cbind(s1, s12, zer, s12, s2, zer, zer, zer, zer)
	S[!is.finite(S)] <- NA
	.get_eig <- function(si){
		sim <- matrix(si, 3)
		navals <- rep(NA, 3)
		vals <- if (any(is.na(sim) | !is.finite(sim))){
			navals
		} else {
			ei <- try(eigen(sim))
			if (inherits(ei, 'try-error')){
				navals
			} else {
				zapsmall(ei[['values']])
			}
		}
		vals <- matrix(vals, ncol=3)
		
		vals
	}
	s1s2s3 <- t(apply(S, 1, .get_eig))
	colnames(s1s2s3) <- c('PS1','PS2','PS3')
	cbind(s11=s1, s12=s12, s22=s2, s1s2s3)
}

MSMS <- function(x, ...) UseMethod('MSMS')

MSMS.LT.diffusive <- MSMS.LT.uniform <- function(x, ...){
	PS <- planestress_principal_stresses(x)
	S1 <- PS[,'PS1']
	S3 <- PS[,'PS3']
	maxshear <- (S1 - S3) / 2
	meanstress <- (S1 + S3) / 2
	cbind(PS, tau_max = maxshear, Smean = meanstress)
}

fault_stresses <- function(x, theta=0, ...) UseMethod('fault_stresses')

fault_stresses.LT.diffusive <- fault_stresses.LT.uniform <- function(x, Beta=45, is.deg=TRUE, ...){
	# Beta is the angle between the plane and the greatest principal stress (S1, S3 < S2 < S1)
	if (is.deg) Beta <- Beta * pi / 180
	msms <- MSMS(x)
	twopsi <- pi - 2*Beta
	tau <- msms[,'tau_max']
	Tau <- tau * sin(twopsi)
	Sig <- msms[,'Smean'] + tau * cos(twopsi)
	return(cbind(msms, Tau = Tau, Snorm = Sig))
}

# (cfs_isoporo.defaultin beeler)
cfs_isoporo.LT.diffusive <- cfs_isoporo.LT.uniform <- function(x, Beta.deg, fric, Skemp, verbose=TRUE, ...){
	if (verbose) message('CFS {S1-to-fault angle = ', Beta.deg, '°, friction = ', fric, "} for isotropic undrained poroelastic response {B = ", Skemp, "}")
	fs <- fault_stresses(x, Beta=Beta.deg, is.deg=TRUE)
	Tau <- fs[,'Tau']
	Snorm <- fs[,'Snorm']
	Smean <- fs[,'Smean']
	cfsip <- cfs_isoporo.default(tau=Tau, mu=fric, Snormal = Snorm, Smean = Smean, B=Skemp)
	zapsmall(cfsip)
}

