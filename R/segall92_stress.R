library(tidyverse)
library(fields)
library(pracma)
library(cubature)
library(viridis)

.stress_G_segall92 <- function(radial.position, depth.position, 
	source.radius=1.001, source.depth=1.001, nu.=0.25, return.tensor=FALSE, verbose=TRUE){
	##
	## Calculate Green's functions for stresses (i.e., equations 59 - 62)
	##
	#
	if (verbose) message(sprintf("G at r = %s, z = %s ...", radial.position, depth.position))
	
		
	.Bess <- function(k, n=1, m=0, tee=0, cee=0, a=source.radius, r=radial.position, dist.sc=1){
		if (dist.sc){
			k <- k * dist.sc
			r <- r * dist.sc
			a <- a * dist.sc
		}
		#
		# Bessel functions for integration (rhs of eq 58, inside integral)
		#	
		#Jn <- Bessel::BesselJ
		#Jm <- Bessel::BesselJ
		Jn <- Bessel::BesselJ(k * a, n)
		Jm <- Bessel::BesselJ(k * r, m)
		Jn * Jm * (k ^ tee) * exp(-1 * cee * k)
	}
	
	.Inmt_integ <- function(n., m., tee., cee., nsub=4000L, tol_pow, upper){
		#
		# Integration to get I^c(n,m;t) (equation 58)
		#
		# default is ** 0.25: ~1e-4
		tol <- .Machine$double.eps ** ifelse(missing(tol_pow), 0.4, tol_pow)
		if (missing(upper)) upper <- Inf
		
		res <- try(stats::integrate(.Bess, n=n., m=m., tee=tee., cee=cee., 
			lower=0, upper=upper, subdivisions = nsub, rel.tol=tol, stop.on.error=FALSE))
		
		is.bad <- inherits(res, 'try-error')
		if (!is.bad){
			if (!res[['message']]=="OK"){
			#print(c(n., m., tee., cee., val))
				NA
			} else {
				res[['value']]
			}
		} else {
			NA
		}
	}
	
	.Inmt_trapz <- function(n., m., tee., cee., nrad = 20){
		require(pracma)
		
		
		#x <- seq(0, source.radius*nrad, length.out=nrad*41)
		#y <- .Bess(x, n=n., m=m., tee=tee., cee=cee.)
		#res <- pracma::trapz(x, y)
		#res #[['value']]
		
		res <- pracma::trapzfun(.Bess, a = 0, b = source.radius * nrad, n=n., m=m., tee=tee., cee=cee.) 
		res[['value']]
	}
	
	.Inmt_cube <- function(n., m., tee., cee.){
		require(cubature)
		
		res <- try(cubature::cubintegrate(.Bess,  
			lower=0, upper=100*source.radius, method = "vegas",
			n=n., m=m., tee=tee., cee=cee.)
		)
		
		is.bad <- inherits(res, 'try-error')
		if (!is.bad){
			res[['integral']]
		} else {
			NA
		}
	}
		
	#.Inmt <- .Inmt_trapz
	.Inmt <- .Inmt_integ
	#.Inmt <- .Inmt_cube
	
	zdiff <- depth.position - source.depth
	zsum <- depth.position + source.depth
	
	#epsilon: -1 for z > d, and 1 for z < d
	epsilon <- ifelse(zdiff >= 0, -1, 1)
	espdiff <- -1*epsilon*zdiff
	
	nu.c.1 <- (3 + 4 * nu.)
	nu.c.2 <- (3 - 4 * nu.)
	two.z <- 2 * depth.position
	
	#message("---> integrating Bessel functions...")
	I1 <- .Inmt(1,0,1,espdiff)	
	I2 <- .Inmt(1,0,1,   zsum) * nu.c.1
	I3 <- .Inmt(1,0,2,   zsum) * two.z
	
	I4 <- .Inmt(1,2,1,espdiff)
	I5 <- .Inmt(1,2,1,   zsum) * nu.c.2
	I6 <- .Inmt(1,2,2,   zsum) * two.z
	
	I7 <- .Inmt(1,1,1,espdiff)
	I8 <- .Inmt(1,1,1,   zsum)
	I9 <- .Inmt(1,1,2,   zsum) * two.z

	
	#message("---> calculating radial stress...")
	#Grr --  equation 59: (rho/2)*{ I-eps(z-d)(1,0;1) + (3+4v)I(z+d)(1,0;1) - 2zI(z+d)(1,0;2) - I-eps(z-d)(1,2;1) - (3-4v)I(z+d)(1,2;1) + 2zI(z+d)(1,2;2) }
	Grr <- (I1 + I2 - I3 - I4 - I5 + I6) * source.radius / 2
	
	#message("---> calculating hoop stress...")
	#Goo -- equation 60: (rho/2)*{ I-eps(z-d)(1,0;1) + (3+4v)I(z+d)(1,0;1) - 2zI(z+d)(1,0;2) + I-eps(z-d)(1,2;1) + (3-4v)I(z+d)(1,2;1) - 2zI(z+d)(1,2;2) }
	Gtt <- (I1 + I2 - I3 + I4 + I5 - I6) * source.radius / 2
	
	#message("---> calculating vertical stress...")
	#Gzz -- equation 61: -rho*{ I-eps(z-d)(1,0;1) - I(z+d)(1,0;1) - 2zI(z+d)(1,0;2)}
	Gzz <- (I1 - I2 / nu.c.1 - I3) * source.radius * -1
	
	#message("---> calculating radial/vertical shear stress...")
	#Grz -- equation 62: rho*{ eps*I-eps(z-d)(1,1;1) - I(z+d)(1,1;1) -2zI(z+d)(1,1;2)}
	Grz <- (I7 - I8 - I9) * source.radius
	
	# Stresses (extension positive)
	g0 <- c(Grr=Grr, Gtt=Gtt, Gzz=Gzz)
	if (return.tensor){
		G <- diag(g0)
		colnames(G) <- rownames(G) <- c('r','t','z')
		# stress tensor is symmetric, add shear strains
		G[1,3] <- Grz
		G[3,1] <- Grz
	} else {
		G <- cbind(data.frame(r=radial.position, z=depth.position), as.data.frame(t(g0)), data.frame(Grz=Grz))
	}
	
	return(G)
	
}

gaussian_reservoir_depth <- function(r, d0, A=1, l_c=8){
	# The depth to the reservoir layer d(r) is taken to be
	# (expression below)
	# where A is the amplitude of the dome, and 
	# lc is the characteristic length of the structure. 
	# In these calculations it is assumed that the thickness 
	# of the producing zone, T, is considerably less than 
	# the reservoir depth and the other characteristic 
	# dimensions of the reservoir.
	d0 - A*(exp(-(r/l_c)^2) - 1/2) # equation 63
}

uniform_declining_pressure <- function(r, z, pressure.change, source.radius, source.thickness, source.depth, inside.R){
	#The Green's functions, Gij(r,z; rho,d), correspond to a 
	# ring of dilatation at radius rho and depth d and are 
	# given by Segall [1992]
	if (missing(inside.R)) inside.R <- source.radius
	stopifnot(inside.R <= source.radius)
	
	# zero if outside the source, otherwise = pressure.change
	radial.heav <- r <= inside.R
	inside.source <- (z >= (source.depth - source.thickness/2)) | (z <= (source.depth + source.thickness/2))

	pressure.change * as.numeric(inside.source) * as.numeric(radial.heav)
	
}

constant_pressure <- function(p0=0){
	approxfun(c(0,Inf), rep(p0, 2))
}

gaussian_pressure_distribution <- function(r, times, P0=constant_pressure(5.1), r_c=8){
	# For the purposes of computation the pressure is constant_pressure(5.1)
	# taken to vary with radius according to 
	# where po(t) is the maximum pressure decline
	# at r = 0, and rc is the characteristic length 
	# of the pressure drop with radial distance.
	p0t <- if (missing(P0)){
		1
	} else {
		if (inherits(P0, 'function')){
			P0(times)
		} else {
			P0[times]
		}
	} 
	p0t * exp(-(r / r_c)^4)
}

#rs <- seq(0,20, length.out=101)
#plot(rs, gaussian_reservoir_depth(rs, d0=2, A=-2), type='l', ylim=c(8,0))
#plot(rs, gaussian_pressure_distribution(rs, times=0), type='l', ylim=c(8,0))
#plot(rs, uniform_declining_pressure(rs, z=2, 1, 2, 0.1, 2), type='l')

isotropic_stresses <- function(S) UseMethod('isotropic_stresses')
isotropic_stresses.default  <- function(S){
	Siso <- diag(S)
	return(Siso)
}

# Stress invariants
# http://www.continuummechanics.org/principalstrain.html
.sinvar_I <- function(S){
	# S11 + S22 + S33
	inv <- sum(isotropic_stresses(S))
	return(inv)
}
.sinvar_II <- function(S){
	#
	# Second: S11 S22 + S22 S33 + S33 S11 - S12^2 - S23^2 - S21^2
	#
	# S11 S22 S33
	D <- diag(S)
	#S11 <- D[1]
	#S22 <- D[2]
	#S33 <- D[3]
	# S12 S13 S23
	#off.D <- S[upper.tri(S)]
	#S12 <- off.D[1]
	#S13 <- off.D[2]
	#S23 <- off.D[3]
	#inv <- S11*S22 + S22*S33 + S33*S11 - S12^2 - S13^2 - S23^2
	#
	# II = (tr(S*S) - tr(S)^2)/2
	inv <- (sum(diag(t(S) %*% S)) - sum(D)^2) / 2
	return(inv)
}
.sinvar_III <- function(S, log.=FALSE){
	# S11 S22 S33 - S11 S23^2 - S22 S13^2 - S33 S12^2 + 2 S12 S13 S23
	inv <- det(S)
	if (log.) inv <- log10(inv)
	return(inv)
}
.sinvar <- function(S, log.=FALSE){
	iI <- .sinvar_I(S)
	iII <- .sinvar_II(S)
	liIII <- .sinvar_III(S, log.=log.)
	allinv <- c(`I`=iI, II=iII, III=liIII)
	if (log.) names(allinv)[3] <- 'logIII'
	return(allinv)
}

stress_invariant <- function(S, invariant=NULL, ...){
	invar <- match.arg(invariant, c('All','I','II','III'))
	FUN <- switch(invar, `I`=.sinvar_I, II=.sinvar_II, III=.sinvar_III, All=.sinvar)
	Si <- FUN(S, ...)
	attr(Si, 'invariant.number') <- invar
	attr(Si, 'units') <- 'Pa'
	if (invar == "All") Si <- as.data.frame(t(Si))
	return(Si)
}

principal_stresses <- function(S) {
	
	if (any(is.na(S)) | any(!is.finite(S))){
		
		return(data.frame(S1=NA, S2=NA, S3=NA))
		
	} else {
		Se <- eigen(S, symmetric = TRUE)
		Sp <- Se[['values']]
		names(Sp) <- paste0("S",1:3)
	
		# Spv <- Se[['vectors']]
		# colnames(Spv) <- paste0("S",1:3)
		# rownames(Spv) <- c('n1','n2','n3')
		# attr(Spm, 'eigen.vectors') <- Spv
	
		Sp <- as.data.frame(t(Sp))
		return(Sp)
	}
}

.calc_inv <- function(S11=1, S22=2, S33=3, S13=-1){
	S <- isotropic_stresses(c(S11, S22, S33))
	S[1,3] <- S13
	S[3,1] <- S13
	#print(S)
	Si <- stress_invariant(S)
	Sp <- principal_stresses(S)

	cbind(Sp, Si)
}

#.calc_inv()
#
#stop()

stress_field <- function(r, z, 
	source.radius=8.001, source.depth=3.501, source.thickness=0.3,
	pressure.drop=-5.1, Biot.coefficient=0.25, Poissons.ratio=0.25, verbose=FALSE, ...){

	expand.grid(r=r, z=z) %>% dplyr::arrange(desc(r), desc(z)) -> rz
	
	message("Calculate Green's functions for stresses...")
	
	rz %>% 
		dplyr::rowwise(.) %>% 
		dplyr::do(.stress_G_segall92(.$r, .$z, source.radius=source.radius, source.depth=source.depth, nu.=Poissons.ratio, verbose=verbose, ...)) -> Gij
	
	Gij %>% 
		dplyr::rowwise(.) %>% 
		dplyr::do(.calc_inv(.$Grr, .$Gtt, .$Gzz, .$Grz)) %>%
		dplyr::mutate(DiffStress = (S1 - S3), MaxShear = DiffStress/2) -> invariants
	
	dplyr::bind_cols(Gij, invariants) %>%
		tidyr::gather(., 'Stress','Value', -r, -z) -> Gij.tidy

	#print(head(Gij.tidy))
	
	stress_lvls <- unique(Gij.tidy$Stress) # if you wanted to change the order
	
	split(Gij.tidy, factor(Gij.tidy$Stress, levels=stress_lvls)) -> Gij.list
	
	Res <- list(Params=list(r=r, z=z, 
							src=list(radius=source.radius, depth=source.depth, thickness=source.thickness, p = pressure.drop), 
							alpha = Biot.coefficient, nu = Poissons.ratio), 
		Results=Gij.list)
	class(Res) <- c("stress_field_G", class(Res))
	
	return(Res)
}

plot.stress_field_G <- function(x, finer=2.2, ...){
	
	Gij <- x[['Results']]
	Params <- x[['Params']]
	src <- Params[['src']]

	source.radius <- src[['radius']]
	res <- with(Params, list(x=r, y=z))

	pp <- src[['p']]
	nu <- Params[['nu']]
	alph <- Params[['alpha']]
	
	scaling <- pp * (1 - 2*nu) / (1 - nu) / 2
	
	print(scaling)
	
	.stress_matrix <- function(S, plot.it=TRUE){
		
		stress_type <- unique(S$Stress)
		message("Working on ", stress_type, "...")
		
		Z <- tidyr::spread(dplyr::select(S, -Stress), r, Value)
		
		z <- scaling * t(as.matrix(dplyr::select(Z, -z)))
		
		ccol <- 'black'
		pal <- if (stress_type %in% c('MaxShear','DiffStress')){
			zl <- NULL #abs(scaling)*c(0,1)
			viridis::magma(32, direction=sign(-scaling))
		} else if (stress_type %in% c("II")){
			zl <- NULL #abs(scaling)*c(0,1)
			ccol <- 'white'
			viridis::plasma(32, direction=sign(scaling))
		} else if (stress_type %in% c("III")){
			zl <- NULL #abs(scaling)*c(0,1)
			viridis::plasma(32, direction=sign(scaling))
		} else {
			zl <- abs(scaling)*c(-1,1)/4
			kook::brewerRamp(num=16) #fields::tim.colors(16)
		}
		
		if (stress_type == "I"){
			z <- z/3
			stress_type <- "Mean Stress"
		}
		
		if (stress_type == "II"){
			z <- sqrt(abs(z))
			stress_type <- "J2^(1/2)"
		}
	
		if (stress_type == "III"){
			z <- abs(z)^(1/3)
			stress_type <- "J3^(1/3)"
		}
		
		ttl <- if (stress_type == "MaxShear"){
			"Maximum Shear Stress"
		} else {
			stress_type
		}
		
		res[['z']] <- zoo::na.approx(z, na.rm = FALSE)
		res[['stress']] <- stress_type

		
		if (plot.it){
			require(fields)
			
			finer <- as.integer(finer)

			if (finer > 1){
				message("upsampling by ", finer, "x...")
				finer.grid <- with(res,{
					nx <- length(x)
					rx <- range(x)
					ny <- length(y)
					ry <- range(y)
					
					list(x = seq(rx[1], rx[2], length.out=ceiling(finer*nx)), y = seq(ry[1], ry[2], length.out=ceiling(finer*ny)))
				})
				
				res <- fields::interp.surface.grid(res, finer.grid)
			}
			
			yl <- rev(range(res[['y']]))
			fields::image.plot(res, 
				#asp=1,
				main=ttl, 
				ylim=yl, ylab="Depth, km",
				zlim=zl, 
				#nlevel=8, 
				col=pal,
				xaxs='i', 
				xlim=rev(range(res[['x']])), 
				xlab="Distance West of source, km"
				)
			levs <- seq(-1,1,by=0.05)
			ltys <- as.numeric(levs < 0)*2 + 1
			contour(res, levels=levs, lty=ltys, col=ccol, labcex=1, vfont=c("sans serif", "bold"), add=TRUE)
			with(Params[['src']],{
				rect(0, depth - thickness/2, source.radius, depth + thickness/2, col='lightgrey', border='red', lty=2)
				points(source.radius, depth, pch=3, lwd=2, col='red')
			})
		}
		return(res)
	}
	
	invisible(lapply(Gij, .stress_matrix))
	
}

# 
# #.stress_G_segall92(5,5)
# 
# area.delaware <- 6843 #km^2
# radius.delaware <- sqrt(area.delaware/pi)
# depth.delaware <- 2.5 
# 
# depths <- seq(0, 10, by=0.3)
# dists <- seq(0, ceiling(2.0*radius.delaware), by=2)
# 
# redo.calc <- FALSE
# if (!exists("Scalc") | redo.calc) Scalc <- stress_field(dists, depths, source.radius=radius.delaware, source.depth=depth.delaware, source.thickness=0.3)
# 
# pdf("~/Desktop/figs/Delaware_stresses_%02d.pdf", height=3, width=7, onefile=FALSE)
# par(las=1, mar=c(3,3,2,1), tcl=-0.25, mgp=c(2,0.5,0))
# try(plot(Scalc))
# dev.off()
# 
# #save(Scalc, file="files/Scalc.rda", compress='xz')

