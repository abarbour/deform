pp_isotropic <- function(Smean, B){ 
	stopifnot(all((B < 1) & (B >= 0)))
	B * Smean
}

effective_normal_stress <- function(Snormal, pp) Snormal - pp

cfs_basic <- function(tau, mu, Snormal, pp){
	sn_eff <- effective_normal_stress(Snormal, pp)
	tau - mu*(sn_eff)
}
cfs_isoporo <- function(tau, mu, Snormal, Smean, B){
	pp.iso <- pp_isotropic(Smean, B)
	cfs_basic(tau, mu, Snormal, pp.iso)
}

apparent_friction_basic <- function(mu, pp, Snormal){
	mu * (1 - pp / Snormal)
}

isoporo_apparent_friction <- function(mu, B, Smean, Snormal){
	pp.iso <- pp_isotropic(Smean, B)
	apparent_friction_basic(mu, pp.iso, Snormal)
}
constant_apparent_friction <- function(mu, B){
	isoporo_apparent_friction(mu, B, 1, 1)
}
psi_angle <- function(Beta, is.deg=TRUE){
	if (is.deg) Beta <- Beta * pi / 180
	pi / 2 - Beta
}
beta_angle <- function(psi, fric, is.deg=TRUE){
	if (missing(fric)){
		psi_angle(psi, is.deg)
	} else {
		pi /4 + atan(fric)/2
	}
}

faultShearNorm_principal_stress <- function(S1, S3, Beta.deg){
	# Beta is the angle between the plane and the greatest principal stress (S1, S3 < S2 < S1)
	# (these are changes in S1 S2 S3)
	psi <- psi_angle(Beta.deg, is.deg=TRUE)
	shear <- (S1 - S3) * sin(2*psi) / 2
	normal <- ((S1 + S3) + (S1 - S3) * cos(2*psi)) / 2
	ret <- cbind(Shear = shear, Normal = normal)
	attr(ret, 'psi') <- psi
	ret
}

stresses_on_reverse_fault <- function(S1, nu_u, Beta.deg){
	# Assumes change in S3 (overburden) is zero
	# and plane strain, i.e., S2 = nu_u * S1, nu_u is undrained poissons ratio
	SN <- faultShearNorm_principal_stress(S1, S3=0, Beta.deg)
	psi <- attr(SN, 'psi') #psi_angle(Beta.deg, is.deg=TRUE)
	shear <- S1 * sin(2*psi) / 2
	Snormal <- S1 * (1 + cos(2*psi)) / 2
	stopifnot(all.equal(SN, cbind(shear, Snormal), check.attributes = FALSE))
	Smean <- S1 * (1 + nu_u) / 3
	ret <- cbind(SN, Mean = Smean)
	attr(ret, 'psi') <- psi
	ret
}

.cfs_rev_flt <- function(s2psi, c2psi, mu, B, nu_u){
	message('reverse cfs -- mu: {', mu, "} B: {", B, "} nu_u: {", nu_u, "} ")
	s1cfs_iso <- s2psi/2 - mu*(1 + c2psi)/2 + mu*B*(1 + nu_u)/3
	s1cfs_caf <- s2psi/2 - mu*(1 - B)*(1 + c2psi)/2
	cbind(Poroiso = s1cfs_iso, CAF=s1cfs_caf)
}

cfs_reverse <- function(Beta.deg, mu, B, nu_u, S1=1){
	#
	# 2D cfs on a reverse fault
	#
	# if S1 == 1 the result is effectively the cfs/S1 ratio -- equations 13a/b
	psi <- psi_angle(Beta.deg, is.deg=TRUE)
	s2psi <- sin(2*psi)
	c2psi <- cos(2*psi)
	
	IC <- .cfs_rev_flt(s2psi, c2psi, mu, B, nu_u)
	
	ret <- S1 * IC
	attr(ret, 'psi') <- psi
	ret
}

cfs_reverse_optimal <- function(mu, B, nu_u, S1=1){
	#
	# 2D cfs on a reverse fault that's optimally oriented for failure
	#
	s2psi_opt <- 1 / sqrt(mu^2 + 1)
	c2psi_opt <- mu / sqrt(mu^2 + 1) # beeler has -1 * [] here, which gives inconsistent results
		
	IC <- .cfs_rev_flt(s2psi_opt, c2psi_opt, mu, B, nu_u)
	
	ret <- S1 * IC
	psi <- asin(s2psi_opt)/2
	attr(ret, 'psi') <- psi
	ret
}

.cfs_ss_flt <- function(s2psi, c2psi, mu, B, nu_u){
	message('ss cfs -- mu: {', mu, "} B: {", B, "} nu_u: {", nu_u, "} ")
	s1cfs_iso <- (1 - nu_u)*s2psi/2 - mu*(1 + nu_u + (1 - nu_u)*c2psi)/2 + mu*B*(1 + nu_u)/3
	s1cfs_caf <- (1 - nu_u)*s2psi/2 - mu*(1 - B)*(1 + nu_u + (1 - nu_u)*c2psi)/2
	cbind(Poroiso = s1cfs_iso, CAF=s1cfs_caf)
}

cfs_strikeslip <- function(Beta.deg, mu, B, nu_u, S1=1){
	#
	# 2D cfs on a strike slip fault
	#
	# if S1 == 1 the result is effectively the cfs/S1 ratio -- equations 13a/b
	psi <- psi_angle(Beta.deg, is.deg=TRUE)
	s2psi <- sin(2*psi)
	c2psi <- cos(2*psi)
	
	IC <- .cfs_ss_flt(s2psi, c2psi, mu, B, nu_u)
	
	ret <- S1 * IC
	attr(ret, 'psi') <- psi
	ret
}


stresses_on_reverse_fault(10, 0.33, 30)

fric <- 0.65
skemp <- 0.6
poiss <- 0.33
ropt <- cfs_reverse_optimal(fric, skemp, poiss)
psifric <- asin(1 / sqrt(fric^2 + 1)) / 2
stopifnot(all.equal(psifric, attr(ropt, 'psi')))

betafric <- beta_angle(psifric, is.deg=FALSE)
#print(c(psifric, betafric)*180/pi)
r <- cfs_reverse(betafric*180/pi, fric, skemp, poiss)

stopifnot(all.equal(ropt,r))

rs <- cfs_strikeslip(betafric*180/pi, fric, skemp, poiss)

rs
