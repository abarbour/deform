library(pbapply)
#pboptions(type = "timer")
library(viridis)


source('LambertCodes.R')
source('beeler_cfs.R')

#-------------------------------------------------------------- Beel tests

stresses_on_reverse_fault(10, 0.33, 30)

fric <- 0.65
opt_ang <- atan(1/fric)*180/pi/2
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



#-------------------------------------------------------------- Lamb tests

n. <- 101
h. <- seq(-5,5,length.out=n.) * 1e3
v. <- h.

n.t <- ceiling(1.1*n.)
t. <- 10**seq(-1,2,length.out=n.t)

td <- 1 * 1e3
hw <- td
thk <- hw/2

ltd <- LambertTsai2020_Diffusive(x1=0 * 1e3, x2=1 * 1e3, y1=hw, resTopDepth=td, resThickness=thk, resDiffusivity=1, Time=t.)
ltu <- LambertTsai2020_Uniform(x1=h., x2=1 * 1e3, resTopDepth=td, resThickness=thk, resHalfWidth=hw)

#> colnames(ltu[['Def']])
#[1] "U1"  "U2"  "S11" "S12" "S22" "E11"
with(ltd, matplot(pars$Time, Def, type='l', lty=c(1,2,1,3,5,1,3,1), col=c(1,1,2,2,2,3,3,4)))
with(ltu, matplot(pars$x1, Def, ylim=c(-3,3), type='l', lty=c(1,2,1,3,5,1,3,1), col=c(1,1,2,2,2,3,3,4)))
abline(v=hw*c(-1,1), col='grey')

#library(magrittr)

examp <- function(){

	fs <- fault_stresses(ltu, Beta=30, is.deg=TRUE)
	cu.nf <- cfs_isoporo(ltu, Beta.deg=30, fric=fric, Skemp=skemp)
	cu.ss <- cfs_isoporo(ltu, Beta.deg=5, fric=fric, Skemp=skemp)

	matplot(h., abs(fs[,7:10]), type='l', log='xy', ylim=range(abs(c(cu.ss,cu.nf)), na.rm=TRUE))
	lines(h., abs(cu.nf), lwd=1, type='h', col=ifelse(cu<0,4,2))
	lines(h., abs(cu.nf), lwd=1.2)
	lines(h., abs(cu.ss), lty=2, lwd=2)

	plot(h., cu.nf / cu.ss, log='')
}

# do by depth
do_calc_by_depth <- function(d){
	ltud <- LambertTsai2020_Uniform(x1=h., x2=d, resTopDepth=td, resThickness=thk, resHalfWidth=hw)
}

d. <- h.[h. >= 0]
d. <- d.[d. < max(d.)]
Du <- lapply(d., do_calc_by_depth)
#str(Du,1)

recomp <- FALSE

if (!exists('cu.o') | recomp){
cu.o <- pbapply::pbsapply(Du, cfs_isoporo, 
	Beta.deg=opt_ang, fric=fric, Skemp=skemp, 
	Shearmod = 10e9, nu_u = poiss, 
	Qt = resQt(m_dot=1e6/(12*30*86400), resRadius=hw, timespan=12*30*86400),
	verbose=FALSE)
}

cu <- cu.o / 1e6 # to [Mk]Pa

#cu <- sapply(Du, mean_stress)
#cu <- sapply(Du, max_shear)
#cu <- sapply(Du, max_shear_principal)
#cu <- sapply(Du, ext1) # OK
#cu <- sapply(Du, ext2) # OK
#cu <- sapply(Du, shear) # out of wack for negative distances -- added kluge
#str(cu,1)

#cu[!is.finite(cu)] <- NA
norm_by <- 1 #td

xy <- list(x = h./norm_by, y = d./norm_by)
Cu <- c(xy, list(z = cu))
laCu <- c(xy, list(z = log10(abs(cu))))
Cupos <- c(xy, list(z = 1*(cu>0)))
lCup <- c(xy, list(z = log10(cu)))
lCun <- c(xy, list(z = log10(-cu)))

fields::image.plot(Cu, asp=1, ylim=rev(range(Cu[['y']], na.rm=TRUE)),
	col = c(NA,viridis::cividis(64, direction=-1)), 
	zlim=max(abs(cu),na.rm=TRUE)*c(0,1))
#image(Cupos, add=TRUE, col=adjustcolor(c('white',NA), alpha=0.5))
#contour(lCup, add=TRUE)
#contour(lCun, lty=3, add=TRUE)
#contour(Cu, add=TRUE)
contour(laCu, add=TRUE, levels=seq(-3,3))
contour(laCu, add=TRUE, levels=seq(-3,3)+log10(2), lty=2)
rect(-hw/norm_by, (td + thk)/norm_by, hw/norm_by, td/norm_by, col='white', border='grey')
