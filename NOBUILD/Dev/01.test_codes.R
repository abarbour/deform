source('LambertCodes.R')
source('beeler_cfs.R')

#-------------------------------------------------------------- Beel tests

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



#-------------------------------------------------------------- Lamb tests

n. <- 251
h. <- seq(-10,10,length.out=n.)
v. <- h.

n.t <- ceiling(1.1*n.)
t. <- 10**seq(-1,2,length.out=n.t)

hw <- max(h.)/5
td <- 1
thk <- td/2

ltd <- LambertTsai2020_Diffusive(x1=0, x2=1, y1=hw, resTopDepth=td, resThickness=thk, resDiffusivity=1, Time=t.)
ltu <- LambertTsai2020_Uniform(x1=h., x2=1, resTopDepth=td, resThickness=thk, resHalfWidth=hw)

#> colnames(ltu[['Def']])
#[1] "U1"  "U2"  "S11" "S12" "S22" "E11"
with(ltd, matplot(pars$Time, Def, type='l', lty=c(1,2,1,3,5,1,3,1), col=c(1,1,2,2,2,3,3,4)))
with(ltu, matplot(pars$x1, Def, ylim=c(-3,3), type='l', lty=c(1,2,1,3,5,1,3,1), col=c(1,1,2,2,2,3,3,4)))
abline(v=hw*c(-1,1), col='grey')


#plane stress principal stresses
PSPS(ltd)
PSPS(ltu)

rbind(
fault_normal_stress(ltd),
fault_shear_stress(ltd),
mean_stress(ltd),
max_shear_principal(ltd),
max_shear(ltd)
)

rbind(
fault_normal_stress(ltu),
fault_shear_stress(ltu),
mean_stress(ltu),
max_shear_principal(ltu),
max_shear(ltu)
)

cd0 <- cfs_isoporo(ltd)
plot(t., cd0, type='l', ylim=max(abs(cd0), na.rm=TRUE)*c(-1,1))
lines(t., cfs_isoporo(ltd, theta=45), lty=5)
lines(t., cfs_isoporo(ltd, theta=90), lty=2)
cu0 <- cfs_isoporo(ltu)
plot(h., cu0, type='l', ylim=max(abs(cu0), na.rm=TRUE)*c(-1,1))
lines(h., cfs_isoporo(ltu, theta=45), lty=5)
lines(h., cfs_isoporo(ltu, theta=90), lty=2)

# do by depth
do_calc_by_depth <- function(d){
	ltud <- LambertTsai2020_Uniform(x1=h., x2=d, resTopDepth=td, resThickness=thk, resHalfWidth=hw)
}

d. <- h.[h. >= 0]
d. <- d.[d. < max(d.)/1.5]
Du <- lapply(d., do_calc_by_depth)

#str(Du,1)

cu <- sapply(Du, cfs_isoporo, theta=60, fric=fric, Skemp=skemp) #, verbose=FALSE)
#cu <- sapply(Du, mean_stress)
#cu <- sapply(Du, max_shear)
#cu <- sapply(Du, ext1) # OK
#cu <- sapply(Du, ext2) # OK
#cu <- sapply(Du, shear) # out of wack for negative distances -- added kluge
#str(cu,1)

xy <- list(x = h./td, y = d./td)
Cu <- c(xy, list(z = cu))
Cupos <- c(xy, list(z = 1*(cu>0)))
lCup <- c(xy, list(z = log10(cu)))
lCun <- c(xy, list(z = log10(-cu)))

fields::image.plot(Cu, asp=1, ylim=rev(range(Cu[['y']])), zlim=max(abs(cu),na.rm=TRUE)*c(-1,1))
image(Cupos, add=TRUE, col=adjustcolor(c('white',NA), alpha=0.5))
#contour(lCup, add=TRUE)
#contour(lCun, lty=3, add=TRUE)
contour(Cu, add=TRUE)
rect(-hw, (td + thk)/td, hw, td/td, col='white', border='grey')
