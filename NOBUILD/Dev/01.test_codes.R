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

n. <- 101
h. <- seq(-10,10,length.out=n.)
v. <- h.
t. <- 10**seq(-1,2,length.out=n.-1)

hw <- max(h.)/3
td <- 1
thk <- td/10

ltd <- LambertTsai2020_Diffusive(x1=0, x2=1, y1=hw, resTopDepth=td, resThickness=thk, resDiffusivity=1, Time=t.)
ltu <- LambertTsai2020_Uniform(x1=h., x2=1, resTopDepth=td, resThickness=thk, resHalfWidth=hw)

with(ltd, matplot(pars$Time, Def, type='l', lty=c(1,2,1,3,5,1,3,1), col=c(1,1,2,2,2,3,3,4)))
with(ltu, matplot(pars$x1, Def, ylim=c(-3,3), type='l', lty=c(1,2,1,3,5,1,3,1), col=c(1,1,2,2,2,3,3,4)))
abline(v=hw*c(-1,1), col='grey')


#plane stress principal stresses
PSPS(ltd)
PSPS(ltu)
