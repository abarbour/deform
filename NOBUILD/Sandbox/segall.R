library(deform)
library(kook)
library(TeachingDemos)

.x. <- sort(unique(c(-7:-3, seq(-3.90,3.90,by=0.15), 3:7)))
su <- surface_displacement(.x.*1e3, C.=1e13, z_src=0.7e3)

sut <- with(su, Tilt(x, z=uz))
sue <- with(su, Uniaxial_extension(x, X=ux))

F1 <- function(){
  plot(c(NA,diff(ux)/diff(x),NA) ~ c(NA,x), su, type="s")
  abline(v=0,h=0,col="grey") 
}
#F1()

F2 <- function(){
  plot(uz ~ x, su, col=NA, ylim=c(-1,1)*20, xlim=c(-1,1)*7*1e3)
  abline(v=0,h=0,col="grey")
  lines(uz ~ x, su, type="l", pch=16, cex=1, lwd=2, col="grey"); text(2e3, -5, "U_z")
  lines(ux ~ x, su, type="l", pch=16, col="blue", cex=1, lwd=2); text(0, 6, "U_x", col="blue", pos=2)
  points(uxz.mag ~ x, su, col="red", pch=16, cex=0.6); text(-3e3, 8, "|U_xz|", col="red")
  suppressWarnings(my.symbols(x=su$x, y=su$uxz.mag, 
             ms.arrows, 
             #r=(su$uxz.mag),
             adj=0, col="red", inches=0.8, add=TRUE, angle=su$uxz.ang*pi/180))
  lines(ztilt*1e2 ~ x, sut, type="h", lwd=5, col="lightgreen"); text(1.1e3, 2, "Tilt = 2*dUz/dx", col="dark green", pos=4)
  lines(ztilt*1e2 ~ x, sut, lwd=3,col="dark green")
  lines(dXdx*1e2 ~ x, sue, type="h", lwd=4, col="grey60"); text(-1e2, 2, "Eee = dUx/dx", pos=2)
  lines(dXdx*1e2 ~ x, sue, lwd=3)
}
#F2()

# Figs 7,8
mxx <- 50
.x.km. <- sort(unique(c((-1*mxx):-3, seq(-2.90,2.90,by=0.1), 3:mxx)))
.z.km. <- sort(unique(c(seq(0,3,by=0.25),seq(3,12,by=0.75))))
yr <- 365*86400
.time. <- seq(2,10,by=2)*10*yr
.Vdot. <- 2e6/yr # volume rate m^3/yr to m^3/s
.D. <- 1e3 # depth of burial
.L. <- 10e3 # length (Vdot/L is the average rate of fluid extraction per unit length)
.B. <- 0.6 # Skemptons coeff
.c. <- 0.1 # hydraulic diffusivity m^2/s
.Sources.x. <- 1e3*c(0)
.TwoSources.x. <- 1e3*c(0,20)
# for mass computation
.t. <- 100 # thickness
.phi. <- 0.2 #porosity
# for pressure computation
.mu. <- 5.6 #GPa -- shear modulus

zz2 <- timevarying_surface_displacement(.x.km.*1e3, .time., .Vdot., .B., .L., .D., .c., Pt.Sources.x=.Sources.x.)
zz2t <- apply(zz2, 2, function(.z.) matrix(Tilt(.x.km.*1e3, z=.z.)$ztilt))

F3 <- function(){
  #matplot(.time./yr, t(zz2)*1e3, type="l", main="Subsidence, mm, Segall 1985, Fig 8B")
  matplot(.x.km., zz2*1e3, type="l", col="black", main="Subsidence, mm, Segall 1985, Fig 8B", sub=Sys.time())
}

F3t <- function(){
  #matplot(.time./yr, t(zz2t), type="l")
  matplot(.x.km., zz2t, type="l", col="black")
}

#try(F3())
#try(F3t())


zz3 <- timevarying_fluidmass(.x.km.*1e3, .time., .Vdot., .L., .t., .c., phi.=.phi.)

F4 <- function(){
  #matplot(.time./yr, t(zz3)*1e2, type="l")
  matplot(.x.km., zz3*1e2, type="l", col="black")
}

#try(F4())

zzp <- timevarying_porepressure(.x.km.*1e3, .z.km.*1e3, .time., .Vdot., .B., .L., .D., .c., .t., .mu., Pt.Sources.x=.TwoSources.x.)

F5 <- function(do.log=FALSE){
  #matplot(.time./yr, t(zz3)*1e2, type="l")
  X<- zzp[,,length(.time.)]
  if (do.log) X <- log10(abs(X))
  matplot(x=.x.km., X, col=NA)
  aaply(zzp, 3, .fun = function(X) {
    if (do.log) X <- log10(abs(X))
    matplot(x=.x.km., X, type="l", add=TRUE)
    return("x")
    })
  invisible()
}

try(F5())


F5c <- function(){
  #matplot(.time./yr, t(zz3)*1e2, type="l")
  #contour(x=.x.km., y=.z.km., zzp[,,length(.time.)], col=NA)
  aaply(zzp, 3, .fun = function(X) {
    image(x=.x.km., y=.z.km., X, ylim=c(12,0), col = brewerRamp())
    contour(x=.x.km., y=.z.km., X, ylim=c(12,0), add = TRUE)
    return("x")
    })
  invisible()
}

#try(F5c())
