library(deform)
library(TeachingDemos)

x. <- sort(unique(c(-7:-3, seq(-3.90,3.90,by=0.15), 3:7)))
su <- surface_displacement(x.*1e3, C.=1e13, z_src=0.7e3)

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
  try(my.symbols(x=su$x, y=su$uxz.mag, 
             ms.arrows, 
             #r=(su$uxz.mag),
             adj=0, col="red", inches=0.8, add=TRUE, angle=su$uxz.ang*pi/180))
  lines(ztilt*1e2 ~ x, sut, type="h", lwd=5, col="lightgreen"); text(1.1e3, 2, "Tilt = 2*dUz/dx", col="dark green", pos=4)
  lines(ztilt*1e2 ~ x, sut, lwd=3,col="dark green")
  lines(dXdx*1e2 ~ x, sue, type="h", lwd=4, col="grey60"); text(-1e2, 2, "Eee = dUx/dx", pos=2)
  lines(dXdx*1e2 ~ x, sue, lwd=3)
}
#F2()

# Fig 8
#y. <- 0:10
t. <- c(0.1,0.5,1:5)
Vdot <- 2e6
D. <- 1e3
L. <- 10e3
B. <- 0.6
c. <- 0.1
zz2 <- timevarying_surface_displacement(x., t., Vdot, B., L., D., c., x_src.=1)

F3 <- function(){
  matplot(x., zz2, type="l")
  matplot(t., t(zz2), type="l")
}

try(F3())

