library(deform)
library(TeachingDemos)

x. <- c(-10:-2, seq(-1.9,1.9,by=0.1), 2:10)
su <- surface_displacement(x.*1e3, C.=1e13, z_src=0.7e3)

plot(uz ~ x, su, type="b", ylim=c(-1,1)*15); text(4e3, -3, "U_z")
lines(ux ~ x, su, type="b", col="red"); text(-1e3, 3, "U_x", col="red")
lines(uxz.mag ~ x, su, type="b", col="blue", pch=16, cex=0.6); text(-3e3, 8, "|U_xz|", col="blue")
my.symbols(x=su$x, y=su$uxz.mag, ms.arrows, angle=su$uxz.ang, adj=0, col="blue", inches=0.4, add=TRUE)

plot(c(NA,diff(ux)/diff(x),NA) ~ c(NA,x), su, type="s"); abline(v=0,h=0,col="grey")

# Fig 8
#y. <- 0:10
t. <- c(0.1,0.5,1:5)
Vdot <- 2e6
D. <- 1e3
L. <- 10e3
B. <- 0.6
c. <- 0.1
zz2 <- timevarying_surface_displacement(x., t., Vdot, B., L., D., c., x_src.=1)
#matplot(x., zz2, type="l")
