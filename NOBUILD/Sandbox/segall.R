library(deform)

y. <- c(-10:-2, seq(-1.9,1.9,by=0.1), 2:10)
su <- surface_displacement(y.*1e3, C.=1e13, z_src=1e3)
#plot(uz ~ y, su, type="b", ylim=c(-2,1)*10)
#lines(uy ~ y, su, type="b", col="red")
#lines((uy+uz) ~ y, su, type="b", col="blue")

# Fig 8
#y. <- 0:10
t. <- c(0.1,0.5,1:5)
Vdot <- 2e6
D. <- 1e3
L. <- 10e3
B. <- 0.6
c. <- 0.1
matplot(y.,t.,
  timevarying_surface_displacement(y., t., Vdot, B., L., D., c., y_src.=1),
  type="l")
