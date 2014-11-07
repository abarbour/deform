x <- 1:3
z <- c(10,20,30,40)
grd <- as.matrix(expand.grid(VarX=x,VarZ=z))

testfun <- function(XY){
  sum(XY)
}

Res <- apply(grd, 1, testfun)

matrix(Res, ncol=length(z))
