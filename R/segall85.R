#' Surface deformation associated with fluid withdrawl
#' @name segall85
#' @export
#' @examples
#' \dontrun{
#' x. <- c(-10:-2, seq(-1.9,1.9,by=0.1), 2:10)
#' su <- surface_displacement(x.*1e3, z_src=1e3)
#' plot(uz ~ x, su, type="b")
#' }
segall85 <- function(){
}

#  spatial coordinates in segall
#   x is positive downwards (here z)
#   y is positive away from source (here x)
#   z is along line source (here y)
#  source coordinates:
#   z_src is zeta (depth)
#   y_src is xi (distance away)
#  parameters:
#   mu shear modulus
#   C units of force, proportional to the source strength
#' @rdname segall85
#' @export
surface_displacement <- function(x, C.=1, mu.=1e9, ...){
  sg <- .surface_g(x, ...)
  c. <- C./mu.
  sg <- mutate(sg,
               ux = gx * c.,
               uz = gz * c.,
               uxz.mag = sqrt(ux^2 + uz^2),
               uxz.ang = atan2(uz,ux) * 180/pi )
  return(sg)
}

#' @rdname segall85
#' @export
.surface_g <- function(x=0, x_src=0, z_src=0, nuu=1/3){
  # segall85 C9
  gx <-  2*(1 - nuu)*(x - x_src)/(z_src^2 + (x - x_src)^2)
  gz <- -2*(1 - nuu)*z_src/(z_src^2 + (x - x_src)^2)
  data.frame(x, gx, gz)
}

# Time varying deformation associatedw with fluid extraction

#' @rdname segall85
#' @export
timevarying_surface_displacement <- function(x, Time, Vdot., B., L., D., HD., nuu.=1/3, x_src.=0){
  # segall85 eq 26
  # at each time slice, calculate a profile
  .tvsd <- function(xi, ti, .Vdot, .B, .L, .D, .HD, .nuu, .x_src){
    message(paste(.Vdot, .B, .L, .D, .HD, .nuu, .x_src))
    # source
    mod <- 2*.B*(1 + .nuu)*.Vdot*.D*sqrt(ti/.HD)/(3*pi*.L)
    # line source correction
    sc <- if (.x_src==0){
      1
    } else {
      Time.src <- sqrt(.x_src^2 / 4 / .HD / ti)
      err <- qnorm(Time.src/2, lower = FALSE)/sqrt(2)  #ierfc
      sum(err/(.D^2 + (xi - .x_src)^2), na.rm=TRUE)
    }
    return(mod*sc)
  }
  outer(X=x, Y=Time, FUN=.tvsd, .Vdot=Vdot., .B=B., .L=L., .D=D., .HD=HD., .nuu=nuu., .x_src=x_src.)
}
