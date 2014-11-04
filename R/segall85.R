#' Surface deformation associated with fluid withdrawl
#' @name segall85
#' @export
#' @examples
#' \dontrun{
#' y. <- c(-10:-2, seq(-1.9,1.9,by=0.1), 2:10)
#' su <- surface_displacement(y.*1e3, z_src=1e3)
#' plot(uz ~ y, su, type="b")
#' }
segall85 <- function(){
}

#  spatial coordinates in segall
#   x is positive downwards
#   y is positive away from source
#   z is along line source
#  source coordinates:
#   z_src is zeta (depth)
#   y_src is xi (distance away)
#  parameters:
#   mu shear modulus
#   C units of force, proportional to the source strength
#' @rdname segall85
#' @export
surface_displacement <- function(y, C.=1, mu.=1e9, ...){
  sg <- .surface_g(y, ...)
  c. <- C./mu.
  sg$uz <- sg$gz * c.
  sg$uy <- sg$gy * c.
  return(sg)
}

#' @rdname segall85
#' @export
.surface_g <- function(y=0, z_src=0, y_src=0, nuu=1/3){
  # segall85 C9
  gz <- -2*(1 - nuu)*z_src/(z_src^2 + (y - y_src)^2)
  gy <-  2*(1 - nuu)*(y - y_src)/(z_src^2 + (y - y_src)^2)
  data.frame(y, gz, gy)
}

# Time varying deformation associatedw with fluid extraction

#' @rdname segall85
#' @export
timevarying_surface_displacement <- function(y, Time, Vdot., B., L., D., HD., nuu.=1/3, y_src.=0){
  # segall85 eq 26
  # at each time slice, calculate a profile
  .tvsd <- function(yi, ti, .Vdot, .B, .L, .D, .HD, .nuu, .y_src){
    message(paste(.Vdot, .B, .L, .D, .HD, .nuu, .y_src))
    # source
    mod <- 2*.B*(1 + .nuu)*.Vdot*.D*sqrt(ti/.HD)/(3*pi*.L)
    # line source correction
    sc <- if (.y_src==0){
      1
    } else {
      Time.src <- sqrt(.y_src^2 / 4 / .HD / ti)
      err <- qnorm(Time.src/2, lower = FALSE)/sqrt(2)  #ierfc
      sum(err/(.D^2 + (yi - .y_src)^2), na.rm=TRUE)
    }
    return(mod*sc)
  }
  outer(X=y, Y=Time, FUN=.tvsd, .Vdot=Vdot., .B=B., .L=L., .D=D., .HD=HD., .nuu=nuu., .y_src=y_src.)
}
