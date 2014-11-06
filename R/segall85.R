#' Surface deformation associated with fluid withdrawl
#' @name segall85
#' @export
#' @seealso \code{\link{Simple-deformation}}
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

#' @rdname segall85
#' @export
timevarying_fluidmass <- function(x, Time, Vdot., L., t., HD., phi.){
  .tvmd <- function(xi, ti, .Vdot, .L, .t, .HD, .phi){
    -1 * .Vdot * sqrt(ti/.HD) * ierfc2(xi^2/4/.HD/.t) / .L / .t / .phi
  }
  outer(X=x, Y=Time, FUN=.tvmd, .Vdot=Vdot., .L=L., .t=t., .HD=HD., .phi=phi.)
}

# Time varying deformation associated with fluid extraction

#' @rdname segall85
#' @export
timevarying_surface_displacement <- function(x, Time, Vdot., B., L., D., HD., nuu.=1/3, x_src.=0){
  # segall85 eq 26
  # at each time slice, calculate a profile
  .tvsd <- function(xi, .t, .Vdot, .B, .L, .D, .HD, .nuu, .x_src){
    message(paste(.Vdot, .B, .L, .D, .HD, .nuu, .x_src))
    # source
    tt <- .t
    xx <- xi
    #
    mod <- 2 * .B * (1 + .nuu) * .Vdot * .D * sqrt(tt / .HD) / (3 * pi * .L)
    # line source correction
    Time.src <- sqrt(.x_src^2 / 4 / .HD / tt)
    err <- ierfc2(Time.src)
    #print(err)
    sc <- (err / (.D^2 + (xx - .x_src)^2)) #sum(, na.rm=TRUE)
    return(mod*sc)
  }
  t(apply(X=matrix(x), MARGIN=1, FUN=.tvsd, .t=Time, .Vdot=Vdot., .B=B., .L=L., .D=D., .HD=HD., .nuu=nuu., .x_src=x_src.))
  #outer(X=x, Y=Time, FUN=.tvsd, .Vdot=Vdot., .B=B., .L=L., .D=D., .HD=HD., .nuu=nuu., .x_src=x_src.)
}

#' Simple numerical deformation estimates
#' @name Simple-deformation
#' @param x numeric; the spatial vector in real units (not an index!)
#' @param X numeric; the quantity varying in the direction of \code{x}
#' @param y numeric; the spatial vector perpendicular to \code{x} in real units
#' @param z numeric; values which vary in the direction perpendicular to the \code{x}-\code{y} plane
#' @param left logical; should padding be on the "left" side of the vector (i.e., index number 1)?
#' @export
#' @seealso \code{\link{segall85}}
#' @examples
#' # Some x values
#' xval. <- sort(unique(c(-7:-3, seq(-3.90,3.90,by=0.15), 3:7)))
#' 
#' #  surface displacements
#' su <- surface_displacement(xval.*1e3, C.=1e13, z_src=0.7e3)
#' 
#' # Vertical tilt -- assumes axial symmetry
#' sut <- with(su, Tilt(x, z=uz))
#' # including anisotropic effects
#' sut <- with(su, Tilt(x, x=x+rnorm(length(x)), z=uz))
#' 
#' # Calculate strain in the 'x' direction
#' sue <- with(su, Uniaxial_extension(x, X=ux))
#' 
#' # see how the 'left' argument affects things:
#' .setleft(1,TRUE) # NA  1
#' .setleft(1,FALSE) # 1  NA
Uniaxial_extension <- function(x, X, left=TRUE){
  deriv <- diff(X)/diff(x)
  dXdx <- .setleft(deriv, left)
  data.frame(x, dXdx)
}
#' @rdname Simple-deformation
#' @export
Tilt <- function(x, y=NULL, z, left=TRUE){
  deriv1 <- diff(z)/diff(x)
  dzdx <- .setleft(deriv1, left)
  dzdy <- if (is.null(y)){
    dzdx
  } else {
    # this needs to be more intelligent:
    # what if z is a matrix!?
    deriv2 <- diff(z)/diff(x)
    .setleft(deriv2, left)
  }
  tlt <- dzdx + dzdy
  #ang <- 90 - tlt * 180 /pi
  data.frame(x=x, ztilt = tlt)
}
#' @rdname Simple-deformation
#' @export
.setleft <- function(x, left=TRUE){
  if (left){
    c(NA, x)
  } else {
    c(x, NA)
  }
}