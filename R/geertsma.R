#' Surface deformation associated with fluid withdrawl from a symmetric, rectangular reservoir
#' @name geertsma73
#' @aliases geertsma1973
#' @export
#' @param help logical; load documentation for \code{\link{geertsma77}}
#' @seealso \code{\link{vasco98}}; \code{\link{Simple-deformation}}
#' @export
#' @examples 
#' 
#' # Reproduce xxx
#' r <- seq(-10,10,by=0.2) # km
#' r.m <- r*1e3
#' xx <- surface_displacement_pressuredisk(r.m, depth=2e3, thickness=100, radius=1e3, delP.=1e3)
#' plot(uz ~ x, xx, type='l')
#' 
geertsma73 <- function(help=FALSE){
  cat("\nThis function is simply a placeholder. See the documentation ( `?geertsma73` ).\n")
  if (help) ?geertsma73
}

#' @rdname geertsma73
#' @export
.surface_g_pdisk <- function(x=0, x_src=0, z_src=0, nuu=1/3, mu., cm){
  if (missing(mu.)){
    message('assuming mu = 10GPa')
    mu. <- 10e9
  }
  #E <- 2*mu.*(1 + nuu) = youngs modulus
  if (missing(cm)) cm <- (1 - 2*nuu) / ((1 - nuu)*2*mu.)
  ##
  sc <- -cm * (1 - nuu)  / pi
  xrel <- (x - x_src)
  g <- sc * z_src / (xrel^2 + z_src^2)**(3/2)
  data.frame(x, g, xz = x/z_src)
}

#' @rdname geertsma73
#' @inheritParams surface_displacement_reservoir
#' @param depth numeric; the depth below the surface to the source
#' @param thickness numeric;
#' @param radius numeric;
#' @param delP. numeric; uniform change in pore-fluid pressure (positive = increase)
#' @param mu. numeric; the shear modulus
#' @param ... additional arguments passed to \code{\link{.surface_g_pdisk}}
surface_displacement_pressuredisk <- function(x, depth, thickness, radius, delP., shear.mod=5e9, ...){
  x <- as.vector(x)
  depth <- depth[1]
  if (depth <= 0) stop('source must be below the surface (positive depth)')
  delP. <- delP.[1]
  # Source reservoir volume
  Vol. <- thickness * pi * radius **2
  sg <- deform::.surface_g_pdisk(x, z_src = depth, mu.=shear.mod, ...)
  sg <- plyr::mutate(sg, uz = Vol. * delP. * g)
  return(sg)               
}
