#' Surface deformation associated with fluid withdrawl from a disk-shaped reservoir
#' @name segall89
#' @aliases segall1989
#' @export
#' @param help logical; load documentation for \code{\link{segall89}}
#' @seealso \code{\link{segall85}}; \code{\link{Simple-deformation}}
#' @export
#' @examples 
#' 
#' # Reproduce Figure 5 from Segall (1989)
#' r <- seq(-4,4,by=0.1)
#' xx <- surface_displacement_reservoir(r)
#' plot(Xi_z ~ x, xx, type='l', ylim=c(-1.6,0.8))
#' lines(Xi_x ~ x, xx)
#' lines(Xi_xx ~ x, xx)
#' 
segall89 <- function(help=FALSE){
  cat("\nThis function is simply a placeholder. See the documentation ( `?segall89` ).\n")
  if (help) ?segall89
}

#' @rdname segall89
#' @export
.surface_g_xi <- function(x, depth, halfwidth, thickness){
  # segall89 7abc
  xi_n <- (x - halfwidth)/depth
  xi_p <- (x + halfwidth)/depth
  # 7b: radial displacements (fixing a typo, and flipping ratio)
  a <- 1 + xi_n**2
  b <- 1 + xi_p**2
  Xi_x <- log10(a / b)
  Xi_x_T <- Xi_x * thickness
  # 7a: vertical displacements
  Xi_z <- atan(xi_n) - atan(xi_p)
  Xi_z_T <- Xi_z * thickness
  # 7c: radial extension (flipping difference)
  Xi_xx <- (xi_n/(1 + xi_n**2)) - (xi_p/(1 + xi_p**2))
  Xi_xx_T <- Xi_xx * thickness / depth
  data.frame(x,  xi_n, xi_p,  Xi_x, Xi_z, Xi_xx,  Xi_x_T, Xi_z_T, Xi_xx_T)
}

#' @rdname segall89
#' @export
#' @param x numeric; the horizontal (or radial) distance from the source
#' @param depth numeric; the depth below the surface to the top of the reservoir
#' @param halfwidth numeric; the half-length of the reservoir
#' @param thickness numeric; the thickness of the reservoir
#' @param DelM. numeric; the change in mass per unit volume in the reservoir
#' @param nuu. numeric; the undrained Poisson's ratio for the reservoir material
#' @param B. numeric; Skempton's ratio for the reservoir material
#' @param rho_f. numeric; the density of the fluid in the reservoir material
#' @param ... additional parameters (currently unused)
surface_displacement_reservoir <- function(x, depth=1, halfwidth=depth, thickness=depth/10, DelM.=1, nuu.=1/3, B.=1, rho_f.=1000, ...){
  sxi <- .surface_g_xi(x, depth, halfwidth, thickness)
  c. <- 2 * (1 + nuu.) * B. * DelM. / (3 * pi * rho_f.)
  sxi <- plyr::mutate(sxi,
               ux = Xi_x_T * c.,
               uz = Xi_z_T * c.,
               exx = Xi_xx_T * c.)
  return(sxi)
}
