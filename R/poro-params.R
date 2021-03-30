#' Poroelastic parameters
#' 
#' @details 
#' \code{\link{effstress}} Wang (2000), equation 2.52 
#' \code{\link{poroelastic_stress_coefficient}} Wang (2000), equation 3.47
#' 
#' @param K numeric; the bulk modulus
#' @param Ks numeric; the constituent bulk modulus (e.g., one from an unjacketed lab test)
#' @param nu numeric; the drained Poisson's ratio  [0,0.5]
#' @param nuu numeric; the undrained Poisson's ratio  [0,0.5] (\code{nu} < \code{nuu})
#' @param B numeric; Skempton's coefficient [0,1]
#' @param mu numeric; the isotropic elastic shear modulus (keep at unity for 'normalized')
#' @param diffusiv numeric; the isotropic hydraulic diffusivity (keep at unity for 'normalized')
#' 
#' @references Wang, H.F. (2000)
#' 
#' @aliases poroelastic-parameters

#' @export
effstress <- function(nu, nuu, B){
  # e.g., Wang (2000) equation 2.52
  3*(nuu - nu) / (1 - 2*nu) / (1 + nuu) / B
}
#' @rdname effstress
#' @export
poroelastic_stress_coefficient <- function(nu, nuu, B){
  alpha <- effstress(nu, nuu, B)
  (1 - 2*nu) * alpha / 2 / (1 - nu)
}
#' @rdname effstress
#' @export
biot_willis_effstress <- function(K, Ks){
  1 - K / Ks
}
#' @rdname effstress
#' @export
biot_compressibility <- function(nu, nuu, B, mu=1){
  9 * (1 - 2*nuu) * (nuu - nu) / 2 / mu / (1 - 2*nu) / (1 + nuu)^2 / B^2
}
#' @rdname effstress
#' @export
darcy_conductivity <- function(nu, nuu, B, mu=1, diffusiv=1){
  9 * diffusiv * (1 - nuu) * (nuu - nu) / 2 / mu / (1 - nu) / (1 + nuu)^2 / B^2
}
#' @rdname effstress
#' @export
lame_first <- function(nu, mu=1){
  2 * mu * nu / (1 - 2*nu)
}