#' Calculate hydraulic diffusivity of a fluid saturated medium
#' @inheritParams segall85
#' @param k. numeric; the permeability
#' @param eta. numeric; the XXX
#' @param phi. numeric; the XXX
#' @param Beta. numeric; the XXX
#' @export
hydraulic_diffusivity <- function(k., eta., mu., B., nu.=1/4, nuu.=1/3){
  # segall85 eq 11
  n1 <- 1 - nu.
  n1u <- 1 - nuu.
  n2 <- 1 - 2*nu.
  k. * B.^2 * (2*mu.*n1/n2) * (n2 * (1 + nuu.)^2 / n1u / (nuu. - nu.)) / eta. / 9
}
#' @rdname hydraulic_diffusivity
#' @export
rigid_hydraulic_diffusivity <- function(k., eta., phi., Beta.){
  # segall85 eq 12
  k./eta./phi./Beta.
}