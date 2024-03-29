#' Rice and Cleary solutions
#'
#' @details \code{\link{rc76_edge_disloc_poro}} is the 2D solution
#' to an edge dislocation (e.g., fault) embedded in a poroelastic medium
#' 
#' @param Times numeric; times to calculate at \code{depth} [s]
#' @param slip numeric; relative slip along the edge dislocation [m]
#' @param radial_dist numeric; distance from xxx to xxx [m]
#' @param flt_angle numeric; angle counter clockwise from fault strike
#' @param diffusiv numeric; in plane hydraulic diffusivity [m^2/s]
#' @param nuu numeric; undrained Poisson's ratio [-]
#' @param nu numeric; drained Poisson's ratio [-]
#' @param B numeric; Skempton's coefficient [-]
#' @param mu numeric; elastic shear modulus [Pa]
#' @param ... additional arguments
#'
#' @return tibble
#' @export
#'
#' @examples
#' 
#' # Pore pressure transient associated with near-field fault slip
#' 
#' ti <- c(-100, -1, 0, 10**seq(-4, log10(6598), length.out=301))
#' Mw <- 4.33
#' M0 <- Mw2M0(Mw) #3.935501e+15 N*m for Mw 4.33
#' bx <- 10^(-3.22 + 0.69 * Mw) #0.5
#' rd <- 250
#' mu <- 2e8
#' Diffu <- 60
#' rc <- rc76_edge_disloc_poro(ti, bx, rd, diffusiv=Diffu, mu=mu)
#' nrc <- rc76_edge_disloc_poro(ti, bx, -rd, diffusiv=Diffu, mu=mu)
#' all.equal(rc$pp, -nrc$pp)
#'
#' \dontrun{
#' plot(pp / 1e3 + 5 ~ I(Times / 3600), rc, type='l', lwd=2, ylab='kPa', 
#'      ylim=c(-8.8,36.9), yaxs='i', xaxs='i')
#' axis(4)
#' abline(v=1*60/3600, lty=3)
#' rect(0,0,100,100)
#' SrcArea <- M0 / (mu * bx)
#' SrcRad <- sqrt(SrcArea / pi)
#' print(c(Area=SrcArea, Rad=SrcRad))
#' }
#' 
rc76_edge_disloc_poro <- function(Times, slip, radial_dist, flt_angle=30, 
                                       diffusiv=1, 
                                       nuu=0.33, nu=0.25, B=0.6,
                                       mu=30e9,
                                       ...){
  # Equation X of Rice and Cleary 1976, or equations 7.149 - 7.152 of Wang (2000)
  data.table::inrange(radial_dist, 0, Inf)
  stopifnot(length(slip) == 1)
  stopifnot(length(radial_dist) == 1)
  stopifnot(length(flt_angle) == 1)
  
  zer <- Times <= 0
  calc_times <- Times
  calc_times[zer] <- NA
  
  theta <- flt_angle * pi / 180
  r <- radial_dist
  geom_fac <- slip / r
  nu_fac <- (nuu - nu) / (1 - nu) / (1 - nuu)
  
  # scaling and temporal changes for all quantities
  common_scaling <- (mu / 2 / pi) * nu_fac * geom_fac 
  rsqtimes <- r^2 / 4 / diffusiv / calc_times
  response <- 1 - exp(-1 * rsqtimes)
  
  # pore pressure
  eta <- poroelastic_stress_coefficient(nu, nuu, B)
  #print(eta)
  pp <- common_scaling * (sin(theta) / eta) * response
  
  # stresses
  srr <- common_scaling * sin(theta) * (-(1 - nu)/(nuu - nu) + rsqtimes * response)
  srt <- common_scaling * cos(theta) * ((1 - nu)/(nuu - nu) - rsqtimes * response)
  stt <- common_scaling * sin(theta) * (2*exp(-1 * rsqtimes) - (1 - nu)/(nuu - nu) - rsqtimes * response)
  
  result <- tibble::tibble(Times=Times, sigma_tt = stt, sigma_rt = srt, sigma_rr = srr, pp = pp)
  
  return(result)
}
