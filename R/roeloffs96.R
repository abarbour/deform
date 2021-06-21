#' Earthquake-related poroelastic solutions in Roeloffs 1996
#' 
#' @details \code{r96_drainage_pressure} refers to the solution
#' given in equation 50 of Roeloffs (1996), the pore pressure
#' response to a static strain. The temporal change in pore
#' pressure represents water table drainage as a function of time.
#' 
#' @param Times numeric; times to calculate at \code{depth} [s]
#' @param depth numeric; observation depth [m]
#' @param depth_water numeric; depth to water table [m]
#' @param diffusiv numeric; hydraulic diffusivity [m^2/s]
#' @param strain0 numeric; static strain change [unscaled strain]
#' @param B numeric; Skempton's coefficient B [-]
#' @param Ku numeric; undrained bulk modulus [Pa]
#' @param ... additional parameters
#' 
#' @export
#' 
#' @references Roeloffs, E. (1996),
#' Poroelastic techniques in the study of earthquake-related hydrologic phenomena,
#' \emph{Advances in Geophysics}, 37, 135-195, \url{https://doi.org/10.1016/S0065-2687(08)60270-8}
#' 
#' @examples 
#' # Remake Fig 16 in Roeloffs (1996)
#' d2s <- 86400
#' lti <- 10^seq(-3,0, length.out=21)
#' ti <- unique(sort(c(0, lti, seq(-10, 50, length.out=201)) * d2s))
#' zsq4c <- c(1, 24, 24*30) * d2s / 24 # timescales in hours
#' z <- 1
#' vert_c <- z^2 / 4 / zsq4c # vertical diffusivity
#' r1  <- r96_drainage_pressure(Times=ti, depth=z, diffusiv=vert_c[1]) # 1 hour
#' r2  <- r96_drainage_pressure(Times=ti, depth=z, diffusiv=vert_c[2]) # 1 day
#' r3  <- r96_drainage_pressure(Times=ti, depth=z, diffusiv=vert_c[3]) # 30 days
#' 
#' \dontrun{
#' plot(c(-10, 50), c(-0.2,1.2), col=NA, xaxs='i', yaxs='i', xlab='TIME, IN DAYS', ylab='P(t)/BKu')
#' lines(ti/d2s, r1, lty=2)
#' lines(ti/d2s, r2, lty=3)
#' lines(ti/d2s, r3)
#' text(c(20,22,32), c(0.12,0.35,0.9), parse(text=sprintf("z^2/4*c == %s", c('1 ~ "hour"','1 ~ "day"','30 ~ "days"'))))
#' 
#' plot(range(10**lti), c(5e-2, 1), col=NA, xaxs='i', log='y', xlab='TIME, IN DAYS', ylab='P(t)/BKu')
#' lines(ti/d2s, r1, lty=2)
#' lines(ti/d2s, r2, lty=3)
#' lines(ti/d2s, r3)
#' }
r96_drainage_pressure <- function(Times, depth, depth_water=0, diffusiv=1, strain0=-1, B=1, Ku=1, ...){
  # Equation 50
  neg <- Times < 0
  Times[neg] <- NA
  zsq <- (depth - depth_water)^2
  evol <- sqrt(zsq / 4 / diffusiv / Times)
  resp <- erf_re(evol)
  scaling <- -1 * B * Ku * strain0
  result <- resp * scaling
  result[neg] <- 0
  return(result)
}